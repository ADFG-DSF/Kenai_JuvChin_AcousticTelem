# ── Load libraries ──
library(tidyverse)
library(lubridate)
library(data.table)
library(vroom)

# ── Define file paths ──
# input_dir  <- "/Users/dakotarygh/Desktop/JSATS Kenai Chinook Smolt Project/Code and Data/Marine Data/Raw Data"
# output_dir <- "/Users/dakotarygh/Desktop/JSATS Kenai Chinook Smolt Project/Code and Data/Marine Data/Cleaned Data"
# unfiltered_clean_dir <- file.path(output_dir, "Unfiltered Clean Files")
# taglist_path <- "/Users/dakotarygh/Desktop/JSATS Kenai Chinook Smolt Project/Code and Data/Marine Data/Tagged Fish Data/chinook_tagged.csv"

# input_dir  <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\lots_of_raw_data\\Stationary Raw Data"
# output_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\lots_of_raw_data\\Cleaned Data"
# unfiltered_clean_dir <- file.path(output_dir, "Unfiltered Clean Files")
# taglist_path <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\metadata\\chinook_tagged.csv"
# 
# stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\metadata\\StationaryMetadata.xlsx"

input_dir  <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\Raw Data"
output_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\Cleaned Data"
unfiltered_clean_dir <- file.path(output_dir, "Unfiltered Clean Files")
taglist_path <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\chinook_tagged.csv"

# stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\StationaryMetadata.csv"
stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\StationaryMetadata.xlsx"




# ── Load tagged fish IDs ──
tagcodes_to_keep <- read_csv(taglist_path, show_col_types = FALSE) %>%
  filter(FATE == "R") %>%  ## this is new
  select(TAG_ID) %>%
  distinct() %>%
  pull(TAG_ID) %>%
  as.character()

### date-time associated with release of each tag
tagcodes_release <- read_csv(taglist_path, show_col_types = FALSE) %>%
  # mutate(release_datetime = as.POSIXlt(paste(PULL_DATE, END)), format="%Y-%m-%d %H:%M:%OS") %>%
  # mutate(release_datetime = parse_date_time(paste(PULL_DT, END), orders = "Ymd HMS", tz = "UTC")) %>%
  mutate(release_datetime = END_DT) %>%
  filter(FATE == "R") %>%
  select(TAG_ID, release_datetime) %>%
  rename(TagCode=TAG_ID) %>%
  distinct(TagCode, .keep_all = TRUE)

### Load stationary metadata
StationaryMetadata <- readxl::read_xlsx(stationary_metadata_dir,
                                        sheet="Metadata"
                                        ) %>% #skip=1
  mutate(ArraySiteName_date = paste(str_replace_all(paste0(Array, SiteID), " ", ""),
                                    paste0(25, substr(as.character(DeploymentDate),6,7), 
                                           substr(as.character(DeploymentDate),9,10)), sep="_"))  %>%
  select(ArraySiteName_date, `River Kilometer`, Latitude, Longitude)   ### can keep additional columns too





# ── Read and clean receiver-reset-enabled CSV ──
read_clean_vroom <- function(file_path) {
  raw_lines <- read_lines(file_path)
  
  header_indices <- which(str_starts(raw_lines, "Internal"))
  if (length(header_indices) == 0) return(NULL)
  
  blocks <- map(header_indices, function(idx) {
    end <- ifelse(idx == tail(header_indices, 1), length(raw_lines), header_indices[which(header_indices > idx)[1]] - 1)
    lines <- raw_lines[idx:end]
    
    header_line <- lines[1]
    header_names <- str_split(header_line, ",")[[1]] %>% str_trim() %>% make.names(unique = TRUE)
    
    data_only <- lines[-1]
    data_only <- data_only[!data_only %in% lines[1]]  # Remove repeated header lines
    
    if (length(data_only) == 0) return(NULL)
    
    tryCatch({
      vroom::vroom(I(c(str_c(header_names, collapse = ","), data_only)), delim = ",", col_types = cols(.default = "c"))
    }, error = function(e) {
      message("⚠️ Failed to parse block in ", basename(file_path), ": ", e$message)
      return(NULL)
    })
  })
  
  blocks <- discard(blocks, is.null)
  if (length(blocks) == 0) return(NULL)
  
  ## is this where the warning is coming from??
  # combined <- bind_rows(blocks)
  suppressWarnings(combined <- bind_rows(blocks))
  
  combined %>%
    filter(!is.na(TagCode)) %>%  # Only drop rows with missing tag values (real errors)
    
    ### detecting Array from filename
    mutate(Array = substr(strsplit(basename(file_path), "_")[[1]][1], 1, nchar(basename(file_path)) - 18 - nchar(SiteName))) %>%
    
    ### forcing SiteName to have 3 characters
    mutate(SiteName = ifelse(nchar(SiteName)==1,
                             paste0("00", SiteName),
                             ifelse(nchar(SiteName)==2,
                                    paste0("0", SiteName), 
                                    SiteName))) %>%
    
    ### inserting site metadata
    mutate(ArraySiteName_date = paste0(Array, SiteName, substr(basename(file_path),
                                                               nchar(basename(file_path)) - 17,
                                                               nchar(basename(file_path)) - 11))) %>%
    left_join(StationaryMetadata, by="ArraySiteName_date")  %>%
    
    ### selecting a subset of columns that are of interest
    select(Array, SiteName, DateTime, TagCode, 
           Tilt, SigStr, 
           ArraySiteName_date, 
           `River Kilometer`, Latitude, Longitude)
}

### quick function to summarize number of individuals
n_indiv <- function(df) length(unique(df$TagCode))

### generalized function to print summary message
print_message <- function(df, the_message) {
  message(the_message, ": ", format(nrow(df), big.mark = ","),
          " (", format(n_indiv(df), big.mark = ","), " individuals)")
}

# ── Clean but unfiltered detection dataframe ──
clean_but_unfiltered_data <- function(df, message=FALSE) {
  initial_n <- nrow(df)
  
  df <- df %>%
    # select(any_of(c("DateTime", "TagCode", "Tilt", "SiteName"))) %>%
    mutate(
      TagCode = str_sub(TagCode, 4, 7),
      Tilt = as.numeric(Tilt),
      DateTime = parse_date_time(DateTime, orders = "mdY HMS", tz = "UTC")
    ) %>%   ############################# took this out to preserve 
    # filter(complete.cases(.))   ############################# took this out to preserve 
    left_join(tagcodes_release, by="TagCode")  ####  
  
  # after_tag_filter_n <- df %>% filter(TagCode %in% tagcodes_to_keep) %>% nrow()
  
  if(message) {
  # message("🧾 Initial rows: ", format(initial_n, big.mark = ","),
  #         " (", format(n_indiv(df), big.mark = ","), " individuals)")
    print_message(df, "✍️ Initial rows")
  }
  
  # # df <- df %>% filter(TagCode %in% tagcodes_to_keep)
  # df <- df %>% filter(DateTime >= release_datetime) %>%
  #   select(-release_datetime)
  # 
  # after_tag_filter_n <- nrow(df)
  # 
  # message("🎣 After TagCode filter: ", format(after_tag_filter_n, big.mark = ","),
  #         " (", format(n_indiv(df), big.mark = ","), " individuals)")
  
  df
}


### pre-filter to only allow individuals with n hits
apply_prefilter <- function(df, n=2, message=FALSE) {
  thetable <- table(df$TagCode)
  to_keep <- names(thetable[thetable >= n])
  df_out <- subset(df, TagCode %in% to_keep)
  
  if(message) {
    print_message(df_out, "🧹 After prefilter")
  }
  df_out
}

apply_tagcode_filter <- function(df, tagcodes_to_keep=tagcodes_to_keep, message=FALSE) {
  
  ### Version 1: check whether TagCode is in the list 
  df_out <- df %>% filter(TagCode %in% tagcodes_to_keep)
  
  # ### Version 2: check whether TagCode has been released
  # df <- df %>% filter(DateTime >= release_datetime) %>%
  #   select(-release_datetime)
  
  # after_tag_filter_n <- nrow(df)
  # 
  # if(message) {
  # message("🎣 After TagCode filter: ", format(after_tag_filter_n, big.mark = ","),
  #         " (", format(n_indiv(df), big.mark = ","), " individuals)")
  # }
  # df
  
  if(message) {
    print_message(df_out, "🎣 After TagCode filter")
  }
  df_out
}

# # ── Apply all filtering steps ──
# apply_all_filters <- function(df) {
#   after_tagcode <- df %>%
#     apply_prefilter(n=2)  %>%  #### change this as needed!!!
#     apply_tagcode_filter
#   
#   after_multipath <- after_tagcode %>%
#     apply_multipath_filter() #threshold_sec = multipath_threshold)
#   message("🎯 After multipath filter: ", format(nrow(after_multipath), big.mark = ","),
#           " (", format(n_indiv(after_multipath), big.mark = ","), " individuals)")
#   
#   after_3sec <- after_multipath %>%
#     filter_tags_with_consistent_3s_intervals()
#   message("⏱ After 3-sec ping filter: ", format(nrow(after_3sec), big.mark = ","),
#           " (", format(n_indiv(after_3sec), big.mark = ","), " individuals)")
#   
#   after_3in30 <- after_3sec %>%
#     # filter_tags_with_3_in_30s()
#     filter_tags_with_3_in_30s_v2()
#   message("🔁 After 3-in-30s burst filter: ", format(nrow(after_3in30), big.mark = ","),
#           " (", format(n_indiv(after_3in30), big.mark = ","), " individuals)")
#   
#   final <- after_3in30 %>%
#     arrange(DateTime) %>%
#     select(SiteName, everything()) %>%
#     mutate(DateTime = format(DateTime, "%m/%d/%Y %H:%M:%S")) %>%
#     select(SiteName, everything())
#   
#   message("✅ Final rows: ", format(nrow(final), big.mark = ","),
#           " (", format(n_indiv(final), big.mark = ","), " individuals)")
#   
#   final
# }

# ── Multipath threshold (seconds) ──
# multipath_threshold <- 0.3

# ── Filter rapid multipath detections ──
apply_multipath_filter <- function(df, threshold_sec=0.3, message=FALSE) {
  df_out <- df %>%
    arrange(SiteName, TagCode, DateTime) %>%
    group_by(SiteName, TagCode) %>%
    mutate(TimeDiff = as.numeric(difftime(DateTime, lag(DateTime), units = "secs"))) %>%
    filter(is.na(TimeDiff) | TimeDiff >= threshold_sec) %>%
    ungroup() %>%
    select(-TimeDiff)
  
  if(message) {
    print_message(df_out, "〽️ After multipath filter")
  }
  df_out
}

# ── Keep tags with ~3-second ping intervals ──
filter_tags_with_consistent_3s_intervals <- function(df, period=3, delta=0.5, message=FALSE) {
  df_out <- df %>%
    arrange(SiteName, TagCode, DateTime) %>%
    group_by(SiteName, TagCode) %>%
    mutate(
      TimeDiff = as.numeric(difftime(DateTime, lag(DateTime), units = "secs")),
      ValidInterval = TimeDiff >= (period-delta) & TimeDiff <= (period+delta),
      ValidCount = sum(ValidInterval, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(ValidCount >= 2) %>%
    select(-TimeDiff, -ValidInterval, -ValidCount)
  
  if(message) {
    print_message(df_out, "⏰ After interval filter")
  }
  df_out
}

# ── Filter for at least 3 detections in 30 seconds ──
filter_tags_with_3_in_30s <- function(df, n_events=3, t_sec=30, message=FALSE) {
  df_out <- df %>%
    arrange(SiteName, TagCode, DateTime) %>%
    group_by(SiteName, TagCode) %>%
    mutate(
      DetectionCount = map_int(seq_along(DateTime), function(i) {
        sum(abs(as.numeric(difftime(DateTime[i], DateTime, units = "secs"))) <= t_sec)
      })
    ) %>%
    ungroup() %>%
    filter(DetectionCount >= n_events) %>%
    select(-DetectionCount)
  
  if(message) {
    print_message(df_out, "*️⃣ After event filter")
  }
  df_out
}

# ── Filter for at least 3 detections in 30 seconds (much faster) ──
filter_tags_with_3_in_30s_v2 <- function(df, n_events=3, t_sec=30, message=FALSE) {
  
  if(nrow(df) == 0) return(df) 
  
  bytag_df <- tapply(df, df$TagCode, \(x) x)
  bytag_toreturn <- list()
  
  if(n_events > 100) n_events <- 100  # maxing this out at 100!
  for(itag in seq_along(bytag_df)) {
    
    # reorder to be sure
    bytag_df[[itag]] <- bytag_df[[itag]][order(bytag_df[[itag]]$DateTime),]
    
    # create a staggered matrix of numeric datetimes
    stagmat <- matrix(nrow=nrow(bytag_df[[itag]])+n_events-1, ncol=n_events)
    
    # populating it
    for(j in 1:n_events) {
      stagmat[(1:nrow(bytag_df[[itag]]))+j-1, j] <- as.numeric(bytag_df[[itag]]$DateTime)
    }
    
    # sum of successive differences
    sumdiffs <- apply(stagmat, 1, \(x) sum(abs(diff(x))))
    theserows <- sumdiffs <= t_sec
    theserows[is.na(theserows)] <- FALSE
    
    # reusing stagmat, maybe it saves on memory!
    stagmat <- matrix(FALSE, nrow=nrow(bytag_df[[itag]])+n_events-1, ncol=n_events)
    
    # populating it
    for(j in 1:n_events) {
      stagmat[(1:nrow(bytag_df[[itag]]))+j-1, j] <- theserows[1:nrow(bytag_df[[itag]])]
    }
    
    # which rows to use for each tag?
    # returnthese <- apply(stagmat[1:nrow(bytag_df[[itag]]),], 1, any, na.rm=TRUE)
    returnthese <- apply(stagmat[-(1:(n_events-1)),], 1, any, na.rm=TRUE)
    bytag_toreturn[[itag]] <- bytag_df[[itag]][returnthese,]
  }
  
  # smash it all into one dataframe
  # return(tibble(do.call(rbind, bytag_toreturn)))
  df_out <- do.call(rbind, bytag_toreturn)
  
  if(message) {
    print_message(df_out, "*️⃣ After event filter")
  }
  df_out
}


# ### left over from testing the new event filter function - i would like to test it further
# df0 <- df
# df <- df0
# df1 <- subset(df, TagCode=="B5A5")
# 
# {
#   tstart <- Sys.time()
#   aaaa <- filter_tags_with_3_in_30s(df1)
#   Sys.time() - tstart
# } # Time difference of 40.84626 secs
# {
#   tstart <- Sys.time()
#   bbbb <- filter_tags_with_3_in_30s_v2(df1)
#   Sys.time() - tstart
# } # Time difference of 0.6515019 secs
# dim(aaaa)
# dim(bbbb)
# dim(df1)
# 
# aaaa %>% ggplot(aes(x=DateTime, y=TagCode)) + geom_point()
# bbbb %>% ggplot(aes(x=DateTime, y=TagCode)) + geom_point()
# df1 %>% ggplot(aes(x=DateTime, y=TagCode)) + geom_point()
# 
# table(df0$TagCode)
# table(aaaa$TagCode)
# table(bbbb$TagCode)




######## interactive Shiny app to investigate processing flow
library(shiny)

server <- shinyServer(function(input, output) {
  
  options(shiny.maxRequestSize=500*1024^2)
  
  printsize <- function(df) {
    paste0(format(nrow(df), big.mark = ","),
           " rows (", format(n_indiv(df), big.mark = ","), " tags)")
  }
  
  df_in <- reactive({
    req(input$file1)
    clean_but_unfiltered_data(read_clean_vroom(input$file1$datapath))
  })
  
  df1 <- reactive({
    apply_prefilter(df_in(), n=input$n_prefilter)
  })
  
  df2 <- reactive({
    if(input$f1=="None") df2 <- df1()
    if(input$f1=="Event") df2 <- filter_tags_with_3_in_30s_v2(df1(), input$n_events, input$t_sec)
    if(input$f1=="Tag Code") df2 <- apply_tagcode_filter(df1(), tagcodes_to_keep=tagcodes_to_keep, message=FALSE)
    if(input$f1=="Multipath") df2 <- apply_multipath_filter(df1(), input$thresh)
    if(input$f1=="Interval") df2 <- filter_tags_with_consistent_3s_intervals(df1(), input$period, input$delta)
    df2
  })
  
  df3 <- reactive({
    if(input$f2=="None") df3 <- df2()
    if(input$f2=="Event") df3 <- filter_tags_with_3_in_30s_v2(df2(), input$n_events, input$t_sec)
    if(input$f2=="Tag Code") df3 <- apply_tagcode_filter(df2(), tagcodes_to_keep=tagcodes_to_keep, message=FALSE)
    if(input$f2=="Multipath") df3 <- apply_multipath_filter(df2(), input$thresh)
    if(input$f2=="Interval") df3 <- filter_tags_with_consistent_3s_intervals(df2(), input$period, input$delta)
    df3
  })
  
  df4 <- reactive({
    if(input$f3=="None") df4 <- df3()
    if(input$f3=="Event") df4 <- filter_tags_with_3_in_30s_v2(df3(), input$n_events, input$t_sec)
    if(input$f3=="Tag Code") df4 <- apply_tagcode_filter(df3(), tagcodes_to_keep=tagcodes_to_keep, message=FALSE)
    if(input$f3=="Multipath") df4 <- apply_multipath_filter(df3(), input$thresh)
    if(input$f3=="Interval") df4 <- filter_tags_with_consistent_3s_intervals(df3(), input$period, input$delta)
    df4
  })
  
  df_out <- reactive({
    if(input$f4=="None") df5 <- df4()
    if(input$f4=="Event") df5 <- filter_tags_with_3_in_30s_v2(df4(), input$n_events, input$t_sec)
    if(input$f4=="Tag Code") df5 <- apply_tagcode_filter(df4(), tagcodes_to_keep=tagcodes_to_keep, message=FALSE)
    if(input$f4=="Multipath") df5 <- apply_multipath_filter(df4(), input$thresh)
    if(input$f4=="Interval") df5 <- filter_tags_with_consistent_3s_intervals(df4(), input$period, input$delta)
    return(df5)
  })
  
  output$thePlot <- renderPlot(
    if(length(unique(df_out()$TagCode)) <= 150) {
      df_out() %>%
        mutate(`Tags in library` = factor(ifelse(TagCode %in% tagcodes_to_keep, "Yes","No"),
                                          levels = c("Yes", "No"))) %>%
        ggplot(aes(x=DateTime, y=TagCode, col=`Tags in library`)) + 
        geom_point(show.legend = TRUE) +
        scale_color_discrete(drop=FALSE) +
        theme_bw()
    }
  )
  
  output$size_in <- renderText(paste(printsize(df_in()), "   -- -- -- Initial"))
  output$size1 <- renderText(paste(printsize(df1()), "   -- -- -- after Prefilter"))
  output$size2 <- renderText(paste(printsize(df2()), "   -- -- -- after", input$f1))
  output$size3 <- renderText(paste(printsize(df3()), "   -- -- -- after", input$f2))
  output$size4 <- renderText(paste(printsize(df4()), "   -- -- -- after", input$f3))
  output$size_out <- renderText(paste(printsize(df_out()), "   -- -- -- after", input$f4))
})
ui <- shinyUI(fluidPage(
  titlePanel("Add Title Here"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose .csv File",
                accept = c(".csv")),
      sliderInput("n_prefilter", "Prefilter n (tags with at least this many hits)",
                  min=2, max=100, value=2),
      
      radioButtons("f1", "Filter 1", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      radioButtons("f2", "Filter 2", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      radioButtons("f3", "Filter 3", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      radioButtons("f4", "Filter 4", choices=c("None", "Event", "Multipath", "Tag Code", "Interval"), inline=TRUE),
      
      hr(),
      h5("Event filter:"),
      sliderInput("n_events", "Number of events:", min=1, max=100, value=3),
      sliderInput("t_sec", "in time interval (sec):", min=1, max=100, value=30),
      
      hr(),
      h5("Multipath filter:"),
      sliderInput("thresh", "Threshold:", min=0.01, max=1, value=0.3, step=0.01),
      
      hr(),
      h5("Interval filter:"),
      sliderInput("period", "Period between pings (sec):", min=1, max=100, value=3),
      sliderInput("delta", "tolerance (sec):", min=0, max=3, value=0.5, step=0.05)
      
    ),
    #beginning of main section
    mainPanel(
      textOutput("size_in"),
      textOutput("size1"),
      textOutput("size2"),
      textOutput("size3"),
      textOutput("size4"),
      textOutput("size_out"),
      
      plotOutput("thePlot", height="600px", width="700px")
    )
  )
))
shinyApp(ui = ui, server = server)






# ──----------- Main processing loop -----------──

### --- LOOP CONTROLS --- ###

# Set filtration order & parameters here!!
apply_all_filters_reorder <- function(df, message=FALSE) {
  df %>%
    apply_prefilter(n=2, 
                    message=message) %>%
    
    apply_tagcode_filter(tagcodes_to_keep=tagcodes_to_keep, 
                         message=message) %>%
    filter_tags_with_consistent_3s_intervals(period=3, delta=0.5, 
                                             message=message) %>%
    filter_tags_with_3_in_30s_v2(n_events=3, t_sec=30, 
                                 message=message) %>%
    apply_multipath_filter(threshold_sec=0.3, 
                           message=message)
}

# Whether to save output to external files, or just create summary table
save_output <- FALSE

# Whether to print summary messages of numbers of entries & individuals
message <- TRUE




### --- initializing processing loop --- ###

csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

##### taking a subset to test
# csv_files <- csv_files[1:6]

#### removing a bad file??
csv_files <- csv_files[basename(csv_files) != "Mile33078_250623_143501.csv"]

### creating a summary table for processing
the_tbl <- data.frame(file = basename(csv_files), 
                      rows_init = NA,
                      rows_unfiltered = NA,
                      rows_filtered = NA,
                      indiv_init = NA,
                      indiv_unfiltered = NA,
                      indiv_filtered = NA,
                      metadata_found = NA)

### storing a list of filtered data.frames for easy checking
all_filtered <- list()

for (file_path in csv_files) {
  
  ### this is silly but it's kinda nice to know how far it's gotten
  cat("file", which(csv_files==file_path), "of", length(csv_files),'\n')
  
  message("📄 Processing: ", basename(file_path))
  
  df <- tryCatch(read_clean_vroom(file_path), error = function(e) {
    message("⚠️ Failed to read: ", basename(file_path), ": ", e$message)
    return(NULL)
  })
  the_tbl$rows_init[which(csv_files==file_path)] <- nrow(df)
  the_tbl$indiv_init[which(csv_files==file_path)] <- n_indiv(df)
  the_tbl$metadata_found[which(csv_files==file_path)] <- all(!is.na(df$Latitude))
  if (is.null(df)) next
  
  if(save_output) {
    # 💾 Save cleaned-but-unfiltered data (before tag filtering)
    unfiltered_vroom_out <- file.path(unfiltered_clean_dir, paste0("Unfiltered_", basename(file_path)))
    fwrite(df, unfiltered_vroom_out)
    message("💾 Saved unfiltered vroom-cleaned file: ", unfiltered_vroom_out)
  }
  
  
  unfiltered <- tryCatch(clean_but_unfiltered_data(df, message=message), error = function(e) {
    message("⚠️ Unfiltered cleaning failed: ", e$message)
    return(NULL)
  })
  the_tbl$rows_unfiltered[which(csv_files==file_path)] <- nrow(unfiltered)
  the_tbl$indiv_unfiltered[which(csv_files==file_path)] <- n_indiv(unfiltered)
  if (is.null(unfiltered)) next
  
  filtered <- tryCatch(apply_all_filters_reorder(unfiltered, message=message), error = function(e) {
    message("⚠️ Filtering failed: ", e$message)
    return(NULL)
  })
  the_tbl$rows_filtered[which(csv_files==file_path)] <- nrow(filtered)
  the_tbl$indiv_filtered[which(csv_files==file_path)] <- n_indiv(filtered)
  all_filtered[[which(csv_files==file_path)]] <- filtered
  if (is.null(filtered)) next
  
  if(save_output) {
    filtered_out <- file.path(output_dir, paste0("Cleaned_", basename(file_path)))
    fwrite(filtered, filtered_out)
    message("✔ Saved: ", filtered_out)
  }
}
the_tbl
sum(the_tbl$rows_filtered) / sum(the_tbl$rows_init)  # global acceptance rate
mean(the_tbl$rows_filtered / the_tbl$rows_init)  # mean acceptance rate across files





# # ── Main processing loop ──
# csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
# 
# ##### test case
# # csv_files <- csv_files[1:6]
# 
# #### removing a bad file??
# csv_files <- csv_files[basename(csv_files) != "Mile33078_250623_143501.csv"]
# 
# ### creating a summary table for processing
# the_tbl <- data.frame(file = basename(csv_files), 
#                       rows_init = NA,
#                       rows_unfiltered = NA,
#                       rows_filtered = NA,
#                       indiv_init = NA,
#                       indiv_unfiltered = NA,
#                       indiv_filtered = NA,
#                       metadata_found = NA)
# 
# ### storing a list of filtered data.frames for easy checking
# all_filtered <- list()
# 
# for (file_path in csv_files) {
#   
#   ### this is silly but it's kinda nice to know how far it's gotten
#   cat("file", which(csv_files==file_path), "of", length(csv_files),'\n')
#   
#   message("📄 Processing: ", basename(file_path))
#   
#   df <- tryCatch(read_clean_vroom(file_path), error = function(e) {
#     message("⚠️ Failed to read: ", basename(file_path), ": ", e$message)
#     return(NULL)
#   })
#   the_tbl$rows_init[which(csv_files==file_path)] <- nrow(df)
#   the_tbl$indiv_init[which(csv_files==file_path)] <- n_indiv(df)
#   the_tbl$metadata_found[which(csv_files==file_path)] <- all(!is.na(df$Latitude))
#   if (is.null(df)) next
#   
#   # 💾 Save cleaned-but-unfiltered data (before tag filtering)
#   unfiltered_vroom_out <- file.path(unfiltered_clean_dir, paste0("Unfiltered_", basename(file_path)))
#   fwrite(df, unfiltered_vroom_out)
#   message("💾 Saved unfiltered vroom-cleaned file: ", unfiltered_vroom_out)
#   
#   
#   unfiltered <- tryCatch(clean_but_unfiltered_data(df), error = function(e) {
#     message("⚠️ Unfiltered cleaning failed: ", e$message)
#     return(NULL)
#   })
#   the_tbl$rows_unfiltered[which(csv_files==file_path)] <- nrow(unfiltered)
#   the_tbl$indiv_unfiltered[which(csv_files==file_path)] <- n_indiv(unfiltered)
#   if (is.null(unfiltered)) next
#   
#   filtered <- tryCatch(apply_all_filters(unfiltered), error = function(e) {
#     message("⚠️ Filtering failed: ", e$message)
#     return(NULL)
#   })
#   the_tbl$rows_filtered[which(csv_files==file_path)] <- nrow(filtered)
#   the_tbl$indiv_filtered[which(csv_files==file_path)] <- n_indiv(filtered)
#   all_filtered[[which(csv_files==file_path)]] <- filtered
#   if (is.null(filtered)) next
#   
#   filtered_out <- file.path(output_dir, paste0("Cleaned_", basename(file_path)))
#   fwrite(filtered, filtered_out)
#   message("✔ Saved: ", filtered_out)
# }
# the_tbl



# # was SiteName something weird?
# the_tbl$SiteName <- sapply(all_filtered, \(x) x$SiteName[1])
# the_tbl

## which receiver data files were not found in the metadata?
receiver_files <- sapply(strsplit(the_tbl$file,"_"), \(x) paste(x[1:2], collapse="_"))

# these were found
receiver_files[receiver_files %in% StationaryMetadata$ArraySiteName_date] %>% 
  sort

# these were NOT found
receiver_files[!(receiver_files %in% StationaryMetadata$ArraySiteName_date)] %>% 
  sort

# these metadata entries had associated files
StationaryMetadata$ArraySiteName_date[StationaryMetadata$ArraySiteName_date %in% receiver_files] %>% 
  sort

# these metadata entries did NOT have associated files
StationaryMetadata$ArraySiteName_date[!(StationaryMetadata$ArraySiteName_date %in% receiver_files)] %>% 
  sort



# smashing all data together into one data.frame
all_filtered_df <- do.call(rbind, all_filtered)
### this should ultimately be saved as an external .csv file

write_csv(all_filtered_df, file=paste0(output_dir, "\\all_filtered.csv"))



