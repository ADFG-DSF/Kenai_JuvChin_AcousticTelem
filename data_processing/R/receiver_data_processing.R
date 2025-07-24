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
  select(TAG_ID) %>%
  distinct() %>%
  pull(TAG_ID) %>%
  as.character()

### date-time associated with release of each tag
tagcodes_release <- read_csv(taglist_path, show_col_types = FALSE) %>%
  # mutate(release_datetime = as.POSIXlt(paste(PULL_DATE, END)), format="%Y-%m-%d %H:%M:%OS") %>%
  mutate(release_datetime = parse_date_time(paste(PULL_DATE, END), orders = "Ymd HMS", tz = "UTC")) %>%
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


# StationaryMetadata <- read_csv(stationary_metadata_dir)
# # add leading zeroes to sitename
# # date to Date to character to substr?? or pull date elements from Date



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
    # mutate(ArraySiteName_date = paste(strsplit(basename(file_path), "_")[[1]][1:2], collapse ="_")) %>% 
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

# ── Clean but unfiltered detection dataframe ──
clean_but_unfiltered_data <- function(df) {
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
  
  message("🧾 Initial rows: ", format(initial_n, big.mark = ","),
          " (", format(n_indiv(df), big.mark = ","), " individuals)")
  
  # df <- df %>% filter(TagCode %in% tagcodes_to_keep)
  df <- df %>% filter(DateTime >= release_datetime) %>%
    select(-release_datetime)
  
  after_tag_filter_n <- nrow(df)
  
  message("🎣 After TagCode filter: ", format(after_tag_filter_n, big.mark = ","),
          " (", format(n_indiv(df), big.mark = ","), " individuals)")
  
  df
}

# ── Apply all filtering steps ──
apply_all_filters <- function(df) {
  after_multipath <- df %>%
    apply_multipath_filter(threshold_sec = multipath_threshold)
  message("🎯 After multipath filter: ", format(nrow(after_multipath), big.mark = ","),
          " (", format(n_indiv(after_multipath), big.mark = ","), " individuals)")
  
  after_3sec <- after_multipath %>%
    filter_tags_with_consistent_3s_intervals()
  message("⏱ After 3-sec ping filter: ", format(nrow(after_3sec), big.mark = ","),
          " (", format(n_indiv(after_3sec), big.mark = ","), " individuals)")
  
  after_3in30 <- after_3sec %>%
    filter_tags_with_3_in_30s()
  message("🔁 After 3-in-30s burst filter: ", format(nrow(after_3in30), big.mark = ","),
          " (", format(n_indiv(after_3in30), big.mark = ","), " individuals)")
  
  final <- after_3in30 %>%
    arrange(DateTime) %>%
    select(SiteName, everything()) %>%
    mutate(DateTime = format(DateTime, "%m/%d/%Y %H:%M:%S")) %>%
    select(SiteName, everything())
  
  message("✅ Final rows: ", format(nrow(final), big.mark = ","),
          " (", format(n_indiv(final), big.mark = ","), " individuals)")
  
  final
}

# ── Multipath threshold (seconds) ──
multipath_threshold <- 0.3

# ── Filter rapid multipath detections ──
apply_multipath_filter <- function(df, threshold_sec) {
  df %>%
    arrange(SiteName, TagCode, DateTime) %>%
    group_by(SiteName, TagCode) %>%
    mutate(TimeDiff = as.numeric(difftime(DateTime, lag(DateTime), units = "secs"))) %>%
    filter(is.na(TimeDiff) | TimeDiff >= threshold_sec) %>%
    ungroup() %>%
    select(-TimeDiff)
}

# ── Keep tags with ~3-second ping intervals ──
filter_tags_with_consistent_3s_intervals <- function(df) {
  df %>%
    arrange(SiteName, TagCode, DateTime) %>%
    group_by(SiteName, TagCode) %>%
    mutate(
      TimeDiff = as.numeric(difftime(DateTime, lag(DateTime), units = "secs")),
      ValidInterval = TimeDiff >= 2.5 & TimeDiff <= 3.5,
      ValidCount = sum(ValidInterval, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    filter(ValidCount >= 2) %>%
    select(-TimeDiff, -ValidInterval, -ValidCount)
}

# ── Filter for at least 3 detections in 30 seconds ──
filter_tags_with_3_in_30s <- function(df) {
  df %>%
    arrange(SiteName, TagCode, DateTime) %>%
    group_by(SiteName, TagCode) %>%
    mutate(
      DetectionCount = map_int(seq_along(DateTime), function(i) {
        sum(abs(as.numeric(difftime(DateTime[i], DateTime, units = "secs"))) <= 30)
      })
    ) %>%
    ungroup() %>%
    filter(DetectionCount >= 3) %>%
    select(-DetectionCount)
}

# ── Main processing loop ──
csv_files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

##### test case
# csv_files <- csv_files[1:6]

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
  
  # 💾 Save cleaned-but-unfiltered data (before tag filtering)
  unfiltered_vroom_out <- file.path(unfiltered_clean_dir, paste0("Unfiltered_", basename(file_path)))
  fwrite(df, unfiltered_vroom_out)
  message("💾 Saved unfiltered vroom-cleaned file: ", unfiltered_vroom_out)
  
  
  unfiltered <- tryCatch(clean_but_unfiltered_data(df), error = function(e) {
    message("⚠️ Unfiltered cleaning failed: ", e$message)
    return(NULL)
  })
  the_tbl$rows_unfiltered[which(csv_files==file_path)] <- nrow(unfiltered)
  the_tbl$indiv_unfiltered[which(csv_files==file_path)] <- n_indiv(unfiltered)
  if (is.null(unfiltered)) next
  
  filtered <- tryCatch(apply_all_filters(unfiltered), error = function(e) {
    message("⚠️ Filtering failed: ", e$message)
    return(NULL)
  })
  the_tbl$rows_filtered[which(csv_files==file_path)] <- nrow(filtered)
  the_tbl$indiv_filtered[which(csv_files==file_path)] <- n_indiv(filtered)
  all_filtered[[which(csv_files==file_path)]] <- filtered
  if (is.null(filtered)) next
  
  filtered_out <- file.path(output_dir, paste0("Cleaned_", basename(file_path)))
  fwrite(filtered, filtered_out)
  message("✔ Saved: ", filtered_out)
}
the_tbl



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


### summarizing passage date for each individual

# orig metadata
metadata_rivermile <- readxl::read_xlsx(stationary_metadata_dir,
                                        sheet="Metadata") %>%
  mutate(Array = str_replace_all(Array, "[ ']", "")) %>%
  mutate(Array = ifelse(Array %in% c("MarineNorth", "MarineSouth"), "Marine", Array)) %>% #### take this out as needed
  select(Array, `River Mile`, `River Kilometer`) %>%
  group_by(Array) %>%
  summarise(rivermile = mean(`River Mile`, na.rm=TRUE),
            river_km = mean(`River Kilometer`, na.rm=TRUE)) %>%
  mutate(Array_order = paste(str_pad(rank(rivermile, ties.method="first"), 2, "left", "0"), Array)) %>%
  mutate(Array_order_rev = paste(str_pad(1+nrow(.)-rank(rivermile, ties.method="first"), 2, "left", "0"), Array))

# summarizing as the DateTime of greatest signal strength, per TagCode & Array
all_filtered_passage <- all_filtered_df %>%
  mutate(Array = ifelse(Array %in% c("MarineNorth", "MarineSouth"), "Marine", Array)) %>% #### take this out as needed
  group_by(TagCode, Array) %>%
  summarise(passage = DateTime[which.max(SigStr)]) %>%
  mutate(Array=ifelse(Array=="Mile19SB", "Mile19", Array)) %>%
  mutate(Array=ifelse(Array=="KeyesPropertySB", "KeyesProperty", Array)) %>%
  left_join(metadata_rivermile) %>%
  mutate(passage=as.POSIXlt(passage, format="%m/%d/%Y %H:%M:%OS"))
### this should ultimately be saved as an external .csv file



# plotting!!!!!!
library(patchwork)

plot1 <- all_filtered_passage %>%
  ggplot(aes(x=passage, y=Array_order, colour=TagCode, group=TagCode)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  # xlab("Passage Date (max sig strength)") +
  xlab("") +
  ylab("Array (ordered)")

plot2 <- all_filtered_passage %>%
  ggplot(aes(x=passage, y=river_km, colour=TagCode, group=TagCode)) +
  geom_point() +
  geom_line() +
  theme_bw()  +
  xlab("Passage Date (max sig strength)") +
  ylab("River km")

plot1 / plot2 + plot_layout(guides="collect")



## making survival matrix
surv_mat_tibble <- pivot_wider(all_filtered_passage,
                        id_cols=TagCode,
                        names_from = Array_order_rev,
                        values_from = passage)
surv_mat <- surv_mat_tibble[,order(colnames(surv_mat_tibble))] %>%
  (\(x) x[, colnames(x) != "TagCode"]) %>%
  as.matrix %>%
  is.na %>% `!`
rownames(surv_mat) <- surv_mat_tibble$TagCode
  
