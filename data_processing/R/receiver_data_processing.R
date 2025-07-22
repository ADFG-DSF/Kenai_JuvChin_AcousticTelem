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
input_dir  <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\lots_of_raw_data\\Stationary Raw Data"
output_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\lots_of_raw_data\\Cleaned Data"
unfiltered_clean_dir <- file.path(output_dir, "Unfiltered Clean Files")
taglist_path <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\metadata\\chinook_tagged.csv"

stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\JSATS_processing_experiments\\metadata\\StationaryMetadata.xlsx"

# ── Load tagged fish IDs ──
tagcodes_to_keep <- read_csv(taglist_path, show_col_types = FALSE) %>%
  select(TAG_ID) %>%
  distinct() %>%
  pull(TAG_ID) %>%
  as.character()

### Load stationary metadata
StationaryMetadata <- readxl::read_xlsx(stationary_metadata_dir,
                                        sheet="Metadata",
                                        skip=1) %>% 
  mutate(ArraySiteName_date = paste(str_replace_all(paste0(Array, SiteName), " ", ""),
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
    
    ### inserting site metadata
    mutate(ArraySiteName_date = paste(strsplit(basename(file_path), "_")[[1]][1:2], collapse ="_")) %>% 
    left_join(StationaryMetadata, by="ArraySiteName_date")  %>%
    
    ### selecting a subset of columns that are of interest
    select(SiteName, DateTime, TagCode, 
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
    select(any_of(c("DateTime", "TagCode", "Tilt", "SiteName"))) %>%
    mutate(
      TagCode = str_sub(TagCode, 4, 7),
      Tilt = as.numeric(Tilt),
      DateTime = parse_date_time(DateTime, orders = "mdY HMS", tz = "UTC")
    ) %>%
    filter(complete.cases(.))
  
  after_tag_filter_n <- df %>% filter(TagCode %in% tagcodes_to_keep) %>% nrow()
  
  message("🧾 Initial rows: ", format(initial_n, big.mark = ","),
          " (", format(n_indiv(df), big.mark = ","), " individuals)")
  
  df <- df %>% filter(TagCode %in% tagcodes_to_keep)
  
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
csv_files <- csv_files[1:3]

for (file_path in csv_files) {
  
  ### this is silly but it's kinda nice to know how far it's gotten
  cat("file", which(csv_files==file_path), "of", length(csv_files),'\n')
  
  message("📄 Processing: ", basename(file_path))
  
  df <- tryCatch(read_clean_vroom(file_path), error = function(e) {
    message("⚠️ Failed to read: ", basename(file_path), ": ", e$message)
    return(NULL)
  })
  if (is.null(df)) next
  
  # 💾 Save cleaned-but-unfiltered data (before tag filtering)
  unfiltered_vroom_out <- file.path(unfiltered_clean_dir, paste0("Unfiltered_", basename(file_path)))
  fwrite(df, unfiltered_vroom_out)
  message("💾 Saved unfiltered vroom-cleaned file: ", unfiltered_vroom_out)
  
  
  unfiltered <- tryCatch(clean_but_unfiltered_data(df), error = function(e) {
    message("⚠️ Unfiltered cleaning failed: ", e$message)
    return(NULL)
  })
  if (is.null(unfiltered)) next
  
  filtered <- tryCatch(apply_all_filters(unfiltered), error = function(e) {
    message("⚠️ Filtering failed: ", e$message)
    return(NULL)
  })
  if (is.null(filtered)) next
  
  filtered_out <- file.path(output_dir, paste0("Cleaned_", basename(file_path)))
  fwrite(filtered, filtered_out)
  message("✔ Saved: ", filtered_out)
}
