# ── Load libraries ──
library(tidyverse)
library(riverdist)

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

filtered_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\Cleaned Data\\all_filtered.csv"
# unfiltered_clean_dir <- file.path(output_dir, "Unfiltered Clean Files")
taglist_path <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\chinook_tagged.csv"

# stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\StationaryMetadata.csv"
stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\StationaryMetadata.xlsx"




# # ── Load tagged fish IDs ──
# tagcodes_to_keep <- read_csv(taglist_path, show_col_types = FALSE) %>%
#   select(TAG_ID) %>%
#   distinct() %>%
#   pull(TAG_ID) %>%
#   as.character()

### loading filtered receiver data
all_filtered_df <- read_csv(filtered_dir)


### date-time associated with release of each tag
tagging_data <- read_csv(taglist_path, show_col_types = FALSE) %>%
  # mutate(release_datetime = as.POSIXlt(paste(PULL_DATE, END)), format="%Y-%m-%d %H:%M:%OS") %>%
  # mutate(release_datetime = parse_date_time(paste(PULL_DT, END), orders = "Ymd HMS", tz = "UTC")) %>%
  mutate(release_datetime = END_DT) %>%
  filter(FATE == "R") %>%
  select(TAG_ID, release_datetime, LAT, LONG) %>%
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








### summarizing passage date for each individual

# orig metadata
metadata_rivermile <- readxl::read_xlsx(stationary_metadata_dir,
                                        sheet="Metadata") %>%
  mutate(Array = str_replace_all(Array, "[ ']", "")) %>%
  mutate(Array = ifelse(Array %in% c("MarineNorth", "MarineSouth"), "Marine", Array)) %>% #### take this out as needed
  select(Array, `River Mile`, `River Kilometer`,Latitude, Longitude) %>%
  group_by(Array) %>%
  summarise(rivermile = mean(`River Mile`, na.rm=TRUE),
            river_km = mean(`River Kilometer`, na.rm=TRUE),
            Latitude = mean(Latitude, na.rm=TRUE),
            Longitude = mean(Longitude, na.rm=TRUE)) %>%
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



## establishing river km location of tagging

# this one is a lot simpler - only named streams
load(file="spatial/kenai1_rivernetwork.Rdata")

## I do this calculation an annoying amount of times!
latlong_to_rkm <- function(latlong) {
  pts_albers <- sf::sf_project(pts=latlong[2:1], to="+proj=aea +lat_1=55 +lat_2=65 
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
    +ellps=GRS80")
  pts_segvert <- xy2segvert(x=pts_albers[,1],
                            y=pts_albers[,2],
                            rivers=kenai1)
  mouthdist(seg=pts_segvert$seg,
            vert=pts_segvert$vert,
            rivers=kenai1)/1000
} 

# StationaryMetadata_albers <- sf::sf_project(pts=StationaryMetadata[,4:3], to="+proj=aea +lat_1=55 +lat_2=65 
#     +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
#     +ellps=GRS80")
# 
# StationaryMetadata_segvert <- xy2segvert(x=StationaryMetadata_albers[,1],
#                                          y=StationaryMetadata_albers[,2],
#                                          rivers=kenai1)
# StationaryMetadata <- cbind(StationaryMetadata,
#                             StationaryMetadata_segvert) 
# StationaryMetadata$seg[substr(StationaryMetadata$ArraySiteName_date,1,6)=="Marine"] <- NA
# StationaryMetadata$vert[substr(StationaryMetadata$ArraySiteName_date,1,6)=="Marine"] <- NA
# StationaryMetadata$rkm_riverdist <- mouthdist(seg=StationaryMetadata_segvert$seg,
#                                               vert=StationaryMetadata_segvert$vert,
#                                               rivers=kenai1)/1000
StationaryMetadata$rkm_riverdist <- latlong_to_rkm(latlong=StationaryMetadata[,3:4])


StationaryMetadata$rkm_riverdist[substr(StationaryMetadata$ArraySiteName_date,1,6)=="Marine"] <- -2


# tagging_albers <- sf::sf_project(pts=tagging_data[,4:3], to="+proj=aea +lat_1=55 +lat_2=65 
#     +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
#     +ellps=GRS80")
# tagging_segvert <- xy2segvert(x=tagging_albers[,1],
#                                          y=tagging_albers[,2],
#                                          rivers=kenai1)
# tagging_data$rkm_riverdist <- mouthdist(seg=tagging_segvert$seg,
#                                               vert=tagging_segvert$vert,
#                                               rivers=kenai1)/1000
tagging_data$rkm_riverdist <- latlong_to_rkm(tagging_data[,3:4])

all_filtered_passage$rkm_riverdist <- latlong_to_rkm(all_filtered_passage[6:7])
all_filtered_passage$rkm_riverdist[all_filtered_passage$Array=="Marine"] <- -2




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

plot2 <- rbind(data.frame(TagCode = tagging_data$TagCode,
                          Array = "tagging",
                          passage = tagging_data$release_datetime,
                          rivermile = NA,
                          river_km = tagging_data$rkm_riverdist,
                          Array_order = NA,
                          Array_order_rev = "0 tagging"), all_filtered_passage) %>%
  ggplot(aes(x=passage, y=river_km, colour=TagCode, group=TagCode)) +
  geom_point() +
  geom_line() +
  theme_bw()  +
  xlab("Passage Date (max sig strength)") +
  ylab("River km")

plot1 / plot2 + plot_layout(guides="collect")



## making detection matrix
detect_mat_tibble <- pivot_wider(all_filtered_passage,
                               id_cols=TagCode,
                               names_from = Array_order_rev,
                               values_from = passage)
detect_mat <- detect_mat_tibble[,order(colnames(detect_mat_tibble))] %>%
  (\(x) x[, colnames(x) != "TagCode"]) %>%
  as.matrix %>%
  is.na %>% `!`
rownames(detect_mat) <- detect_mat_tibble$TagCode



## establishing array of entry
# compare tagging_data to all_filtered_passage

# vector of rkm_riverdist for each array
rkm_arrays <- tapply(all_filtered_passage$rkm_riverdist,
       all_filtered_passage$Array_order_rev,
       mean)

# df of rkm_riverdist from tagging_data, corresponding to detect_mat
detect_mat_tagging <- data.frame(TagCode=rownames(detect_mat)) %>%
  left_join(tagging_data)

detect_mat_tagging$array_entry_vec <- NA
for(i in 1:nrow(detect_mat_tagging)) {
  diffs <- detect_mat_tagging$rkm_riverdist[i] - rkm_arrays
  diffs[diffs < 0] <- 1000
  detect_mat_tagging$array_entry_vec[i] <- which.min(diffs)
}

# making sure that detection only happens after tagging!
par(mfrow=c(2,3))
for(check in sort(unique(detect_mat_tagging$array_entry_vec))) {
  plot(colMeans(detect_mat[detect_mat_tagging$array_entry_vec==check,,drop=F]),
       main=check)
}
