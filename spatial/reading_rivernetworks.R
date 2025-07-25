library(riverdist)
library(tidyverse)

## loading shapefiles and creating two versions of
## rivernetwork objects for use within riverdist.

## These have been stored in loadable .Rdata files.

## -- intentionally commented out
# kenai1_raw <- line2network(path="C:\\Users\\mbtyers\\Documents\\ArcGIS\\",
#                       layer="Kenai_albers")
# kenai1 <- cleanup(kenai1_raw)
# 
# kenai2_raw <- line2network(path="C:\\Users\\mbtyers\\Documents\\ArcGIS\\",
#                            layer="Kenai_all_dissolve")
# kenai2 <- cleanup(kenai2_raw)
# 
# plot(kenai1)
# plot(kenai2)
# save(kenai1, file="spatial/kenai1_rivernetwork.Rdata")
# save(kenai2, file="spatial/kenai2_rivernetwork.Rdata")
## -- intentionally commented out

# this one is a lot simpler - only named streams
load(file="spatial/kenai1_rivernetwork.Rdata")

# # this one is much more complex: don't use it unless we have to!
# load(file="spatial/kenai2_rivernetwork.Rdata")


plot(kenai1)

# stationary_metadata_dir <- "C:\\Users\\mbtyers\\Documents\\Current Projects\\Kenai_JuvChin_AcousticTelem\\data_processing\\TESTING\\StationaryMetadata.xlsx"
stationary_metadata_dir <- "data_processing\\TESTING\\StationaryMetadata.xlsx"
StationaryMetadata <- readxl::read_xlsx(stationary_metadata_dir,
                                        sheet="Metadata") %>% #skip=1
  mutate(ArraySiteName_date = paste(str_replace_all(paste0(Array, SiteID), " ", ""),
                                    paste0(25, substr(as.character(DeploymentDate),6,7), 
                                           substr(as.character(DeploymentDate),9,10)), sep="_"))  %>%
  select(ArraySiteName_date, `River Kilometer`, Latitude, Longitude)   ### can keep additional columns too


StationaryMetadata_albers <- sf::sf_project(pts=StationaryMetadata[,4:3], to="+proj=aea +lat_1=55 +lat_2=65 
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs 
    +ellps=GRS80")
points(StationaryMetadata_albers)

StationaryMetadata_segvert <- xy2segvert(x=StationaryMetadata_albers[,1],
                                         y=StationaryMetadata_albers[,2],
                                         rivers=kenai1)
StationaryMetadata <- cbind(StationaryMetadata,
                            StationaryMetadata_segvert) 
# StationaryMetadata$seg[substr(StationaryMetadata$ArraySiteName_date,1,6)=="Marine"] <- NA
# StationaryMetadata$vert[substr(StationaryMetadata$ArraySiteName_date,1,6)=="Marine"] <- NA
StationaryMetadata$rkm_riverdist <- mouthdist(seg=StationaryMetadata$seg,
                                              vert=StationaryMetadata$vert,
                                              rivers=kenai1)/1000

StationaryMetadata$rkm_riverdist[substr(StationaryMetadata$ArraySiteName_date,1,6)=="Marine"] <- -2
