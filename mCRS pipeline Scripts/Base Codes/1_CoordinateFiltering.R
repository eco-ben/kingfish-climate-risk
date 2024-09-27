########################################## RUN IF Occurrence.csv Gets Updated! ###########################

setwd("/home/sun/Analyses/Kingfish/Global/")
library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(sf)
library(ggplot2)
library(lwgeom)
library(tidyr)

## 1) Filter Only Australia --------------

## OBIS Data 
OBISData <- read.csv("Data/Occurrence.csv")
OBISData <- OBISData[,c('decimallongitude','decimallatitude','date_mid')]
OBISData$longitude_360 <- ifelse(OBISData$decimallongitude < 0, 
                                 OBISData$decimallongitude + 360, OBISData$decimallongitude)
OBISData <- OBISData[!duplicated(OBISData[,c(1,2)]),]

OBISData$PopID <- 0
OBISData$PopID <- ifelse(OBISData$longitude_360 >= 105 & OBISData$longitude_360 <= 295 & OBISData$decimallatitude < -10, 'South_Pacific', OBISData$PopID)
OBISData$PopID <- ifelse(OBISData$longitude_360 <= 265 & OBISData$longitude_360 >= 220 & OBISData$decimallatitude > 0, 'NE_Pacific', OBISData$PopID)
OBISData$PopID <- ifelse(OBISData$longitude_360 >= 120 & OBISData$longitude_360 <= 150 & OBISData$decimallatitude > 20, 'NW_Pacific', OBISData$PopID)
OBISData$PopID <- ifelse(OBISData$longitude_360 >= 0 & OBISData$longitude_360 <= 45 & OBISData$decimallatitude <= -20, 'South_Africa', OBISData$PopID)
OBISData <- OBISData[!OBISData$PopID == 0,]

colnames(OBISData)[c(1,2)] <- c('Longitude','Latitude')
OBISData <- OBISData[!duplicated(OBISData[,c("Longitude","Latitude")]),]

GBIF <- read_delim('Data/GBIF_Occurrence.csv', delim = "\t", escape_double = FALSE, trim_ws = TRUE)
GBIF <- GBIF[,c('decimalLongitude','decimalLatitude','eventDate')]
GBIF$longitude_360 <- ifelse(GBIF$decimalLongitude < 0, 
                             GBIF$decimalLongitude + 360, GBIF$decimalLongitude)

## GBIF Data
GBIF <- GBIF[!duplicated(GBIF[,c(1,2)]),]
GBIF <- GBIF[complete.cases(GBIF[,c(1,2)]),]

GBIF$PopID <- 0
GBIF$PopID <- ifelse(GBIF$longitude_360 >= 105 & GBIF$longitude_360 <= 295 & GBIF$decimalLatitude < -10, 'South_Pacific', GBIF$PopID)
GBIF$PopID <- ifelse(GBIF$longitude_360 <= 265 & GBIF$longitude_360 >= 220 & GBIF$decimalLatitude > 0, 'NE_Pacific', GBIF$PopID)
GBIF$PopID <- ifelse(GBIF$longitude_360 >= 120 & GBIF$longitude_360 <= 150 & GBIF$decimalLatitude > 20, 'NW_Pacific', GBIF$PopID)
GBIF$PopID <- ifelse(GBIF$longitude_360 >= 0 & GBIF$longitude_360 <= 45 & GBIF$decimalLatitude <= -20, 'South_Africa', GBIF$PopID)
GBIF <- GBIF[!GBIF$PopID == 0,]
colnames(GBIF)[c(1,2)] <- c('Longitude','Latitude')
GBIF <- GBIF[!duplicated(GBIF[,c("Longitude","Latitude")]),]

## Combining Data Sources
OBISData$Date.Mid <- as.POSIXct(OBISData$date_mid/1000, origin = '1970-01-01')
OBISData$eventDate <- OBISData$Date.Mid
OBISData <- OBISData[,c(1,2,7,4,5)]
CombData <- rbind(GBIF, OBISData)
CombData <- CombData[!duplicated(CombData[,c(1,2)]),]
CombData$Month.Observed <- months(CombData$eventDate)

for(i in 1:nrow(CombData)){
  if(CombData[i,]$Latitude < 0){
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('January','February','December'), 'Summer',NA)
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('March','April','May'), 'Autumn', CombData$Season)
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('June','July','August'),'Winter', CombData$Season)
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('September','October','November'),'Spring', CombData$Season)
  }else{
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('January','February','December'), 'Winter',NA)
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('March','April','May'), 'Spring', CombData$Season)
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('June','July','August'),'Summer', CombData$Season)
    CombData$Season <- ifelse(CombData$Month.Observed %in% c('September','October','November'),'Autumn', CombData$Season)
  }
}


## 2) Remove coordinates on land -------------------
world <- rgdal::readOGR('Data/GIS/ne_10m_land.shp') # Natural Earth 10m resolution land shapefile

CombData.pts <- CombData[,c(1,2)]
colnames(CombData.pts) <- c('Longitude','Latitude')

pts <- sp::SpatialPoints(CombData.pts, proj4string=sp::CRS(sp::proj4string(world))) # CRS calibration for the generated grids
land.coords <- sp::over(pts, world) # Find coordinates over land
CombData.pts <- droplevels(CombData.pts[which(is.na(land.coords$featurecla)),]) # Restrict data to marine coordinates
row.names(CombData.pts) <- NULL

pts.id <- CombData.pts
pts.id$ID <- row.names(CombData.pts)

CombData <- sqldf::sqldf("SELECT * from 'CombData' LEFT JOIN
                                'pts.id' USING (Longitude, Latitude)")
CombData <- CombData[!is.na(CombData$ID),]
rm(pts.id,pts,world,land.coords)

CombDataSP <- CombData[CombData$PopID == 'South_Pacific',]
rm(CombData, OBISData, GBIF)
## 3) Add labels for occurrences as Population ID, PopID ID, Date Observed, Month and Season --------------

#change timestamp date into date format

#add season to the df


## 4) Moving coordinates that are too close to land -----

#### 4.1) Create Ocean Polygon -----
# climate_data <- brick('/media/sun/Data1/CMIP6/CESM2-WACCM/thetao/GeoTiff/historical/historical_5m_thetao.tif')
climate_data <- brick('/media/sun/Data1/CMIP6/CESM2-WACCM/thetao/GeoTiff/historical/historical_5m_thetao.tif')
climate_data <- climate_data[[1]] #only use the first time layer of the climate data
data_present <- !is.na(climate_data) #should identify squares with data in them
Data_polygon <- rasterToPolygons(data_present,dissolve = T)

Data_polygon_ocean <- st_as_sf(Data_polygon, crs=4326)
Data_polygon_ocean <- Data_polygon_ocean[Data_polygon_ocean$layer==1,]

#### 4.2) Test which points are too close to land ----
# load('/home/sun/Analyses/Kingfish/RData/OBIS/historical_summary_CESM2_WACCM.RData')
# OBISDataANZ$ID <- as.numeric(OBISDataANZ$ID)
# AllIDs <- levels(as.factor(OBISDataANZ$ID))
# ClimatePosIDs <- levels(as.factor(Historical.CESM2_WACCM.stack[complete.cases(Historical.CESM2_WACCM.stack),]$ID))
# rm(Historical.CESM2_WACCM.stack)
# noDataIDs <- setdiff(AllIDs,ClimatePosIDs)
# OBISDataANZ$Extracted <- ifelse(OBISDataANZ$ID %in% ClimatePosIDs, 'Yes', 'No')

#load coordinates and data sf
coordinates <- CombDataSP[,c(8,1,2)]
coordinates_sf <- st_as_sf(coordinates,coords = c(2,3),crs =4326) %>%
  dplyr::mutate(ID = coordinates[,1])

#identify points outside ocean
land.IDs <- sapply(st_intersects(coordinates_sf,Data_polygon_ocean),function(x){length(x)==1})

#### 4.3) Move points too close to land ------------
land.points <- coordinates_sf[which(land.IDs==T),]
edge.points <- st_nearest_points(land.points, Data_polygon_ocean) 
edge.linestring <- st_cast(edge.points, 'MULTILINESTRING')

extend_linestring <- function(linestring, additional_length) {
  # Extract the coordinates of the linestring
  coords <- st_coordinates(linestring)
  
  # Calculate the total length of the linestring
  total_length <- sum(sqrt(diff(coords[, "X"])^2 + diff(coords[, "Y"])^2))
  
  # Calculate the extension factor based on the desired additional length
  extension_factor <- 1 + (additional_length / total_length)
  
  # Interpolate new points along the linestring
  new_coords <- lapply(1:(nrow(coords)-1), function(i) {
    p1 <- coords[i, c("X", "Y")]
    p2 <- coords[i+1, c("X", "Y")]
    diff_vec <- p2 - p1
    new_p2 <- p1 + extension_factor * diff_vec
    rbind(p1, new_p2)
  })
  
  new_coords <- new_coords[seq(1, length(new_coords), 2)]
  
  return(new_coords)
}

# Specify the desired increase as a fixed additional length
increase_length <- 2500  # Replace with desired value

sfc_geom_extended <- extend_linestring(linestring=edge.linestring, additional_length=increase_length)
sfc_geom_extended <- lapply(sfc_geom_extended, st_linestring)
sfc_geom_extended <- lapply(sfc_geom_extended, st_sfc, crs=4326)

###
extended.linestring <- lapply(sfc_geom_extended, st_endpoint)
extended.coords <- as.data.frame(do.call(rbind, extended.linestring))
extended.coords$Longitude <- sapply(strsplit(as.character(extended.coords$V1), ','),'[', 1)
extended.coords$Longitude <- as.numeric(substr(extended.coords$Longitude, 3, nchar(extended.coords$Longitude)))
extended.coords$Latitude <- sapply(strsplit(as.character(extended.coords$V1), ','),'[', 2)
extended.coords$Latitude <- as.numeric(substr(extended.coords$Latitude, 1, (nchar(extended.coords$Latitude)-1)))
extended.coords <- extended.coords[,-1L]
extended.coords$ID <- land.points$ID
colnames(extended.coords) <- c('Alt.Longitude','Alt.Latitude','ID')

#  test.extract <- raster::extract(climate_data, extended.coords[,c(1,2)])
#  extended.coords$Value <- test.extract
# #

All.points <- dplyr::left_join(coordinates,extended.coords, by = c('ID'))
All.points$Alt.Longitude <- ifelse(is.na(All.points$Alt.Longitude),All.points$Longitude,All.points$Alt.Longitude)
All.points$Alt.Latitude <- ifelse(is.na(All.points$Alt.Latitude), All.points$Latitude, All.points$Alt.Latitude)

# new.points <- st_as_sf(All.points, coords = c(4,5),crs = 4326)
# old.points <- st_as_sf(All.points, coords = c(2,3), crs = 4326)
# 
# ggplot() +
#   geom_sf(data = Data_polygon_ocean)+
#   geom_sf(data = new.points, color = 'red') #+
#   # coord_sf(xlim = c(120,150), ylim = c(0,50))

CombDataSP <- dplyr::left_join(CombDataSP,All.points[,c('ID','Alt.Longitude','Alt.Latitude')], by = c('ID'))

#csv for filtered coordinates with likely useful information
write.csv(CombDataSP[,c(9,10,4,5,3,6,7,8)], file = "Data/SouthPacific_Coordinates.csv", row.names = F)
