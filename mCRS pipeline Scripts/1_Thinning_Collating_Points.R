setwd('/home/sun/Analyses/Kingfish/Global/')
library(sf)
library(ggplot2)
library(spThin)
library(dplyr)

#### South Pacific ----
SP_Points <- read.csv('Data/SouthPacific_Coordinates.csv')
thin(SP_Points, lat.col = 'Alt.Latitude', long.col = 'Alt.Longitude', spec.col = 'PopID', thin.par = 99, reps = 1, locs.thinned.list.return = FALSE, 
     write.files = T, max.files = 1, out.dir = 'Data/', out.base = 'SP',write.log.file = F, verbose = T)
SP_thin <- read.csv('Data/SP_thin1.csv')

SP_Points_sf <- st_as_sf(SP_Points, coords = c(1,2), crs = 4326)
SP_thin_sf <- st_as_sf(SP_thin, coords = c(2,3), crs = 4326)

ggplot() + geom_sf(data = SP_Points_sf) + geom_sf(data = SP_thin_sf, col = 'red', size = 0.5)

SP_Data <- left_join(SP_thin[c(2,3)], SP_Points, by = c('Alt.Longitude','Alt.Latitude'))
rm(SP_Points, SP_Points_sf, SP_thin, SP_thin_sf)


#### South Africa ----
SA_Points <- read.csv('Data/SouthAfrica_Coordinates.csv')
thin(SA_Points, lat.col = 'Alt.Latitude', long.col = 'Alt.Longitude', spec.col = 'PopID', thin.par = 99, reps = 1, locs.thinned.list.return = FALSE, 
     write.files = T, max.files = 1, out.dir = 'Data/', out.base = 'SA',write.log.file = F, verbose = T)
SA_thin <- read.csv('Data/SA_thin1.csv')

SA_Points_sf <- st_as_sf(SA_Points, coords = c(1,2), crs = 4326)
SA_thin_sf <- st_as_sf(SA_thin, coords = c(2,3), crs = 4326)

ggplot() + geom_sf(data = SA_Points_sf) + geom_sf(data = SA_thin_sf, col = 'red', size = 0.5)  + coord_sf(xlim = c(0,70), ylim = c(0,-50))

SA_Data <- left_join(SA_thin[c(2,3)], SA_Points, by = c('Alt.Longitude','Alt.Latitude'))
rm(SA_Points, SA_Points_sf, SA_thin, SA_thin_sf)


#### NE Pacific
NE_Points <- read.csv('Data/NEPacific_Coordinates.csv')
thin(NE_Points, lat.col = 'Alt.Latitude', long.col = 'Alt.Longitude', spec.col = 'PopID', thin.par = 99, reps = 1, locs.thinned.list.return = FALSE, 
     write.files = T, max.files = 1, out.dir = 'Data/', out.base = 'NE',write.log.file = F, verbose = T)
NE_thin <- read.csv('Data/NE_thin1.csv')

NE_Points_sf <- st_as_sf(NE_Points, coords = c(1,2), crs = 4326)
NE_thin_sf <- st_as_sf(NE_thin, coords = c(2,3), crs = 4326)

ggplot() + geom_sf(data = NE_Points_sf) + geom_sf(data = NE_thin_sf, col = 'red', size = 0.5)  + coord_sf(xlim = c(-140,-100), ylim = c(0,50))

NE_Data <- left_join(NE_thin[c(2,3)], NE_Points, by = c('Alt.Longitude','Alt.Latitude'))
rm(NE_Points, NE_Points_sf, NE_thin, NE_thin_sf)


#### NW Pacific ----
NW_Points <- read.csv('Data/NWPacific_Coordinates.csv')
thin(NW_Points, lat.col = 'Alt.Latitude', long.col = 'Alt.Longitude', spec.col = 'PopID', thin.par = 99, reps = 1, locs.thinned.list.return = FALSE, 
     write.files = T, max.files = 1, out.dir = 'Data/', out.base = 'NW',write.log.file = F, verbose = T)
NW_thin <- read.csv('Data/NW_thin1.csv')

NW_Points_sf <- st_as_sf(NW_Points, coords = c(1,2), crs = 4326)
NW_thin_sf <- st_as_sf(NW_thin, coords = c(2,3), crs = 4326)

ggplot() + geom_sf(data = NW_Points_sf) + geom_sf(data = NW_thin_sf, col = 'red', size = 0.5)  + coord_sf(xlim = c(120,150), ylim = c(20,40))

NW_Data <- left_join(NW_thin[c(2,3)], NW_Points, by = c('Alt.Longitude','Alt.Latitude'))
rm(NW_Points, NW_Points_sf, NW_thin, NW_thin_sf)


### Combining ----
Total <- rbind(NE_Data, NW_Data, SA_Data, SP_Data)
Total$Extract.ID <- seq(1,217,1)
write.csv(Total, 'Data/Global_Coordinates.csv', row.names = F)
