netcdf.to.raster <- function(input.path, time.var, var.name.loc, var.name.loc2) {
  library(ncdf4)
  library(raster)
  
  ipath <- input.path
  
  dir.nc <- paste(list.dirs(path = ipath, full.names = TRUE, recursive = FALSE)) #list directories within model folder
  opath <- paste(list.dirs(path = ipath, full.names = TRUE, recursive = FALSE)) 
  
  lyr <- list()
  
  for(i in 1:length(dir.nc)) {
    dir.nc[i] <- paste0(dir.nc[i], '/Merged/', time.var) #replaces values in dir.nc with the folder for the right time var
    opath[i] <- paste0(opath[i], '/GeoTiff/', time.var) #replaces values in output path witht the folder for the right time var
    files.nc <- paste(dir.nc[i], list.files(path = paste(dir.nc[i], sep = "/"), pattern = "*.nc"), 
                      sep = "/")
  # files.nc <- paste(list.files(path = input.path, full.names = TRUE, recursive = FALSE))
    
  # lyr <- list()
  
  for (j in 1:length(files.nc)) {
    # var.name <- sapply(strsplit(as.character(files.nc[[j]]), '/'),'[', var.name.loc)
    # var.name <- sapply(strsplit(as.character(var.name), '[.]'),'[', 1)
    var.name <- sapply(strsplit(as.character(files.nc[[j]]), '/'),'[', var.name.loc)
    var.name <- sapply(strsplit(as.character(var.name), '_'),'[', var.name.loc2)
    
    nc.info <- nc_open(files.nc[[j]])
    depth.levels <- nc.info$dim[[5]]
    
    if (depth.levels$units %in% c('m','meters','ms')) {
      depth.levels$vals <- depth.levels$vals
    } else if (depth.levels$units %in% c('cm', 'centimeters', 'cms')) {
      depth.levels$vals <- depth.levels$vals / 100
      depth.levels$units <- 'm'
    }
    
    depth.ids <- depth.levels$vals[depth.levels$vals <= 150]
    
    for (k in 1:length(depth.ids)) {
      lyr[[j]] <- brick(files.nc[[j]],  level=k)
      writeRaster(x=lyr[[j]],
                  # filename=paste0(opath[i],'/',time.var,'_',k*10,'m_',var.name,'.tif'),
                  filename=paste0(opath[i],'/',time.var,'_',depth.ids[[k]],'m_',var.name,'.tif'),
                  format='GTiff', overwrite=T)
      }
    }
  }
}

# CESM2-WACCM
netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', time.var = 'historical',
          var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)

# CMCC-ESM2
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', time.var = 'historical',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)

# GFDL-ESM4.gr
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', time.var = 'historical',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)

# IPSL-CM6A-LR
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', time.var = 'historical',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)

# MPI-ESM1-2-HR
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', time.var = 'historical',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)

# NorESM2-LM.gr
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', time.var = 'historical',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)

# NorESM2-MM.gr
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', time.var = 'historical',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', time.var = 'ssp26',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', time.var = 'ssp45',
#                  var.name.loc=10, var.name.loc2 = 2)
# netcdf.to.raster(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', time.var = 'ssp85',
#                  var.name.loc=10, var.name.loc2 = 2)
