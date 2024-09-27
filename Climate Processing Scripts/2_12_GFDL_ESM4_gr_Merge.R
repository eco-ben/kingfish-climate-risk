### CMIP6 data regridding 
merge.nc <- function(input.path, cores, time.loc, var.name.loc, var.name.loc2, model.loc, time.var) {
  library(doParallel)
  library(parallel)
  
  ipath <- input.path
  
  dir.nc <- paste(list.dirs(path = ipath, full.names = TRUE, recursive = FALSE)) #list directories within model folder
  opath <- paste(list.dirs(path = ipath, full.names = TRUE, recursive = FALSE)) 
  
  for(i in 1:length(dir.nc)) {
    dir.nc[i] <- paste0(dir.nc[i], '/Regridded/', time.var) #replaces values in dir.nc with the folder for the right time var
    opath[i] <- paste0(opath[i], '/Merged/', time.var) #replaces values in output path witht the folder for the right time var
    files.nc <- paste(dir.nc[i], list.files(path = paste(dir.nc[i], sep = "/"), pattern = "*.nc"), 
                      sep = "/")
    time.period <- sapply(strsplit(as.character(files.nc[1]), '/'),'[', time.loc)
    model.name <- sapply(strsplit(as.character(files.nc[1]), '/'),'[', model.loc)
    
    cl <- makeCluster(cores)  
    registerDoParallel(cl)
    
    
    parameters <- c()
    
    for (j in 1:length(files.nc)) {
      var.name <- sapply(strsplit(as.character(files.nc[j]), '/'),'[', var.name.loc)
      parameter.name <- sapply(strsplit(as.character(var.name), '_'),'[', var.name.loc2) 
      parameters <- c(parameters, parameter.name)
    }
    
    parameters <- unique(parameters)
    
    parameter.list <- vector('list', length=length(parameters))
    
    for (j in 1:length(parameters)) {
      dummy.list <- subset(files.nc, grepl(parameters[[j]], files.nc))
      parameter.list[[j]] <- dummy.list
    }
    
    for (j in 1:length(parameter.list)) {
      parameter.name <- parameters[[j]]
      system(paste(paste("cdo -P 10 mergetime", paste(parameter.list[[j]], collapse=' ')),
                   paste0(opath[i], '/', 'Merged_',parameter.name,'_',model.name,'_',time.period,'.nc')))
    }
    
    # system(paste(paste("cdo -P 10 mergetime ", paste(files.nc, collapse=' ')),
    #                paste0(opath, 'Merged_',parameter.name, '_', model.name, '.nc')))
    
    
    # foreach(j = 1:length(files.nc)) %dopar% {
    #   var.name <- sapply(strsplit(as.character(files.nc[j]), '/'),'[', var.name.loc)
    #   parameter.name <- sapply(strsplit(as.character(var.name), '_'),'[', var.name.loc2)
    #   model.name <- sapply(strsplit(as.character(var.name), '_'),'[', model.loc)
    #   system(paste(paste("cdo -P 10 ", parameter.name, "*.nc "),
    #                paste0(opath, 'Merged_',parameter.name, '_', model.name, '.nc')))
    # } 
    stopCluster(cl)
  }
}

# CESM2-WACCM
# merge.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', cores=10, time.loc=9, var.name.loc=10,
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
# merge.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp85')

#CMCC-ESM2
# merge.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', cores=10, time.loc=9, var.name.loc=10,
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
# merge.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp85')

#GFDL-ESM4.gr
# merge.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', cores=10, time.loc=9, var.name.loc=10,
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', cores=10, time.loc=9, var.name.loc=10,
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
merge.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr', cores=10, time.loc=9, var.name.loc=10,
         var.name.loc2=2, model.loc=6, time.var='ssp85')

#IPSL-CM6A-LR
# merge.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', cores=10, time.loc=9, var.name.loc=10, 
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
# merge.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp85')

#MPI-ESM1-2-HR
# merge.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', cores=10, time.loc=9, var.name.loc=10, 
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
# merge.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp85')

#NorESM2-LM.gr
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', cores=10, time.loc=9, var.name.loc=10, 
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp85')

#NorESM2-MM.gr
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', cores=10, time.loc=9, var.name.loc=10, 
#              var.name.loc2=2, model.loc=6, time.var='historical')
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp26')
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp45')
# merge.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr', cores=10, time.loc=9, var.name.loc=10, 
#          var.name.loc2=2, model.loc=6, time.var='ssp85')






