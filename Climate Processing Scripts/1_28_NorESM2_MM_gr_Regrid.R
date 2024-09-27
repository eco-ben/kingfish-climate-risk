### CMIP6 data regridding 
regrid.nc <- function(input.path, cores, var.name.loc, var.name.loc2,
                      time.var, grd.path) {
  library(doParallel)
  library(parallel)
  
  ipath <- input.path
  
  dir.nc <- paste(list.dirs(path = ipath, full.names = TRUE, recursive = FALSE)) #lists folders inside model folder (variables)
  opath <- paste(list.dirs(path = ipath, full.names = TRUE, recursive = FALSE)) 
  file.grd <- paste(grd.path, list.files(path = paste(grd.path, sep = "/"),
                                         pattern = "*.grd"), sep = "/") #grid file location
  
  for(i in 1:length(dir.nc)) { #for each folder inside model folder:
    dir.nc[i] <- paste0(dir.nc[i], '/RAW/', time.var) #lists folders within raw folder?
    opath[i] <- paste0(opath[i], '/Regridded/', time.var) #output into regridded folder
    files.nc <- paste(dir.nc[i], list.files(path = paste(dir.nc[i], sep = "/"), pattern = "*.nc"), 
                      sep = "/") #lists the files within each raw/time folder 
    
    cl <- makeCluster(cores)  
    registerDoParallel(cl)
    
    
    foreach(j = 1:length(files.nc)) %dopar% { #for each folder within the raw/time folder (e.g. RAW/ssp26)
      var.name <- sapply(strsplit(as.character(files.nc[j]), '/'),'[', var.name.loc)
      var.name <- sapply(strsplit(as.character(var.name), '_'),'[', var.name.loc2)
      system(paste(paste("cdo -P 10 -remapbil,", file.grd, sep = ""), #creates the cdo command (-P 10 = use 10 cores) (remapbil = bilinear interpolation)
                   files.nc[j], paste0(opath[i], '/', 'Regridded_', basename(files.nc[j])), sep = (" "))) #output file name
    } 
    stopCluster(cl)
  }
}

## CESM2-WACCM
## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp85',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# 
# ## CMCC-ESM2
# ## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
# regrid.nc(input.path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp85',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# 
# ## GFDL-ESM4.gr
# ## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
# regrid.nc(input.path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp85',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# 
# ## IPSL-CM6A-LR
# ## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
# regrid.nc(input.path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp85',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# 
# ## MPI-ESM1-2-HR
# ## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
# regrid.nc(input.path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp85',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# 
# ## NorESM2-LM.gr
# ## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp85',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')

# ## NorESM2-MM.gr
# ## 1) historical
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='historical',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 2) ssp26
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp26',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp45
# regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr',
#           cores=10, var.name.loc=10, var.name.loc2=1,
#           time.var='ssp45',
#           grd.path='/media/sun/Data1/CMIP6/Process_Codes')
# ## 4) ssp85
regrid.nc(input.path='/media/sun/Data1/CMIP6/NorESM2-MM.gr',
          cores=10, var.name.loc=10, var.name.loc2=1,
          time.var='ssp85',
          grd.path='/media/sun/Data1/CMIP6/Process_Codes')



