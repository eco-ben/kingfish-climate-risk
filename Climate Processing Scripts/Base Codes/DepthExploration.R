setwd('/media/sun/Data1/CMIP6/')
library(ncdf4)
library(ggplot2)
depths <- data.frame()



CESM2.info <- nc_open('CESM2-WACCM/no3/Merged/ssp26/Merged_no3_CESM2-WACCM_ssp26.nc')
CMCC.info <- nc_open('CMCC-ESM2/no3/Merged/ssp26/Merged_no3_CMCC-ESM2_ssp26.nc')
GFDL.info <- nc_open('GFDL-ESM4.gr/no3/Merged/ssp26/Merged_no3_GFDL-ESM4.gr_ssp26.nc')
IPSL.info <- nc_open('IPSL-CM6A-LR/no3/Merged/ssp26/Merged_no3_IPSL-CM6A-LR_ssp26.nc')
MPI.info <- nc_open('MPI-ESM1-2-HR/no3/Merged/ssp26/Merged_no3_MPI-ESM1-2-HR_ssp26.nc')
NorLM.info <- nc_open('NorESM2-LM.gr/no3/Merged/ssp26/Merged_no3_NorESM2-LM.gr_ssp26.nc')
NorMM.info <- nc_open('NorESM2-MM.gr/no3/Merged/ssp26/Merged_no3_NorESM2-MM.gr_ssp26.nc') 
infos <- list('CESM2.info','CMCC.info','GFDL.info','IPSL.info','MPI.info','NorLM.info','NorMM.info')

for(i in 1:length(infos)){
  info <- get(infos[[i]])
  levels <- info$dim[[5]]
  if (levels$units %in% c('m','meters','ms')) {
    levels$vals <- levels$vals
  } else if (levels$units %in% c('cm', 'centimeters', 'cms')) {
    levels$vals <- levels$vals / 100
    levels$units <- 'm'
  }
  ids <- levels$vals[levels$vals <= 150]
  model <- stringr::str_extract(infos[[i]],'^(?<=^).*?(?=\\.info$)')
  depth <- data.frame(Model <- character(length(ids)))
  depth$Levels <- ids
  for(i in nrow(depth)){
    depth$Model <- model
  }
  depth <- depth[,-1]
  depths <- rbind(depths,depth)
}

ggplot() +
  geom_point(data = depths, aes(x = Model, y = Levels)) + scale_y_continuous(breaks = seq(0,150,by = 10))


