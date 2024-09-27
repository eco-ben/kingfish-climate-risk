setwd('/home/sun/Analyses/Kingfish/Global/')
library(data.table)
library(ggplot2)


RData <- list.files(path = 'RData/CRS_results/')
RData <- RData[grep('ensemble',RData)]
CRS.dt <- data.table()
PremeanDT <- data.table()
stacks <- list('CRS.ssp26.100m.ensemble','CRS.ssp26.10m.ensemble','CRS.ssp26.150m.ensemble','CRS.ssp26.50m.ensemble','CRS.ssp45.100m.ensemble','CRS.ssp45.10m.ensemble',
               'CRS.ssp45.150m.ensemble','CRS.ssp45.50m.ensemble','CRS.ssp85.100m.ensemble','CRS.ssp85.10m.ensemble',
               'CRS.ssp85.150m.ensemble','CRS.ssp85.50m.ensemble') 

for(i in 1:length(RData)){
  RData[[i]] <- paste0('RData/CRS_results/',RData[[i]])
}

new_dts <- list()
PrePop <- list()

for(i in 1:length(RData)){
  load(RData[[i]])
  dt.name <- stacks[[i]]
  dt <- get(stacks[[i]])
  # dt$Depth <- 0
  # dt$Depth <- stringr::str_extract(dt.name,"(?<=\\.)\\d+(?=m)")
  CRS.refmean <- dt[dt$Scenario == 'historical',.(CRS.refmean = mean(CRS)),by =.(PopID,Extract.ID)]
  dt$Date <- as.Date(dt$Date)
  dt$Year <- as.numeric(format(dt$Date,'%Y'))
  dt <- dt[CRS.refmean, CRS.refmean := CRS.refmean, on = 'Extract.ID']
  dt <- dt[, .(DeltaCRS = (CRS - CRS.refmean)), by = .(Extract.ID, PopID,Date, Year, Scenario, CRS)]
  Popmeans <- dt[,.(DeltaCRS.Popmean = mean(DeltaCRS),CRS.Popmean = mean(CRS)),by=.(PopID,Extract.ID,Year, Scenario)]
  Popmeans <- Popmeans[,.(DeltaCRS.Popmean = mean(DeltaCRS.Popmean),CRS.Popmean = mean(CRS.Popmean),Sample.size = nrow(Popmeans)),by=.(PopID,Year, Scenario)]
  Popmeans$Depth <- stringr::str_extract(dt.name,"(?<=\\.)\\d+(?=m)")
  dt$Depth <- stringr::str_extract(dt.name,"(?<=\\.)\\d+(?=m)")
  #Popmeans$Scenario <- stringr::str_extract(dt.name,"ssp\\d+")
  Popmeans$PopID_Depth <- factor(paste(Popmeans$PopID, Popmeans$Depth, sep = "-"))
  new_dts[[i]] <- Popmeans
  PrePop[[i]] <- dt
  PremeanDT <- rbindlist(list(PremeanDT, PrePop[[i]])) # Store results for each ID and each Date
  PrePop[[i]] <- NULL
  CRS.dt <- rbindlist(list(CRS.dt, new_dts[[i]])) # Store results for the population means
  new_dts[[i]] <- NULL
  rm(dt,Popmeans,CRS.refmean)
  rm(list = c(dt.name))
}

write.csv(PremeanDT, file = 'Results/Ensemble/DeltaCRS/DeltaCRS_mCRS_RAW.csv', row.names = F)
write.csv(CRS.dt, file = 'Results/Ensemble/DeltaCRS/DeltaCRS_PopValues.csv', row.names = F)




