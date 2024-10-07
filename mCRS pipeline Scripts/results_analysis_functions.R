preprocess_results_df = function(x, coordinates, level){
  # Calculate annual average of CRS and DeltaCRS
  x_dt = data.table(x)
  names(x_dt)[names(x_dt) == "Extract_ID"] <- "Extract.ID"
  x_dt <- x_dt[,.(DeltaCRS = mean(DeltaCRS),CRS = mean(CRS)),by=.(PopID,Depth,Extract.ID,Year, Scenario)]
  x_dt <- as.data.frame(x_dt)
  
  # Attach coordinates and analysis levels
  x_dt <- dplyr::left_join(x_dt, coordinates[,c('Extract.ID','Alt.Longitude',"Alt.Latitude","longitude_360","PopID")], by = c('Extract.ID','PopID'))
  x_dt$Distance_from_equator <- ifelse(x_dt$Alt.Latitude < 0, -x_dt$Alt.Latitude, x_dt$Alt.Latitude)
  x_dt$Depth <- as.factor(x_dt$Depth)
  x_dt$PopID <- as.factor(x_dt$PopID)
  x_dt$Level = level
  
  # Add Population label instead of full ID
  x_dt$PopLab <- as.character(x_dt$PopID)
  x_dt$PopLab <- ifelse(x_dt$PopLab == 'NE_Pacific', 'NEP', x_dt$PopLab)
  x_dt$PopLab <- ifelse(x_dt$PopLab == 'NW_Pacific', 'NWP', x_dt$PopLab)
  x_dt$PopLab <- ifelse(x_dt$PopLab == 'South_Africa', 'SA', x_dt$PopLab)
  x_dt$PopLab <- ifelse(x_dt$PopLab == 'South_Pacific', 'SIP', x_dt$PopLab)
  
  return(x_dt)
}

PopID_add_PopLab = function(x){
  x$PopLab = as.character(x$PopLab)
  x$PopLab <- ifelse(x$PopLab == 'NE_Pacific', 'NEP', x$PopLab)
  x$PopLab <- ifelse(x$PopLab == 'NW_Pacific', 'NWP', x$PopLab)
  x$PopLab <- ifelse(x$PopLab == 'South_Africa', 'SA', x$PopLab)
  x$PopLab <- ifelse(x$PopLab == 'South_Pacific', 'SIP', x$PopLab)
  
  return(x)
}
  
mean_Delta_temperature = function(x){
  # Calculate the mean historical temperature each population and depth is exposed to. 
  refmeanTemp <- x[x$Scenario == 'historical',] %>% group_by(PopID,Depth) %>% summarise(TempMean = mean(thetao.mean))
  temp_and_ref <- left_join(x, refmeanTemp, by = c('PopID', 'Depth'))
  
  # Calculate the change in temperature from historical mean
  temp_and_ref <- temp_and_ref %>% group_by(PopID,Depth,Year,Scenario) %>% reframe(DeltaTemp = (thetao.mean - TempMean))
  
  # Attach historical data to each scenario and relabel
  ssp26 <- temp_and_ref[temp_and_ref$Scenario %in% c('historical','ssp26'),]
  ssp45 <- temp_and_ref[temp_and_ref$Scenario %in% c('historical','ssp45'),]
  ssp85 <- temp_and_ref[temp_and_ref$Scenario %in% c('historical','ssp85'),]
  ssp26[ssp26$Scenario == 'historical',]$Scenario <- 'ssp26'
  ssp45[ssp45$Scenario == 'historical',]$Scenario <- 'ssp45'
  ssp85[ssp85$Scenario == 'historical',]$Scenario <- 'ssp85'
  
  # Combine all scenarios' data 
  all_scens <- rbind(ssp26,ssp45,ssp85)
  all_scens <- all_scens[!all_scens$Year <= 2000,]
  all_scens$PopID <- as.factor(all_scens$PopID)
  all_scens$Depth <- as.factor(all_scens$Depth)
  all_scens$Scenario <- as.factor(all_scens$Scenario)
  
  return(all_scens)
}