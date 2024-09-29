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
  