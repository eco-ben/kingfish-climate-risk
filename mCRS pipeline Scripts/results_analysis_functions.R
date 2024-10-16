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

PCA_plotting = function(dist, env, scenario, env_startcol, env_endcol, refyear_start, refyear_end, plot_title, method=""){
  all_pca = data.frame()
  all_centroids = data.frame()
  all_polys = data.frame()
  
  env <- as.data.frame(env)
  env <- env[!env$Extract.ID %in% c('103', '204', '200', '111', '112','34','52'),] # ID locations removed as of thesis writing (not including NWP loc)
  env$Scenario <- ifelse(env$Date >= '2001-01-01' & env$Date <= '2014-12-31', scenario, env$Scenario)
  env <- left_join(env, coordinates[,c(9,1,2)], by = c('Extract.ID'))
  env.prcomp <- env[env$Date >= refyear_start &
                      env$Date <= refyear_end, 
                    c(env_startcol:env_endcol)]#environmental conditions for the reference period
  env.prcomp <- prcomp(env.prcomp, scale=T) # summary prcomp - 71+% - Ensemble proportion of variance explained by the first 2 PC axes
  env.prcomp.data <- as.data.frame(predict(object=env.prcomp,
                                           newdata=env[,c(env_startcol:env_endcol)]))# rotated environmental data
  env.prcomp.data <- env.prcomp.data[,c(1:2)] #only use PC1 and PC2
  env.prcomp.data$Extract.ID <- env$Extract.ID 
  env.prcomp.data$Date <- env$Date
  env.prcomp.data$Scenario <- unlist(env$Scenario)
  env.prcomp.data$Month <- as.numeric(format(env.prcomp.data$Date,'%m'))
  env.prcomp.data$Depth <- unlist(env$Depth)
  env.prcomp.data <- env.prcomp.data[,c(3:5,7,1,2)]
  env.prcomp.data <- env.prcomp.data[complete.cases(env.prcomp.data),] #remove invomplete rows 
  sp <- unique(dist$PopID)
  all_pca <- env.prcomp.data
  all_pca = left_join(all_pca, dist, by='Extract.ID')
  polygons <- vector('list',length=length(sp))
  pc.centr <- vector('list',length = length(sp))
  for (m in 1:length(sp)) {
    tmp <- unique(dist[dist$PopID == sp[[m]],]) #IDs for relevant population
    tmp.data <- sqldf::sqldf("SELECT * from tmp LEFT JOIN
                               'env.prcomp.data' USING (`Extract.ID`)") #join tmp to PC results -> tmp.data
    tmp.data <- tmp.data[complete.cases(tmp.data),] #remove incomplete results
    env.coords <- droplevels(tmp.data[tmp.data$Date >= refyear_start & tmp.data$Date <= refyear_end, ]) #include only the reference period and relevant population
    #env.coords <- env.coords[chull(env.coords[,5], env.coords[,6]),] #create convex hull around condition points
    
    env.poly <- distfree.cr(env.coords[,c("PC1","PC2")],alpha = 0.05, draw = F)
    env.poly.pts <- data.frame(env.poly$polygon)
    FL <- env.poly.pts[1,]
    env.poly.pts <- rbind(env.poly.pts, FL)
    env.poly <- sp::SpatialPolygons(list(
      sp::Polygons(list(sp::Polygon(env.poly.pts)), #the polygon used is the convex hull of environmental points
                   'historical')))
    env.polygon_sf <- sf::st_as_sf(env.poly, coords = env.polygon@polygons[[1]]@Polygons[[1]]@coords)
    env.polygon_sf$PopID <- unlist(sp[m])
    
    if(method == "surface"){
      spat.pcs = Full[Full$PopID == sp[m] & Full$Date <= refyear_end & Full$Depth == unique(env$Depth),]
    }
    else{
      spat.pcs <- Full[Full$PopID == sp[m] & Full$Date <= refyear_end,]
    }

    spat.pcs <- spat.pcs[complete.cases(spat.pcs),]
    env.centroid <- c(mean(spat.pcs$pred.PC1),mean(spat.pcs$pred.PC2))
    env.centroid <- data.frame(env.centroid[1],env.centroid[2])
    env.centroid$PopID = sp[m]
    all_polys = rbind(all_polys, env.polygon_sf)
    all_centroids = rbind(all_centroids, env.centroid)
  }
  
  inset <- ggplot() +
    geom_segment(data=as.data.frame(env.prcomp$rotation),
                 aes(x=0, y=0, xend=(PC1), yend=(PC2)),
                 arrow=arrow(length=unit(0.01, 'cm')), color='coral1') +
    annotate(geom='text',
             x=(as.data.frame(env.prcomp$rotation)[1,1]),
             y=(as.data.frame(env.prcomp$rotation)[1,2]),
             size=4,
             label=paste0((row.names(env.prcomp$rotation)[1]))) +
    annotate(geom='text',
             x=(as.data.frame(env.prcomp$rotation)[2,1]),
             y=(as.data.frame(env.prcomp$rotation)[2,2]),
             size=4,
             label=paste0((row.names(env.prcomp$rotation)[2]))) +
    annotate(geom='text',
             x=(as.data.frame(env.prcomp$rotation)[3,1]),
             y=(as.data.frame(env.prcomp$rotation)[3,2]),
             size=4,
             label=paste0((row.names(env.prcomp$rotation)[3]))) +
    annotate(geom='text',
             x=(as.data.frame(env.prcomp$rotation)[4,1]),
             y=(as.data.frame(env.prcomp$rotation)[4,2]),
             size=4,
             label=paste0((row.names(env.prcomp$rotation)[4]))) +
    theme_void() + scale_x_continuous(limits = c((min(as.data.frame(env.prcomp$rotation)$PC1)-1.5), 
                                                 (max(as.data.frame(env.prcomp$rotation)$PC2)+1.5))) + 
    scale_y_continuous(limits = c((min(as.data.frame(env.prcomp$rotation)$PC2)-1.5), 
                                  (max(as.data.frame(env.prcomp$rotation)$PC2)+1.5))) 
  
  all_pca$Year = as.numeric(format(as.Date(all_pca$Date,'%Y-m-%d'),'%Y'))
  year_pca = all_pca %>% group_by(Extract.ID,Depth, PopID, Scenario, Year) %>%
    summarise(PC1 = mean(PC1), PC2 = mean(PC2))
  
  pca_plot = ggplot() + 
    geom_point(data = year_pca,
               aes(x = PC1, y = PC2, color = as.factor(Depth)),alpha = 0.4) + 
    geom_point(data = all_centroids,aes(x = env.centroid.1., y = env.centroid.2.), color = 'black', size = 1, alpha = 0.7) +
    geom_sf(data = all_polys, alpha =0.3) +
    annotation_custom(grob = ggplotGrob(inset))+
    theme_classic() +
    theme(legend.position = 'none') + 
    facet_wrap(~PopID) + ggtitle(paste0(plot_title)) + theme(axis.text = element_text(size = 14), 
                                                                                axis.title = element_text(size = 17, face = 'bold'), 
                                                                                strip.text.x = element_text(size = 12)) 
  return(pca_plot)
}

DeltamCRSMap <- function(Scenario){
  Pacific <- ggplot() + 
    geom_sf(data =  Warmlevel.data_360[Warmlevel.data_360$Scenario == Scenario,], 
            aes(color = DeltaCRS), size = 1.2, alpha = 0.6) +
    geom_sf(data = World_360) +
    geom_sf(data = Boxes, fill = NA, linewidth = 0.2) +
    #ggtitle('DeltaCRS - 10m') + 
    theme_void() + theme(legend.position = "none") + theme(panel.background = element_rect(fill = 'white')) +
    scale_color_scico(palette = 'batlow', midpoint = 0, limits = c(minDcrs, maxDcrs)) + 
    coord_sf(xlim = c(100,300),ylim = c(50,-60)) + facet_wrap(~Depth)
  SA <- ggplot() + 
    geom_sf(data =  all_depths_sf[all_depths_sf$Scenario == Scenario,], 
            aes(color = DeltaCRS), size = 1.2, alpha = 0.6) +
    geom_sf(data = World) +
    theme_void() + theme(legend.position = "none") + theme(panel.background = element_rect(fill = 'white')) + 
    scale_color_scico(palette = 'batlow', midpoint = 0, limits = c(minDcrs, maxDcrs)) + 
    coord_sf(xlim = c(10,40),ylim = c(-20,-40)) + facet_wrap(~Depth)
  final <- ggarrange(Pacific, SA, heights = c(3,1), widths = c(2.5,1), labels = c('(a)','(b)'), vjust = 1)
  final <- ggarrange(final,leg, ncol = 1, heights = c(6,1), widths = c(6,1))

  return(final)
}

append_historical = function(x){
  # Attach historical data to each scenario and relabel
  ssp26 <- x[x$Scenario %in% c('historical','ssp26'),]
  ssp45 <- x[x$Scenario %in% c('historical','ssp45'),]
  ssp85 <- x[x$Scenario %in% c('historical','ssp85'),]
  ssp26[ssp26$Scenario == 'historical',]$Scenario <- 'ssp26'
  ssp45[ssp45$Scenario == 'historical',]$Scenario <- 'ssp45'
  ssp85[ssp85$Scenario == 'historical',]$Scenario <- 'ssp85'
  
  # Combine all scenarios' data 
  all_scens <- rbind(ssp26,ssp45,ssp85)
  all_scens$PopID <- as.factor(all_scens$PopID)
  all_scens$Depth <- as.factor(all_scens$Depth)
  all_scens$Scenario <- as.factor(all_scens$Scenario)
  
  return(all_scens)
}

normalise = function(x){
  x = (x-min(x)) / (max(x) - min(x))

  return(x)  
}