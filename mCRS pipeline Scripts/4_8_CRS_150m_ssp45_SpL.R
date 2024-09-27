### Working directory
setwd('/home/sun/Analyses/Kingfish/Global/')
library(data.table)
############ Call data and data manipulation #######################
####################################################################

load('RData/Compiled/historical_stack_ensemble.RData')
load('RData/Compiled/ssp26_stack_ensemble.RData')
load('RData/Compiled/ssp45_stack_ensemble.RData')
load('RData/Compiled/ssp85_stack_ensemble.RData')
ENV <- list(historical.stack.ensemble, ssp26.stack.ensemble, ssp45.stack.ensemble, ssp85.stack.ensemble)

##########################################################################################
# Environmental parameter summary
ENV.summary <- rbindlist(ENV) #combine all data tables into one ENV
ENV.summary <- as.data.frame((ENV.summary))

rm(ENV,historical.stack.ensemble, ssp26.stack.ensemble, ssp45.stack.ensemble, ssp85.stack.ensemble, DataCollation)

#Loading Population ID into the workspace
Populations <- read.csv('Data/Global_Coordinates.csv')
Populations <- Populations[,c(9,4)]
Populations$PopID <- as.factor(Populations$PopID)

#use only species level
Populations$SP <- as.factor('Seriola.l')
Populations <- Populations[,-2L]

### Climate risk score calculation ####
CRS.computation <- function(dist, env, scenario, env_startcol, env_endcol,
                            refyear_start, refyear_end, depth) {
  
  #filter env to only contain values for chosen depth level
  env <- env[env$Depth == depth,]
  
  #Remove  South Pacific and USA outlier IDs
  env <- env[!env$Extract.ID %in% c('103', '204', '200', '111', '112','34'),]
  
  #change the scenario label for dates between the reference period and actual ssp model data period 
  env$Scenario <- ifelse(env$Date >= '2001-01-01' & env$Date <= '2014-12-31', scenario, env$Scenario) 
  
  dat <- droplevels(dist)
  sp <- droplevels(unique(dat$SP)) #Identify unique populations
  results <- vector('list', length=length(sp)) #have results for both populations
  
  marginality <- vector('list', length=length(sp)) #marginality for both populations
  
  for (i in 1:length(marginality)) {
    marginality[[i]] <- vector('list', length=nlevels(as.factor(unique(env$Date)))) #Makes marginality[[i]] a list with the same number of entries as unique dates
  }
  
  radiallines <- vector('list', length=length(sp)) #radiallines for both populations
  
  for (i in 1:length(radiallines)) {
    radiallines[[i]] <-vector('list', length=nlevels(as.factor(unique(env$Date)))) #Makes radiallines[[i]] a list with the same entries as unique dates
  }
  
  env.dist <- vector('list', length=length(sp)) #env.dist for each population
  
  for (i in 1:length(env.dist)) {
    env.dist[[i]] <- vector('list', length=nlevels(as.factor(unique(env$Date)))) #makes env.dist[[i]] a list with the same entries as unique dates
  }
  
  env.prcomp <- env[env$Date >= refyear_start &
                      env$Date <= refyear_end, 
                    c(env_startcol:env_endcol)]#environmental conditions for the reference period
  env.prcomp <- prcomp(env.prcomp, scale=T) # summary prcomp - 80.60% - Ensemble proportion of variance explained by the first 2 PC axes
  
  env.prcomp.data <- as.data.frame(predict(object=env.prcomp,
                                           newdata=env[,c(env_startcol:env_endcol)]))# rotated environmental data
  
  env.prcomp.data <- env.prcomp.data[,c(1:2)] #only use PC1 and PC2
  env.prcomp.data$Extract.ID <- env$Extract.ID 
  env.prcomp.data$Date <- env$Date
  env.prcomp.data$Scenario <- env$Scenario
  env.prcomp.data <- env.prcomp.data[,c(3:5,1,2)]
  env.prcomp.data <- env.prcomp.data[complete.cases(env.prcomp.data),] #remove invomplete rows 
  
  for (i in 1:length(sp)) {
    for (j in 1:nlevels(as.factor(unique(env$Date)))) {
      cat('\nAnalysing', i, 'of', length(sp), 'populations -', as.character(sp[[i]]), 'in Date', (levels(as.factor(unique(env$Date)))[[j]])) #print analysis progress
      tmp <- unique(dat[dat$SP == sp[[i]],]) #IDs for relevant population
      
      tmp.data <- sqldf::sqldf("SELECT * from tmp LEFT JOIN
                               'env.prcomp.data' USING (`Extract.ID`)") #join tmp to PC results -> tmp.data
      tmp.data <- tmp.data[complete.cases(tmp.data),] #remove incomplete results
      
      env.coords <- droplevels(tmp.data[tmp.data$Date >= refyear_start & tmp.data$Date <= refyear_end, ]) #include only the reference period and relevant population
      env.coords <- env.coords[chull(env.coords[,5], env.coords[,6]),] #create convex hull around condition points
      
      env.centroid <- c(mean(tmp.data[tmp.data$Scenario=='historical' & #create climate centroid points for reference period
                                        tmp.data$Date >= refyear_start &
                                        tmp.data$Date <= refyear_end, 5]),
                        mean(tmp.data[tmp.data$Scenario=='historical' &
                                        tmp.data$Date >= refyear_start &
                                        tmp.data$Date <= refyear_end, 6]))
      
      
      # Finding the nearest point of the convex hull via a point (site)
      # https://gis.stackexchange.com/questions/154479/how-to-find-the-distance-from-a-point-to-a-polygon-via-another-point-in-r
      marginality[[i]][[j]] <- tmp.data[tmp.data$Date==levels(as.factor(unique(env$Date)))[[j]], c(5,6)] #assign PC coordinate values to marginality
      
      
      env.polygon <- sp::SpatialPolygons(list(
        sp::Polygons(list(sp::Polygon(env.coords[,c(5,6)])), #the polygon used is the convex hull of environmental points
                     'historical')
      ))
      
      
      marginality[[i]][[j]]$PC1_centroid <- env.centroid[[1]] 
      marginality[[i]][[j]]$PC2_centroid <- env.centroid[[2]]
      
      # placeholder values for coordinates outside hull
      marginality[[i]][[j]]$PC1_outside <- 0
      marginality[[i]][[j]]$PC2_outside <- 0
      
      for (k in 1:nrow(marginality[[i]][[j]])) {
        # Outside convex hull points that are an extension of the line that connects the historical envelope centroid and sample point
        marginality[[i]][[j]]$PC1_outside[[k]] <- marginality[[i]][[j]][,1L][[k]] +
          (marginality[[i]][[j]][,1L][[k]] - (env.centroid[1])) / sqrt((env.centroid[1] - marginality[[i]][[j]][,1L][[k]])^2 +
                                                                         (env.centroid[2] - marginality[[i]][[j]][,2L][[k]])^2) * 20
        marginality[[i]][[j]]$PC2_outside[[k]] <- marginality[[i]][[j]][,2L][[k]] +
          (marginality[[i]][[j]][,2L][[k]] - (env.centroid[2])) / sqrt((env.centroid[1] - marginality[[i]][[j]][,1L][[k]])^2 +
                                                                         (env.centroid[2] - marginality[[i]][[j]][,2L][[k]])^2) * 20
      }
      
      # Add radial lines from the historical envelope centroid to each sample point
      RadialLines <- vector('list', length=nrow(marginality[[i]][[j]]))
      
      for (k in 1:nrow(marginality[[i]][[j]])) {
        RadialLines[[k]] <- sp::Lines(list(
          sp::Line(rbind(env.centroid, marginality[[i]][[j]][k ,c(5,6)]))# centroids are the same regardless of climate trajectory b/c they are based on historical conditions in which conditions are unaffected by climate trajectories yet
        ), ID=tmp.data[tmp.data$Date==levels(as.factor(unique(env$Date)))[[j]],]$Extract.ID[[k]])
      }
      
      RadialLines <- sp::SpatialLines(RadialLines) 
      
      # Cut radial lines at polygon edge
      Centre_edge_lines <- rgeos::gIntersection(RadialLines, env.polygon, byid=T)
      
      # Extract edge points from lines
      # placeholder values for edge points
      marginality[[i]][[j]]$PC1_edge <- 0
      marginality[[i]][[j]]$PC2_edge <- 0
      marginality[[i]][[j]]$Date <- levels(as.factor(unique(env$Date)))[[j]]
      
      # Extract points
      for (k in 1:nrow(marginality[[i]][[j]])) {   # i,j-th list element, 2 and 4 are edge coords (1 and 3 are centre points)
        marginality[[i]][[j]]$PC1_edge[[k]] <- sp::coordinates(Centre_edge_lines)[[k]][[1]][2]
        marginality[[i]][[j]]$PC2_edge[[k]] <- sp::coordinates(Centre_edge_lines)[[k]][[1]][4]
      }
    }
    
    marginality[[i]] <- data.table::rbindlist(marginality[[i]])
    
    tmp <- marginality[[i]]
    tmp$Date <- as.Date(tmp$Date)
    tmp <- unique(tmp)
    tmp.data <- sqldf::sqldf("SELECT * from 'tmp.data' LEFT JOIN
                             tmp USING (PC1, PC2, Date)")
    
    # placeholder values for niche marginality and centrality
    tmp.data$dM <- 0
    tmp.data$dC <- 0
    tmp.data$dE <- 0
    tmp.data$CRS <- 0
    
    tmp.data$dM <- as.numeric(unlist(lapply(1:nrow(tmp.data), function(x) {
      edge <- tmp.data[x, c('PC1_edge','PC2_edge')]
      obs <- tmp.data[x, c('PC1','PC2')]
      
      dist.obj <- rbind(as.numeric(edge),
                        as.numeric(obs))
      
      dM.val <- as.numeric(dist(dist.obj, method='euclidean'))
      
      return(dM.val)
    })))
    
    tmp.data$dC <- as.numeric(unlist(lapply(1:nrow(tmp.data), function(x) {
      centre <- tmp.data[x, c('PC1_centroid','PC2_centroid')]
      obs <- tmp.data[x, c('PC1','PC2')]
      
      dist.obj <- rbind(as.numeric(centre),
                        as.numeric(obs))
      
      dC.val <- as.numeric(dist(dist.obj, method='euclidean'))
      
      return(dC.val)
    })))
    
    tmp.data$dE <- as.numeric(unlist(lapply(1:nrow(tmp.data), function(x) {
      centre <- tmp.data[x, c('PC1_centroid','PC2_centroid')]
      edge <- tmp.data[x, c('PC1_edge','PC2_edge')]
      
      dist.obj <- rbind(as.numeric(centre),
                        as.numeric(edge))
      
      dE.val <- as.numeric(dist(dist.obj, method='euclidean'))
      
      return(dE.val)
    })))
    
    tmp.data$CRS <- -log(tmp.data$dE/tmp.data$dC)
    
    results[[i]] <- tmp.data
    
    write.csv(tmp.data, file=paste0('Results/SpL/CRS/',
                                    scenario,'/',scenario,'_',
                                    as.character(sp[[i]]),'_',
                                    depth,'m_SpL.csv'),row.names=F)
    
  }
  
  results <- rbindlist(results)
  
  return(results)
}

#### ENSEMBLE ----
### ssp26depth = 10
# # 1)
# CRS.ssp26.10m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp26'),],
#                                  scenario='ssp26',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 10)
# save.image('RData/CRS_results/CRScomputation_ssp26_10m_SpL.RData')

# # 2)
# CRS.ssp26.50m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp26'),],
#                                  scenario='ssp26',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 50)
# save.image('RData/CRS_results/CRScomputation_ssp26_50m_SpL.RData')

##3)
# CRS.ssp26.100m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp26'),],
#                                  scenario='ssp26',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 100)
# save.image('RData/CRS_results/CRScomputation_ssp26_100m_SpL.RData')

##4)
# CRS.ssp26.150m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp26'),],
#                                  scenario='ssp26',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 150)
# save.image('RData/CRS_results/CRScomputation_ssp26_150m_SpL.RData')

# ### ssp45depth = 10
# ## 1)
# CRS.ssp45.10m.ensemble <- CRS.computation(dist=Populations,
#                                           env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp45'),],
#                                           scenario='ssp45',
#                                           env_startcol=9, env_endcol=12,
#                                           refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 10)
# save.image('RData/CRS_results/CRScomputation_ssp45_10m_SpL.RData')
# 
# ## 2)
# CRS.ssp45.50m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp45'),],
#                                  scenario='ssp45',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 50)
# save.image('RData/CRS_results/CRScomputation_ssp45_50m_SpL.RData')
# 
# ##3)
# CRS.ssp45.100m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp45'),],
#                                  scenario='ssp45',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 100)
# save.image('RData/CRS_results/CRScomputation_ssp45_100m_SpL.RData')
# 
# ##4)
CRS.ssp45.150m.ensemble <- CRS.computation(dist=Populations,
                                 env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp45'),],
                                 scenario='ssp45',
                                 env_startcol=9, env_endcol=12,
                                 refyear_start='1850-01-01', refyear_end='2000-12-31', depth = 150)

save.image('RData/CRS_results/CRScomputation_ssp45_150m_SpL.RData')

### ssp85
# ## 5
# CRS.ssp85.10m.ensemble <- CRS.computation(dist=Populations,
#                                           env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp85'),],
#                                           scenario='ssp85',
#                                           env_startcol=9, env_endcol=12,
#                                           refyear_start='1850-01-01', refyear_end= '2000-12-31', depth = 10)
# save.image('RData/CRS_results/CRScomputation_ssp85_10m_SpL.RData')

## 6)
# CRS.ssp85.50m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp85'),],
#                                  scenario='ssp85',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end= '2000-12-31', depth = 50)
#save.image('RData/CRS_results/CRScomputation_ssp85_50m_SpL.RData')

##7)
# CRS.ssp85.100m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp85'),],
#                                  scenario='ssp85',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end= '2000-12-31', depth = 100)
#save.image('RData/CRS_results/CRScomputation_ssp85_100m_SpL.RData')

##8)
# CRS.ssp85.150m.ensemble <- CRS.computation(dist=Populations,
#                                  env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp85'),],
#                                  scenario='ssp85',
#                                  env_startcol=9, env_endcol=12,
#                                  refyear_start='1850-01-01', refyear_end= '2000-12-31', depth = 150)
#save.image('RData/CRS_results/CRScomputation_ssp85_150m_SpL.RData')







