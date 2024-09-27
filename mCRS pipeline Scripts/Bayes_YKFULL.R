library("leaps")
library("lmtest")
library("nlme")
library("ape")
library("broom")
library("FRK")
library("purrr")
library("lattice")
library("ggplot2")
library("RColorBrewer")
library("dplyr")
library("gstat")
library("sp")
library("spacetime")
library('INLA')
library(data.table)
library("tidyr")
library(sf)


setwd('/home/sun/Analyses/Kingfish/Global')
load('RData/Compiled/historical_stack_ensemble.RData')
load('RData/Compiled/ssp26_stack_ensemble.RData')
load('RData/Compiled/ssp45_stack_ensemble.RData')
load('RData/Compiled/ssp85_stack_ensemble.RData')

ENV <- list(historical.stack.ensemble, ssp26.stack.ensemble, ssp45.stack.ensemble, ssp85.stack.ensemble)
ENV.summary <- rbindlist(ENV) #combine all data tables into one ENV
rm(ENV,historical.stack.ensemble, ssp26.stack.ensemble, ssp45.stack.ensemble,ssp85.stack.ensemble, DataCollation)

env_startcol=9
env_endcol=12
refyear_start='1850-01-01'
refyear_end='2000-12-31'
Coordinates <- read.csv('Data/Global_Coordinates.csv')
dist <- Coordinates[,c(9,4)]
dist$PopID <- as.factor(dist$PopID)

scenario = 'ssp26'
env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp26'),]


env <- as.data.frame(env)
env <- env[!env$Extract.ID %in% c('103', '204', '200', '111', '112','34','52'),]
env$Scenario <- ifelse(env$Date >= '2001-01-01' & env$Date <= '2014-12-31', scenario, env$Scenario)
env <- left_join(env, Coordinates[,c(9,1,2)], by = c('Extract.ID'))
env.pr <- env[env$Date >= refyear_start &
                env$Date <= refyear_end,]#environmental conditions for the reference period
env.prcomp <- env.pr[,c(env_startcol:env_endcol)]
env.prcomp <- prcomp(env.prcomp, scale=T) # summary prcomp - 80.60% - Ensemble proportion of variance explained by the first 2 PC axes
env.prcomp.data <- as.data.frame(predict(object=env.prcomp,
                                         newdata=env[,c(env_startcol:env_endcol)]))# rotated environmental data
env.prcomp.data <- env.prcomp.data[,c(1:2)] #only use PC1 and PC2
env.prcomp.data$Extract.ID <- env$Extract.ID 
env.prcomp.data$Date <- env$Date
env.prcomp.data$Scenario <- env$Scenario
env.prcomp.data$Depth <- env$Depth
env.prcomp.data$Month <- as.numeric(format(env.prcomp.data$Date,'%m'))
env.prcomp.data <- env.prcomp.data[,c(3:6,1,2)]
env.prcomp.data <- env.prcomp.data[complete.cases(env.prcomp.data),] #remove invomplete rows 
env.prcomp.data <- left_join(env.prcomp.data, Coordinates[,c(9,1,2,4)], by = c('Extract.ID'))
env.sptcor <- env.prcomp.data


# env.sptcor <- env.sptcor[env.sptcor$Date < '2001-01-01',]
# env.sptcor <- data.table(env.sptcor)
# env.sptcor$Year <- format(as.Date(env.sptcor$Date, format="%Y/%m/%d"),"%Y")
# env.sptcor <- env.sptcor[,.(PC1 = mean(PC1), PC2 = mean(PC2)), by = .(Extract.ID, Year, PopID, Depth, Alt.Longitude, Alt.Latitude)]

env.sptcor$MonthT <- as.numeric(env.sptcor$Date)
Populations <- c('NE_Pacific','NW_Pacific','South_Africa','South_Pacific')
Depths <- c(10, 50, 100, 150)
INLA.Res.depths <- vector('list', length = 4)
INLA.Res <- vector('list', length = 4)

for(i in 1:length(Populations)){
  for(j in 1:length(Depths)){
    cat('\nAnalysing Population', i, 'of', length(Populations), 'and Depth', j, 'of', length(Depths))
    tst <- env.sptcor[env.sptcor$PopID == Populations[i] & env.sptcor$Depth == Depths[j],]
    
    
    coords <- unique(tst[c("Extract.ID", "Alt.Longitude", "Alt.Latitude")])
    boundary <- inla.nonconvex.hull(as.matrix(coords[, 2:3]))
    tst_coords <- st_as_sf(tst, coords = c('Alt.Longitude','Alt.Latitude'))
    max.edge = diff(range(st_coordinates(tst_coords)[,1]))/(3*5)
    YKmesh <- inla.mesh.2d(boundary = boundary,
                           max.edge = c(1,2)*max.edge, # max. edge length
                           cutoff = max.edge/5) 
    # plot(YKmesh, asp = 1, main = "")
    # lines(coords[c("Alt.Longitude", "Alt.Latitude")], col = "red", type = "p")
    # spde <- inla.spde2.pcmatern(mesh = YKmesh,
    #                             alpha = 2,
    #                             prior.range = c(1, 0.01),
    #                             prior.sigma = c(4, 0.01)
    # )
    
    spde <- inla.spde2.matern(mesh = YKmesh)
    
    n_months <- length(unique(tst$MonthT)) #number of unique years
    n_spatial <- YKmesh$n #number of cells in the grid 
    s_index <- inla.spde.make.index(name = "spatial.field",
                                    n.spde = n_spatial,
                                    n.group = n_months)
    
    coords.allyear <- tst[c("Alt.Longitude", "Alt.Latitude")] %>%
      as.matrix() #make a matrix of all coordinates
    PHI <- inla.spde.make.A(mesh = YKmesh,
                            loc = coords.allyear,
                            group = as.numeric(as.factor(tst$MonthT)),
                            n.group = n_months)
    dim(PHI)
    nrow(tst)
    length(s_index$spatial.field)
    tst$Depth <- as.factor(tst$Depth)
    tst$BL.PC1 <- NA
    tst$BL.PC2 <- NA
    stack.est.PC1 = inla.stack(data=list(PC1=tst$PC1),
                               A=list(PHI, 1),
                               effects=list(s_index,
                                            list(Intercept = rep(1, nrow(tst)))),
                               tag="est")
    
    stack.pred.PC1 = inla.stack(data = list(PC1 = tst$BL.PC1),
                                A = list(PHI, 1),
                                effects=list(s_index,
                                             list(Intercept = rep(1, nrow(tst)))),
                                tag='pred')
    
    stack.PC1 = inla.stack(stack.est.PC1, stack.pred.PC1)
    
    formula.PC1 <- (PC1 ~ -1 + Intercept + f(spatial.field,
                                             model = spde,
                                             group = spatial.field.group,
                                             control.group=list(model="ar1")))
    
    mod.mode.PC1 = inla(formula.PC1,
                        data=inla.stack.data(stack.est.PC1, spde=spde),
                        family="gaussian",
                        control.predictor=list(A=inla.stack.A(stack.est.PC1), compute=FALSE), 
                        control.inla = list(int.strategy = 'eb'),
                        num.threads = 1)
    
    mod.PC1 = inla(formula.PC1, data=inla.stack.data(stack.PC1, spde=spde),family="gaussian",
                   control.predictor=list(A=inla.stack.A(stack.PC1), compute=TRUE),
                   control.mode=list(theta=mod.mode.PC1$mode$theta, restart=FALSE), 
                   control.inla = list(int.strategy = 'eb'),
                   num.threads = 1)
    
    index.pred.PC1 = inla.stack.index(stack.PC1,"pred")$data
    tst$pred.PC1 <- mod.PC1$summary.linear.predictor[index.pred.PC1,'mean']
    
    
    stack.est.PC2 = inla.stack(data=list(PC2=tst$PC2),
                               A=list(PHI, 1),
                               effects=list(s_index,
                                            list(Intercept = rep(1, nrow(tst)))),
                               tag="est")
    
    stack.pred.PC2 = inla.stack(data = list(PC2 = tst$BL.PC2),
                                A = list(PHI, 1),
                                effects=list(s_index,
                                             list(Intercept = rep(1, nrow(tst)))),
                                tag='pred')
    
    stack.PC2 = inla.stack(stack.est.PC2, stack.pred.PC2)
    
    formula.PC2 <- (PC2 ~ -1 + Intercept + f(spatial.field,
                                             model = spde,
                                             group = spatial.field.group,
                                             control.group=list(model="ar1")))
    
    mod.mode.PC2 = inla(formula.PC2,
                        data=inla.stack.data(stack.est.PC2, spde=spde),
                        family="gaussian",
                        control.predictor=list(A=inla.stack.A(stack.est.PC2), compute=FALSE), 
                        control.inla = list(int.strategy = 'eb'),
                        num.threads = 1)
    
    mod.PC2 = inla(formula.PC2, data=inla.stack.data(stack.PC2, spde=spde),family="gaussian",
                   control.predictor=list(A=inla.stack.A(stack.PC2), compute=TRUE),
                   control.mode=list(theta=mod.mode.PC2$mode$theta, restart=FALSE), 
                   control.inla = list(int.strategy = 'eb'),
                   num.threads = 1)
    
    index.pred.PC2 = inla.stack.index(stack.PC2,"pred")$data
    tst$pred.PC2 <- mod.PC2$summary.linear.predictor[index.pred.PC2,'mean']
    
    save(tst, file = paste0('Results/INLA/INLA_Res_',Populations[i],'_',
                                                  Depths[j],'m.RData'))
  }
}






