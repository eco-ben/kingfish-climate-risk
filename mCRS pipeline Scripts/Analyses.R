setwd('c://Users/grier/OneDrive/Documents/Honours/Lab Records/')
library(data.table)
library(ggplot2)
library(sf)
library(ggpubr)
library(maps)
library(dplyr)
library(scico)
library(nlme)
library(emmeans)
library(SIBER)
#### A) Load DeltaCRS Results ----
PremeanDT <- read.csv('Results/mCRS-DeltamCRS RAW/DeltaCRS_mCRS_RAW.csv')
PremeanDT$Depth <- as.factor(PremeanDT$Depth)
# CRS.dt <- read.csv('Results/DeltaCRS_PopValues.csv')
# CRS.dt$PopID <- as.factor(CRS.dt$PopID)
# CRS.dt$Depth <- as.factor(CRS.dt$Depth)
Coordinates <- read.csv('Required Files/Global_Coordinates.csv')

## Load Species level analysis results 
SpL <- read.csv('Results/mCRS-DeltamCRS RAW/DeltaCRS_mCRS_RAW_SpL.csv')

#Calculate annual averages for mCRS and DeltaCRS ID values
PremeanDT.Year <- data.table(PremeanDT)
PremeanDT.Year <- PremeanDT.Year[,.(DeltaCRS = mean(DeltaCRS),CRS = mean(CRS)),by=.(PopID,Depth,Extract.ID,Year, Scenario)]
PremeanDT.Year <- as.data.frame(PremeanDT.Year)
PremeanDT.Year <- dplyr::left_join(PremeanDT.Year, Coordinates[,c(9,1,2,3,4)], by = c('Extract.ID','PopID'))
PremeanDT.Year$Distance_from_equator <- ifelse(PremeanDT.Year$Alt.Latitude < 0, -PremeanDT.Year$Alt.Latitude, PremeanDT.Year$Alt.Latitude)
PremeanDT.Year$Depth <- as.factor(PremeanDT.Year$Depth)
PremeanDT.Year$PopID <- as.factor(PremeanDT.Year$PopID)

SpL <- data.table(SpL)
SpL <- SpL[,.(DeltaCRS = mean(DeltaCRS),CRS = mean(CRS)),by=.(PopID,Depth,Extract.ID,Year, Scenario)]
SpL <- as.data.frame(SpL)
SpL <- dplyr::left_join(SpL, Coordinates[,c(9,1,2,3,4)], by = c('Extract.ID','PopID'))
SpL$Distance_from_equator <- ifelse(SpL$Alt.Latitude < 0, -SpL$Alt.Latitude, SpL$Alt.Latitude)
SpL$Depth <- as.factor(SpL$Depth)
SpL$PopID <- as.factor(SpL$PopID)
PremeanDT.Year$Level <- 'Population'
SpL$Level <- 'Species'
Comb.Data <- rbind(PremeanDT.Year, SpL)

#### B) Load Environmental Data for IDs, Dates and Scenarios (From CRS Scripts) ----

load('Required Files/historical_stack_ensemble.RData')
load('Required Files/ssp26_stack_ensemble.RData')
load('Required Files/ssp45_stack_ensemble.RData')
load('Required Files/ssp85_stack_ensemble.RData')
ENV <- list(historical.stack.ensemble, ssp26.stack.ensemble, ssp45.stack.ensemble, ssp85.stack.ensemble)
ENV.summary <- rbindlist(ENV) #combine all data tables into one ENV
rm(ENV,historical.stack.ensemble, ssp26.stack.ensemble, ssp45.stack.ensemble,ssp85.stack.ensemble, DataCollation)
ENV.summary <- ENV.summary[!ENV.summary$Extract.ID %in% c('103', '204', '200', '111', '112','34'),]
ENV.summary.Y <- ENV.summary[,.(no3.mean = mean(no3), ph.mean = mean(ph), so.mean= mean(so), thetao.mean = mean(thetao)),by =.(Extract.ID,PopID,Year,Scenario,Depth)]
ENV.PopSummary <- ENV.summary.Y[,.(no3.mean = mean(no3.mean), ph.mean = mean(ph.mean), so.mean= mean(so.mean), thetao.mean = mean(thetao.mean)),by =.(PopID,Year,Scenario,Depth)]
ENV.PopSummary <- data.frame(ENV.PopSummary)
ENV.PopSummary_long <- tidyr::pivot_longer(ENV.PopSummary, cols = c(5:8), names_to = c('Variable'), values_to = 'Values')

#### C) Calculate Warming levels ----
#Calculating mean warming levels from environmental data
meanTemp <- ENV.PopSummary[,c(1,2,3,4,8)]
refmeanTemp <- meanTemp[meanTemp$Scenario == 'historical',] %>% group_by(PopID,Depth) %>% summarise(TempMean = mean(thetao.mean))
meanTemp <- left_join(meanTemp, refmeanTemp, by = c('PopID', 'Depth'))
meanTemp <- meanTemp %>% group_by(PopID,Depth,Year,Scenario) %>% reframe(DeltaTemp = (thetao.mean - TempMean))
meanTemp26 <- meanTemp[meanTemp$Scenario %in% c('historical','ssp26'),]
meanTemp45 <- meanTemp[meanTemp$Scenario %in% c('historical','ssp45'),]
meanTemp85 <- meanTemp[meanTemp$Scenario %in% c('historical','ssp85'),]
meanTemp26[meanTemp26$Scenario == 'historical',]$Scenario <- 'ssp26'
meanTemp45[meanTemp45$Scenario == 'historical',]$Scenario <- 'ssp45'
meanTemp85[meanTemp85$Scenario == 'historical',]$Scenario <- 'ssp85'
meanTemp <- rbind(meanTemp26,meanTemp45,meanTemp85)
meanTemp <- meanTemp[!meanTemp$Year <= 2000,]

##### Modelling mean Temperatures ----
meanTemp$PopID <- as.factor(meanTemp$PopID)
meanTemp$Depth <- as.factor(meanTemp$Depth)
meanTemp$Scenario <- as.factor(meanTemp$Scenario)
newTempGAM <- meanTemp[,c(1,2,3,4)]
newTempGAM$DeltaTModelled <- 0
ssps = list('ssp26', 'ssp45', 'ssp85')
depths = list(10,50,100,150)
pops = list('NE_Pacific','NW_Pacific','South_Africa','South_Pacific')

## Population Specific Warming levels --
for(i in 1:length(ssps)){
  for(j in 1:length(depths)){
    for(z in 1:length(pops)){
      model <- mgcv::gam(DeltaTemp ~ s(Year), data = meanTemp[meanTemp$PopID == pops[z]&
                                                                meanTemp$Depth == depths[j] &
                                                                meanTemp$Scenario == ssps[i],])
      newTempGAM[newTempGAM$Scenario == ssps[i]&
                   newTempGAM$PopID == pops[z]&
                   newTempGAM$Depth == depths[j],]$DeltaTModelled <- predict(model,
                                                                             newTempGAM[newTempGAM$Scenario == ssps[i]&
                                                                                       newTempGAM$PopID == pops[z]&
                                                                                       newTempGAM$Depth == depths[j],])
    }
  }
}

## Species level warming levels -- 
# ENV.summary.Y$Depth <- as.factor(ENV.summary.Y$Depth)
# Global.refTemp <- ENV.summary.Y[ENV.summary.Y$Scenario == 'historical',c(1,3,4,5,9)] %>% group_by(Depth) %>% 
#   summarise(TempMean = mean(thetao.mean))
# DeltaT.Glob <- ENV.summary.Y[,c(1,3,4,5,9)] %>% group_by(Depth, Year, Scenario) %>% summarise(thetao.mean = mean(thetao.mean))
# DeltaT.Glob <- left_join(DeltaT.Glob, Global.refTemp, by = c('Depth'))
# DeltaT.Glob <- DeltaT.Glob %>% group_by(Depth,Year,Scenario) %>% reframe(DeltaTemp = (thetao.mean - TempMean))
# DeltaT.Glob26 <- DeltaT.Glob[DeltaT.Glob$Scenario %in% c('historical','ssp26'),]
# DeltaT.Glob45 <- DeltaT.Glob[DeltaT.Glob$Scenario %in% c('historical','ssp45'),]
# DeltaT.Glob85 <- DeltaT.Glob[DeltaT.Glob$Scenario %in% c('historical','ssp85'),]
# DeltaT.Glob26[DeltaT.Glob26$Scenario == 'historical',]$Scenario <- 'ssp26'
# DeltaT.Glob45[DeltaT.Glob45$Scenario == 'historical',]$Scenario <- 'ssp45'
# DeltaT.Glob85[DeltaT.Glob85$Scenario == 'historical',]$Scenario <- 'ssp85'
# DeltaT.Glob <- rbind(DeltaT.Glob26, DeltaT.Glob45, DeltaT.Glob85)
# DeltaT.Glob <- DeltaT.Glob[DeltaT.Glob$Year > 2000,]
# 
# DeltaT.Glob$Depth <- as.factor(DeltaT.Glob$Depth)
# DeltaT.Glob$Scenario <- as.factor(DeltaT.Glob$Scenario)
# newTempGlob <- DeltaT.Glob[,c(1,2,3,4)]
# newTempGlob$DeltaTModelled <- 0
# ssps = list('ssp26', 'ssp45', 'ssp85')
# depths = list(10,50,100,150)
# 
# for(i in 1:length(ssps)){
#   for(j in 1:length(depths)){
#     model <- mgcv::gam(DeltaTemp ~ s(Year), data = DeltaT.Glob[DeltaT.Glob$Depth == depths[j] &
#                                                                 DeltaT.Glob$Scenario == ssps[i],])
#       newTempGAM.Glob[newTempGAM.Glob$Scenario == ssps[i]&
#                    newTempGAM.Glob$Depth == depths[j],]$DeltaTModelled <- predict(model,
#                                                                              newTempGAM.Glob[newTempGAM.Glob$Scenario == ssps[i]&
#                                                                                           newTempGAM.Glob$Depth == depths[j],])
#   }
# }

###### modelled warming years and results data ----
#Find years when warming levels are reached for a given PopID and Depth 
ssp26.0.75deg.Year <- as.data.frame(newTempGAM[newTempGAM$Scenario == 'ssp26',]) %>% filter(DeltaTModelled >= 0.75) %>% 
  group_by(PopID, Depth, Scenario) %>% 
  summarize(Warming.year = min(Year))
ssp26.0.75deg.Year$Warming.level <- 0.75
ssp45.1deg.Year <- as.data.frame(newTempGAM[newTempGAM$Scenario == 'ssp45',]) %>% filter(DeltaTModelled >= 1) %>% 
  group_by(PopID, Depth, Scenario) %>% 
  summarize(Warming.year = min(Year))
ssp45.1deg.Year$Warming.level <- 1
ssp85.1.75deg.Year <- as.data.frame(newTempGAM[newTempGAM$Scenario == 'ssp85',]) %>% filter(DeltaTModelled >= 1.75) %>% 
  group_by(PopID, Depth, Scenario) %>% 
  summarize(Warming.year = min(Year))
ssp85.1.75deg.Year$Warming.level = 1.75
Warming.levels <- rbind(ssp26.0.75deg.Year, ssp45.1deg.Year, ssp85.1.75deg.Year)
Warming.levels$Depth <- as.factor(Warming.levels$Depth)

Warmlevel.data <- left_join(PremeanDT.Year[!PremeanDT.Year$Scenario == 'historical',], Warming.levels, by = c('PopID','Depth','Scenario'))
Warmlevel.data <- Warmlevel.data[Warmlevel.data$Year == Warmlevel.data$Warming.year,] #df containing DeltaCRS for every coordinateat the calculated PopID/Depth Warming levels

# ### Analysis level comparison data 
# Comb.Data$Extract.ID <- as.factor(Comb.Data$Extract.ID)
# Comb.Data_WL <- left_join(Comb.Data_WL[!Comb.Data_WL$Scenario == 'historical',], Warming.levels, by = c('PopID','Depth','Scenario'))
# Comb.Data_WL <- Comb.Data[SpL.WY$Year == SpL.WY$Warming.year,] #df containing DeltaCRS for every coordinateat the calculated PopID/Depth Warming levels

#### D) Regional Risk Emergence Calculation ----
Novel.Year <- as.data.frame(PremeanDT.Year) %>% filter(CRS >= 0) %>% 
  group_by(Extract.ID, PopID, Depth, Scenario) %>% 
  summarize(Novel.Year = min(Year))
dummy.df <- as.data.frame(PremeanDT.Year) %>%  
  group_by(Extract.ID, PopID, Depth, Scenario) %>% 
  summarize(Novel.Year = min(Year))
dummy.df <- dummy.df[,-5]
Novel.Year <- left_join(dummy.df, Novel.Year, by = c('Extract.ID','PopID','Depth','Scenario'))
Novel.Year <- left_join(Novel.Year, Coordinates[,c(9,4,1,2)], by = c('Extract.ID','PopID'))
Novel.Year <- Novel.Year[!Novel.Year$Scenario == 'historical',]
NovelY.med <- Novel.Year %>% group_by(Scenario, Depth, PopID) %>% summarise(NY.med = mean(Novel.Year, na.rm = T))
write.csv(NovelY.med, file = 'Results/Statistical Results/Regional_R_emergence_med.csv',row.names = F)

Combined.NY <- as.data.frame(Comb.Data) %>% filter(CRS >= 0) %>% 
  group_by(Extract.ID, Level, Depth, Scenario) %>% 
  summarize(Novel.Year = min(Year))
dummy.df <- as.data.frame(Comb.Data) %>%  
  group_by(Extract.ID, Level, Depth, Scenario) %>% 
  summarize(Novel.Year = min(Year))
dummy.df <- dummy.df[,-5]
#### E) Statistical Analyses ----
##### regional risk analyses ----
###### regional risk emergence ----
# calculating the year when regional risk emerges at a location
# regional risk emergence linear regression
NY.lm.ssp26 <-nlme::gls(Novel.Year ~ Depth * PopID,
                        data = Novel.Year[Novel.Year$Scenario == 'ssp26',], na.action = na.omit)
NY.lm.ssp45 <-nlme::gls(Novel.Year ~ Depth * PopID,
                        data = Novel.Year[Novel.Year$Scenario == 'ssp45',], na.action = na.omit)
NY.lm.ssp85 <-nlme::gls(Novel.Year ~ Depth * PopID,
                        data = Novel.Year[Novel.Year$Scenario == 'ssp85',], na.action = na.omit)
ssp26_NYemm <- emmeans(NY.lm.ssp26, pairwise ~ Depth*PopID)
ssp45_NYemm <- emmeans(NY.lm.ssp45, pairwise ~ Depth*PopID)
ssp85_NYemm <- emmeans(NY.lm.ssp85, pairwise ~ Depth*PopID)
ssp26_NYpair <- pairs(ssp26_NYemm, simple = list('Depth','PopID'), adjust = 'bonferroni')
ssp45_NYpair <- pairs(ssp45_NYemm, simple = list('Depth','PopID'), adjust = 'bonferroni')
ssp85_NYpair <- pairs(ssp85_NYemm, simple = list('Depth','PopID'), adjust = 'bonferroni')

NY.AP.ssp26 <- as.data.frame(ssp26_NYpair$`simple contrasts for PopID`)
NY.AP.ssp26$Scenario <- 'ssp26'
NY.AP.ssp45 <- as.data.frame(ssp45_NYpair$`simple contrasts for PopID`)
NY.AP.ssp45$Scenario <- 'ssp45'
NY.AP.ssp85 <- as.data.frame(ssp85_NYpair$`simple contrasts for PopID`)
NY.AP.ssp85$Scenario <- 'ssp85'
NY.Ac.Pops <- rbind(NY.AP.ssp26,NY.AP.ssp45,NY.AP.ssp85)
NY.Ac.Pops$Significance <- ifelse(NY.Ac.Pops$p.value < 0.05, 'significant','not significant')
write.csv(NY.Ac.Pops, file = 'Results/Statistical Results/NYResults_Across_Populations_CustWL.csv', row.names = F)

NY.AD.ssp26 <- as.data.frame(ssp26_NYpair$`simple contrasts for Depth`)
NY.AD.ssp26$Scenario <- 'ssp26'
NY.AD.ssp45 <- as.data.frame(ssp45_NYpair$`simple contrasts for Depth`)
NY.AD.ssp45$Scenario <- 'ssp45'
NY.AD.ssp85 <- as.data.frame(ssp85_NYpair$`simple contrasts for Depth`)
NY.AD.ssp85$Scenario <- 'ssp85'
NY.Ac.Depths <- rbind(NY.AD.ssp26,NY.AD.ssp45,NY.AD.ssp85)
NY.Ac.Depths$Significance <- ifelse(NY.Ac.Depths$p.value < 0.05, 'significant','not significant')
write.csv(NY.Ac.Depths, file = 'Results/Statistical Results/NYResults_Across_Depths_CustWL.csv', row.names = F)

###### regional risk rate of change ----
# mCRS ~ Year linear regression
Year_lm.ssp26 <- nlme::gls(CRS ~ Year * Depth * PopID,
                           data = PremeanDT.Year[PremeanDT.Year$Scenario == 'ssp26'&PremeanDT.Year$Year>2000,],
                           correlation = corAR1(form = ~Year|PopID/Depth/Extract.ID)) #time series models for DeltaCRS
Year_lm.ssp45 <- nlme::gls(CRS ~ Year * Depth * PopID,
                           data = PremeanDT.Year[PremeanDT.Year$Scenario == 'ssp45'&PremeanDT.Year$Year>2000,],
                           correlation = corAR1(form = ~Year|PopID/Depth/Extract.ID))
Year_lm.ssp85 <- nlme::gls(CRS ~ Year * Depth * PopID,
                           data = PremeanDT.Year[PremeanDT.Year$Scenario == 'ssp85'&PremeanDT.Year$Year>2000,],
                           correlation = corAR1(form = ~Year|PopID/Depth/Extract.ID))
YearT.ssp26 <- emtrends(Year_lm.ssp26, ~Depth:PopID, var = 'Year') #Trends for Year by Depth/PopID
YearT.ssp45 <- emtrends(Year_lm.ssp45, ~Depth:PopID, var = 'Year')
YearT.ssp85 <- emtrends(Year_lm.ssp85, ~Depth:PopID, var = 'Year')
Y.eff.ssp26 <- as.data.frame(test(YearT.ssp26))
Y.eff.ssp26$Scenario <- 'ssp26'
Y.eff.ssp45 <- as.data.frame(test(YearT.ssp45))
Y.eff.ssp45$Scenario <- 'ssp45'
Y.eff.ssp85 <- as.data.frame(test(YearT.ssp85))
Y.eff.ssp85$Scenario <- 'ssp85'

Year.effect <- rbind(Y.eff.ssp26, Y.eff.ssp45, Y.eff.ssp85)
write.csv(Year.effect, file = 'Results/Statistical Results/Results_Yeartimeseries.csv', row.names = F)

# Year trend plots
sjPlot::plot_model(Year_lm.ssp26, type = 'pred', terms = c('Year','Depth','PopID'))
ggsave(filename = 'Results/Plots/YearModelplot_ssp26.png',height = 100, width = 300,dpi = 300, units = 'mm')
sjPlot::plot_model(Year_lm.ssp45, type = 'pred', terms = c('Year','Depth','PopID'))
ggsave(filename = 'Results/Plots/YearModelplot_ssp45.png',height = 100, width = 300,dpi = 300, units = 'mm')
sjPlot::plot_model(Year_lm.ssp85, type = 'pred', terms = c('Year','Depth','PopID'))
ggsave(filename = 'Results/Plots/YearModelplot_ssp85.png',height = 100, width = 300,dpi = 300, units = 'mm')

##### local risk analyses ----
Warmlevel.data$PopLab <- as.character(Warmlevel.data$PopID)
Warmlevel.data$PopLab <- ifelse(Warmlevel.data$PopLab == 'NE_Pacific', 'NEP', Warmlevel.data$PopLab)
Warmlevel.data$PopLab <- ifelse(Warmlevel.data$PopLab == 'NW_Pacific', 'NWP', Warmlevel.data$PopLab)
Warmlevel.data$PopLab <- ifelse(Warmlevel.data$PopLab == 'South_Africa', 'SA', Warmlevel.data$PopLab)
Warmlevel.data$PopLab <- ifelse(Warmlevel.data$PopLab == 'South_Pacific', 'SIP', Warmlevel.data$PopLab)
Warmlevel.data$PopLab <- as.factor(Warmlevel.data$PopLab)
Warmlevel.data$randoff <- sample(1:2448,2448, replace = F)
#adding random number to Long to create slightly unique points for each depth (found small amount to not change p values much)
Warmlevel.data$new.Long <- Warmlevel.data$Alt.Longitude + (0.00000000001*Warmlevel.data$randoff) 

ssp26_WL.lm <- nlme::gls(DeltaCRS ~ new.Long + Distance_from_equator * Depth * PopLab, #ssp26 regression including interaction with Latitude
                         data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp26',], 
                         correlation = corExp(form=~new.Long + Distance_from_equator))
ssp45_WL.lm <- nlme::gls(DeltaCRS ~ new.Long + Distance_from_equator * Depth * PopLab, #ssp45 regression including interaction with Latitude
                         data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp45',], 
                         correlation = corExp(form=~new.Long + Distance_from_equator))
ssp85_WL.lm <- nlme::gls(DeltaCRS ~ new.Long + Distance_from_equator * Depth * PopLab, #ssp85 regression including interaction with Latitude
                         data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp85',], 
                         correlation = corExp(form=~new.Long + Distance_from_equator))

###### Differences Across Pops and Depths ----
ssp26_emm <- emmeans(ssp26_WL.lm, pairwise ~ Depth*PopLab)
ssp45_emm <- emmeans(ssp45_WL.lm, pairwise ~ Depth*PopLab)
ssp85_emm <- emmeans(ssp85_WL.lm, pairwise ~ Depth*PopLab)
ssp26_pair <- pairs(ssp26_emm, simple = list('Depth','PopLab'), adjust = 'bonferroni')
ssp45_pair <- pairs(ssp45_emm, simple = list('Depth','PopLab'),adjust = 'bonferroni')
ssp85_pair <- pairs(ssp85_emm, simple = list('Depth','PopLab'),adjust='bonferroni')

AP.ssp26 <- as.data.frame(ssp26_pair$`simple contrasts for PopLab`)
AP.ssp26$Scenario <- 'ssp26'
AP.ssp45 <- as.data.frame(ssp45_pair$`simple contrasts for PopLab`)
AP.ssp45$Scenario <- 'ssp45'
AP.ssp85 <- as.data.frame(ssp85_pair$`simple contrasts for PopLab`)
AP.ssp85$Scenario <- 'ssp85'
Ac.Pops <- rbind(AP.ssp26,AP.ssp45,AP.ssp85)
Ac.Pops$Significance <- ifelse(Ac.Pops$p.value < 0.05, 'significant','not significant')
write.csv(Ac.Pops, file = 'Results/Statistical Results/Results_Across_Populations_CustWL.csv', row.names = F)

AD.ssp26 <- as.data.frame(ssp26_pair$`simple contrasts for Depth`)
AD.ssp26$Scenario <- 'ssp26'
AD.ssp45 <- as.data.frame(ssp45_pair$`simple contrasts for Depth`)
AD.ssp45$Scenario <- 'ssp45'
AD.ssp85 <- as.data.frame(ssp85_pair$`simple contrasts for Depth`)
AD.ssp85$Scenario <- 'ssp85'
Ac.Depths <- rbind(AD.ssp26,AD.ssp45,AD.ssp85)
Ac.Depths$Significance <- ifelse(Ac.Depths$p.value < 0.05, 'significant','not significant')
write.csv(Ac.Depths, file = 'Results/Statistical Results/Results_Across_Depths_CustWL.csv', row.names = F)

###### Latitude effect ----
Lat.emt.ssp26 <- (emtrends(ssp26_WL.lm, ~ Depth:PopLab, var = 'Distance_from_equator'))
Lat.eff.ssp26 <- as.data.frame(test(Lat.emt.ssp26))
Lat.emt.ssp45 <- emtrends(ssp45_WL.lm, ~ Depth:PopLab, var = 'Distance_from_equator')
Lat.eff.ssp45 <- as.data.frame(test(Lat.emt.ssp45))
Lat.emt.ssp85 <- emtrends(ssp85_WL.lm, ~ Depth:PopLab, var = 'Distance_from_equator')
Lat.eff.ssp85 <- as.data.frame(test(Lat.emt.ssp85))
Lat.eff.ssp26$Scenario <- 'ssp26'
Lat.eff.ssp45$Scenario <- 'ssp45'
Lat.eff.ssp85$Scenario <- 'ssp85'

Lat.Trends <- rbind(Lat.eff.ssp26, Lat.eff.ssp45, Lat.eff.ssp85)
Lat.Trends$Significance <- ifelse(Lat.Trends$p.value < 0.05, 'significant','not significant')
write.csv(Lat.Trends, file = 'Results/Statistical Results/Results_Latitude_CustWL.csv', row.names = F)

#Latitude trend plots
Lat_26 <- ggeffects::ggpredict(ssp26_WL.lm, terms = c('Distance_from_equator','Depth','PopLab'))
Lat_26 <- plot(Lat_26, rawdata = TRUE) + theme_classic() + theme(plot.title = element_blank(), 
                                                                 axis.text = element_text(size = 13), 
                                                                 axis.title = element_text(size = 14, 
                                                                                           face = 'bold'), 
                                                                 strip.text.x = element_text(size = 10)) + 
  ylab('\u0394mCRS') + xlab('Distance from equator (\u00b0)') + ggtitle('\u0394mCRS across Latitude - SSP1-2.6')
ggsave(plot = Lat_26, filename = 'Results/Plots/LatitudeModelplot_ssp26.png',height = 150, width = 250,dpi = 300, units = 'mm')
Lat_45 <- ggeffects::ggpredict(ssp45_WL.lm, terms = c('Distance_from_equator','Depth','PopLab'))
Lat_45 <- plot(Lat_45, rawdata = TRUE) + theme_classic() + theme(plot.title = element_blank(), 
                                                                 axis.text = element_text(size = 13), 
                                                                 axis.title = element_text(size = 14, 
                                                                                           face = 'bold'), 
                                                                 strip.text.x = element_text(size = 10)) + 
  ylab('\u0394mCRS') + xlab('Distance from equator (\u00b0)') + ggtitle('\u0394mCRS across Latitude - SSP2-4.5')
ggsave(plot = Lat_45, filename = 'Results/Plots/LatitudeModelplot_ssp45.png',height = 150, width = 250,dpi = 300, units = 'mm')
Lat_85 <- ggeffects::ggpredict(ssp85_WL.lm, terms = c('Distance_from_equator','Depth','PopLab'))
Lat_85 <- plot(Lat_85, rawdata = TRUE) + theme_classic() + theme(plot.title = element_blank(), 
                                                                 axis.text = element_text(size = 13), 
                                                                 axis.title = element_text(size = 14, 
                                                                                           face = 'bold'), 
                                                                 strip.text.x = element_text(size = 10)) + 
  ylab('\u0394mCRS') + xlab('Distance from equator (\u00b0)') + ggtitle('\u0394mCRS across Latitude - SSP5-8.5')
ggsave(plot = Lat_85, filename = 'Results/Plots/LatitudeModelplot_ssp85.png',height = 150, width = 250,dpi = 300, units = 'mm')

###### Species versus Population level analyses ----

ssp26.SpL.lm <- nlme::gls(CRS ~ Depth * Level,
                          data = Comb.Data[Comb.Data$Scenario == 'ssp26'&Comb.Data$Year == 2100,],)
ssp45.SpL.lm <- nlme::gls(CRS ~ Depth * Level,
                          data = Comb.Data[Comb.Data$Scenario == 'ssp45'&Comb.Data$Year == 2100,],)
ssp85.SpL.lm <- nlme::gls(CRS ~ Depth * Level,
                          data = Comb.Data[Comb.Data$Scenario == 'ssp85'&Comb.Data$Year == 2100,],)
ssp26.SpL_emm <- emmeans(ssp26.SpL.lm, ~ Depth:Level)
ssp45.SpL_emm <- emmeans(ssp45.SpL.lm, ~ Depth:Level)
ssp85.SpL_emm <- emmeans(ssp85.SpL.lm, ~ Depth:Level)
ssp26.SpL_pair <- as.data.frame(pairs(ssp26.SpL_emm, simple = list('Level')))
ssp26.SpL_pair$Scenario <- 'ssp26'
ssp45.SpL_pair <- as.data.frame(pairs(ssp45.SpL_emm, simple = list('Level')))
ssp45.SpL_pair$Scenario <- 'ssp45'
ssp85.SpL_pair <- as.data.frame(pairs(ssp85.SpL_emm, simple = list('Level')))
ssp85.SpL_pair$Scenario <- 'ssp85'
Analysis.Lvls <- rbind(ssp26.SpL_pair, ssp45.SpL_pair, ssp85.SpL_pair)
write.csv(Analysis.Lvls, file = 'Results/Statistical Results/Species_PopulationlvlAN.csv',row.names = F)

### F) Figures ----
##### modelled temperature plots ----
meanTemp$Depth <- as.factor(meanTemp$Depth)
# axes ranges
minTemp <- min(newTempGAM$DeltaTModelled)
maxTemp <- max(newTempGAM$DeltaTModelled)
#Pop names to labels for plotting
Warming.levels$PopLab <- as.character(Warming.levels$PopID)
Warming.levels$PopLab <- ifelse(Warming.levels$PopLab == 'NE_Pacific', 'NEP', Warming.levels$PopLab)
Warming.levels$PopLab <- ifelse(Warming.levels$PopLab == 'NW_Pacific', 'NWP', Warming.levels$PopLab)
Warming.levels$PopLab <- ifelse(Warming.levels$PopLab == 'South_Africa', 'SA', Warming.levels$PopLab)
Warming.levels$PopLab <- ifelse(Warming.levels$PopLab == 'South_Pacific', 'SIP', Warming.levels$PopLab)
Warming.levels$PopLab <- as.factor(Warming.levels$PopLab)
newTempGAM$PopLab <- as.character(newTempGAM$PopID)
newTempGAM$PopLab <- ifelse(newTempGAM$PopLab == 'NE_Pacific', 'NEP', newTempGAM$PopLab)
newTempGAM$PopLab <- ifelse(newTempGAM$PopLab == 'NW_Pacific', 'NWP', newTempGAM$PopLab)
newTempGAM$PopLab <- ifelse(newTempGAM$PopLab == 'South_Africa', 'SA', newTempGAM$PopLab)
newTempGAM$PopLab <- ifelse(newTempGAM$PopLab == 'South_Pacific', 'SIP', newTempGAM$PopLab)
newTempGAM$PopLab <- as.factor(newTempGAM$PopLab)
meanTemp$PopLab <- as.character(meanTemp$PopID)
meanTemp$PopLab <- ifelse(meanTemp$PopLab == 'NE_Pacific', 'NEP', meanTemp$PopLab)
meanTemp$PopLab <- ifelse(meanTemp$PopLab == 'NW_Pacific', 'NWP', meanTemp$PopLab)
meanTemp$PopLab <- ifelse(meanTemp$PopLab == 'South_Africa', 'SA', meanTemp$PopLab)
meanTemp$PopLab <- ifelse(meanTemp$PopLab == 'South_Pacific', 'SIP', meanTemp$PopLab)
meanTemp$PopLab <- as.factor(meanTemp$PopLab)

ggplot() +
  geom_line(data = newTempGAM[newTempGAM$Scenario == 'ssp26',], aes(x = Year, y = DeltaTModelled, color = as.factor(Depth))) +
  geom_line(data = meanTemp[meanTemp$Scenario == 'ssp26',],aes(x = Year, y = DeltaTemp, color = as.factor(Depth)), alpha = 0.2) +
  geom_segment(data = Warming.levels[Warming.levels$Scenario == 'ssp26',],aes(x = Warming.year,y = minTemp,xend = Warming.year, yend = Warming.level), alpha = 0.4) + 
  geom_segment(data = Warming.levels[Warming.levels$Scenario == 'ssp26'&Warming.levels$Depth ==150,],aes(x = 2001, y = Warming.level, xend = Warming.year, yend = Warming.level), alpha = 0.4) +
  geom_point(data = Warming.levels[Warming.levels$Scenario == 'ssp26',], aes(x = Warming.year, y = Warming.level), alpha = 0.5, size = 2) +
  scale_y_continuous(limits = c(minTemp, maxTemp)) + labs(color = 'Depth') + ylab('\u0394Temp (\u00b0C)') +
  facet_wrap(~PopLab, nrow = 1) + theme_classic() + ggtitle('Population Temperature change - SSP1-2.6') + 
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/Modelled_TempChange_ssp26.png',height = 150, width = 300, dpi = 300, units = 'mm')
ggplot() +
  geom_line(data = newTempGAM[newTempGAM$Scenario == 'ssp45',], aes(x = Year, y = DeltaTModelled, color = as.factor(Depth))) +
  geom_line(data = meanTemp[meanTemp$Scenario == 'ssp45',],aes(x = Year, y = DeltaTemp, color = as.factor(Depth)), alpha = 0.2) +
  geom_segment(data = Warming.levels[Warming.levels$Scenario == 'ssp45',],aes(x = Warming.year,y = minTemp,xend = Warming.year, yend = Warming.level), alpha = 0.4) + 
  geom_segment(data = Warming.levels[Warming.levels$Scenario == 'ssp45'&Warming.levels$Depth ==150,],aes(x = 2001, y = Warming.level, xend = Warming.year, yend = Warming.level), alpha = 0.4) +
  geom_point(data = Warming.levels[Warming.levels$Scenario == 'ssp45',], aes(x = Warming.year, y = Warming.level), alpha = 0.5, size = 2) +
  scale_y_continuous(limits = c(minTemp, maxTemp)) + labs(color = 'Depth') + ylab('\u0394Temp (\u00b0C)') +
  facet_wrap(~PopLab, nrow = 1) + theme_classic() + ggtitle('Population Temperature change - SSP2-4.5') + 
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/Modelled_TempChange_ssp45.png',height = 150, width = 300, dpi = 300, units = 'mm')
ggplot() +
  geom_line(data = newTempGAM[newTempGAM$Scenario == 'ssp85',], aes(x = Year, y = DeltaTModelled, color = as.factor(Depth))) +
  geom_line(data = meanTemp[meanTemp$Scenario == 'ssp85',],aes(x = Year, y = DeltaTemp, color = as.factor(Depth)), alpha = 0.2) +
  geom_segment(data = Warming.levels[Warming.levels$Scenario == 'ssp85',],aes(x = Warming.year,y = minTemp,xend = Warming.year, yend = Warming.level), alpha = 0.4) + 
  geom_segment(data = Warming.levels[Warming.levels$Scenario == 'ssp85'&Warming.levels$Depth ==150,],aes(x = 2001, y = Warming.level, xend = Warming.year, yend = Warming.level), alpha = 0.4) +
  geom_point(data = Warming.levels[Warming.levels$Scenario == 'ssp85',], aes(x = Warming.year, y = Warming.level), alpha = 0.5, size = 2) +
  scale_y_continuous(limits = c(minTemp, maxTemp)) + labs(color = 'Depth') + ylab('\u0394Temp (\u00b0C)') +
  facet_wrap(~PopLab, nrow = 1) + theme_classic() + ggtitle('Population Temperature change - SSP5-8.5') + 
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/Modelled_TempChange_ssp85.png',height = 150, width = 300, dpi = 300, units = 'mm')

##### Regional risk emergence plots ----
Novel.Year$PopLab <- Novel.Year$PopID
Novel.Year$PopLab <- ifelse(Novel.Year$PopLab == 'NE_Pacific', 'NEP', Novel.Year$PopLab)
Novel.Year$PopLab <- ifelse(Novel.Year$PopLab == 'NW_Pacific', 'NWP', Novel.Year$PopLab)
Novel.Year$PopLab <- ifelse(Novel.Year$PopLab == 'South_Africa', 'SA', Novel.Year$PopLab)
Novel.Year$PopLab <- ifelse(Novel.Year$PopLab == 'South_Pacific', 'SIP', Novel.Year$PopLab)
Novel.Year$PopLab <- as.factor(Novel.Year$PopLab)
#Novel.Year$Novel.Year <- ifelse(is.na(Novel.Year$Novel.Year),2100,Novel.Year$Novel.Year)
###### Across Depth Boxplots ----
####### 1) ssp26 ----
ggplot() +
  geom_boxplot(data = Novel.Year[Novel.Year$Scenario == 'ssp26',], 
               aes(x = Depth, y = Novel.Year,fill = Depth)) +
  geom_jitter(data = Novel.Year[Novel.Year$Scenario == 'ssp26',],
              aes(x = Depth, y = Novel.Year),color = 'grey', width = 0.3) +
  scale_y_continuous(limits = c(2000,2101))+ylab('Novel Year') +
  facet_wrap(~PopLab, nrow = 1) + theme_classic() + ggtitle('Novel Year Across Depths - SSP1-2.6') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NovelYearAD_ssp26.png', height = 150, width = 250, dpi = 300, units = 'mm')
####### 2) ssp45 ----
ggplot() +
  geom_boxplot(data = Novel.Year[Novel.Year$Scenario == 'ssp45',], 
               aes(x = Depth, y = Novel.Year,fill = Depth)) +
  geom_jitter(data = Novel.Year[Novel.Year$Scenario == 'ssp45',],
              aes(x = Depth, y = Novel.Year),color = 'grey', width = 0.3) +
  scale_y_continuous(limits = c(2000,2101))+ylab('Novel Year') +
  facet_wrap(~PopLab, nrow = 1) + theme_classic() + ggtitle('Novel Year Across Depths - SSP2-4.5') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NovelYearAD_ssp45.png', height = 150, width = 250, dpi = 300, units = 'mm')
####### 3) ssp85 ----
ggplot() +
  geom_boxplot(data = Novel.Year[Novel.Year$Scenario == 'ssp85',], 
               aes(x = Depth, y = Novel.Year,fill = Depth)) +
  geom_jitter(data = Novel.Year[Novel.Year$Scenario == 'ssp85',],
              aes(x = Depth, y = Novel.Year),color = 'grey', width = 0.3) +
  scale_y_continuous(limits = c(2000,2101))+ylab('Novel Year') +
  facet_wrap(~PopLab, nrow = 1) + theme_classic() + ggtitle('Novel Year Across Depths - SSP5-8.5') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NovelYearAD_ssp85.png', height = 150, width = 250, dpi = 300, units = 'mm')
###### Across Population Boxplots ----
####### 1) ssp26 ----
ggplot() +
  geom_boxplot(data = Novel.Year[Novel.Year$Scenario == 'ssp26',], 
               aes(x = PopLab, y = Novel.Year,fill = PopLab)) +
  geom_jitter(data = Novel.Year[Novel.Year$Scenario == 'ssp26',],
              aes(x = PopLab, y = Novel.Year),color = 'grey', width = 0.3) +
  scale_y_continuous(limits = c(2000,2101))+ylab('Novel Year') +
  facet_wrap(~Depth, nrow = 1) + theme_classic() + xlab('Population') + ggtitle('Novel Year Across Populations - SSP1-2.6') + 
  labs(fill = 'Population') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NovelYearAP_ssp26.png', height = 150, width = 250, dpi = 300, units = 'mm')
####### 2) ssp45 ----
ggplot() +
  geom_boxplot(data = Novel.Year[Novel.Year$Scenario == 'ssp45',], 
               aes(x = PopLab, y = Novel.Year,fill = PopLab)) +
  geom_jitter(data = Novel.Year[Novel.Year$Scenario == 'ssp45',],
              aes(x = PopLab, y = Novel.Year),color = 'grey', width = 0.3) +
  scale_y_continuous(limits = c(2000,2101))+ylab('Novel Year') +
  facet_wrap(~Depth, nrow = 1) + theme_classic() + xlab('Population') + ggtitle('Novel Year Across Populations - SSP2-4.5') + 
  labs(fill = 'Population') +  
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NovelYearAP_ssp45.png', height = 150, width = 250, dpi = 300, units = 'mm')
####### 3) ssp85 ----
ggplot() +
  geom_boxplot(data = Novel.Year[Novel.Year$Scenario == 'ssp85',], 
               aes(x = PopLab, y = Novel.Year,fill = PopLab)) +
  geom_jitter(data = Novel.Year[Novel.Year$Scenario == 'ssp85',],
              aes(x = PopLab, y = Novel.Year),color = 'grey', width = 0.3) +
  scale_y_continuous(limits = c(2000,2101))+ ylab('Novel Year') +
  facet_wrap(~Depth, nrow = 1) + theme_classic() + xlab('Population') + ggtitle('Novel Year Across Populations - SSP5-8.5') + 
  labs(fill = 'Population') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NovelYearAP_ssp85.png', height = 150, width = 250, dpi = 300, units = 'mm')
###### Regional risk maps ----
BB_list = list(list(rbind(c(50,-15), c(50,-40), c(0,-40), c(0,-15), c(50,-15))),
               list(rbind(c(285,10), c(200,10), c(200, 40), c(285,40), c(285,10))),
               list(rbind(c(115,25), c(150,25), c(150, 45), c(115,45), c(115,25))),
               list(rbind(c(100,-5),c(100,-60),c(300,-60),c(300,-5),c(100,-5))))
Boxes <- st_multipolygon(BB_list)
Boxes <- st_sfc(Boxes, crs = 4326)
Boxes_360 <- st_shift_longitude(Boxes)
World <- st_read('Required Files/GIS/ne_10m_land.shp')
World_360 <- sf::st_as_sf(map('world2', plot = FALSE, fill = TRUE))
Novel.Year_sf <- st_as_sf(Novel.Year, coords = c(6,7), crs = 4326)
Novel.Year_sf <- st_shift_longitude(Novel.Year_sf)

Test2 <- ggplot()+
  geom_sf(data = World_360) + 
  geom_sf(data = Novel.Year_sf[Novel.Year_sf$Depth == 10 & Novel.Year_sf$Scenario == 'ssp85',], aes(color = Novel.Year)) + 
  scale_color_viridis_c(limits = c(2001,2100)) + labs(color = 'Novel Year') + theme(legend.position = 'bottom')
leg2 <- ggpubr::get_legend(Test2)

mCRSMap <- function(Scenario){
  Pacific <- ggplot() + 
    geom_sf(data =  Novel.Year_sf[Novel.Year_sf$Scenario == Scenario,], 
            aes(color = Novel.Year), size = 1.2, alpha = 0.7) +
    geom_sf(data = World_360) +
    geom_sf(data = Boxes, fill = NA, linewidth = 0.2) +
    theme_void() + theme(legend.position = "none") + theme(panel.background = element_rect(fill = 'white')) +
    scale_color_viridis_c(limits = c(2001, 2100)) + 
    coord_sf(xlim = c(100,300),ylim = c(50,-60)) + facet_wrap(~Depth)
  SA <- ggplot() + 
    geom_sf(data =  Novel.Year_sf[Novel.Year_sf$Scenario == Scenario,], 
            aes(color = Novel.Year), size = 1.2, alpha = 0.7) +
    geom_sf(data = World) +
    theme_void() + theme(legend.position = "none") + theme(panel.background = element_rect(fill = 'white')) + 
    scale_color_viridis_c(limits = c(2001, 2100)) + 
    coord_sf(xlim = c(10,40),ylim = c(-20,-40)) + facet_wrap(~Depth)
  final <- ggarrange(Pacific, SA, heights = c(3,1), widths = c(2.5,1), labels = c('(a)','(b)'), vjust = 1)
  final <- ggarrange(final,leg2, ncol = 1, heights = c(6,1), widths = c(6,1))
  ggsave(final, filename = paste0('Results/Plots/NovelYear_Map_ ',Scenario,'_custWL.png'), bg = 'white',height = 150, 
         width = 300, dpi = 300, units = 'mm')
}
mCRSMap(Scenario = 'ssp26')
mCRSMap(Scenario = 'ssp45')
mCRSMap(Scenario = 'ssp85')
##### Local Risk magnitude plots ----
minDcrs <- min(Warmlevel.data$DeltaCRS)
maxDcrs <- max(Warmlevel.data$DeltaCRS)
###### Across Depth Boxplots ----
####### 1) 0.75deg warming (ssp26) ----
ggplot() +
  geom_boxplot(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp26',], aes(x = Depth, y = DeltaCRS, fill = Depth)) + 
  geom_jitter(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp26',], aes(x = Depth, y = DeltaCRS),width = 0.25, color = 'grey') + 
  facet_wrap(~PopLab,nrow = 1) + ggtitle('\u0394mCRS - SSP1-2.6 - 0.75\u00b0C - Across Depths') + theme_classic() + 
  scale_y_continuous(limits = c(minDcrs,maxDcrs),breaks = c(-2,-1,0,1,2)) + ylab('\u0394mCRS') + xlab('Depth') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/DeltaCRS_BoxPlot_acDepth_ssp26_custWL.png',height = 150, width = 250, dpi = 300, units = 'mm')
####### 2) 1deg warming (ssp45) ----
ggplot() + 
  geom_boxplot(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp45',], aes(x = Depth, y = DeltaCRS, fill = Depth)) + 
  geom_jitter(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp45',], aes(x = Depth, y = DeltaCRS),width = 0.25, color = 'grey') + 
  facet_wrap(~PopLab,nrow = 1) + ggtitle('\u0394mCRS - SSP2-4.5 - 1\u00b0C - Across Depths') + theme_classic() + 
  scale_y_continuous(limits = c(minDcrs,maxDcrs),breaks = c(-2,-1,0,1,2)) + ylab('\u0394mCRS') + xlab('Depth') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/DeltaCRS_BoxPlot_acDepth_ssp45_custWL.png',height = 150, width = 250, dpi = 300, units = 'mm')
####### 3) 1.75deg warming (ssp85) ----
ggplot() + 
  geom_boxplot(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp85',], aes(x = Depth, y = DeltaCRS, fill = Depth)) + 
  geom_jitter(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp85',], aes(x = Depth, y = DeltaCRS),width = 0.25, color = 'grey') + 
  facet_wrap(~PopLab,nrow = 1) + ggtitle('\u0394mCRS - SSP5-8.5 - 1.75\u00b0C - Across Depths') + theme_classic() + 
  scale_y_continuous(limits = c(minDcrs,maxDcrs),breaks = c(-2,-1,0,1,2)) + ylab('\u0394mCRS') + xlab('Depth') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10))  
ggsave('Results/Plots/DeltaCRS_BoxPlot_acDepth_ssp85_custWL.png',height = 150, width = 250, dpi = 300, units = 'mm')
###### Across Population Boxplots ----
####### 1) 0.75deg warming (ssp26) ----
ggplot() +
  geom_boxplot(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp26',], aes(x = PopLab, y = DeltaCRS, fill = PopLab)) +
  geom_jitter(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp26',], aes(x = PopLab, y = DeltaCRS),width = 0.25, color = 'grey') +
  facet_wrap(~Depth,nrow = 1) + ggtitle('\u0394mCRS - SSP1-2.6 - 0.75\u00b0C - across Populations') + theme_classic() + 
  scale_y_continuous(limits = c(minDcrs,maxDcrs),breaks = c(-2,-1,0,1,2)) + ylab('\u0394mCRS') + xlab('Population') + 
  labs(fill = 'Population') +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/DeltaCRS_BoxPlot_acPop_ssp26_custWL.png',height = 150, width = 250, dpi = 300, units = 'mm')
####### 2) 1deg warming (ssp45) ----
ggplot() +
  geom_boxplot(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp45',], aes(x = PopLab, y = DeltaCRS, fill = PopLab)) +
  geom_jitter(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp45',], aes(x = PopLab, y = DeltaCRS),width = 0.25, color = 'grey') +
  facet_wrap(~Depth,nrow = 1) + ggtitle('\u0394mCRS - SSP2-4.5 - 1\u00b0C - across Populations') + theme_classic() + 
  scale_y_continuous(limits = c(minDcrs,maxDcrs),breaks = c(-2,-1,0,1,2)) + ylab('\u0394mCRS') + xlab('Population') + 
  labs(fill = 'Population') +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/DeltaCRS_BoxPlot_acPop_ssp45_custWL.png',height = 150, width = 250, dpi = 300, units = 'mm')
####### 3) 1.75deg warming (ssp85) ----
ggplot() +
  geom_boxplot(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp85',], aes(x = PopLab, y = DeltaCRS, fill = PopLab)) +
  geom_jitter(data = Warmlevel.data[Warmlevel.data$Scenario == 'ssp85',], aes(x = PopLab, y = DeltaCRS),width = 0.25, color = 'grey') +
  facet_wrap(~Depth,nrow = 1) + ggtitle('\u0394mCRS - SSP5-8.5 - 1.75\u00b0C warming - across Populations') + theme_classic() + 
  scale_y_continuous(limits = c(minDcrs,maxDcrs),breaks = c(-2,-1,0,1,2)) + ylab('\u0394mCRS') + xlab('Population') + 
  labs(fill = 'Population') +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave('Results/Plots/DeltaCRS_BoxPlot_acPop_ssp85_custWL.png',height = 150, width = 250, dpi = 300, units = 'mm')
###### local risk maps ----
Warmlevel.data_sf <- st_as_sf(Warmlevel.data, coords = c(8,9), crs = 4326) # make sf object to use in DeltaCRS mapping 


BB_list = list(list(rbind(c(50,-15), c(50,-40), c(0,-40), c(0,-15), c(50,-15))),
               list(rbind(c(285,10), c(200,10), c(200, 40), c(285,40), c(285,10))),
               list(rbind(c(115,25), c(150,25), c(150, 45), c(115,45), c(115,25))),
               list(rbind(c(100,-5),c(100,-60),c(300,-60),c(300,-5),c(100,-5))))
Boxes <- st_multipolygon(BB_list)
Boxes <- st_sfc(Boxes, crs = 4326)
Boxes_360 <- st_shift_longitude(Boxes)
Warmlevel.data_360 <- st_shift_longitude(Warmlevel.data_sf)
World <- st_read('Required Files/GIS/ne_10m_land.shp')
World_360 <- sf::st_as_sf(map('world2', plot = FALSE, fill = TRUE))

Test <- ggplot() + 
  geom_sf(data = World_360) +
  geom_sf(data =  Warmlevel.data_360[Warmlevel.data_360$Depth == 10 & Warmlevel.data_360$Scenario == 'ssp85',], aes(color = DeltaCRS), size = 2, alpha = 0.6) +
  theme_void() + theme(legend.position = "right") + 
  scale_color_scico(palette = 'batlow', midpoint = 0, limits = c(minDcrs, maxDcrs)) + labs(color = '\u0394mCRS') + 
  theme(legend.position = 'bottom')
leg <- ggpubr::get_legend(Test)
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
    geom_sf(data =  Warmlevel.data_sf[Warmlevel.data_sf$Scenario == Scenario,], 
            aes(color = DeltaCRS), size = 1.2, alpha = 0.6) +
    geom_sf(data = World) +
    # geom_sf(data = Boxes, fill = NA, linewidth = 0.2) + geom_sf(data = SP_lines, linewidth = 0.2) +
    #ggtitle('DeltaCRS - 10m') + 
    theme_void() + theme(legend.position = "none") + theme(panel.background = element_rect(fill = 'white')) + 
    scale_color_scico(palette = 'batlow', midpoint = 0, limits = c(minDcrs, maxDcrs)) + 
    coord_sf(xlim = c(10,40),ylim = c(-20,-40)) + facet_wrap(~Depth)
  final <- ggarrange(Pacific, SA, heights = c(3,1), widths = c(2.5,1), labels = c('(a)','(b)'), vjust = 1)
  final <- ggarrange(final,leg, ncol = 1, heights = c(6,1), widths = c(6,1))
  ggsave(final, filename = paste0('Results/Plots/DeltaCRS_Map_ ',Scenario,'_custWL.png'), bg = 'white',height = 150, 
         width = 300, dpi = 300, units = 'mm')
}
DeltamCRSMap(Scenario = 'ssp26')
DeltamCRSMap(Scenario = 'ssp45')
DeltamCRSMap(Scenario = 'ssp85')

##### Species Level Figures ----
minComb.CRS <- min(Comb.Data[!Comb.Data$Scenario == 'historical',]$CRS)
maxComb.CRS <- max(Comb.Data[!Comb.Data$Scenario == 'historical',]$CRS)

ggplot() + geom_boxplot(data = Comb.Data[Comb.Data$Scenario == 'ssp26'&Comb.Data$Year == 2100,], aes(x = Depth, y = CRS, fill = Level)) + 
  scale_y_continuous(limits = c(minComb.CRS, maxComb.CRS), breaks = c(-4,-3,-2,-1,0,1)) +
  theme_classic() + ggtitle('Niche Analysis Level - SSP1-2.6') + ylab('mCRS') +
  labs(fill = 'Analysis Level') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NicheLevelAnalysis_ssp26.png', height = 150, width = 250, dpi = 300, units = 'mm')

ggplot() + geom_boxplot(data = Comb.Data[Comb.Data$Scenario == 'ssp45'&Comb.Data$Year == 2100,], aes(x = Depth, y = CRS, fill = Level)) + 
  scale_y_continuous(limits = c(minComb.CRS, maxComb.CRS), breaks = c(-4,-3,-2,-1,0,1)) +
  theme_classic() + ggtitle('Niche Analysis Level - SSP2-4.5') + ylab('mCRS') +
  labs(fill = 'Analysis Level') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NicheLevelAnalysis_ssp45.png', height = 150, width = 250, dpi = 300, units = 'mm')

ggplot() + geom_boxplot(data = Comb.Data[Comb.Data$Scenario == 'ssp85'&Comb.Data$Year == 2100,], aes(x = Depth, y = CRS, fill = Level)) + 
  scale_y_continuous(limits = c(minComb.CRS, maxComb.CRS), breaks = c(-4,-3,-2,-1,0,1)) +
  theme_classic() + ggtitle('Niche Analysis Level - SSP5-8.5') + ylab('mCRS') +
  labs(fill = 'Analysis Level') + 
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 14, face = 'bold'), 
        strip.text.x = element_text(size = 10)) 
ggsave(filename = 'Results/Plots/NicheLevelAnalysis_ssp85.png', height = 150, width = 250, dpi = 300, units = 'mm')

##### PCA Results Plots ---- 
# Copied code from 4_CRS computation scripts to create datasets of PCA results for each scenario and depth
env_startcol=9
env_endcol=12
refyear_start='1850-01-01'
refyear_end='2000-12-31'
# Coordinates <- read.csv('Results/Global_Coordinates.csv')
Coordinates$PopLab <- as.character(Coordinates$PopID)
Coordinates$PopLab <- ifelse(Coordinates$PopLab == 'NE_Pacific', 'NEP', Coordinates$PopLab)
Coordinates$PopLab <- ifelse(Coordinates$PopLab == 'NW_Pacific', 'NWP', Coordinates$PopLab)
Coordinates$PopLab <- ifelse(Coordinates$PopLab == 'South_Africa', 'SA', Coordinates$PopLab)
Coordinates$PopLab <- ifelse(Coordinates$PopLab == 'South_Pacific', 'SIP', Coordinates$PopLab)
Coordinates$PopLab <- as.factor(Coordinates$PopLab)
dist <- Coordinates[,c(9,10)]
dist$PopLab <- as.factor(dist$PopLab)
depths = list(10,50,100,150)
pops = list('NEP','NWP','SA','SIP')
Colours <- Polychrome::createPalette(217, c("#ff0000", "#00ff00", "#0000ff"))

for(j in 1:length(depths)){
      depth.poly <- data.frame()
      depth.pca <- data.frame()
      depth.centroids <- data.frame()
      env=ENV.summary[ENV.summary$Scenario %in% c('historical','ssp26'),]
      env <- as.data.frame(env)
      env <- env[env$Depth == depths[j],]
      env <- env[!env$Extract.ID %in% c('103', '204', '200', '111', '112','34'),] # ID locations removed as of thesis writing (not including NWP loc)
      env$Scenario <- ifelse(env$Date >= '2001-01-01' & env$Date <= '2014-12-31', 'ssp26', env$Scenario)
      env <- left_join(env, Coordinates[,c(9,1,2)], by = c('Extract.ID'))
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
      env.prcomp.data <- env.prcomp.data[,c(3:5,1,2)]
      env.prcomp.data$Depth <- unlist(depths[j])
      env.prcomp.data <- env.prcomp.data[complete.cases(env.prcomp.data),] #remove invomplete rows 
      dat <- droplevels(dist)
      sp <- droplevels(unique(dat$PopLab))
      depth.pca <- rbind(depth.pca, env.prcomp.data)
      results <- vector('list', length=length(sp))
      polygons <- vector('list',length=length(sp))
      pc.centr <- vector('list',length = length(sp))
      insets <- vector('list',length = length(depths))
      for (m in 1:length(sp)) {
        tmp <- unique(dat[dat$PopLab == sp[[m]],]) #IDs for relevant population
        tmp.data <- sqldf::sqldf("SELECT * from tmp LEFT JOIN
                               'env.prcomp.data' USING (`Extract.ID`)") #join tmp to PC results -> tmp.data
        tmp.data <- tmp.data[complete.cases(tmp.data),] #remove incomplete results
        env.coords <- droplevels(tmp.data[tmp.data$Date >= refyear_start & tmp.data$Date <= refyear_end, ]) #include only the reference period and relevant population
        env.coords <- env.coords[chull(env.coords[,5], env.coords[,6]),] #create convex hull around condition points
        ### regular centroid
        env.centroid <- c(mean(tmp.data[tmp.data$Scenario=='historical' & #create climate centroid points for reference period
                                          tmp.data$Date >= refyear_start &
                                          tmp.data$Date <= refyear_end, 5]),
                          mean(tmp.data[tmp.data$Scenario=='historical' &
                                          tmp.data$Date >= refyear_start &
                                          tmp.data$Date <= refyear_end, 6]))
        env.centroid <- data.frame(env.centroid[1],env.centroid[2])
        env.polygon <- sp::SpatialPolygons(list(
          sp::Polygons(list(sp::Polygon(env.coords[,c(5,6)])), #the polygon used is the convex hull of environmental points
                       'historical')
        ))
        env.polygon_sf <- sf::st_as_sf(env.polygon, coords = env.polygon@polygons[[1]]@Polygons[[1]]@coords)
        env.polygon_sf$PopLab <- unlist(sp[m])
        env.polygon_sf$Depth <- unlist(depths[j])
        env.polygon_sf$Scenario <- unlist(ssps[i])
        depth.poly <- rbind(depth.poly, env.polygon_sf)
        env.centroid$PopLab <- unlist(sp[m])
        depth.centroids <- rbind(depth.centroids, env.centroid)
      }
      depth.pca <- left_join(depth.pca, dist, by = c('Extract.ID'))
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
      
      depthplot <- ggplot() + geom_sf(data = depth.poly) +
        geom_point(data = depth.pca[depth.pca$Date <= refyear_end,],
                   aes(x = PC1, y = PC2, color = as.factor(Extract.ID)),alpha = 0.4) + 
        geom_point(data = depth.centroids,aes(x = env.centroid.1., y = env.centroid.2.), color = 'black', size = 1, alpha = 0.7) +
        scale_color_manual(values = Colours) +
        annotation_custom(grob = ggplotGrob(inset))+
        theme_classic() +
        theme(legend.position = 'none') + 
        facet_wrap(~PopLab) + ggtitle(paste0('Historical Niche - ',depths[j],'m')) + theme(axis.text = element_text(size = 14), 
                                                                        axis.title = element_text(size = 17, face = 'bold'), 
                                                                        strip.text.x = element_text(size = 12)) 
      ggsave(plot = depthplot, filename=paste0('Results/Plots/',depths[j],'m_Reference.png'), width = 300, height = 200, dpi = 300, units = 'mm')
}
