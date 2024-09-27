setwd("/home/sun/Analyses/Kingfish/")

#Loading Coordinates
Coordinates <- read.csv('Global/Data/Global_Coordinates.csv')

#### Environmental data merger ####
### Data preparation
source('Base_Codes/Functions/env_extract_ecoregionSummary.R')

### Historical period
# # 1) CESM2-WACCM
# historical.CESM2_WACCM.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#                                                            coords=Coordinates[,c(1,2)],
#                                                            time.period='historical',
#                                                            time1='1850-01-01',
#                                                            time2='2014-12-31')
# historical.CESM2_WACCM.stack$Scenario <- 'historical'
# 
# save.image('Global/RData/historical_summary_CESM2_WACCM.RData')
# rm(historical.CESM2_WACCM.stack)

# # # 2) CMCC-ESM2
# historical.CMCC_ESM2.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#                                                            coords=Coordinates[,c(1,2)],
#                                                            time.period='historical',
#                                                            time1='1850-01-01',
#                                                            time2='2014-12-31')
# historical.CMCC_ESM2.stack$Scenario <- 'historical'
# 
# save.image('Global/RData/historical_summary_CMCC_ESM2.RData')
# rm(historical.CMCC_ESM2.stack)

# # 3) GFDL-ESM4.gr
# historical.GFDL_ESM4_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#                                                          coords=Coordinates[,c(1,2)],
#                                                          time.period='historical',
#                                                          time1='1850-01-01',
#                                                          time2='2014-12-31')
# historical.GFDL_ESM4_gr.stack$Scenario <- 'historical'
# 
# save.image('Global/RData/historical_summary_GFDL_ESM4_gr.RData')
# rm(historical.GFDL_ESM4_gr.stack)

# # 4) IPSL-CM6A-LR
# historical.IPSL_CM6A_LR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#                                                          coords=Coordinates[,c(1,2)],
#                                                          time.period='historical',
#                                                          time1='1850-01-01',
#                                                          time2='2014-12-31')
# historical.IPSL_CM6A_LR.stack$Scenario <- 'historical'
# 
# save.image('Global/RData/historical_summary_IPSL_CM6A_LR.RData')
# rm(historical.IPSL_CM6A_LR.stack)

# # 5) MPI-ESM1-2-HR
# historical.MPI_ESM1_2_HR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#                                                          coords=Coordinates[,c(1,2)],
#                                                          time.period='historical',
#                                                          time1='1850-01-01',
#                                                          time2='2014-12-31')
# historical.MPI_ESM1_2_HR.stack$Scenario <- 'historical'
# 
# save.image('Global/RData/historical_summary_MPI_ESM1_2_HR.RData')
# rm(historical.MPI_ESM1_2_HR.stack)

# # 6) NorESM2-LM.gr
# historical.NorESM2_LM_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#                                                          coords=Coordinates[,c(1,2)],
#                                                          time.period='historical',
#                                                          time1='1850-01-01',
#                                                          time2='2014-12-31')
# historical.NorESM2_LM_gr.stack$Scenario <- 'historical'
# 
# save.image('Global/RData/historical_summary_NorESM2_LM_gr.RData')
# rm(historical.NorESM2_LM_gr.stack)


# ### SSP26
# # 7) CESM2-WACCM
# ssp26.CESM2_WACCM.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#                                                            coords=Coordinates[,c(1,2)],
#                                                            time.period='ssp26',
#                                                            time1='2015-01-01',
#                                                            time2='2100-12-31')
# ssp26.CESM2_WACCM.stack$Scenario <- 'ssp26'
# 
# save.image('Global/RData/ssp26_summary_CESM2_WACCM.RData')
# rm(ssp26.CESM2_WACCM.stack)

# # 8) CMCC-ESM2
# ssp26.CMCC_ESM2.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#                                                          coords=Coordinates[,c(1,2)],
#                                                          time.period='ssp26',
#                                                          time1='2015-01-01',
#                                                          time2='2100-12-31')
# ssp26.CMCC_ESM2.stack$Scenario <- 'ssp26'
# 
# save.image('Global/RData/ssp26_summary_CMCC_ESM2.RData')
# rm(ssp26.CMCC_ESM2.stack)

# # 9) GFDL-ESM4.gr
# ssp26.GFDL_ESM4_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#                                                             coords=Coordinates[,c(1,2)],
#                                                             time.period='ssp26',
#                                                             time1='2015-01-01',
#                                                             time2='2100-12-31')
# ssp26.GFDL_ESM4_gr.stack$Scenario <- 'ssp26'
# 
# save.image('Global/RData/ssp26_summary_GFDL_ESM4_gr.RData')
# rm(ssp26.GFDL_ESM4_gr.stack)

# # 10) IPSL-CM6A-LR
# ssp26.IPSL_CM6A_LR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#                                                             coords=Coordinates[,c(1,2)],
#                                                             time.period='ssp26',
#                                                             time1='2015-01-01',
#                                                             time2='2100-12-31')
# ssp26.IPSL_CM6A_LR.stack$Scenario <- 'ssp26'
# 
# save.image('Global/RData/ssp26_summary_IPSL_CM6A_LR.RData')
# rm(ssp26.IPSL_CM6A_LR.stack)

# # 11) MPI-ESM1-2-HR
# ssp26.MPI_ESM1_2_HR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#                                                              coords=Coordinates[,c(1,2)],
#                                                              time.period='ssp26',
#                                                              time1='2015-01-01',
#                                                              time2='2100-12-31')
# ssp26.MPI_ESM1_2_HR.stack$Scenario <- 'ssp26'
# 
# save.image('Global/RData/ssp26_summary_MPI_ESM1_2_HR.RData')
# rm(ssp26.MPI_ESM1_2_HR.stack)

# # 12) NorESM2-LM.gr
# ssp26.NorESM2_LM_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#                                                              coords=Coordinates[,c(1,2)],
#                                                              time.period='ssp26',
#                                                              time1='2015-01-01',
#                                                              time2='2100-12-31')
# ssp26.NorESM2_LM_gr.stack$Scenario <- 'ssp26'
# 
# save.image('Global/RData/ssp26_summary_NorESM2_LM_gr.RData')
# rm(ssp26.NorESM2_LM_gr.stack)

# ### SSP45
# # 13) CESM2-WACCM
# ssp45.CESM2_WACCM.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#                                                            coords=Coordinates[,c(1,2)],
#                                                            time.period='ssp45',
#                                                            time1='2015-01-01',
#                                                            time2='2100-12-31')
# ssp45.CESM2_WACCM.stack$Scenario <- 'ssp45'
# 
# save.image('Global/RData/ssp45_summary_CESM2_WACCM.RData')
# rm(ssp45.CESM2_WACCM.stack)

# # 14) CMCC-ESM2
# ssp45.CMCC_ESM2.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#                                                          coords=Coordinates[,c(1,2)],
#                                                          time.period='ssp45',
#                                                          time1='2015-01-01',
#                                                          time2='2100-12-31')
# ssp45.CMCC_ESM2.stack$Scenario <- 'ssp45'
# 
# save.image('Global/RData/ssp45_summary_CMCC_ESM2.RData')
# rm(ssp45.CMCC_ESM2.stack)

# # 15) GFDL-ESM4.gr
# ssp45.GFDL_ESM4_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#                                                             coords=Coordinates[,c(1,2)],
#                                                             time.period='ssp45',
#                                                             time1='2015-01-01',
#                                                             time2='2100-12-31')
# ssp45.GFDL_ESM4_gr.stack$Scenario <- 'ssp45'
# 
# save.image('Global/RData/ssp45_summary_GFDL_ESM4_gr.RData')
# rm(ssp45.GFDL_ESM4_gr.stack)

# # 16) IPSL-CM6A-LR
# ssp45.IPSL_CM6A_LR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
#                                                             coords=Coordinates[,c(1,2)],
#                                                             time.period='ssp45',
#                                                             time1='2015-01-01',
#                                                             time2='2100-12-31')
# ssp45.IPSL_CM6A_LR.stack$Scenario <- 'ssp45'
# 
# save.image('Global/RData/ssp45_summary_IPSL_CM6A_LR.RData')
# rm(ssp45.IPSL_CM6A_LR.stack)

# # 17) MPI-ESM1-2-HR
# ssp45.MPI_ESM1_2_HR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#                                                              coords=Coordinates[,c(1,2)],
#                                                              time.period='ssp45',
#                                                              time1='2015-01-01',
#                                                              time2='2100-12-31')
# ssp45.MPI_ESM1_2_HR.stack$Scenario <- 'ssp45'
# 
# save.image('Global/RData/ssp45_summary_MPI_ESM1_2_HR.RData')
# rm(ssp45.MPI_ESM1_2_HR.stack)

# # 18) NorESM2-LM.gr
# ssp45.NorESM2_LM_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#                                                              coords=Coordinates[,c(1,2)],
#                                                              time.period='ssp45',
#                                                              time1='2015-01-01',
#                                                              time2='2100-12-31')
# ssp45.NorESM2_LM_gr.stack$Scenario <- 'ssp45'
# 
# save.image('Global/RData/ssp45_summary_NorESM2_LM_gr.RData')
# rm(ssp45.NorESM2_LM_gr.stack)

 
# ### SSP85
# # 19) CESM2-WACCM
# ssp85.CESM2_WACCM.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CESM2-WACCM',
#                                                       coords=Coordinates[,c(1,2)],
#                                                       time.period='ssp85',
#                                                       time1='2015-01-01',
#                                                       time2='2100-12-31')
# ssp85.CESM2_WACCM.stack$Scenario <- 'ssp85'
# 
# save.image('Global/RData/ssp85_summary_CESM2_WACCM.RData')
# rm(ssp85.CESM2_WACCM.stack)

# # 20) CMCC-ESM2
# ssp85.CMCC_ESM2.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/CMCC-ESM2',
#                                                     coords=Coordinates[,c(1,2)],
#                                                     time.period='ssp85',
#                                                     time1='2015-01-01',
#                                                     time2='2100-12-31')
# ssp85.CMCC_ESM2.stack$Scenario <- 'ssp85'
# 
# save.image('Global/RData/ssp85_summary_CMCC_ESM2.RData')
# rm(ssp85.CMCC_ESM2.stack)

# # 21) GFDL-ESM4.gr
# ssp85.GFDL_ESM4_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/GFDL-ESM4.gr',
#                                                        coords=Coordinates[,c(1,2)],
#                                                        time.period='ssp85',
#                                                        time1='2015-01-01',
#                                                        time2='2100-12-31')
# ssp85.GFDL_ESM4_gr.stack$Scenario <- 'ssp85'
# 
# save.image('Global/RData/ssp85_summary_GFDL_ESM4_gr.RData')
# rm(ssp85.GFDL_ESM4_gr.stack)
# 
# 22) IPSL-CM6A-LR
ssp85.IPSL_CM6A_LR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/IPSL-CM6A-LR',
                                                       coords=Coordinates[,c(1,2)],
                                                       time.period='ssp85',
                                                       time1='2015-01-01',
                                                       time2='2100-12-31')
ssp85.IPSL_CM6A_LR.stack$Scenario <- 'ssp85'

save.image('Global/RData/ssp85_summary_IPSL_CM6A_LR.RData')
rm(ssp85.IPSL_CM6A_LR.stack)
# 
# # 23) MPI-ESM1-2-HR
# ssp85.MPI_ESM1_2_HR.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/MPI-ESM1-2-HR',
#                                                         coords=Coordinates[,c(1,2)],
#                                                         time.period='ssp85',
#                                                         time1='2015-01-01',
#                                                         time2='2100-12-31')
# ssp85.MPI_ESM1_2_HR.stack$Scenario <- 'ssp85'
# 
# save.image('Global/RData/ssp85_summary_MPI_ESM1_2_HR.RData')
# rm(ssp85.MPI_ESM1_2_HR.stack)
# 
# # 24) NorESM2-LM.gr
# ssp85.NorESM2_LM_gr.stack <- env.extract.bulk.ecoregion(path='/media/sun/Data1/CMIP6/NorESM2-LM.gr',
#                                                         coords=Coordinates[,c(1,2)],
#                                                         time.period='ssp85',
#                                                         time1='2015-01-01',
#                                                         time2='2100-12-31')
# ssp85.NorESM2_LM_gr.stack$Scenario <- 'ssp85'
# 
# save.image('Global/RData/ssp85_summary_NorESM2_LM_gr.RData')
# rm(ssp85.NorESM2_LM_gr.stack)
# 
