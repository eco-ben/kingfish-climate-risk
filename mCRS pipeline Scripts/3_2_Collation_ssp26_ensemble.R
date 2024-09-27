### Working directory
setwd('/home/sun/Analyses/Kingfish/Global/')
library(data.table)
####################################################################
############ Call data and data manipulation #######################
####################################################################
####################################################################

DataCollation <- function(time.var,summary.method){  #define the DataCollation function, summary.method refers to whether enemble or each model data is needed
  RData <- list.files(path = "RData/", pattern = "*.RData")
  for(i in 1:length(RData)){
    RData[[i]] <- paste0("RData/", RData[[i]])
  }
  RData <- RData[grep(time.var,RData)] #filter out other time variant files
  RData <- RData[grep('IPSL',RData,invert = T)]
  models <- list()
  stacks <- list()
  for(i in 1:length(RData)){
    models[[i]] <- stringr::str_extract(RData[[i]], paste0('(?<=',time.var,'_summary_).*?(?=\\.RData)'))
    stacks[[i]] <- paste0(time.var,'.',models[[i]],'.stack')                              
  }
  
  new_dfs <- list()
  
  Populations <- read.csv('Data/Global_Coordinates.csv')
  Populations <- Populations[,c(9,4)]
  Populations <- data.table(Populations)
  new.stack <- data.table()
  # rm(list = c(df_name), envir= parent.frame())
  for(i in 1:length(RData)){
    load(RData[[i]])
    print(paste0('Adding ',stacks[[i]]))
    df_name <- paste0(stacks[[i]]) #name of df within RData file
    df <- get(df_name) #get only relevant df in RData file
    rm(list = c(df_name)) #remove original data to free memory
    gc()
    df <- df[!is.na(df$Value),] #remove any na values
    df <- data.table(df)
        ###categorise depth values into bins###
    df$Depth.cats <- 0 # placeholder value
    df[,Depth.cats:= fifelse(df$Depth <= 10, 10, df$Depth.cats)]
    df[,Depth.cats:= fifelse(df$Depth > 10 & df$Depth <= 50, 50, df$Depth.cats)]
    df[,Depth.cats:= fifelse(df$Depth > 50 & df$Depth <= 100, 100, df$Depth.cats)]
    df[,Depth.cats:= fifelse(df$Depth > 100 & df$Depth <= 150, 150, df$Depth.cats)]
    df <- df[, -6L]
    colnames(df)[8] <- 'Depth'
    colnames(df)[1] <- 'Extract.ID'
    gc()
    df <- df[Populations, PopID := PopID, on = 'Extract.ID']
    # df <- df[OBISpopulations, RegID := RegID, on = 'ID']
    # df <- df[df$RegID == Region,] #filter to only AUS data, not NZ
    new_name <- paste0(stacks[[i]],"_filtered") 
    new_dfs[[new_name]] <- df #add new df to list for combining
    rm(df, df_name)
    gc()
    print("combining datatables")
    new.stack <- rbindlist(list(new.stack, new_dfs[[new_name]]))
    new_dfs[[new_name]] <- NULL
    rm(new_name)
    gc()
  }
 

  #new.stack <- rbindlist(new_dfs) #combine all the new dfs together
  rm(new_dfs, RData, stacks,models, time.var, OBISpopulations)
  gc()
  
  print("data table summary")
  
  if(summary.method == 'ensemble'){
    new.stack <- new.stack[, .(Value.mean = mean(Value)), by =.(Extract.ID,Model,Scenario,Variable,PopID,Year,Month,Depth)] #this may be unnecessary if the following line gives the same result anyway
    gc()
    #model ensemble mean
    new.stack <- new.stack[, .(Value.mean = mean(Value.mean)), by=.(Extract.ID,Scenario,Variable,PopID,Year,Month,Depth)]
    gc()
    new.stack$Day <- '01'
    new.stack$Date <- with(new.stack, sprintf('%d-%02d', Year, Month))
    new.stack$Date <- paste0(new.stack$Date,'-',new.stack$Day)
    new.stack$Date <- as.Date(new.stack$Date)
    rm(env.extract.bulk.ecoregion)
    new.stack <- dcast(new.stack, Extract.ID+PopID+Scenario+Year+Month+Day+Date+Depth ~ Variable, value.var = c('Value.mean'))
  } else{
    new.stack <- new.stack[, .(Value.mean = mean(Value)), by =.(Extract.ID,Model,Scenario,Variable,PopID,Year,Month,Depth)] #this may be unnecessary if the following line gives the same result anyway
    gc()
    new.stack$Day <- '01'
    new.stack$Date <- with(new.stack, sprintf('%d-%02d', Year, Month))
    new.stack$Date <- paste0(new.stack$Date,'-',new.stack$Day)
    new.stack$Date <- as.Date(new.stack$Date)
    rm(env.extract.bulk.ecoregion)
    new.stack <- dcast(new.stack, Model+Extract.ID+PopID+Scenario+Year+Month+Day+Date+Depth ~ Variable, value.var = c('Value.mean'))
  }
 return(new.stack)
}

ssp26.stack.ensemble <- DataCollation(time.var = 'ssp26', summary.method = 'ensemble')

save.image("RData/Compiled/ssp26_stack_ensemble.RData")
