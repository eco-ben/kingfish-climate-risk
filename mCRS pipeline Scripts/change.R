library(ggnewscale)
mean_change_future = environmental_change(ENV.summary, coordinates) # need to have Warming.levels data in env
mean_change_future$Longitude_360 = ifelse(mean_change_future$Alt.Longitude < 0, mean_change_future$Alt.Longitude + 360, mean_change_future$Alt.Longitude)
# ggplot() + geom_point(
#   data = mean_change_future[mean_change_future$PopID == "NE_Pacific" & mean_change_future$Scenario == "ssp45",], 
#   aes(x = Variable, y = value, color = as.factor(Extract.ID))) +
#   theme(legend.position = "none") + facet_wrap(~Depth, nrow=1)
# 
# data_breaks = data.frame(start = c(0.5,1.5,2.5,3.5),
#                          end = c(1.5,2.5,3.5,4.5),
#                          colors= factor(1:4),
#                          Variable = c("mean_Delta_thetao", "mean_Delta_so", "mean_Delta_ph", "mean_Delta_no3"),
#                          y = rep(max(plotting_data$value), 4))
# 
# mean_change_future$Longitude_360 = ifelse(
#   mean_change_future$Alt.Longitude < 0, 
#   mean_change_future$Alt.Longitude + 360, 
#   mean_change_future$Alt.Longitude
# )

plotting_data = mean_change_future[mean_change_future$Depth == 10 & mean_change_future$Scenario == "ssp85",]
plotting_data = left_join(plotting_data, warmlevel_data_all_depths[, c('Extract.ID', 'PopID', 'Depth', 'Scenario', 'CRS', 'DeltaCRS')], by=c('Extract.ID', 'PopID', 'Depth', 'Scenario'))








# After normalising the data to compare changes in variables, the largest change appears to be in pH (ocean acidification)
# These environmental variable plots do not account for the what the historical ranges of the populations were. 
# So in a population with a wider range in a certain variable the same amount of environmental change [as another population] 
# may cause less increase in climate risk.

