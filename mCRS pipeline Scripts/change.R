ENV_summary_normal = append_historical(ENV.summary) %>%
  group_by(Scenario) %>% # Can also groupb_by(Scenario, Depth, PopID, Extract.ID) as a comparison
  mutate(
    ., 
    norm_no3 = normalise(no3), 
    norm_ph = normalise(ph), 
    norm_so = normalise(so),
    norm_thetao = normalise(thetao)
  )  
  
mean_historical_conditions = ENV_summary_normal[ENV_summary_normal$Year < 2000, ] %>%
  group_by(Extract.ID, PopID, Depth) %>% # Can also group_by(PopID, Depth) as a comparison
  summarise(no3_mean = mean(norm_no3), ph_mean = mean(norm_ph), so_mean = mean(norm_so), thetao_mean = mean(norm_thetao))

change_future = ENV_summary_normal[ENV_summary_normal$Year >= 2000, ] %>%
  group_by(Extract.ID, Depth, PopID, Scenario) %>%
  left_join(., mean_historical_conditions, by=c("Extract.ID", "PopID", "Depth")) %>%
  mutate(., Delta_no3 = (norm_no3-no3_mean), Delta_ph = (norm_ph-ph_mean), Delta_so = (norm_so-so_mean), Delta_thetao = (norm_thetao-thetao_mean))

mean_change_future = change_future %>%
  group_by(Extract.ID, Depth, PopID, Scenario) %>%
  summarise(
    mean_Delta_no3 = mean(Delta_no3), 
    mean_Delta_ph = mean(Delta_ph), 
    mean_Delta_so = mean(Delta_so), 
    mean_Delta_thetao = mean(Delta_thetao)
  )

mean_change_future = tidyr::pivot_longer(
  mean_change_future, 
  c(mean_Delta_no3, mean_Delta_ph, mean_Delta_so, mean_Delta_thetao),
  names_to="Variable")

ggplot() + geom_point(
  data = mean_change_future[mean_change_future$Depth == 10 & mean_change_future$Scenario == "ssp45",], 
  aes(x = Variable, y = value, color = as.factor(Extract.ID))) +
  theme(legend.position = "none") + facet_wrap(~PopID, nrow=1)

# After normalising the data to compare changes in variables, the largest change appears to be in pH (ocean acidification)
# These environmental variable plots do not account for the what the historical ranges of the populations were. 
# So in a population with a wider range in a certain variable the same amount of environmental change [as another population] 
# may cause less increase in climate risk.

