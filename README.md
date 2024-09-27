# kingfish-climate
Honours project assessing the extent of climate risk for Seriola lalandi

CMIP6 data is first downloaded and processed, occurrence data are then downloaded from OBIS and GBIF. 
mCRS pipeline scripts are then run on the processed CMIP6 data and occurrence data. This produces the raw mCRS and DeltamCRS results. Analyses script then produces statistical and visual analyses. 
For more details see 'Lab Notes' folder.

For analyses labelled SpL, the PopID variable has been replaced by "Seriola lalandi" for every ID because these analyses use the species level as a single population. The original populations are then attached to the IDs again for further analysis. 

# Data Requirements
Data files required to run scripts 4_CRS_calculation onwards, including statistical analysis are held on OneDrive to avoid exceeding GitHub size limits. Please contact `grierben777@gmail.com` for access to the current and updated data files required to run scripts. 
Datafiles required to run climate processing scripts and mCRS processing scripts (1-3 scripts) can be accessed through various websites. See `Lab Notes` folder for more details on data locations. 
