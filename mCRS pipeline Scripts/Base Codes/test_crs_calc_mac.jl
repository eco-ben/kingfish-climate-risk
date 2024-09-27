include("kingfish_crs_functions.jl")

historical_stack_ensemble = RData.load("Required Files/historical_stack_ensemble.RData")["historical.stack.ensemble"]
ssp26_stack_ensemble = RData.load("Required Files/ssp26_stack_ensemble.RData")["ssp26.stack.ensemble"]
ssp45_stack_ensemble = RData.load("Required Files/ssp45_stack_ensemble.RData")["ssp45.stack.ensemble"]
ssp85_stack_ensemble = RData.load("Required Files/ssp85_stack_ensemble.RData")["ssp85.stack.ensemble"]
spat_pcs = RData.load("Required Files/INLA_Res_Combined.RData")["Full"]

ENV = vcat(historical_stack_ensemble, ssp26_stack_ensemble, ssp45_stack_ensemble, ssp85_stack_ensemble)
historical_stack_ensemble, ssp26_stack_ensemble, ssp45_stack_ensemble, ssp85_stack_ensemble = nothing, nothing, nothing, nothing

Populations = CSV.read("Required Files/Global_Coordinates.csv", DataFrame)
select!(Populations, ["Extract.ID", "PopID"])
rename!(Populations, "Extract.ID" => :Extract_ID)
Populations.Extract_ID = string.(Populations.Extract_ID)

unwanted_IDs = ["103", "204", "200", "111", "112","34", "52"]

#
# env = ENV[ENV.Scenario .∈ [["historical", "ssp26"]], :]
# dist=Populations
# scenario="ssp26"
# env_startcol=9
# env_endcol=12
# refyear_start="1850-01-01"
# refyear_end="2000-12-31"
# Depth = "10"
#



CRS_computation_depth_together(
    Populations, 
    ENV[ENV.Scenario .∈ [["historical", "ssp26"]], :];
    scenario="ssp26",
    env_startcol=9,
    env_endcol=12,
    refyear_start="1850-01-01",
    refyear_end="2000-12-31",
    Depth="10"
)
