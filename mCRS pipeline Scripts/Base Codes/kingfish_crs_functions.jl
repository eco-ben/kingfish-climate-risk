using CSV, DataFrames, ArchGDAL, MultivariateStats, RData, Dates, Statistics, RCall, Distances, FLoops

function CRS_computation_depth_together(
    dist, 
    env; 
    scenario, 
    env_startcol, 
    env_endcol, 
    refyear_start, 
    refyear_end, 
    Depth
    )

    env.Extract_ID = string.(Int64.(env.Extract_ID))
    env.Depth = string.(Int64.(env.Depth))

    sp = ["NE_Pacific", "NW_Pacific", "South_Africa", "South_Pacific"]
    env = env[(env.Extract_ID .∉ [unwanted_IDs]), :]
    env.Scenario = ifelse.(Date("2001-01-01") .<= env.Date .<= Date("2014-12-31"), scenario, env.Scenario)

    env_prcomp = env[(Date(refyear_start) .<= env.Date .<= Date(refyear_end)), env_startcol:env_endcol]

    env_prcomp_norm = Matrix(env_prcomp)

    env_train_mean = mean(env_prcomp_norm, dims=1)
    env_train_std = std(env_prcomp_norm, dims=1)

    env_prcomp_norm = (env_prcomp_norm .- env_train_mean) ./ env_train_std
    env_prcomp_norm = Matrix{Float64}(env_prcomp_norm')

    env_prcomp = fit(PCA, env_prcomp_norm; method=:svd, maxoutdim=4)

    env_res_norm = Matrix(env[:,env_startcol:env_endcol])
    env_res_norm = (env_res_norm .- env_train_mean) ./ env_train_std

    env_prcomp_data = predict(env_prcomp, env_res_norm')
    env_prcomp_data = DataFrame(env_prcomp_data', [:PC1, :PC2, :PC3, :PC4]) 
    env_prcomp_data = env_prcomp_data[:,1:2]
    env_prcomp_data.Extract_ID = env.Extract_ID
    env_prcomp_data.Depth = env.Depth
    env_prcomp_data.Date = env.Date
    env_prcomp_data.Scenario = env.Scenario

    select!(env_prcomp_data, Not(:PC1, :PC2), [:PC1, :PC2])
    env_prcomp_data.PC1 = -1 .* env_prcomp_data.PC1
    Dates_list = unique(env.Date)

    function calculate_crs_for_population(pop_name::String)::DataFrame

        tmp = dist[dist.PopID .== pop_name, :]
        tmp_data = leftjoin(env_prcomp_data, tmp; on=:Extract_ID, order=:left)
        dropmissing!(tmp_data)
    
        env_coords = tmp_data[(Date(refyear_start) .<= tmp_data.Date .<= Date(refyear_end)),:]
        sort!(env_coords, [:Depth])
        env_coords_pc1_pc2 = env_coords[:, [:PC1, :PC2]]
        env_polygon_points = R"data.frame(distfree.cr::distfree.cr($(env_coords_pc1_pc2), alpha = 0.05, draw = F)$polygon)"
        env_polygon_points = rcopy(env_polygon_points)
    
        env_polygon = ArchGDAL.createpolygon()
        coord_tuples = [(row[1], row[2]) for row in eachrow(env_polygon_points)]
        push!(coord_tuples, (env_polygon_points[1,1], env_polygon_points[1,2]))
        reverse!(coord_tuples)
        ring = ArchGDAL.unsafe_createlinearring(coord_tuples)
        ArchGDAL.addgeom!(env_polygon, ring)
    
        pop_spat_pcs = spat_pcs[(spat_pcs.PopID .== pop_name) .& (spat_pcs.Date .<= Date(refyear_end)), :]
        dropmissing!(pop_spat_pcs)
        env_centroid = (mean(pop_spat_pcs[:, "pred_PC1"]), mean(pop_spat_pcs[:, "pred_PC2"]))
    
        tmp_data = tmp_data[tmp_data.Depth .== Depth, :]
    
        tmp_data.PC1_centroid .= 0.0
        tmp_data.PC2_centroid .= 0.0
        tmp_data.PC1_outside .= 0.0
        tmp_data.PC2_outside .= 0.0
        tmp_data.PC1_edge .= 0.0
        tmp_data.PC2_edge .= 0.0
        tmp_data.DM .= 0.0
        tmp_data.DC .= 0.0
        tmp_data.DE .= 0.0
        tmp_data.CRS .= 0.0
    
        function pc_outside(x_1, x_2, centroid_1, centroid_2; pc::String)
            if pc == "PC1"
                return x_1 + (x_1 - centroid_1) / sqrt((centroid_1 - x_1)^2 + (centroid_2 - x_2)^2) * 20
            elseif pc == "PC2"
                return x_2 + (x_2 - centroid_2) / sqrt((centroid_1 - x_1)^2 + (centroid_2 - x_2)^2) * 20
            end
        end
        
        function pc_edge(centroid_1, centroid_2, pc1_outside, pc2_outside, polygon)
            outside_line = ArchGDAL.createlinestring()
            ArchGDAL.addpoint!(outside_line, centroid_1, centroid_2)
            ArchGDAL.addpoint!(outside_line, pc1_outside, pc2_outside)
                
            centre_edge_line = ArchGDAL.intersection(polygon, outside_line)
            PC1_edge = ArchGDAL.getpoint(centre_edge_line, 1)[1]
            PC2_edge = ArchGDAL.getpoint(centre_edge_line, 1)[2]
        
            return PC1_edge, PC2_edge
        end
        
        function distance_and_crs(PC1, PC2, PC1_edge, PC2_edge, centroid_1, centroid_2)
            DM = Distances.Euclidean()((PC1, PC2), (PC1_edge, PC2_edge))
            DC = Distances.Euclidean()((PC1, PC2), (centroid_1, centroid_2))
            DE = Distances.Euclidean()((PC1_edge, PC2_edge), (centroid_1, centroid_2)) 
        
            CRS = (DC / DE) - 1
        
            return DM, DC, DE, CRS
        end
        
        function calculations(centroid_1, centroid_2, data, env_polygon)
            data[:, :PC1_centroid] .= centroid_1
            data[:, :PC2_centroid] .= centroid_2
            
            data[:, :PC1_outside] .= pc_outside(data[:, :PC1][1], data[:, :PC2][1], centroid_1, centroid_2; pc = "PC1")
            data[:, :PC2_outside] .= pc_outside(data[:, :PC1][1], data[:, :PC2][1], centroid_1, centroid_2; pc = "PC2")
            
            data[:, :PC1_edge] .= pc_edge(centroid_1, centroid_2, data[:, :PC1_outside][1], data[:, :PC2_outside][1], env_polygon)[1]
            data[:, :PC2_edge] .= pc_edge(centroid_1, centroid_2, data[:, :PC1_outside][1], data[:, :PC2_outside][1], env_polygon)[2]
        
            data[:, :DM] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[1]
            data[:, :DC] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[2]
            data[:, :DE] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[3]
            data[:, :CRS] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[4]
            
            return data
        end
        
        for month in Dates_list
            println("Analysing $pop_name of $(length(sp)) populations - in Date $month")
            @floop for ID in unique(tmp_data.Extract_ID)
                tmp_data[(tmp_data.Date .== month) .& (tmp_data.Extract_ID .== ID), :] = calculations(env_centroid[1], env_centroid[2], tmp_data[(tmp_data.Date .== month) .& (tmp_data.Extract_ID .== ID), :], env_polygon)
            end
        end
    
        return tmp_data
    end

    results_NE_Pacific = calculate_crs_for_population("NE_Pacific")
    results_NW_Pacific = calculate_crs_for_population("NW_Pacific")
    results_South_Africa = calculate_crs_for_population("South_Africa")
    results_South_Pacific = calculate_crs_for_population("South_Pacific")
    
    results = vcat(results_NE_Pacific, results_NW_Pacific, results_South_Africa, results_South_Pacific)

    results_NE_Pacific = nothing
    results_NW_Pacific = nothing
    results_South_Africa = nothing
    results_South_Pacific = nothing

    CSV.write("./Results/mCRS-DeltamCRS RAW/depth-together/mCRS_raw_populations_$(scenario)_$(Depth)m.csv", results)
    results = nothing
    GC.gc()
end

function CRS_computation_depth_together_SpL(
    dist, 
    env; 
    scenario, 
    env_startcol, 
    env_endcol, 
    refyear_start, 
    refyear_end, 
    Depth
    )

    env.Extract_ID = string.(Int64.(env.Extract_ID))
    env.Depth = string.(Int64.(env.Depth))

    sp = ["NE_Pacific", "NW_Pacific", "South_Africa", "South_Pacific"]
    env = env[(env.Extract_ID .∉ [unwanted_IDs]), :]
    env.Scenario = ifelse.(Date("2001-01-01") .<= env.Date .<= Date("2014-12-31"), scenario, env.Scenario)

    env_prcomp = env[(Date(refyear_start) .<= env.Date .<= Date(refyear_end)), env_startcol:env_endcol]

    env_prcomp_norm = Matrix(env_prcomp)

    env_train_mean = mean(env_prcomp_norm, dims=1)
    env_train_std = std(env_prcomp_norm, dims=1)

    env_prcomp_norm = (env_prcomp_norm .- env_train_mean) ./ env_train_std
    env_prcomp_norm = Matrix{Float64}(env_prcomp_norm')

    env_prcomp = fit(PCA, env_prcomp_norm; method=:svd, maxoutdim=4)

    env_res_norm = Matrix(env[:,env_startcol:env_endcol])
    env_res_norm = (env_res_norm .- env_train_mean) ./ env_train_std

    env_prcomp_data = predict(env_prcomp, env_res_norm')
    env_prcomp_data = DataFrame(env_prcomp_data', [:PC1, :PC2, :PC3, :PC4]) 
    env_prcomp_data = env_prcomp_data[:,1:2]
    env_prcomp_data.Extract_ID = env.Extract_ID
    env_prcomp_data.Depth = env.Depth
    env_prcomp_data.Date = env.Date
    env_prcomp_data.Scenario = env.Scenario

    select!(env_prcomp_data, Not(:PC1, :PC2), [:PC1, :PC2])
    env_prcomp_data.PC1 = -1 .* env_prcomp_data.PC1
    Dates_list = unique(env.Date)

    function calculate_crs_for_population(pop_name::String)::DataFrame

        tmp = dist[dist.PopID .== pop_name, :]
        tmp_data = leftjoin(env_prcomp_data, tmp; on=:Extract_ID, order=:left)
        dropmissing!(tmp_data)
    
        env_coords = tmp_data[(Date(refyear_start) .<= tmp_data.Date .<= Date(refyear_end)),:]
        sort!(env_coords, [:Depth])
        env_coords_pc1_pc2 = env_coords[:, [:PC1, :PC2]]
        env_polygon_points = R"data.frame(distfree.cr::distfree.cr($(env_coords_pc1_pc2), alpha = 0.05, draw = F)$polygon)"
        env_polygon_points = rcopy(env_polygon_points)
    
        env_polygon = ArchGDAL.createpolygon()
        coord_tuples = [(row[1], row[2]) for row in eachrow(env_polygon_points)]
        push!(coord_tuples, (env_polygon_points[1,1], env_polygon_points[1,2]))
        reverse!(coord_tuples)
        ring = ArchGDAL.unsafe_createlinearring(coord_tuples)
        ArchGDAL.addgeom!(env_polygon, ring)
    
        spat_pcs.species .= "Seriola lalandi SpL"
        pop_spat_pcs = spat_pcs[(spat_pcs.species .== pop_name) .& (spat_pcs.Date .<= Date(refyear_end)), :]
        dropmissing!(pop_spat_pcs)
        env_centroid = (mean(pop_spat_pcs[:, "pred_PC1"]), mean(pop_spat_pcs[:, "pred_PC2"]))
    
        tmp_data = tmp_data[tmp_data.Depth .== Depth, :]
    
        tmp_data.PC1_centroid .= 0.0
        tmp_data.PC2_centroid .= 0.0
        tmp_data.PC1_outside .= 0.0
        tmp_data.PC2_outside .= 0.0
        tmp_data.PC1_edge .= 0.0
        tmp_data.PC2_edge .= 0.0
        tmp_data.DM .= 0.0
        tmp_data.DC .= 0.0
        tmp_data.DE .= 0.0
        tmp_data.CRS .= 0.0
    
        function pc_outside(x_1, x_2, centroid_1, centroid_2; pc::String)
            if pc == "PC1"
                return x_1 + (x_1 - centroid_1) / sqrt((centroid_1 - x_1)^2 + (centroid_2 - x_2)^2) * 20
            elseif pc == "PC2"
                return x_2 + (x_2 - centroid_2) / sqrt((centroid_1 - x_1)^2 + (centroid_2 - x_2)^2) * 20
            end
        end
        
        function pc_edge(centroid_1, centroid_2, pc1_outside, pc2_outside, polygon)
            outside_line = ArchGDAL.createlinestring()
            ArchGDAL.addpoint!(outside_line, centroid_1, centroid_2)
            ArchGDAL.addpoint!(outside_line, pc1_outside, pc2_outside)
                
            centre_edge_line = ArchGDAL.intersection(polygon, outside_line)
            PC1_edge = ArchGDAL.getpoint(centre_edge_line, 1)[1]
            PC2_edge = ArchGDAL.getpoint(centre_edge_line, 1)[2]
        
            return PC1_edge, PC2_edge
        end
        
        function distance_and_crs(PC1, PC2, PC1_edge, PC2_edge, centroid_1, centroid_2)
            DM = Distances.Euclidean()((PC1, PC2), (PC1_edge, PC2_edge))
            DC = Distances.Euclidean()((PC1, PC2), (centroid_1, centroid_2))
            DE = Distances.Euclidean()((PC1_edge, PC2_edge), (centroid_1, centroid_2)) 
        
            CRS = (DC / DE) - 1
        
            return DM, DC, DE, CRS
        end
        
        function calculations(centroid_1, centroid_2, data, env_polygon)
            data[:, :PC1_centroid] .= centroid_1
            data[:, :PC2_centroid] .= centroid_2
            
            data[:, :PC1_outside] .= pc_outside(data[:, :PC1][1], data[:, :PC2][1], centroid_1, centroid_2; pc = "PC1")
            data[:, :PC2_outside] .= pc_outside(data[:, :PC1][1], data[:, :PC2][1], centroid_1, centroid_2; pc = "PC2")
            
            data[:, :PC1_edge] .= pc_edge(centroid_1, centroid_2, data[:, :PC1_outside][1], data[:, :PC2_outside][1], env_polygon)[1]
            data[:, :PC2_edge] .= pc_edge(centroid_1, centroid_2, data[:, :PC1_outside][1], data[:, :PC2_outside][1], env_polygon)[2]
        
            data[:, :DM] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[1]
            data[:, :DC] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[2]
            data[:, :DE] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[3]
            data[:, :CRS] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[4]
            
            return data
        end
        
        for month in Dates_list
            println("Analysing $pop_name of $(length(sp)) populations - in Date $month")
            @floop for ID in unique(tmp_data.Extract_ID)
                tmp_data[(tmp_data.Date .== month) .& (tmp_data.Extract_ID .== ID), :] = calculations(env_centroid[1], env_centroid[2], tmp_data[(tmp_data.Date .== month) .& (tmp_data.Extract_ID .== ID), :], env_polygon)
            end
        end
    
        return tmp_data
    end

    results_SpL = calculate_crs_for_population("Seriola lalandi SpL")

    CSV.write("./Results/mCRS-DeltamCRS RAW/depth-together/mCRS_raw_SpL_$(scenario)_$(Depth)m.csv", results_SpL)
    results_SpL = nothing
    GC.gc()
end

function CRS_computation_surface(
    dist, 
    env; 
    scenario, 
    env_startcol, 
    env_endcol, 
    refyear_start, 
    refyear_end, 
    Depth
    )

    env.Extract_ID = string.(Int64.(env.Extract_ID))
    env.Depth = string.(Int64.(env.Depth))

    sp = ["NE_Pacific", "NW_Pacific", "South_Africa", "South_Pacific"]
    env = env[(env.Extract_ID .∉ [unwanted_IDs]), :]
    env.Scenario = ifelse.(Date("2001-01-01") .<= env.Date .<= Date("2014-12-31"), scenario, env.Scenario)

    env_prcomp = env[(Date(refyear_start) .<= env.Date .<= Date(refyear_end)), env_startcol:env_endcol]

    env_prcomp_norm = Matrix(env_prcomp)

    env_train_mean = mean(env_prcomp_norm, dims=1)
    env_train_std = std(env_prcomp_norm, dims=1)

    env_prcomp_norm = (env_prcomp_norm .- env_train_mean) ./ env_train_std
    env_prcomp_norm = Matrix{Float64}(env_prcomp_norm')

    env_prcomp = fit(PCA, env_prcomp_norm; method=:svd, maxoutdim=4)

    env_res_norm = Matrix(env[:,env_startcol:env_endcol])
    env_res_norm = (env_res_norm .- env_train_mean) ./ env_train_std

    env_prcomp_data = predict(env_prcomp, env_res_norm')
    env_prcomp_data = DataFrame(env_prcomp_data', [:PC1, :PC2, :PC3, :PC4]) 
    env_prcomp_data = env_prcomp_data[:,1:2]
    env_prcomp_data.Extract_ID = env.Extract_ID
    env_prcomp_data.Depth = env.Depth
    env_prcomp_data.Date = env.Date
    env_prcomp_data.Scenario = env.Scenario

    select!(env_prcomp_data, Not(:PC1, :PC2), [:PC1, :PC2])
    env_prcomp_data.PC1 = -1 .* env_prcomp_data.PC1
    Dates_list = unique(env.Date)

    function calculate_crs_for_population(pop_name::String)::DataFrame

        tmp = dist[dist.PopID .== pop_name, :]
        tmp_data = leftjoin(env_prcomp_data, tmp; on=:Extract_ID, order=:left)
        dropmissing!(tmp_data)
    
        env_coords = tmp_data[(Date(refyear_start) .<= tmp_data.Date .<= Date(refyear_end)),:]
        sort!(env_coords, [:Depth])
        env_coords_pc1_pc2 = env_coords[:, [:PC1, :PC2]]
        env_polygon_points = R"data.frame(distfree.cr::distfree.cr($(env_coords_pc1_pc2), alpha = 0.05, draw = F)$polygon)"
        env_polygon_points = rcopy(env_polygon_points)
    
        env_polygon = ArchGDAL.createpolygon()
        coord_tuples = [(row[1], row[2]) for row in eachrow(env_polygon_points)]
        push!(coord_tuples, (env_polygon_points[1,1], env_polygon_points[1,2]))
        reverse!(coord_tuples)
        ring = ArchGDAL.unsafe_createlinearring(coord_tuples)
        ArchGDAL.addgeom!(env_polygon, ring)
    
        pop_spat_pcs = spat_pcs[(spat_pcs.PopID .== pop_name) .& (spat_pcs.Date .<= Date(refyear_end)) .& (spat_pcs.Depth .== Depth), :]
        dropmissing!(pop_spat_pcs)
        env_centroid = (mean(pop_spat_pcs[:, "pred_PC1"]), mean(pop_spat_pcs[:, "pred_PC2"]))
    
        tmp_data = tmp_data[tmp_data.Depth .== Depth, :]
    
        tmp_data.PC1_centroid .= 0.0
        tmp_data.PC2_centroid .= 0.0
        tmp_data.PC1_outside .= 0.0
        tmp_data.PC2_outside .= 0.0
        tmp_data.PC1_edge .= 0.0
        tmp_data.PC2_edge .= 0.0
        tmp_data.DM .= 0.0
        tmp_data.DC .= 0.0
        tmp_data.DE .= 0.0
        tmp_data.CRS .= 0.0
    
        function pc_outside(x_1, x_2, centroid_1, centroid_2; pc::String)
            if pc == "PC1"
                return x_1 + (x_1 - centroid_1) / sqrt((centroid_1 - x_1)^2 + (centroid_2 - x_2)^2) * 20
            elseif pc == "PC2"
                return x_2 + (x_2 - centroid_2) / sqrt((centroid_1 - x_1)^2 + (centroid_2 - x_2)^2) * 20
            end
        end
        
        function pc_edge(centroid_1, centroid_2, pc1_outside, pc2_outside, polygon)
            outside_line = ArchGDAL.createlinestring()
            ArchGDAL.addpoint!(outside_line, centroid_1, centroid_2)
            ArchGDAL.addpoint!(outside_line, pc1_outside, pc2_outside)
                
            centre_edge_line = ArchGDAL.intersection(polygon, outside_line)
            PC1_edge = ArchGDAL.getpoint(centre_edge_line, 1)[1]
            PC2_edge = ArchGDAL.getpoint(centre_edge_line, 1)[2]
        
            return PC1_edge, PC2_edge
        end
        
        function distance_and_crs(PC1, PC2, PC1_edge, PC2_edge, centroid_1, centroid_2)
            DM = Distances.Euclidean()((PC1, PC2), (PC1_edge, PC2_edge))
            DC = Distances.Euclidean()((PC1, PC2), (centroid_1, centroid_2))
            DE = Distances.Euclidean()((PC1_edge, PC2_edge), (centroid_1, centroid_2)) 
        
            CRS = (DC / DE) - 1
        
            return DM, DC, DE, CRS
        end
        
        function calculations(centroid_1, centroid_2, data, env_polygon)
            data[:, :PC1_centroid] .= centroid_1
            data[:, :PC2_centroid] .= centroid_2
            
            data[:, :PC1_outside] .= pc_outside(data[:, :PC1][1], data[:, :PC2][1], centroid_1, centroid_2; pc = "PC1")
            data[:, :PC2_outside] .= pc_outside(data[:, :PC1][1], data[:, :PC2][1], centroid_1, centroid_2; pc = "PC2")
            
            data[:, :PC1_edge] .= pc_edge(centroid_1, centroid_2, data[:, :PC1_outside][1], data[:, :PC2_outside][1], env_polygon)[1]
            data[:, :PC2_edge] .= pc_edge(centroid_1, centroid_2, data[:, :PC1_outside][1], data[:, :PC2_outside][1], env_polygon)[2]
        
            data[:, :DM] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[1]
            data[:, :DC] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[2]
            data[:, :DE] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[3]
            data[:, :CRS] .= distance_and_crs(data[:, :PC1][1], data[:, :PC2][1], data[:, :PC1_edge][1], data[:, :PC2_edge][1], centroid_1, centroid_2)[4]
            
            return data
        end
        
        for month in Dates_list
            println("Analysing $pop_name of $(length(sp)) populations - in Date $month")
            @floop for ID in unique(tmp_data.Extract_ID)
                tmp_data[(tmp_data.Date .== month) .& (tmp_data.Extract_ID .== ID), :] = calculations(env_centroid[1], env_centroid[2], tmp_data[(tmp_data.Date .== month) .& (tmp_data.Extract_ID .== ID), :], env_polygon)
            end
        end
    
        return tmp_data
    end

    results_NE_Pacific = calculate_crs_for_population("NE_Pacific")
    results_NW_Pacific = calculate_crs_for_population("NW_Pacific")
    results_South_Africa = calculate_crs_for_population("South_Africa")
    results_South_Pacific = calculate_crs_for_population("South_Pacific")

    results = vcat(results_NE_Pacific, results_NW_Pacific, results_South_Africa, results_South_Pacific)

    results_NE_Pacific = nothing
    results_NW_Pacific = nothing
    results_South_Africa = nothing
    results_South_Pacific = nothing

    CSV.write("./Results/mCRS-DeltamCRS RAW/surface_depth/mCRS_raw_populations_$(scenario)_$(Depth)m.csv", results)
    results = nothing
    GC.gc()
end
