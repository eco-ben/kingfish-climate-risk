library(ggnewscale)
library(cowplot)

preprocess_results_df <- function(x, coordinates, level) {
    # Calculate annual average of CRS and DeltaCRS
    x_dt <- data.table(x)
    names(x_dt)[names(x_dt) == "Extract_ID"] <- "Extract.ID"
    x_dt <- x_dt[, .(DeltaCRS = mean(DeltaCRS), CRS = mean(CRS)), by = .(PopID, Depth, Extract.ID, Year, Scenario)]
    x_dt <- as.data.frame(x_dt)

    # Attach coordinates and analysis levels
    x_dt <- dplyr::left_join(x_dt, coordinates[, c("Extract.ID", "Alt.Longitude", "Alt.Latitude", "longitude_360")], by = c("Extract.ID"))
    x_dt$Distance_from_equator <- ifelse(x_dt$Alt.Latitude < 0, -x_dt$Alt.Latitude, x_dt$Alt.Latitude)
    x_dt$Depth <- as.factor(x_dt$Depth)
    x_dt$PopID <- as.factor(x_dt$PopID)
    x_dt$Level <- level

    # Add Population label instead of full ID
    x_dt$PopLab <- as.character(x_dt$PopID)
    x_dt$PopLab <- ifelse(x_dt$PopLab == "NE_Pacific", "NEP", x_dt$PopLab)
    x_dt$PopLab <- ifelse(x_dt$PopLab == "NW_Pacific", "NWP", x_dt$PopLab)
    x_dt$PopLab <- ifelse(x_dt$PopLab == "South_Africa", "SA", x_dt$PopLab)
    x_dt$PopLab <- ifelse(x_dt$PopLab == "South_Pacific", "SIP", x_dt$PopLab)

    return(x_dt)
}

PopID_add_PopLab <- function(x) {
    x$PopLab <- as.character(x$PopLab)
    x$PopLab <- ifelse(x$PopLab == "NE_Pacific", "NEP", x$PopLab)
    x$PopLab <- ifelse(x$PopLab == "NW_Pacific", "NWP", x$PopLab)
    x$PopLab <- ifelse(x$PopLab == "South_Africa", "SA", x$PopLab)
    x$PopLab <- ifelse(x$PopLab == "South_Pacific", "SIP", x$PopLab)

    return(x)
}

mean_Delta_temperature <- function(x) {
    # Calculate the mean historical temperature each population and depth is exposed to.
    refmeanTemp <- x[x$Scenario == "historical", ] %>%
        group_by(PopID, Depth) %>%
        summarise(TempMean = mean(thetao.mean))
    temp_and_ref <- left_join(x, refmeanTemp, by = c("PopID", "Depth"))

    # Calculate the change in temperature from historical mean
    temp_and_ref <- temp_and_ref %>%
        group_by(PopID, Depth, Year, Scenario) %>%
        reframe(DeltaTemp = (thetao.mean - TempMean))

    # Attach historical data to each scenario and relabel
    ssp26 <- temp_and_ref[temp_and_ref$Scenario %in% c("historical", "ssp26"), ]
    ssp45 <- temp_and_ref[temp_and_ref$Scenario %in% c("historical", "ssp45"), ]
    ssp85 <- temp_and_ref[temp_and_ref$Scenario %in% c("historical", "ssp85"), ]
    ssp26[ssp26$Scenario == "historical", ]$Scenario <- "ssp26"
    ssp45[ssp45$Scenario == "historical", ]$Scenario <- "ssp45"
    ssp85[ssp85$Scenario == "historical", ]$Scenario <- "ssp85"

    # Combine all scenarios' data
    all_scens <- rbind(ssp26, ssp45, ssp85)
    all_scens <- all_scens[!all_scens$Year <= 2000, ]
    all_scens$PopID <- as.factor(all_scens$PopID)
    all_scens$Depth <- as.factor(all_scens$Depth)
    all_scens$Scenario <- as.factor(all_scens$Scenario)

    return(all_scens)
}

pca_calculation <- function(dist, env, scenario = "ssp45", env_startcol = 9, env_endcol = 12, refyear_start = "1850-01-01", refyear_end = "2000-01-01", method = "") {
    all_pca <- data.frame()
    all_centroids <- data.frame()
    all_polys <- data.frame()

    env <- as.data.frame(env)
    env <- env[!env$Extract.ID %in% c("103", "204", "200", "111", "112", "34", "52"), ] # ID locations removed as of thesis writing (not including NWP loc)
    env$Scenario <- ifelse(env$Date >= "2001-01-01" & env$Date <= "2014-12-31", scenario, env$Scenario)
    env <- left_join(env, coordinates[, c(9, 1, 2)], by = c("Extract.ID"))
    env.prcomp <- env[
        env$Date >= refyear_start &
            env$Date <= refyear_end,
        c(env_startcol:env_endcol)
    ] # environmental conditions for the reference period
    env.prcomp <- prcomp(env.prcomp, scale = T) # summary prcomp - 71+% - Ensemble proportion of variance explained by the first 2 PC axes
    env.prcomp.data <- as.data.frame(predict(
        object = env.prcomp,
        newdata = env[, c(env_startcol:env_endcol)]
    )) # rotated environmental data
    env.prcomp.data <- env.prcomp.data[, c(1:2)] # only use PC1 and PC2
    env.prcomp.data$Extract.ID <- env$Extract.ID
    env.prcomp.data$Date <- env$Date
    env.prcomp.data$Scenario <- unlist(env$Scenario)
    env.prcomp.data$Month <- as.numeric(format(env.prcomp.data$Date, "%m"))
    env.prcomp.data$Depth <- unlist(env$Depth)
    env.prcomp.data <- env.prcomp.data[, c(3:5, 7, 1, 2)]
    env.prcomp.data <- env.prcomp.data[complete.cases(env.prcomp.data), ] # remove invomplete rows
    sp <- unique(dist$PopID)
    all_pca <- env.prcomp.data
    all_pca <- left_join(all_pca, dist, by = "Extract.ID")
    polygons <- vector("list", length = length(sp))
    pc.centr <- vector("list", length = length(sp))
    for (m in 1:length(sp)) {
        tmp <- unique(dist[dist$PopID == sp[[m]], ]) # IDs for relevant population
        tmp.data <- sqldf::sqldf("SELECT * from tmp LEFT JOIN
                               'env.prcomp.data' USING (`Extract.ID`)") # join tmp to PC results -> tmp.data
        tmp.data <- tmp.data[complete.cases(tmp.data), ] # remove incomplete results
        env.coords <- droplevels(tmp.data[tmp.data$Date >= refyear_start & tmp.data$Date <= refyear_end, ]) # include only the reference period and relevant population
        # env.coords <- env.coords[chull(env.coords[,5], env.coords[,6]),] #create convex hull around condition points

        env.poly <- distfree.cr(env.coords[, c("PC1", "PC2")], alpha = 0.05, draw = F)
        env.poly.pts <- data.frame(env.poly$polygon)
        FL <- env.poly.pts[1, ]
        env.poly.pts <- rbind(env.poly.pts, FL)
        env.poly <- sp::SpatialPolygons(list(
            sp::Polygons(
                list(sp::Polygon(env.poly.pts)), # the polygon used is the convex hull of environmental points
                "historical"
            )
        ))
        env.polygon_sf <- sf::st_as_sf(env.poly, coords = env.polygon@polygons[[1]]@Polygons[[1]]@coords)
        env.polygon_sf$PopID <- unlist(sp[m])

        if (method == "surface") {
            spat.pcs <- Full[Full$PopID == sp[m] & Full$Date <= refyear_end & Full$Depth == unique(env$Depth), ]
        } else {
            spat.pcs <- Full[Full$PopID == sp[m] & Full$Date <= refyear_end, ]
        }

        spat.pcs <- spat.pcs[complete.cases(spat.pcs), ]
        env.centroid <- c(mean(spat.pcs$pred.PC1), mean(spat.pcs$pred.PC2))
        env.centroid <- data.frame(env.centroid[1], env.centroid[2])
        env.centroid$PopID <- sp[m]
        all_polys <- rbind(all_polys, env.polygon_sf)
        all_centroids <- rbind(all_centroids, env.centroid)
    }

    inset <- ggplot() +
        geom_segment(
            data = as.data.frame(env.prcomp$rotation),
            aes(x = 0, y = 0, xend = (PC1), yend = (PC2)),
            arrow = arrow(length = unit(0.01, "cm")), color = "coral1"
        ) +
        annotate(
            geom = "text",
            x = (as.data.frame(env.prcomp$rotation)[1, 1]),
            y = (as.data.frame(env.prcomp$rotation)[1, 2]),
            size = 4,
            label = paste0((row.names(env.prcomp$rotation)[1]))
        ) +
        annotate(
            geom = "text",
            x = (as.data.frame(env.prcomp$rotation)[2, 1]),
            y = (as.data.frame(env.prcomp$rotation)[2, 2]),
            size = 4,
            label = paste0((row.names(env.prcomp$rotation)[2]))
        ) +
        annotate(
            geom = "text",
            x = (as.data.frame(env.prcomp$rotation)[3, 1]),
            y = (as.data.frame(env.prcomp$rotation)[3, 2]),
            size = 4,
            label = paste0((row.names(env.prcomp$rotation)[3]))
        ) +
        annotate(
            geom = "text",
            x = (as.data.frame(env.prcomp$rotation)[4, 1]),
            y = (as.data.frame(env.prcomp$rotation)[4, 2]),
            size = 4,
            label = paste0((row.names(env.prcomp$rotation)[4]))
        ) +
        theme_void() +
        scale_x_continuous(limits = c(
            (min(as.data.frame(env.prcomp$rotation)$PC1) - 1.5),
            (max(as.data.frame(env.prcomp$rotation)$PC2) + 1.5)
        )) +
        scale_y_continuous(limits = c(
            (min(as.data.frame(env.prcomp$rotation)$PC2) - 1.5),
            (max(as.data.frame(env.prcomp$rotation)$PC2) + 1.5)
        ))

    all_pca$Year <- as.numeric(format(as.Date(all_pca$Date, "%Y-m-%d"), "%Y"))
    mean_pca <- all_pca %>%
        group_by(Extract.ID, Depth, PopID, Scenario) %>%
        summarise(PC1 = mean(PC1), PC2 = mean(PC2))
    mean_pca$Depth <- as.factor(mean_pca$Depth)

    pca_results <- list("mean_pca" = mean_pca, "all_centroids" = all_centroids, "all_polys" = all_polys, "inset" = inset)

    return(pca_results)
}

# PCA_plotting = function(dist, env, scenario="ssp45", env_startcol=9, env_endcol=12, refyear_start="1850-01-01", refyear_end="2000-01-01",plot_title, method=""){
#   pca_results = pca_calculation(dist, env, scenario, env_startcol, env_endcol, refyear_start, refyear_end, plot_title, method)
#   year_pca = pca_results$year_pca
#   all_centroids = pca_results$all_centroids
#   all_polys = pca_results$all_polys
#   inset = pca_results$inset
#
#   pca_plot = ggplot() +
#     #geom_point(data = year_pca,
#     #           aes(x = PC1, y = PC2, color = as.factor(Depth)),alpha = 0.4) +
#
#     geom_point(
#       data = year_pca[year_pca$PopID == "South_Pacific",],
#       aes(x = PC1, y = PC2, color= DeltaCRS),
#       shape = 1, stroke = 2, size = 2) +
#     scale_color_scico(palette = "batlow", midpoint=0) + new_scale("color") +
#     geom_point(
#       data = year_pca[year_pca$PopID == "NE_Pacific",],
#       aes(x = PC1, y = PC2, color= DeltaCRS),
#       shape = 1, stroke = 2, size = 2) +
#     scale_color_scico(palette = "batlow", midpoint=0) + new_scale("color") +
#     geom_point(
#       data = year_pca[year_pca$PopID == "NW_Pacific",],
#       aes(x = PC1, y = PC2, color= DeltaCRS),
#       shape = 1, stroke = 2, size = 2) +
#     scale_color_scico(palette = "batlow", midpoint=0) + new_scale("color") +
#     geom_point(
#       data = year_pca[year_pca$PopID == "South_Africa",],
#       aes(x = PC1, y = PC2, color= DeltaCRS),
#       shape = 1, stroke = 2, size = 2) +
#     scale_color_scico(palette = "batlow", midpoint=0) + new_scale("color") +
#
#     geom_point(data = all_centroids,aes(x = env.centroid.1., y = env.centroid.2.), color = 'black', size = 1, alpha = 0.7) +
#     geom_sf(data = all_polys, alpha =0.3) +
#     annotation_custom(grob = ggplotGrob(inset))+
#     theme_classic() +
#     facet_wrap(~PopID) + ggtitle(paste0(plot_title)) + theme(axis.text = element_text(size = 14),
#                                                                                 axis.title = element_text(size = 17, face = 'bold'),
#                                                                                 strip.text.x = element_text(size = 12))
#   return(pca_plot)
# }

DeltamCRSMap <- function(Scenario) {
    Pacific <- ggplot() +
        geom_sf(
            data = Warmlevel.data_360[Warmlevel.data_360$Scenario == Scenario, ],
            aes(color = DeltaCRS_log), shape = 1, stroke = 2, size = 2, alpha = 0.6
        ) +
        geom_sf(data = World_360) +
        geom_sf(data = Boxes, fill = NA, linewidth = 0.2) +
        # ggtitle('DeltaCRS - 10m') +
        theme_void() +
        theme(legend.position = "none") +
        theme(panel.background = element_rect(fill = "white")) +
        scale_color_scico(palette = "batlow", limits = c(minDcrs, maxDcrs)) +
        coord_sf(xlim = c(100, 300), ylim = c(50, -60)) +
        facet_wrap(~Depth)
    SA <- ggplot() +
        geom_sf(
            data = all_depths_sf[all_depths_sf$Scenario == Scenario, ],
            aes(color = DeltaCRS_log), shape = 1, stroke = 2, size = 2, alpha = 0.6
        ) +
        geom_sf(data = World) +
        theme_void() +
        theme(legend.position = "none") +
        theme(panel.background = element_rect(fill = "white")) +
        scale_color_scico(palette = "batlow", limits = c(minDcrs, maxDcrs)) +
        coord_sf(xlim = c(10, 40), ylim = c(-20, -40)) +
        facet_wrap(~Depth)
    final <- ggarrange(Pacific, SA, heights = c(3, 1), widths = c(2.5, 1), labels = c("(a)", "(b)"), vjust = 1)
    final <- ggarrange(final, leg, ncol = 1, heights = c(6, 1), widths = c(6, 1))

    return(final)
}

append_historical <- function(x) {
    # Attach historical data to each scenario and relabel
    ssp26 <- x[x$Scenario %in% c("historical", "ssp26"), ]
    ssp45 <- x[x$Scenario %in% c("historical", "ssp45"), ]
    ssp85 <- x[x$Scenario %in% c("historical", "ssp85"), ]
    ssp26[ssp26$Scenario == "historical", ]$Scenario <- "ssp26"
    ssp45[ssp45$Scenario == "historical", ]$Scenario <- "ssp45"
    ssp85[ssp85$Scenario == "historical", ]$Scenario <- "ssp85"

    # Combine all scenarios' data
    all_scens <- rbind(ssp26, ssp45, ssp85)
    all_scens$PopID <- as.factor(all_scens$PopID)
    all_scens$Depth <- as.factor(all_scens$Depth)
    all_scens$Scenario <- as.factor(all_scens$Scenario)

    return(all_scens)
}

normalise <- function(x) {
    x <- (x - min(x)) / (max(x) - min(x))

    return(x)
}

environmental_change <- function(env_variable_data, coordinates) {
    env_data_normal <- append_historical(env_variable_data) %>%
        group_by(Scenario) %>% # Can also groupb_by(Scenario, Depth, PopID, Extract.ID) as a comparison
        mutate(
            .,
            norm_no3 = normalise(no3),
            norm_ph = normalise(ph),
            norm_so = normalise(so),
            norm_thetao = normalise(thetao)
        )

    mean_historical_conditions <- env_data_normal[env_data_normal$Year < 2000, ] %>%
        group_by(Extract.ID, PopID, Depth) %>% # Can also group_by(PopID, Depth) as a comparison
        summarise(
            no3_mean = mean(norm_no3),
            ph_mean = mean(norm_ph),
            so_mean = mean(norm_so),
            thetao_mean = mean(norm_thetao)
        )

    # Extract future change values at years when warming levels are reached
    change_future <- env_data_normal[env_data_normal$Year >= 2000, ] %>%
        group_by(Extract.ID, Depth, PopID, Scenario) %>%
        left_join(., mean_historical_conditions, by = c("Extract.ID", "PopID", "Depth")) %>%
        left_join(., Warming.levels, by = c("PopID", "Depth", "Scenario")) %>%
        filter(., Year == Warming.year) %>%
        mutate(
            .,
            Delta_no3 = (norm_no3 - no3_mean),
            Delta_ph = (norm_ph - ph_mean),
            Delta_so = (norm_so - so_mean),
            Delta_thetao = (norm_thetao - thetao_mean)
        )

    mean_change_future <- change_future %>%
        group_by(Extract.ID, Depth, PopID, Scenario) %>%
        summarise(
            mean_Delta_no3 = mean(Delta_no3),
            mean_Delta_ph = mean(Delta_ph),
            mean_Delta_so = mean(Delta_so),
            mean_Delta_thetao = mean(Delta_thetao)
        )

    mean_change_future <- tidyr::pivot_longer(
        mean_change_future,
        c(mean_Delta_no3, mean_Delta_ph, mean_Delta_so, mean_Delta_thetao),
        names_to = "Variable"
    )
    mean_change_future <- left_join(
        mean_change_future,
        coordinates[, c("Extract.ID", "Alt.Longitude", "Alt.Latitude")],
        by = "Extract.ID"
    )

    return(mean_change_future)
}

horizontal_env_jitter <- function(plotting_data, population, ylims, x_var, y_var, colouring) {
    pop_data <- plotting_data[plotting_data$PopID == population, ]
    # colouring = as.factor(pop_data$Alt.Longitude + pop_data$Alt.Latitude)
    # colouring = pop_data$DeltaCRS

    plot_x <- ggplot() +
        annotate("rect",
            xmin = c(0.5, 1.5, 2.5, 3.5),
            xmax = c(1.5, 2.5, 3.5, 4.5),
            ymin = rep(ylims[1], 4),
            ymax = rep(ylims[2], 4),
            alpha = c(0.5, 0.3, 0.1, 1),
            fill = c("black", "black", "black", "white"),
            color = rep("black", 4),
            linewidth = rep(1, 4)
        ) +
        geom_jitter(
            data = pop_data,
            aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[colouring]]),
            shape = 1, size = 3, stroke = 2
        ) +
        geom_hline(yintercept = 0) +
        scale_y_continuous(limits = ylims) +
        theme_classic() +
        theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "mm")) +
        coord_flip() +
        # scale_color_viridis_c(option="turbo")
        scale_color_scico(palette = "batlow", midpoint = 0)

    return(plot_x)
}

# see this link for how to select columns in ggplot() tidy way https://stackoverflow.com/questions/13260626/selecting-data-frame-columns-to-plot-in-ggplot2

base_points <- function(plotting_data, colouring) {
    points <- ggplot() +
        geom_point(
            data = plotting_data[plotting_data$PopID == "South_Pacific", ],
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(palette = "batlow", midpoint = 0) +
        new_scale("color") +
        geom_point(
            data = plotting_data[plotting_data$PopID == "South_Africa", ],
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(palette = "batlow", midpoint = 0) +
        new_scale("color") +
        geom_point(
            data = plotting_data[plotting_data$PopID == "NW_Pacific", ],
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(palette = "batlow", midpoint = 0) +
        new_scale("color") +
        geom_point(
            data = plotting_data[plotting_data$PopID == "NE_Pacific", ],
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(palette = "batlow", midpoint = 0) +
        theme(legend.position = "none", panel.background = NULL)

    return(points)
}

env_change_map <- function(plotting_data, colouring, legend_title, title) {
    world_360 <- sf::st_as_sf(map("world2", plot = FALSE, fill = TRUE))
    base_map <- geom_sf(data = world_360) + coord_sf(xlim = c(0, 275), ylim = c(-52, 52))

    # Creating legend for plot
    legend <- ggplot() +
        geom_point(
            data = plotting_data,
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(
            palette = "batlow",
            midpoint = 0
            # )
        ) +
        guides(color = guide_colorbar(
            barwidth = unit(40, "mm"), # bar width
            barheight = unit(5, "mm"), # bar height
            # frame.colour = "black",  # bar frame color
            # frame.linewidth = 1,  # bar frame width
            # frame.linetype = "solid",  # bar frame linetype
            ticks = T, # show the ticks on the bar
            ticks.linewidth = 1, # tick width
            ticks.colour = "black", # tick color
            draw.ulim = T, # show the tick at the upper limit
            draw.llim = T # show the tick at the lower limit
        )) +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
        ) +
        labs(color = legend_title)
    legend <- ggpubr::get_legend(legend)

    # Creating horizontal jitter plots
    S_P_grob <- horizontal_env_jitter(plotting_data, "South_Pacific", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    S_P_grob_coords <- c(170, 260, -80, -50) # These coordintates are in order: xmin, xmax, ymin, ymax

    NEP_grob <- horizontal_env_jitter(plotting_data, "NE_Pacific", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    NEP_grob_coords <- c(150, 230, 0, 30)

    NWP_grob <- horizontal_env_jitter(plotting_data, "NW_Pacific", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    NWP_grob_coords <- c(40, 110, 25, 55)

    SA_grob <- horizontal_env_jitter(plotting_data, "South_Africa", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    SA_grob_coords <- c(0, 70, -15, 15)

    # Assembling final map
    final_map <- base_points(plotting_data, colouring) +
        geom_sf(data = world_360) +
        coord_sf(xlim = c(0, 275), ylim = c(-52, 52)) +
        ggtitle(title) +
        annotation_custom(grob = ggplotGrob(S_P_grob), xmin = S_P_grob_coords[1], xmax = S_P_grob_coords[2], ymin = S_P_grob_coords[3], ymax = S_P_grob_coords[4]) +
        annotation_custom(grob = ggplotGrob(NEP_grob), xmin = NEP_grob_coords[1], xmax = NEP_grob_coords[2], ymin = NEP_grob_coords[3], ymax = NEP_grob_coords[4]) +
        annotation_custom(grob = ggplotGrob(NWP_grob), xmin = NWP_grob_coords[1], xmax = NWP_grob_coords[2], ymin = NWP_grob_coords[3], ymax = NWP_grob_coords[4]) +
        annotation_custom(grob = ggplotGrob(SA_grob), xmin = SA_grob_coords[1], xmax = SA_grob_coords[2], ymin = SA_grob_coords[3], ymax = SA_grob_coords[4])

    final <- ggarrange(final_map, legend, ncol = 1, heights = c(6, 1), widths = c(6, 1))

    return(final)
}

pca_plot <- function(pca_results, population, centroids, polys, inset, colouring) {
    pca_plot_pop <- ggplot() +
        # plot pca results
        geom_point(
            data = pca_results[pca_results$PopID == population & pca_results$Scenario == "historical", ],
            aes(x = PC1, y = PC2, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        geom_point(
            data = pca_results[pca_results$PopID == population & pca_results$Scenario != "historical", ],
            aes(x = PC1, y = PC2, color = .data[[colouring]]),
            shape = 2, stroke = 2, size = 2
        ) +
        scale_color_scico(palette = "batlow", midpoint = 0) +
        theme_classic() +
        theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "mm")) +

        # plot centroids
        geom_point(data = centroids[centroids$PopID == population, ], aes(x = env.centroid.1., y = env.centroid.2.), color = "black", size = 2) +
        # xlim(range(pca_results[pca_results$PopID == population,]$PC1)) +
        # plot polygon outlines
        geom_sf(data = polys[polys$PopID == population, ], alpha = 0.3) +

        # plot inset diagram
        annotation_custom(grob = ggplotGrob(inset))

    return(pca_plot_pop)
}

# pca_plot <- function(pca_results, centroids, polys, inset, colouring) {
#     pca_plot_x <- ggplot() +
#         # plot pca results
#         geom_point(
#             data = pca_results[pca_results$Scenario == "historical", ],
#             aes(x = PC1, y = PC2, color = .data[[colouring]]),
#             shape = 1, stroke = 2, size = 2
#         ) +
#         theme_classic() +
#         theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +

#         # plot centroids
#         geom_point(data = centroids, aes(x = env.centroid.1., y = env.centroid.2.), color = "black", size = 2) +
#         # xlim(range(pca_results[pca_results$PopID == population,]$PC1)) +
#         # plot polygon outlines
#         geom_sf(data = polys, alpha = 0.3) +

#         # plot inset diagram
#         annotation_custom(grob = ggplotGrob(inset)) +
#         facet_wrap(~PopID)

#     return(pca_plot_x)
# }

results_and_pca_plot <- function(coordinates,
                                 environmental_data,
                                 crs_results,
                                 depth,
                                 scenario,
                                 colouring,
                                 legend_title,
                                 title) {
    # getting pca results and plot ready
    env <- environmental_data
    pca_results <- pca_calculation(coordinates, env)
    mean_pca <- pca_results$mean_pca
    centroids <- pca_results$all_centroids
    polys <- pca_results$all_polys
    inset <- pca_results$inset

    crs_results <- crs_results[crs_results$Depth == depth, ]

    mean_pca <- left_join(mean_pca, crs_results[, c("Extract.ID", "PopID", "DeltaCRS", "CRS")], by = c("Extract.ID", "PopID"))
    mean_pca <- mean_pca[mean_pca$Depth == depth, ]

    coordinates_crs <- left_join(coordinates, crs_results[, c("Extract.ID", "PopID", "DeltaCRS", "CRS")], by = c("Extract.ID", "PopID"))
    coordinates_crs <- rename(coordinates_crs, Longitude_360 = longitude_360)

    ##
    world_360 <- sf::st_as_sf(map("world2", plot = FALSE, fill = TRUE))
    base_map <- geom_sf(data = world_360) + coord_sf(xlim = c(0, 275), ylim = c(-52, 52))


    # Creating legend for plot
    legend <- ggplot() +
        geom_point(
            data = coordinates_crs,
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(
            palette = "batlow",
            midpoint = 0
            # )
        ) +
        guides(color = guide_colorbar(
            barwidth = unit(40, "mm"), # bar width
            barheight = unit(5, "mm"), # bar height
            # frame.colour = "black",  # bar frame color
            # frame.linewidth = 1,  # bar frame width
            # frame.linetype = "solid",  # bar frame linetype
            ticks = T, # show the ticks on the bar
            ticks.linewidth = 1, # tick width
            ticks.colour = "black", # tick color
            draw.ulim = T, # show the tick at the upper limit
            draw.llim = T # show the tick at the lower limit
        )) +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
        ) +
        labs(color = legend_title)
    legend <- ggpubr::get_legend(legend)

    pca_plot_legend <- ggplot() +
        geom_point(aes(x = 1, y = 1, shape = "Historical")) +
        geom_point(aes(x = 1, y = 2, shape = scenario)) +
        scale_shape_manual(
            name = "Time Period",
            breaks = c("Historical", scenario),
            values = c(1, 2)
        )
    pca_plot_legend <- ggpubr::get_legend(pca_plot_legend)

    # Creating horizontal jitter plots
    S_P_grob <- pca_plot(mean_pca, "South_Pacific", centroids, polys, inset, colouring)
    # S_P_grob_coords = c(170, 260, -80, -50) # These coordintates are in order: xmin, xmax, ymin, ymax

    NEP_grob <- pca_plot(mean_pca, "NE_Pacific", centroids, polys, inset, colouring)
    # NEP_grob_coords = c(150, 230, 0, 30)

    NWP_grob <- pca_plot(mean_pca, "NW_Pacific", centroids, polys, inset, colouring)
    # NWP_grob_coords = c(50, 60, 45, 55)

    SA_grob <- pca_plot(mean_pca, "South_Africa", centroids, polys, inset, colouring)
    # SA_grob_coords = c(0, 70, -15, 15)

    # Assembling final map
    final_map <- base_points(coordinates_crs, colouring) +
        geom_sf(data = world_360) +
        coord_sf(xlim = c(0, 275), ylim = c(-52, 52))
    # annotation_custom(grob = ggplotGrob(S_P_grob), xmin=S_P_grob_coords[1], xmax=S_P_grob_coords[2], ymin=S_P_grob_coords[3], ymax=S_P_grob_coords[4]) +
    # annotation_custom(grob = ggplotGrob(NEP_grob), xmin=NEP_grob_coords[1], xmax=NEP_grob_coords[2], ymin=NEP_grob_coords[3], ymax=NEP_grob_coords[4]) +
    # annotation_custom(grob = test, xmin=NWP_grob_coords[1], xmax=NWP_grob_coords[2], ymin=NWP_grob_coords[3], ymax=NWP_grob_coords[4]) +
    # annotation_custom(grob = ggplotGrob(SA_grob), xmin=SA_grob_coords[1], xmax=SA_grob_coords[2], ymin=SA_grob_coords[3], ymax=SA_grob_coords[4])
    final_map <- ggdraw(final_map) +
        draw_plot(S_P_grob, x = 0.675, y = 0.05, width = 0.3, height = 0.275) +
        draw_plot(NEP_grob, x = 0.55, y = 0.5, width = 0.3, height = 0.25) +
        draw_plot(NWP_grob, x = 0.215, y = 0.6, width = 0.3, height = 0.3) +
        draw_plot(SA_grob, x = 0.05, y = 0.41, width = 0.3, height = 0.3) +
        draw_label(title, x = 0.13, y = 0.96)

    final <- ggarrange(final_map, legend, pca_plot_legend, ncol = 1, heights = c(6, 1, 1), widths = c(6, 1, 1))

    return(final)
}

env_change_map <- function(plotting_data, colouring, legend_title, title) {
    world_360 <- sf::st_as_sf(map("world2", plot = FALSE, fill = TRUE))
    base_map <- geom_sf(data = world_360) + coord_sf(xlim = c(0, 275), ylim = c(-52, 52))

    # Creating legend for plot
    legend <- ggplot() +
        geom_point(
            data = plotting_data,
            aes(x = Longitude_360, y = Alt.Latitude, color = .data[[colouring]]),
            shape = 1, stroke = 2, size = 2
        ) +
        scale_color_scico(
            palette = "batlow",
            midpoint = 0
            # )
        ) +
        guides(color = guide_colorbar(
            barwidth = unit(40, "mm"), # bar width
            barheight = unit(5, "mm"), # bar height
            # frame.colour = "black",  # bar frame color
            # frame.linewidth = 1,  # bar frame width
            # frame.linetype = "solid",  # bar frame linetype
            ticks = T, # show the ticks on the bar
            ticks.linewidth = 1, # tick width
            ticks.colour = "black", # tick color
            draw.ulim = T, # show the tick at the upper limit
            draw.llim = T # show the tick at the lower limit
        )) +
        theme(
            legend.position = "bottom",
            legend.direction = "horizontal",
        ) +
        labs(color = legend_title)
    legend <- ggpubr::get_legend(legend)

    # Creating horizontal jitter plots
    S_P_grob <- horizontal_env_jitter(plotting_data, "South_Pacific", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    S_P_grob_coords <- c(170, 260, -80, -50) # These coordintates are in order: xmin, xmax, ymin, ymax

    NEP_grob <- horizontal_env_jitter(plotting_data, "NE_Pacific", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    NEP_grob_coords <- c(150, 230, 0, 30)

    NWP_grob <- horizontal_env_jitter(plotting_data, "NW_Pacific", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    NWP_grob_coords <- c(40, 110, 25, 55)

    SA_grob <- horizontal_env_jitter(plotting_data, "South_Africa", c(min(plotting_data$value) - 0.01, max(plotting_data$value) + 0.01), "Variable", "value", colouring)
    SA_grob_coords <- c(0, 70, -15, 15)

    # Assembling final map
    final_map <- base_points(plotting_data, colouring) +
        geom_sf(data = world_360) +
        coord_sf(xlim = c(0, 275), ylim = c(-52, 52)) +
        ggtitle(title) +
        annotation_custom(grob = ggplotGrob(S_P_grob), xmin = S_P_grob_coords[1], xmax = S_P_grob_coords[2], ymin = S_P_grob_coords[3], ymax = S_P_grob_coords[4]) +
        annotation_custom(grob = ggplotGrob(NEP_grob), xmin = NEP_grob_coords[1], xmax = NEP_grob_coords[2], ymin = NEP_grob_coords[3], ymax = NEP_grob_coords[4]) +
        annotation_custom(grob = ggplotGrob(NWP_grob), xmin = NWP_grob_coords[1], xmax = NWP_grob_coords[2], ymin = NWP_grob_coords[3], ymax = NWP_grob_coords[4]) +
        annotation_custom(grob = ggplotGrob(SA_grob), xmin = SA_grob_coords[1], xmax = SA_grob_coords[2], ymin = SA_grob_coords[3], ymax = SA_grob_coords[4])

    final <- ggarrange(final_map, legend, ncol = 1, heights = c(6, 1), widths = c(6, 1))

    return(final)
}
