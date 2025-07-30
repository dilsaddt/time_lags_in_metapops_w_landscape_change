############################################################################
# This script is to plot
# proportion of occupied sites over time
# with both scenarios
#
# Date: 20/05/2025
# Author: Dilsad Dagtekin

############################################################################
rm(list=ls())

library(ggplot2)
library(gridExtra)
library(dplyr)

# load fixed parameters to be used in plotting, if saved
load("./path_to_output/land_restoration/fixed_params/sim_landrest_1_a500_80p_fixed_params.RData")
rm(list=setdiff(ls(), c("ntime_first", "ntime_deg", "ntime_rest", "ntime_stop"))) # keep only ntime_phases to be used in plotting to begin with
gc()

year_deg_start <- ntime_first # end of initial phase, start of degradation phase
year_rest_start <- ntime_first + ntime_deg # end of degradation phase, start of restoration phase
year_rest_end <- ntime_first + ntime_deg + ntime_rest # end of restoration phase
year_sim_end <- ntime_first + ntime_deg + ntime_rest + ntime_stop # end of simulation

### load stored calculated values to makes this faster:
load(file=("./path_to_output/results_inspection/Propoccs_wLandscapeReps.RData"))
load(file=("./path_to_output/results_inspection/PropoccsMean.RData"))

str(propocc)
str(propocc.mean)

# add the phase marks -----------------------------------------------------
# phases are:
# start of degradation
# end of degradation / start of restoration
# end of restoration

year_phase <- data.frame(NA)

year_phase <- rbind(propocc.mean[which(propocc.mean$year == year_deg_start),],
                    propocc.mean[which(propocc.mean$year == year_rest_start),],
                    propocc.mean[which(propocc.mean$year == year_rest_end),])

year_phase$phase <- NA


year_phase$phase <- ifelse(year_phase$year == year_deg_start, "start of degradation",
                           ifelse(year_phase$year == year_rest_start, "end of degradation / start of restoration", "end of restoration"))

year_phase$phase <- as.factor(year_phase$phase)
year_phase$phase <- factor(year_phase$phase, levels = c("start of degradation", "end of degradation / start of restoration", "end of restoration"))


# plots -------------------------------------------------------------------

# only take after year 50, assuming all reached stable state after then
## subset to one alpha, e, c, degradation percentage first to check

thisalpha <- "a500"
thise <- "0.2"
thisc <- "1"
thisdp <- 80


temp <- propocc[propocc$alphaCat==thisalpha & 
                  propocc$e == thise & 
                  propocc$c == thisc & 
                  propocc$deg_percent == thisdp & 
                  propocc$year>50,] 

temp2 <- propocc.mean[propocc.mean$alphaCat==thisalpha & 
                        propocc.mean$e == thise & 
                        propocc.mean$c == thisc & 
                        propocc.mean$deg_percent == thisdp & 
                        propocc.mean$year>50,] 


temp3 <- year_phase %>%
  filter(alphaCat == thisalpha,
         e == thise,
         c == thisc,
         deg_percent == thisdp,
         year_phase$scenario == "Restoration",
         year > 50)

# this plot is to check if all is right, also get the legend
p_legend <- ggplot() +
  geom_vline(xintercept = c(100, 200, 250), col="grey")+
  geom_line(data = temp2, aes(x = year, y = propocc.mean, group = scenario, col = scenario, lty = scenario), lwd = 1, alpha = 1) +
  # geom_ribbon(data = temp, aes(x = year, ymin =ribbon_low_corrected, ymax = ribbon_up, group = scenario, fill = scenario), alpha = 0.1) +
  geom_point(data=temp3, aes(x=year, y= -0.1, group = scenario, pch=phase), col ="black", fill="grey", alpha = 0.9, cex = 1.5, stroke = 0.2) + 
  
  scale_color_brewer(palette = "Dark2", direction = -1) +  # Color-blind friendly categorical palette
  scale_linetype_manual(
    values = c("Baseline degradation" = "solid", "Restoration" = "dashed"),
    breaks = c("Baseline degradation", "Restoration"),  # order in the legend
    name = "scenario") +
  # scale_x_continuous(name = "Year", breaks = c(50, seq(0, max(unique(temp$year)), by = 50)), limits = c(50, max(unique(temp$year)))) +
  scale_x_continuous(name = "Year", breaks = c(50, seq(0, 300, by = 50)), limits = c(50, 300)) +
  scale_shape_manual(values = c(
    "start of degradation" = 21,
    "end of degradation / start of restoration" = 22,
    "end of restoration" = 23,
    "minimum occupancy" = 8  # or whatever shape you used
  )) +
  scale_y_continuous(name = "Proportion of Occupied Sites", breaks = seq(0, 1, by = 0.25), limits = c(-0.1, 1)) +
  
 # ylim(c(0, 1)) +
  #ggtitle(plot_title) +
  ylab("Proportion of Occupied Sites") +
  xlab("Year") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x=element_text(size=8.5),
        axis.text.y=element_text(size=8.5),
        axis.title.x=element_text(size=10), # margin = margin(10,0,0,0)),
        axis.title.y=element_text(size=10), #margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12),
        plot.title=element_text(size=14),
        #plot.subtitle = element_text(size = 11),
        legend.position = "top",
        legend.direction = "vertical",
        legend.title =  element_text(size = 16),
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 14))

out_dir <- "./path_to_output/results_inspection/plots/propocc_timeseries_plots/"  # Change this to your desired output folder

filename_base <- paste0("traj_", thisalpha, "_e", thise, "_c", thisc, "_", thisdp, "p")

# Save as PDF
ggsave(filename = file.path(out_dir, paste0("legend_", filename_base, ".pdf")),
       plot = p_legend, width = 8, height = 6)


# plotting all variable combinations in a loop ----------------------------
# all alphas, e values, c values, and degradation percentages

alphas <- unique(propocc.mean$alphaCat)
es <- unique(propocc.mean$e)
cs <- unique(propocc.mean$c)
deg_percents <- unique(propocc.mean$deg_percent)

# Output folder
out_dir <- "./path_to_output/results_inspection/plots/propocc_timeseries_plots/"  # Change this to your desired output folder

# Loop over all combinations
for (thisalpha in alphas) {
  for (thise in es) {
    for (thisc in cs) {
      for (thisdp in deg_percents) {
        
        temp <- propocc.mean %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year > 50)
        
        
        temp3 <- year_phase %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year_phase$scenario == "Restoration",
                 year > 50)
        
        if (nrow(temp) == 0) next  # Skip if no data
        if (nrow(temp3) == 0) next  # Skip if no data
        
        temp$ribbon_up <- temp$propocc.mean + temp$sd
        temp$ribbon_low <- temp$propocc.mean - temp$sd
        temp$ribbon_low_corrected <- ifelse(temp$ribbon_low < 0, 0, temp$ribbon_low)
        
        
        plot_title <- paste0("Trajectories, ", thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp, "% degradation")
        filename_base <- paste0("traj_", thisalpha, "_e", thise, "_c", thisc, "_", thisdp, "p")
        
        p <- ggplot() +
          geom_vline(xintercept = c(100, 200, 250), col="grey")+
          geom_line(data = temp, aes(x = year, y = propocc.mean, group = scenario, col = scenario, lty = scenario), lwd = 1, alpha = 1) +
          # geom_ribbon(data = temp, aes(x = year, ymin =ribbon_low_corrected, ymax = ribbon_up, group = scenario, fill = scenario), alpha = 0.1) +
          geom_point(data=temp3, aes(x=year, y= -0.1, group = scenario, pch=phase), col ="black", fill="grey", alpha = 0.9, cex = 1.5, stroke = 0.2) + 
          
          scale_color_brewer(palette = "Dark2", direction = -1) +  # Color-blind friendly categorical palette
          scale_linetype_manual(
            values = c("Baseline degradation" = "solid", "Restoration" = "dashed"),
            breaks = c("Baseline degradation", "Restoration"),  # order in the legend
            name = "scenario") +
          # scale_x_continuous(name = "Year", breaks = c(50, seq(0, max(unique(temp$year)), by = 50)), limits = c(50, max(unique(temp$year)))) +
          scale_x_continuous(name = "Year", breaks = c(50, seq(0, 300, by = 50)), limits = c(50, 300)) +
          scale_shape_manual(values = c(
            "start of degradation" = 21,
            "end of degradation / start of restoration" = 22,
            "end of restoration" = 23,
            "minimum occupancy" = 8  # or whatever shape you used
          )) +
          scale_y_continuous(name = "Proportion of Occupied Sites", breaks = seq(0, 1, by = 0.25), limits = c(-0.1, 1)) +
          
          #ggtitle(plot_title) +
          ylab("Proportion of Occupied Sites") +
          xlab("Year") +
          theme_bw() +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.text.x=element_text(size=8.5),
                axis.text.y=element_text(size=8.5),
                axis.title.x=element_text(size=10), # margin = margin(10,0,0,0)),
                axis.title.y=element_text(size=10), #margin = margin(0,10,0,0)),
                strip.text = element_text(size = 12),
                plot.title=element_text(size=10),
                #plot.subtitle = element_text(size = 11),
                legend.position = "none",
                legend.direction = "vertical",
                legend.title =  element_text(size = 12),
                legend.text = element_text(size = 12))
        
        # Save as PDF
        ggsave(filename = file.path(out_dir, paste0(filename_base, ".pdf")),
               plot = p, width = 2.5, height = 2.5)
        
        # Save as PNG
        # ggsave(filename = file.path(out_dir, paste0(filename_base, ".png")),
        #       plot = p, width = 8, height = 6, dpi = 300)
      }
    }
  }
}

