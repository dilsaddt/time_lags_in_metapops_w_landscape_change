############################################################################
# This script is to plot
# percentage of extinct populations over time
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
rm(list=setdiff(ls(), c("ntime_first", "ntime_deg", "ntime_rest", "ntime_stop"))) # keep only out to begin with
gc()

year_deg_start <- ntime_first # end of initial phase, start of degradation phase
year_rest_start <- ntime_first + ntime_deg # end of degradation phase, start of restoration phase
year_rest_end <- ntime_first + ntime_deg + ntime_rest # end of restoration phase
year_sim_end <- ntime_first + ntime_deg + ntime_rest + ntime_stop # end of simulation

### load stored calculated values to makes this faster:
load("./path_to_output/results_inspection/ForMetapopExtPercentage.RData")

str(out.sim)
str(cumulative.extinct.by.year)

cumulative.extinct.by.year$cumulative_extinct_perc <- (cumulative.extinct.by.year$cumulative_extinct*100)/50
# calculate mean over landscapeIDs
mean.cumulative.extinct.perc <- aggregate(cumulative.extinct.by.year$cumulative_extinct_perc, 
                                          by = list(cumulative.extinct.by.year$year, cumulative.extinct.by.year$alphaCat, 
                                                    cumulative.extinct.by.year$e, cumulative.extinct.by.year$c, 
                                                    cumulative.extinct.by.year$deg_percent, cumulative.extinct.by.year$scenario), 
                                          FUN = mean, na.rm = T)

names(mean.cumulative.extinct.perc) <- c("year", "alphaCat", "e", "c", "deg_percent", "scenario", "mean_cumulative_extinct_perc")


# add the phase marks -----------------------------------------------------
# phases are:
# start of degradation
# end of degradation / start of restoration
# end of restoration

year_phase <- data.frame(NA)

year_phase <- rbind(mean.cumulative.extinct.perc[which(mean.cumulative.extinct.perc$year == year_deg_start),],
                    mean.cumulative.extinct.perc[which(mean.cumulative.extinct.perc$year == year_rest_start),],
                    mean.cumulative.extinct.perc[which(mean.cumulative.extinct.perc$year == year_rest_end),])

year_phase$phase <- NA


year_phase$phase <- ifelse(year_phase$year == year_deg_start, "start of degradation",
                           ifelse(year_phase$year == year_rest_start, "end of degradation / start of restoration", "end of restoration"))

year_phase$phase <- as.factor(year_phase$phase)
year_phase$phase <- factor(year_phase$phase, levels = c("start of degradation", "end of degradation / start of restoration", "end of restoration"))



# plots -------------------------------------------------------------------

# only take after year 50, assuming all reached stable state after then
## subset to one alpha, e, c, degradation percentage first to check

thisalpha <- "a500"
thise <- "0.8"
thisc <- "1"
thisdp <- 80
# thisland <- 1 # if need to be checked for only one simulated landscape (out of 20)

temp <- cumulative.extinct.by.year[cumulative.extinct.by.year$alphaCat==thisalpha & 
                                     cumulative.extinct.by.year$e == thise & 
                                     cumulative.extinct.by.year$c == thisc & 
                                     cumulative.extinct.by.year$deg_percent == thisdp & 
                                     #cumulative.extinct.by.year$landscapeID == thisland &
                                     cumulative.extinct.by.year$year>50,] 

temp2 <- mean.cumulative.extinct.perc[mean.cumulative.extinct.perc$alphaCat==thisalpha & 
                                        mean.cumulative.extinct.perc$e == thise & 
                                        mean.cumulative.extinct.perc$c == thisc & 
                                        mean.cumulative.extinct.perc$deg_percent == thisdp & 
                                        mean.cumulative.extinct.perc$year>50,] 


temp3 <- year_phase %>%
  filter(alphaCat == thisalpha,
         e == thise,
         c == thisc,
         deg_percent == thisdp,
         year_phase$scenario == "Restoration",
         year > 50)

# this plot is just to check

ggplot() +
  geom_vline(xintercept = c(100, 200, 250), col="grey")+
  # geom_line(data=temp, aes(x=year,  y = cumulative_extinct_perc, group = interaction(landscapeID,scenario),lty = scenario),  col = 'gray',  alpha = 0.3, lwd = 0.5) +
  geom_line(data=temp2, aes(x=year, y = mean_cumulative_extinct_perc, group = scenario, col=scenario, lty = scenario), lwd = 1, alpha = 1) + 
  geom_point(data=temp3, aes(x=year, y= (mean_cumulative_extinct_perc-3), group = scenario, pch=phase), col ="black", fill="grey", alpha = 0.7, cex = 2, stroke = 1) + 
  
  scale_color_brewer(palette = "Dark2", direction = -1) +  # Color-blind friendly categorical palette
  scale_linetype_manual(
    values = c("Baseline degradation" = "solid", "Restoration" = "dashed"),
    breaks = c("Baseline degradation", "Restoration"),  # order in the legend
    name = "scenario") +
  # scale_x_continuous(name = "Year", breaks = c(50, seq(0, max(unique(temp$year)), by = 50)), limits = c(50, max(unique(temp$year)))) +
  scale_x_continuous(name = "Year", breaks = c(50, seq(0, 300, by = 50)), limits = c(50, 300)) +
  scale_y_continuous(name = "Metapopulation Extinction %", breaks = seq(0, 100, by = 25), limits = c(-3, 100)) +
  scale_shape_manual(values = c(
    "start of degradation" = 21,
    "end of degradation / start of restoration" = 22,
    "end of restoration" = 23,
    "minimum occupancy" = 8  # or whatever shape you used
  )) +
  # ylim(c(-0.5, 100)) +
  #ggtitle(paste0("Trajectories, ",  thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp," % degradation"))+
  ylab("Metapopulation Extinction %") +
  xlab("Year") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=12, margin = margin(10,0,0,0)),
        axis.title.y=element_text(size=12, margin = margin(0,10,0,0)),
        strip.text = element_text(size = 12),
        plot.title=element_text(size=14),
        #plot.subtitle = element_text(size = 11),
        legend.position = "top",
        legend.direction = "vertical",
        legend.title =  element_text(size = 12),
        legend.text = element_text(size = 12))

# plotting all variable combinations in a loop ----------------------------
# all alphas, e values, c values, and degradation percentages

alphas <- unique(cumulative.extinct.by.year$alphaCat)
es <- unique(cumulative.extinct.by.year$e)
cs <- unique(cumulative.extinct.by.year$c)
deg_percents <- unique(cumulative.extinct.by.year$deg_percent)

# Output folder
out_dir <- "./path_to_output/results_inspection/plots/metapop_ext_perc_plots/"  # Change this to your desired output folder

# Loop over all combinations
for (thisalpha in alphas) {
  for (thise in es) {
    for (thisc in cs) {
      for (thisdp in deg_percents) {
        
        temp <- cumulative.extinct.by.year %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year > 50)
        
        temp2 <- mean.cumulative.extinct.perc %>%
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
        if (nrow(temp2) == 0) next  # Skip if no data
        if (nrow(temp3) == 0) next  # Skip if no data
        
        plot_title <- paste0("Trajectories, ", thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp, "% degradation")
        filename_base <- paste0("metapop_ext_", thisalpha, "_e", thise, "_c", thisc, "_", thisdp, "p")
        
        p <- ggplot() +
          geom_vline(xintercept = c(100, 200, 250), col="grey")+
          # geom_line(data=temp, aes(x=year,  y = cumulative_extinct_perc, group = interaction(landscapeID,scenario),lty = scenario),  col = 'gray',  alpha = 0.3, lwd = 0.5) +
          geom_line(data=temp2, aes(x=year, y = mean_cumulative_extinct_perc, group = scenario, col=scenario, lty = scenario), lwd = 1, alpha = 1) + 
          geom_point(data=temp3, aes(x=year, y= -3, group = scenario, pch=phase), col ="black", fill="grey", alpha = 0.9, cex = 1.5, stroke = 0.2) + 
          
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
          scale_y_continuous(name = "Metapopulation Extinction %", breaks = seq(0, 100, by = 25), limits = c(-3, 100)) +
          # ylim(c(-0.5, 100)) +
          #ggtitle(paste0("Trajectories, ",  thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp," % degradation"))+
          ylab("Metapopulation Extinction %") +
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


