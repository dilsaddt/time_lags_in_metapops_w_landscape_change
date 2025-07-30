############################################################################
# This script is to plot the proportion of occupied sites over npatch change
# aka hysteresis plots
# 
# Date: 20/05/2025
# Author: Dilsad Dagtekin

############################################################################
rm(list=ls())

library(ggplot2)
library(gridExtra)
library(dplyr)

########################################################################################################################


# load fixed parameters to be used in plotting, if saved
load("./path_to_output/land_restoration/fixed_params/sim_landrest_1_a500_80p_fixed_params.RData")
rm(list=setdiff(ls(), c("ntime_first", "ntime_deg", "ntime_rest", "ntime_stop"))) # keep only out to begin with
gc()

year_deg_start <- ntime_first # end of initial phase, start of degradation phase
year_rest_start <- ntime_first + ntime_deg # end of degradation phase, start of restoration phase
year_rest_end <- ntime_first + ntime_deg + ntime_rest # end of restoration phase
year_sim_end <- ntime_first + ntime_deg + ntime_rest + ntime_stop # end of simulation


### load stored calculated values to makes this faster:
load("./path_to_output/results_inspection/Propoccs_wLandscapeReps.RData") # with metapop repetitions averaged over landscape repetitions
str(propocc)
names(propocc)[9] <- "propocc.sd"

propocc.base <- propocc[propocc$scenario == "Baseline degradation",]
propocc.restor <- propocc[propocc$scenario == "Restoration",]

load("./path_to_output/results_inspection/PropoccsMean.RData") # averaged over both metapop and lanscape repetitions
str(propocc.mean)
names(propocc.mean)[8] <- "propocc.sd"

mean.propocc.base <- propocc.mean[propocc.mean$scenario == "Baseline degradation",]
mean.propocc.restor <- propocc.mean[propocc.mean$scenario == "Restoration",]

npatch.base <- read.csv("./path_to_output/results_inspection/npatch_base.csv")
str(npatch.base)
npatch.restor <- read.csv("./path_to_output/results_inspection/npatch_restor.csv")
str(npatch.restor)


# degradation scenario ----------------------------------------------------
# combine prop.occ.mean:

# only take after year 50, assuming all reached stable state after then
hysteresis.propocc.base <- propocc.base[propocc.base$year >=50,]

# add the number of available sites:
hysteresis.propocc.base <- hysteresis.propocc.base %>%
  left_join(npatch.base, by = c("year", "deg_percent"))

hysteresis.propocc.base$uniquecombi <- paste(hysteresis.propocc.base$deg_percent, hysteresis.propocc.base$e, hysteresis.propocc.base$c, hysteresis.propocc.base$alphaCat, sep="_")

# plot(hysteresis.propocc.base$year, hysteresis.propocc.base$npatch, col = hysteresis.propocc.base$deg_percent)

hysteresis.mean.propocc.base <- mean.propocc.base[mean.propocc.base$year >=50,]

# add the number of available sites:
hysteresis.mean.propocc.base <- hysteresis.mean.propocc.base %>%
  left_join(npatch.base, by = c("year", "deg_percent"))

hysteresis.mean.propocc.base$uniquecombi <- paste(hysteresis.mean.propocc.base$deg_percent, hysteresis.mean.propocc.base$e, hysteresis.mean.propocc.base$c, hysteresis.mean.propocc.base$alphaCat, sep="_")

hysteresis.propocc.base$scenario <- "Baseline degradation"
hysteresis.mean.propocc.base$scenario <- "Baseline degradation"

min.propocc.base <- hysteresis.mean.propocc.base %>%
  group_by(e, c, uniquecombi, deg_percent) %>%
  slice_min(propocc.mean, with_ties = FALSE) %>%  # Only the row with the true min
  ungroup()
min.propocc.base$phase <- "minimum occupancy"


# restoration scenario ----------------------------------------------------
# add the number of available sites:

# only take after year 50, assuming all reached stable state after then
hysteresis.propocc.restor <- propocc.restor[propocc.restor$year >=50,]

hysteresis.propocc.restor <- hysteresis.propocc.restor %>%
  left_join(npatch.restor, by = c("year", "deg_percent"))

hysteresis.propocc.restor$uniquecombi <- paste(hysteresis.propocc.restor$deg_percent, hysteresis.propocc.restor$e, hysteresis.propocc.restor$c, hysteresis.propocc.restor$alphaCat, sep="_")

# plot(hysteresis.propocc.restor$year, hysteresis.propocc.restor$npatch, col = hysteresis.propocc.restor$deg_percent)

hysteresis.mean.propocc.restor <- mean.propocc.restor[mean.propocc.restor$year >=50,]

# add the number of available sites:
hysteresis.mean.propocc.restor <- hysteresis.mean.propocc.restor %>%
  left_join(npatch.restor, by = c("year", "deg_percent"))

hysteresis.mean.propocc.restor$uniquecombi <- paste(hysteresis.mean.propocc.restor$deg_percent, hysteresis.mean.propocc.restor$e, hysteresis.mean.propocc.restor$c, hysteresis.mean.propocc.restor$alphaCat, sep="_")

hysteresis.propocc.restor$scenario <- "Restoration"
hysteresis.mean.propocc.restor$scenario <- "Restoration"

min.propocc.restor <- hysteresis.mean.propocc.restor %>%
  group_by(e, c, uniquecombi, deg_percent) %>%
  slice_min(propocc.mean, with_ties = FALSE) %>%  # Only the row with the true min
  ungroup()
min.propocc.restor$phase <- "minimum occupancy"


# put both scenarios together ---------------------------------------------

hysteresis.propocc <- rbind(hysteresis.propocc.base, hysteresis.propocc.restor)
hysteresis.mean.propocc <- rbind(hysteresis.mean.propocc.base, hysteresis.mean.propocc.restor)
min.propocc <- rbind(min.propocc.base, min.propocc.restor)


# add the phase marks to merged data frame
# phases are:
# start of degradation
# end of degradation / start of restoration
# end of restoration

year_phase <- data.frame(NA)

year_phase <- rbind(hysteresis.mean.propocc[which(hysteresis.mean.propocc$year == year_deg_start),],
                    hysteresis.mean.propocc[which(hysteresis.mean.propocc$year == year_rest_start),],
                    hysteresis.mean.propocc[which(hysteresis.mean.propocc$year == year_rest_end),])

year_phase$phase <- NA


year_phase$phase <- ifelse(year_phase$year == year_deg_start, "start of degradation",
                           ifelse(year_phase$year == year_rest_start, "end of degradation / start of restoration", "end of restoration"))

year_phase$phase <- as.factor(year_phase$phase)
year_phase$phase <- factor(year_phase$phase, levels = c("start of degradation", "end of degradation / start of restoration", "end of restoration"))


# plot -------------------------------------------------------------

## remove first 90 years to get to the start of degradation closer (year 100):
## subset to one alpha, e, c, degradation percentage first to check

thisalpha <- "a500"
thise <- "0.8"
thisc <- "0.4"
thisdp <- 80

# also not plotting c = 1.6 and deg_perc = 20 below

temp <- hysteresis.propocc[hysteresis.propocc$alphaCat==thisalpha & 
                             hysteresis.propocc$e == thise & 
                             hysteresis.propocc$c == thisc & 
                             hysteresis.propocc$deg_percent == thisdp & 
                             hysteresis.propocc$year>90,] 

temp2 <- hysteresis.mean.propocc[hysteresis.mean.propocc$alphaCat==thisalpha & 
                                   hysteresis.mean.propocc$e == thise & 
                                   hysteresis.mean.propocc$c == thisc & 
                                   hysteresis.mean.propocc$deg_percent == thisdp & 
                                   hysteresis.mean.propocc$year>90,] 

temp3 <- year_phase[year_phase$alphaCat==thisalpha & 
                      year_phase$e == thise & 
                      year_phase$c == thisc & 
                      year_phase$deg_percent == thisdp & 
                      year_phase$scenario == "Restoration" &
                      year_phase$year>90,]

temp4 <- min.propocc[min.propocc$alphaCat==thisalpha & 
                       min.propocc$e == thise & 
                       min.propocc$c == thisc & 
                       min.propocc$deg_percent == thisdp & 
                       min.propocc$year>90,]

# this plot is just to check and get the legend

p_legend <- ggplot() +
  #  geom_vline(xintercept = c(100, 200), col="grey")+
  
  #geom_path(data=temp, aes(x=npatch,  y = propocc.mean, group = scenario), col = 'grey',  alpha = 0.5, lwd = 1) + # lwd = 0.5  , linetype=c
  geom_path(data=temp2, aes(x=npatch,  y = propocc.mean, group = scenario, col=scenario, lty = scenario), alpha = 1, lwd = 1) + # lwd = 0.5  , linetype=c
  
  geom_point(data=temp4, aes(x=npatch, y=propocc.mean, group = scenario, col=scenario, shape=phase),  cex = 3, stroke = 2) + #show.legend = FALSE) + 
  
  scale_color_brewer(palette = "Dark2", direction = -1) +  # Color-blind friendly categorical palette
  scale_fill_brewer(palette = "Dark2") +  # Color-blind friendly categorical palette
  scale_linetype_manual(
    values = c("Baseline degradation" = "solid", "Restoration" = "dashed"),
    breaks = c("Baseline degradation", "Restoration"),  # order in the legend
    name = "scenario"
  ) +
  geom_point(data=temp3, aes(x=npatch, y=propocc.mean, group = scenario, pch=phase), col ="black", fill="grey", alpha = 0.7, cex = 3, stroke = 1) + 
  scale_shape_manual(values = c(
    "start of degradation" = 21,
    "end of degradation / start of restoration" = 22,
    "end of restoration" = 23,
    "minimum occupancy" = 8  # or whatever shape you used
  )) +
  ylim(c(0, 0.6)) +
  ggtitle(paste0("Trajectories, ", thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp," % degradation"))+
  ylab("Proportion of Occupied Sites") +
  xlab("Number of Patches") +
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
        legend.title =  element_text(size = 12),
        legend.text = element_text(size = 12))


# save plot legend
out_dir <- "./path_to_output/results_inspection/plots/hysteresis_plots/"  # Change this to your desired output folder

filename_base <- paste0("hysteresis_", thisalpha, "_e", thise, "_c", thisc, "_", thisdp, "p")

# Save as PDF
ggsave(filename = file.path(out_dir, paste0("legend_", filename_base, ".pdf")),
       plot = p_legend, width = 8, height = 6)


# plotting all variable combinations in a loop ----------------------------
# all alphas, e values, c values, and degradation percentages

alphas <- unique(hysteresis.propocc$alphaCat)
es <- unique(hysteresis.propocc$e)
cs <- unique(hysteresis.propocc$c)
deg_percents <- unique(hysteresis.propocc$deg_percent)

# Output folder
out_dir <- "./path_to_output/results_inspection/plots/hysteresis_plots/"  # Change this to your desired output folder

# Loop over all combinations
for (thisalpha in alphas) {
  for (thise in es) {
    for (thisc in cs) {
      for (thisdp in deg_percents) {
        
        
        temp <- hysteresis.propocc %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year > 90)
        
        temp2 <- hysteresis.mean.propocc %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year > 90)
        
        temp3 <- year_phase %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year_phase$scenario == "Restoration",
                 year > 90)
        
        temp4 <- min.propocc %>%
          filter(alphaCat == thisalpha,
                 e == thise,
                 c == thisc,
                 deg_percent == thisdp,
                 year > 90)
        
        temp <- temp[order(temp[,"year"]), ]
        temp2 <- temp2[order(temp2[,"year"]), ]
        temp3 <- temp3[order(temp3[,"year"]), ]
        
        #if (nrow(temp) == 0) next  # Skip if no data
        #if (nrow(temp2) == 0) next  # Skip if no data
        #if (nrow(temp3) == 0) next  # Skip if no data
        #if (nrow(temp4) == 0) next  # Skip if no data
        
        plot_title <- paste0("Trajectories, ", thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp," % degradation")
        filename_base <- paste0("hysteresis_", thisalpha, "_e", thise, "_c", thisc, "_", thisdp, "p")
        
        p <- ggplot() +
          geom_path(data=temp2, aes(x=npatch,  y = propocc.mean, group = scenario, col=scenario, lty = scenario), alpha = 1, lwd = 1) + # lwd = 0.5  , linetype=c
          scale_color_brewer(palette = "Dark2", direction = -1) +  # Color-blind friendly categorical palette
          scale_linetype_manual(
            values = c("Baseline degradation" = "solid", "Restoration" = "dashed"),
            breaks = c("Baseline degradation", "Restoration"),  # order in the legend
            name = "scenario"
          ) +
          geom_point(data=temp3, aes(x=npatch, y=propocc.mean, group = scenario, pch=phase), col ="black", fill="grey", alpha = 0.7, cex = 2, stroke = 1) + 
          geom_point(data=temp4, aes(x=npatch, y=propocc.mean, group = scenario, col=scenario, shape=phase),  cex = 2, stroke = 1, alpha = 0.7) + #show.legend = FALSE) + 
          
          scale_shape_manual(values = c(
            "start of degradation" = 21,
            "end of degradation / start of restoration" = 22,
            "end of restoration" = 23,
            "minimum occupancy" = 8  # or whatever shape you used
          )) +
          ylim(c(0, 1)) +
          # scale_x_continuous(name = "Number of Patches",
          #                    breaks = seq(100, 600, by = 100),
          #                    limits = c(100, 600)) +
          scale_x_reverse(name = "Number of Patches", breaks = seq(100, 600, by = 100), limits = c(600,100)) +
          # ggtitle(paste0("Trajectories, ", thisalpha, " | e = ", thise, " | c = ", thisc, " | ", thisdp," % degradation"))+
          ylab("Proportion of Occupied Sites") +
          #xlab("Number of Patches") +
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

