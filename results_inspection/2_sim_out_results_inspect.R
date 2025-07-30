############################################################################
# Data simulation script for the theoretical simulation models based on
# spatially explicit version of Levins' metapopulation model 
# using landscape_simulation_functions.R and sim_obsDat_function.R
#
# investigating simulated data for two scenarios
# land degradation (landdeg)
# land degradation + restoration (landrest)
# checking time series for simulated populations
#
# This script is to prep and calculate derived parameters
# to investigate and plot the simulation output 
#
# Date: 20/05/2025
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data ----
###########################################################################

rm(list=ls())

## 1.1. Loading libraries ----
##############################

library(ggplot2)
library(ggthemes)
library(gridExtra)
library(dplyr)
library(viridis)

## 1.2. Load data ----
##################################################################
load("./path_to_output/results_inspection/sim_out.RData")
str(out)

###########################################################################
# 2. Calculate values for summarizing ----
###########################################################################

## 2.1. Proportion of occupied sites ----
#########################################

# average across metapopulations but keep the landscape repetitions -------

propocc <- aggregate(out$n.occ.prop, list(out$landscapeID, out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=mean, na.rm=T)
names(propocc) <- c("landscapeID", "deg_percent", "year", "e", "c", "alphaCat", "scenario", "propocc.mean")

## add standard deviation for a simple measure of spread 
sd.propocc <- aggregate(out$n.occ.prop, list(out$landscapeID, out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=sd, na.rm=T)
propocc$sd <- sd.propocc$x
rm(sd.propocc); gc()

## add e-c combination for easier plotting:
propocc$e_c <- paste(propocc$e, propocc$c, sep="_")

# factors for plotting
propocc$e <- factor(propocc$e )
propocc$c <- factor(propocc$c )

# round propocc.mean to 3 digits
propocc$propocc.mean <- round(propocc$propocc.mean, digits = 3)

### store these dataframes and load instead later (it takes a couple of minutes to generate them):
save(propocc, file=("./path_to_output/results_inspection/Propoccs_wLandscapeReps_200525.RData"))


# average across both metapopulations and landscapes ----------------------

propocc.mean <- aggregate(out$n.occ.prop, list(out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=mean, na.rm=T)
names(propocc.mean) <- c("deg_percent", "year", "e", "c", "alphaCat", "scenario", "propocc.mean")

## add standard deviation for a simple measure of spread 
sd.propocc.mean <- aggregate(out$n.occ.prop, list(out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=sd, na.rm=T)
propocc.mean$sd <- sd.propocc.mean$x
rm(sd.propocc.mean); gc()

## add e-c combination for easier plotting:
propocc.mean$e_c <- paste(propocc.mean$e, propocc.mean$c, sep="_")

# factors for plotting
propocc.mean$e <- factor(propocc.mean$e )
propocc.mean$c <- factor(propocc.mean$c )

# round propocc.mean to 3 digits
propocc.mean$propocc.mean <- round(propocc.mean$propocc.mean, digits = 3)

### store these dataframes and load instead later (it takes a couple of minutes to generate them):
save(propocc.mean, file=("./path_to_output/results_inspection/PropoccsMean_200525.RData"))

## 2.2. Number of occupied sites ----
#########################################

# average across metapopulations but keep the landscape repetitions -------

nocc <- aggregate(out$n.occ, list(out$landscapeID, out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=mean, na.rm=T)
names(nocc) <- c("landscapeID", "deg_percent", "year", "e", "c", "alphaCat", "scenario", "nocc.mean")

## add standard deviation for a simple measure of spread 
sd.nocc <- aggregate(out$n.occ, list(out$landscapeID, out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=sd, na.rm=T)
nocc$sd <- sd.nocc$x
rm(sd.nocc); gc()

## add e-c combination for easier plotting:
nocc$e_c <- paste(nocc$e, nocc$c, sep="_")

# factors for plotting
nocc$e <- factor(nocc$e )
nocc$c <- factor(nocc$c )

# round propocc.mean to 3 digits
nocc$nocc.mean <- round(nocc$nocc.mean) # to make them integers

### store these dataframes and load instead later (it takes a couple of minutes to generate them):
save(nocc, file=("./path_to_output/results_inspection/Noccs_wLandscapeReps_200525.RData"))

# average across both metapopulations and landscapes ----------------------

nocc.mean <- aggregate(out$n.occ, list(out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=mean, na.rm=T)
names(nocc.mean) <- c("deg_percent", "year", "e", "c", "alphaCat", "scenario", "nocc.mean")

## add standard deviation for a simple measure of spread 
sd.nocc.mean <- aggregate(out$n.occ, list(out$deg_percent, out$year, out$e, out$c, out$alphaCat, out$scenario), FUN=sd, na.rm=T)
nocc.mean$sd <- sd.nocc.mean$x
rm(sd.nocc.mean); gc()

## add e-c combination for easier plotting:
nocc.mean$e_c <- paste(nocc.mean$e, nocc.mean$c, sep="_")

# factors for plotting
nocc.mean$e <- factor(nocc.mean$e )
nocc.mean$c <- factor(nocc.mean$c )

# round propocc.mean to 3 digits
nocc.mean$nocc.mean <- round(nocc.mean$nocc.mean) # to make them integers

### store these dataframes and load instead later (it takes a couple of minutes to generate them):
save(nocc.mean, file=("./path_to_output/results_inspection/NoccsMean_200525.RData"))

## 2.3. For percengtage of extinct metapopulations ----
########################################################

# year, metapopID, e, c, alpha, npatch, n.occ, n.occ.prop, landscapeID, deg_percent, scenario, alphaCat
str(out)

out.sim <- out[, c(1:8,12,14:16)]
str(out.sim)


# add simulation nr -------------------------------------------------------
out.sim$sim.nr.name <- paste0("land ", out.sim$landscapeID, " | metapop ",  out.sim$metapopID)
length(unique(out.sim$sim.nr.name))
out.sim$sim.nr <- rep(1:length(unique(out.sim$sim.nr.name)), each = table(out.sim$sim.nr.name)[1] )

# add extinct column ------------------------------------------------------

out.sim$uniquecombi <- paste(out.sim$deg_percent, out.sim$e, out.sim$c, out.sim$alphaCat, sep="_")
out.sim$extinct <- ifelse(out.sim$n.occ == 0, 1, 0) # 1 means extinct
table(out.sim$extinct)

# calculate cumulative metapopulations that went extinct by year ----------

# STEP 1: Find first extinction year for each metapopID (within each group)
extinction.years <- out.sim %>%
  filter(extinct == 1) %>%
  group_by(landscapeID, alphaCat, e, c, deg_percent, metapopID, scenario) %>%
  summarise(first_extinction_year = min(year), .groups = "drop")

# STEP 2: Create all years per group to support cumulative counting
all.years <- out.sim %>%
  distinct(year, landscapeID, alphaCat, e, c, deg_percent, scenario)

# STEP 3: Count how many metapopIDs went extinct *by or before* each year
cumulative.extinct.by.year <- all.years %>%
  group_by(landscapeID, alphaCat, e, c, deg_percent, year, scenario) %>%
  summarise(
    cumulative_extinct = sum(
      extinction.years$first_extinction_year <= year &
        extinction.years$landscapeID == unique(landscapeID) &
        extinction.years$alphaCat == unique(alphaCat) &
        extinction.years$e == unique(e) &
        extinction.years$c == unique(c) &
        extinction.years$deg_percent == unique(deg_percent) &
        extinction.years$scenario == unique(scenario)
    ),
    .groups = "drop"
  )


str(cumulative.extinct.by.year)

# save these --------------------------------------------------------------

save(out.sim, cumulative.extinct.by.year, file=("./path_to_output/results_inspection/ForMetapopExtPercentage_200525.RData"))


## 2.4. npatch per year calculation ----
########################################################

# prep and save npatch per year depending on deg_percent

# degradation only scenario -----------------------------------------------
out.base <- out[out$scenario == "Baseline degradation",]
dp20_temp <- out.base[out.base$deg_percent == 20,]
dp20_temp <- dp20_temp[!duplicated(dp20_temp$year),]
dp20_temp <- select(dp20_temp, c("year", "deg_percent", "npatch"))

dp50_temp <- out.base[out.base$deg_percent == 50,]
dp50_temp <- dp50_temp[!duplicated(dp50_temp$year),]
dp50_temp <- select(dp50_temp, c("year", "deg_percent", "npatch"))

dp80_temp <- out.base[out.base$deg_percent == 80,]
dp80_temp <- dp80_temp[!duplicated(dp80_temp$year),]
dp80_temp <- select(dp80_temp, c("year", "deg_percent", "npatch"))


plot(dp20_temp$year, dp20_temp$npatch, type = 'p', ylim = c(0, 600))
points(dp50_temp$year, dp50_temp$npatch, col = 'red')
points(dp80_temp$year, dp80_temp$npatch, col = 'blue')


npatch.base <- rbind(dp20_temp, dp50_temp, dp80_temp)
plot(npatch.base$year, npatch.base$npatch, col = npatch.base$deg_percent)

write.csv(npatch.base, file="./path_to_output/results_inspection/npatch_base_200525.csv", row.names = F)

# restoration scenario ----------------------------------------------------
out.restor <- out[out$scenario == "Restoration",]
dp20_temp <- out.restor[out.restor$deg_percent == 20,]
dp20_temp <- dp20_temp[!duplicated(dp20_temp$year),]
dp20_temp <- select(dp20_temp, c("year", "deg_percent", "npatch"))

dp50_temp <- out.restor[out.restor$deg_percent == 50,]
dp50_temp <- dp50_temp[!duplicated(dp50_temp$year),]
dp50_temp <- select(dp50_temp, c("year", "deg_percent", "npatch"))

dp80_temp <- out.restor[out.restor$deg_percent == 80,]
dp80_temp <- dp80_temp[!duplicated(dp80_temp$year),]
dp80_temp <- select(dp80_temp, c("year", "deg_percent", "npatch"))


plot(dp20_temp$year, dp20_temp$npatch, type = 'p', ylim = c(0, 600))
points(dp50_temp$year, dp50_temp$npatch, col = 'red')
points(dp80_temp$year, dp80_temp$npatch, col = 'blue')


npatch.restor <- rbind(dp20_temp, dp50_temp, dp80_temp)
plot(npatch.restor$year, npatch.restor$npatch, col = npatch.restor$deg_percent)

write.csv(npatch.restor, file="./path_to_output/results_inspection/npatch_restor_200525.csv", row.names = F)

