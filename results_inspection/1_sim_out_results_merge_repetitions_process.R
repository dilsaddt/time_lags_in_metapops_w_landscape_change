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
# This script is to process the simulation output 
# for further inspection
#
# Date: 20/05/2025
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data ----
###########################################################################

rm(list=ls())

## 1.2. Loading libraries ----
##############################

library(ggplot2)
library(ggthemes)
library(gridExtra)
library(dplyr)
library(viridis)

###########################################################################
# 2. Baseline degradation ----
###########################################################################

###########################################################################
# 2.1. csv output of the simulations ----
###########################################################################
outdir <- "./path_to_output/land_degradation/sim_out/"
thesefiles <- list.files(outdir)
thesefiles

# first script of 10 landscape repetitions
simfile <- grep(pattern = "_1_", thesefiles) # simulation
simfiles <- thesefiles[simfile]
simfiles

temp <- list(NA)
for(i in 1:length(simfiles)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles[i]), header=T, sep=",")
}

temp_out_1 <- do.call(rbind, temp)
str(temp_out_1)

# second script of 10 landscape repetitions
simfile_temp <- grep(pattern = "_2_", thesefiles) # simulation
simfiles_temp <- thesefiles[simfile_temp]
simfiles_temp

temp <- list(NA)
for(i in 1:length(simfiles_temp)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles_temp[i]), header=T, sep=",")
}

temp_out_2 <- do.call(rbind, temp)
temp_out_2$landscapeID <- temp_out_2$landscapeID + 10
table(temp_out_2$landscapeID)
str(temp_out_2)

# merge them
out_base <- rbind(temp_out_1, temp_out_2)
out_base <- out_base[order(out_base$landscapeID),]

# add other columns
out_base$scenario <- "Baseline degradation"

out_base$alphaCat <- ifelse(out_base$alpha == 1/200, "a200", ifelse(out_base$alpha == 1/500, "a500", "a1000"))
out_base$alphaCat <- as.factor(out_base$alphaCat)
out_base$alphaCat <- factor(out_base$alphaCat, levels =  c("a200", "a500", "a1000"))
levels(out_base$alphaCat)

# check
str(out_base)
#table(out_base$metapopID)


###########################################################################
# 2.2. lambdaM outputs of the simulations ----
###########################################################################
# load the fixed parameters if saved
load("./path_to_output/land_degradation/fixed_params/sim_landdeg_1_a500_80p_fixed_params.RData")
outdir <- "./path_to_output/land_degradation/lambdaM/"
thesefiles <- list.files(outdir)
thesefiles

# first script of 10 landscape repetitions
simfile <- grep(pattern = "_1_", thesefiles) # simulation
simfiles <- thesefiles[simfile]
simfiles

temp <- list(NA)
for(i in 1:length(simfiles)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles[i]), header=T, sep=",")
  temp[[i]]$year <- rep(1:ntime, times = nlands)
}

temp_out_1 <- do.call(rbind, temp)
temp_out_1$alphaCat <- rep(c("a1000", "a200", "a500"), each = (3*nrow(temp[[1]])))
temp_out_1$deg_percent <- rep(c(20, 50, 80), each = nrow(temp[[2]]))

str(temp_out_1)

# second script of 10 landscape repetitions
simfile_temp <- grep(pattern = "_2_", thesefiles) # simulation
simfiles_temp <- thesefiles[simfile_temp]
simfiles_temp

temp <- list(NA)
for(i in 1:length(simfiles_temp)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles_temp[i]), header=T, sep=",")
  temp[[i]]$year <- rep(1:ntime, times = nlands)
  
}

temp_out_2 <- do.call(rbind, temp)
temp_out_2$landscapeID <- temp_out_2$landscapeID + 10
table(temp_out_2$landscapeID)

temp_out_2$alphaCat <- rep(c("a1000", "a200", "a500"), each = (3*nrow(temp[[1]])))
temp_out_2$deg_percent <- rep(c(20, 50, 80), each = nrow(temp[[1]]))

str(temp_out_2)

# merge them
lambdaM_base <- rbind(temp_out_1, temp_out_2)

# add other columns
lambdaM_base$scenario <- "Baseline degradation"

# check
str(lambdaM_base)

###########################################################################
# 3. Restoration ----
###########################################################################

###########################################################################
# 3.1. csv output of the simulations ----
###########################################################################
outdir <- "./path_to_output/land_restoration/sim_out/"
thesefiles <- list.files(outdir)
thesefiles

# first script of 10 landscape repetitions
simfile <- grep(pattern = "_1_", thesefiles) # simulation
simfiles <- thesefiles[simfile]
simfiles

temp <- list(NA)
for(i in 1:length(simfiles)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles[i]), header=T, sep=",")
}

temp_out_1 <- do.call(rbind, temp)
str(temp_out_1)

# second script of 10 landscape repetitions
simfile_temp <- grep(pattern = "_2_", thesefiles) # simulation
simfiles_temp <- thesefiles[simfile_temp]
simfiles_temp

temp <- list(NA)
for(i in 1:length(simfiles_temp)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles_temp[i]), header=T, sep=",")
}

temp_out_2 <- do.call(rbind, temp)
temp_out_2$landscapeID <- temp_out_2$landscapeID + 10
table(temp_out_2$landscapeID)
str(temp_out_2)

# merge them
out_rest <- rbind(temp_out_1, temp_out_2)
out_rest <- out_rest[order(out_rest$landscapeID),]

# add other columns
out_rest$scenario <- "Restoration"

out_rest$alphaCat <- ifelse(out_rest$alpha == 1/200, "a200", ifelse(out_rest$alpha == 1/500, "a500", "a1000"))
out_rest$alphaCat <- as.factor(out_rest$alphaCat)
out_rest$alphaCat <- factor(out_rest$alphaCat, levels =  c("a200", "a500", "a1000"))
levels(out_rest$alphaCat)

# check
str(out_rest)
#table(out_rest$metapopID)


###########################################################################
# 3.2. lambdaM outputs of the simulations ----
###########################################################################
load("./path_to_output/sim_landrest_1_a500_80p_fixed_params.RData")
outdir <- "./path_to_output/land_restoration/lambdaM/"
thesefiles <- list.files(outdir)
thesefiles

# first script of 10 landscape repetitions
simfile <- grep(pattern = "_1_", thesefiles) # simulation
simfiles <- thesefiles[simfile]
simfiles

temp <- list(NA)
for(i in 1:length(simfiles)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles[i]), header=T, sep=",")
  temp[[i]]$year <- rep(1:ntime, times = nlands)
}

temp_out_1 <- do.call(rbind, temp)
temp_out_1$alphaCat <- rep(c("a1000", "a200", "a500"), each = (3*nrow(temp[[1]])))
temp_out_1$deg_percent <- rep(c(20, 50, 80), each = nrow(temp[[2]]))

str(temp_out_1)

# second script of 10 landscape repetitions
simfile_temp <- grep(pattern = "_2_", thesefiles) # simulation
simfiles_temp <- thesefiles[simfile_temp]
simfiles_temp

temp <- list(NA)
for(i in 1:length(simfiles_temp)){
  temp[[i]] <- read.table(file=paste0(outdir, simfiles_temp[i]), header=T, sep=",")
  temp[[i]]$year <- rep(1:ntime, times = nlands)
  
}

temp_out_2 <- do.call(rbind, temp)
temp_out_2$landscapeID <- temp_out_2$landscapeID + 10
table(temp_out_2$landscapeID)

temp_out_2$alphaCat <- rep(c("a1000", "a200", "a500"), each = (3*nrow(temp[[1]])))
temp_out_2$deg_percent <- rep(c(20, 50, 80), each = nrow(temp[[1]]))

str(temp_out_2)

# merge them
lambdaM_rest <- rbind(temp_out_1, temp_out_2)

# add other columns
lambdaM_rest$scenario <- "Restoration"

# check
str(lambdaM_rest)

###########################################################################
# 3. MERGE TWO SCENARIOS ----
###########################################################################
out <- rbind(out_base, out_rest)
lambdaM_df <- rbind(lambdaM_base, lambdaM_rest)

###########################################################################
# 4. SAVE ----
###########################################################################

save(out, file=("/path_to_output/results_inspection/sim_out.RData"))
save(lambdaM_df, file=("path_to_output/results_inspection/sim_lambdaM.RData"))

