# time_lags_in_metapops_w_landscape_change
Descriptipn: Simulation models using metapopulation modeling framework to investigate relationship between species-traits (colonisation, extinction, dipsersal) and time lags (extinction debt and colonisation lag) in response to landscape degradation and restoration.


**Manuscript in preparation:**  
Dagtekin, D., Moor, H., & Sutherland, C. (2025).  
_**Extinction debts and colonisation lags in the response of spatially structured populations to degradation and restoration**_.  
_To be submitted to: Ecology Letters_. [OSF Project Link](https://osf.io/xxxxx)

---

## Summary

This repository contains all code and supporting files for simulating extinction debts and colonisation lags in metapopulations exposed to habitat degradation and restoration.

We use a spatially explicit, stochastic patch-occupancy model to explore how species traits (dispersal, extinction, colonisation) and landscape change dynamics interact to produce time-delayed responses in biodiversity.

---

## Repository Structure

### Simulation functions: R functions used in simulations scripts

_landscape_simulation_functions_290724.R:_ Simulating landscapes with random patch distribution with log-normal distributed patch areas, adapted from Walker & Gilbert, 2023, Ecology, doi: 10.1002/ecy.3840

_PStarFunction_290724.R:_ Initial occupancy probability calculator for initial, stable state in simulations, adapted from Walker & Gilbert, 2023, Ecology, doi: 10.1002/ecy.3840

_sim_obsDat_landdeg_function_faster_100425.R:_ Simulation function for metapopulation under landscape degradation, adapted from Sutherland et al., 2014, Ecology, doi: 10.1890/14-0384.1 & Bertasello et al., 2021, Royal Soc. Open Sci. doi: 10.1098/rsos.201309

_sim_obsDat_landrest_function_faster_100425.R:_ Simulation function for metapopulation under landscape degradation and restoration, adapted from Sutherland et al., 2014, Ecology, doi: 10.1890/14-0384.1 & Bertasello et al., 2021, Royal Soc. Open Sci. doi: 10.1098/rsos.201309

### Simulation scripts: R functions for simulation

The scripts here are examples for only intermediate dispersal (alpha = 1/500) and degradation percentage = 80. All of them are used to simulate combinations with different species traits and degradation percentages. This includes low-, mid-, high-values of these parameters:
alpha = {1/200, 1/500, 1/1000}
e = {0.2, 0.5, 0.8}
c = {0.4, 1, 1.6}
degradation percentage = {20, 50, 80}

These scripts include 10 landscape repetitions and 50 metapopulation repetitions in each landscape. For the whole study in total 20 landscape repetitions and 50 metapopulation repetitions in each landscape were simulated.

_sim_landdeg_a500_80p.R:_ Simulation script for landscape degradation only scenario.
Example script for only intermediate dispersal (alpha = 1/500) and degradation percentage = 80.

### Results inspection: Scripts for processing simulation ouputs

_1_sim_out_results_merge_repetitions_process.R:_ Script for checking simulation output and merging landscape repetitions (2 scripts, each with 10 reps) for each scenario (degradation only & degradation + restoration) 

_2_sim_out_results_inspect:_ Script for calculating values for results and visualization for each scenario (degradation only & degradation + restoration), e.g., proportion of occupied sites, number of occupied sites, number of patches in each year

### Visualization: Scripts for visualizing ouputs
All plots were done fore each species and scenarios. Each species meaning: all combinations or e, c, and alpha.

_3_propocc_time_series.R:_ Script for plotting proportion of occupied sites over time.
_4_percentage_of_extinct_metapops.R:_ Script for metapopulation extinction percentage over time.
_5_hysteresis.R:_ Script for hysteresis plots - proportion of occupied sites over number of available patches.
_6_landscape_change_lambdaM_plots.R:_ Script for checking metapopulation capacity change over time with landscape change (degradation and restoration). This was used in methods to show the change in landscape and its relationship with metapopulation capacity, not in results.




