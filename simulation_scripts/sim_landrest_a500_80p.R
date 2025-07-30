############################################################################
# Data simulation script for the theoretical simulation models based on
# spatially explicit version of Levins' metapopulation model 
# using landscape_simulation_functions.R and sim_obsDat_function.R
#
# simulating data for landscape restoration
# 
# nlands = 10 + 10 (2 scripts)
# nmetapops = 50 (in each land)
#
# Date: 15/05/2024
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data ----
###########################################################################

rm(list=ls())

set.seed(193)

## 1.1. Loading libraries ----
##############################

library(dplyr)

## 1.2. Set directory ----
##################################################################

# using relative paths in an Rstudio project makes it clear where stuff is and makes scripts portable
fundir <- "path_to_folder/simulation_functions/"

outdir <- "path_to_folder/simulation_out_150525/"
# dir.create(outdir) # run each time you create a new output folder

## 1.3. Load functions ----
##################################################################

## Landscape:
source(paste0(fundir, "landscape_simulation_functions_290724.R"))
# currently only using createLandscape(npatch, landscapeLimit, landscapeType, niter)
rm(createMultipleLandscapes)

## Initial occupancy function, adapted from Walker & Gilbert, 2023:
source(paste0(fundir, "PStarFunction_290724.R"))

## Data simulation:
source(paste0(fundir, "sim_obsDat_landrest_function_faster_100425.R"))

###########################################################################
# 2. Simulation ----
###########################################################################


## 2.1. Set path to output ----
###############################

## output: set a filename using the relative output path:
run_name <- "sim_landrest_a500_80p"
filename <- paste0(outdir, run_name, ".csv")
filename2 <- paste0(outdir, run_name, "_lambdaM.csv")


## 2.2. Fixed parameters ----
##############################

nlands <- 10 # number of landscape replications
nmetapops <- 50 # number of metapopulation replications per landscape
nsim <- nmetapops

npatch <- 600 # number of patches in a landscape
ntime <- 400 # number of total years, 400 for degradation + restoration scenario


alpha <- 1/500 # dispersal distance, change for low = 1/200, mid = 1/500, high = 1/1000
S_int <- 0 # intercept for suitability function
S_b1 <- 1 # coefficient for area in suitability function
x <- 1 # exponent of area effect on extinction equation

# for landscape simulation function:
landscapeLimit <- 10000 # meters, extent of the landscape along the x and y axes (ie the space within which patches will be distributed)
landscapeType <- "clustered" # wanted landscape type

# e and c parameters:
# 3 values each: low, mid, and high
e.values <- c(0.2, 0.5, 0.8) # scaling parameter for extinction prob
c.values <- c(0.4, 1, 1.6) # scaling parameter for colonization prob

## 2.3. Land degradation and restoration parameters ----
########################################################
ntime_first <- 100 # 100 years of stable landscape
ntime_deg <- 100 # 100 years of degradation
ntime_rest <- 50 # 50 years of restoration right after degradation
ntime_stop <- 150 # 150 years of no interference after stopping everything
# the total ntime should be equal to the sum of ntime_X above

deg_percent <- 80 # degradation severity, % of the patches will go in ntime_deg years, change for {low = 20, mid = 50, high = 80}

## 2.4. Empty data frame to collect results ----
################################################

# save the fixed and prpared parameters
save.image(paste0(outdir, run_name, "_fixed_params.RData"))

# simulations are done like: npatch x nlandscapes (nland) x nmetapops (nsim) x e-c combinations

results <- data.frame(  year = numeric(),
                        metapopID = integer(),
                        e = numeric(),  # species settings, e: extinction vulnerability
                        c = numeric(),  # species settings, c: colonization capability
                        alpha = numeric(),  # species settings, alpha: dispersal capability
                        npatch = numeric(),  # number of patches
                        n.occ = numeric(), # number of occupied patches
                        n.occ.prop = numeric(), # proportion of occupied patches
                        avg.psi = numeric(), # average occupancy probability
                        avg.gamma = numeric(), # average colonization probability
                        avg.eps = numeric(), # average extinction probability
                        landscapeID = integer(), # metapop replications per landscape and e-c combination (collect from nsim counter d, 3rd dim in out array)
                        lambdaM = numeric(), # metapopulation capacity
                        deg_percent = integer()) # degradation percent

l# save metapopulation capacity over time too
lambdaM_out <- data.frame(lambdaM = numeric(),
                          landscapeID = integer()) # degradation percent


# write the file with headers, still empty:
write.table(results, file = filename, append = F, sep = ",", col.names = T)
write.table(lambdaM_out, file = filename2, append = F, sep = ",", col.names = T)

# empty holders for probability results
psi_out_list <- gamma_out_list <- eps_out_list <- vector(mode = "list", length = nlands)

psi_out <- gamma_out <- eps_out <- NA


## 2.5. Simulate  ----
################################
start.time <- Sys.time()

for(land in 1:nlands){ # 10 different landscapes
  
  print(paste("Running landscape repetition", land, ":"))
  
  # calculate new lands
  lands_new <- patchID_new <- vector(mode = "list", length = (1 + ntime_deg + ntime_rest))
  dmat_new <- S_new <- vector(mode = "list", length = (1 + ntime_deg + ntime_rest))
  
  # naming based on years
  names(lands_new) <- (ntime_first):(ntime_first + ntime_deg + ntime_rest)
  names(dmat_new) <- (ntime_first):(ntime_first + ntime_deg + ntime_rest)
  names(S_new) <- (ntime_first):(ntime_first + ntime_deg + ntime_rest)
  names(patchID_new) <- (ntime_first):(ntime_first + ntime_deg + ntime_rest)
  
  ### initial landscape
  sim_land <- createLandscape(npatch = npatch, landscapeLimit = landscapeLimit, landscapeType = landscapeType, niter = 50)
  
  initial_land <- sim_land
  initial_dmat <- as.matrix(initial_land[, -c(1:4)])
  initial_area <- initial_land$A
  initial_A <- scale(initial_area)[,1]
  initial_S <- plogis(S_int + S_b1*initial_A) 
  initial_S_mat <- as.matrix(initial_S)
  rownames(initial_S_mat) <- initial_land$patchID
  
  ################################
  ## INITIAL LANDSCAPE  ----
  ################################
  
  #### fill the lists for stable phase ####
  ##########################################
  
  lands_new[[1]] <- initial_land
  dmat_new[[1]] <- initial_dmat
  S_new[[1]] <- initial_S_mat
  patchID_new[[1]] <- initial_land$patchID
  
  ################################
  ## DEGRADATION  ----
  ################################
  
  # how many years of degradation is ntime_deg
  # years of degradation
  deg_years <- (ntime_first + 1):(ntime_first + ntime_deg)
  
  # we choose which patch to remove randomly during degradation years 
  # out of total npatch (600), we lose deg_percent % of patches over ntime_deg years
  nr_patch_removed_per_year <- ceiling(npatch*(deg_percent/100)/ntime_deg)
  
  # which patches are removed, with patch IDs:
  to_remove <- rep(NA, length = ntime_deg)
  to_remove <- sample(initial_land$patchID, size = nr_patch_removed_per_year*ntime_deg, replace = FALSE) 
  
  # save the removed patch IDs in each year
  which_patch_removed_per_year <- split(to_remove, ceiling(seq_along(to_remove) / nr_patch_removed_per_year))
  names(which_patch_removed_per_year) <- deg_years
  which_patch_removed_per_year
  
  #### fill the lists for degradation phase ####
  ##############################################
  
  for(t in 1:ntime_deg){
    
    lands_new[[t+1]] <- lands_new[[t]][!rownames(lands_new[[t]]) %in% as.character(which_patch_removed_per_year[[t]]), 
                                       !colnames(lands_new[[t]]) %in% paste("X", which_patch_removed_per_year[[t]], sep = "")]
    patchID_new[[t+1]] <-  lands_new[[t+1]]$patchID
  }
  
  for(t in 2:(1+ntime_deg)){
    
    dmat_new[[t]] <- as.matrix(lands_new[[t]][, -c(1:4)])
    
    S_new[[t]] <- as.matrix(S_new[[t-1]][!rownames(S_new[[t-1]]) %in% as.character(which_patch_removed_per_year[[t-1]]),])
  }
  
  # # checks
  # duplicated(to_remove)
  # str(lands_new)
  # str(patchID_new)
  # 
  # str(dmat_new)
  # str(S_new)
  
  
  # plot(1:npatch, S_new[[1]], type = 'l', ylim = c(0,1))
  # for(i in 2:length(S_new)){
  #   lines(S_new[[i]], col = i)
  # }
  
  ################################
  ## RESTORATION  ----
  ################################
  
  # how many years of restoration is ntime_rest
  ntime_rest
  
  # years of restoration
  rest_years <- (ntime_first + ntime_deg + 1):(ntime_first + ntime_deg + ntime_rest)
  
  # we choose which patch to restore randomly during restoration years
  # but we try to restore the same patch IDs we removed before
  # still trying to keep the rate the same
  # so if we removed X number of patches per year, we will restore X number of patches per year
  # so our deg_percent = rest_percent
  # so here essentially nr_patch_removed_per_year = nr_patch_restored_per_year
  # but we write the code anyway for reproducibility
  # below line calculates this number based on our input
  nr_patch_restored_per_year <- nr_patch_removed_per_year
  
  
  to_remove # all of the removed patch IDs
  to_restore <- rep(NA, length = ntime_rest)
  # since we have all removed patch IDs in the to_remove vector
  # we can sample from that to randomly decide which removed patch should be restored
  # because we do not want to restore in the same order we destroyed
  to_restore <- sample(to_remove, size = nr_patch_restored_per_year*ntime_rest, replace = FALSE)
  # we can keep the nr_patch_removed_per_year same, because we are trying to keep the rate the same for degradation and restoration
  to_restore
  
  # save the removed patch IDs in each year
  which_patch_restored_per_year <- split(to_restore, ceiling(seq_along(to_restore) / nr_patch_restored_per_year))
  names(which_patch_restored_per_year) <- rest_years
  which_patch_restored_per_year
  
  #### fill the lists for restoration phase ####
  ##############################################
  
  for(t in (ntime_first + ntime_deg):((ntime_first + ntime_deg + ntime_rest)-1)){
    
    # get the landscape from the previous time step
    last_land <- lands_new[[as.character(t)]]
    
    # add columns 
    colbind <- as.matrix(initial_land[last_land$patchID,(which_patch_restored_per_year[[as.character(t+1)]] + 4)])
    colnames(colbind) <- colnames(initial_land[which_patch_restored_per_year[[as.character(t+1)]] + 4])
    
    temp_land <- cbind(last_land, colbind)
    
    # order col names
    fixed_cols <- colnames(temp_land)[c(1:4)]
    
    to_order_cols <- colnames(temp_land)[-c(1:4)]
    
    sorted_cols <- to_order_cols[order(as.numeric(gsub("X", "", to_order_cols)))]
    
    new_colnames <- c(fixed_cols, sorted_cols)
    
    temp_land <- temp_land[, new_colnames]
    
    # add rows
    rowbind <- initial_land[which_patch_restored_per_year[[as.character(t+1)]], ]
    rowbind <- rowbind[,which(colnames(rowbind) %in% colnames(temp_land))]
    
    temp_land <- rbind(temp_land, rowbind)  
    
    # order row names
    temp_land <- temp_land[order(as.numeric(temp_land$patchID)),]
    
    # final outcomes after restoration
    lands_new[[as.character(t+1)]] <- temp_land
    
    patchID_new[[as.character(t+1)]] <-  lands_new[[as.character(t+1)]]$patchID
    
  }
  
  for(t in (ntime_first + ntime_deg):(ntime_first + ntime_deg + ntime_rest)){
    
    dmat_new[[as.character(t)]] <- as.matrix(lands_new[[as.character(t)]][, -c(1:4)])
    
    S_new[[as.character(t)]] <- as.matrix(initial_S_mat[rownames(initial_S_mat) %in% rownames(lands_new[[as.character(t)]]),])
    
  }
  
  
  ####################################################
  ## AFTER YEARS (no degradation or restoration) ----
  ####################################################
  
  after_years <- (ntime_first + ntime_deg + ntime_rest + 1):(ntime_first + ntime_deg + ntime_rest + ntime_stop)
  
  after_npatch <- nrow(lands_new[[as.character(ntime_first + ntime_deg + ntime_rest)]])
  after_dmat <- dmat_new[[as.character(ntime_first + ntime_deg + ntime_rest)]]
  after_S <- c(S_new[[as.character(ntime_first + ntime_deg + ntime_rest)]])
  
  #####################################
  ## TRANSITION MATRIX - M ----
  #####################################
  
  # npatch x npatch matrix
  # diagonal elements = 0
  # non-diagonal elements = S[i] * S[j] * exp(-alpha * dmat[i,j]) * (1 - ifelse(IND[i]==IND[j],1,0)) 
  # transition matrix = M
  
  ################################
  ## INITIAL LANDSCAPE  ----
  ################################
  
  M_init = matrix(NA, ncol = nrow(initial_land), nrow = nrow(initial_land))
  for (i in 1:nrow(initial_land)) {
    for (j in 1:nrow(initial_land)) {
      if (i == j) { # diagnoal elements
        M_init[i,j] = 0
      } else { # non-diagonal elements
        M_init[i,j] = initial_S[i] * initial_S[j] * exp(-alpha * initial_dmat[i,j])
      }
    }
  }
  
  rownames(M_init) <- rownames(initial_dmat)
  colnames(M_init) <- colnames(initial_dmat)
  
  # let's put them in a list together
  # deg_years # degradation years
  # rest_years # restoration years
  # after_years # after stopping
  
  ###########################################
  ## DEGRADATION  & RESTORATION TOGETHER ----
  ### for lambdaM calculation
  ###########################################
  
  # + 1 is for the initial land and then the degradation and then the restoration
  dimensions <- vector(mode = "list", length = 1 + ntime_deg + ntime_rest)
  
  for(t in 1:length(dimensions)){
    dimensions[[t]] <- c(nrow(lands_new[[t]]), nrow(lands_new[[t]]))
  }
  
  M_all <- lapply(dimensions, function(dim) matrix(NA, nrow = dim[1], ncol = dim[2]))
  names(M_all) <- (ntime_first):(ntime_first + ntime_deg + ntime_rest)
  M_all[[as.character(ntime_first)]] <- M_init
  
  # during degradation years
  
  for(t in deg_years){
    M_all[[as.character(t)]] <- outer(S_new[[as.character(t)]][,1], S_new[[as.character(t)]][,1]) * exp(-alpha * dmat_new[[as.character(t)]])
    diag(M_all[[as.character(t)]]) <- 0 
  }
  
  # during restoration years
  
  for(t in rest_years){
    M_all[[as.character(t)]] <- outer(S_new[[as.character(t)]][,1], S_new[[as.character(t)]][,1]) * exp(-alpha * dmat_new[[as.character(t)]])
    diag(M_all[[as.character(t)]]) <- 0 
  }
  
  rownames(M_all[[as.character(ntime_first)]]) <- rownames(dmat_new[[as.character(ntime_first)]])
  colnames(M_all[[as.character(ntime_first)]]) <- colnames(dmat_new[[as.character(ntime_first)]])
  
  for(t in deg_years){
    rownames(M_all[[as.character(t)]]) <- rownames(dmat_new[[as.character(t)]])
    colnames(M_all[[as.character(t)]]) <- colnames(dmat_new[[as.character(t)]])
    
  }
  
  for(t in rest_years){
    rownames(M_all[[as.character(t)]]) <- rownames(dmat_new[[as.character(t)]])
    colnames(M_all[[as.character(t)]]) <- colnames(dmat_new[[as.character(t)]])
    
  }
  
  # lambdaM is the metapopulation capacity
  # it is the leading eigenvalue of the transition matrix
  initial_lambdaM <- eigen(M_init)$values[1] # same with max(M_eigen$values)
  
  #lambdaM <- c(rep(NA, length = ntime))
  lambdaM <- vector(mode = "numeric", length = ntime)
  lambdaM[1:ntime_first] <- initial_lambdaM
  
  for(t in (ntime_first+1):(ntime_first+ntime_deg+ntime_rest)){
    lambdaM[t] <- Re(eigen(M_all[[as.character(t)]])$values[1])
  }
  
  lambdaM[max(rest_years):ntime] <- lambdaM[max(rest_years)]
  
  lambdaM_df <- NA
  lambdaM_df <- data.frame(lambdaM = lambdaM, landscapeID = land)
  
 
  rm(M_init, M_all,  initial_A, initial_lambdaM)
  
  ## Calculate probabilities ----
  ###################################################
  
  for(e in e.values){
    
    print(paste("Running e value", e, ":"))
    
    for(c in c.values){
      
      print(paste("Running c value", c, ":"))
      
      ## initial occupancy
      ## use pstar function as starting values
      init_psi1_mean <- pstar.function(landscape=initial_land, alpha, delta=e/c, iterations=1000, S_int, S_b1)
      
      
      out <- sim_obsDat_landrest_100425(nsim = nsim, npatch = npatch, ntime = ntime,
                                        ntime_first = ntime_first, ntime_deg = ntime_deg, ntime_rest = ntime_rest, ntime_stop = ntime_stop,
                                        lands_new = lands_new, dmat_new = dmat_new, S_new = S_new,
                                        init_psi1_mean = init_psi1_mean, alpha = alpha, c = c, e = e, x = x)
      
      ######## save output before simulating the next metapopulation ######## 
      out_list <- list(NA)
      prob_list <- list(NA)
      
      for(i in 1:nsim){
        out_list[[i]] <- out[[i]][[1]]
        prob_list[[i]] <- out[[i]][[2]]
        
      }
      
      my.results <- bind_rows(out_list, .id = "metapopID") # from dplyr package
      my.results$landscapeID <- land
      my.results$lambdaM <- round(lambdaM, digits = 2)
      my.results$deg_percent <- deg_percent
      
    
      # write to table results:
      write.table(my.results, file = filename, append=T, sep=",", col.names = F, row.names = F)
      
      # remove for the next repetition
      rm(out_list, prob_list, my.results, out); gc()
      
    } # end c.values
  } # end e.values
  
  # write to table landscape related results - lambdaM:
  write.table(lambdaM_df, file = filename2, append=T, sep=",", col.names = F, row.names = F)
  
  # remove current landscape before creating new
  rm(sim_land, initial_land, initial_dmat, initial_S,  dmat_new, S_new, 
     patchID_new, to_remove, which_patch_removed_per_year, to_restore, which_patch_restored_per_year,
     lambdaM, lambdaM_df);gc()
  
} # end nlandscapes

end.time <- Sys.time()

print(end.time-start.time)
