############################################################################
# Observation data simulation for the theoretical models
# In this script we create a function that can simulate observation data 
# based on formulations adapted from:
# Sutherland et al., 2014, Ecology, doi: 10.1890/14-0384.1
# AND 
# Bertasello et al., 2021, Royal Soc. Open Sci. doi: 10.1098/rsos.201309
#
# UPDATE: separating the landscape simulation and species simulation
# so that we can simulate different species in the same landscape
# Thus, removing landscape function outside.
#
# ### #####################################################
# ### This function involves the land degradation part ####
# ### #####################################################
#
# Date: 07/11/2024
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data ----
###########################################################################

# rm(list = ls())

## 2. Observation Data Simulation Function ----
##########################################################

# TO RUN THIS FUNCTION YOU NEED THE LANDSCAPE CHARACTERISTICS ALREADY
# SO RUN THE ~/landscape_simulation_functions.R FIRST THEN THIS

###########################################################################
# INPUT:
# nsim = how many metapopulations
# npatch = number of patches in the initial landscape
# ntime = total nummber of time (years) spanning entire simulation - should be equal to sum of ntime_first, ntime_deg, ntime_stop
# 
# ntime_first = length of equilibrium years
# ntime_deg = length of degradation years
# ntime_stop = length of rest years
#
# lands_new = list of different landscapes - as we degrade the initial landscape
# dmat_new = list of distance matrices for different landscapes - as we degrade the initial landscape
## contains one dmat for equilibrium timespan, all of degradation timespan, last dmat in the degradation timespan is considered as dmat for rest timespan
# S_new = list of habitat suitability for different landscapes - as we degrade the initial landscape
# 
# init_psi1_mean = first year occupancy probabilitt based on initial landscape
# alpha = species-specific dispersal rate
# c = species-specific colonization scaling parameter
# e = species-specific extinction scaling parameter
# x = indicator for the steepness of the relationship vetween e and suitability
###########################################################################


sim_obsDat_landdeg_100425 <- function(nsim, npatch, ntime,
                                      ntime_first, ntime_deg, ntime_stop,
                                      lands_new, dmat_new, S_new, # landscape characteristics
                                      init_psi1_mean, alpha, c, e, x){
  
  
  ## collect results ----
  ################################################################
  
  # what outputs we want to save
  # outputs_to_save <- c("year", "metapopID", "npatch", "e", "c", "alpha", "x", "e_over_c", "lambdaM", 
  #                      "Persistence.pred", "Persistence.obs", "t.extinct",
  #                      "n.occ", "avg.psi", "avg.gamma", "avg.eps") # removed "psi.fs"# removed "psi.fs"
  # 
  # vector_outputs_to_save <- c("year", "metapopID", "e", "c", "alpha", 
  #                             "npatch", "n.occ", "n.occ.prop", "avg.psi", "avg.gamma", "avg.eps")
  
  ###########################################################################
  # PARAMETERS:
  # psi1 = initial occupancy probability
  # psi = occupancy probability
  # muZ = true occupancy state placeholder
  # z = true occupancy state
  # eps = extinction probability
  # E = extinction probability placeholder
  # e = extinction probability scaling parameter (species-specific)
  # S = patch suitability (can be logit-function of patch characteristics and change through time)
  # con = patch connectivity
  # C = connectivity
  # gamma = colonization probability
  # c = colonization probability scaling parameter (species-specific)
  ###########################################################################
  
  ## provide empty holders  ----
  ################################################################
  
  vector_sim_out <- vector("list", nsim)
  
  probs_outputs_to_save <- c("psi", "gamma", "eps")
  
  probs_sim_out <- vector("list", nsim)
  
  sim_out <- list(NA)
  
  muZ <- z <- array(dim = c(npatch, ntime)) # Occupancy, occurrence
  con <- array(dim = c(npatch, npatch, ntime)) # Connectivity
  eps <- gamma <- E <- C <- array(NA, dim = c(npatch, (ntime-1))) # Extinction
  
  ## provide out-of-loop fixed parametes ----
  ################################################################
  
  ### initial equilibrium state years ----
  ################################################################
  init_npatch <- nrow(dmat_new[[1]])
  init_dmat <- dmat_new[[1]]
  init_S <- c(S_new[[1]])
  
  ### degradation years ----
  ################################################################
  deg_years <- (ntime_first + 1):(ntime_first + ntime_deg)
  
  ### after years ----
  ################################################################
  after_years <- (ntime_first + ntime_deg + 1):(ntime_first + ntime_deg + ntime_stop)
  
  after_npatch <- nrow(lands_new[[as.character(ntime_first + ntime_deg)]])
  after_dmat <- dmat_new[[as.character(ntime_first + ntime_deg)]]
  after_S <- c(S_new[[as.character(ntime_first + ntime_deg)]])
  
  ## simulation start  ----
  ################################################################
  
  library(lubridate)
  
  print("Running simulation function:")
  
  # progress_bar = txtProgressBar(min=0, max=nsim, style = 3, width = nsim, char="=")
  
  # init <- numeric(nsim)
  # end <- numeric(nsim)
  
  
  for(d in 1:nsim){
    
    # init[d] <- Sys.time()
    
    
    #############################
    ##### EQUILIBRIUM YEARS #####
    #############################
    
    ## Ecological process  ----
    ##################################
    
    # initial occupancy, no covariates
    # psi[,1] <- init_psi1_mean
    muZ[,1] <- init_psi1_mean
    
    # z[,1] <- rbinom(init_npatch, 1, psi[,1])
    z[,1] <- rbinom(init_npatch, 1, muZ[,1])
    
    # extinction probability, changes with years and patch area
    
    for(i in 1:init_npatch){ # so far suitability does not change with time, only patches
      
      # S[i] <- plogis(S_int + S_b1*A[i]) # calculated outside
      
      for(t in 2:ntime_first){
        E[i,t-1] <- e/(init_S[i])^x
        eps[i,t-1] <- 1-exp(-E[i,t-1])
      }
    }
    
    # connectivity
    for(i in 1:init_npatch){
      for(j in 1:init_npatch){
        for(t in 1:ntime_first){
          
          con[i,j,t] <- exp(-alpha * init_dmat[i,j]) * init_S[j]
        }
      }
    }
    
    # colonization probability
    for(t in 2:ntime_first){
      for(i in 1:init_npatch){
        C[i,t-1] <- c * (z[1:init_npatch,t-1] %*% con[i,1:init_npatch,t-1])                           # measure of connectivity between patches
        
        gamma[i,t-1] <- 1 - exp(-C[i,t-1])                                               # colonization rate
        
        muZ[i,t] <- z[i,t-1]*(1-eps[i,t-1]) + (1-z[i,t-1])*gamma[i,t-1]
        
        # True state z
        z[i,t] <- rbinom(1, 1, muZ[i,t])
        
      }
    }
    
    #############################
    ##### DEGRADATION YEARS #####
    #############################
    
    deg_years <- (ntime_first + 1):(ntime_first + ntime_deg)
    
    ## Ecological process  ----
    ##################################
    
    # # extinction probability, changes with years and patch area
    
    for (t in deg_years) {
      S_vec <- S_new[[as.character(t)]][,1]  # Extract vector of suitability values
      n <- length(S_vec)
      
      # Extinction rate and probability for all patches
      E[1:n, t-1] <- e / (S_vec^x)
      eps[1:n, t-1] <- 1 - exp(-E[1:n, t-1])
      
      # Vectorized connectivity matrix
      con[1:n, 1:n, t] <- exp(-alpha * dmat_new[[as.character(t)]]) * matrix(rep(S_vec, each = n), nrow = n)
    }
    
    
    # colonization probability
    for(t in deg_years){
      for(i in 1:nrow(lands_new[[as.character(t)]])){ 
        C[i,t-1] <- c * (z[1:nrow(lands_new[[as.character(t)]]),t-1] %*% con[i,1:nrow(lands_new[[as.character(t)]]),t-1]) # measure of connectivity between patches
        
        gamma[i,t-1] <- 1 - exp(-C[i,t-1])                                               # colonization rate
        
        muZ[i,t] <- z[i,t-1]*(1-eps[i,t-1]) + (1-z[i,t-1])*gamma[i,t-1]
        
        # True state z
        z[i,t] <- rbinom(1, 1, muZ[i,t])
        
      }
    }
    
    ###########################
    ##### AFTER DEG YEARS #####
    ###########################
    
    after_years <- (ntime_first + ntime_deg + 1):(ntime_first + ntime_deg + ntime_stop)
    
    after_npatch <- nrow(lands_new[[as.character(ntime_first + ntime_deg)]])
    after_dmat <- dmat_new[[as.character(ntime_first + ntime_deg)]]
    after_S <- c(S_new[[as.character(ntime_first + ntime_deg)]])
    
    ## Ecological process  ----
    ##################################
    
    # extinction probability, changes with years and patch area
    
    for(t in after_years){
      for(i in 1:after_npatch){
        
        E[i,t-1] <- e/(after_S[i])^x
        eps[i,t-1] <- 1-exp(-E[i,t-1])
        
        # connectivity
        for(j in 1:after_npatch){
          
          con[i,j,t] <- exp(-alpha * after_dmat[i,j]) * after_S[j]
        }
      }
    }
    
    # colonization probability
    for(t in after_years){
      for(i in 1:after_npatch){
        
        # measure of connectivity between patches
        C[i,t-1] <- c * (z[1:after_npatch,t-1] %*% con[i,1:after_npatch,t-1])
        
        # colonization rate
        gamma[i,t-1] <- 1 - exp(-C[i,t-1])
        
        muZ[i,t] <- z[i,t-1]*(1-eps[i,t-1]) + (1-z[i,t-1])*gamma[i,t-1]
        
        # True state z
        z[i,t] <- rbinom(1, 1, muZ[i,t])
        
      }
    }
    
    ## Transition matrix and lambda calculated outside
    ###################################################
    
    ## Derived parameters ----
    ###################################
    
    # Compute annual population occupancy
    # same with muZ, thus not doing it again (11/11/2024)
    # for(i in 1:npatch){
    #   for (t in 2:ntime){
    #     psi[i,t] <- psi[i,t-1]*(1-eps[i,t-1]) + (1-psi[i,t-1])*gamma[i,t-1]
    #   }
    # }
    # psi <- muZ
    
    # Number of occupied sites and proportion of occupied sites
    tot.sites <-  n.occ <- n.occ.prop <- NA
    
    tot.sites[1] <- length(which(!is.na(z[1:npatch,1])))
    n.occ[1] <- sum(z[1:npatch,1], na.rm = T)
    n.occ.prop[1] = n.occ[1]/tot.sites[1]
    
    for(t in 2:ntime){
      tot.sites[t] <- length(which(!is.na(z[1:npatch,t])))
      
      n.occ[t] <- sum(z[1:npatch,t], na.rm = T) # true number of occupied sites
      
      n.occ.prop[t] = n.occ[t]/tot.sites[t]
      
    }
    
    # Compute annual average of psi, gamma and eps
    avg.psi <- round(colMeans(muZ, na.rm = T),3)
    avg.gamma <- round(colMeans(gamma, na.rm = T),3)
    avg.eps <- round(colMeans(eps, na.rm = T),3)
    
    ## Take out derived parameters  ----
    ###########################################
    
    ###########################################################################
    # OUTPUT:
    # year = 1:ntime
    # c = species-specific colonization scaling parameter
    # e = species-specific extinction scaling parameter
    # alpha = species-specific dispersal rate
    # npatch = number of patches in the landscape over the years - as we degrade
    # n.occ = number of occupied sites
    # n.occ.prop = proportion of occupied sites
    ###########################################################################
    # psi = occupancy probability
    # gamma = colonization probability
    # eps = extinction probability
    ###########################################################################
    
    
    # if(n.occ[ntime]>0){
    #   t.extinct <- NA
    # } else {
    #   t.extinct <- min(which(n.occ==0))
    # } 
    
    vector.data.out <- data.frame(year = 1:ntime, 
                                  metapopID = d,
                                  e = e, 
                                  c = c,
                                  alpha = alpha, 
                                  npatch = tot.sites,
                                  n.occ = n.occ,
                                  n.occ.prop = round(n.occ.prop, digits = 2),
                                  avg.psi, 
                                  avg.gamma=c(NA, avg.gamma), 
                                  avg.eps=c(NA, avg.eps))
    
    psi_df <- as.data.frame(muZ)
    psi_df$e <- e
    psi_df$c <- c
    psi_df$metapopID <- d
    
    gamma_df <- as.data.frame(gamma)
    gamma_df$e <- e
    gamma_df$c <- c
    gamma_df$metapopID <- d
    
    eps_df <- as.data.frame(eps)
    eps_df$e <- e
    eps_df$c <- c
    eps_df$metapopID <- d    
    
    probs.data.out <- list(psi_df, gamma_df, eps_df)
    names(probs.data.out) <- c("psi", "gamma", "eps")
    
    # data.out <- as.matrix(data.out) # store in an array instead (TRUE=1, FALSE=0)
    
    # sim_out[,,d] <- data.out
    
    vector_sim_out[[d]] <- vector.data.out
    probs_sim_out[[d]] <- probs.data.out
    sim_out[[d]] <- list(vector_sim_out[[d]], probs_sim_out[[d]])
    
    rm(vector.data.out, probs.data.out)
    
    # end[d] <- Sys.time()
    
    # ################# PROGRESS BAR SETTINGS #################
    # setTxtProgressBar(progress_bar, value = d)
    # time <- round(seconds_to_period(sum(end - init)), 0)
    # 
    # # Estimated remaining time based on the
    # # mean time that took to run the previous iterations
    # est <- nsim * (mean(end[end != 0] - init[init != 0])) - time
    # remainining <- round(seconds_to_period(est), 0)
    # 
    # cat(paste(" // Execution time:", time,
    #           " // Estimated time remaining:", remainining), "")
    # ################# END PROGRESS BAR SETTINGS #################
    
  }
  
  # close(progress_bar)
  
  return(sim_out)
}

###########################################################################
# trying the function
###########################################################################

# nmetapops <- 3
# nsim <- nmetapops
# npatch <- 600 
# ntime <- 300
# ntime_first <- 100 
# ntime_deg <- 100
# ntime_stop <- 100
# 
# # need to run landscape function and init psi function for these
# lands_new = lands_new
# dmat_new = dmat_new
# S_new = S_new
# init_psi1_mean = init_psi1_mean
# 
# 
# alpha <- 1/500 
# e <- 0.2
# c <- 0.4
# x <- 1
# 
# start <- Sys.time()
# 
# out2 <- sim_obsDat_landdeg_100425(nsim = nsim, npatch = npatch, ntime = ntime,
#                                   ntime_first = ntime_first, ntime_deg = ntime_deg, ntime_stop = ntime_stop,
#                                   lands_new = lands_new, dmat_new = dmat_new, S_new = S_new,
#                                   init_psi1_mean = init_psi1_mean, alpha = alpha, c = c, e = e, x = x)
# 
# end <- Sys.time()
# 
# end - start
# 
# str(out2)

