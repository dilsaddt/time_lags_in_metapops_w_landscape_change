############################################################################
# Landscape simulation for the theoretical models
# In this script we create a function that can can simulate landscapes
# Code adapted from Walker & Gilbert, 2023, Ecology, doi: 10.1002/ecy.3840
#
# Landscape simulation:
# instead of habitat quality (capacity) metric is K, we used area
#
# We create either randomly distributed or clustered landscapes
#
# Date: 29/07/2024
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data ----
###########################################################################

# rm(list = ls())

library(EnvStats)

## 2. Landscape Function ----
#############################

###########################################################################
# Input:
# *npatch* = number of patches in desired landscape
# *landscapeLimit* = size of desired landscape (max possible x and y coordinate of a patch)
# *landscapeType* = specify whether the desired landscape's patches should be randomly distributed or clustered
# accepts either "clustered" or "random"
# *niter* = number of iterations for which the clustering algorithm is run for, specifying the level of clustering
# higher the value, the more patches are clustered or regularly spaced
# if it is a "clustered" landscape versus a "random" landscape
###########################################################################

createLandscape<-function(npatch, landscapeLimit, landscapeType, niter){
  
  npatch <- npatch
  
  # assign ID numbers to each patch (This is so that we can track which patches are
  # removed throughout destruction of patches from the landscape)
  patchID <- c(1:npatch)
  
  ################
  ## PATCH AREA ##
  ################
  
  # A <- 100 # for starters all areas are equal, 100 m2
  
  # log normal distrubted areas
  # A <- rlnorm(npatch, meanlog = 4, sdlog=1) # mean ~ log(100)
  
  # truncated log normal, min area = 10 m2
  A <- rlnormTrunc(npatch, meanlog = log(150), sdlog = 1, min = 10, max = 1000)

  # A <- rlnorm(npatch, meanlog = 2, sdlog=1 )
  # patch areas lognormally distributed with a mean of 2 and standard deviation of 1
  radii <- sqrt(A/pi) # calculate radii of hypothetically circular patches
  
  radii <- rep(radii, times = npatch) # we need a vector for that 
  
  #################
  ## COORDINATES ##
  #################
  
  # pick a random number for the x and y coordinate of each patch between 0 and the extent of the landscape
  xCoord <- runif(npatch, min = 0, max = landscapeLimit) 
  yCoord <- runif(npatch, min = 0, max = landscapeLimit) 
  
  # the landscape of patch locations is given by the randomly chosen x and y coordinates paired together
  coordinates <- data.frame(xCoord, yCoord) 
  
  #####################
  ## DISTANCE MATRIX ##
  #####################
  
  # Create distance matrix of patches - dmat
  # calculate distances between patches based on euclidean distances
  
  dmat  <- as.matrix(dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE))
  # dii=0 and dij=dji, which is indicated by diag=TRUE and upper=TRUE
  
  # LOOP TO PREVENT PATCH OVERLAP
  for (i in 1:npatch){ # for every patch
    for (j in 1:npatch){ # and every other patch
      while (j!=i & dmat[i,j] < (radii[i]+radii[j])) { 
        # while the distance between the two patches is less than the sum of their radii
        xCoord[i] <- runif(1, min = 0, max = landscapeLimit) 
        # pick a new x coordinate for that patch
        yCoord[i] <- runif(1, min = 0, max = landscapeLimit) 
        # pick a new y coordinate for that patch
        coordinates[i,]<-c(xCoord[i], yCoord[i]) # update the coordinates 
        
        # update the distances
        dmat <- as.matrix(dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE))
      }
    }
  }
  
  ####################
  ## LANDSCAPE TYPE ##
  ####################
  
  if (landscapeType == "clustered"){ 
    # If we want a non-random landscape (clustered)...
    for (i in 1:niter){
      
      # STEP 1: PICK A RANDOM PATCH UP FROM THE LANDSCAPE
      ######################################################################################
      
      r <- sample(1:npatch, 1, replace=T) # pick a random number between 1 and npatch
      rCoord <- coordinates[-r,] # remove that patch from the landscape (pick it up to be potentially relocated)
      
      # STEP 2: CALCULATE THE DISTANCE AND CONNECTIVITY OF THE PATCH CHOSEN TO OTHER PATCHES
      ######################################################################################
      
      d_r_to_j <- rep(NA, (npatch-1))
      con_r_to_j <- rep(NA, (npatch-1))
      for (j in 1:(npatch-1)) { # for every patch in the network calculate...
        d_r_to_j[j] <- sqrt((xCoord[r]-rCoord[j,1])^2+(yCoord[r]-rCoord[j,2])^2) # the distance between point r and all patches
        con_r_to_j[j] <- exp(-d_r_to_j[j]) # the connectivity of point r and all others
      }
      
      con_r <- sum(con_r_to_j) # the connectivity of point r is the sum of it's connectivity to all patches
      
      # STEP 3: PICK A RANDOM POINT
      ######################################################################################
      
      p_xCoord <- runif(1, min=0, max=landscapeLimit) 
      # pick a random x coordinate for the point between 0 and the landscape extent
      p_yCoord <- runif(1, min=0, max=landscapeLimit) 
      # pick a random y coordinate for the point between 0 and the landscape extent
      
      # STEP 4: CALCULATE THE DISTANCE AND CONNECTIVITY OF THE NEW POINT TO OTHER PATCHES
      ######################################################################################
      
      d_x_to_j<-rep(NA, (npatch-1))
      con_x_to_j<-rep(NA, (npatch-1))
      for (j in 1:(npatch-1)) { # for every patch in the network calculate...
        d_x_to_j[j] <- sqrt((p_xCoord-rCoord[j,1])^2+(p_yCoord-rCoord[j,2])^2) # the distance between point x and all patches
        con_x_to_j[j] <- exp(-d_x_to_j[j]) # the connectivity of point x and all others
      }
      
      con_x <- sum(con_x_to_j) # the connectivity of point x is the sum of it's connectivity to all patches
      
      # STEP 5: LOOP TO ENSURE NEW POINT WILL NOT RESULT IN PATCH OVERLAP
      ######################################################################################
      for (g in 1:(npatch-1)) { #for each patch
        while (g!=r & d_x_to_j[g]<(radii[r]+radii[g])){ 
          # while any patch other than the chosen patch r has a distance between it's 
          # radius and r's radius greater than the distance between that patches center 
          # and the proposed new location for r
          p_xCoord <- runif(1, min=0, max=landscapeLimit) 
          # pick a new x coordinate for the new point
          p_yCoord <- runif(1, min=0, max=landscapeLimit) 
          # pick a new y coordinate for the new point
          
          d_x_to_j<-rep(NA, (npatch-1))
          
          for (j in 1:(npatch-1)){ 
            # for every patch in the network update calculations for...
            d_x_to_j[j] <- sqrt((p_xCoord-rCoord[j,1])^2+(p_yCoord-rCoord[j,2])^2) # the distance between point x and all patches
            con_x_to_j[j] <- exp(-d_x_to_j[j]) # the connectivity of point x and all others
          }
          
          con_x <- sum(con_x_to_j) # the connectivity of point x is the sum of it's connectivity to all patches
          
        }
      }
      
      
      # STEP 6: COMPARE THE CONNECTIVITY OF THE PATCH AT ITS ORIGINAL LOCATION TO ITS 
      ######################################################################################
      # NEW LOCATION
      # to generate clustered landscape: increase the pr(r relocated) 
      # if con_x > con_r
      if (con_x > con_r){ 
        #if the new location has a higher connectivity...
        xCoord[r] <- p_xCoord # set the x coordinate of the patch to be the new location's
        yCoord[r] <- p_yCoord # set the y coordinate of the patch to be the new location's
        coordinates[r,]<-c(xCoord[r],yCoord[r]) # enter the patch's new coordinates into the landscape
      }
    }
  }
      
  # STEP 7: UPDATE MATRIX OF DISTANCES BETWEEN PATCHES
  ######################################################################################
  dmat <- as.matrix(dist(coordinates[,1:2], method="euclidean", diag=TRUE, upper=TRUE))
  # calculate distances between patches based on euclidean distances. 
  # dii=0 and dij=dji, which is indicated by diag=TRUE and upper=TRUE

  landscape <- data.frame(patchID, A, coordinates, dmat)
  return(landscape)
}

## 3. Multiple Landscape Function ----
######################################
# creates a multiple landscapes of all possible landscape types of a specified size with a specified number
# of patches, with a patch distribution of specified level of clustering

###########################################################################
# Input:
# *npatch* = number of patches in desired landscape
# *landscapeLimit* = size of desired landscape (max possible x and y coordinate of a patch)
# *landscapeType* = may be a list of landscape types to be created eg. "clustered", "random"
# *niter* = number of iterations for which the clustering algorithm is run for, specifying the level of clustering
# higher the value, the more patches are clustered or regularly spaced
# if it is a "clustered" landscape versus a "random" landscape
# *nrep* = number of landscapes of each type to create
# directory = where the output data should be saved as csv
###########################################################################

createMultipleLandscapes <- function(npatch, landscapeLimit, landscapeType, niter, nrep, directory){
  setwd(directory)
  for (j in 1:length(landscapeType)){
    for (i in 1:nrep){
      landscape_i_data <- createLandscape(npatch=npatch,
                                        landscapeLimit=landscapeLimit,
                                        niter=niter,
                                        landscapeType = landscapeType[j])
      
      landscape_i_data$landscapeNo <- i # add a column indicating the individual landscape's number
      
      # stitching all the individual landscapes into one data frame
      if (i==1){type_j_data <- landscape_i_data} 
      else {type_j_data <- rbind(type_j_data, landscape_i_data)}
    }
    type_j_data$landscapeType <- landscapeType[j] # add a column indicating the type of landscape's generated
    
    # stitching the data on all landscapes of each type together
    if (j==1){data <- type_j_data} 
    else {data <- rbind(data, type_j_data)}
  }
  # output file to be written with all the landscape's data and metadata
  file.name<-paste0(npatch, "Patch_", nrep, "xLandscapes", niter, "ClusterLevel", ".csv")
  write.csv(data, file=file.name)
  return(data)
}







