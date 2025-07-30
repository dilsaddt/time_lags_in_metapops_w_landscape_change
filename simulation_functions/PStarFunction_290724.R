############################################################################
# Observation data simulation for the theoretical models
# In this script we create a function that 
# calculates the probability individual patches are occupied at equilibrium
# P* function
# Adapted from PStar Function.R from Walker & Gilbert 2023
#
# Date: 29/07/2024
# Author: Dilsad Dagtekin
###########################################################################

###########################################################################
# 1. House keeping and loading libraries and data ----
###########################################################################

# rm(list = ls())

## 2. P* Function ----
##########################################################
# calculates the probability of individual patches are occupied at equilibrium

###########################################################################
# INPUT:
# n.patch = number of patches
# landscape = a dataframe created by the create landscape function which specifies patch locations, areas, 
#  and an interpatch distance matrix (dmat - yes, checked it, its dmat)
# alpha = 1/(average dispersal distance) for a species
# delta = the ratio of the extinction to colonization rate of a species  (HM e/c  ?)
# iterations = the number of iterations for which the iterative function f is iterated, 
## greater iterations = greater accuracy, less iterations = lower accuracy

# NOTE: delta can be factored out of this calculation by letting delta = 1
###########################################################################

pstar.function<-function(landscape, alpha, delta, iterations, S_int, S_b1){
  
  n.patch <- length(landscape$patchID) ## here patch.ID, in our sim_land$patchID
  A <- landscape$A
  d <- as.matrix(landscape[,5:(dim(landscape)[[2]])])
  
  # standardize log Areas:
  # A <- scale(log(landscape$A))[,1]
  ## we changes this to scale(area), dropped the log:
  A <- scale(landscape$A)[,1]
  
  # calculate suitability S:
  S <-  plogis(S_int + S_b1*A) 
  
  # STEP 1: SET UP OF MATRIX M
  M <- rep(NA, n.patch*n.patch)
  dim(M) <- c(n.patch,n.patch) 
  
  for (i in 1:n.patch){ #for each row of matrix M
    for (j in 1:n.patch){ #for each column of matrix M
      if (i==j){
        M[i,j] <- 0} #let the connectivity of a patch to itself be 0
      else{
        M[i,j] <- S[i]*S[j]*exp(-alpha*d[i,j])} 
      
      # let the connectivity of other patches to that patch be AiAje^(-adij)Pj
      
      ### is this necessary?:
      #if(M[i,j]<1*(10^-8)){M[i,j]<-1*(10^-8)}#assume that there is always some very very small
      ##probability that a patch may be colonized because there is never 0 chance and because
      ##R makes a rounding error making this 0, if this is too close to 0
      
    }}
  
  # ITERATING FUNCTION TO FIND P*
  p <- rep(NA, n.patch*iterations)
  dim(p) <- c(iterations,n.patch)
  
  p[1,1:n.patch] <- rep(0.1,n.patch)
  
  for (t in 1:(iterations-1)){ #iterate the following for a given number of iterations
    P <- p[t,] #use p of previous iteration
    g <- (M%*%P) #SRLM g function from Ovaskainen and Hanski, 2001
    f <- (g/(g+delta)) #SRLM f function from Ovaskainen and Hanski, 2001, 
    #which when iterated many times will give p*
    p[(t+1),]<-t(f)} 
  #set p for the next iteration = to the output value of p for this iteration
  p.star<-p[iterations,] #p* = the final iterations value of p after many iterations
  ## HM: expected equilibrium value for each patch i?
  return(p.star)
}
