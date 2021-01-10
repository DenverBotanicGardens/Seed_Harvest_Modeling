#' @title Seed Harvest Simulation
#' @description Simulation of climate and seed harvest frequency and intensity
#' 
#' @export
#' 
#' @param Mx_list A list of annual transition matrices indicating survival and growth from col to row; fecundity in the first row (square matrix with N stages).
#' @param TMx_list A list of annual transition matrices without fecundity (square matrix with N stages).
#' @param Nx A population size vector by stage (length N). Stable stage distribution or most recent count by stage. 
#' @param Nx_list A list of annual size by stage vectors
#' @param StartPopSize A number for the starting population size of the simulation. 
#' @param TotYrs A number indicating how many years to project in the simulation.
#' 
#' @return 

# Test
Nx_list <- list(Nxs_annual)
StartPopSize <- 100
TotYrs <- 100
Nx_start <- stable.stage(mean(Mx_sample))
# rm(Nx_list); rm(StPopSz);rm(StartPopSize);rm(projlength);rm(TotYrs);rm(PopSize);rm(nx);rm(mats);rm(Mx_sample);rm(Nx_scale);rm(vec)

# starting population of only above ground individuals - pop size is above ground
SeedHarvestSimple <- function(Mx_list, TMx_list, Nx, Nx_list, StartPopSize,TotYrs){
  # take the distribution of classes (the stable stage distribution) and scale to starting population size with integers
  # Scale to starting population size with only the aboveground, detectable stages (no seedbank or dormant); undetectable stages have NA in Nx_list
  Nx_scale <- Nx[!is.na(Nx_list[[1]])]
  # vec <- matrix(floor(Nx*(StartPopSize/sum(Nx))), ncol = 1)
  vec <- matrix(round(Nx*(StartPopSize/sum(Nx_scale)),0), ncol = 1)
  
  # Initialize 
  # PopSize <- rep(NA, TotYrs)
  # lambdas <- rep(NA, TotYrs)
  PopSize <- c()
  lambdas <- c()
  Extant <- c()
  mats <- list()
  Tulapprox <- c()
  LogGrowthRate <- c()
  Yr <- c()
  
  ######################## Annual Transitions without clusters ######################
  for(i in 1:TotYrs){
    nx <- sample(Mx_list, 1)
    # Multiply N_t+1 = Matrix*N_t
    vec <- floor(nx[[1]]%*%vec)
    if(sum(vec)<1) break
    
    Yr[i] <- i
    PopSize[i] <- sum(vec)
    mats[[i]] <- nx[[1]]
    lambdas[i] <- lambda(nx[[1]])
    if(sum(vec)<1) Extant[i] <- 0 else Extant[i] <- 1
    Tulapprox[i] <- stoch.growth.rate(mats[[1]], verbose=FALSE)$approx
    # if the matrix is singular (the determinant is zero) then $sim will be NaN
    # det(mats[[1]])
    LogGrowthRate[i] <- stoch.growth.rate(mats[[1]], verbose=FALSE)$sim
  }
  if(any(Extant==0)) t2e <- min(which(Extant==0)) else t2e <- NA 
  data.frame(PopSize = PopSize, Year = Yr, lambdas = lambdas, StPopSize = StartPopSize, 
             Time2Ext = t2e, Tulapprox = Tulapprox, LogGrowthRate = LogGrowthRate) # $approx is log stochastic growth rate by Tuljapukar's approximation; $sim and CI
}
