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
#' @param GoodBadTm A transition matrix (row to col) of good and bad years. Annual transition matrices are defined as bad < mean fecundity <= good 
#' @param Freq A numerical vector for the the probability of collecting seed given the year c(good, bad).
#' @param Int A numerical vector for the intensity (the percentage decrease in fecundity) if seed harvest for c(good, bad).
#' @param TotYrs A number indicating how many years to project in the simulation.
#' @param ClusteredColl A numerical value indicating the number of years of collection to be clustered. Results in increased likelihood of collection in t+1 when collecting in year t
#' @param ddceiling TRUE or FALSE to revert to the dominant eigenvalue with logarithmic term when within 90% of the ceiling, the ceiling set at 10*starting population size 
#' 
#' @return 


# starting population of only above ground individuals - pop size is above ground
SeedHarvestSim <- function(Mx_list, TMx_list, Nx, Nx_list, StartPopSize, GoodBadTm, Freq = c(0,0), Int = c(0,0), 
                           TotYrs, ClusteredColl = 1, ddceiling = FALSE){
  # take the distribution of classes (likely the stable stage distribution) and scale to starting population size with integers
  # Scale to starting population size with only the aboveground, detectable stages (no seedbank or dormant); undetectable stages have NA in Nx_list
  Nx_scale <- Nx[!is.na(Nx_list[[1]])]
  # vec <- matrix(floor(Nx*(StartPopSize/sum(Nx))), ncol = 1)
  vec <- matrix(round(Nx*(StartPopSize/sum(Nx_scale)),0), ncol = 1)
  # Annual transition function needs local environment variables
  environment(AnnualTransition) <- environment() 
  
  # Seed Harvest Rates, frequency, probability of harvest in a given year and intensity or percent of seed harvested
  #     Freq
  #Good X
  #Bad  Y
  # Good vs. bad - either determine the percent per good and bad after given clusters or just keep not in clusters and
  Freq_m <- matrix(Freq, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Freq")))
  Inten_m <- matrix(Int, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Inten")))
  
  # fecundity is the cutoff (practical sense of a good year with lots of seed vs few seed produced) but this doesn't really tell me if there is high production of seed, it's the ratio of recruitment to seed production
  f_all <- apply(mapply(function(x,y) x-y, Mx_list, TMx_list), 2, sum) # subtracting the matrices with fecundity from the matrices with all transitions leaves only the fecundity transitions, adding all remaining fecundity values
  mn_f <- mean(f_all) 
  #The index of above/equal and below average
  above_f <- which(f_all >= mn_f)
  below_f <- which(f_all < mn_f)
  
  # Initialize 
  ceil <- 10*StartPopSize
  Yr <- "Good"
  GBs <- rep(NA,TotYrs)
  Freqs <- rep(NA,TotYrs)
  Intens <- rep(NA,TotYrs)
  PopSize <- rep(NA, TotYrs)
  lambdas <- rep(NA, TotYrs)
  mats <- list()
  yrs <- 1
  
  # if clustered collection then does not depend on good or bad, will collect during a bad year, all match the larger frequency
  if(ClusteredColl > 1){
    numCollections <- max(Freq)*TotYrs
    numclusters <- floor(numCollections/ClusteredColl) # might underestimate
    i <- 0
    collectionyr_ix <- c()
    notavailable <- c()
    while(i < numclusters){
      available <- setdiff((1:(TotYrs-ClusteredColl)),notavailable) # but need to substract ClusteredColl from each start (clst) 
      clst <- sample(size = 1, x = available) # select from possible non-overlapping starts
      notavailable <- c(notavailable, (clst-ClusteredColl):(clst+ClusteredColl-1))
      collectionyr_ix <- c(collectionyr_ix,c(clst:(clst+ClusteredColl-1)))
      i <- i+1
      if(length(available)<ClusteredColl){
        warning("Cannot find non-overlapping cluster period, set at frequency")
        collectionyr_ix <- c(collectionyr_ix, sample(size = (numCollections-length(collectionyr_ix)), setdiff((1:TotYrs),collectionyr_ix)))
        break
      }
    } # allow adding overlapping cluster periods
    Freqs <- rep(0,TotYrs)
    Freqs[collectionyr_ix] <- 1
    
    while(yrs <= TotYrs){
      # Transition between good and bad years
      Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
        "Good"
      } else {
        "Bad"
      } # Make next a good or bad year
      
      CollectYN <- Freqs[yrs] # this was in the Annual transition
      onetransition <- AnnualTransition()
      GBs[yrs] <- Yr
      if(CollectYN==1){
        Intens[yrs] <- Inten_m[Yr,]
      }
      
      PopSize[yrs] <- floor(sum(vec))
      lambdas[yrs] <- onetransition[[2]]
      mats[[yrs]] <- onetransition[[3]]
      yrs <- yrs+1
    }
  } else {
    ######################## Annual Transitions without clusters ######################
    while(yrs <= TotYrs){
      # Transition between good and bad years
      Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
        "Good"
      } else {
        "Bad"
      } 
      CollectYN <- rbinom(1,1,Freq_m[Yr,]) # this was in the Annual transition
      onetransition <- AnnualTransition()  ### CHECK BACK HERE
      GBs[yrs] <- Yr
      Freqs[yrs] <- CollectYN
      if(CollectYN==1){
        Intens[yrs] <- Inten_m[Yr,]
      }
      
      vec <- onetransition[[1]] 
      PopSize[yrs] <- sum(vec)
      lambdas[yrs] <- onetransition[[2]]
      mats[[yrs]] <- onetransition[[3]]
      yrs <- yrs+1
    }
  } # ends all annual transitions without clusters
  SimDF <- data.frame(GB_yrs = GBs, Frequency = Freqs, Intensity = Intens, PopulationSize = PopSize, Year = 1:TotYrs, 
                      lambdas = lambdas)
  list(SimDF, FreqGB = Freq, IntGB = Int, GBMx = GoodBadTm, StartingPopSize = StartPopSize, Mx = mats,
       LogGrowthRate = stoch.growth.rate(mats, verbose=FALSE)) # $approx is log stochastic growth rate by Tuljapukar's approximation; $sim and CI
}




# starting population of only above ground individuals - pop size is above ground
# SeedHarvestSim <- function(Mx_list, TMx_list, Nx, Nx_list, StartPopSize, GoodBadTm, Freq = c(0,0), Int = c(0,0), 
#                            TotYrs, ClusteredColl = 1, ddceiling = FALSE){
#   # take the distribution of classes (likely the stable stage distribution) and scale to starting population size with integers
#   # Scale to starting population size with only the aboveground, detectable stages (no seedbank or dormant); undetectable stages have NA in Nx_list
#   Nx_scale <- Nx[!is.na(Nx_list[[1]])]
#   # vec <- matrix(floor(Nx*(StartPopSize/sum(Nx))), ncol = 1)
#   vec <- matrix(round(Nx*(StartPopSize/sum(Nx_scale)),0), ncol = 1)
#   # Annual transition function needs local environment variables
#   environment(AnnualTransition) <- environment() 
#   
#   # Number of seed
#   
#   
#   # Seed Harvest Rates, frequency, probability of harvest in a given year and intensity or percent of seed harvested
#   #     Freq
#   #Good X
#   #Bad  Y
#   # Good vs. bad - either determine the percent per good and bad after given clusters or just keep not in clusters and
#   Freq_m <- matrix(Freq, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Freq")))
#   Inten_m <- matrix(Int, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Inten")))
#   
#   # fecundity is the cutoff (practical sense of a good year with lots of seed vs few seed produced) but this doesn't really tell me if there is high production of seed, it's the ratio of recruitment to seed production
#   f_all <- apply(mapply(function(x,y) x-y, Mx_list, TMx_list), 2, sum) # subtracting the matrices with fecundity from the matrices with all transitions leaves only the fecundity transitions
#   mn_f <- mean(f_all) 
#   #The index of above/equal and below average
#   above_f <- which(f_all >= mn_f)
#   below_f <- which(f_all < mn_f)
#   
#   # Initialize 
#   ceil <- 10*StartPopSize
#   Yr <- "Good"
#   GBs <- rep(NA,TotYrs)
#   Freqs <- rep(NA,TotYrs)
#   Intens <- rep(NA,TotYrs)
#   PopSize <- rep(NA, TotYrs)
#   lambdas <- rep(NA, TotYrs)
#   mats <- list()
#   yrs <- 1
#   
#   # if clustered collection then does not depend on good or bad, will collect during a bad year, all match the larger frequency
#   if(ClusteredColl > 1){
#     numCollections <- max(Freq)*TotYrs
#     numclusters <- floor(numCollections/ClusteredColl) # might underestimate
#     i <- 0
#     collectionyr_ix <- c()
#     notavailable <- c()
#     while(i < numclusters){
#       available <- setdiff((1:(TotYrs-ClusteredColl)),notavailable) # but need to substract ClusteredColl from each start (clst) 
#       clst <- sample(size = 1, x = available) # select from possible non-overlapping starts
#       notavailable <- c(notavailable, (clst-ClusteredColl):(clst+ClusteredColl-1))
#       collectionyr_ix <- c(collectionyr_ix,c(clst:(clst+ClusteredColl-1)))
#       i <- i+1
#       if(length(available)<ClusteredColl){
#         warning("Cannot find non-overlapping cluster period, set at frequency")
#         collectionyr_ix <- c(collectionyr_ix, sample(size = (numCollections-length(collectionyr_ix)), setdiff((1:TotYrs),collectionyr_ix)))
#         break
#       }
#     } # allow adding overlapping cluster periods
#     Freqs <- rep(0,TotYrs)
#     Freqs[collectionyr_ix] <- 1
#     
#     while(yrs <= TotYrs){
#       # Transition between good and bad years
#       Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
#         "Good"
#       } else {
#         "Bad"
#       } # Make next a good or bad year
#       
#       CollectYN <- Freqs[yrs] # this was in the Annual transition
#       onetransition <- AnnualTransition()
#       GBs[yrs] <- Yr
#       if(CollectYN==1){
#         Intens[yrs] <- Inten_m[Yr,]
#       }
#       
#       PopSize[yrs] <- floor(sum(vec))
#       lambdas[yrs] <- onetransition[[2]]
#       mats[[yrs]] <- onetransition[[3]]
#       yrs <- yrs+1
#     }
#   } else {
#     ######################## Annual Transitions without clusters ######################
#     while(yrs <= TotYrs){
#       # Transition between good and bad years
#       Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
#         "Good"
#       } else {
#         "Bad"
#       } 
#       CollectYN <- rbinom(1,1,Freq_m[Yr,]) # this was in the Annual transition
#       onetransition <- AnnualTransition()
#       GBs[yrs] <- Yr
#       Freqs[yrs] <- CollectYN
#       if(CollectYN==1){
#         Intens[yrs] <- Inten_m[Yr,]
#       }
#       if(ddceiling == TRUE){
#         if(sum(onetransition[[1]]) < ceil){
#           vec <- onetransition[[1]] # output is a matrix with one column of stage class numbers, needs to be a vector
#         }
#       } else {
#         vec <- onetransition[[1]] # output is a matrix with one column of stage class numbers, needs to be a vector
#       }
#       
#       PopSize[yrs] <- sum(vec)
#       lambdas[yrs] <- onetransition[[2]]
#       mats[[yrs]] <- onetransition[[3]]
#       yrs <- yrs+1
#     }
#   } # ends all annual transitions without clusters
#   SimDF <- data.frame(GB_yrs = GBs, Frequency = Freqs, Intensity = Intens, PopulationSize = PopSize, Year = 1:TotYrs, 
#                       lambdas = lambdas)
#   list(SimDF, FreqGB = Freq, IntGB = Int, GBMx = GoodBadTm, StartingPopSize = StartPopSize, 
#        LogGrowthRate = stoch.growth.rate(mats, verbose=FALSE)) # $approx is log stochastic growth rate by Tuljapukar's approximation; $sim and CI
# }
