#' @title Seed Harvest Simulation
#' @description Simulation of climate and seed harvest frequency and intensity
#' 
#' @export
#' 
#' @param Mx_list A list of annual transition matrices indicating survival and growth from col to row; fecundity in the first row (square matrix with N stages).
#' @param TMx_list A list of annual transition matrices without fecundity (square matrix with N stages).
#' @param Nx A population size vector by stage (length N). Stable stage distriubiton or most recent count by stage. 
#' @param StartPopSize A number for the starting population size of the simulation. 
#' @param GoodBadTm A transition matrix (row to col) of good and bad years. Annual transition matrices are defined as bad < mean fecundity <= good 
#' @param Freq A numerical vector for the the probability of collecting seed given the year c(good, bad).
#' @param Int A numerical vector for the intensity (the percentage decrease in fecundity) if seed harvest for c(good, bad).
#' @param TotYrs A number indicating how many years to project in the simulation.
#' @param ClusteredColl "Y" or "N" indicating if collection years should be clustered. Results in increased likelihood of collection in t+1 when collecting in year t
#' @param seedbank TRUE or FALSE to skip the first element in the first row when reducing fecundity. 
#' 
#' @return 

SeedHarvestSim <- function(Mx_list, TMx_list, Nx, StartPopSize, GoodBadTm, Freq = c(0,0), Int = c(0,0), 
                           TotYrs, ClusteredColl = "N", seedbank = FALSE){

  vec <- floor(Nx*(StartPopSize/sum(Nx)))
  
  # Set a maximum population size of 4*StartPopSize or set flexible as function input??? Or 10 times, wishfull thinking??
  dendep <- 10*StartPopSize
  
  # Seed Harvest Rates, frequency, probability of harvest in a given year and intensity or percent of seed harvested
  Freq_m <- matrix(Freq, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Freq")))
  Inten_m <- matrix(Int, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Inten")))
  
  # if Lambda is the cutoff for good vs. bad year
  l_all <- sapply(Mx_list, function(x)  lambda(t(x)), simplify = FALSE, USE.NAMES = TRUE) 
  mn_l <- mean(unlist(l_all))
  above_l <- Mx_list[which(l_all >= mn_l)]
  below_l <- Mx_list[which(l_all < mn_l)]
  # If fecundity is the cutoff (practical sense of a good year with lots of seed vs few seed produced) but this doesn't really tell me if there is high production of seed, it's the ratio of seed production to recruitment
  f_all <- mapply(function(x,y){
    fec <- x-y
    rowSums(fec)[1] #Only the first row is funcunity across classes (subtracted out survival by Mx-TMx)
  }, Mx_list, TMx_list)
  mn_f <- mean(unlist(f_all)) 
  above_f <- Mx_list[which(f_all >= mn_f)]
  below_f <- Mx_list[which(f_all < mn_f)]
  
  # Start with good! positive
  Yr <- "Good"
  GBs <- rep(NA,TotYrs)
  Freqs <- rep(NA,TotYrs)
  Intens <- rep(NA,TotYrs)
  PopSize <- c(sum(vec),rep(NA, TotYrs))
  yrs <- 1
  while(yrs <= TotYrs){
    if(Yr == "Good"){
      nx <- above_f[[sample(1:length(above_f), 1)]] # randomly sample from the high fecundty transition matrices
    } else {
      nx <- below_f[[sample(1:length(below_f), 1)]]
    }
    # Freq likely to collect seed in good or bad years
    CollectYN <- rbinom(1,1,Freq_m[Yr,])
    if(CollectYN == 1){
      # yes collect, reduce fecuntity by how much?
      if(seedbank == TRUE){ 
        nx[1,-1] <- nx[1,-1]*Inten_m[Yr,] 
      } else {
        nx[1,] <- nx[1,]*Inten_m[Yr,]
      }
    }
    
    # Multiple N_t+1 = M*N_t
    if(sum(nx%*%vec) < dendep) vec <- nx%*%vec
    
    GBs[yrs] <- Yr
    Freqs[yrs] <- CollectYN
    if(CollectYN==1){
      Intens[yrs] <- Inten_m[Yr,]
    }
    PopSize[(yrs+1)] <- sum(vec)
    # Transition between good and bad years
    Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
      "Good"
    } else {
      "Bad"
    } 
    Yr
    yrs <- yrs+1
  }  
  list(GB_yrs = GBs, Frequency = Freqs, Intensity = Intens, PopulationSize = PopSize)
}
