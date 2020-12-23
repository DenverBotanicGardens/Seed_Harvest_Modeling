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


AnnualTransition <-function(){
  # above_f <- Mx_list[which(f_all >= mn_f)]  from the SeedHarvestSim function get index of one above 
  # below_f <- Mx_list[which(f_all < mn_f)]   from the SeedHarvestSim function get index of one below
  if(Yr == "Good"){
    Mx_i <- sample(above_f, 1) # randomly sample from the high fecundity transition matrices
  } else {
    Mx_i <- sample(below_f, 1)
  }
  
  # Freq[] is the likelihood to collect seed in good or bad years
  if(CollectYN == 1){  # yes collect, reduce fecuntity by the year good or bad
    fecundMx <- Mx_list[[Mx_i]] - TMx_list[[Mx_i]] # Subtract the transition matrix (Mx) by the matrix without fecundity rate (TMx) to get only fecundity rate
    fecundMx <- fecundMx - fecundMx*Inten_m[Yr,] # Reduce all fecundity rates by intensity of harvest
    nx <- fecundMx + TMx_list[[Mx_i]] # add the reduced fecundity rates to the transition matrix that lacks the fecundity rate to get new transition matrix
  } else {
    nx <- Mx_list[[Mx_i]]
  }
  K <- 100*StartPopSize  # maximum population size
  if(ddceiling == TRUE & sum(vec[!is.na(Nx_list[[1]])]) >= K*0.9){
    vec.t1 <- floor(nx%*%vec)
    vec <- floor(vec.t1/(lambda(nx))*(1+(lambda(nx)-1)*(1-(sum(vec[!is.na(Nx_list[[1]])])/K ))))  # dominant eigenvalue is lambda and the same as eigen(nx)$values[1]; what now? 1+(lambda(nx)-1) is just lambda(nx)
  } else {
    # Multiply N_t+1 = Matrix*N_t; vec comes from the SeedHarvestSim function
    vec <- floor(nx%*%vec)
  }
  list(vec, lambda(nx), projmat = nx)
}

# AnnualTransition <-function(){
#   # above_f <- Mx_list[which(f_all >= mn_f)]  from the SeedHarvestSim function get index of one above 
#   # below_f <- Mx_list[which(f_all < mn_f)]   from the SeedHarvestSim function get index of one below
#   if(Yr == "Good"){
#     Mx_i <- sample(above_f, 1) # randomly sample from the high fecundity transition matrices
#   } else {
#     Mx_i <- sample(below_f, 1)
#   }
#   
#   # Freq[] is the likelihood to collect seed in good or bad years
#   if(CollectYN == 1){  # yes collect, reduce fecuntity by the year good or bad
#     fecundMx <- Mx_list[[Mx_i]] - TMx_list[[Mx_i]] # Subtract the transition matrix (Mx) by the matrix without fecundity rate (TMx) to get only fecundity rate
#     fecundMx <- fecundMx - fecundMx*Inten_m[Yr,] # Reduce all fecundity rates by intensity of harvest
#     nx <- fecundMx + TMx_list[[Mx_i]] # add the reduced fecundity rates to the transition matrix that lacks the fecundity rate to get new transition matrix
#   } else {
#     nx <- Mx_list[[Mx_i]]
#   }
#   K <- 100*StartPopSize  # maximum population size
#   if(ddceiling == TRUE & sum(vec[!is.na(Nx_list[[1]])]) >= K*0.9){
#     vec.t1 <- floor(nx%*%vec)
#     vec <- floor(vec.t1/(lambda(nx))*(1+(lambda(nx)-1)*(1-(sum(vec[!is.na(Nx_list[[1]])])/K ))))  # dominant eigenvalue is lambda and the same as eigen(nx)$values[1]; what now? 1+(lambda(nx)-1) is just lambda(nx)
#   } else {
#     # Multiply N_t+1 = Matrix*N_t; vec comes from the SeedHarvestSim function
#     vec <- floor(nx%*%vec)
#   }
#   list(vec, lambda(nx), projmat = nx)
# }