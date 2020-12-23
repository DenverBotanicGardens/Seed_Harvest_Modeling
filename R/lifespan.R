#' @title lifespan
#' @description Take a projection matrix, a population of 100 seedlings (recruits), remove the reproductive transitions (first row all zero) and measure years until the 100 individuals are dead
#' @author Dan Doak
#' @export
#' @param nx A Matrix population model or an Integral projection model square matrix with growth, survival, and fecundity (in the top row only) as transitions in each kernel or cell
#' @param seedbank TRUE/FALSE to feed 100 individuals to the smallest stage/age/size class instead of the seed bank.


lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))
  nx[1,]=0
  jj=1
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
  }
  return(jj)
}

# Seeds a population with 100 individuals in the smallest size class, removes the top row (so no more reproduction), and follows the population until all have died. Based on a list of transition matrices (with fecundity rate: Mx, and without fecundity rate included: TMx)
lifespanMPM <- function(Mx, TMx, seedlingBin = 1){
  nclasses=dim(Mx[[1]])[1]
  # add 100 individuals to the smallest non-seedbank category
  vec <- rep(0, nclasses)
  vec[seedlingBin] <- 100
  # remove all fecundity from matrices
  nx <- mapply(function(x,y) x - (x-y), Mx, TMx, SIMPLIFY = FALSE) # keep output as matrices
  # if any transition is 100%
  nx[[1]][nx[[1]] == 1] <- 0.9
  jj=1
  while (sum(vec)>1){
    vec=nx[[sample(1:length(Mx),1)]]%*%(floor(vec)) # randomly sample from the list of matrices with fecundity rate removed
    jj=jj+1
  }
  return(jj)
}

# lifespan <- function(nx, seedbank = FALSE){
#   nclasses=dim(nx)[1]
#   if(seedbank){
#     vec <- c(0, 100, rep(0,(nclasses-2))) # so that 'life' starts at germination 
#   } else {
#     vec=c(100,rep(0,(nclasses-1)))
#   }
#   nx[1,]=0
#   jj=1
#   while (sum(vec)>1){
#     vec=nx%*%vec
#     jj=jj+1
#   }
#   return(jj)
# }
