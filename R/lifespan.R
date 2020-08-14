#' @title lifespan
#' @description Take a projection matrix, a population of 100 seedlings (recruits), remove the reproductive transitions (first row all zero) and measure years until the 100 individuals are dead
#' @author Dan Doak
#' @export
#' @param nx A Matrix population model or an Integral projection model square matrix with growth, survival, and fecundity (in the top row only) as transitions in each kernel or cell
#' @param seedbank TRUE/FALSE to feed 100 individuals to the smallest stage/age/size class instead of the seed bank.

lifespan <- function(nx, seedbank = FALSE){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))
  nx[1,]=0
  jj=1
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
    #print(sum(vec))
  }
  return(jj)
}
