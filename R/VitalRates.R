#' Vital rates for creating fast and slow, itero and semelparous, and annual matrices
#' 
#' A type I survival curve was most common among listed plants (Salguero-Gomez 2017) 
#'   equation 2 from Fujiwara and Diaz-Lopez was used to estimate growth from stage 1 to stage 2 
#'   and stasis within stage 1 and stage 2 
#'   An fertility is based on the survival of the reproductive stage with a concave slope for 
#'   iteroparous and covex slope for semelparous species where survival is limited to [0,0.8] and fecundity is [0,10]
#' 
#' @description Calculate stasis (along the diagonal), growth (below the diagonal), and fecundity for a virtual species matrix population model
#' @author Michelle DePrenger-Levin
#' @export
#' @section Survival
#' @describeIn equation 2 from Fujiwara and Diaz-Lopez 2017 
#'   The x is age Fujiwara and Diaz-Lopez 2017; hazard is h(x) = alpha2*exp(beta2*x); 
#'   exponentially increasing risk of mortality with age, risk due to aging
survivalTypeI <- function(alpha2, beta2, x){
  exp((alpha2/beta2)*(1-(exp(beta2*x))))
}

#' @describeIn paramsAll This creates a data frame of the parameter values and survival at ages 1:30 years

paramsAll <- do.call(rbind, lapply(seq(0.01,1, by=0.01), function(a2){
  outb2 <- do.call(rbind, lapply(seq(0.01,1, by=0.01), function(b2){
    surv <- survivalTypeI(a2, b2, 1:30)
    data.frame(a2, b2, age = 1:30, survival = surv, a2b2 = paste(a2,b2,sep=""))
  }))
}))

#' @section Vital Rates (S)tasis, (G)rowth, and (R)eproduction
#' @describeIn itero_fecundsurv Concave for iteroparous from Takada and Kawai 2020; represents the output of seed
#'   but does not account for survival or recruitment to time t+1
#' @param s survival (mean) of the reproductive stage

itero_fecundsurv <- function(s){
  -12.5*((s + 0.1)^2) + 10.125
}

#' @describeIn semel_fecundsurv Concave for semelparous from Takada and Kawai 2020; represents the output of seed
#'   but does not account for survival or recruitment to time t+1
#' @param s survival of the reproductive stage 

semel_fecundsurv <- function(s){
  12.5*((s - 0.9)^2) - 0.125
}


#' @describeIn Stasis The product of survival at age x+1/x for ages within a stage class 
#' @param Age_first The first age of stage x
#' @param Age_last The last age of stage x
#' @param alpha2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)
#' @param beta2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)

Stasis <- function(Age_first, Age_last, alpha2, beta2){
  if(Age_first < (Age_last)){
    stasis <- prod(unlist(lapply(Age_first:(Age_last), function(i){
      survivalTypeI(alpha2, beta2,i+1)/survivalTypeI(alpha2, beta2,i)
    })))
  } else stasis <- 0 
  stasis
}

#' @describeIn Growth The proportion of individuals that mature (survive) from the beginning of stage x 
#'   to the end of stage x  
#' @param Age_first The first age of stage x
#' @param Age_last The last age of stage x
#' @param alpha2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)
#' @param beta2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)
Growth <- function(Age_first, Age_last, alpha2, beta2){
  matching_proportion <- survivalTypeI(alpha2, beta2, Age_last)/
    sum(survivalTypeI(alpha2,beta2,Age_first:Age_last))
  matching_proportion * sum(survivalTypeI(alpha2, beta2, 
                                          Age_first:Age_last))/(length(Age_first:Age_last))
}

#' @describeIn matrix_elast Make a matrix of specified dimensions to sum growth, stasis, and reproductive vital rates
#'   from matrix elements
#' @param Mx_dim Dimensions of an NxN square matrix 
matrix_elast <- function(Mx_dim = 2){
  x <- matrix(rep("L",Mx_dim^2), ncol = Mx_dim)
  x[lower.tri(x)] <- "G"
  x[upper.tri(x)] <- "R" # retrogressive growth
  x[1,Mx_dim] <- "F" # fecundity
  x
}

#' @examples  
matrix_elast(4)
generic_mat <- matrix_elast(2)


