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

#' @describeIn paramsAll This creates a data frame of the parameter values and survival at ages 1:30 years for a type I curve

paramsAll <- do.call(rbind, lapply(seq(0.01,1, by=0.01), function(a2){
  outb2 <- do.call(rbind, lapply(seq(0.01,1, by=0.01), function(b2){
    surv <- survivalTypeI(a2, b2, 1:10)
    data.frame(a2, b2, age = 1:30, survival = surv, a2b2 = paste(a2,b2,sep=""))
  }))
}))

# constant risk of mortality with age
#' @describeIn paramsAll This creates a data frame of the parameter values and survival at ages 1:30 years
survivalTypeIII <- function(alpha1, x){
  exp(-alpha1 * x)
}

# exponentially declining survivorship curve with age with constant risk; use with slow
paramsAll_typeIII <- do.call(rbind, lapply(seq(0.001, 1, by = 0.001), function(a1){
  data.frame(a1,age=1:30, survival = survivalTypeIII(a1, 1:30))
}))

#' @section Vital Rates (S)tasis, (G)rowth, and (R)eproduction
#' @describeIn itero_fecundsurv Concave for iteroparous from Takada and Kawai 2020; represents the output of seed
#'   but does not account for survival or recruitment to time t+1; upwardly convex and concave down ??? (0,0.8); 
#'   Farrell 2020: h(x) : [0, 1] ??? [0, 1] is twice differentiable with h(0) = 1, h(1) = 0 and h(x) ??? 0.
#'   Further, h is either strictly concave (i.e. h''(x) < 0 for x ???(0, 1), also called 'concave down'), 
#'   linear (i.e. h''(x) ??? 0 for x ??? (0, 1)) or strictly convex (i.e. h''(x) < 0 for x ??? (0, 1), 
#'   also called 'concave up').
#' @param s survival (mean) of the reproductive stage

itero_fecundsurv <- function(s){
  -12.5*((s + 0.1)^2) + 10.125
}

# A negative second derivative implies a concave down function
it_fecsur <- expression(-12.5*((s + 0.1)^2) + 10.125)
D(D(it_fecsur, 's'),'s')


#' @describeIn semel_fecundsurv Concave for semelparous from Takada and Kawai 2020; represents the output of seed
#'   but does not account for survival or recruitment to time t+1; convex (concave up or convex down)
#' @param s survival of the reproductive stage 

semel_fecundsurv <- function(s){
  12.5*((s - 0.9)^2) - 0.125
}

# A positive second derivative implies a concave up or convex (convex down) function
se_fecsur <- expression(12.5*((s - 0.9)^2) - 0.125)
D(D(se_fecsur,'s'),'s')

#' @describeIn Stasis The product of survival at age x+1/x for ages within a stage class 
#' @param Age_first The first age of stage x
#' @param Age_last The last age of stage x
#' @param alpha2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)
#' @param beta2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)

#Test
# alpha2 <- params[i,1]
# beta2 <- params[i,2]
# Age_first <- 2
# Age_last <- 3

# beta beta(alpha, beta) - continuous random variables [0-1] a proportion, survival, transition
# method of moments to parameterize matrix elements from survival curves (mu = mean and sig2 = variance within stage) that keep lambda ca. 1
# alpha1 <- (mu^2 - mu^3 - mu*sig2)/sig2
# beta1 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2


library(DirichletReg)
# Test
# Age_first_stage2 = 6
# Age_last_stage2 = 7
# Age_last_stage2 <- 20

SlowVitalRates <- function(Age_first_stage2, Age_last_stage2, alpha1 = 0.201){
  # Weighted arithmetic mean where weight (proportion of individuals at that age) is survival at that age
  muSurv1 <- survivalTypeIII(alpha1, 1)
  # age specific survival to next age from last one
  s_x <- unlist(lapply(1:(Age_first_stage2 - 1), function(x){ 
    survivalTypeIII(alpha1,x+1)/survivalTypeIII(alpha1,x)
    }))
  w_vegx <- survivalTypeIII(alpha1,1:(Age_first_stage2 -1))
  muStasis1 <- sum(apply(cbind(s_x,w_vegx), 1, function(x) prod(x)))/
    (sum(survivalTypeIII(alpha1, 1:(Age_first_stage2 -1)))) 
  
  s_xrep <- unlist(lapply(Age_first_stage2:Age_last_stage2, function(x){
    survivalTypeIII(alpha1, x+1)/survivalTypeIII(alpha1, x)
  }))
  w_vegxrep <- survivalTypeIII(alpha1, Age_first_stage2:Age_last_stage2)
  muStasis2 <- sum(apply(cbind(s_xrep,w_vegxrep), 1, function(x) prod(x)))/
                        (sum(survivalTypeIII(alpha1, Age_first_stage2:Age_last_stage2))) 
  # Growth should be of surviving at final age, what percent of that grows? 
  muGrowth <- survivalTypeIII(alpha1, Age_first_stage2-1)/
                 (sum(survivalTypeIII(alpha1, 1:(Age_first_stage2-1)))) 
  # if(muGrowth > 0.99) muGrowth <- survivalTypeIII(alpha1, Age_first_stage2-1)
  # if(muStasis1 > 0.99) muStasis1 <- survivalTypeIII(alpha1, 1)
  muDeath1 <- 1 - (muStasis1 + muGrowth)
  if(muDeath1 < 0){
    Dirichletvital <- rdirichlet(1, c(muStasis1, muGrowth, 0.01))
    muStasis1 <- Dirichletvital[1]
    muGrowth <- Dirichletvital[2]
    muDeath1 <- Dirichletvital[3]
  }
  
  data.frame(muSurv1,muStasis1,muStasis2,muGrowth,muDeath1)
}

# Age_first_stage2 <- 1
# Age_last_stage2 <- 3


FastVitalRates <- function(Age_first_stage2, Age_last_stage2, alpha2 = 0.05, beta2 = 0.91){
  # Weighted arithmetic mean where weight (proportion of individuals at that age) is survival at that age
  muSurv1 <- sum(survivalTypeI(alpha2, beta2, 1)*survivalTypeI(alpha2, beta2, 1))/(sum(survivalTypeI(alpha2, beta2, 1)))
  s_x <- unlist(lapply(1:(Age_first_stage2 - 1), function(x){ 
    survivalTypeIII(alpha2, beta2, x+1)/survivalTypeIII(alpha2, beta2, x)
  }))
  w_vegx <- survivalTypeIII(alpha2, beta2,1:(Age_first_stage2 -1))
  
  muStasis1 <- sum( survivalTypeI(alpha2, beta2, 1:(Age_first_stage2 -1)) * survivalTypeI(alpha2, beta2, 1:(Age_first_stage2 -1)) )/
    (sum(survivalTypeI(alpha2, beta2, 1:(Age_first_stage2 -1)))) 
  muStasis2 <- sum(survivalTypeI(alpha2, beta2, Age_first_stage2:Age_last_stage2)* survivalTypeI(alpha2, beta2, Age_first_stage2:Age_last_stage2 ))/
    (sum(survivalTypeI(alpha2, beta2, Age_first_stage2:Age_last_stage2))) 
  muGrowth <- (survivalTypeI(alpha2, beta2, Age_first_stage2-1)/(sum(survivalTypeI(alpha2, beta2, 1:(Age_first_stage2-1)))))
  if(muGrowth > 0.99) muGrowth <- survivalTypeI(alpha2, beta2, Age_first_stage2-1)
  if(muStasis1 > 0.99) muStasis1 <- survivalTypeI(alpha2, beta2, Age_first_stage2-1)
  muDeath1 <- 1 - (muStasis1 + muGrowth)
  if(muDeath1 < 0){
    Dirichletvital <- rdirichlet(1, c(muStasis1, muGrowth, 0.01))
    muStasis1 <- Dirichletvital[1]
    muGrowth <- Dirichletvital[2]
    muDeath1 <- Dirichletvital[3]
  }
  data.frame(muSurv1,muStasis1,muStasis2,muGrowth,muDeath1)
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


