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


# Test
# alpha2 <- 0.76
alpha2 <- 0.94
beta2 < 0.27
# beta2 <- 0.23
# Age_first <- 2
# Age_last <- 10
Age_last <- 5

library(mvtnorm)
library(MultiRNG)
library(MASS)
library(NonNorMvtDist)
# library(Countr)
# library(dplyr)
# library(xtable)
library(rWishart)
library(DirichletReg)
# These functions provide information about the multivariate normal distribution 
# with mean equal to mean and covariance matrix sigma. 
# dmvnorm gives the density and rmvnorm generates random deviates.

beta_mom <- function(mean, var){
  term <- mean * (1 - mean)/var - 1
  alpha <- mean * term
  beta <- (1 - mean) * term
  if (var >= mean * (1 - mean)) stop("var must be less than mean * (1 - mean)")
  return(list(alpha = alpha, beta = beta))
}
beta.pars <- mapply(function(mu,v) beta_mom(mean = mu, var = v),
                                            apply(varcovdf, 2, mean), apply(varcovdf, 2, var))


# Do both stasis and growth for each surival curve - bivariate distribution 
# sig2 -> the variance and covariance matrix 
# pull from bivariate normal (just indexing...)
# Variance from var_cov matrix 
# promotion 1 ->2; stasis at 1; stasis 2
# three dimensional for iteroparous, semel is 2D
# params1 across all curves
# jointly distributed values
Age_first_stage2 <- 3
Age_last_stage2 <- 9

plot(0:15, survivalTypeI(params1[1,1], params1[1,2], 0:15), type = "l", col = rainbow(nrow(params_itero_fast), alpha = 1)[1],
     xlab = "Age in years",
     ylab = "Survival",
     ylim = c(0,1))
for(i in 2:nrow(params1)){
  lines(0:15, survivalTypeI(params1[i,1], params1[i,2], 0:15), type = "l", col = rainbow(nrow(params1), alpha = 0.4)[i])
}

# assuming equal proportion of individuals within each age within a stage class 
IteroStasisGrowth <- function(Age_first_stage2, Age_last_stage2, params1 = params1){
  varcovdf <- do.call(rbind,lapply(1:nrow(params1), function(i){
    alpha2 <- params1$a2[i]
    beta2 <- params1$b2[i]
    survivalcurve <- function(x){
      survivalTypeI(alpha2, beta2, x)
    }
    muSurv1 <- survivalTypeI(alpha2, beta2, 1)
    muStasis1 <- integrate(survivalcurve, lower = 1, upper = (Age_first_stage2 -1))$value/length(1:(Age_first_stage2-1))
    muStasis2 <- integrate(survivalcurve, lower = Age_first_stage2, upper = Age_last_stage2)$value/
      length(Age_first_stage2:Age_last_stage2)
    muGrowth <- survivalTypeI(alpha2, beta2, Age_first_stage2 -1)
    muDeath1 <- 1 - (muStasis1 + muGrowth)
    data.frame(muSurv1,muStasis1,muStasis2,muGrowth,muDeath1)
    }))
  # just bootstrap; copula
  varcovdf
}

hist(varcovdf$muStasis1)
hist(varcovdf$muStasis2)
hist(varcovdf$muGrowth)
hist(varcovdf$muDeath1)
# variance <- apply(varcovdf,2,var)
# covariance <- cov(varcovdf)
# vcov(lm(rep(1,nrow(varcovdf)) ~ muStasis1 + muStasis2 + muGrowth, data = varcovdf))
# varcovdf$Curve <- 1:nrow(varcovdf)
# form <- Curve ~ muSurv1 + muStasis1 + muStasis2 + muGrowth 
# gam <- renewalCount(formula = form, data = varcovdf, dist = "gamma", computeHessian = TRUE)
# mvnmle::mlest(varcovdf)$sigmahat
# draw.d.variate.normal(1, 3, mlest(varcovdf)$muhat, mlest(varcovdf)$sigmahat)
# NonNorMvtDist::rmvinvbeta(2, parm1 = mlest(varcovdf)$muhat, parm2 = mlest(varcovdf)$sigmahat)
# round(MASS::mvrnorm(2, mu = mlest(varcovdf)$muhat, Sigma =  mlest(varcovdf)$sigmahat),4)
# MultiRNG::draw.multivariate.laplace(2,4,gamma = 2, mu = mlest(varcovdf)$muhat, Sigma =  mlest(varcovdf)$sigmahat)
# rWishart::rSingularWishart(n = 2, df = nrow(varcovdf), Sigma = )
# Icmix::mvgamma(1, mu = mlest(varcovdf)$muhat, Sigma =  mlest(varcovdf)$sigmahat)
## Dirichlet
# mean
alpha_stasis <- varcovdf$muStasis1 
# precision
rdirichlet(100, alpha = c(mean(varcovdf$muStasis1), mean(varcovdf$muGrowth), mean(varcovdf$muDeath1)))
dircl <- rdirichlet(100, alpha = c(mean(varcovdf$muStasis1), mean(varcovdf$muGrowth), mean(varcovdf$muDeath1)))
plot(1:100,dircl[,1],type = "l")
lines(1:100, dircl[,2], type = "l", col="red")

SemelStasisGrowth <- function(Age_first_stage2, Age_last_stage2, params1 = params1){
  varcovdf <- do.call(rbind,lapply(1:nrow(params1), function(i){
    alpha2 <- params1$a2[i]
    beta2 <- params1$b2[i]
    survivalcurve <- function(x){
      survivalTypeI(alpha2, beta2, x)
    }
    muStasis1 <- integrate(survivalcurve, lower = 1, upper = (Age_first_stage2 -1))$value
    muGrowth <- survivalTypeI(alpha2, beta2, Age_first_stage2 -1)
    data.frame(muStasis1,muGrowth)
  }))
  round(MASS::mvrnorm(1, mu = mlest(varcovdf)$muhat, Sigma =  mlest(varcovdf)$sigmahat),4)
}



Stasis <- function(Age_first, Age_last, alpha2, beta2){
  if(Age_first < Age_last){
    survivalcurve <- function(x){
      survivalTypeI(alpha2, beta2, x)
    }
    # df <- get stasis1 mu, growth mu, and stasis 2 mu - paired with survival curve for each survival curve
    # or Nx3 matrix 
    # get 3D vector and return 
    # mls
    # survival curve shapes -> integrals for stasis or fixed point for growth (age_last) mlest --> 3D means
    # variance will be matrix - diagnol is var of each    "sighat" 
    # rmvnorm multivariate norm to pick from the variance/covariance that will pick the three values based on mus and var/cov structure
    
    
    mu <- (integrate(survivalcurve, lower = Age_first, upper = Age_last)$value) # /((Age_last+1) - Age_first)
      # survivalcurve(integrate(survivalcurve, lower = Age_first, upper = Age_last)$value)
      # DescTools::AUC(x = 1:20, y = survivalTypeI(alpha2, beta2, 1:20), from = Age_first, to = Age_last)
    # sig2 <- var(survivalTypeI(alpha2, beta2, Age_first:Age_last))
    # sig2 <- (sum(survivalTypeI(alpha2, beta2, Age_first:Age_last)-mu)^2)/((Age_last+1)-Age_first)
      # sigma(lm(survivalTypeI(alpha2, beta2, Age_first:Age_last)~ c(Age_first:Age_last)))^2
    alpha1 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
    beta1 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
    stasis <- rbeta(1,alpha1, beta1)
      # rbeta(100,alpha1/(alpha1 + beta1), (alpha1 * beta1)/((alpha1 + beta1)^2 * (alpha1 + beta1 + 1)) )
      # beta(alpha1/(alpha1 + beta1), (alpha1 * beta1)/((alpha1 + beta1)^2 * (alpha1 + beta1 + 1)))
      # rbeta(1,shape1 = alpha1, shape2 = beta1)

    } else stasis <- 0 
    stasis
}



#' @describeIn Growth The proportion of individuals that mature (survive) from the beginning of stage x 
#'   to the end of stage x  
#' @param Age_first The first age of stage x
#' @param Age_last The last age of stage x
#' @param alpha2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)
#' @param beta2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)

# alpha2 <- .82
# beta2 <- .35
Age_first <- stage1
Age_last <- stage2-1
# params1 
alpha2 <- params1$a2[1000]
beta2 <- params1$b2[1000]

Growth <- function(Age_first, Age_last, params1, alpha2, beta2){
  # mu 
  if(Age_last > Age_first){
  mu <- matching_proportion <- survivalTypeI(alpha2, beta2, Age_last) #/ # probability of promotion
    # sum(survivalTypeI(alpha2,beta2,Age_first:Age_last))
  sig2 <- sum((sapply(1:nrow(params1), function(i){
      survivalTypeI(params1[i,1], params1[i,2], Age_last)/
        sum(survivalTypeI(params1[i,1], params1[i,2], Age_first:Age_last))
      }) - mu)^2 )/nrow(params1)
  alpha1 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
  beta1 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
  rbeta(1,alpha1,beta1)
  } else {
    survivalTypeI(alpha2, beta2, Age_last)
  }
  
  # matching_proportion
  # matching_proportion * (integrate(survivalcurve, lower = Age_first, upper = Age_last)$value)/(Age_last - Age_first)
  # This or add variance (what is the variance among all the curves)
  # these are beta - proportions 
  # matching_proportion * sum(survivalTypeI(alpha2, beta2, 
  #                                         Age_first:Age_last))/(length(Age_first:Age_last))
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


