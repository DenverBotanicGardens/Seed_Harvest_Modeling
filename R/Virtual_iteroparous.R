#' @title Matrix model for an iteroparous species.
#' @description  Lefkovich stage based population, two by two matrix with one vegetative stage and one reproductive stage
#' @author Michelle DePrenger-Levin
#' @export
#' 
#' @param params A data frame of survival curve parameters (alpha2 and beta2) from a Type I survival curve
#' @param stage1 The oldest age of stage 1
#' @param stage2 The oldest age of stage 2
#' @param Age_first The first age of stage x
#' @param Age_last The last age of stage x
#' @param alpha2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)
#' @param beta2 parameter for a Type I survival curve (from Fujiwara and Diaz-Lopez 2017)

# params <- paramsAll
maxSurvRep <- 0.8
maxSurvVeg <- 0.8
minSurv <- 0.1
stage1 <- 2
stage2 <- 9

MPM_itero <- function(params, stage1, stage2, maxSurvVeg = 0.99, maxSurvRep = 0.8, minSurv = 0.001){
  possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurvRep & 
                                                          params$survival > minSurv & params$age == stage1+1],]
  possibleParams <- possibleParams[possibleParams$survival < maxSurvVeg & possibleParams$age == 1,]
  params1 <- possibleParams
  
  plot(0:15, survivalTypeI(params1[1,1], params1[1,2], 0:15), type = "l", col = rainbow(nrow(params_itero_fast), alpha = 1)[1],
       xlab = "Age in years",
       ylab = "Survival",
       ylim = c(0,1))
  for(i in 2:nrow(params1)){
    lines(0:15, survivalTypeI(params1[i,1], params1[i,2], 0:15), type = "l", col = rainbow(nrow(params1), alpha = 0.4)[i])
  }
  # i <- sample(1:nrow(params1),1)
  # S <- Stasis(alpha2 = params1[i,1], beta2 = params1[i,2], Age_first = 1, Age_last = stage1)
  # S2 <- Stasis(alpha2 = params1[i,1], beta2 = params1[i,2], Age_first = stage1+1, Age_last = stage2)
  # Kendall et al 2019: fecundity multiplied by the survival of the offspring for prebreeding census 
  # (i.e. f = b_x * sigma_0 where parent in class x produces b_x seed after the census which survive 
  # to the end of the timestep at rate sigma_0)
  SS1S2G <- IteroStasisGrowth(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, params1 = params1)
  S <- SS1S2G[2]
  S2 <- SS1S2G[3]
  G <- SS1S2G[4]
  f <- itero_fecundsurv(S2) * SS1S2G[1] # survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = 1)
  # f <- itero_fecundsurv(S2) * survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = 1)
  # (G <- Growth(alpha2 = params1[i,1], beta2 = params1[i,2], params1 = params1, Age_first = 1, Age_last = stage1))
  # if(S+G > maxSurvVeg){
  #   S <- S*(maxSurvVeg/(S+G))
  #   G <- G*(maxSurvVeg/(S+G))
  # }
  t_ij <- matrix(c(S, G,
                   f, S2),
                 nrow = 2)
  # print(t_ij)
  if(colSums(t_ij)[1] > 1){ 
    print(colSums(t_ij))
    print(lambda(t_ij))
    }
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                        gentim = generation.time(t_ij), 
                                        netrep = net.reproductive.rate(t_ij),
                                        AgeRep = stage1+1, lifelength = stage2), params1)
  return(listout)
}

#' @examples 
#' MPM_itero()
#' 
#' @describeIn MPM_semel

# test
params <- paramsAll
maxSurv <- 0.4
minSurv <- 0.0001
stage1 <- 1 # max(1,rpois(1,1))
stage2 <- stage1+1

MPM_semel <- function(params, stage1, stage2, maxSurv = 0.8, minSurv = 0.0001){
  possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurv & params$survival > minSurv & 
                                                          params$age == stage2],] 
  possibleParams <- possibleParams[possibleParams$survival < 0.99 & possibleParams$age == 1,]
  possibleParams <- possibleParams[possibleParams$survival > 0.9 & possibleParams$age == 1,]
  params1 <- possibleParams
  i <- sample(1:nrow(params1),1)
  (S <- Stasis(alpha2 = params1[i,1], beta2 = params1[i,2], Age_first = 1, Age_last = stage1))
  # Kendall et al 2019: fecundity multiplied by the survival of the offspring for prebreeding census 
  # (i.e. f = b_x * sigma_0 where parent in class x produces b_x seed after the census which survive 
  # to the end of the timestep at rate sigma_0)
  (f <- semel_fecundsurv(survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = stage2)) * 
      survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = 1))
  # (f <- itero_fecundsurv(survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = stage2)) * 
  # survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = 1))
  (G <- Growth(alpha2 = params1[i,1], beta2 = params1[i,2],params1 = params1, Age_first = 1, Age_last = stage1))
  # if(S+G > maxSurv){
  #   S <- S*(maxSurv/(S+G))
  #   G <- G*(maxSurv/(S+G))
  # }
  t_ij <- matrix(c(S, G,
                   f, 0),
                 nrow = 2)
  if(colSums(t_ij)[1] > 1) print(colSums(t_ij))
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                        gentim = generation.time(t_ij), 
                                        netrep = net.reproductive.rate(t_ij),
                                        AgeRep = stage1+1, lifelength = stage2), params1)
  return(listout)
}


#' @describeIn MPM_annual One seed bank stage and one above ground reproductive stage


MPM_annual <- function(params = paramsAll, stage2 = 1, maxSurv = 0.8, minSurv = 0.0001, minGerm = 0.0001, maxGerm = 0.8, 
                       persSeedBank = 3){
  S <- runif(1,minSurv,maxSurv) # in the soil seed bank
  (f <- semel_fecundsurv(runif(1, 0.001,0.1))  )
  f2 <- f * runif(1, minSurv, maxSurv) # recruitment to next year adults is reduced by servival to rep of seed at t-1
  (G <- S/(sum(rep(S,persSeedBank))) * # rate individuals reach the end of the first stage
      S) # the average survival of the stage (soil seed bank) 
  t_ij <- matrix(c(S, G,
                   f, f2),
                 nrow = 2)
  print(t_ij)
  print(colSums(t_ij))
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[1])
  growthElast <- sum(e_ij[2])
  fecundElast <- sum(e_ij[c(3,4)])
  listout <- list(t_ij,e_ij, 
                  data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                             gentim = generation.time(t_ij, r=c(1,2), c=2), 
                             netrep = net.reproductive.rate(t_ij, r=c(1,2), c=2),
                             AgeRep = stage2, lifelength = stage2))
  return(listout)
}
