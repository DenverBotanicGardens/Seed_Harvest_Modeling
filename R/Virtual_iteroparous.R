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



MPM_itero <- function(params, stage1, stage2, maxSurv = 0.8, minSurv = 0.001){
  possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < 0.8 & 
                                                          params$survival > minSurv & params$age == stage1+1 &
                                                          !is.na(params$survival)],]
  possibleParams <- possibleParams[possibleParams$survival < maxSurv & possibleParams$age == stage1,]
  params1 <- possibleParams
  i <- sample(1:nrow(params1),1)
  S <- Stasis(alpha2 = params1[i,1], beta2 = params1[i,2], Age_first = 1, Age_last = stage1)
  S2 <- Stasis(alpha2 = params1[i,1], beta2 = params1[i,2], Age_first = stage1+1, Age_last = stage2)
   # Kendall et al 2019: fecundity multiplied by the survival of the offspring for prebreeding census (i.e. f = b_x * sigma_0 where parent in class x produces b_x seed after the census which surviv to the end of the timestep at rate sigma_0)
  f <- itero_fecundsurv(S2) * survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = 1)
  (G <- Growth(alpha2 = params1[i,1], beta2 = params1[i,2],Age_first = 1, Age_last = stage1))
  t_ij <- matrix(c(S, G,
                   f, S2),
                 nrow = 2)
  print(t_ij)
  print(colSums(t_ij))
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                        gentim = generation.time(t_ij), 
                                        AgeRep = stage1+1, lifelength = stage2), params1)
  return(listout)
}

#' @examples 
#' MPM_itero()
#' 
#' @describeIn MPM_semel


MPM_semel <- function(params, stage1, stage2, maxSurv = 0.8, minSurv = 0.0001){
  possibleParams <- params[params$survival < maxSurv & params$survival > minSurv & params$age == stage1+1 &
                             !is.na(params$survival),]
  params1 <- possibleParams
  i <- sample(1:nrow(params1),1)
  S <- Stasis(alpha2 = params1[i,1], beta2 = params1[i,2], Age_first = 1, Age_last = stage1)
  # Kendall et al 2019: fecundity multiplied by the survival of the offspring for prebreeding census (i.e. f = b_x * sigma_0 where parent in class x produces b_x seed after the census which surviv to the end of the timestep at rate sigma_0)
  f <- semel_fecundsurv(survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = stage2)) * survivalTypeI(alpha2 = params1[i,1], beta2 = params1[i,2], x = 1)
  (G <- Growth(alpha2 = params1[i,1], beta2 = params1[i,2],Age_first = 1, Age_last = stage1))
  t_ij <- matrix(c(S, G,
                   f, 0),
                 nrow = 2)
  print(t_ij)
  print(colSums(t_ij))
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                        gentim = generation.time(t_ij), 
                                        AgeRep = stage1+1, lifelength = stage2), params1)
  return(listout)
}


#' @describeIn MPM_annual One seed bank stage and one above ground reproductive stage


MPM_annual <- function(params = paramsAll, stage2 = 1, maxSurv = 0.8, minSurv = 0.0001, minGerm = 0.0001, maxGerm = 0.8){
  S <- runif(1,minSurv,maxSurv) # in the soil seed bank
  (f <- semel_fecundsurv(runif(1, 0.001,0.1))  )
  f2 <- f * runif(1, minGerm, maxGerm) # recruitment to next year adults is reduced by servival to rep of seed at t-1
  
  (G <- runif(1,minGerm, maxGerm))
  t_ij <- matrix(c(S, G,
                   f, S2),
                 nrow = 2)
  print(t_ij)
  print(colSums(t_ij))
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                        gentim = generation.time(t_ij, r=c(1,2), c=2), 
                                        AgeRep = stage2, lifelength = stage2), params1)
  return(listout)
}
