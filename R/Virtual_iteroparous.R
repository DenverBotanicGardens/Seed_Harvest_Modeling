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
# maxSurvRep <- 0.8
# maxSurvVeg <- 0.95
# minSurv <- 0.1
# stage1 <- 10
# stage2 <- 20
# NumMats <- 100

MPM_itero <- function(params, stage1, stage2, maxSurvVeg = 0.99, maxSurvRep = 0.8, minSurv = 0.001, NumMats = 100){
  possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurvRep & 
                                                          params$survival > minSurv & params$age == stage1+1],]
  possibleParams <- possibleParams[possibleParams$a2b2 %in% 
                                     possibleParams$a2b2[possibleParams$survival < maxSurvVeg & possibleParams$age == 1],]
  params1 <- possibleParams
  vitalrates <- IteroStasisGrowth(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, params1 = params1)
  mats <- lapply(sample(1:nrow(vitalrates), size = NumMats, replace = TRUE), function(i){
    S <- vitalrates[i,2]
    G <- vitalrates[i,4]
    S2 <- vitalrates[i,3]
    f <- itero_fecundsurv(S2) * vitalrates[i,1]
    t_ij <- matrix(c(S, G,
                     f, S2),
                   nrow = 2)
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
  })
  return(mats)
}

#' @examples 
#' MPM_itero()
#' 
#' @describeIn MPM_semel

# test
# params <- paramsAll
# maxSurv <- 0.4
# minSurv <- 0.0001
# stage1 <- 1 # max(1,rpois(1,1))
# stage2 <- stage1+1

MPM_semel <- function(params, stage1, stage2, maxSurv = 0.8, minSurv = 0.0001, vegminSurv = 0.9, vegmaxSurv = 0.99, NumMats = 100){
  possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurv & params$survival > minSurv & 
                                                          params$age == stage2],] 
  possibleParams <- possibleParams[possibleParams$survival < vegmaxSurv & possibleParams$age == 1,]
  possibleParams <- possibleParams[possibleParams$survival > vegminSurv & possibleParams$age == 1,]
  params1 <- possibleParams
  vitalrates <- IteroStasisGrowth(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, params1 = params1)
  mats <- lapply(sample(1:nrow(vitalrates), size = NumMats, replace = TRUE), function(i){
    S <- vitalrates[i,2]
    G <- vitalrates[i,4]
    S2 <- vitalrates[i,3]
    f <- semel_fecundsurv(S2) * vitalrates[i,1]
    t_ij <- matrix(c(S, G,
                     f, 0),
                   nrow = 2)
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
  })
  return(mats)
}


#' @describeIn MPM_annual One seed bank stage and one above ground reproductive stage

# test
# params <- paramsAll
# seedminSurv <- 0.001
# seedmaxSurv <- 0.99

MPM_annual <- function(params, stage1 = 3, stage2 = 4, maxSurv = 0.8, minSurv = 0.0001, seedminSurv = 0.001, 
                       seedmaxSurv = 0.99, NumMats = 100){
  # possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurv & params$survival > minSurv & 
  #                                                         params$age == stage2],] 
  # possibleParams <- possibleParams[possibleParams$survival < seedmaxSurv & possibleParams$age == 1,]
  # possibleParams <- possibleParams[possibleParams$survival > seedminSurv & possibleParams$age == 1,]
  # params1 <- possibleParams
  # vitalrates <- IteroStasisGrowth(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, params1 = params1)
  
  
  Stasis <- 
  
  mats <- lapply(sample(1:nrow(vitalrates), size = NumMats, replace = TRUE), function(i){
    S <- vitalrates[i,2]
    G <- vitalrates[i,4] # should be germination instead of growth
    S2 <- vitalrates[i,3]
    f <- itero_fecundsurv(S2) # Seed production that either goes into the soil seed bank or to rep # loss to dispersal and death
    Fecund <- rdirichlet(1, alpha = c(f,f,(1-S))) # a proportion goes to the seed bank, a proportion germinates, a portion dies
    # Seed produced either enters the soil seed bank, dies, or germinates 
    # Germ <- rdirichlet(1, alpha = c(f-(S+G),f-(S+G), (1-S), G)) # seed production to seed bank, seed to rep, death and germination)
    Germ <- rdirichlet(1, alpha = c((1-S), G*S)) # Seeds either die (1-S) or germinate immediately or from each age of seed bank
    # f2 <- itero_fecundsurv(S2) * Germ[2]
    t_ij <- matrix(c(S, Germ[2],
                     f*Fecund[1], f*Fecund[2]),
                   nrow = 2)
    (lambda1 <- lambda(t_ij))
    (e_ij <- popbio::elasticity(t_ij))
    survivalElast <- sum(e_ij[1])
    growthElast <- sum(e_ij[2])
    fecundElast <- sum(e_ij[c(3,4)])
    listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                          gentim = generation.time(t_ij, r=c(1,2), c=2),
                                          netrep = net.reproductive.rate(t_ij, r=c(1,2), c=2),
                                          AgeRep = stage1+1, lifelength = stage2), params1)
  })
  return(mats)
}

# MPM_annual <- function(maxSurv = 0.8, minSurv = 0.0001, minGerm = 0.01, maxGerm = 0.8, 
#                        persSeedBank = 3, NumMats = 100){
#   S <- runif(1,minGerm,maxGerm) # in the soil seed bank
#   G <- runif(1, minSurv, maxSurv)
#   D <- runif(1, 1-maxSurv, 1-minSurv)
#   SGD <- rdirichlet(NumMats, alpha = c(S*persSeedBank, G, D))
#   S2 <- runif(NumMats, minSurv, maxSurv)
#   (f <- semel_fecundsurv(S2))
#   f2 <- f * S2 # recruitment to next year adults is reduced by survival to rep of seed at t-1
#   mats <- lapply(sample(1:nrow(SGD), size = NumMats, replace = TRUE), function(i){
#     S <- SGD[i,1]
#     G <- SGD[i,2]
#     S2 <- S2[i]
#     f <- f[i]
#     f2 <- f2[i]
#     
#     t_ij <- matrix(c(S, G,
#                      f, f2),
#                    nrow = 2)
#   (lambda1 <- lambda(t_ij))
#   (e_ij <- popbio::elasticity(t_ij))
#   survivalElast <- sum(e_ij[1])
#   growthElast <- sum(e_ij[2])
#   fecundElast <- sum(e_ij[c(3,4)])
#   listout <- list(t_ij,e_ij, 
#                   data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
#                              gentim = generation.time(t_ij, r=c(1,2), c=2), 
#                              netrep = net.reproductive.rate(t_ij, r=c(1,2), c=2),
#                              AgeRep = stage2, lifelength = stage2))
#   
#   })
#   
#   return(mats)
# }
