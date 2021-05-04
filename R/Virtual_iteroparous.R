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

# MPM_itero <- function(params, stage1, stage2, maxSurvVeg = 0.99, maxSurvRep = 0.8, minSurv = 0.001, NumMats = 100){
#   possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurvRep & 
#                                                           params$survival > minSurv & params$age == stage1+1],]
#   possibleParams <- possibleParams[possibleParams$a2b2 %in% 
#                                      possibleParams$a2b2[possibleParams$survival < maxSurvVeg & possibleParams$age == 1],]
#   params1 <- possibleParams
#   vitalrates <- IteroStasisGrowth(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, params1 = params1)
#   mats <- lapply(sample(1:nrow(vitalrates), size = NumMats, replace = TRUE), function(i){
#     S <- vitalrates[i,2]
#     G <- vitalrates[i,4]
#     S2 <- vitalrates[i,3]
#     f <- itero_fecundsurv(S2) * vitalrates[i,1]
#     t_ij <- matrix(c(S, G,
#                      f, S2),
#                    nrow = 2)
#     if(colSums(t_ij)[1] > 1){ 
#       print(colSums(t_ij))
#       print(lambda(t_ij))
#       }
#     (lambda1 <- lambda(t_ij))
#     (e_ij <- popbio::elasticity(t_ij))
#     survivalElast <- sum(e_ij[which(generic_mat == "L")])
#     growthElast <- sum(e_ij[which(generic_mat == "G")])
#     fecundElast <- sum(e_ij[which(generic_mat == "F")])
#     listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
#                                         gentim = generation.time(t_ij), 
#                                         netrep = net.reproductive.rate(t_ij),
#                                         AgeRep = stage1+1, lifelength = stage2), params1)
#   })
#   return(mats)
# }

#' @examples 
#' MPM_itero()
#' 
#' @describeIn MPM_semel

# test
# params <- paramsAll
# maxSurv <- 0.4
# minSurv <- 0.0001
# stage1 <- 1 
# stage2 <- stage1+1

# MPM_semel <- function(params, stage1, stage2, maxSurv = 0.8, minSurv = 0.0001, vegminSurv = 0.9, vegmaxSurv = 0.99, NumMats = 100){
#   possibleParams <- params[params$a2b2 %in% params$a2b2[params$survival < maxSurv & params$survival > minSurv & 
#                                                           params$age == stage2],] 
#   possibleParams <- possibleParams[possibleParams$survival < vegmaxSurv & possibleParams$age == 1,]
#   possibleParams <- possibleParams[possibleParams$survival > vegminSurv & possibleParams$age == 1,]
#   params1 <- possibleParams
#   vitalrates <- IteroStasisGrowth(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, params1 = params1)
#   mats <- lapply(sample(1:nrow(vitalrates), size = NumMats, replace = TRUE), function(i){
#     S <- vitalrates[i,2]
#     G <- vitalrates[i,4]
#     S2 <- vitalrates[i,3]
#     f <- semel_fecundsurv(S2) * vitalrates[i,1]
#     t_ij <- matrix(c(S, G,
#                      f, 0),
#                    nrow = 2)
#     if(colSums(t_ij)[1] > 1){ 
#       print(colSums(t_ij))
#       print(lambda(t_ij))
#     }
#     (lambda1 <- lambda(t_ij))
#     (e_ij <- popbio::elasticity(t_ij))
#     survivalElast <- sum(e_ij[which(generic_mat == "L")])
#     growthElast <- sum(e_ij[which(generic_mat == "G")])
#     fecundElast <- sum(e_ij[which(generic_mat == "F")])
#     listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
#                                           gentim = generation.time(t_ij, r = 1, c = 2), 
#                                           lifeexpectancy = sum(fundamental.matrix(t_ij, r = 1, c = 2)$meaneta),
#                                           netrep = net.reproductive.rate(t_ij),
#                                           AgeRep = stage1+1, lifelength = stage2), params1)
#   })
#   return(mats)
# }


#' @describeIn MPM_annual One seed bank stage and one above ground reproductive stage

# germRatemu = 0.1; germRatesig2 = 0.01; seedSurvmu = 0.1; seedSurvsig2 = 0.01; NumMats = 100; survS2mu = 0.5; survS2sig2 = 0.01
# rm(germRatemu);rm(germRatesig2);rm(seedSurvmu);rm(seedSurvsig2);rm(survS2mu);rm(survS2sig2)

MPM_annual <- function(germRatemu = 0.01, germRatesig2 = 0.01, seedSurvmu = 0.01, seedSurvsig2 = 0.01,
                       survS2mu = 0.5, survS2sig2 = 0.01, NumMats = 100){
  alpha <- function(mu, sig2){
    ((mu^2) - (mu^3) - (mu*sig2))/sig2
  }
  beta <- function(mu, sig2){
    (mu - (2 * mu^2) + (mu^3) -(sig2) + (mu * sig2))/sig2
  }
  seedPersist <- rbeta(NumMats, shape1 = alpha(seedSurvmu, seedSurvsig2), shape2 = beta(seedSurvmu, seedSurvsig2))
  seedGerm <- rbeta(NumMats, shape1 = alpha(germRatemu, germRatesig2), shape2 = beta(germRatemu, germRatesig2))
  seedDeath <- rbeta(NumMats, shape1 = alpha((1-seedSurvmu), seedSurvsig2), shape2 = beta((1-seedSurvmu), seedSurvsig2))
  # soil seed bank, death, germination
  Stasis <- rdirichlet(NumMats, alpha = matrix(c(seedPersist,
                                                 seedDeath,
                                                 seedGerm), nrow = NumMats))
  
  mats <- lapply(1:NumMats, function(i){
    S <- Stasis[i,1]
    G <- Stasis[i,ncol(Stasis)] # germination
    S2 <- rbeta(1, shape1 = alpha(survS2mu,survS2sig2), shape2 = beta(survS2mu,survS2sig2))
    f <- itero_fecundsurv(S2) # Seed production that either goes into the soil seed bank or to rep # loss to dispersal and death
    Fecund <- rdirichlet(1, alpha = c(f,f,(1-S))) # a proportion goes to the seed bank, a proportion germinates, a portion dies
    t_ij <- matrix(c(S, G,
                     f*Fecund[1], f*Fecund[2]),
                   nrow = 2)
    (lambda1 <- lambda(t_ij))
    (e_ij <- popbio::elasticity(t_ij))
    survivalElast <- sum(e_ij[1])
    growthElast <- sum(e_ij[2])
    fecundElast <- sum(e_ij[c(3,4)])
    listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                          gentim = generation.time(t_ij, r=c(1,2), c=2),
                                          lifeexpectancy = fundamental.matrix(t_ij, r = c(1,2), c=2)$meaneta[2],
                                          netrep = net.reproductive.rate(t_ij, r=c(1,2), c=2),
                                          AgeRep = 1, lifelength = 1))
  })
  return(mats)
}



# Ternary::TernaryPlot(point = "up",
#                      atip = "Seed",
#                      btip = "Death",
#                      ctip = "Germination",
#                      alab = "Seed bank \u2192", blab = "Death \u2192", clab = "\u2190 germination", # with arrows
#                      grid.minor.lines = 0) 
# AddToTernary(points, Stasis, col = "red", pch = 16, cex = 1) 



#### Try with constant risk survival curves
# stage1 <- 3
# stage2 <- 5
# params <- paramsAll_typeIII
# params <- paramsAll_sub

MPM_itero <- function(params = paramsAll_typeIII, stage1, stage2, lambdarange = c(0.99,1.05)){
  possibleParams <- params[params$a1 %in% params$a1[params$survival < 0.8 & params$age == stage1+1],]
  # possibleParams <- possibleParams[possibleParams$a1 %in%
  #                                    possibleParams$a1[possibleParams$survival < maxSurvVeg & possibleParams$age == 1],]
  a1 <- unique(possibleParams$a1)

  vitalrates <- do.call(rbind, lapply(a1, function(a1){
    IteroStasisGrowth_constant(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, alpha1 = a1)
  }))
  mats <- unlist(lapply(1:nrow(vitalrates), function(i){
  # mats <- lapply(sample(1:nrow(vitalrates), size = NumMats, replace = TRUE), function(i){
  # for(i in 1:nrow(vitalrates)){
    S <- vitalrates$muStasis1[i]
    G <- vitalrates$muGrowth[i]
    S2 <- vitalrates$muStasis2[i]
    f <- itero_fecundsurv(S2) * vitalrates$muSurv1[i]
    t_ij <- matrix(c(S, G,
                     f, S2),
                   nrow = 2)
    (lambda1 <- lambda(t_ij))
    if(lambda1 > lambdarange[1] & lambda1 < lambdarange[2] & generation.time(t_ij) < stage2){
      D <- max((1-(S+G)), 0.01)
      diri <- rdirichlet(10, c(S,G,D))
      mu <- S2
      sig2 <- 0.01
      alpha1 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
      beta1 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
      S2beta <- sapply(rbeta(10,alpha1, beta1), function(rb) min(rb,0.79))
      mu <- vitalrates$muSurv1[i]
      alpha2 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
      beta2 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
      f <- itero_fecundsurv(S2beta) * rbeta(10, alpha2, beta2)
      
      listout <- lapply(1:10, function(x){
        t_ij <- matrix(c(diri[x,1], diri[x,2],
                         f[x], S2beta[x]), nrow = 2)
        if(generation.time(t_ij) < stage2 & generation.time(t_ij) > stage1){
          (e_ij <- popbio::elasticity(t_ij))
          survivalElast <- sum(e_ij[which(generic_mat == "L")])
          growthElast <- sum(e_ij[which(generic_mat == "G")])
          fecundElast <- sum(e_ij[which(generic_mat == "F")])
          
          list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                              gentim = generation.time(t_ij),
                                              lifeexpectancy = sum(fundamental.matrix(t_ij, r = 1, c=2)$meaneta),
                                              netrep = net.reproductive.rate(t_ij),
                                              AgeRep = stage1+1, lifelength = stage2, alpha1 = a1[i]))
        } # end if generation time less than stage 2
        }) # end making list for stochasticity
        listout
      } # end if lambda and gen time before adding stochasticity
    }), recursive=FALSE) # end mats
  
  return(mats[!sapply(mats, is.null)])
}

# test
# params <- paramsAll_typeIII
# head(params)
# stage1 <- 1
# stage2 <- 2

MPM_semel <- function(params, stage1, stage2, lambdarange = c(0.99, 1.05)){
  possibleParams <- params[params$a1 %in% params$a1[params$survival < 0.8 & params$age == stage1+1],]
  a1 <- unique(possibleParams$a1)
  
  vitalrates <- do.call(rbind, lapply(a1, function(a1){
    IteroStasisGrowth_constant(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, alpha1 = a1)
  }))
  mats <- unlist(lapply(1:nrow(vitalrates), function(i){
    S <- vitalrates$muStasis1[i]
    G <- vitalrates$muGrowth[i]
    S2 <- vitalrates$muStasis2[i]
    f <- semel_fecundsurv(S2) * vitalrates$muSurv1[i]
    t_ij <- matrix(c(S, G,
                     f, 0),
                   nrow = 2)
    (lambda1 <- lambda(t_ij))
    if(lambda1 > lambdarange[1] & lambda1 < lambdarange[2] & generation.time(t_ij) < stage2){
      D <- max((1-(S+G)), 0.01)
      diri <- rdirichlet(10, c(S,G,D))
      mu <- S2
      sig2 <- 0.01
      alpha1 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
      beta1 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
      S2beta <- sapply(rbeta(10,alpha1, beta1), function(rb) min(rb,0.79))
      mu <- vitalrates$muSurv1[i]
      alpha2 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
      beta2 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
      f <- semel_fecundsurv(S2beta) * rbeta(10, alpha2, beta2)
      
      listout <- lapply(1:10, function(x){
        t_ij <- matrix(c(diri[x,1], diri[x,2],
                         f[x], 0), nrow = 2)
        if(generation.time(t_ij) < stage2 & generation.time(t_ij) > stage1){
          (e_ij <- popbio::elasticity(t_ij))
          survivalElast <- sum(e_ij[which(generic_mat == "L")])
          growthElast <- sum(e_ij[which(generic_mat == "G")])
          fecundElast <- sum(e_ij[which(generic_mat == "F")])
          
          list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                     gentim = generation.time(t_ij),
                                     lifeexpectancy = sum(fundamental.matrix(t_ij, r = 1, c=2)$meaneta),
                                     netrep = net.reproductive.rate(t_ij),
                                     AgeRep = stage1+1, lifelength = stage2, alpha1 = a1[i]))
        } # end if generation time less than stage 2
      }) # end making list for stochasticity
      listout
    } # end if lambda and gen time before adding stochasticity
  }), recursive=FALSE) # end mats
  
  return(mats[!sapply(mats, is.null)])
}

 