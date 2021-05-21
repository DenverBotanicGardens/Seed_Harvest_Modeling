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

#' @examples 
#' @describeIn MPM_itero
 
# Type 1 for short lived/fast, type III for long lived/slow via Salguero-Gomez et al. 2016
# params = paramsAll_typeIII,
stage1 <- 10
stage2 <- stage1+1
itero = FALSE
fast = FALSE
 
MPM_iterosemel <- function(itero = TRUE, fast = TRUE, stage1, stage2, lambdarange = c(0.8,1.2), nummats = 100){
  if(fast){
    # paramsAll type I: x is age Fujiwara and Diaz-Lopez 2017; hazard is h(x) = alpha2*exp(beta2*x); 
    # exponentially increasing risk of mortality with age, risk due to aging
    possibleParams <- paramsAll[paramsAll$survival < 0.8 & paramsAll$age == stage1+1,]
    vitalrates <- do.call(rbind, mapply(function(a2, b2) FastVitalRates(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, alpha2 = a2, beta2 = b2),
                                        possibleParams$a2, possibleParams$b2, SIMPLIFY = FALSE))
  } else {
    # Slow paramsAll_typeIII
    if(itero){
      possibleParams <- paramsAll_typeIII[paramsAll_typeIII$a1 %in% paramsAll_typeIII$a1[paramsAll_typeIII$survival < 0.8 & paramsAll_typeIII$age == stage1+1],]
      a1 <- unique(possibleParams$a1)
    } else {
        possibleParams <- paramsAll_typeIII[paramsAll_typeIII$a1 %in% paramsAll_typeIII$a1[paramsAll_typeIII$survival < 0.99 & 
                                                                                             paramsAll_typeIII$age == 1],]
        possibleParams <- paramsAll_typeIII[paramsAll_typeIII$a1 %in% paramsAll_typeIII$a1[paramsAll_typeIII$survival < 0.8 & paramsAll_typeIII$age == stage1+1],]
      }
    vitalrates <- do.call(rbind, lapply(a1, function(a1){
      SlowVitalRates(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, alpha1 = a1)
    }))
  }
  
  start <- 1
  end <- nrow(vitalrates)
  lambdavec <- 1
  i_list <- c()
  
  while(length(i_list)<=nummats){
    i <- sample(start:end, 1)
    S <- vitalrates$muStasis1[i]
    G <- vitalrates$muGrowth[i]
    S2 <- vitalrates$muStasis2[i]
    if(itero){
      f <- itero_fecundsurv(S2) * vitalrates$muSurv1[i]
      t_ij <- matrix(c(S, G,
                       f, S2),
                     nrow = 2)
    } else {
      f <- semel_fecundsurv(S2) * vitalrates$muSurv1[i]
      t_ij <- matrix(c(S, G,
                       f, 0),
                     nrow = 2)
    }

    (lambda1 <- lambda(t_ij))
    if(lambda1 >= lambdarange[1] & lambda1 <= lambdarange[2]){
      i_list <- c(i_list,i)
      }
    } # end while building index list
  
  mats <- unlist(lapply(i_list, function(i){
      if(lambda1 > lambdarange[1] & lambda1 < lambdarange[2]){ # } & generation.time(t_ij) <= stage2){
        S <- vitalrates$muStasis1[i]
        G <- vitalrates$muGrowth[i]
        S2 <- vitalrates$muStasis2[i]
        if(itero){
          f <- itero_fecundsurv(S2) * vitalrates$muSurv1[i]
        } else {
          f <- semel_fecundsurv(S2) * vitalrates$muSurv1[i]
        }
        if(itero){
        t_ij <- matrix(c(S, G,
                         f, S2), nrow = 2)
        } else {
          t_ij <- matrix(c(S, G,
                           f, 0), nrow = 2)
        }
        (lambda1 <- lambda(t_ij))
        
          D <- max((1-(S+G)), 0.1)
          diri <- rdirichlet(nummats, c(S,G,D))
          mu <- S2
          sig2 <- 0.01
          alpha1 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
          beta1 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
          # when alpha1 and or beta1 is negative, get NA for rbeta
          if(alpha1 < 0 && beta1 < 0){
            alpha1 <- 0.01
            beta1 <- 0.01
          }
            S2beta <- sapply(rbeta(nummats,alpha1, beta1), function(rb) min(rb,0.79))
            mu <- vitalrates$muSurv1[i]
            alpha2 <- ((mu^2) - (mu^3) - (mu*sig2))/sig2
            beta2 <- (mu - 2*mu^2 + mu^3 - sig2 + mu*sig2)/sig2
            if(itero){
              f <- itero_fecundsurv(S2beta) * rbeta(nummats, alpha2, beta2)
            } else {
              f <- semel_fecundsurv(S2beta) * rbeta(nummats, alpha2, beta2)
            } 
            listout <- lapply(1:nummats, function(x){
              if(itero){
                t_ij <- matrix(c(diri[x,1], diri[x,2],
                                 f[x], S2beta[x]), nrow = 2)
              } else {
                t_ij <- matrix(c(diri[x,1], diri[x,2],
                                 f[x], 0), nrow = 2)
              }
              if(is.na(generation.time(t_ij))) print(t_ij) # print(paste("alpha1",alpha1, "beta1", beta1, "mu",mu))
              if(!is.na(generation.time(t_ij))){
                # if(generation.time(t_ij) <= stage2 & generation.time(t_ij) > stage1){
                if(!is.infinite(generation.time(t_ij)) & generation.time(t_ij) < stage2*3){
                  (e_ij <- popbio::elasticity(t_ij))
                  survivalElast <- sum(e_ij[which(generic_mat == "L")])
                  growthElast <- sum(e_ij[which(generic_mat == "G")])
                  fecundElast <- sum(e_ij[which(generic_mat == "F")])
                  
                  list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1,
                                             gentim = generation.time(t_ij),
                                             lifeexpectancy = sum(fundamental.matrix(t_ij, r = 1, c=2)$meaneta),
                                             netrep = net.reproductive.rate(t_ij),
                                             AgeRep = stage1+1, lifelength = stage2, 
                                             survparams = if(fast){
                                               possibleParams$a2b2[i]
                                             } else {
                                               a1[i]}
                                             ))
                  } # end if generation time is not infinite and not huge; from between stage 1 and stage 2
                } # end if generation.time worked or not
              }) # end making list for stochasticity
          listout
          } # else { # end if lambda and gen time before adding stochasticity
      
      # } # end first if generation.time worked before taking and adding stochasticity
    }), recursive=FALSE) # end mats
  return(mats[!sapply(mats, is.null)])
}



#' @describeIn MPM_annual One seed bank stage and one above ground reproductive stage

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
                                          AgeRep = 1, lifelength = 1, survparams = NA ))
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

