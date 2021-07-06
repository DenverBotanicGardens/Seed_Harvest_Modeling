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
# stage1 <- 10
# stage2 <- stage1+1
# itero = FALSE
# fast = FALSE
# 
# stage1 <- 25
# stage2 <- 26
# itero <- FALSE
# fast <- FALSE 
# if slow, then for iteroparous need to restrict survival to 0.8 at transition to reproductive (stage 1 + 1) while for semelparous need to , I'd say do same for both
# Slow paramsAll_typeIII
# if(itero){
#   possibleParams <- paramsAll_typeIII[paramsAll_typeIII$a1 %in% paramsAll_typeIII$a1[paramsAll_typeIII$survival < 0.8 & paramsAll_typeIII$age == stage1+1],]
#   a1 <- unique(possibleParams$a1)
# } else {
# In Fujiwara and Diaz-Lopez early maturation is this curve so should be fast... 
# under constant risk, exponentially declining survivorship curve Type III

MPM_iterosemel <- function(itero = TRUE, fast = TRUE, stage1, stage2, lambdarange = c(0.8,1.2), nummats = 100, fastslowcutoff = 5){
  if(fast){
    # paramsAll type I: x is age Fujiwara and Diaz-Lopez 2017; hazard is h(x) = alpha2*exp(beta2*x); 
    # exponentially increasing risk of mortality with age, risk due to aging
    possibleParams <- paramsAll[paramsAll$survival < 0.8 & paramsAll$age == stage1+1,]
    vitalrates <- do.call(rbind, mapply(function(a2, b2) FastVitalRates(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, 
                                                                        alpha2 = a2, beta2 = b2),
                                        possibleParams$a2, possibleParams$b2, SIMPLIFY = FALSE))
  } else {
# moved above
    if(itero){
        possibleParams <- paramsAll_typeIII[paramsAll_typeIII$a1 %in% paramsAll_typeIII$a1[paramsAll_typeIII$survival < 0.8 & 
                                                                                             paramsAll_typeIII$age == 1],]
    } else {
      #semel
      possibleParams <- paramsAll_typeIII[paramsAll_typeIII$a1 %in% paramsAll_typeIII$a1[paramsAll_typeIII$survival < 0.93 & 
                                                                                           paramsAll_typeIII$age == 1],]
        }
        # possibleParams <- possibleParams[possibleParams$a1 %in% possibleParams$a1[possibleParams$survival < 0.8 & 
        #                                                                             possibleParams$age == stage1+1],]
        if(itero){
        possibleParams <- possibleParams[possibleParams$survival < 0.2 & # 0.8 & 
                                           possibleParams$age == stage1 + 1,]
        a1 <- possibleParams$a1
        } else {
          possibleParams <- possibleParams[possibleParams$survival < 0.8 & possibleParams$age == stage1 + 1,]
          a1 <- possibleParams$a1
        }
      # }
    vitalrates <- do.call(rbind, lapply(a1, function(a1){
      SlowVitalRates(Age_first_stage2 = stage1+1, Age_last_stage2 = stage2, alpha1 = a1)
    }))
  }
  
  start <- 1
  end <- nrow(vitalrates)
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
      if(!fast) S2 <- runif(1, 0.001, 0.01)
      f <- semel_fecundsurv(S2) * vitalrates$muSurv1[i]
      t_ij <- matrix(c(S, G,
                       f, 0),
                     nrow = 2)
    }

    (lambda1 <- lambda(t_ij))
    if(lambda1 >= lambdarange[1] & lambda1 <= lambdarange[2]){
      if(fast){
        if(generation.time(t_ij) < fastslowcutoff){
          i_list <- c(i_list,i)
        }
        } else {
          # Slow
          if(generation.time(t_ij) >= fastslowcutoff) 
            i_list <- c(i_list,i)
        }
      }
    } # end while building index list

  mats <- #unlist(
    lapply(i_list, function(i){
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
      
      e_ij <- popbio::elasticity(t_ij)
      survivalElast <- sum(e_ij[which(generic_mat == "L")])
      growthElast <- sum(e_ij[which(generic_mat == "G")])
      fecundElast <- sum(e_ij[which(generic_mat == "F")])
      
      listout <- list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, 
                                            lam = lambda(t_ij), gentim = generation.time(t_ij),
                                            lifeexpectancy = sum(fundamental.matrix(t_ij, r = 1, c=2)$meaneta),
                                            netrep = net.reproductive.rate(t_ij),
                                            AgeRep = stage1+1, lifelength = stage2, 
                                            survparams = if(fast){
                                              possibleParams$a2b2[i]
                                              } else {
                                                a1[i]
                                                }
                                            ))
      listout
      }) # , recursive=FALSE) # end mats
  # return(mats[!sapply(mats, is.null)])
  return(mats)
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

