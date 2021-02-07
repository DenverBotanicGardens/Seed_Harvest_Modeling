library(popbio)
rm(list=ls())

# Rare plants from Salguero-Gomez et al. 2016
# Keyfitz' entropy < 1, Type I, K-selected species, low juvenile mortality with most individuals living to an old age. 



#Trade-off between survival and fecundity from Takada and Kawai (2020) 
# f=-1.25(s+4.6)^2+36.45+ error

# ------------------------------------------------ constrain lambda ---------------------------
# functions 
# Takada and Kawai 2020; s would be s_44
fecunditysurvival <- function(s, errorSD = 0.01){ 
  out <- -1.25*((s+4.6)^2) + 36.45 + rnorm(1, 0 ,errorSD)  
  out
}

# Convex for semelparous - see Farrell 2020 and Takada and Kawai 2020
semel_fecundsurv <- function(s){
  f <- -1.25*((s + 4.6)^2) + 36.45
}

# Concave for iteroparous from Takada and Kawai 2020
itero_fecundsurv <- function(s){
  f = 12.5*((s - 0.9)^2) - 0.125
}

#Hazard
# h(x) <- alpha3*exp((-beta3 * x))

# Fujiwara and Diaz-Lopez 2017 
# Do the x's represent the median age of the stage class??? I think yes
survivalTypeIII <- function(alpha1,x,alpha3,beta3){

  exp(-alpha1*x) * exp((alpha3/beta3)*(1-exp(beta3 * x)))
}

# not really
survivalTypeI <- function(alpha2, beta2, x){
  exp((alpha2/beta2)*(1-(exp(beta2*x))))
}

hazard <- function(alpha2,beta2,x){
  alpha2*exp(beta2*x)
}

## ------------------ Fast as years 1:4 where reproductive when age 4------------------
parametersfast <- do.call(rbind, lapply(seq(0.0,0.6,by=0.1), function(a1){ # a1 < 0.2
  outa3 <- do.call(rbind, lapply(seq(0.05,0.6,by=0.1), function(a3){ # a3 < 0.3
    outb3 <- do.call(rbind, lapply(seq(0.05,0.3,by=0.1), function(b3){ # b3 < 0.6
      surv <- survivalTypeIII(alpha1 = a1, 1:4, alpha3 = a3, beta3 = b3)
      data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
    }))
    outb3
  }))
  outa3
}))

parametersslow <- do.call(rbind, lapply(seq(0.0,0.6,by=0.1), function(a1){ # a1 < 0.2
  outa3 <- do.call(rbind, lapply(seq(0.05,0.6,by=0.1), function(a3){ # a3 < 0.3
    outb3 <- do.call(rbind, lapply(seq(0.05,0.3,by=0.1), function(b3){ # b3 < 0.6
      surv <- survivalTypeIII(alpha1 = a1, c(1,4,8,12), alpha3 = a3, beta3 = b3)
      data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
    }))
    outb3
  }))
  outa3
}))

paramlist_slow <- parametersslow[parametersslow$stage == "x4" & parametersslow$survival > 0.1,]
slowmin <- aggregate(survival ~ stage, min, data = parametersslow[interaction(parametersslow[,1:3]) %in% interaction(paramlist_slow[,1:3]),])
slowmax <- aggregate(survival ~ stage, max, data = parametersslow[interaction(parametersslow[,1:3]) %in% interaction(paramlist_slow[,1:3]),])

# Fast should have lower adult survival to increase fecundity
paramlist_fast <- parametersfast[parametersfast$stage == "x4" & parametersfast$survival > 0.01 & !is.na(parametersfast$survival) &
                                   parametersfast$survival < 0.2,]
fastmin <- aggregate(survival ~ stage, min, data = parametersfast[interaction(parametersfast[,1:3]) %in% interaction(paramlist_fast[,1:3]),])
fastmax <- aggregate(survival ~ stage, max, data = parametersfast[interaction(parametersfast[,1:3]) %in% interaction(paramlist_fast[,1:3]),])

fecunditysurvival(slowmin$survival[slowmin$stage == "x4"])
fecunditysurvival(slowmax$survival[slowmax$stage == "x4"]) 
fecunditysurvival(fastmin$survival[fastmin$stage == "x4"]) 
fecunditysurvival(fastmax$survival[fastmax$stage == "x4"]) 


# ------------------- parameters for slow where its 1:3 and 7 representing the median time to reproductive matruity
parametersslow2 <- do.call(rbind, lapply(seq(0.0,0.6,by=0.1), function(a1){ # a1 < 0.2
  outa3 <- do.call(rbind, lapply(seq(0.05,0.6,by=0.1), function(a3){ # a3 < 0.3
    outb3 <- do.call(rbind, lapply(seq(0.05,0.3,by=0.1), function(b3){ # b3 < 0.6
      surv <- survivalTypeIII(alpha1 = a1, c(1,2,3,7), alpha3 = a3, beta3 = b3)
      data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
    }))
    outb3
  }))
  outa3
}))


paramlist_slow2 <- parametersslow2[parametersslow2$stage == "x4" & parametersslow2$survival > 0.2,]
slowmin <- aggregate(survival ~ stage, min, data = parametersslow2[interaction(parametersslow2[,1:3]) %in% interaction(paramlist_slow2[,1:3]),])
slowmax <- aggregate(survival ~ stage, max, data = parametersslow2[interaction(parametersslow2[,1:3]) %in% interaction(paramlist_slow2[,1:3]),])

# ------------------ Iteroparous fast ------------------------
# only progressive, no retrogressive growth; r_ij = 0
MPM_iterofast <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_fast),1)
  s_s <- survivalTypeIII(alpha1 = paramlist_fast[i,1],1:4, alpha3 = paramlist_fast[i,2], beta3 = paramlist_fast[i,3] )
  f <- fecunditysurvival(s_s[4])
  t_t <- survivalTypeI(alpha2 = 0.3,beta2 = 0.1, 1:4)
  t_t <- t_t/sum(t_t)
  t_ij <- matrix(c(0, s_s[1]*sum(t_t[1:2]), s_s[1]*sum(t_t[3]),  s_s[1]*sum(t_t[4]),
                   0, s_s[2]*sum(t_t[1:2]), s_s[2]*sum(t_t[3]),  s_s[2]*sum(t_t[4]),
                   0, 0,            s_s[3]*sum(t_t[1:3]),  s_s[3]*sum(t_t[4]),
                   f, 0,            0,             s_s[4]), 
                 nrow = 4)
  lambda1 <- lambda(t_ij)
  lambdarange <- abs(1-rnorm(1, 1, 0.1))
  # plot(seq(0.51,1.5,by=0.01), dnorm(seq(0.51,1.5,by=0.01), 1, 0.1))
  if(lambda1 > (1+lambdarange)) { # allow more or less variability by speed of life, more variable for fastest
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(0.1,0.5, by = 0.1)){
      t_ij[mostelastic] <- (1-k) * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 <= (1+lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda(t_ij) < (1-lambdarange)){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(1.1,1.5, by=0.1)){
      t_ij[mostelastic] <- k * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 >= (1-lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  return(t_ij)
}

MPMs_itfast <- lapply(1:100, function(x) MPM_iterofast())
generation.time(mean(MPMs_itfast))
hist(unlist(lapply(MPMs_itfast, function(x) lambda(x))), xlab = "lambda", main = "Iteroparous Fast")


# --------------------- Iteroparous slow -------------------
MPM_iteroslow <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_slow2),1)
  s_s <- survivalTypeIII(alpha1 = paramlist_slow2[i,1],c(1:3,7), alpha3 = paramlist_slow2[i,2], beta3 = paramlist_slow2[i,3] )
  t_t <- survivalTypeI(alpha2 = 0.2, beta2 = .01, 1:4)
  # plot(1:4, t_t)
  t_t <- t_t/sum(t_t)
  # t_t <- s_s/sum(s_s)
  f <- fecunditysurvival(s_s[4])
  t_ij <- matrix(c(0, s_s[1]*(sum(t_t[1:3])),0,                    0,
                   0, s_s[2]*sum(t_t[1:2]),  s_s[2]*sum(t_t[3:4]),0,
                   0, 0,                      s_s[3]*sum(t_t[1:2]),s_s[3]*sum(t_t[3:4]),
                   f, 0,            0,             s_s[4]), 
                 nrow = 4)
  lambda1 <- lambda(t_ij)
  lambdarange <- abs(1-rnorm(1, 1, 0.05))
  if(lambda1 > (1+lambdarange)) { # allow more or less variability by speed of life, more variable for fastest
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(0.1,0.5, by = 0.01)){
      t_ij[mostelastic] <- (1-k) * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 <= (1+lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda(t_ij) < (1-lambdarange)){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    print(i)
    for(k in seq(1.1,1.5, by=0.01)){
      t_ij[mostelastic] <- k * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 >= (1-lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  return(t_ij)
}

MPMs_itslow <- lapply(1:100, function(x) MPM_iteroslow())
generation.time(mean(MPMs_itslow))
hist(unlist(lapply(MPMs_itslow, function(x) lambda(x))), xlab = "lambda", main = "Iteroparous Slow")


# --------------------- semelparous fast  ------------------------
MPM_semelfast <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_fast),1)
  s_s <- survivalTypeIII(alpha1 = paramlist_fast[i,1],1:4, alpha3 = paramlist_fast[i,2], beta3 = paramlist_fast[i,3] )
  # let the same Type III decay happen in transitions; first is most likely, sharp drop off
  t_t <- survivalTypeI(alpha2 = 0.2, beta2 = .01, 1:4)
  # plot(1:4, t_t)
  t_t <- t_t/sum(t_t)
  f <- fecunditysurvival(s_s[4])
  t_ij <- matrix(c(0, s_s[1]*(sum(t_t[1:2])),s_s[1]*t_t[3],       s_s[1]*t_t[4],
                   0, s_s[2]*sum(t_t[1:2]),  s_s[2]*sum(t_t[3]),  s_s[2]*t_t[4],
                   0, 0,                      s_s[3]*sum(t_t[1:3]),s_s[3]*sum(t_t[4]),
                   f, 0,            0,             0), 
                 nrow = 4)
  lambda1 <- lambda(t_ij)
  lambdarange <- abs(1-rnorm(1, 1, 0.1))
  # plot(seq(0.51,1.5,by=0.01), dnorm(seq(0.51,1.5,by=0.01), 1, 0.1))
  if(lambda1 > (1+lambdarange)) { # allow more or less variability by speed of life, more variable for fastest
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    # secondmostelastic <- which(popbio::elasticity(t_ij) == sort(popbio::elasticity(t_ij), TRUE)[2])
    for(k in seq(0.1,0.5, by = 0.1)){
      t_ij[mostelastic] <- (1-k) * t_ij[mostelastic]
      # t_ij[secondmostelastic] <- (1-k) * t_ij[secondmostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 <= (1+lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda(t_ij) < (1-lambdarange)){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(1.1,1.5, by=0.1)){
      t_ij[mostelastic] <- k * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 >= (1-lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  return(t_ij)
}

MPMs_sefast <- lapply(1:100, function(x) MPM_semelfast())
generation.time(mean(MPMs_sefast))
hist(unlist(lapply(MPMs_sefast, function(x) lambda(x))), xlab = "lambda", main = "Semelparous Fast")


# --------------------- Semelparous slow -------------------
MPM_semelslow <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_slow2),1)
  s_s <- survivalTypeIII(alpha1 = paramlist_slow2[i,1],c(1:3,7), alpha3 = paramlist_slow2[i,2], beta3 = paramlist_slow2[i,3] )
  f <- fecunditysurvival(s_s[4]) # remove retrogressive growth
  t_t <- survivalTypeI(alpha2 = 0.001, beta2 = .001, 1:4)
  # plot(1:4, t_t)
  t_t <- t_t/sum(t_t)
  t_ij <- matrix(c(0, s_s[1],                     0,                    0,
                   0, s_s[2]*sum(t_t[c(1:2,4)]),  s_s[2]*sum(t_t[3]),   0,
                   0, 0,                          s_s[3]*sum(t_t[1:3]), s_s[3]*sum(t_t[4]),
                   f, 0,            0,             0),
                 nrow = 4)
  (lambda1 <- lambda(t_ij))
  lambdarange <- abs(1-rnorm(1, 1, 0.05))
  if(lambda1 > (1+lambdarange)) { # allow more or less variability by speed of life, more variable for fastest
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(0.1,0.5, by = 0.1)){
      t_ij[mostelastic] <- (1-k) * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 <= (1+lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda(t_ij) < (1-lambdarange)){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(1.1,1.5, by=0.1)){
      t_ij[mostelastic] <- k * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 >= (1-lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  return(t_ij)
}

MPMs_seslow <- lapply(1:100, function(x) MPM_semelslow())
generation.time(mean(MPMs_seslow))
hist(unlist(lapply(MPMs_seslow, function(x) lambda(x))), xlab = "lambda", main = "Semelparous Slow")


#--------------------- Parameters for annuals --------------------
paramlist_annuals <- do.call(rbind, lapply(seq(0.0,0.6,by=0.1), function(a1){ # a1 < 0.2
  outa3 <- do.call(rbind, lapply(seq(0.05,0.6,by=0.1), function(a3){ # a3 < 0.3
    outb3 <- do.call(rbind, lapply(seq(0.05,0.3,by=0.1), function(b3){ # b3 < 0.6
      surv <- survivalTypeIII(alpha1 = a1, c(1,2,3,7), alpha3 = a3, beta3 = b3)
      data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
    }))
    outb3
  }))
  outa3
}))
# --------------------- Annuals ------------------------
MPM_annual <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_fast),1)
  s_s <- survivalTypeIII(alpha1 = paramlist_fast[i,1],1:4, alpha3 = paramlist_fast[i,2], beta3 = paramlist_fast[i,3] )
  f <- fecunditysurvival(s_s[4])
  t_t <- survivalTypeI(alpha2 = 0.3,beta2 = 0.1, 1:4)
  t_t <- t_t/sum(t_t)
  t_ij <- matrix(c(0, s_s[1]*(2/3), 0,             s_s[1]*(1/3),
                   0, 0,            s_s[2]*(2/3),  s_s[2]*(1/3),
                   0, 0,            0,             s_s[3],
                   f, 0,            0,             0), 
                 nrow = 4)
  lambda1 <- lambda(t_ij)
  lambdarange <- abs(1-rnorm(1, 1, 0.1))
  # plot(seq(0.51,1.5,by=0.01), dnorm(seq(0.51,1.5,by=0.01), 1, 0.1))
  if(lambda1 > (1+lambdarange)) { # allow more or less variability by speed of life, more variable for fastest
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    # secondmostelastic <- which(popbio::elasticity(t_ij) == sort(popbio::elasticity(t_ij), TRUE)[2])
    for(k in seq(0.1,0.5, by = 0.1)){
      t_ij[mostelastic] <- (1-k) * t_ij[mostelastic]
      # t_ij[secondmostelastic] <- (1-k) * t_ij[secondmostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 <= (1+lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda(t_ij) < (1-lambdarange)){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    for(k in seq(1.1,1.5, by=0.1)){
      t_ij[mostelastic] <- k * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 >= (1-lambdarange)) break
      mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  return(t_ij)
}

MPMs_annual <- lapply(1:100, function(x) MPM_annual())
generation.time(mean(MPMs_annual))
hist(unlist(lapply(MPMs_annual, function(x) lambda(x))), xlab = "lambda", main = "Annual")




MPMs <- list(annual = MPM_annual, semelfast = MPM_semelfast, iterofast = MPM_iterofast,
             semelslow = MPM_semelslow, iteroslow = MPM_iteroslow)


# Testing
# Nxstages <- c(rep(NA,3),1)
# StartPopSize <- 10
# MatrixType <- "annual"
# repl <- 1
# rm(Nxstages); rm(StartPopSize); rm(MatrixType); rm(repl)

# ------------------------ Simulation of virtual species -----------------------------
# Modified from VirtualSpeciesSimulationLength_v2
virtualPVA <- function(Nxstages = rep(1,4), StartPopSize, MatrixType){
  dfout <- do.call(rbind,lapply(1:100, function(repl){
    # Take a sample of size at least 5, remove fecundity for TMx_list 
    Mx_sample <- lapply(1:100, function(i) MPMs[[grep(MatrixType, c("annual","semelfast","iterofast","semelslow","iteroslow"))]]())
    # calculate generation time to scale by generation time or set number of years for projection length
    gentim <- popbio::generation.time(mean(Mx_sample), r=1, c=4)
    projlength <- max(gentim*3, 100)
    # scale the above ground growth to the starting population size
    Nx <- stable.stage(mean(Mx_sample))
    Nx_scale <- Nx[!is.na(Nxstages)]
    vec1 <- matrix(floor(Nx*(StartPopSize/sum(Nx_scale))), ncol = 1)
    # Initilize for loop
    Extant <- c()
    yr <- c()
    mats <- list()
    popsz <- c()
    Time2Ext <- NA
    for(i in 1:projlength){
      yr[i] <- i
      mats[[i]] <- (sample(Mx_sample,1))[[1]]
      vec1 <- floor(mats[[i]]%*%vec1)
      popsz[i] <- sum(vec1[!is.na(Nxstages)])
      Extant[i] <- if(popsz[i]<1) 0 else 1
      if(Extant[i] == 0) Time2Ext <- i
      if(popsz[i] == 0) break
    }
    # A list of 2 or more, will be NA until second year
    if(length(popsz) > 2){
      if(length(popsz) < 6){
        lambdas <- do.call(rbind,lapply(seq(2,length(mats)), function(x){
          stochLambda <- popbio::stoch.growth.rate(mats[1:x], verbose = FALSE)
          data.frame(Year = x, Tulapprox = stochLambda$approx, LogGrowthsim = stochLambda$sim,
                     LogGrowthSimLowerCI = stochLambda$sim.CI[1], LogGrowthSimUpperCI = stochLambda$sim.CI[2])
        }))
      } else {
        lambdas <- do.call(rbind,lapply(seq(5,length(popsz),by = 5), function(x){
          stochLambda <- popbio::stoch.growth.rate(mats[1:x], verbose = FALSE)
          data.frame(Year = x, Tulapprox = stochLambda$approx, LogGrowthsim = stochLambda$sim,
                     LogGrowthSimLowerCI = stochLambda$sim.CI[1], LogGrowthSimUpperCI = stochLambda$sim.CI[2])
        }))
      }
    } else {
      lambdas <- data.frame(Year = yr, Tulapprox = unlist(lapply(mats, function(x) lambda(x))),
                            LogGrowthsim = NA, LogGrowthSimLowerCI = NA, LogGrowthSimUpperCI = NA)
    }
    dfMPM <- data.frame(PopSize = popsz, Year = yr, MatType = MatrixType, GenTime = gentim,
                        Extant = Extant, StPopSz = StartPopSize, Time2Extinction = Time2Ext, 
                        ProjectionLength = projlength, Replicate = repl)
    out <- merge(dfMPM, lambdas, by = "Year", all = TRUE)
    out
  }))
  dfout
}

# test 
Nxstages <- c(NA,NA,NA,1)
StartPopSize <- 10
MatrixType <- "annual"
# ------------------- Simulation runs ----------------------------
pathstartVirtual <-  "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Simulation_output_virt_lambda_1/"

annuals <- do.call(rbind,lapply(c(10,50,100,500), function(StPopSz){
  out <- virtualPVA(Nxstages = c(NA,NA,NA,1), StartPopSize = StPopSz, MatrixType = "annual" )
  save(out, file =  paste(pathstartVirtual,"annuals",StPopSz,".Rdata", sep=""))
  out
}))

iterofast <- do.call(rbind,lapply(c(10,50,100,500), function(StPopSz){
  out <- virtualPVA(StartPopSize = StPopSz, MatrixType = "iterofast")
  save(out, file =  paste(pathstartVirtual,"iterofast",StPopSz,".Rdata", sep=""))
  out
}))

# error eigen(Abar) : non-square matrix in 'eigen'; Abar likely means mean(A) 
# When first thing went extinct, then mats[2:1] gave the one matrix and NULL 
iteroslow <- do.call(rbind,lapply(c(10,50,100,500), function(StPopSz){
  out <- virtualPVA(StartPopSize = StPopSz, MatrixType = "iteroslow")
  save(out, file =  paste(pathstartVirtual,"iteroslow",StPopSz,".Rdata", sep=""))
  out
}))

semelfast <- do.call(rbind,lapply(c(10,50,100,500), function(StPopSz){
  out <- virtualPVA(StartPopSize = StPopSz, MatrixType = "semelfast")
  save(out, file =  paste(pathstartVirtual,"semelfast",StPopSz,".Rdata", sep=""))
  out
}))

semelslow <- do.call(rbind,lapply(c(10,50,100,500), function(StPopSz){
  out <- virtualPVA(StartPopSize = StPopSz, MatrixType = "semelslow")
  save(out, file =  paste(pathstartVirtual,"semelslow",StPopSz,".Rdata", sep=""))
  out
}))
