rm(list=ls())
library(dplyr)
library(popbio)
library(MuMIn)
library(binr)
library(matrixStats)
require(AICcmodavg)
library(prism)
library(raster)
library(lme4)
library(ggplot2)
library(patchwork)
# install.packages("roxygen2")

library(stringr)
setwd("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/")

lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))
  nx[1,]=0
  jj=1
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
  }
  return(jj)
}

# Seeds a population with 100 individuals in the smallest size class, removes the top row (so no more reproduction), and follows the population until all have died. Based on a list of transition matrices (with fecundity rate: Mx, and without fecundity rate included: TMx)
lifespanMPM <- function(Mx, TMx, seedlingBin = 1){
  nclasses=dim(Mx[[1]])[1]
  # add 100 individuals to the smallest non-seedbank category
  vec <- rep(0, nclasses)
  vec[seedlingBin] <- 100
  # remove all fecundity from matrices
  nx <- mapply(function(x,y) x - (x-y), Mx, TMx, SIMPLIFY = FALSE) # keep output as matrices
  # if any transition is 100%
  nx[[1]][nx[[1]] == 1] <- 0.9
  jj=1
  while (sum(vec)>1){
    vec=nx[[sample(1:length(Mx),1)]]%*%(floor(vec)) # randomly sample from the list of matrices with fecundity rate removed
    jj=jj+1
  }
  return(jj)
}
AnnualTransition <-function(){
  # above_f <- Mx_list[which(f_all >= mn_f)]  from the SeedHarvestSim function get index of one above 
  # below_f <- Mx_list[which(f_all < mn_f)]   from the SeedHarvestSim function get index of one below
  if(Yr == "Good"){
    Mx_i <- sample(above_f, 1) # randomly sample from the high fecundity transition matrices
  } else {
    Mx_i <- sample(below_f, 1)
  }
  
  # Freq[] is the likelihood to collect seed in good or bad years
  if(CollectYN == 1){  # yes collect, reduce fecuntity by the year good or bad
    fecundMx <- Mx_list[[Mx_i]] - TMx_list[[Mx_i]] # Subtract the transition matrix (Mx) by the matrix without fecundity rate (TMx) to get only fecundity rate
    fecundMx <- fecundMx - fecundMx*Inten_m[Yr,] # Reduce all fecundity rates by intensity of harvest
    nx <- fecundMx + TMx_list[[Mx_i]] # add the reduced fecundity rates to the transition matrix that lacks the fecundity rate to get new transition matrix
  } else {
    nx <- Mx_list[[Mx_i]]
  }
  K <- 100*StartPopSize  # maximum population size
  if(ddceiling == TRUE & sum(vec[!is.na(Nx_list[[1]])]) >= K*0.9){
    vec.t1 <- floor(nx%*%vec)
    vec <- floor(vec.t1/(lambda(nx))*(1+(lambda(nx)-1)*(1-(sum(vec[!is.na(Nx_list[[1]])])/K ))))  # dominant eigenvalue is lambda and the same as eigen(nx)$values[1]; what now? 1+(lambda(nx)-1) is just lambda(nx)
  } else {
    # Multiply N_t+1 = Matrix*N_t; vec comes from the SeedHarvestSim function
    vec <- floor(nx%*%vec)
  }
  list(vec, lambda(nx), projmat = nx)
}

# starting population of only above ground individuals - pop size is above ground
SeedHarvestSim <- function(Mx_list, TMx_list, Nx, Nx_list, StartPopSize, GoodBadTm, Freq = c(0,0), Int = c(0,0), 
                           TotYrs, ClusteredColl = 1, ddceiling = FALSE){
  # take the distribution of classes (likely the stable stage distribution) and scale to starting population size with integers
  # Scale to starting population size with only the aboveground, detectable stages (no seedbank or dormant); undetectable stages have NA in Nx_list
  Nx_scale <- Nx[!is.na(Nx_list[[1]])]
  # vec <- matrix(floor(Nx*(StartPopSize/sum(Nx))), ncol = 1)
  vec <- matrix(round(Nx*(StartPopSize/sum(Nx_scale)),0), ncol = 1)
  # Annual transition function needs local environment variables
  environment(AnnualTransition) <- environment() 
  
  # Seed Harvest Rates, frequency, probability of harvest in a given year and intensity or percent of seed harvested
  #     Freq
  #Good X
  #Bad  Y
  # Good vs. bad - either determine the percent per good and bad after given clusters or just keep not in clusters and
  Freq_m <- matrix(Freq, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Freq")))
  Inten_m <- matrix(Int, nrow=2, byrow = TRUE, dimnames = list(c("Good","Bad"), c("Inten")))
  
  # fecundity is the cutoff (practical sense of a good year with lots of seed vs few seed produced) but this doesn't really tell me if there is high production of seed, it's the ratio of recruitment to seed production
  f_all <- apply(mapply(function(x,y) x-y, Mx_list, TMx_list), 2, sum) # subtracting the matrices with fecundity from the matrices with all transitions leaves only the fecundity transitions, adding all remaining fecundity values
  mn_f <- mean(f_all) 
  #The index of above/equal and below average
  above_f <- which(f_all >= mn_f)
  below_f <- which(f_all < mn_f)
  
  # Initialize 
  ceil <- 10*StartPopSize
  Yr <- "Good"
  GBs <- rep(NA,TotYrs)
  Freqs <- rep(NA,TotYrs)
  Intens <- rep(NA,TotYrs)
  PopSize <- rep(NA, TotYrs)
  lambdas <- rep(NA, TotYrs)
  mats <- list()
  yrs <- 1
  
  # if clustered collection then does not depend on good or bad, will collect during a bad year, all match the larger frequency
  if(ClusteredColl > 1){
    numCollections <- max(Freq)*TotYrs
    numclusters <- floor(numCollections/ClusteredColl) # might underestimate
    i <- 0
    collectionyr_ix <- c()
    notavailable <- c()
    while(i < numclusters){
      available <- setdiff((1:(TotYrs-ClusteredColl)),notavailable) # but need to substract ClusteredColl from each start (clst) 
      clst <- sample(size = 1, x = available) # select from possible non-overlapping starts
      notavailable <- c(notavailable, (clst-ClusteredColl):(clst+ClusteredColl-1))
      collectionyr_ix <- c(collectionyr_ix,c(clst:(clst+ClusteredColl-1)))
      i <- i+1
      if(length(available)<ClusteredColl){
        warning("Cannot find non-overlapping cluster period, set at frequency")
        collectionyr_ix <- c(collectionyr_ix, sample(size = (numCollections-length(collectionyr_ix)), setdiff((1:TotYrs),collectionyr_ix)))
        break
      }
    } # allow adding overlapping cluster periods
    Freqs <- rep(0,TotYrs)
    Freqs[collectionyr_ix] <- 1
    
    while(yrs <= TotYrs){
      # Transition between good and bad years
      Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
        "Good"
      } else {
        "Bad"
      } # Make next a good or bad year
      
      CollectYN <- Freqs[yrs] # this was in the Annual transition
      onetransition <- AnnualTransition()
      GBs[yrs] <- Yr
      if(CollectYN==1){
        Intens[yrs] <- Inten_m[Yr,]
      }
      
      PopSize[yrs] <- floor(sum(vec))
      lambdas[yrs] <- onetransition[[2]]
      mats[[yrs]] <- onetransition[[3]]
      yrs <- yrs+1
    }
  } else {
    ######################## Annual Transitions without clusters ######################
    while(yrs <= TotYrs){
      # Transition between good and bad years
      Yr <- if(rbinom(1,1,GoodBadTm[Yr,1])==1){
        "Good"
      } else {
        "Bad"
      } 
      CollectYN <- rbinom(1,1,Freq_m[Yr,]) # this was in the Annual transition
      onetransition <- AnnualTransition()  ### CHECK BACK HERE
      GBs[yrs] <- Yr
      Freqs[yrs] <- CollectYN
      if(CollectYN==1){
        Intens[yrs] <- Inten_m[Yr,]
      }
      
      vec <- onetransition[[1]] 
      PopSize[yrs] <- sum(vec)
      lambdas[yrs] <- onetransition[[2]]
      mats[[yrs]] <- onetransition[[3]]
      yrs <- yrs+1
    }
  } # ends all annual transitions without clusters
  SimDF <- data.frame(GB_yrs = GBs, Frequency = Freqs, Intensity = Intens, PopulationSize = PopSize, Year = 1:TotYrs, 
                      lambdas = lambdas)
  list(SimDF, FreqGB = Freq, IntGB = Int, GBMx = GoodBadTm, StartingPopSize = StartPopSize, Mx = mats,
       LogGrowthRate = stoch.growth.rate(mats, verbose=FALSE)) # $approx is log stochastic growth rate by Tuljapukar's approximation; $sim and CI
}


# ----------------------- Data -----------------------
# Ellis et al. 2012 paper
tm <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Kaye_transitionmatrices.csv")

sp_info <- read.csv("C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Kaye_speciesinfo.csv")

m1 <- gsub("\\[|\\]", "", tm$Mx[100]) 
m2 <- strsplit(m1, ";")  # remove the ";" that separates rows
(m3 <- lapply(m2, function(x) strsplit(trimws(x,"l"), " "))) # split each list element after removing leading whitespace
(m4 <- matrix(as.numeric(unlist(m3)), nrow = length(m3[[1]][[1]]), byrow = TRUE))
# Mx: Annual transition matrix
# Tmx: Annual trnasition matrix with transition probabilities, no fecundity
# Nx: Vector of observed stage structures (counts in each stage); must be at the start

# NA for unobserved stages
Nx_all <- lapply(1:nrow(tm), function(i){
  v1 <- gsub("\\[|\\]", "", tm$Nx[i])
  v2 <- sapply(v1, function(x) strsplit(trimws(x,"l"), " ")) 
  v3 <- as.numeric(unlist(v2[[1]]))
  v3
})

#CIPI get rid of the matrix with different number of stages, do that on its own, most have 3, one population has 6, rename species CIPI3 and CIPI6
Nx_names <- lapply(1:nrow(tm), function(i){
  spp <- tm$SPP[i]
  pop <- tm$POP[i]
  yr <- tm$YR[i]
  paste(spp, pop, yr)
})

# Nx_all[grep("CIPI", unlist(Nx_names))][which(sapply(Nx_all[grep("CIPI", unlist(Nx_names))], function(x) length(x))==3)]

# which(sapply(Nx_all[grep("CIPI", unlist(Nx_names))], function(x) length(x))==3)
# which(sapply(Nx_all[grep("CIPI", unlist(Nx_names))], function(x) length(x))==6)

Nx_names[grep("CIPI", Nx_names)][1:30] <- str_replace_all(Nx_names[grep("CIPI", Nx_names)][1:30], "CIPI", "CIPI3")
Nx_names[grep("CIPI", Nx_names)][31:35] <- str_replace_all(Nx_names[grep("CIPI", Nx_names)][31:35], "CIPI", "CIPI6")

Nx_all[grep("CIPI3", unlist(Nx_names))]
Nx_all[grep("CIPI6", unlist(Nx_names))]

tm$SPP <- as.character(tm$SPP)
tm$SPP[grep("CIPI", tm$SPP)][1:30] <- "CIPI3"
tm$SPP[grep("CIPI", tm$SPP)][31:35] <- "CIPI6"
tm$SPP <- as.factor(tm$SPP)

names(Nx_all) <- Nx_names

# Transitions without fecundity
TMx_all <- sapply(1:nrow(tm), function(i){
  m1 <- gsub("\\[|\\]", "", tm$Tmx[i])
  m2 <- strsplit(m1, ";")  # remove the ";" that separates rows
  m3 <- lapply(m2, function(x) strsplit(trimws(x,"l"), " ")) # split each list element after removing leading whitespace
  m4 <- matrix(as.numeric(unlist(m3)), nrow = length(m3[[1]][[1]]), byrow = TRUE)
  m4
}, simplify =FALSE, USE.NAMES = TRUE)

names(TMx_all) <- Nx_names

sp_info$SPP <- as.character(sp_info$SPP)
sp_info <- rbind(sp_info, sp_info[sp_info$SPP == "CIPI",])
sp_info$SPP[sp_info$SPP == "CIPI"][1] <- "CIPI3" 
sp_info$SPP[sp_info$SPP == "CIPI"] <- "CIPI6" 



# Transitions and fecundity
Mx_all <- sapply(1:nrow(tm), function(i){
  m1 <- gsub("\\[|\\]", "", tm$Mx[i])
  m2 <- strsplit(m1, ";")  # remove the ";" that separates rows
  m3 <- lapply(m2, function(x) strsplit(trimws(x,"l"), " ")) # split each list element after removing leading whitespace
  m4 <- matrix(as.numeric(unlist(m3)), nrow = length(m3[[1]][[1]]), byrow = TRUE)
  m4
}, simplify =FALSE, USE.NAMES = TRUE)

names(Mx_all) <- Nx_names
Mx_all[grep("CIPI", names(Mx_all))]

lambdas <- sapply(Mx_all, function(x) lambda(x), simplify = FALSE, USE.NAMES = TRUE)

fecundity <- mapply(function(x,y){ 
  fec <- x-y
  rowSums(fec)[1] #Only the first row is funcunity across classes (subtracted out survival by Mx-TMx)
}, Mx_all, TMx_all)

x <- unique(tm$SPP)[3]

# --------------------------- Annuals ---------------------------------------------
load("C:/Users/deprengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/PVA/Paniw etal/patterns_temporal_autocorrelation-master/COMPADRE_v.4.0.0.RData")
AEVI <- compadre$mat[compadre$metadata$Genus == "Aeschynomene"]
Nx_names_annuals <- list("AEVI AEVI_1", "AEVI AEVI_2", "AEVI AEVI_2")
Nx_all_annuals <- list(c(NA,NA,1,1,1), c(NA,NA,1,1,1), c(NA,NA,1,1,1)) # Don't know but only need which are NA - unobserved seedbank
Mx_all_annuals <- lapply(AEVI, function(x) x$matA) # they are all singular, can't solve, where does it invert? 
TMx_all_annuals <- lapply(AEVI, function(x) x$matU)

Nx_names_SEER <- list("SEER 2004", "SEER 2005", "SEER 2006", "SEER 2007") 
Nx_all_SEER <- list(c(NA,1),c(NA,1),c(NA,1),c(NA,1)) # I have the number per population and the average but don't need
Mx_all_SEER <- list(matrix(c(0.409,10.667,0.026,1.47), ncol = 2, byrow = TRUE),
                    matrix(c(0.409,2.532,0.006,0.25), ncol = 2, byrow = TRUE),
                    matrix(c(0.409,1.244,0.011,2.14), ncol = 2, byrow = TRUE),
                    matrix(c(0.409,4.95,0.001,0.03), ncol = 2, byrow = TRUE))
TMx_all_SEER <- list(matrix(c(0.409,0,0.026,0), ncol = 2, byrow = TRUE),
                     matrix(c(0.409,0,0.006,0), ncol = 2, byrow = TRUE),
                     matrix(c(0.409,0,0.011,0), ncol = 2, byrow = TRUE),
                     matrix(c(0.409,0,0.001,0), ncol = 2, byrow = TRUE))


# --------------- Lifespan --------------------
# Enter plants in first observed stage
seedlingbins <- data.frame(SPP = unique(tm$SPP), 
                           seedlingBin = do.call(rbind,lapply(Nx_all[!duplicated(tm$SPP)], function(x){
                             if(length(which(is.na(x[1])))==0){
                               1
                             } else {
                               which(is.na(x[1]))+1
                             }
                           })))

lifelengthsXmatrix <- do.call(rbind,lapply(1:length(Mx_all), function(i){
  data.frame(Lifespan = lifespanMPM(Mx_all[i], TMx_all[i], 
                                    seedlingbins$seedlingBin[seedlingbins$SPP == strsplit(Nx_names[[i]][1]," ")[[1]][1] ]),
             SPP = strsplit(Nx_names[[i]][1]," ")[[1]][1], Site = strsplit(Nx_names[[i]][1]," ")[[1]][2],
             Yr =  strsplit(Nx_names[[i]][1]," ")[[1]][3])
}))

generationtimeXspecies <- do.call(rbind,lapply(sp_info$SPP, function(sp){
  matF <- lapply(1:length(grep(sp, Nx_names)), function(x){ 
    Mx_all[grep(sp, Nx_names)][[x]] - TMx_all[grep(sp, Nx_names)][[x]]
  })
  data.frame(generationTime = generation.time(mean(Mx_all[grep(sp, Nx_names)]),
                                              r = mean(matF)), 
             SPP = sp)
}))

avg.lifespan.sp <- aggregate(Lifespan~SPP,data = lifelengthsXmatrix, median)
sd.lifespan.sp <- aggregate(Lifespan~SPP,data = lifelengthsXmatrix, sd)

lifespan.sp <- cbind(avg.lifespan.sp, SDlifespan = sd.lifespan.sp[,2])
lifespan.sp[order(lifespan.sp$Lifespan),]

lifelengthsXmatrix[lifelengthsXmatrix$SPP == "NEMA",]

lifespan.sp <- merge(lifespan.sp, generationtimeXspecies, by = "SPP")

# ---------------- Basic Simulation -----------------------

pathstart <- "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Simulation_output/"

# lifespan is really generation time
SimSimple <- function(HarvestType = "No", sppVector = unique(tm$SPP), FrG = 0, FrB = 0, intG = 0, intB = 0, 
                      AvgInt = 0, Cgg_gb_bg_bb = c(0.5,0.5,0.5,0.5), reps = 100, Mxs, TMxs, Nxs,
                      Simlength = 50*(floor(lifespan.sp$generationTime[lifespan.sp$SPP==sp])), 
                      lifespan = floor(lifespan.sp$generationTime[lifespan.sp$SPP==sp])){
  lapply(sppVector, function(sp){
    bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
      outdf <- do.call(rbind,lapply(1:reps, function(repNow){
        outSimulation <- SeedHarvestSim(Mx_list = Mxs,
                                        TMx_list = TMxs,
                                        Nx_list = Nxs,
                                        Nx = stable.stage(mean(Mxs)),
                                        StartPopSize = popsz,
                                        GoodBadTm = matrix(Cgg_gb_bg_bb, 
                                                           nrow=2, byrow = TRUE, 
                                                           dimnames = list(c("Good","Bad"),
                                                                           c("Good","Bad"))),
                                        Freq = c(FrG,FrB), Int = c(intG,intB), 
                                        TotYrs = Simlength,
                                        ClusteredColl = 1, # clust,
                                        ddceiling = FALSE)
        # Damping ratio
        AMx <- mean(outSimulation[[6]])
        dampR <- eigen(AMx)$values[1]/abs(eigen(AMx)$values[2])
        
        outSimulation[[1]]$Rep <- repNow
        outSimulation[[1]]$SPP <- sp
        if(any(outSimulation[[1]]$PopulationSize==0)){
          Exterpated <- 1
          Time2Ext <- which(outSimulation[[1]]$PopulationSize == 0)[1]
        } else {
          Exterpated <- 0
          Time2Ext <- NA
        }
        
        # Extinction probability, stochastic lambda, and damping ratio
        out <- data.frame(Tulapprox = outSimulation$LogGrowthRate$approx,
                          LogGrowthSim = outSimulation$LogGrowthRate$sim,
                          Exterpated, Time2Ext, 
                          DampingRatio = as.numeric(dampR),
                          FreqG = FrG, FreqB = FrB, 
                          IntG = intG, IntB = intB,
                          SPP = sp, Clust = 1, StPopSz = popsz, 
                          lifespan = lifespan,
                          Frequency = 0, IntRatio = 0)# Int[1]/Int[2])
        rm(outSimulation)
        out
      })) # end replicates
      outdf
    })) # end bypopsz 
    # bypopsz
    save(bypopsz, file = paste(pathstart, HarvestType,"Harvest",
                               "Fr",FrG,"_",FrB,"In",intG,"_",intB,"AverageInt",AvgInt,
                               "Climate",paste(Cgg_gb_bg_bb, collapse = "_"),sp,
                               ".Rdata", sep=""))
  }) # end species
}


# ---------------- No Harvest different climate ------------
SimSimple(HarvestType = "No", sppVector = unique(tm$SPP), Cgg_gb_bg_bb = c(0.25,0.75,0.25,0.75))


# ----------------- Harvest across intensity and frequency different climate -------------
lapply(seq(0.1,1, by = 0.1), function(intG){
  lapply(seq(0.1,1, by = 0.1), function(FrG){
    SimSimple(HarvestType = "GoodYrOnly", sppVector = unique(tm$SPP), intG = intG, FrG = FrG, Cgg_gb_bg_bb = c(0.25,0.75,0.25,0.75) )
  })
})


# ---------------- Harvest same intensity and Frequency in both good and bad yeras ----------------


# --------------- Only good across levels --------------------------
lapply(c(0.1,0.25,0.5,0.75,0.9), function(int){
  lapply(c(0.1,0.25,0.5,0.75), function(fr){
    SimSimple(HarvestType = "OnlyGB", sppVector = unique(tm$SPP), 
              intG = int, intB = 0, FrG = fr, FrB = 0)
  })
})

# -------------------- Plots on a RStudio not running everything else --------------
# ----------------------- 2020.11.15 --------------------------------
# ----------------------- desktop only ------------------------------
# load to bind and plot
HarvestType <- "SameGB"
files <- list.files(path = pathstart, pattern = paste(HarvestType, "Harvest", sep=""))
toplot <- do.call(rbind, lapply(files, function(x){
  load(file = paste(pathstart,x,sep = ""))
  bypopsz
}))

library(scales)
library(splines)
library(MASS)
library(nlme)

table(toplot$SPP, toplot$IntB)

model <- glm(DampingRatio ~ lifespan + FreqG + FreqB, data = toplot, family = poisson)
summary(model)

# DRpNo <-
  ggplot(toplot, aes(lifespan, DampingRatio, colour= as.factor(FreqG)))+ 
  geom_jitter(height = 0, width = 1, alpha = 0.4, shape = 20)+
  # stat_smooth(method = "glm", family = "quasipoisson", formula = y ~ ns(x, 3))+
  # geom_smooth(method = MASS::glm.nb, aes(color = "NB"), se = FALSE)+
  geom_line(aes(y=fitted(model)))+
    facet_wrap(~as.factor(IntG))+
  theme_bw()+
  ylab("Damping Ratio")+
  xlab("Generation Time")


ggsave(filename = "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/DampingRatioNoharvest.jpg", DRpNo,
       width=100, height=100,units='mm', dpi=300)

StochLamNo <- ggplot(toplot, aes(lifespan, Tulapprox))+
  geom_jitter(height = 0, width = 0.5, alpha = 0.4, shape = 20)+
  stat_smooth(method = "glm")+
  theme_bw()+
  ylab("Log stochastic growth rate \n (Tuljapukar's approximation)")+
  xlab("Generation Time")

ggsave(filename = "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/StochLambdaNoharvest.jpg", StochLamNo,
       width=100, height=100,units='mm', dpi=300)

ggsave(filename = "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/Time2ExtNoharvest.jpg",
       ggplot(toplot, aes(as.factor(StPopSz), Time2Ext))+
         geom_boxplot()+
         theme_bw()+
         facet_wrap(~SPP, scales = "free")+
         geom_jitter(height=0, width = 0.01)+
         ylab("Time to extinction (< 1)")+
         xlab("Starting population size"), width=300, height=300,units='mm', dpi=300)


extinctionrisk <- aggregate(Exterpated ~ SPP + StPopSz + FreqG + FreqB + IntG + IntB + lifespan, function(x) sum(x)/length(x) ,data = toplot)

ggsave(filename = "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/ExtRiskNoharvest.jpg",
       ggplot(extinctionrisk, aes(lifespan, Exterpated))+
         # geom_smooth(method=MASS::glm.nb, se = FALSE)+
         stat_smooth(method="glm", family="binomial", formula=y~log(x),
                     alpha=0.2, size=0.2, se=FALSE)+
         theme_bw()+
         facet_wrap(~as.factor(StPopSz))+
         geom_point()+
         xlab("Generation Time")+
         ylab("Extinction Risk"),
       width=100, height=100,units='mm', dpi=300)



# ---------------- Basic Simulation with 100 years instead of generation time -----------------------

pathstart100 <- "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Simulation_output_100years/"

script <- getURL("https://raw.githubusercontent.com/DenverBotanicGardens/Seed_Harvest_Modeling/master/R/SimulationMPM.R",ssl.verifypeer = FALSE)
# See SimulationMPM.R
# lifespan is really generation time
# SimSimple100 <- function(HarvestType = "No", sppVector = unique(tm$SPP), FrG = 0, FrB = 0, intG = 0, intB = 0, 
#                       AvgInt = 0, Cgg_gb_bg_bb = c(0.5,0.5,0.5,0.5), reps = 100, Mxs, TMxs, Nxs,
#                       Simlength = 100){  
#   lapply(sppVector, function(sp){
#     bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
#       outdf <- do.call(rbind,lapply(1:reps, function(repNow){
#         outSimulation <- SeedHarvestSim(Mx_list = Mx_all[grep(sp, unlist(Nx_names))],
#                                         TMx_list = TMx_all[grep(sp, unlist(Nx_names))],
#                                         Nx_list = Nx_all[grep(sp, unlist(Nx_names))],
#                                         Nx = stable.stage(mean(Mx_all[grep(sp, unlist(Nx_names))])),
#                                         StartPopSize = popsz,
#                                         GoodBadTm = matrix(Cgg_gb_bg_bb, 
#                                                            nrow=2, byrow = TRUE, 
#                                                            dimnames = list(c("Good","Bad"),
#                                                                            c("Good","Bad"))),
#                                         Freq = c(FrG,FrB), Int = c(intG,intB), 
#                                         TotYrs = Simlength,
#                                         ClusteredColl = 1, # clust,
#                                         ddceiling = FALSE)
#         # Damping ratio
#         AMx <- mean(outSimulation[[6]])
#         dampR <- eigen(AMx)$values[1]/abs(eigen(AMx)$values[2])
#         
#         outSimulation[[1]]$Rep <- repNow
#         outSimulation[[1]]$SPP <- sp
#         if(any(outSimulation[[1]]$PopulationSize==0)){
#           Exterpated <- 1
#           Time2Ext <- which(outSimulation[[1]]$PopulationSize == 0)[1]
#         } else {
#           Exterpated <- 0
#           Time2Ext <- NA
#         }
#         
#         # Extinction probability, stochastic lambda, and damping ratio
#         out <- data.frame(Tulapprox = outSimulation$LogGrowthRate$approx,
#                           LogGrowthSim = outSimulation$LogGrowthRate$sim,
#                           Exterpated, Time2Ext, 
#                           DampingRatio = as.numeric(dampR),
#                           FreqG = FrG, FreqB = FrB, 
#                           IntG = intG, IntB = intB,
#                           SPP = sp, Clust = 1, StPopSz = popsz, 
#                           lifespan = floor(lifespan.sp$generationTime[lifespan.sp$SPP==sp]),
#                           Frequency = 0, IntRatio = 0)# Int[1]/Int[2])
#         rm(outSimulation)
#         out
#       })) # end replicates
#       outdf
#     })) # end bypopsz 
#     # bypopsz
#     save(bypopsz, file = paste(pathstart100, HarvestType,"Harvest",
#                                "Fr",FrG,"_",FrB,"In",intG,"_",intB,"AverageInt",AvgInt,
#                                "Climate",paste(Cgg_gb_bg_bb, collapse = "_"),sp,
#                                ".Rdata", sep=""))
#   }) # end species
# }
# 



# ---------------- No Harvest for only 100 years  ------------
SimSimple100(HarvestType = "No100years", sppVector = unique(tm$SPP)[-1],  Mxs = Mx_all, TMxs = TMx_all, Nxs = Nx_all)

SimSimple100(HarvestType = "100yearsRandomNx", sppVector = unique(tm$SPP)[-1],  Mxs = Mx_all, TMxs = TMx_all, Nxs = Nx_all, stablestagestart = FALSE)


SimSimple100(HarvestType = "No10years", sppVector = unique(tm$SPP)[-1],  Mxs = Mx_all, TMxs = TMx_all, Nxs = Nx_all, Simlength = 10)

SimSimple100(HarvestType = "10yearsRandomNx", sppVector = unique(tm$SPP)[-1],  Mxs = Mx_all, TMxs = TMx_all, Nxs = Nx_all, Simlength = 10, stablestagestart = FALSE)


# --------------- Need to adjust simple100 to add new species without the lifespan.sp table ----------------------
# lifespan is really generation time



SimSimple100new <- function(HarvestType = "No", sppVector = sp, FrG = 0, FrB = 0, intG = 0, intB = 0, 
                         AvgInt = 0, Cgg_gb_bg_bb = c(0.5,0.5,0.5,0.5), reps = 100, Mxs, TMxs, Nxs,
                         lifespan = 1){  
  matF <- lapply(1:length(Mxs), function(i){
    Mxs[[i]] - TMxs[[i]]
  })
  generationtime <- popbio::generation.time(mean(Mxs),r = mean(matF))
  
  lapply(sppVector, function(sp){
    bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
      outdf <- do.call(rbind,lapply(1:reps, function(repNow){
        outSimulation <- SeedHarvestSim(Mx_list = Mxs,
                                        TMx_list = TMxs,
                                        Nx_list = Nxs,
                                        Nx = stable.stage(mean(Mxs)),
                                        StartPopSize = popsz,
                                        GoodBadTm = matrix(Cgg_gb_bg_bb, 
                                                           nrow=2, byrow = TRUE, 
                                                           dimnames = list(c("Good","Bad"),
                                                                           c("Good","Bad"))),
                                        Freq = c(FrG,FrB), Int = c(intG,intB), 
                                        TotYrs = 100,
                                        ClusteredColl = 1, 
                                        ddceiling = FALSE)
        # Damping ratio
        # need to make new added by defining matrices into square matrices
        AMx <- mean(lapply(outSimulation[[6]], function(x) matrix(x, ncol = sqrt(nrow(outSimulation[[6]][1][[1]])), byrow = FALSE)))
        dampR <- eigen(AMx)$values[1]/abs(eigen(AMx)$values[2])
        
        outSimulation[[1]]$Rep <- repNow
        outSimulation[[1]]$SPP <- sp
        if(any(outSimulation[[1]]$PopulationSize==0)){
          Exterpated <- 1
          Time2Ext <- which(outSimulation[[1]]$PopulationSize == 0)[1]
        } else {
          Exterpated <- 0
          Time2Ext <- NA
        }
        
        # Extinction probability, stochastic lambda, and damping ratio
        out <- data.frame(Tulapprox = outSimulation$LogGrowthRate$approx,
                          LogGrowthSim = outSimulation$LogGrowthRate$sim,
                          Exterpated, Time2Ext, 
                          DampingRatio = as.numeric(dampR),
                          FreqG = FrG, FreqB = FrB, 
                          IntG = intG, IntB = intB,
                          SPP = sp, Clust = 1, StPopSz = popsz, 
                          lifespan = generationtime,
                          Frequency = 0, IntRatio = 0)# Int[1]/Int[2])
        rm(outSimulation)
        out
      })) # end replicates
      outdf
    })) # end bypopsz 
    # bypopsz
    save(bypopsz, file = paste(pathstart100, HarvestType,"Harvest",
                               "Fr",FrG,"_",FrB,"In",intG,"_",intB,"AverageInt",AvgInt,
                               "Climate",paste(Cgg_gb_bg_bb, collapse = "_"),sp,
                               ".Rdata", sep=""))
  }) # end species
}

# --------------------------- Add in no harvest for only 100 years for AEVI -----------------
SimSimple100new(HarvestType = "No100years", sppVector = "AEVI",  
             Mxs = Mx_all_annuals, TMxs = TMx_all_annuals, Nxs = Nx_all_annuals)

# --------------- Same for good and bad across levels --------------
lapply(c(0.1,0.25,0.5,0.75,0.9), function(int){
  lapply(c(0.1,0.25,0.5,0.75,0.9), function(fr){
    lapply(unique(tm$SPP), function(sp){
      SimSimple100(HarvestType = "SameGB100", sppVector = sp, 
                   Mxs = Mx_all[grep(sp, unlist(Nx_names))],
                   TMxs = TMx_all[grep(sp, unlist(Nx_names))],
                   Nxs = Nx_all[grep(sp, unlist(Nx_names))],
                   intG = int, intB = int, FrG = fr, FrB = fr, 
                   lifespan = lifespan.sp$generationTime[lifespan.sp$SPP==sp])
      }) # by species
    }) # by frequencies
  }) # by intensities




# ---------------------------- go for the whole 50 generations but for additional species --------------------
pathstart <- "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Simulation_output/"

lapply(c(0.1,0.5,0.9), function(Fr){
  lapply(c(0.1,0.25,0.5,0.75,0.9), function(AvgInt){
    x1 <-combn(seq(0,1,by=0.1),2)
    cols <- which(apply(x1,2,function(x) sum(x)/2) == AvgInt)
    ints <- x1[,cols, drop = FALSE] # keep matrix structure even when only one column
    lapply(1:(2*ncol(ints)),function(i){
      if(i <= ncol(ints)){
        intG <- ints[1,i]
        intB <- ints[2,i]
      }  else {
        intG <- ints[2,i-ncol(ints)]
        intB <- ints[1,i-ncol(ints)]
      }
      lapply(1:100, function(repeatAvgInt){
        SimSimple(HarvestType = "BiasAvg", AvgInt = AvgInt, sppVector = unique(tm$SPP), intG = intG, intB = intB,
                  FrG = Fr, FrB = Fr)
      }) # repeat each combination 100 times
    }) # work through all combinations where you can get that average wtih both good and bad having the larger amount
  }) # end average intentsity of 0.1, 0.25, 0.5, 0.75, 0.9
}) # end frequency 0.1, half, 0.9


# Missed 0.9 freq last time
# --------------- Same for good and bad across levels --------------
lapply(c(0.1,0.25,0.5,0.75,0.9), function(int){
  lapply(c(0.9), function(fr){
    SimSimple(HarvestType = "SameGB", sppVector = unique(tm$SPP), 
              intG = int, intB = int, FrG = fr, FrB = fr)
  })
})
