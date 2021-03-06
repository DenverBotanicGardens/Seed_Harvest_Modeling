# testing
# Nx_list <- Nx_all[grep(sp, unlist(Nx_names))]
# Nx <- sample(Nx_all[grep(sp, unlist(Nx_names))],1)[[1]]
# species <- unique(tm$SPP)[1]
# Mxs <- Mx_all[grep(species, unlist(Nx_names))]
# TMxs <- TMx_all[grep(species, unlist(Nx_names))]
# Fm <- (mean(Mxs)-mean(TMxs))
# which(Fm > 0, arr.ind = TRUE)
# rm(Mxs);rm(TMxs);rm(species);rm(Fm);rm(Nx_list)

# stop doing all species at once, just one at a time, sppVector should just be character of what the species is
SimSimple <- function(HarvestType = "No", sppVector, FrG = 0, FrB = 0, intG = 0, intB = 0, 
                         AvgInt = 0, Cgg_gb_bg_bb = c(0.5,0.5,0.5,0.5), reps = 100, Mxs, TMxs, Nxs,
                         Simlength = 100, generationspan = FALSE, stablestagestart = TRUE, ps = pathstart){  
  lapply(sppVector, function(sp){
    gentim <- popbio::generation.time(mean(Mxs))
    if(is.infinite(gentim)){
      Fm <- (mean(Mxs)-mean(TMxs))
      gentim <-  popbio::generation.time(mean(Mxs), 
                                         r = c(unique(which(Fm > 0, arr.ind = TRUE)[,1])), 
                                         c = c(unique(which(Fm > 0, arr.ind = TRUE)[,2])))}
    if(stablestagestart == TRUE) Nx_start <- stable.stage(mean(Mxs))
    if(stablestagestart == FALSE) Nx_start <- sample(Nxs,1)[[1]]
    if(generationspan == TRUE) {
      projlength <- floor(Simlength * gentim)
    } else {
      projlength <- Simlength
    }
    bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
      outdf <- do.call(rbind,lapply(1:reps, function(repNow){
        if(stablestagestart == TRUE) Nx_start <- stable.stage(mean(Mxs))
        if(stablestagestart == FALSE) Nx_start <- sample(Nx_all[grep(sp, unlist(Nx_names))],1)[[1]]
        outSimulation <- SeedHarvestSim(Mx_list = Mxs,
                                        TMx_list = TMxs,
                                        Nx_list = Nxs,
                                        Nx = Nx_start,
                                        StartPopSize = popsz,
                                        GoodBadTm = matrix(Cgg_gb_bg_bb, 
                                                           nrow=2, byrow = TRUE, 
                                                           dimnames = list(c("Good","Bad"),
                                                                           c("Good","Bad"))),
                                        Freq = c(FrG,FrB), Int = c(intG,intB), 
                                        TotYrs = projlength,
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
                          generationtime = gentim,
                          Frequency = 0, IntRatio = 0)# Int[1]/Int[2])
        rm(outSimulation)
        out
      })) # end replicates
      outdf
    })) # end bypopsz 
    # bypopsz
    save(bypopsz, file = paste(ps, HarvestType,"Harvest",
                               "Fr",FrG,"_",FrB,"In",intG,"_",intB,"AverageInt",AvgInt,
                               "Climate",paste(Cgg_gb_bg_bb, collapse = "_"),sp,
                               ".Rdata", sep=""))
  }) # end species
}


# ---------------- Some differences with an additional species not part of Ellis et al. 2012 ---------------

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


# -------------------------------------- Virtual species --------------------------------
# TMx excludes fecundity
# Mx survival, growth, and fecundity
# Nx starting poulation stage distriubtions so should be c(NA,NA,NA,x) for annuals 
# Testing
# Mxs <- Tmx_annual
# hist(unlist(lapply(1:length(Mxs), function(x) Mxs[[x]][1,4])))
# Nxs <- lapply(1:100, function(i){
#   c(max(1,rpois(1,runif(1, 1,500))),NA,NA,max(1,rpois(1,runif(1, 1,500))))
#   })

# HarvestType <- "NoneVirtual"; MatrixType <- "Annual"; Mxs <- Tmx_annual
# Nxs <- Nxs_annual; Simlength <- 100; generationspan <- FALSE 
# stablestagestart <- TRUE; ps <- pathstartVirtual
# popsz <- 50
# rm(HarvestType);rm(MatrixType);rm(Mxs);rm(Nxs);rm(Simlength);rm(generationspan);rm(stablestagestart);rm(ps);rm(popsz)

SimSimpleVirtual <- function(HarvestType = "No", MatrixType, FrG = 0, FrB = 0, intG = 0, intB = 0, 
                             AvgInt = 0, Cgg_gb_bg_bb = c(0.5,0.5,0.5,0.5), reps = 100, Mxs, Nxs,
                             Simlength = 100, generationspan = FALSE, stablestagestart = TRUE, ps = pathstart){  
  lapply(1:100, function(sp){ # repeat for random number of matrices from each type
    Mx_sample <- sample(Mxs, size = max(5,rpois(1,10)))
    gentim <- popbio::generation.time(mean(Mx_sample))
    bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
      outdf <- do.call(rbind,lapply(1:reps, function(repNow){
        if(stablestagestart == TRUE) Nx_start <- stable.stage(mean(Mx_sample))
        if(stablestagestart == FALSE) Nx_start <- sample(Nxs,1)[[1]]
        if(generationspan == TRUE) {
          projlength <- floor(Simlength * gentim)
        } else {
          projlength <- Simlength
        }
        outSimulation <- SeedHarvestSim(Mx_list = Mx_sample,
                                        # The TMx are transition matrices without fecundity, virtual species only have largest as fecundity
                                        TMx_list = lapply(Mx_sample, function(i){ 
                                          i[1,4] <- 0
                                          i
                                        }),
                                        Nx_list = Nxs,
                                        Nx = Nx_start,
                                        StartPopSize = popsz,
                                        GoodBadTm = matrix(Cgg_gb_bg_bb, 
                                                           nrow=2, byrow = TRUE, 
                                                           dimnames = list(c("Good","Bad"),
                                                                           c("Good","Bad"))),
                                        Freq = c(FrG,FrB), Int = c(intG,intB), 
                                        TotYrs = projlength,
                                        ClusteredColl = 1, # clust,
                                        ddceiling = FALSE)
        # save(outSimulation, file = paste(ps, HarvestType,"Harvest",
        #                            "Fr",FrG,"_",FrB,"In",intG,"_",intB,"AverageInt",AvgInt,
        #                            "Climate",paste(Cgg_gb_bg_bb, collapse = "_"),"Virtual",MatrixType, sp,
        #                            ".Rdata", sep=""))
        
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
                          SPP = MatrixType, Clust = 1, StPopSz = popsz, 
                          generationtime = gentim,
                          SimLength = projlength,
                          MatrixType = MatrixType,
                          Frequency = 0, IntRatio = 0)# Int[1]/Int[2])
        rm(outSimulation)
        out
      })) # end replicates
      outdf
    })) # end bypopsz 
    # bypopsz
    save(bypopsz, file = paste(ps, HarvestType,"Harvest",
                               "Fr",FrG,"_",FrB,"In",intG,"_",intB,"AverageInt",AvgInt,
                               "Climate",paste(Cgg_gb_bg_bb, collapse = "_"),MatrixType,sp,
                               ".Rdata", sep=""))
  }) # end species
}


# Test
ps <- pathstart
MatrixType <- "Iteroslow"
Mxs <- Tmx_iteroslow
Nxs <- Nxs_othertypes
Simlength <- 100
generationspan <- FALSE
popsz <- 10
reps <- 100
stablestagestart <- TRUE
rm(ps);rm(MatrixType);rm(Mxs);rm(Nxs);rm(Simlength);rm(generationspan);rm(popsz)

SimVirtualNoHarvest <- function(MatrixType,reps = 100, Mxs, Nxs,
                                Simlength = 100, generationspan = FALSE, stablestagestart = TRUE, ps = pathstart){  
  lapply(1:100, function(sp){ # repeat for random number of matrices from each type
    Mx_sample <- sample(Mxs, size = max(5,rpois(1,10)))
    gentim <- popbio::generation.time(mean(Mx_sample))
    bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
      outdf <- do.call(rbind,lapply(1:reps, function(repNow){
        if(stablestagestart == TRUE) Nx_start <- stable.stage(mean(Mx_sample))
        if(stablestagestart == FALSE) Nx_start <- sample(Nxs,1)[[1]]
        if(generationspan == TRUE) {
          projlength <- floor(Simlength * gentim)
        } else {
          projlength <- Simlength
        }
        outSimulation <- SeedHarvestSimple(Mx_list = Mx_sample,
                                        # The TMx are transition matrices without fecundity, virtual species only have largest as fecundity
                                        TMx_list = lapply(Mx_sample, function(i){ 
                                          i[1,4] <- 0
                                          i
                                        }),
                                        Nx_list = Nxs,
                                        Nx = Nx_start,
                                        StartPopSize = popsz,
                                        TotYrs = projlength)
        
        out <- data.frame(outSimulation, SPP = MatrixType, StPopSz = popsz, 
                          generationtime = gentim,
                          SimLength = projlength,
                          MatrixType = MatrixType)
        rm(outSimulation)
        out
      })) # end replicates
      outdf
    })) # end bypopsz 
    # bypopsz
    save(bypopsz, file = paste(ps,"PVA","SimulationLength",projlength,MatrixType,sp,
                               ".Rdata", sep=""))
  }) # end species
}
