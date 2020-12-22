SimSimple100 <- function(HarvestType = "No", sppVector = unique(tm$SPP), FrG = 0, FrB = 0, intG = 0, intB = 0, 
                         AvgInt = 0, Cgg_gb_bg_bb = c(0.5,0.5,0.5,0.5), reps = 100, Mxs, TMxs, Nxs,
                         Simlength = 100){  
  lapply(sppVector, function(sp){
    bypopsz <- do.call(rbind,lapply(c(10,50,100,500), function(popsz){
      outdf <- do.call(rbind,lapply(1:reps, function(repNow){
        outSimulation <- SeedHarvestSim(Mx_list = Mx_all[grep(sp, unlist(Nx_names))],
                                        TMx_list = TMx_all[grep(sp, unlist(Nx_names))],
                                        Nx_list = Nx_all[grep(sp, unlist(Nx_names))],
                                        Nx = stable.stage(mean(Mx_all[grep(sp, unlist(Nx_names))])),
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
                          lifespan = floor(lifespan.sp$generationTime[lifespan.sp$SPP==sp]),
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

