library(popbio)

# -----------------------Virtual -----------------------------

Tmx_iteroslow <- lapply(1:100, function(Mx){
  s <- runif(1,0,0.8) # adult survival
  f <- -1.25 * ((s + 4.6)^2) + 36.45 + rnorm(1, mean = 0, sd = 0.1) # add some variance around the line
  t_ij <- matrix(0, nrow = 4, ncol = 4)
  t_ij[2,1] <- runif(1, min = 0.01, max = 0.99)
  t_ij[3,2] <- runif(1, min = 0.01, max = 0.99)
  t_ij[4,3] <- runif(1, min = 0.01, max = 0.99)
  # stasis or retrogression
  t_ij[2,2:3] <- t_ij[3,3:4] <- t_ij[4,4] <- s
  if(any(colSums(t_ij)>=1)){
    for(i in which(colSums(t_ij)>=1)){
      t_ij[,i] <- t_ij[,i]/(sum(t_ij[,i]))
    } 
  }
  t_ij[1,4] <- f
  t_ij
})


Tmx_semelslow <- lapply(1:100, function(Mx){
  s <- runif(1,0,0.8) # adult survival
  f <- -1.25 * ((s + 4.6)^2) + 36.45 + rnorm(1, mean = 0, sd = 0.1) # add some variance around the line
  t_ij <- matrix(0, nrow = 4, ncol = 4)
  t_ij[2,1] <- runif(1, min = 0.01, max = 0.99)
  t_ij[3,2] <- runif(1, min = 0.01, max = 0.99)
  t_ij[4,3] <- runif(1, min = 0.01, max = 0.99)
  # stasis or retrogression
  t_ij[2,2:3] <- t_ij[3,3:4] <- s
  if(any(colSums(t_ij)>=1)){
    for(i in which(colSums(t_ij)>=1)){
      t_ij[,i] <- t_ij[,i]/(sum(t_ij[,i]))
    } 
  }
  t_ij[1,4] <- f
  t_ij
})

Tmx_iterofast <- lapply(1:100, function(Mx){
  s <- runif(1,0,0.8) # adult survival
  f <- -1.25 * ((s + 4.6)^2) + 36.45 + rnorm(1, mean = 0, sd = 0.1) # add some variance around the line
  t_ij <- matrix(0, nrow = 4, ncol = 4)
  t_ij[2:4,1:3] <- runif(9, min = 0.01, max = 0.99)
  t_ij[2,3] <- 0
  t_ij[2,2] <- t_ij[3,3] <-t_ij[4,4] <- s
  
  if(any(colSums(t_ij)>=1)){
    for(i in which(colSums(t_ij)>=1)){
      t_ij[,i] <- t_ij[,i]/(sum(t_ij[,i]))
    } 
  }
  t_ij[1,4] <- f
  t_ij
})

Tmx_semelfast <- lapply(1:100, function(Mx){
  s <- runif(1,0,0.8) # adult survival
  f <- -1.25 * ((s + 4.6)^2) + 36.45 + rnorm(1, mean = 0, sd = 0.1) # add some variance around the line
  t_ij <- matrix(0, nrow = 4, ncol = 4)
  t_ij[2:4,1:3] <- runif(9, min = 0.01, max = 0.99)
  t_ij[2,3] <- 0
  t_ij[2,2] <- t_ij[3,3] <- s
  
  if(any(colSums(t_ij)>=1)){
    for(i in which(colSums(t_ij)>=1)){
      t_ij[,i] <- t_ij[,i]/(sum(t_ij[,i]))
    } 
  }
  t_ij[1,4] <- f
  t_ij
})


Tmx_annual <- lapply(1:100, function(rep){
  # three seed stages, 1 reproductive
  # lower adult survival then reproduction is higher
  s <- runif(1,0,0.8) # adult survival
  f <- -1.25 * ((s + 4.6)^2) + 36.45 + rnorm(1, mean = 0, sd = 0.1) # add some variance around the line
  t_ij <- matrix(0, nrow = 4, ncol = 4)
  t_ij[4,1:3] <- rlnorm(1, meanlog = s, sdlog = 0.5)
  t_ij[2,1] <- t_ij[3,2] <- runif(1, min = 0, max = s)
  if(any(colSums(t_ij)>=1)){
    for(i in which(colSums(t_ij)>=1)){
      t_ij[,i] <- t_ij[,i]/ (sum(t_ij[,i]) + runif(1,0.01,0.1)) # make there be some death, 80-99% survival
    } 
  }
  t_ij[1,4] <- f
  t_ij
})

# 2020-12-23
pathstartVirtual <-  "C:/Users/DePrengm/Denver Botanic Gardens/Conservation - Demography/SeedHarvestModelling/Seed_Harvest_Modeling/Simulation_output_virtual50Gen/"

Nxs_annual <- lapply(1:100, function(i){
  c(max(1,rpois(1,runif(1, 1,500))),NA,NA,max(1,rpois(1,runif(1, 1,500))))
  })


SimSimpleVirtual(HarvestType = "NoneVirtual50Gen", MatrixType = "Annual", Mxs = Tmx_annual,
                 Nxs = Nxs_annual, Simlength = 50, generationspan = TRUE,
                 stablestagestart = TRUE, ps = pathstartVirtual)

Nxs_othertypes <- lapply(1:100, function(i){
  rpois(4,runif(4, 1,500))
})

SimSimpleVirtual(HarvestType = "NoneVirtual50Gen", MatrixType = "SemelFast", Mxs = Tmx_semelfast,
                 Nxs = Nxs_othertypes, Simlength = 50, generationspan = TRUE,
                 stablestagestart = TRUE, ps = pathstartVirtual)


SimSimpleVirtual(HarvestType = "NoneVirtual50Gen", MatrixType = "SemelSlow", Mxs = Tmx_semelslow,
                 Nxs = Nxs_othertypes, Simlength = 50, generationspan = TRUE,
                 stablestagestart = TRUE, ps = pathstartVirtual)

SimSimpleVirtual(HarvestType = "NoneVirtual50Gen", MatrixType = "IteroFast", Mxs = Tmx_iterofast,
                 Nxs = Nxs_othertypes, Simlength = 50, generationspan = TRUE,
                 stablestagestart = TRUE, ps = pathstartVirtual)

SimSimpleVirtual(HarvestType = "NoneVirtual50Gen", MatrixType = "IteroSlow", Mxs = Tmx_iteroslow,
                 Nxs = Nxs_othertypes, Simlength = 50, generationspan = TRUE,
                 stablestagestart = TRUE, ps = pathstartVirtual)


