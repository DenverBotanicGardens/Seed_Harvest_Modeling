library(popbio)
library(ggplot2)
library(devtools)
# library(Rcompadre)
library(popdemo)
library(Rage)
library(ggtern)

rm(list=ls())

# Rare plants from Salguero-Gomez et al. 2016
# Keyfitz' entropy < 1, Type I, K-selected species, low juvenile mortality with most individuals living to an old age. 

#Trade-off between survival and fecundity from Takada and Kawai (2020) 

# ------------------------------------------------ constrain lambda ---------------------------
# functions 
# Takada and Kawai 2020; s would be s_44
# Convex for semelparous - see Farrell 2020 and Takada and Kawai 2020
semel_fecundsurv <- function(s){
  -12.5*((s + 0.1)^2) + 10.125
}

# Concave for iteroparous from Takada and Kawai 2020
itero_fecundsurv <- function(s){
  12.5*((s - 0.9)^2) - 0.125
}

#Hazard
# h(x) <- alpha3*exp((-beta3 * x))

# Fujiwara and Diaz-Lopez 2017 
# The x's represent the median age of the stage class
# Salguero-Gomez has signs flipped, type I is the area where rare plants fall; equation 2 from Fujiwara and Diaz-Lopez 2017; hazard is
# h(x) = alpha2*exp(beta2*x); exponentially increasing risk of mortality with age, risk due to aging
survivalTypeI <- function(alpha2, beta2, x){
  exp((alpha2/beta2)*(1-(exp(beta2*x))))
}

plot(1:4, survivalTypeI(0.03,0.60,1:4), xlab = "Stage", ylab = "Survival (S)", type = "l") # expression("Stage"[1-4])

transitions <- function(b1, b2, x){
  b1*exp((x*-b2))
}


#Flat when b1 = 1 and b2 = 0
plot(1:4, transitions(.9,0.1,1:4)/sum(transitions(.9,0.1,1:4)))
plot(1:2, transitions(.9,0.1,1:2)/sum(transitions(.9,0.1,1:2)))

# decreasing 
plot(1:4, transitions(0.1,.9,1:4)/sum(transitions(0.1,.9,1:4)))
plot(1:2, transitions(0.1,.9,1:2)/sum(transitions(0.1,.9,1:2)))


params <- do.call(rbind, lapply(seq(0,.1, by=0.01), function(a2){
  outb2 <- do.call(rbind, lapply(seq(0.5,.7, by=0.01), function(b2){
    surv <- survivalTypeI(a2, b2, 1:4)
    data.frame(a2, b2, stage = c("x1","x2","x3","x4"), survival = surv)
  }))
}))


(paramstransitions <- params[params$stage == "x4" & params$survival > 0.6 & params$survival < 0.8,])

paramlist_type1 <- params[params$stage == "x4" & params$survival > 0.6 & params$survival < 0.8,]
paramsitero <-  params[interaction(params[,1:2]) %in% interaction(paramlist_type1[,1:2]),]
(slowmin <- aggregate(survival ~ stage, min, data = params[interaction(params[,1:2]) %in% interaction(paramlist_type1[,1:2]),]))
(slowmax <- aggregate(survival ~ stage, max, data = params[interaction(params[,1:2]) %in% interaction(paramlist_type1[,1:2]),]))
min(paramsitero$a2)
max(paramsitero$a2)

min(paramsitero$b2)
max(paramsitero$b2)
# ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveTypeIIteroParams.jpg",
       ggplot(paramsitero, aes(stage, survival, colour = as.factor(a2), group = interaction(a2,b2)))+
         geom_line()+
         geom_point()+
         theme_bw()+
         facet_wrap(~b2)
       # ,width=300, height=300,units='mm', dpi=300)
  
paramlist_type1_semel <- params[params$stage == "x4" & params$survival < 0.2,]
paramssemel <-  params[interaction(params[,1:2]) %in% interaction(paramlist_type1_semel[,1:2]),]
(semelmin <- aggregate(survival ~ stage, min, data = params[interaction(params[,1:2]) %in% interaction(paramlist_type1_semel[,1:2]),]))
(semelmax <- aggregate(survival ~ stage, max, data = params[interaction(params[,1:2]) %in% interaction(paramlist_type1_semel[,1:2]),]))
# alpha2
min(paramssemel$a2)
max(paramssemel$a2)

min(paramssemel$b2)
max(paramssemel$b2)

# ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveTypeISemelParams.jpg",
       ggplot(paramssemel, aes(stage, survival, colour = as.factor(a2), group = interaction(a2,b2)))+
         geom_line()+
         geom_point()+
         theme_bw()+
         facet_wrap(~b2)
       # ,width=300, height=300,units='mm', dpi=300)

 semel_fecundsurv(semelmin$survival[semelmin$stage == "x4"])
 semel_fecundsurv(semelmax$survival[semelmax$stage == "x4"]) 
 itero_fecundsurv(slowmin$survival[slowmin$stage == "x4"]) 
 itero_fecundsurv(slowmax$survival[slowmax$stage == "x4"]) 


hazard <- function(alpha2,beta2,x){
  alpha2*exp(beta2*x)
}

plot(1:4, hazard(alpha2 = 1.2, beta2 = -.3, 1:4))
## ------------------ years 1:4 where reproductive when age 4------------------
# Early maturation so all are kept as ages 1:4 but then just stay as a reproductive after 4 years
# conceptual model of model creation and expections on lambda and extinction risk
jpeg("C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/ConceptualMatrixCreation.jpg",
     width = 170, height = 170, units="mm", res = 300)
      
      # Plot layout of plots
      layout(matrix(c(1,2,3),1,3))
      
      # Plot 1
            plot(1:4, survivalTypeI(0.03,0.60,1:4), 
           xaxt = "n",
           xlab = "Stage", ylab = "Survival (S)", type = "l")
      axis(1, at = 1:4)
      mtext("a)", side = 3, adj = 0)
      
      # Plot 2
      plot(seq(0,1,by=0.1), xlim = c(0,0.8), ylim = c(0,10), semel_fecundsurv(seq(0,1,by=0.1)), type = "l", 
           ylab = "Fecundity (R)", xlab = expression("Adult survival (S"[44]~")"))
      lines(seq(0,1,by=0.1), itero_fecundsurv(seq(0,1,by=0.1)))
      text(0.6,7, "semelparous")
      text(0.2,3, "iteroparous")
      mtext("b)", side = 3, adj = 0)
      
      # Plot 3
      plot(seq(0,0.9,0.25),seq(0,0.9,0.25),type = "n", ylim = c(0,.5), xlab = "Stage", ylab = "Growth (G)", xaxt = "n")
      axis(1, at = seq(0,0.9,by = 0.25), labels = 1:4)
      mtext("c)", side = 3, adj = 0)
      for(b1 in c(0,1)){
        for(b2 in c(0,1)){
          lines(seq(0,0.9,0.25), transitions(b1,b2,seq(0,0.91,0.25))/sum(transitions(b1,b2,seq(0,0.91,0.25)))) #, col = rgb((1-b2),(1-b1),b1,1))
        }}
      text(0.65,0.27, "fast")
      text(0.65, 0.165, "slow")

      
dev.off()

# Elasticity space, survival of the reproductive stage
# need to sum to one
## No make this after making lots of matrices
elast_type <- do.call(rbind,lapply(1:20, function(x){
  S <- runif(1,0.61,0.77)
  G <- runif(1,0.3,0.4)
  R <- itero_fecundsurv(S)
  itero <- data.frame(R = R/(R+S+G), G = G/(R+S+G), S = S/(R+S+G), parity = "iteroparous", speed = "fast")
  
  S <- runif(1,0.61, 0.77)
  R <- semel_fecundsurv(S)
  G <- runif(1,0.2,0.3) 
  iteroslow <- data.frame(R = R/(R+S+G), G = G/(R+S+G), S = S/(R+S+G), parity = "iteroparous", speed = "slow")
  
  # Semel
  G <- runif(1,0.3,0.4)
  S <- runif(1,0.11,0.2)
  R <- semel_fecundsurv(S)
  semelfast <- data.frame(R = R/(R+S+G), G = G/(R+S+G), S = S/(R+S+G), parity = "semelparous", speed = "fast")
  
  S <- runif(1,0.11, 0.2)
  R <- semel_fecundsurv(S)
  G <- runif(1,0.2,0.3)
  semelslow <- data.frame(R = R/(R+S+G), G = G/(R+S+G), S = S/(R+S+G), parity = "semelparous", speed = "slow")
  
  out <- rbind(itero,iteroslow,semelfast,semelslow)
  }))

elast_type <- data.frame(R = c(runif(20, 0.01,0.1),
                               runif(20, 0.85, 0.99)),
                         G = c(runif(20, 0.1, 0.12),
                               runif(20, 0.1, 0.12)),
                         S = c(runif(20, 0.71, 0.97),
                               runif(20, 0.01, 0.1)), Type = c(rep("iteroparous",20),rep("semelparous",20)))
elast_type[1:3] <- apply(elast_type[1:3], 1, function(x) x/sum(x))

# ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/TernaryPlotConceptual.jpg",
       ggtern::ggtern(elast_type, aes(R, G, S, colour = Type))+ # colour = parity, shape = speed))+
         geom_point()+
         # scale_color_viridis_c()+
         # geom_polygon_closed()+
         theme_showarrows()+
         theme_clockwise()
       # , width=100, height = 100, units = 'mm', dpi = 300)
      

### Set the SD of rnorm for variability in part of the function
### spit out the most elastic as a check
## ie. Silvertown et al 1996 "Interpretation of Elasticity Matrices as an aid to teh managment of plant populations for conservation"
generic_mat <- matrix(c("L","G","G","G",
         " ","L","G","G",
         " "," ","L","G",
         "F"," "," ","L"), nrow = 4)
# Survival
which(generic_mat == "L")
# Growth
which(generic_mat == "G")
# Fecundity
which(generic_mat == "F")

# Semelparous elasticities:            growth > 0.5, survival < 0.5, fecundity > 0.5
# Iteroparous elasticities:            growth [0,1], survival [0,1], fecundity < 0.5
# trees and shrubs (slow iteroparous): growth < 0.5, survival [0.25,1], fecundity < 0.25

# ------------------ Iteroparous fast ------------------------
# only progressive, no retrogressive growth; r_ij = 0
MPM_iterofast <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_type1),1)
  s_s <- survivalTypeI(alpha2 = paramlist_type1[i,1], beta2 = paramlist_type1[i,2], 1:4)
  f <- itero_fecundsurv(s_s[4])
  # transitions for fast should be fairly flat
  t_t <- transitions(b1 = 0.9,b2 = 0.1, x = 1:3)/sum( transitions(b1 = 0.9,b2 = 0.1, x = 1:3))
  t_ij <- matrix(c(0, s_s[1]*sum(t_t[1]), s_s[1]*sum(t_t[2]),  s_s[1]*sum(t_t[3]),
                   0, s_s[2]*sum(t_t[1]), s_s[2]*sum(t_t[2]),  s_s[2]*sum(t_t[3]),
                   0, 0,                  s_s[3]*sum(t_t2[1:2]),s_s[3]*sum(t_t2[3]),
                   f, 0,                  0,                   s_s[4]),
                 nrow = 4)
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1)))
}

MPMs_itfast <- lapply(1:100, function(x) MPM_iterofast()[[1]])
generation.time(mean(MPMs_itfast))
hist(unlist(lapply(MPMs_itfast, function(x) lambda(x))), xlab = "lambda", main = "Iteroparous Fast")

#Elasticities
Elasts_itfast <- do.call(rbind,lapply(1:100, function(x) MPM_iterofast(StDev = 0.1)[[3]]))
ggtern::ggtern(Elasts_itfast, aes(R, G, S, colour = lam))+
  geom_point()+
  scale_color_viridis_c()+
  theme_showarrows()+
  theme_clockwise()



#### Semel
# ------------------ Iteroparous fast ------------------------
# only progressive, no retrogressive growth; r_ij = 0
MPM_semelfast <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_type1_semel),1)
  s_s <- survivalTypeI(alpha2 = paramlist_type1_semel[i,1], beta2 = paramlist_type1_semel[i,2], 1:4)
  f <- semel_fecundsurv(s_s[4])
  # transitions for fast should be fairly flat
  t_t <- transitions(b1 = 0.9,b2 = 0.1, x = 1:3)/sum( transitions(b1 = 0.9,b2 = 0.1, x = 1:3))
  # t_t2 <- transitions(b1 = 0.9,b2 = 0.1, x = 1:2)/sum( transitions(b1 = 0.9,b2 = 0.1, x = 1:2))
  t_ij <- matrix(c(0, s_s[1]*sum(t_t[1]), s_s[1]*sum(t_t[2]),  s_s[1]*sum(t_t[3]),
                   0, s_s[2]*sum(t_t[1]), s_s[2]*sum(t_t[2]),  s_s[2]*sum(t_t[3]),
                   0, 0,                  s_s[3]*sum(t_t2[1:2]),s_s[3]*sum(t_t2[3]),
                   f, 0,                  0,                   s_s[4]),
                 nrow = 4)
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1)))
}

MPMs_itfast <- lapply(1:100, function(x) MPM_iterofast()[[1]])
generation.time(mean(MPMs_itfast))
hist(unlist(lapply(MPMs_itfast, function(x) lambda(x))), xlab = "lambda", main = "Iteroparous Fast")

#Elasticities
Elasts_itfast <- do.call(rbind,lapply(1:100, function(x) MPM_iterofast(StDev = 0.1)[[3]]))
ggtern::ggtern(Elasts_itfast, aes(R, G, S, colour = lam))+
  geom_point()+
  scale_color_viridis_c()+
  theme_showarrows()+
  theme_clockwise()



# --------------------- Semelparous slow -------------------
MPM_iteroslow <- function(StDev = 0.05){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_type1),1)
  s_s <- survivalTypeI(alpha2 = paramlist_type1[i,1],beta2 = paramlist_type1[i,2],1:4)
  t_t2 <- transitions(b1 = 0.1,b2 = 0.9, x = 1:3)/sum( transitions(b1 = 0.1,b2 = 0.9, x = 1:3))
  f <- semel_fecundsurv(s_s[4])
  t_ij <- matrix(c(0, s_s[1],              0,                  0,
                   0, s_s[2]*(1-t_t2[2]),  s_s[2]*t_t2[2],     0,
                   0, 0,                   s_s[3]*(1-t_t2[3]), s_s[3]*t_t2[3],
                   f, 0,                   0,                  s_s[4]), 
                 nrow = 4)
  (lambda1 <- lambda(t_ij))
  (lambdarange <- abs(1-rnorm(1, 1, StDev)))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1)))
}

MPMs_itslow <- lapply(1:100, function(x) MPM_iteroslow(StDev = 0.05)[[1]])
generation.time(mean(MPMs_itslow))
hist(unlist(lapply(MPMs_itslow, function(x) lambda(x))), xlab = "lambda", main = "Iteroparous Slow")

#Elasticities
Elasts_itslow <- do.call(rbind,lapply(1:100, function(x) MPM_iteroslow(StDev = 0.1)[[3]]))

ggtern::ggtern(Elasts_itslow, aes(R, G, S, colour = lam))+ #, shape = Type))+
  geom_point()+
  scale_color_viridis_c(name = expression(lambda))+
  theme_showarrows()+
  theme_clockwise()

# Elasticities final
ggplot(Elasts_itslowall, aes(lam, fill = Type))+
  geom_density(alpha = 0.5)+
  theme_bw()

str(Elasts_itslowall)

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
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
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

