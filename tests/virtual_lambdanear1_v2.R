library(popbio)
library(ggplot2)
library(devtools)
# library(Rcompadre)
library(popdemo)
library(Rage)
library(ggtern)
library(patchwork)
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

# decreasing 
plot(1:4, transitions(0.1,.9,1:4)/sum(transitions(0.1,.9,1:4)))

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


paramlist_type1_semel <- params[params$stage == "x4" & params$survival > 0 & 
                                             !is.na(params$survival) & params$survival < 0.2,]
paramssemel <-  params[interaction(params[,1:2]) %in% interaction(paramlist_type1_semel[,1:2]),]
ggplot(paramssemel, aes(stage, survival, colour = as.factor(a2), group = interaction(a2,b2)))+
  geom_line()+
  geom_point()+
  theme_bw()+
  facet_wrap(~b2)
(semelmin <- aggregate(survival ~ stage, min, data = params_semel_test[interaction(params_semel_test[,1:2]) %in% interaction(paramlist_type1_semel[,1:2]),]))
(semelmax <- aggregate(survival ~ stage, max, data = params_semel_test[interaction(params_semel_test[,1:2]) %in% interaction(paramlist_type1_semel[,1:2]),]))
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
elast_type <- data.frame(R = c(0.3, 0.2, 0.7, 0.025, 0.2, 0.01, 0.25, 0.05, 0.9,  0.4, 0.05, 0.6, 0.49, 0.01, 0.98, 0.4),
                         G = c(0.5, 0.4, 0.29, 0.025, 0.7, 0.6,  0.7,  0.65,0.05, 0.5, 0.7,  0.5, 0.49, 0.01, 0.01, 0.4),
                         S = c(0.2, 0.5, 0.01, 0.95,  0.1, 0.39, 0.05, 0.3, 0.05, 0.6, 0.25, 0.4, 0.02, 0.98, 0.01, 0.2),
                         Type = c("semelparous","iteroparous",
                                  "semelparous","iteroparous",
                                  "semelparous","iteroparous",
                                  "semelparous","iteroparous",
                                  "semelparous","iteroparous",
                                  "semelparous","iteroparous",
                                  "semelparous","iteroparous",
                                  "semelparous","iteroparous"))

ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/TernaryPlotConceptual.jpg",
       ggtern::ggtern(elast_type, aes(R, G, S, colour = Type, fill = Type))+ # colour = parity, shape = speed))+label = Type, 
         # geom_point(size = 30, shape = 2)+
         # scale_color_viridis_c()+
         # geom_text(size=3.5)+
         # geom_polygon()+
         geom_mean_ellipse()+
         theme_showarrows()+
         theme_clockwise()+
         ggtitle("d)")
         # theme(legend.position = "none")
       , width=100, height = 70, units = 'mm', dpi = 300)
      

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
# keep between 0.8 < lambda < 1.2
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
  t_ij_norm <- t_ij
  
  (targetlam <- rnorm(1, 1, 0.1))
  if(lambda1 > targetlam) { # allow more or less variability by speed of life, more variable for fastest
    (i <- mostelastic <- which(e_ij == max(e_ij)))
    for(k in seq(0.01,0.5, by = 0.01)){
      t_ij_norm[mostelastic] <- (1-k) * t_ij_norm[mostelastic]
      if(t_ij_norm[13] < 0) t_ij_norm[13] <- 0
      lambda_norm <- lambda(t_ij_norm)
      if(lambda_norm <= targetlam) break
      mostelastic <- which(popbio::elasticity(t_ij_norm) == max(popbio::elasticity(t_ij_norm))) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda1 < targetlam){
    i <- mostelastic <- which(e_ij == max(e_ij))
    for(k in seq(1.01,1.5, by=0.01)){
      t_ij_norm[mostelastic] <- k * t_ij_norm[mostelastic]
      lambda_norm <- lambda(t_ij_norm)
      if(lambda_norm >= targetlam) break
      mostelastic <- which(popbio::elasticity(t_ij_norm) == max(popbio::elasticity(t_ij_norm))) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  e_ij_norm <- popbio::elasticity(t_ij_norm)
  survivalElast_norm <- sum(e_ij_norm[which(generic_mat == "L")])
  growthElast_norm <- sum(e_ij_norm[which(generic_mat == "G")])
  fecundElast_norm <- sum(e_ij_norm[which(generic_mat == "F")])
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1), 
              t_ij_norm, e_ij_norm, data.frame(S = survivalElast_norm, G = growthElast_norm, R = fecundElast_norm, lam_norm = lambda_norm)))
}

# get 100 runs
itfast <- lapply(1:100, function(x) MPM_iterofast())

MPMs_itfast <- lapply(itfast, function(x) x[[1]]) # The first item is the MPM
lamitfast <-data.frame(lam = unlist(lapply(MPMs_itfast, function(x) lambda(x))), parity = "itero", speed = "fast") 
gentimitfast <- data.frame(gentim = unlist(lapply(MPMs_itfast, function(x) generation.time(x))), parity = "itero", speed = "fast")

MPMs_itfast_norm <- lapply(itfast, function(x) x[[4]])
lamitfast_norm <-data.frame(lam = unlist(lapply(MPMs_itfast_norm, function(x) lambda(x))), parity = "itero", speed = "fast") 
gentimitfast_norm <- data.frame(gentim = unlist(lapply(MPMs_itfast_norm, function(x) generation.time(x))), parity = "itero", speed = "fast")

generation.time(mean(MPMs_itfast))
hist(gentimitfast$gentim, xlab = "Generation time", main = "Iteroparous Fast")
hist(lamitfast$lam, xlab = "lambda", main = "Iteroparous Fast")


# ---------------------------------------------

layout(matrix(c(1,2,3,4), 2,2),widths = c(6,2), heights = c(2,6))

  # Plot 1 density generation time
    den <- density(gentimitfast_norm$gentim)
    den2 <- density(gentimitfast$gentim)
    par(mar=c(0,4,0,0))
    plot(den$x, den$y, xlab = "",ylab="", main="", xaxt = "n", yaxt = "n", type = "l",bty="n", col = "red", ylim=c(0,max(c(den$y,den2$y))))
    lines(den2$x, den2$y, xlab = "",ylab="", main="", xaxt = "n", yaxt = "n", type = "l",bty="n", col = "lightcoral", lty = 5)
  # Plot 2 scatter plot
    par(mar=c(4,4,0,0))
    plot(gentimitfast_norm$gentim,lamitfast_norm$lam, ylab = expression(lambda), xlab = "generation time",
         main = "", pch = 16, col = "red")
    points(gentimitfast$gentim, lamitfast$lam, col = "lightcoral", pch = 1)
  # Plot 3 blank
    frame()
  # Plot 4 lambda density
    den <- density(lamitfast_norm$lam)
    den2 <- density(lamitfast$lam)
    par(mar=c(4,0,0,0))
    plot(den$y, den$x, xlab = "",ylab="", main="", xaxt = "n", yaxt = "n", type = "l",bty="n", col = "red", xlim=c(0,max(c(den$y,den2$y))))
    lines(den2$y, den2$x, xlab = "",ylab="", main="", xaxt = "n", yaxt = "n", bty="n", col = "lightcoral", lty = 5)
    

    dev.off()
      
# ---------------------------------------------
dev.off() # to reset par
            
plot(density(gentimitfast$gentim),  main = "Iteroparous Fast", xlab = "Generation time")
plot(density(lamitfast$lam), main = "Iteroparous Fast", xlab = expression(lambda))

#Elasticities

Elasts_itfast <- data.frame(do.call(rbind,lapply(1:100, function(x) MPM_iterofast()[[3]])), 
                            GenTime = unlist(lapply(1:100, function(x) generation.time(MPM_iterofast()[[1]]))))

Elasts_itfast_norm <- data.frame(do.call(rbind,lapply(1:100, function(x) MPM_iterofast()[[4]])), 
                            GenTime = unlist(lapply(1:100, function(x) generation.time(MPM_iterofast()[[1]]))))

tern_itfast <- ggtern::ggtern(Elasts_itfast, aes(R, G, S, colour = lam))+ #, size = as.factor(floor(GenTime))))+ #lam))+
                  geom_point()+
                  scale_color_viridis_c(name = expression(lambda))+
                  theme_showarrows()+
                  theme_clockwise()
                  # stat_mean_ellipse()

ggplot(Elasts_itfast)+
  geom_density(aes(GenTime))



# ---------------------------- Iteroparous Slow ---------------------
MPM_iteroslow <- function(){
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
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1)))
}

MPMs_itslow <- lapply(1:100, function(x) MPM_iteroslow()[[1]])
lamitslow <-data.frame(lam = unlist(lapply(MPMs_itslow, function(x) lambda(x))), parity = "itero", speed = "slow") 
gentimitslow <- data.frame(gentim = unlist(lapply(1:100, function(x) generation.time(MPM_iteroslow()[[1]]))), parity = "itero", speed = "slow") 
generation.time(mean(MPMs_itslow))

ggplot(rbind(lamitslow,lamitfast), aes(lam, fill = speed))+
  geom_density(alpha = 0.5)

# Figure 1 -
plot(density(lamitslow$lam), col = "red", xlab = expression(lambda), xlim = c(min(rbind(lamitslow,lamitfast)$lam),
                                                                              max(rbind(lamitslow,lamitfast)$lam)),
     main = "Iteroparous")
lines(density(lamitfast$lam), col = "blue")
text(1.025,30, "slow")
text(0.9,6,"fast")
# ----------

# Figure 2 - 
plot(density(gentimitslow$gentim), col = "red", xlab = "Generation time", xlim = c(min(rbind(gentimitslow,gentimitfast)$gentim)-1,
                                                                              max(rbind(gentimitslow,gentimitfast)$gentim)+1),
     main = "Iteroparous")
lines(density(gentimitfast$gentim), col = "blue")
text(mean(gentimitslow$gentim),0.1, "slow")
text(mean(gentimitfast$gentim),0.1,"fast")
# ------------


# Elasticities
Elasts_itfast <- do.call(rbind,lapply(1:100, function(x) MPM_iterofast()[[3]]))
ggtern::ggtern(Elasts_itfast, aes(R, G, S, colour = lam))+
  geom_point()+
  scale_color_viridis_c(name = expression(lambda))+
  theme_showarrows()+
  theme_clockwise()

IteroAll <- rbind(data.frame(Elasts_itfast, parity = "iteroparous", speed = "fast"),
                  data.frame(Elasts_itslow, parity = "iteroparous", speed = "slow"))

ggtern::ggtern(IteroAll, aes(R, G, S, colour = lam, shape = interaction(speed,parity)))+
            geom_point()+
            scale_color_viridis_c(name = expression(lambda))+
            scale_shape_manual(name = "Speed and\n parity", values = c(1,2))+
            theme_showarrows()+
            theme_clockwise()

# ------------------ Semelparous fast ------------------------
# only progressive, no retrogressive growth; r_ij = 0
MPM_semelfast <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_type1_semel),1)
  s_s <- survivalTypeI(alpha2 = paramlist_type1[i,1], beta2 = paramlist_type1[i,2], 1:4)
  s_s <- s_s - s_s[4]
  f <- semel_fecundsurv(s_s[3])
  # transitions for fast should be fairly flat
  t_t <- transitions(b1 = 0.9,b2 = 0.1, x = 1:3)/sum( transitions(b1 = 0.9,b2 = 0.1, x = 1:3))
  t_ij <- matrix(c(0, s_s[1]*sum(t_t[1]), s_s[1]*sum(t_t[2]),   s_s[1]*sum(t_t[3]),
                   0, s_s[2]*sum(t_t[1]), s_s[2]*sum(t_t[2]),   s_s[2]*sum(t_t[3]),
                   0, 0,                  s_s[3]*sum(t_t2[2]),  s_s[3]*sum(t_t2[c(1,3)]),
                   f, 0,                  0,                    s_s[4]),
                 nrow = 4)
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1)))
}

MPMs_sefast <- lapply(1:100, function(x) MPM_semelfast()[[1]])
lamsefast <- data.frame(lam = unlist(lapply(MPMs_sefast, function(x) lambda(x))), parity = "semel", speed = "fast")
gentimsefast <- data.frame(gentim = unlist(lapply(1:100, function(x) generation.time(MPM_semelfast()[[1]]))), parity = "semel", speed = "fast")
generation.time(mean(MPMs_sefast))
# hist(, xlab = "lambda", main = "Semelparous Fast")

# Figure 1 -
plot(density(lamitslow$lam), col = "red", xlab = expression(lambda), xlim = c(min(rbind(lamitslow,lamitfast,lamsefast)$lam)-0.05,
                                                                              max(rbind(lamitslow,lamitfast,lamsefast)$lam)+0.05),
     main="")
lines(density(lamitfast$lam), col = "blue")
lines(density(lamsefast$lam), col = "blue", lty = 3)
# text(1.025,30, "slow")
# text(0.9,6,"fast")
legend(0.77,35, legend = c("iteroparous","semelparous"),
       lty = c(1,3), cex = 0.8)
legend(0.77,27, legend = c("fast","slow"), col = c("blue","red"), cex = 0.8, lty=1)
# ----------

# Figure 2 - 
plot(density(gentimitslow$gentim), col = "red", xlab = "Generation time", 
     xlim = c(min(rbind(gentimitslow,gentimitfast,gentimsefast)$gentim)-1,
              max(rbind(gentimitslow,gentimitfast,gentimsefast)$gentim)+1),
     ylim = c(0,0.8),
     main = "")
lines(density(gentimitfast$gentim), col = "blue")
polygon(density(gentimsefast$gentim), col = rgb(0,0,1,0.2), lty = 3)
legend(7.27,.8, legend = c("iteroparous","semelparous"),
       lty = c(1,3), cex = 0.8)
legend(10,0.6, legend = c("fast","slow"), col = c("blue","red"), cex = 0.8, lty=1)
# ------------


#Elasticities
Elasts_sefast <- do.call(rbind,lapply(1:100, function(x) MPM_semelfast()[[3]]))
ggtern::ggtern(Elasts_sefast, aes(R, G, S, colour = lam))+
  geom_point()+
  scale_color_viridis_c(name = expression(lambda))+
  theme_showarrows()+
  theme_clockwise()

SemelIteroAll <- rbind(IteroAll, data.frame(Elasts_sefast, parity = "semelparous", speed = "fast"))
ggtern::ggtern(SemelIteroAll, aes(R, G, S, colour = lam, shape = interaction(speed,parity)))+
  geom_point(size = 2)+
  scale_color_viridis_c(name = expression(lambda))+
  scale_shape_manual(name = "Speed and\n parity", values = 1:3)+
  theme_showarrows()+
  theme_clockwise()


# --------------------- Semelparous slow -------------------
MPM_semelslow <- function(){
  # three juvenile stages, one reproductive
  i <- sample(1:nrow(paramlist_type1_semel),1)
  s_s <- survivalTypeI(alpha2 = paramlist_type1_semel[i,1], beta2 = paramlist_type1_semel[i,2], 1:4)
  f <- semel_fecundsurv((s_s[4]-.1))  
  t_t2 <- transitions(b1 = 0.9,b2 = 0.1, x = 1:3)/sum( transitions(b1 = 0.9,b2 = 0.1, x = 1:3))
  t_ij <- matrix(c(0, s_s[1],              0,                  0,
                   0, s_s[2]*sum(t_t2[c(1,3)]),  s_s[2]*t_t2[2],     0,
                   0, 0,                   s_s[3]*sum(t_t2[1:2]), s_s[3]*t_t2[3],
                   f, 0,                   0,                  0), 
                 nrow = 4)
  (lambda1 <- lambda(t_ij))
  (e_ij <- popbio::elasticity(t_ij))
  survivalElast <- sum(e_ij[which(generic_mat == "L")])
  growthElast <- sum(e_ij[which(generic_mat == "G")])
  fecundElast <- sum(e_ij[which(generic_mat == "F")])
  (targetlam <- rnorm(1, 1, 0.1))
  
  if(lambda1 > targetlam) { # allow more or less variability by speed of life, more variable for fastest
    (i <- mostelastic <- which(e_ij == max(e_ij)))
    for(k in seq(0.01,0.5, by = 0.01)){
      t_ij[mostelastic] <- (1-k) * t_ij[mostelastic]
      if(t_ij[13] < 0) t_ij[13] <- 0
      lambda1 <- lambda(t_ij)
      if(lambda1 <= targetlam) break
      mostelastic <- which(e_ij == max(e_ij)) # check if still same element
      if(mostelastic != i) break
    } # end reduction by k for loop of that element
  } # end while lambda is too big
  if(lambda1 < targetlam){
    i <- mostelastic <- which(e_ij == max(e_ij))
    for(k in seq(1.01,1.5, by=0.01)){
      t_ij[mostelastic] <- k * t_ij[mostelastic]
      lambda1 <- lambda(t_ij)
      if(lambda1 >= targetlam) break
      mostelastic <- which(e_ij == max(e_ij)) # check if still same element
      if(mostelastic !=i) break
    } # end reduction by i for loop of that element
  } # end if lambda is too small
  return(list(t_ij,e_ij, data.frame(S = survivalElast, G = growthElast, R = fecundElast, lam = lambda1)))
}

MPMs_seslow <- lapply(1:100, function(x) MPM_semelslow()[[1]])
generation.time(mean(MPMs_seslow))
lamseslow <- data.frame(lam = unlist(lapply(MPMs_seslow, function(x) lambda(x))), parity = "semel", speed = "slow")
gentimseslow <- data.frame(lam = unlist(lapply(MPMs_seslow, function(x) generation.time(x))), parity = "semel", speed = "slow")

#Elasticities
Elasts_seslow <- do.call(rbind,lapply(1:100, function(x) MPM_semelslow()[[3]]))
ggtern::ggtern(Elasts_seslow, aes(R, G, S, colour = lam))+ #, shape = Type))+
  geom_point()+
  scale_color_viridis_c(name = expression(lambda))+
  theme_showarrows()+
  theme_clockwise()

# Figure 1 -
plot(density(lamitslow$lam), col = "red", xlab = expression(lambda), xlim = c(min(rbind(lamitslow,lamitfast,lamsefast,lamseslow)$lam)-0.05,
                                                                              max(rbind(lamitslow,lamitfast,lamsefast,lamseslow)$lam)+0.05),
     main="")
lines(density(lamitfast$lam), col = "blue")
lines(density(lamsefast$lam), col = "blue", lty = 3)
lines(density(lamseslow$lam), col = "red", lty = 3)
# text(1.025,30, "slow")
# text(0.9,6,"fast")
legend(0.73,35, legend = c("iteroparous","semelparous"),
       lty = c(1,3), cex = 0.73)
legend(0.73,27, legend = c("fast","slow"), col = c("blue","red"), cex = 0.73, lty=1)
# ----------

# Figure 2 - 
plot(density(gentimitslow$gentim), col = "red", xlab = "Generation time", 
     xlim = c(min(rbind(gentimitslow,gentimitfast,gentimsefast,gentimseslow)$gentim)-1,
              max(rbind(gentimitslow,gentimitfast,gentimsefast,gentimseslow)$gentim)+1),
     ylim = c(0,0.8),
     main = "")
lines(density(gentimitfast$gentim), col = "blue")
polygon(density(gentimsefast$gentim), col = rgb(0,0,1,0.2), lty = 3)
lines(density(gentimseslow$gentim), col = rgb(1,0,0,0.2), lty = 3, lwd = 2)
legend(7.27,.8, legend = c("iteroparous","semelparous"),
       lty = c(1,3), cex = 0.8)
legend(10,0.6, legend = c("fast","slow"), col = c("blue","red"), cex = 0.8, lty=1)
# ------------

# Elasticities final
SemelIteroAll <- rbind(SemelIteroAll, data.frame(Elasts_seslow, parity = "semelparous", speed = "slow"))
table(SemelIteroAll$parity,SemelIteroAll$speed)
ggtern::ggtern(SemelIteroAll, aes(R, G, S, colour = lam, shape = interaction(speed,parity)))+
  geom_point(size = 2)+
  scale_color_viridis_c(name = expression(lambda))+
  scale_shape_manual(name = "Speed and\n parity", values = 1:4)+
  theme_showarrows()+
  theme_clockwise()

# --------------------- Annuals ------------------------
MPM_annual <- function(){
  # three seed stages, one reproductive
  i <- sample(1:nrow(paramlist_type1_semel),1)
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

