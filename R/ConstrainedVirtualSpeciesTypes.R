library(popbio)
library(ggplot2)
library(dplyr)
rm(list=ls())

#Trade-off between survival and fecundity from Takada and Kawai (2020) 
# f=-1.25(s+4.6)^2+36.45+ error

#   s <- runif(1,0,0.8) for adult survival in Takada and Kawai 2020 
fecundity <- function(s, errorSD = 0.01){ 
  out <- -1.25*((s+4.6)^2) + 36.45 + rnorm(1, 0 ,errorSD)  
  out
  }

survivorship <- function(alpha1, x, alpha3, beta3){
  out <- exp(-alpha1*x) * exp((alpha3/beta3)*(1-exp(beta3 * x)))
  out
}

# just to visualize hazard changes over stage
hazard <- function(alpha2, beta2, x){
  alpha2 * exp(beta2*x)
}


# ------------------------------------------------ constrain lambda ---------------------------
# functions 
# Takada and Kawai 2020; s would be s_44
fecunditysurvival <- function(s, errorSD = 0.01){ 
  out <- -1.25*((s+4.6)^2) + 36.45 + rnorm(1, 0 ,errorSD)  
  out
}

# or fecundity from Fujiwara and Diaz-Lopez 2017: b(x); NO - only have one stage class that is reproductive, not applicable
fecundityProportionalSize <- function(R, L_inf, kapa, x){
  R*(L_inf * (1-exp(-kapa*x))^3)
}

# Fujiwara and Diaz-Lopez 2017 
# Do the x's represent the median age of the stage class???
survivalTypeIII <- function(alpha1,x,alpha3,beta3){
  #Hazard
  # h(x) <- alpha3*exp((-beta3 * x))
  exp(-alpha1*x) * exp((alpha3/beta3)*(1-exp(beta3 * x)))
}

survivalTypeIII(alpha1 = 0.01, seq(1.5,4.5), alpha3 = 0.04, beta3 = 0.3)
survivalTypeIII(alpha1 = 0.01, c(1.5,8.5,10.5,40.5), alpha3 = 0.04, beta3 = 0.3)
survivalTypeIII(alpha1 = 0.01, c(1.5,4.5,8.5,12.5), alpha3 = 0.04, beta3 = 0.3) # slow? longer lived? 


# Visualize what the survivalTypeIII equation is doing
toplot <- do.call(rbind, lapply(seq(0.1,0.9,by=0.1), function(a1){
  outa3 <- do.call(rbind, lapply(seq(0.1,0.9,by=0.1), function(a3){
    outb3 <- do.call(rbind, lapply(seq(0.1,0.9,by=0.1), function(b3){
      surv <- survivalTypeIII(alpha1 = a1, 1:4, alpha3 = a3, beta3 = b3)
       data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
     }))
     outb3
   }))
   outa3
 }))
head(toplot)
#
# ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveParameters.jpg",
      ggplot(toplot, aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
        geom_point()+
        geom_line()+
        facet_wrap(~as.factor(b3)+as.factor(a1))+
        theme_bw()
      # ,width=300, height=300,units='mm', dpi=300)
# a3 0.2 and b3 
# ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveParameters2.jpg",
       ggplot(toplot[toplot$b3 < 0.7 &
                       toplot$a1 < 0.5 &
                       toplot$a3 < 0.3,], aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
         geom_point()+
         geom_line()+
         facet_wrap(~as.factor(b3)+as.factor(a1))+
         theme_bw()
       # ,width=300, height=300,units='mm', dpi=300)

       
parametersfast <- do.call(rbind, lapply(seq(0,0.19,by=0.01), function(a1){ # a1 < 0.2
   outa3 <- do.call(rbind, lapply(seq(0,0.29,by=0.03), function(a3){ # a3 < 0.3
     outb3 <- do.call(rbind, lapply(seq(0,0.59,by=0.05), function(b3){ # b3 < 0.6
       surv <- survivalTypeIII(alpha1 = a1, 1:4, alpha3 = a3, beta3 = b3)
       data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
     }))
     outb3
   }))
   outa3
 }))
 
parametersslow <- do.call(rbind, lapply(seq(0.0,0.6,by=0.1), function(a1){ # a1 < 0.2
   outa3 <- do.call(rbind, lapply(seq(0.05,0.59,by=0.05), function(a3){ # a3 < 0.3
     outb3 <- do.call(rbind, lapply(seq(0.05,0.29,by=0.05), function(b3){ # b3 < 0.6
       surv <- survivalTypeIII(alpha1 = a1, c(1,4,10,16), alpha3 = a3, beta3 = b3)
       data.frame(a1, a3, b3, stage = c("x1","x2","x3","x4"), survival = surv)#x1 = surv[1], x2 = surv[2], x3 = surv[3], x4 = surv[4])
     }))
     outb3
   }))
   outa3
 }))

paramlist_slow <- parametersslow[parametersslow$stage == "x4" & parametersslow$survival > 0.01,]

ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveParametersslow.jpg",
       ggplot(parametersslow[interaction(parametersslow[,1:3]) %in% interaction(paramlist_slow[,1:3]),], 
              aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
         geom_point()+
         geom_line()+
         facet_wrap(~as.factor(b3)+as.factor(a1))+
         theme_bw() ,width=300, height=300,units='mm', dpi=300)

aggregate(survival ~ stage, min, data = parametersslow[interaction(parametersslow[,1:3]) %in% interaction(paramlist_slow[,1:3]),])
aggregate(survival ~ stage, max, data = parametersslow[interaction(parametersslow[,1:3]) %in% interaction(paramlist_slow[,1:3]),])

paramlist_fast <- parametersfast[parametersfast$stage == "x4" & parametersfast$survival > 0.2 & !is.na(parametersfast$survival) &
                                   parametersfast$survival < 0.5,]
aggregate(survival ~ stage, min, data = parametersfast[interaction(parametersfast[,1:3]) %in% interaction(paramlist_fast[,1:3]),])
aggregate(survival ~ stage, max, data = parametersfast[interaction(parametersfast[,1:3]) %in% interaction(paramlist_fast[,1:3]),])
ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveParametersfast.jpg",
       ggplot(parametersfast[interaction(parametersfast[,1:3]) %in% interaction(paramlist_fast[,1:3]),], 
       aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
  geom_point()+
  geom_line()+
  facet_wrap(~as.factor(b3)+as.factor(a1))+
  theme_bw() ,width=300, height=300,units='mm', dpi=300)
# SLOW
# which combinations of parameters needed to keep survival of adults above 0.2 - this would be slow, get to adult and stay adult a while
paramsSurvSlow <- unique(parameters[parameters$stage == "x4" & parameters$survival > 0.2,c("a1","a3","b3")])
plot(1:3, paramsSurvSlow[2,], type = "l", xaxt='n', xlab = "parameters", ylab = "parameter value combinations",
     ylim = c(0,0.6))
axis(1, at = c(1,2,3),labels = c("a1","a3","b3"))
for(i in 3:nrow(paramsSurvSlow)){
  lines(1:3, paramsSurvSlow[i,], col = i)
}
slow <- merge(parameters, paramsSurvSlow, by = c("a1","a3","b3"))
min(slow$survival[slow$stage == "x1"])
max(slow$survival[slow$stage == "x1"])

# FAST
# which combinations of parameters needed to keep survival of adults below 0.2 - this would be fast, get to adult and are likely to die, 
# lower iteroparity, but should have higher output of seed
paramsSurvfast <- unique(parameters[parameters$stage == "x4" & parameters$survival <= 0.2,c("a1","a3","b3")])
# What is the minimum stage 1 can be with these parameters?
plot(1:3, paramsSurvfast[2,], type = "l", xaxt='n', xlab = "parameters", ylab = "parameter value combinations",
     ylim = c(0,1))
axis(1, at = c(1,2,3),labels = c("a1","a3","b3"))
for(i in 3:nrow(paramsSurvfast)){
  lines(1:3, paramsSurvfast[i,], col = i, lty = i)
}
fast <- merge(parameters, paramsSurvfast, by = c("a1","a3","b3"))
min(fast$survival[fast$stage == "x1"])
max(fast$survival[fast$stage == "x1"])


# ------------------ Iteroparous fast ------------------------
# only progressive, no retrogressive growth; r_ij = 0
MPM_iterofast <- function(){
  # three juvenile stages, one reproductive
  # s_s <- survivalTypeIII(alpha1 = runif(1,0,0.1), 1:4, alpha3 = runif(1, 0,0.2), beta3 = runif(1, 0,0.4))
  i <- sample(1:nrow(paramsSurvfast),1)
  s_s <- survivalTypeIII(alpha1 = paramsSurvfast[i,1],1:4, alpha3 = paramsSurvfast[i,2], beta3 = paramsSurvfast[i,3] )
  f <- fecunditysurvival(s_s[4])
  # g_j <- s_s/sum(s_s) # scale to 1, transitions follow same pattern as survival most likely to climb one stage; same for dropping in 
  # r_j <- s_s/sum(s_s) 
  # t_ij <- matrix(c(0, s_s[1]*sum(g_j[1:2]), s_s[1]*g_j[3],       s_s[1]*g_j[4],
  #                  0, s_s[2]*sum(g_j[1:2]), s_s[2]*g_j[3],       s_s[2]*g_j[4],
  #                  0, 0,                    s_s[3]*sum(g_j[1:3]),s_s[3]*g_j[4],
  #                  f, 0,                    0,                   s_s[4]*sum(g_j[1:4])), nrow = 4)
  t_ij <- matrix(c(0, s_s[1]*(1/3), s_s[1]*(1/3),       s_s[1]*(1/3),
                   0, s_s[2]*(1/3), s_s[2]*(1/3),       s_s[2]*(1/3),
                   0, 0,                    s_s[3]*(1/2),s_s[3]*(1/2),
                   f, 0,                    0,                   s_s[4], nrow = 4)
    lambda1 <- lambda(t_ij)
  while(lambda1>1.2){ # allow more or less variability by speed of life, more variable for fastest rnorm(1, 1, 0.1) for limits
  i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    while(i == mostelastic){ # What if the most elastic transition changes? 
      for(i in seq(0.1,0.5, by = 0.1)){
        t_ij[mostelastic] <- (1-i) * t_ij[mostelastic]
        lambda1 <- lambda(t_ij)
        if(lambda1 <= 1.2) break
        mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
        } # end reduction by i for loop of that element
      } # end while same element is most elastic
    } # end while lambda is too big
  while(lambda(t_ij)<0.8){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    while(i == mostelastic){
      for(i in seq(1.1,1.5, by=0.1)){
        t_ij[mostelastic] <- i * t_ij[mostelastic]
        lambda1 <- lambda(t_ij)
        if(lambda1 >= 0.8) break
        mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
        } # end reduction by i for loop of that element
    } # end while same element is most elastic
  } # end while lambda is too small
  return(t_ij)
}

MPMs_itfast <- lapply(1:100, function(x) MPM_iterofast())
generation.time(mean(MPMs_itfast))
hist(unlist(lapply(MPMs_itfast, function(x) lambda(x))))

# ------------------------- just one, cover them all ------------------
MPM_virtual <- function(speed = "fast", parity = "iteroparous", annual = FALSE){
  # three juvenile stages, one reproductive
  
  s_s <- survivalTypeIII(alpha1 = runif(1,0,0.1), 1:4, alpha3 = runif(1, 0,0.2), beta3 = runif(1, 0,0.4))
  f <- fecunditysurvival(s_s[4])
  g_j <- s_s/sum(s_s) # scale to 1, transitions follow same pattern as survival most likely to climb one stage; same for dropping in 
  r_j <- s_s/sum(s_s) 
  # Make the matrix
  if(grepl(parity, "iteroparous")){
    t_ij <- matrix(c(0, s_s[1]*sum(g_j[1:2]), s_s[1]*g_j[3],       s_s[1]*g_j[4],
                     0, s_s[2]*sum(g_j[1:2]), s_s[2]*g_j[3],       s_s[2]*g_j[4],
                     0, 0,                    s_s[3]*sum(g_j[1:3]),s_s[3]*g_j[4],
                     f, 0,                    0,                   s_s[4]*sum(g_j[1:4])), nrow = 4)
    lambda1 <- lambda(t_ij)
  }
  if(grepl(parity, "semelparous")){
    t_ij <- matrix(c(0, s_s[1]*sum(g_j[1:2]), s_s[1]*g_j[3],       s_s[1]*g_j[4],
                     0, s_s[2]*sum(g_j[1:2]), s_s[2]*g_j[3],       s_s[2]*g_j[4],
                     0, 0,                    s_s[3]*sum(g_j[1:3]),s_s[3]*g_j[4],
                     f, 0,                    0,                   s_s[4]*sum(g_j[1:4])), nrow = 4)
    lambda1 <- lambda(t_ij)    
  }
  
  while(lambda1>1.2){ # allow more or less variability by speed of life, more variable for fastest
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    while(i == mostelastic){ # What if the most elastic stage changes? 
      for(i in seq(0.1,0.5, by = 0.1)){
        t_ij[mostelastic] <- (1-i) * t_ij[mostelastic]
        lambda1 <- lambda(t_ij)
        if(lambda1 <= 1.2) break
        mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      } # end reduction by i for loop of that element
    } # end while same element is most elastic
  } # end while lambda is too big
  while(lambda(t_ij)<0.8){
    i <- mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij)))
    while(i == mostelastic){
      for(i in seq(1.1,1.5, by=0.1)){
        t_ij[mostelastic] <- i * t_ij[mostelastic]
        lambda1 <- lambda(t_ij)
        if(lambda1 >= 0.8) break
        mostelastic <- which(popbio::elasticity(t_ij) == max(popbio::elasticity(t_ij))) # check if still same element
      } # end reduction by i for loop of that element
    } # end while same element is most elastic
  } # end while lambda is too small
  return(t_ij)
}

# ----------------------- Notes: ----------------------------------------------------------------

# Calculate Lambda and 
elast <- popbio::elasticity(t_ij)
which(elast == max(elast), arr.ind = TRUE)
maxelast <- which(elast == max(elast))
lamb1 <- lambda(t_ij)

scaledvalues <- rep(1,16) 
scaledvalues[maxelast] <- scaledvalues[maxelast] * 0.5  
lambda(t_ij*scaledvalues) # dot matrix multiplication basically


#-------------------------------------------------------------------------------------------------



# <https://rushinglab.github.io/WILD3810/articles/lab7_tree_harvest.html> 
eigen(t_ij)
# a negative value flips space, inverts space, absolute value is how scaled, first derivative is the 
# determinant is the volume, columns are linearly dependent if squishes to a plane, line, or point
# thumb = khat, pointer = ihat, middle = jhat
det(t_ij) # A - lambda*I where roots of polynomial are the eigenvalues, how much the linear transformation scales the unit area, 
# This makes it very small, if goes onto a line or point, all regions becomes zero, into a smaller dimension
solve(t_ij) # The inverse, so multiplying by the matirx and by the inverse gets back to original 
vec <- matrix(runif(4, 0,1), ncol = 1)
solve(t_ij) %*% (t_ij %*% vec)
solve(t_ij) %*% t_ij


# imaginay numbers means rotation-scaling matrix. change to eigenbasis, transform matrix to eigenbasis 
# eigenvectors
eigenvectors <- eigen(t_ij)$vectors
eigenbasis <- solve(eigenvectors) %*% t_ij %*% eigenvectors

# complex eigenvectors Trick for 2x2 matrices
# https://textbooks.math.gatech.edu/ila/complex-eigenvalues.html

# Determinant is how much the A tranformation matrix changes the unit area, is zero when dimension reduction 
# negative amounts just means flipped orientation, i hat is no longer going right while j hat is to the left and going up 
det(t_ij) # is negative but small -0.000037; rotation past zero by a little so reduced greatly. nearly zero, nearly linearly dependent
eigen(t_ij) # main eigenvalue is 2.96

# <https://www.youtube.com/watch?v=PFDu9oVAE-g&list=PLZHQObOWTQDPD3MizzM2xVFitgF8hE_ab&index=14> 
A <- matrix(c(0,1,1,1),nrow = 2)
# Compute powers of A by hand, sure looks like Fibonacci Sequence 
A2 <- A %*% A
A3 <- A %*% A %*% A 
A4 <- A %*% A %*% A %*% A
A5 <- A %*% A %*% A %*% A %*% A
A6 <- A %*% A %*% A %*% A %*% A %*% A
A7 <- A %*% A %*% A %*% A %*% A %*% A %*% A


eigenvector1 <- matrix(c(2, 1+sqrt(5)), ncol = 1)
eigenvector2 <- matrix(c(2, 1-sqrt(5)), ncol = 1)

# Change to an eigenbasis by transforming Tr^-1 %*% A %*% Tr where Tr = [eigenvector1 eigenvetor2]
Tr <- matrix(c(eigenvector1,eigenvector2), ncol = 2) 
eigenbasis <- solve(Tr) %*% A %*% Tr
lambda <- -1
det(eigenbasis - matrix(c(lambda,0,0,lambda), nrow = 2))
# ???? 
# 

popbio::lambda(t_ij)
# sensitivity change in lambda caused by small change in vital rate but elasticity takes proportional change in matrix element on lambda 
popbio::sensitivity(t_ij)
popbio::elasticity(t_ij) # This is what Wunder suggested to use to figure out what to change
popbio::stable.stage(t_ij)
popbio::generation.time(t_ij)

# Constant with different units depending on the type of fecundity
Rm

# Early maturation types
x1 <- x2 <- 1.5 

f <- expression(x^2 + 3*x)
D(f,'x') # first derivative
D(D(f, 'x'),'x')
# evaluate over values of x
x <- 1:5
eval(f)
eval(D(f,'x'))
# max is where the first derivative crosses the x from positive to negative

# integrate over a function from lower to upper limits
f <- function(x) x^2 + 3*x
integrate(f, 0, 1)


  