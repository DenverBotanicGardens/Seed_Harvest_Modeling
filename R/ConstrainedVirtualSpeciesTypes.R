# install.packages("PACLasso")
# library(PACLasso)
library(popbio)
library(ggplot2)
rm(list=ls())

#Trade-off between survival and fecundity from Takada and Kawai (2020) 
# f=-1.25(s+4.6)^2+36.45+ error

#   s <- runif(1,0,0.8) # adult survival
fecundity <- function(s, errorSD = 0.1){ 
  out <- -1.25*((s+4.6)^2) + 36.45 + rnorm(1, 0 ,errorSD)  
  out
  }

fecundity(0.5,0.1)


survivorship <- function(alpha1, x, alpha3, beta3){
  out <- exp(-alpha1*x) * exp((alpha3/beta3)*(1-exp(beta3 * x)))
  out
}

survivorship(.5,.3,.4,.5)

hazard <- function(alpha2, beta2, x){
  alpha2 * exp(beta2*x)
}


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

# ------------------------------------------------ constrain lambda ---------------------------
# functions
fecunditysurvival <- function(s, errorSD = 0.1){ 
  out <- -1.25*((s+4.6)^2) + 36.45 + rnorm(1, 0 ,errorSD)  
  out
}

# or fecundity from Fujiwara and Diaz-Lopez 2017: b(x)
fecundityProportionalSize <- function(R, L_inf, kapa, x){
  R*(L_inf * (1-exp(-kapa*x))^3)
}

# Fujiwara and Diaz-Lopez 2017
survivalTypeIII <- function(alpha1,x,alpha3,beta3){
  #Hazard
  # h(x) <- alpha3*exp((-beta3 * x))
  exp(-alpha1*x) * exp((alpha3/beta3)*(1-exp(beta3 * x)))
}

# Find lambda ????????
# Euler-Lotka equation; solve for 1!! 
# Euler_Lotka <- function(x, fecund){
#   fun <- l^(-x) * fecund *  exp(-alpha1*x) * exp((alpha3/beta3)*(1-exp(beta3 * x)))
#   f <- expression(fun)
#   D(f, 'l')
# }
# Euler_Lotka(1:4, fecunditysurvival(0.8))

# Constant Hazard and Increasing fecundity
surv <- survivalTypeIII(0.5, 1:4, 0,0.4)
plot(1:4,surv, type="l")
abline(lm(surv ~ c(1:4)), col="red")


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
ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveParameters.jpg",
       ggplot(toplot, aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
         geom_point()+
         geom_line()+
         facet_wrap(~as.factor(b3)+as.factor(a1))+
         theme_bw(), 
       width=300, height=300,units='mm', dpi=300)
# a3 0.2 and b3 
ggsave(filename =  "C:/Users/DePrengm/OneDrive - Denver Botanic Gardens/P drive/My Documents/UCDenver_phd/Dissertation/Chapter3/Figures/SurvivalCurveParameters2.jpg",
       ggplot(toplot[toplot$b3 < 0.7 &
                       toplot$a1 < 0.5 &
                       toplot$a3 < 0.3,], aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
         geom_point()+
         geom_line()+
         facet_wrap(~as.factor(b3)+as.factor(a1))+
         theme_bw(), 
       width=300, height=300,units='mm', dpi=300)
head(toplot)
aggregate(survival ~ a1+a3+b3, sum, data=toplot)
aggregate(survival ~ a1+a3+b3, min, data=toplot)

# range of iteroparaty as in larger survival of older ones when b3 < 0.5 & a1 < 0.4
 ggplot(toplot[toplot$a1 < 0.2 & toplot$a3 < 0.3 & toplot$b3 < 0.5,], aes(stage, survival, colour = as.factor(a3), group = interaction(a1,a3)))+
   geom_point()+
   geom_line()+
   facet_wrap(~as.factor(b3)+as.factor(a1))+
   theme_bw()
 
 ggplot(toplot[toplot$a1==0.1 & toplot$a3<0.3 & toplot$b3 < 0.4,], 
        aes(stage, survival, colour = as.factor(a3)))+
   geom_point()+
   facet_grid(~as.factor(b3))

# which combinations of parameters needed to keep survival of adults above 0.2
paramsSurv <- unique(toplot[toplot$stage == "x4" & toplot$survival > 0.2,c("a1","a3","b3")])


# ------------------ Iteroparous fast ------------------------
MPM_iterofast <- function(){
  # three juvenile stages, one reproductive
  s_s <- survivalTypeIII(alpha1 = runif(1,0,0.1), 1:4, alpha3 = runif(1, 0,0.2), beta3 = runif(1, 0,0.4))
  f <- fecunditysurvival(s_s[4])
  g_j <- s_s/sum(s_s) # scale to 1, transitions follow same pattern as survival most likely to climb one stage; same for dropping in 
  r_j <- s_s/sum(s_s) 
  t_ij <- matrix(c(0, s_s[1]*sum(g_j[1:2]), s_s[1]*g_j[3],       s_s[1]*g_j[4],
                   0, s_s[2]*sum(g_j[1:2]), s_s[2]*g_j[3],       s_s[2]*g_j[4],
                   0, 0,                    s_s[3]*sum(g_j[1:3]),s_s[3]*g_j[4],
                   f, 0,                    0,                   s_s[4]*sum(g_j[1:4])), nrow = 4)
  lambda1 <- lambda(t_ij)
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

MPMs_itfast <- lapply(1:100, function(x) MPM_iterofast())
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



# Calculate Lambda and 
elast <- popbio::elasticity(t_ij)
which(elast == max(elast), arr.ind = TRUE)
maxelast <- which(elast == max(elast))
lamb1 <- lambda(t_ij)

scaledvalues <- rep(1,16) 
scaledvalues[maxelast] <- scaledvalues[maxelast] * 0.5  
lambda(t_ij*scaledvalues) # dot matrix multiplication basically


#-------------------------------------------------------------------------------------------------

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


  