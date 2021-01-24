#First load the functions and data used

set.seed(1111)
source('functions/toolbox_reorganization.R') #new set of tools
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}
#source('analysis/functions/function.R') #functions
#source('analysis/functions/alphantw.R') #Plot aphas ntw
#source('analysis/functions/triangle_projection_3sp_with_pairwise.R') #plotting
load(file = "data/matrices/plants.RData") #plant matrices as calculayted in plant_fit.R
load(file = "data/matrices/polinators.RData") #from fit_pol.R
load(file = "data/matrices/betas.Rdata") #Beta matrix as constructed in strstab.R 

#This overrides toolbox function to allow the null model to run smoothly
Probabilities_all <- function(alpha,S) {
  survival <- rep(0,S)
  sims <- 250  ### number of simualtions to obtain the probabilities of species persistence. In the text we used 100,000
  for(y in 1:sims){
    delta_t <- 0.01 #time step
    time_step <- seq(0,300,by=delta_t)
    N0 <- rep(1,S) #initial conditions for species abundances
    N0 <- N0 / sum(N0)
    blow <- 0
    while(blow==0){
      r <- sphere_sampling(S) ### sampling random K's in the unit sphere
      parms <- list(r=r, alpha = alpha) ##ODE
      model <- function(t,N,parms){ dN <- N * (parms$r - parms$alpha %*% N); list(dN)}
      e <- try(sol <- ode(N0,time_step,model,parms),silent = T)
      check <- (sol[nrow(sol),1+1]>0.000001)*1
      if (is.na(check)==FALSE){
        #sol <- ode(N0,time_step,model,parms)
        for(z in 1:S){
          survival[z] <- survival[z] + (sol[nrow(sol),z+1]>0.000001)*1 ### check if the species survived in the simulation
        }
        blow <- 1
      }
    }
  }
  survival <- survival / sims ### calculate probabilities
  out <- survival
  return(out)
}


#Calculate individual probs PLANTS and POLLINATORS (assuming that r can be positive and negative):

probs_no_link <- Probabilities_all(beta_no_link,6)
probs_link <- Probabilities_all(beta_link,6)
probs_analitical <- Probabilities_all(beta_no_link_analitical2,6)

#print(c(probs_no_link))
#print(c(probs_link))
#print(c(probs_analitical))

original_strength_diff <- beta_link - beta_no_link
hist(original_strength_diff)
mean(original_strength_diff[which(original_strength_diff != 0)])
sd(original_strength_diff[which(original_strength_diff != 0)])

########### This is the random analysis  ##### you can skip it (20 hours in my laptop)

#sims <- 10
sims <- 250
dat <- matrix(0, nrow = sims, ncol = 6)
reference <- beta_no_link_analitical2
reference[reference>0] <- 1
reference[reference<0] <- -1

hist(beta_link)
mean(beta_link) 
sd(beta_link) 

for(i in 1:sims){
  print(i)
  #M <- abs(rnorm(36, mean = 0.1, sd = 0.25)) 
  M <- beta_no_link_analitical2
  for(j in 1:length(M)){
    if(beta_no_link_analitical2[j] == 0){ #preserve zeroes
      M[j] <- 0 
    } else {
      M[j] <- rnorm(1, beta_no_link_analitical2[j], #add a random variation to each cell
                    sd = sd(original_strength_diff[which(original_strength_diff != 0)])) 
    }
  }
  #M <- runif(36)*2
  #M <- matrix(M, ncol = 6, nrow = 6)*reference
  #domain <- 0
  dat[i,] <- Probabilities_all(M,6)
}

dat2 <- dat
#save(dat2, file = "reorganization_random_nacho.RData")
#############################


load("data/reorganization_random_nacho.RData")
dat <- dat2
head(dat)
print(original_diff <- probs_no_link - probs_link) 
print(mean_random_diff <- colMeans(dat) - probs_link) 


#order T           R             H          O            B            F
par(mfrow = c(2,3))
hist(dat[,1], las = 1, main = "Tomato", xlab = "simulated")
abline(v = probs_no_link[1], col = "grey", lwd = 3)
abline(v = probs_link[1], col = "DodgerBlue", lty = 2, lwd = 3)
#abline(v = probs_analitical[1], col = "green")
abline(v = mean(dat[,1]), col = "orangered3", lwd = 3)
hist(dat[,2], las = 1, main = "Radish", xlab = "simulated")
abline(v = probs_no_link[2], col = "grey", lwd = 3)
abline(v = probs_link[2], col = "DodgerBlue", lty = 2, lwd = 3)
#abline(v = probs_analitical[2], col = "green")
abline(v = mean(dat[,2]), col = "orangered3", lwd = 3)
hist(dat[,3], las = 1, main = "Field Bean", xlab = "simulated")
abline(v = probs_no_link[3], col = "grey", lwd = 3)
abline(v = probs_link[3], col = "DodgerBlue", lty = 2, lwd = 3)
#abline(v = probs_analitical[3], col = "green")
abline(v = mean(dat[,3]), col = "orangered3", lwd = 3)
hist(dat[,4], las = 1, main = "Mason bee", xlab = "simulated")
abline(v = probs_no_link[4], col = "grey", lwd = 3)
abline(v = probs_link[4], col = "DodgerBlue", lty = 2, lwd = 3)
#abline(v = probs_analitical[4], col = "green")
abline(v = mean(dat[,4]), col = "orangered3", lwd = 3)
hist(dat[,5], las = 1, main = "Bumblebee", xlab = "simulated")
abline(v = probs_no_link[5], col = "grey", lwd = 3)
abline(v = probs_link[5], col = "DodgerBlue", lty = 2, lwd = 3)
#abline(v = probs_analitical[5], col = "green")
abline(v = mean(dat[,5]), col = "orangered3", lwd = 3)
hist(dat[,6], las = 1, main = "Green bottle fly", xlab = "simulated")
abline(v = probs_no_link[6], col = "grey", lwd = 3)
abline(v = probs_link[6], col = "DodgerBlue", lty = 2, lwd = 3)
#abline(v = probs_analitical[6], col = "green")
abline(v = mean(dat[,6]), col = "orangered3", lwd = 3)
par(mfrow = c(1,1))

#calculate p-val
pnorm(probs_no_link[1], mean = mean(dat[,1]), sd = sd(dat[,1]))
pnorm(probs_no_link[2], mean = mean(dat[,2]), sd = sd(dat[,2]))
pnorm(probs_no_link[3], mean = mean(dat[,3]), sd = sd(dat[,3]))
pnorm(probs_no_link[4], mean = mean(dat[,4]), sd = sd(dat[,4]))
pnorm(probs_no_link[5], mean = mean(dat[,5]), sd = sd(dat[,5]))
pnorm(probs_no_link[6], mean = mean(dat[,6]), sd = sd(dat[,6]))
#alternative for non gausian
(sum(probs_no_link[1] > dat[,1]))/length(dat[,1])
(sum(probs_no_link[2] > dat[,2]))/length(dat[,2])
(sum(probs_no_link[3] > dat[,3]))/length(dat[,3])
(sum(probs_no_link[4] > dat[,4]))/length(dat[,4])
(sum(probs_no_link[5] > dat[,5]))/length(dat[,5])
(sum(probs_no_link[6] > dat[,6]))/length(dat[,6])

#plot

original_diff <- probs_no_link - probs_link
mean_random_diff <- colMeans(dat) - probs_link

plot(x = 1,                 
     xlab = " ", ylab = "Probability",
     xlim = c(1, 12), ylim = c(0, 1),
     type = "n", las = 1)
segments(x0 = c(1,3,5,7,9,11), 
         y0 = probs_link, 
         x1 = c(2,4,6,8,10,12), 
         y1 = probs_analitical, 
         col = "orangered3")
segments(x0 = c(1,3,5,7,9,11), 
         y0 = probs_link, 
         x1 = c(2,4,6,8,10,12), 
         y1 = probs_no_link, 
         col = "darkorange")
segments(x0 = c(1,3,5,7,9,11), 
         y0 = probs_link, 
         x1 = c(2,4,6,8,10,12), 
         y1 = colMeans(dat), 
         col = "grey")
points(x = c(1,3,5,7,9,11), y = probs_link, pch = 16, col = "DodgerBlue")
points(x = c(2,4,6,8,10,12), y = probs_analitical, pch = 16, col = "orangered3")
points(x = c(2,4,6,8,10,12), y = probs_no_link, pch = 16, col = "darkorange")
points(x = c(2,4,6,8,10,12), y = colMeans(dat), pch = 16, col = "grey")

#in %
change_analitical <- (100-(probs_analitical*100)/probs_link)
change_nolink <- (100-(probs_no_link*100)/probs_link)
change_null <- (100-(colMeans(dat)*100)/probs_link)
plot(x = 1,                 
     xlab = " ", ylab = "effect (%)",
     xlim = c(1, 12), ylim = c(-100, 80),
     type = "n", las = 1)
segments(x0 = c(1,3,5,7,9,11), 
         y0 = change_analitical, 
         x1 = c(2,4,6,8,10,12), 
         y1 = change_nolink, 
         col = "grey")
segments(x0 = c(1,3,5,7,9,11), 
         y0 = change_analitical, 
         x1 = c(2,4,6,8,10,12), 
         y1 = change_null, 
         col = "grey")
points(x = c(1,3,5,7,9,11), y = change_analitical, pch = 16, col = "orangered3")
points(x = c(2,4,6,8,10,12), y = change_nolink, pch = 16, col = "darkorange")
points(x = c(2,4,6,8,10,12), y = change_null, pch = 16, col = "grey")
# segments(x0 = c(1,3,5,7,9,11), 
#          y0 = c(0,0,0,0,0,0), 
#          x1 = c(1,3,5,7,9,11), 
#          y1 = change_analitical, 
#          col = "orangered3")
# segments(x0 = c(2,4,6,8,10,12), 
#          y0 = c(0,0,0,0,0,0), 
#          x1 = c(2,4,6,8,10,12), 
#          y1 = change_nolink, 
#          col = "darkorange")
# segments(x0 = c(2,4,6,8,10,12), 
#          y0 = c(0,0,0,0,0,0), 
#          x1 = c(2,4,6,8,10,12), 
#          y1 = change_null, 
#          col = "grey")
abline(h = 0, lty = 3)



