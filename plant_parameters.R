#This script calculates plant lambda and alpha matrices

#load functions and data needed----
rm(list=ls())
library(MASS)
library(DHARMa)

d <- read.csv(file='data/data_plants.csv',header = T)
str(d)
#d$Y <- (d$seeds  - (1-d$g) * d$s) / d$g
d$Y <- d$seeds
levels(d$treatment)

#First for no pollinators----

d2 <- d[d$treatment == 'NO_Pol',]
levels(d2$focal_plant)
str(d2)

d_T <- d2[d2$focal_plant == 'T',]
d_R <- d2[d2$focal_plant == 'R',]
d_H <- d2[d2$focal_plant == 'H',]

X_T <- as.matrix(d_T[4:6])
Y_T <- d_T$Y

X_R <- as.matrix(d_R[4:6])
Y_R <- d_R$Y

X_H <- as.matrix(d_H[4:6])
Y_H <- d_H$Y

m <- lm(log(Y_T) ~ X_T)
#plot(m) #highly skewed, estimates ~0,-0.18 and -0.44, the later two sig.
#summary(m)
m <- glm(round(Y_T, digits = 0) ~ X_T, family = "poisson") 
#plot(simulateResiduals(m)) #Bad fit
out_T <- glm.nb(round(Y_T, digits = 0) ~ X_T)
plot(simulateResiduals(out_T)) #good model.
summary(out_T) #Good fit, estimates ~0,-0.22,-0.51, congruent!

#For consistency, we use neg.bin for all models.
out_R <- glm.nb(round(Y_R, digits = 0) ~ X_R)
plot(simulateResiduals(out_R))
summary(out_R)

out_H <- glm.nb(round(Y_H, digits = 0) ~ X_H)
plot(simulateResiduals(out_H))
summary(out_H)

#Save r and alphas
r <- c(out_T$coefficients[1],out_R$coefficients[1],out_H$coefficients[1])
alpha <- -rbind(out_T$coefficients[2:4],out_R$coefficients[2:4],out_H$coefficients[2:4]) 
# Coefficient has been multiplied by -1 to change from the negative binomial fit
# to a alpha matrix that the sign of inter. coeff. are according to a L-V model. 

rownames(alpha) <- c('T','R','H')
colnames(alpha) <- c('T','R','H')
names(r) <- c('T','R','H')

# Second, with no link-----

d2 <- d[(d$treatment == 'no_link'),]
levels(d2$focal_plant)
str(d2)

d_T <- d2[d2$focal_plant == 'T',]
d_R <- d2[d2$focal_plant == 'R',]
d_H <- d2[d2$focal_plant == 'H',]

X_T_plant <- as.matrix(d_T[4:6])
X_T_pol <- as.matrix(d_T[7:9])
Y_T <- d_T$Y

X_R_plant <- as.matrix(d_R[4:6])
X_R_pol <- as.matrix(d_R[7:9])
Y_R <- d_R$Y

X_H_plant <- as.matrix(d_H[4:6])
X_H_pol <- as.matrix(d_H[7:9])
Y_H <- d_H$Y

#models
out_T <- glm.nb(round(Y_T, digits = 0) ~ X_T_plant + X_T_pol[,2])
plot(simulateResiduals(out_T))
summary(out_T)

out_R <- glm.nb(round(Y_R, digits = 0) ~ X_R_plant + X_R_pol[,c(1,3)]) #with 2 included, promts NA.
plot(simulateResiduals(out_R))
summary(out_R)

out_H <- glm.nb(round(Y_H, digits = 0)~ X_H_plant + X_H_pol[,c(1,2)])
plot(simulateResiduals(out_H))
summary(out_H)

r_no_link <- c(out_T$coefficients[1],out_R$coefficients[1],out_H$coefficients[1])
alpha_no_link <- -rbind(out_T$coefficients[2:4],out_R$coefficients[2:4],out_H$coefficients[2:4])
gamma_no_link <- matrix(0,nrow = 3,ncol = 3)
gamma_no_link[c(4,2,8,3,6)] <- c(out_T$coefficients[5],out_R$coefficients[5:6],out_H$coefficients[5:6])

rownames(alpha_no_link) <- c('T','R','H')
colnames(alpha_no_link) <- c('T','R','H')
rownames(gamma_no_link) <- c('T','R','H')
colnames(gamma_no_link) <- c('O','B','F')
names(r_no_link) <- c('T','R','H')

#Finally, link -----

d2 <- d[(d$treatment == 'link'),] 
d2 <- d2[d2$Y > 0,] # Remove individuals that did not fruits
levels(d2$focal_plant)
str(d2)

d_T <- d2[d2$focal_plant == 'T',]
d_R <- d2[d2$focal_plant == 'R',]
d_H <- d2[d2$focal_plant == 'H',]

X_T_plant <- as.matrix(d_T[4:6])
X_T_pol <- as.matrix(d_T[7:9])
Y_T <- d_T$Y

X_R_plant <- as.matrix(d_R[4:6])
X_R_pol <- as.matrix(d_R[7:9])
Y_R <- d_R$Y

X_H_plant <- as.matrix(d_H[4:6])
X_H_pol <- as.matrix(d_H[7:9])
Y_H <- d_H$Y

#models
out_T <- glm.nb(round(Y_T, digits = 0) ~ X_T_plant + X_T_pol[,2])
plot(simulateResiduals(out_T))
summary(out_T)

out_R <- glm.nb(round(Y_R, digits = 0) ~ X_R_plant + X_R_pol[,c(1,2,3)])
plot(simulateResiduals(out_R))
summary(out_R)

out_H <- glm.nb(round(Y_H, digits = 0) ~ X_H_plant + X_H_pol[,c(1,2)])
plot(simulateResiduals(out_H))
summary(out_H)

#omegas
r_link <- c(out_T$coefficients[1],out_R$coefficients[1],out_H$coefficients[1])
alpha_link <- -rbind(out_T$coefficients[2:4],out_R$coefficients[2:4],out_H$coefficients[2:4])
gamma_link <- matrix(0,nrow = 3,ncol = 3)
gamma_link[c(4,2,5,8,3,6)] <- c(out_T$coefficients[5],out_R$coefficients[5:7],out_H$coefficients[5:6])

rownames(alpha_link) <- c('T','R','H')
colnames(alpha_link) <- c('T','R','H')
rownames(gamma_link) <- c('T','R','H')
colnames(gamma_link) <- c('O','B','F')
names(r_link) <- c('T','R','H')

#wrap up----

r
r_link
r_no_link

alpha
alpha_link
alpha_no_link

gamma_link
gamma_no_link

#Export matrices----

# save(r,
#      r_link,
#      r_no_link,
#      alpha,
#      alpha_link,
#      alpha_no_link,
#      gamma_link,
#      gamma_no_link,
#      file = "data/matrices/plants.RData")
