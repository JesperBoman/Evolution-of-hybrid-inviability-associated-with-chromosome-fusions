#!/usr/bin/env Rscript

# Development time survival model ####

library(nadiv)
library(MCMCglmm)
library(QGglmm)
library(reshape2)

f0_f2_pedigree <- read.csv(file = "f0-f2_pedigree_and_data_EL.csv", header = T, sep = ";", dec =",", stringsAsFactors = F, na.strings=c(""," ", "NA"))


F2data<- f0_f2_pedigree[f0_f2_pedigree$Lifespan != "Egg_harvest",]
F2data$ID <- as.factor(F2data$ID)
F2data$Survival <- ifelse(F2data$Lifespan == "Imago", 1, 0)
F2data[1:13,]$Survival <- NA

data_long <- melt(F2data, id.vars=c("ID", "Sex", "Sire", "Dam", "Lifespan", "Inbreeding_coefficient", "Comment", "Egg_laying_day", "Survival"))
data_long$devtime <- as.double(data_long$value)
data_long$Stage <- as.character(data_long$variable)
data_long$Stage <- ifelse(data_long$Stage == "I_II", "II", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "II_III", "III",data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "III_IV", "IV", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "IV_V", "V", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "V_P", "P", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "P_Imago", "Imago", data_long$Stage)
data_long$Stage <- factor(data_long$Stage , levels=c("Egg_to_emergence", "II", "III", "IV", "V", "P", "Imago"))
data_long$animal <- data_long$ID

ped <- f0_f2_pedigree[f0_f2_pedigree$Lifespan != "Egg_harvest",1:3]
ped$ID <- as.factor(ped$ID)
ped$Sire <- as.factor(ped$Sire)
ped$Dam <- as.factor(ped$Dam)
ped <- prepPed(ped)


Ainv <- inverseA(ped)$Ainv

data_long <- data_long[!is.na(data_long$animal),]
data_long <- data_long[!is.na(data_long$Sex),]
data_long <- data_long[!is.na(data_long$Survival),]

#Skipping egg-to-emergence, there is quite some measurment error for Egg_to_emergence in relation to the length of development up until that point
#E.g. when a certain egg was laid on a particular day was not recorded (imagine the work!), so there will be quite some uncertainty in that measure
data_long$Stage <- as.character(data_long$Stage)
data_long_noEtE  <- data_long[data_long$Stage != "Egg_to_emergence",]
data_long_noEtE$Stage <- factor(data_long_noEtE$Stage , levels=c("II", "III", "IV", "V", "P", "Imago"))
data_long_noEtE<-data_long_noEtE[!is.na(data_long_noEtE$Stage),]

Prior_gauss4 <- list(R = list(V = 1, nu = 0.002),
                     G = list(G1 = list(V = diag(6), nu = 0.002, alpha.mu = rep(0, 6), alpha.V= diag(1, 6, 6)), G2 = list(V = 1, nu = 0.002,  alpha.mu = 0,  alpha.V = 1)))


modeldev4 <-MCMCglmm(scale(devtime) ~ Sex+as.factor(Survival),
                     random = ~ us(1+Stage):animal+animal, ginv = list(animal = Ainv),
                     data = data_long_noEtE, prior = Prior_gauss4, rcov =~ units, family = "gaussian", nitt=10^5, burnin=10^4,thin=100, saveX=T,saveZ=T,saveXL=T, pr=T )

save(modeldev4, file="modeldev4_noEtE_pexpand_mcmcGLMM.rda")

Prior_gauss4 <- list(R = list(V = 1, nu = 1e-6),
                     G = list(G1 = list(V = diag(6), nu = 1e-6), G2 = list(V = 1, nu = 1e-6)))


modeldev4 <-MCMCglmm(scale(devtime) ~ Sex+as.factor(Survival),
                     random = ~ us(1+Stage):animal+animal, ginv = list(animal = Ainv),
                     data = data_long_noEtE, prior = Prior_gauss4, rcov =~ units, family = "gaussian", nitt=10^5, burnin=10^4,thin=100, saveX=T,saveZ=T,saveXL=T, pr=T )

save(modeldev4, file="modeldev4_noEtE_uninformative_mcmcGLMM.rda")
