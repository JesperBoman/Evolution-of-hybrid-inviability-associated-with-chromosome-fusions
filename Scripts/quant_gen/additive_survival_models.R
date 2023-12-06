#!/usr/bin/env Rscript


library(nadiv)
library(MCMCglmm)
library(QGglmm)
library(reshape2)

f0_f2_pedigree <- read.csv(file = "f0-f2_pedigree_and_data_EL.csv", header = T, sep = ";", dec =",", stringsAsFactors = F, na.strings=c(""," ", "NA"))

ped <- f0_f2_pedigree[f0_f2_pedigree$Lifespan != "Egg_harvest",1:3]
F2data<- f0_f2_pedigree[f0_f2_pedigree$Lifespan != "Egg_harvest",]
F2data$ID <- as.factor(F2data$ID)
F2data$Survival <- ifelse(F2data$Lifespan == "Imago", 1, 0)
F2data[1:13,]$Survival <- NA

#F2data<-F2data[F2data$Dam != "22",]
ped <- F2data[,1:3]

ped$ID <- as.factor(ped$ID)
ped$Sire <- as.factor(ped$Sire)
ped$Dam <- as.factor(ped$Dam)
ped <- prepPed(ped)

Ainv <- inverseA(ped)$Ainv

F2data$animal <- F2data$ID

#Uninformative prior
Prior1.1 <- list(R = list(V = 1, fix = 1),
                   G = list(G1 = list(V = 1, nu = 1e-6)))
                   
model1.1 <- MCMCglmm(Survival ~ 1,
                     random = ~ animal, ginv = list(animal = Ainv),
                     data = F2data, prior = Prior1.1, family = "threshold", trunc=TRUE, nitt=10^6, burnin=10^4, thin=50)
                     
model1.1_uninf <- model1.1
save(model1.1_uninf, file="model1.1_uninfprior_survival_mcmcGLMM_1M.rda")


#Parameter-expanded prior
Prior1.1 <- list(R = list(V = 1, fix = 1),
                 G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V =1000)))
#http://cran.nexr.com/web/packages/MCMCglmm/vignettes/CourseNotes.pdf

model1.1 <- MCMCglmm(Survival ~ 1,
                     random = ~ animal, ginv = list(animal = Ainv),
                     data = F2data, prior = Prior1.1, family = "threshold", trunc=TRUE, nitt=10^6, burnin=10^4, thin=50)


model1.1_pexpand <- model1.1
save(model1.1_pexpand, file="model1.1_pexpandprior_survival_mcmcGLMM_1M.rda")
