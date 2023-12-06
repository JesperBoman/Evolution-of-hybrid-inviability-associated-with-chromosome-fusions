# MCMCglmm - Survival model ####


library(nadiv)
library(MCMCglmm)
library(QGglmm)

#Additive survival model ####

#Example for the uninformative prior
load("model1.1_uninfprior_survival_mcmcGLMM.rda")

model1.1 <- model1.1_uninf
summary(model1.1)

mu <- mean(model1.1[["Sol"]][ , "(Intercept)"])
va <- mean(model1.1[["VCV"]][ , "animal"])
vp <- mean(rowSums(model1.1[["VCV"]]))

QGparams(mu = mu, var.a = va, var.p = vp, model = "binom1.probit")


dfQG <- data.frame(mu = as.vector(model1.1[["Sol"]][, "(Intercept)"]),
                 va = as.vector(model1.1[["VCV"]][, "animal"]),
                 vp = rowSums(model1.1[["VCV"]]))

post <- do.call("rbind", apply(dfQG, 1, function(row){
  QGparams(mu = row[["mu"]], var.a = row[["va"]], var.p = row[["vp"]], model = "binom1.probit", verbose = FALSE)
}))

ggplot(post, aes(x=h2.obs))+geom_density(fill="yellow", alpha=0.4)+
  xlab(expression(italic("h"^"2")))+
  theme_classic()+
  ylab("Density")+
  scale_colour_manual(values = c("#882255", "#DDCC77"))+
  scale_x_continuous(expand = c(0, 0), limits=c(0,0.55))+
  scale_y_continuous(expand = c(0, 0), limits=c(0,6.5))+
  theme(aspect.ratio=1, legend.position = "right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=16))




# Development time survival model ####
load("modeldev4_noEtE_pexpand_mcmcGLMM.rda")
summary(modeldev4)
