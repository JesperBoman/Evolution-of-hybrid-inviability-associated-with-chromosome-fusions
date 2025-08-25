#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
args<-as.numeric(args)

#args[1] = constant pool fraction
#args[2] = number of loci
#args[3] = recombination modifier

#Sample sizes of pools

#Alive
nAF=80*2
nAM=76*2

#Dead
nDE=298*2
nDLDP=72*2


#Constant pool fraction
pool_fraction<-args[1]

#Mean read depth
mean_nb <- 115

#Variance in read depth
var_nb<-8000

#Number of loci
loci=args[2]
physL=30

#Recombination probability
recP<-(args[3]*physL)/(100*loci) #Implicitly physL=genL, so roughly 60 cM with chiasmatic male meiosis. Sex-averaged 30 cM. This gives 1 cM/Mb, which is conservative




FDres=data.frame()
resdataFULL=data.frame()
cutoff=0.05


for(j in 1:1000){
  resdata<-data.frame()
  
  for(i in 1:loci){
    set.seed(i*j)
    if(i == 1){
      #Adult females
      seg_freq_AF<-mean(rbinom(nAF, 1, p=0.5))
      pool <- rbinom(nAF*pool_fraction, 1, p=seg_freq_AF)
      pool_freq<-mean(pool)
      #Code below is from PoolHelper (https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14185)
      # calculate the parameters for the negative binomial
      pnb <- mean_nb/var_nb
      rnb <- (mean_nb^2)/(var_nb - mean_nb)
      AF_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
      ###
      
      
      #Adult males
      seg_freq_AM<-mean(rbinom(nAM, 1, p=0.5))
      pool <- rbinom(nAM*pool_fraction, 1, p=seg_freq_AM)
      pnb <- mean_nb/var_nb
      rnb <- (mean_nb^2)/(var_nb - mean_nb)
      AM_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
      ###
      
      pAw<-((AF_freq*80)+(AM_freq*76))/(80+76)
      
      
      #Dead embryos
      seg_freq_DE<-mean(rbinom(nDE, 1, p=0.5))
      pool <- rbinom(nDE*pool_fraction, 1, p=seg_freq_DE)
      pnb <- mean_nb/var_nb
      rnb <- (mean_nb^2)/(var_nb - mean_nb)
      DE_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
      ###
      
      
      #Dead larvae + dead pupae
      seg_freq_DLDP<-mean(rbinom(nDLDP, 1, p=0.5))
      pool <- rbinom(nDLDP*pool_fraction, 1, p=seg_freq_DLDP)
      pnb <- mean_nb/var_nb
      rnb <- (mean_nb^2)/(var_nb - mean_nb)
      DLDP_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
      ###
      
      pDw<-((DE_freq*298)+(DLDP_freq*72))/(298+72)
      
      
    }
    else{
      n<-sum(runif(nAF)<=recP)
      if(n == 0){
        #Adult females
        pool <- rbinom(nAF*pool_fraction, 1, p=seg_freq_AF)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        AF_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      else{
        seg_freq_AF_new<-mean(rbinom(n, 1, p=0.5))
        seg_freq_AF<-((seg_freq_AF*(nAF-n))+(seg_freq_AF_new*n))/nAF
        #Adult females
        pool <- rbinom(nAF*pool_fraction, 1, p=seg_freq_AF)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        AF_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      
      
      n=sum(runif(nAM)<=recP)
      if(n == 0){
        #Adult males
        pool <- rbinom(nAM*pool_fraction, 1, p=seg_freq_AM)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        AM_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      else{
        seg_freq_AM_new<-mean(rbinom(n, 1, p=0.5))
        seg_freq_AM<-((seg_freq_AM*(nAM-n))+(seg_freq_AM_new*n))/nAM
        #Adult males
        pool <- rbinom(nAM*pool_fraction, 1, p=seg_freq_AM)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        AM_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      
      
      pAw<-((AF_freq*80)+(AM_freq*76))/(80+76)
      
      
      n=sum(runif(nDE)<=recP)
      if(n == 0){
        #Dead embryos
        pool <- rbinom(nDE*pool_fraction, 1, p=seg_freq_DE)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        DE_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      else{
        seg_freq_DE_new<-mean(rbinom(n, 1, p=0.5))
        seg_freq_DE<-((seg_freq_DE*(nDE-n))+(seg_freq_DE_new*n))/nDE
        #Dead embryos
        pool <- rbinom(nDE*pool_fraction, 1, p=seg_freq_DE)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        DE_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      
      
      n=sum(runif(nDLDP)<=recP)
      if(n == 0){
        #Dead larvae and dead pupae
        pool <- rbinom(nDLDP*pool_fraction, 1, p=seg_freq_DLDP)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        DLDP_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      else{
        seg_freq_DLDP_new<-mean(rbinom(n, 1, p=0.5))
        seg_freq_DLDP<-((seg_freq_DLDP*(nDLDP-n))+(seg_freq_DLDP_new*n))/nDLDP
        #Dead larvae and dead pupae
        pool <- rbinom(nDLDP*pool_fraction, 1, p=seg_freq_DLDP)
        pnb <- mean_nb/var_nb
        rnb <- (mean_nb^2)/(var_nb - mean_nb)
        DLDP_freq<-mean(sample(pool, rnbinom(n=1, size=rnb, prob=pnb), replace=T))
        ###
        
      }
      
      pDw<-((DE_freq*298)+(DLDP_freq*72))/(298+72)
      
    }
    
    #Save stuff here
    resdata<-rbind(resdata, cbind("pAw"=pAw, "pDw"=pDw))
    
  }
  
  resdata$iteration <- j
  resdataFULL <- rbind(resdataFULL, resdata)
  
  
  p<-ggplot(resdata, aes(x=seq(from=1, to=physL*10^6, length.out=loci), y=pAw-pDw))+geom_point()+
    geom_smooth(formula=y ~ s(x, bs = "cs"), method='gam')
  p2<-ggplot_build(p)
  GAM_data<-p2$data[[2]]
  
  FD<-as.data.frame(ifelse(GAM_data$ymin > cutoff | GAM_data$ymax < -cutoff, "Y", "N"))
  FD$iteration <- j
  
  
  FDres<-rbind(FDres, FD)
}

FDres.df<-as.data.frame(table(FDres[,1], FDres$iteration))
FDres.df$FP<-ifelse(FDres.df$Var1 == "Y" & FDres.df$Freq>0, "FP", "TN" )

print(noquote(paste(as.integer(table(FDres.df[FDres.df$Var1 == "Y",]$FP)[1]), args[1], args[2], args[3], sep =",")))



