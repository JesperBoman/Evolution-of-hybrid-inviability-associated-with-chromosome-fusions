library(ggplot2)
library(reshape2)

#Reading and formatting data
setwd("Downloads")
dataF2_frq <- data.frame()
for(sample in c("BE", "DE", "DLDP", "F1ad", "F2ad_female", "F2ad_male", "F2RP")){
tmp<- read.table(paste(sample, ".Swe.fixed.diff_noDel.mq20.af6.SweFreq_AncFreq.pol", sep=""))
tmp$Sample <- paste(sample, "Swe", sep=".")
dataF2_frq<-rbind(dataF2_frq, tmp)
}

for(sample in c("BE", "DE", "DLDP", "F1ad", "F2ad_female", "F2ad_male", "F2RP")){
  tmp<- read.table(paste(sample, ".Cat.fixed.diff_noDel.mq20.af6.SweFreq_AncFreq.pol", sep=""))
  tmp$Sample <- paste(sample, "Cat", sep=".")
  dataF2_frq<-rbind(dataF2_frq, tmp)
}

colnames(dataF2_frq) <- c("Chromosome", "Position", "Coverage", "SWE_Frequency", "Swedish_allele", "Catalan_allele", "Ancestral_allele","ANC_Frequency", "Sample")


faidx <- read.table("Lsin_DToL.fasta.fai")
chr_file <- faidx[,1:2]
chr_file$Start <- 1 
colnames(chr_file) <- c("Chromosome", "End", "Start")
chr_file$Chr_num <- as.numeric(gsub("Chr_", "", chr_file$Chromosome))
chr_file$Chr_type <- ifelse(chr_file$Chr_num  == 2 | chr_file$Chr_num  == 3 | chr_file$Chr_num == 48, "Z", "A" )


dataF2_frq$region_ID <- paste(dataF2_frq$Chromosome, dataF2_frq$Position, sep="_")
dataF2_frq_wide <- dcast(dataF2_frq, Chromosome+Position+region_ID+Ancestral_allele~Sample, value.var = "SWE_Frequency")
dataF2_frq_wide_cov <- dcast(dataF2_frq, Chromosome+Position+region_ID~paste(Sample, ".Cov", sep=""), value.var = "Coverage")
dataF2_frq_wide <- cbind(dataF2_frq_wide, dataF2_frq_wide_cov[,4:17])
dataF2_frq_wide$Chr_num <- as.numeric(gsub("Chr_", "", dataF2_frq_wide$Chromosome))


#Kolmogorov-Smirnov test for the Alive vs Dead comparison
ks.results<-data.frame()
for(Chromosome in unique(dataF2_frq_wide$Chromosome)){
kd<-dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chromosome,]
ks<-ks.test(((kd$F2ad_female.Swe+kd$F2ad_female.Cat)*80+(kd$F2ad_male.Swe+kd$F2ad_male.Cat)*76)/(2*(76+80)),
  ((kd$DLDP.Swe+kd$DLDP.Cat)*72+(kd$DE.Swe+kd$DE.Cat)*298)/(2*(72+298)))
ks.results <- rbind(ks.results, cbind(Chromosome, p.val=ks$p.value, D.stat=ks$statistic))
}

ks.results$p.val <- as.numeric(ks.results$p.val)
ks.results$p.val.bonferroni <- ks.results$p.val*48
ks.results$bonferroni.sig <- ifelse(ks.results$p.val.bonferroni < 0.05, "SIG", "NOT_SIG")


#Plotting allele frequency distributions per chromosome for exploratory data analysis
ggplot(dataF2_frq_wide, aes(x=((F2ad_female.Swe+F2ad_female.Cat)*80+(F2ad_male.Swe+F2ad_male.Cat)*76)/(2*(76+80))))+geom_histogram(col="red", alpha=0.3)+
  geom_histogram(aes(x=((DLDP.Swe+DLDP.Cat)*72+(DE.Swe+DE.Cat)*298)/(2*(72+298))), col="black", alpha=0.3)+
  facet_wrap(~Chr_num, scales = 'free_y', nrow=8, ncol=6)+
  ylab("Count")+
  xlab("SWE allele frequency")
