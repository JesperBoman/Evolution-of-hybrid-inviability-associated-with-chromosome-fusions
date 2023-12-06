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


sign_regions <- read.table(file=file.choose(), header=T) #Use candidate regions from Inferring_candidate_regions.R

AvD_results_QTLseq <- read.table(file=file.choose(), header=T)

numChr=25

ggplot(dataF2_frq_wide[dataF2_frq_wide$Chr_num < numChr,], aes(x=Position, y=(((F2ad_female.Swe+F2ad_female.Cat)*83+(F2ad_male.Swe+F2ad_male.Cat)*77)/(2*(77+83)) - ((DLDP.Swe+DLDP.Cat)*72+(DE.Swe+DE.Cat)*298)/(2*(72+298))), col=as.factor(Chr_num)))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept=cutoff, lty=2)+
  geom_hline(yintercept=-cutoff, lty=2)+
  geom_smooth(formula=y ~ s(x, bs = "cs"), col="purple", fill="purple", method='gam')+
  ylab("Alive - Dead")+
  geom_hline(yintercept = 0)+
  ylim(-0.3,0.3)+
  geom_rect(aes(xmin = start, xmax = end, ymin = cutoff, ymax = 0.3),
            alpha = 0.3, fill="orange",
            data = sign_regions[sign_regions$Chr_num < numChr & sign_regions$pop == "SWE" ,],
            inherit.aes = FALSE) +
  
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.3, ymax = -cutoff, fill=pop),
            alpha = 0.3, fill = "red",
            data = sign_regions[sign_regions$Chr_num < numChr & sign_regions$pop == "CAT" ,],
            inherit.aes = FALSE) +
  geom_point(data=AvD_results_QTLseq[AvD_results_QTLseq$Chr_num < numChr & AvD_results_QTLseq$nSNPs > 1,], aes(x=start+(end-start), y=0.29, group=as.factor(Chr_num)), shape=8, col="red", size=3)+
  theme_classic()+
  theme(legend.position="none")+
  scale_color_manual(values = rep(c("grey", "black"), 48 )) +
  facet_grid(~Chr_num, scales = 'free_x', space = 'free_x', switch='x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 0))+
  theme(strip.text.x = element_text(size = 14), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14, colour="black"))


