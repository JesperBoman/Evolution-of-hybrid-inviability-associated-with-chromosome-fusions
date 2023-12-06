# install devtools first to download packages from github
install.packages("devtools")
library(devtools)
#install.packages("bioconductor")

# use devtools to install QTLseqr
devtools::install_github("bmansfeld/QTLseqr")
#install.packages("QTLseqr")
library("QTLseqr")

library(ggplot2)
library(reshape2)

#Reading and formatting data

setwd("./Downloads")
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


F2ad_female_df <- dataF2_frq_wide[, c(1:2, 13:14, 27:28)]
F2ad_female_df$AD_REF <- round((F2ad_female_df$F2ad_female.Cat * F2ad_female_df$F2ad_female.Cat.Cov + F2ad_female_df$F2ad_female.Swe * F2ad_female_df$F2ad_female.Swe.Cov)/2)
F2ad_female_df$AD_ALT <- round(( (1-F2ad_female_df$F2ad_female.Cat) * F2ad_female_df$F2ad_female.Cat.Cov + (1-F2ad_female_df$F2ad_female.Swe) * F2ad_female_df$F2ad_female.Swe.Cov)/2)
F2ad_female_df$DP <- F2ad_female_df$AD_REF+F2ad_female_df$AD_ALT

F2ad_male_df <- dataF2_frq_wide[, c(1:2, 15:16, 29:30)]
F2ad_male_df$AD_REF <- round((F2ad_male_df$F2ad_male.Cat * F2ad_male_df$F2ad_male.Cat.Cov + F2ad_male_df$F2ad_male.Swe * F2ad_male_df$F2ad_male.Swe.Cov)/2)
F2ad_male_df$AD_ALT <- round(( (1-F2ad_male_df$F2ad_male.Cat) * F2ad_male_df$F2ad_male.Cat.Cov + (1-F2ad_male_df$F2ad_male.Swe) * F2ad_male_df$F2ad_male.Swe.Cov)/2)


F2ad_male_df$DP <- F2ad_male_df$AD_REF+F2ad_male_df$AD_ALT

DE_df <- dataF2_frq_wide[, c(1:2, 7:8, 21:22)]
DE_df$AD_REF <- round((DE_df$DE.Cat * DE_df$DE.Cat.Cov + DE_df$DE.Swe * DE_df$DE.Swe.Cov)/2)
DE_df$AD_ALT <- round(( (1-DE_df$DE.Cat) * DE_df$DE.Cat.Cov + (1-DE_df$DE.Swe) * DE_df$DE.Swe.Cov)/2)
DE_df$DP <- DE_df$AD_REF+DE_df$AD_ALT

DLDP_df <- dataF2_frq_wide[, c(1:2, 9:10, 23:24)]
DLDP_df$AD_REF <- round((DLDP_df$DLDP.Cat * DLDP_df$DLDP.Cat.Cov + DLDP_df$DLDP.Swe * DLDP_df$DLDP.Swe.Cov)/2)
DLDP_df$AD_ALT <- round(( (1-DLDP_df$DLDP.Cat) * DLDP_df$DLDP.Cat.Cov + (1-DLDP_df$DLDP.Swe) * DLDP_df$DLDP.Swe.Cov)/2)
DLDP_df$DP <- DLDP_df$AD_REF+DLDP_df$AD_ALT

alive_df<-F2ad_female_df[, c(1:2)]
alive_df$AD_REF.HIGH <- F2ad_female_df$AD_REF+F2ad_male_df$AD_REF
alive_df$AD_ALT.HIGH <- F2ad_female_df$AD_ALT+F2ad_male_df$AD_ALT
alive_df$DP.HIGH <- F2ad_female_df$DP+F2ad_male_df$DP

dead_df<-F2ad_female_df[, c(1:2)]
dead_df$AD_REF.LOW <- DE_df$AD_REF + DLDP_df$AD_REF
dead_df$AD_ALT.LOW <- DE_df$AD_ALT + DLDP_df$AD_ALT
dead_df$DP.LOW <- DE_df$DP+DLDP_df$DP

dead_df$AD_REF.LOW <- DE_df$AD_REF 
dead_df$AD_ALT.LOW <- DE_df$AD_ALT
dead_df$DP.LOW <- DE_df$DP

AvD_df<-cbind(alive_df, dead_df[, -(1:2)])

AvD_df$deltaSNP <- AvD_df$AD_REF.HIGH/AvD_df$DP.HIGH - AvD_df$AD_REF.LOW/AvD_df$DP.LOW


#Filtering used in paper
#minTotalDepth = 100,
#maxTotalDepth = 800,
#depthDifference = 400,
#minSampleDepth = 80,

AvD_df_filt <-
  filterSNPs(
    SNPset = AvD_df,
    minTotalDepth = 100,
    maxTotalDepth = 800,
    depthDifference = 400,
    minSampleDepth = 80,
    verbose = TRUE
  )
colnames(AvD_df_filt) <- c("CHROM", "POS", "AD_REF.HIGH", "AD_ALT.HIGH", "DP.HIGH", "AD_REF.LOW", "AD_ALT.LOW", "DP.LOW", "deltaSNP")

AvD_df_QTLseq <-runQTLseqAnalysis(AvD_df_filt,
                                  windowSize = 1.5e6,
                                  popStruc = "F2",
                                  bulkSize = c(156, 370),
                                  replications = 10000,
                                  intervals = c(90, 95, 99)
)


AvD_results_QTLseq <- getQTLTable(SNPset = AvD_df_QTLseq, method = "QTLseq", interval = 95)

str(AvD_results_QTLseq[AvD_results_QTLseq$nSNPs > 1,])
write.table(AvD_results_QTLseq[AvD_results_QTLseq$nSNPs > 1,], file="QTLseq_sim_CI95.list", sep = "\t", quote = F, row.names=F)


AvD_df_QTLseq$Chr_num <- as.numeric(gsub("Chr_", "", AvD_df_QTLseq$CHROM))
AvD_results_QTLseq$Chr_num <- as.numeric(gsub("Chr_", "", AvD_results_QTLseq$CHROM))


numChr=24
ggplot(AvD_df_QTLseq[AvD_df_QTLseq$Chr_num > numChr,], aes(x=POS, y=tricubeDeltaSNP, col=as.factor(Chr_num)))+geom_line()+
  geom_line(aes(y=CI_95), col="blue", alpha=0.4)+
  geom_line(aes(y=CI_99), col="blue", alpha=0.2)+
  geom_line(aes(y=-CI_90), col="blue", alpha=0.6)+
  geom_line(aes(y=CI_90), col="blue", alpha=0.6)+
  geom_line(aes(y=-CI_95), col="blue", alpha=0.4)+
  geom_line(aes(y=-CI_99), col="blue", alpha=0.2)+
  ylab("Alive - Dead")+
  xlab("Position")+
  theme_classic()+
  ylim(-0.2, 0.2)+
  geom_hline(yintercept=0, lty=2)+
  theme(legend.position="none")+
  scale_color_manual(values = rep(c("grey", "black"), 48 )) +
  facet_grid(~Chr_num, scales = 'free_x', space = 'free_x', switch='x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 0))+
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 0.2),
            alpha = 0.5, fill="orange",
            data = AvD_results_QTLseq[AvD_results_QTLseq$Chr_num > numChr & AvD_results_QTLseq$avgDeltaSNP>0 ,],
            inherit.aes = FALSE) +
  
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.2, ymax =0, fill=pop),
            alpha = 0.5, fill = "red",
            data = AvD_results_QTLseq[AvD_results_QTLseq$Chr_num > numChr  & AvD_results_QTLseq$avgDeltaSNP<0,],
            inherit.aes = FALSE) +
  theme(strip.text.x = element_text(size = 14), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14, colour="black"))

