
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


#Inferring candidate regions: ALIVE vs DEAD

p <- ggplot(dataF2_frq_wide, aes(x=Position, y=(((F2ad_female.Swe+F2ad_female.Cat)*83+(F2ad_male.Swe+F2ad_male.Cat)*77)/(2*(77+83)) - ((DLDP.Swe+DLDP.Cat)*72+(DE.Swe+DE.Cat)*298)/(2*(72+298))), col=as.factor(Chr_num)))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept=0.05, lty=2)+
  geom_hline(yintercept=-0.05, lty=2)+
  geom_smooth(formula=y ~ s(x, bs = "cs"), col="purple", fill="purple", method='gam')+
  ylab("Alive - Dead freq")+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(legend.position="none")+
  scale_color_manual(values = rep(c("grey", "black"), 48 )) +
  facet_grid(~Chr_num, scales = 'free_x', space = 'free_x', switch='x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 0))+
  theme(strip.text.x = element_text(size = 14), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14, colour="black"))

p2 <- ggplot_build(p)

#We extract the results from the generalized additive model (GAM) called by geom_smooth
#The GAM is a fitted curve based on the observed allele frequency differences along the chromosomes
GAM_data <-p2$data[[4]]

#Here is the cutoff in allele frequency difference for which the 95% CI shouldn't overlap (i.e. the mean allele frequency difference will be greater than this cutoff)
cutoff=0.05
GAM_data$Cutoff_sign <- ifelse(GAM_data$ymin > cutoff | GAM_data$ymax < -cutoff, "Y", "N")

sign_regions<-data.frame()
prevSign<-"N"
prevChr <-""


for(i in 1:length(GAM_data$Cutoff_sign)){
  
  if(GAM_data$Cutoff_sign[i] == "Y" & (prevSign=="N" | prevChr != GAM_data$PANEL[i])){
    
    x_start <-  ifelse(prevChr != GAM_data$PANEL[i], 1, (GAM_data$x[i]+ GAM_data$x[i-1])/2)
    
    pop <- ifelse(GAM_data$y[i] >0, "SWE", "CAT")
  }
  
  stop_criterion <- ifelse( is.na(GAM_data$Cutoff_sign[i+1] == "N" | GAM_data$PANEL[i+1] != GAM_data$PANEL[i]), T, (GAM_data$Cutoff_sign[i+1] == "N" | GAM_data$PANEL[i+1] != GAM_data$PANEL[i]))
  
  if( stop_criterion == T & GAM_data$Cutoff_sign[i] == "Y") {
    
    x_end <- ifelse(GAM_data$PANEL[i+1] == GAM_data$PANEL[i], (GAM_data$x[i]+GAM_data$x[i+1])/2, chr_file[chr_file$Chr_num == GAM_data$PANEL[i], ]$End )
    x_end <- ifelse(is.na(x_end), chr_file[chr_file$Chr_num == GAM_data$PANEL[i], ]$End, x_end)
    
    sign_regions <- rbind(sign_regions, cbind(Chr=GAM_data$PANEL[i],start=x_start, end=x_end, pop=pop))
  }
  prevChr <- GAM_data$PANEL[i]
  prevSign <- GAM_data$Cutoff_sign[i]
}

#Sorry for the confusion here with sign_regions instead of e.g. cand_regions. At first I called them significant regions but that name could be confusing.
sign_regions$start <- round(as.numeric(sign_regions$start))
sign_regions$end <- round(as.numeric(sign_regions$end))
sign_regions$sizeMB <- (sign_regions$end - sign_regions$start)/1000000
sign_regions$Chr_num <- as.numeric(sign_regions$Chr)
sign_regions$Chr_type <- ifelse(sign_regions$Chr_num == 2 | sign_regions$Chr_num == 3 | sign_regions$Chr_num == 48, "Z", "A" )


table(sign_regions$pop)
range(sign_regions$sizeMB)

write.table(sign_regions, file = paste("sign_regions_ADvDEAD_dedup_wave", cutoff, ".list", sep=""), sep = "\t", quote = F, row.names=F)

sg_prop <- aggregate(sizeMB ~ pop, sign_regions, sum )
sg_prop$Prop <- sg_prop$sizeMB/ (sum(chr_file$End-chr_file$Start)/10^6)
sg_prop

sg_prop <- aggregate(sizeMB ~ Chr_type, sign_regions, sum )
sg_prop$Chr_type_size <- ifelse(sg_prop$Chr_type == "A", (sum(chr_file[chr_file$Chr_type == "A",]$End-chr_file[chr_file$Chr_type == "A",]$Start)/10^6), (sum(chr_file[chr_file$Chr_type == "Z",]$End-chr_file[chr_file$Chr_type == "Z",]$Start)/10^6))
sg_prop$Prop <- sg_prop$sizeMB / sg_prop$Chr_type_size
sg_prop


