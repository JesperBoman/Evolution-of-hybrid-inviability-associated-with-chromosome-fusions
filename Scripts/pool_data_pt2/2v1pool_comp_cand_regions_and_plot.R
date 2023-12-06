#This script can be used to obtain candidate regions when using pools of different sex ratio and you want to combine two pools in one group and one pool in another
#This was important for us since due to our pedigree, our expected allele frequency at Z sex chromosomes were 75 % SWE for males but 50 % SWE for females


# Sex ratio testing ####
Chr="Chr_48"
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F2ad_female.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F2ad_female.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F2ad_female.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F2ad_female.Cat.Cov, na.rm=T)-1

2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F2ad_male.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F2ad_male.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F2ad_male.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F2ad_male.Cat.Cov, na.rm=T)-1

2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F2RP.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F2RP.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F2RP.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F2RP.Cat.Cov, na.rm=T)-1

2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$DE.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$DE.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$DE.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$DE.Cat.Cov, na.rm=T)-1

2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$DLDP.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$DLDP.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$DLDP.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$DLDP.Cat.Cov, na.rm=T)-1

2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F1ad.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F1ad.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$F1ad.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$F1ad.Cat.Cov, na.rm=T)-1

2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$BE.Swe.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$BE.Swe.Cov, na.rm=T)-1
2*mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome == Chr,]$BE.Cat.Cov, na.rm=T)/mean(dataF2_frq_wide[dataF2_frq_wide$Chromosome != "Chr_48" & dataF2_frq_wide$Chromosome != "Chr_2" & dataF2_frq_wide$Chromosome != "Chr_3",]$BE.Cat.Cov, na.rm=T)-1




data_sex_rat_adj <- dataF2_frq_wide[,-grep("Cov", colnames(dataF2_frq_wide))]

group1 <- "F2ad_male"
group2 <- "F2ad_female"
group3 <- "F2RP"



group1_male_prop <- 1
group2_male_prop <- 0
group3_male_prop <- (0.322404+0.3219253)/2


group1_normalizer <- (1+group1_male_prop*0.5)
group2_normalizer <- (1+group2_male_prop*0.5)
group3_normalizer <- (1+group3_male_prop*0.5)

data_sex_rat_adj[, grep(group1, colnames(data_sex_rat_adj))[1]] <- ifelse(data_sex_rat_adj$Chr_num == 2 |data_sex_rat_adj$Chr_num == 3 | data_sex_rat_adj$Chr_num == 48,
                                                                          data_sex_rat_adj[, grep(group1, colnames(data_sex_rat_adj))[1]]/group1_normalizer,
                                                                          data_sex_rat_adj[, grep(group1, colnames(data_sex_rat_adj))[1]] )

data_sex_rat_adj[, grep(group1, colnames(data_sex_rat_adj))[2]] <- ifelse(data_sex_rat_adj$Chr_num == 2 |data_sex_rat_adj$Chr_num == 3 | data_sex_rat_adj$Chr_num == 48,
                                                                          data_sex_rat_adj[, grep(group1, colnames(data_sex_rat_adj))[2]]/group1_normalizer,
                                                                          data_sex_rat_adj[, grep(group1, colnames(data_sex_rat_adj))[2]] )

data_sex_rat_adj[, grep(group2, colnames(data_sex_rat_adj))[1]] <- ifelse(data_sex_rat_adj$Chr_num == 2 |data_sex_rat_adj$Chr_num == 3 | data_sex_rat_adj$Chr_num == 48,
                                                                          data_sex_rat_adj[, grep(group2, colnames(data_sex_rat_adj))[1]]/group2_normalizer,
                                                                          data_sex_rat_adj[, grep(group2, colnames(data_sex_rat_adj))[1]] )

data_sex_rat_adj[, grep(group2, colnames(data_sex_rat_adj))[2]] <- ifelse(data_sex_rat_adj$Chr_num == 2 |data_sex_rat_adj$Chr_num == 3 | data_sex_rat_adj$Chr_num == 48,
                                                                          data_sex_rat_adj[, grep(group2, colnames(data_sex_rat_adj))[2]]/group2_normalizer,
                                                                          data_sex_rat_adj[, grep(group2, colnames(data_sex_rat_adj))[2]] )

data_sex_rat_adj[, grep(group3, colnames(data_sex_rat_adj))[1]] <- ifelse(data_sex_rat_adj$Chr_num == 2 |data_sex_rat_adj$Chr_num == 3 | data_sex_rat_adj$Chr_num == 48,
                                                                          data_sex_rat_adj[, grep(group3, colnames(data_sex_rat_adj))[1]]/group3_normalizer,
                                                                          data_sex_rat_adj[, grep(group3, colnames(data_sex_rat_adj))[1]] )

data_sex_rat_adj[, grep(group3, colnames(data_sex_rat_adj))[2]] <- ifelse(data_sex_rat_adj$Chr_num == 2 |data_sex_rat_adj$Chr_num == 3 | data_sex_rat_adj$Chr_num == 48,
                                                                          data_sex_rat_adj[, grep(group3, colnames(data_sex_rat_adj))[2]]/group3_normalizer,
                                                                          data_sex_rat_adj[, grep(group3, colnames(data_sex_rat_adj))[2]] )



data_sex_rat_adj<-data_sex_rat_adj[, c("Chromosome", "Chr_num", "region_ID", "Position", colnames(data_sex_rat_adj)[grep(group1, colnames(data_sex_rat_adj))], colnames(data_sex_rat_adj)[grep(group2, colnames(data_sex_rat_adj))], colnames(data_sex_rat_adj)[grep(group3, colnames(data_sex_rat_adj))]) ]

colnames(data_sex_rat_adj) <- c("Chromosome", "Chr_num", "region_ID", "Position", "g1.1", "g1.2", "g2.1", "g2.2", "g3.1", "g3.2")

n.g1=77
n.g2=83

p <- ggplot(data_sex_rat_adj, aes(x=Position, y=((g1.1+g1.2)*n.g1+(g2.1+g2.2)*n.g2)/(2*(n.g1+n.g2)) - (g3.1+g3.2)/2, col=as.factor(Chr_num)))+
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
GAM_data <-p2$data[[4]]



cutoff=0.075
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

sign_regions$start <- round(as.numeric(sign_regions$start))
sign_regions$end <- round(as.numeric(sign_regions$end))
sign_regions$sizeMB <- (sign_regions$end - sign_regions$start)/1000000
sign_regions$Chr_num <- as.numeric(sign_regions$Chr)

sign_regions$Chr_type <- ifelse(sign_regions$Chr_num == 2 | sign_regions$Chr_num == 3 | sign_regions$Chr_num == 48, "Z", "A" )
sign_regions
write.table(sign_regions, file = paste("sign_regions_F2ADvF2RP_dedup_wave", cutoff, ".list", sep=""), sep = "\t", quote = F, row.names=F)


sg_prop <- aggregate(sizeMB ~ pop, sign_regions, sum )
sg_prop$Prop <- sg_prop$sizeMB/ (sum(chr_file$End-chr_file$Start)/10^6)

sg_prop <- aggregate(sizeMB ~ Chr_type, sign_regions, sum )
sg_prop$Chr_type_size <- ifelse(sg_prop$Chr_type == "A", (sum(chr_file[chr_file$Chr_type == "A",]$End-chr_file[chr_file$Chr_type == "A",]$Start)/10^6), (sum(chr_file[chr_file$Chr_type == "Z",]$End-chr_file[chr_file$Chr_type == "Z",]$Start)/10^6))
sg_prop$Prop <- sg_prop$sizeMB / sg_prop$Chr_type_size

table(sign_regions$pop)
sg_prop


numChr=24

ggplot(data_sex_rat_adj[data_sex_rat_adj$Chr_num > numChr,], aes(x=Position, y=(g1.1+g1.2+g2.1+g2.2)/4 - (g3.1+g3.2)/2, col=as.factor(Chr_num)))+
  geom_point(alpha=0.3)+
  geom_hline(yintercept=cutoff, lty=2)+
  geom_hline(yintercept=-cutoff, lty=2)+
  geom_smooth(col="purple", fill="purple")+
  ylab("Alive - Egg pool")+
  geom_hline(yintercept = 0)+
  ylim(-0.3,0.3)+
  geom_rect(aes(xmin = start, xmax = end, ymin = cutoff, ymax = 0.3),
            alpha = 0.3, fill="orange",
            data = sign_regions[sign_regions$Chr_num > numChr & sign_regions$pop == "SWE" ,],
            inherit.aes = FALSE) +
  
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.3, ymax = -cutoff, fill=pop),
            alpha = 0.3, fill = "red",
            data = sign_regions[sign_regions$Chr_num > numChr & sign_regions$pop == "CAT" ,],
            inherit.aes = FALSE) +
  theme_classic()+
  theme(legend.position="none")+
  scale_color_manual(values = rep(c("grey", "black"), 48 )) +
  facet_grid(~Chr_num, scales = 'free_x', space = 'free_x', switch='x')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 0))+
  theme(strip.text.x = element_text(size = 14), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=14, colour="black"))
