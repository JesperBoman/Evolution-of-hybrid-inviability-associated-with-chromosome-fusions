#Data restructuring and visualization of phenotypic data

library(ggplot2)
library(reshape2)

f0_f2_pedigree <- read.csv(file = file.choose() , header = T, sep = ";", dec =",", stringsAsFactors = F, na.strings=c(""," ", "NA"))


data <- f0_f2_pedigree[f0_f2_pedigree$Dam == 20 | f0_f2_pedigree$Dam == 22 | f0_f2_pedigree$Dam == 23 | 
                         f0_f2_pedigree$Dam == 26 | f0_f2_pedigree$Dam == 27 | f0_f2_pedigree$Dam == 28 | f0_f2_pedigree$Dam == 29 | f0_f2_pedigree$Dam == 34, ]

data$Sex <- gsub("^F$", "Female", data$Sex)
data$Sex <- gsub("^M$", "Male", data$Sex)
data$Sex <- gsub("^M ", "Male", data$Sex)

lh_data <- as.data.frame(table(data$Lifespan, data$Dam))
colnames(lh_data) <- c("Stage", "Dam", "Freq")

lh_data$Stage <- as.character(lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_I", "I", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_II", "II", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_III", "III", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_IV", "IV", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_V", "V", lh_data$Stage)

#I also recorded transitional stages, but for simplicity we group them with the post-moult phenotype
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_I/II", "II", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_II/III", "III",lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_III/IV", "IV", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_IV/V", "V", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Instar_V/P", "Pupa", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Pupa/Imago", "Imago", lh_data$Stage)
lh_data$Stage <- ifelse(lh_data$Stage == "Egg", "Embryo", lh_data$Stage)

lh_data <- aggregate(Freq~Stage+Dam, lh_data, sum)
lh_data<-lh_data[lh_data$Stage != "Egg_harvest",]

lh_data$Stage <- factor(lh_data$Stage , levels=c("Embryo", "I", "II", "III","IV","V",  "Pupa", "Imago"))
lh_data <- lh_data[order(lh_data$Stage),]

#Get total number of offspring per female
tot_fit <- aggregate(Freq~Dam, lh_data, sum)


lh_data$PropFreq <- NA

for (i in 1:length(tot_fit$Dam)){
  matches <- tot_fit$Dam[i] == lh_data$Dam
  m_list <- (1:length(matches))[matches]
  lh_data$PropFreq[m_list] <-  lh_data$Freq[m_list]/tot_fit$Freq[i]
}

lh_data$Alive <- NA
lh_data$PropAlive <- NA

for(dam in tot_fit$Dam){
  alive <- tot_fit[tot_fit$Dam == dam,]$Freq
  for(stage in unique(lh_data$Stage)){
      lh_data[lh_data$Dam == dam & lh_data$Stage == stage, ]$Alive <- alive
      lh_data[lh_data$Dam == dam & lh_data$Stage == stage, ]$PropAlive <-  alive/tot_fit[tot_fit$Dam == dam,]$Freq
      alive <- alive - lh_data[lh_data$Dam == dam & lh_data$Stage == stage, ]$Freq
    }

}

#Proportion of offspring alive per female and stage
ggplot(data=lh_data, aes(y=PropAlive, x=Stage, group=Dam), col) +
  theme_classic()+
  ylab("Proportion offspring alive")+
  ylim(0,1)+
  stat_summary(fun=sum, geom="line", lwd=1)+
  theme(axis.text.x = element_text(angle = 45, hjust=1), aspect.ratio=1, legend.position = "right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=16))


# Proportion of offspring alive per sex and major developmental transition. 
#Sex of imagines have been phenotypically determined while sex of embryos and dead larvae + dead pupae were determined from resequencing data of pools####
male_survival<-c((1-(149/253))*100, (1-(27/(253-149)))*100)
female_survival<-c((1-(149/277))*100, (1-(45/(277-149)))*100)

survData <- as.data.frame(c(male_survival, female_survival))
survData <- cbind(survData, c("Embryo to Larva", "Larva to Imago"), c("M", "M", "F", "F"))
colnames(survData)<-c("Survival", "Stage", "Sex")

ggplot(survData, aes(y=Survival, x=Stage, col=Sex))+geom_point(size=5)+
  geom_segment(x="Embryo to Larva", xend="Larva to Imago", y=41.10672, yend=74.03846, col="#5D3A9B")+
  geom_segment(x="Embryo to Larva", xend="Larva to Imago", y=46.20939, yend=64.84375, col="#E66100")+
  ylim(0,100)+
  scale_color_manual(values=c("#E66100", "#5D3A9B"))+
  ylab("Survival (%)")+
  theme_classic()+
  theme(aspect.ratio=1, legend.position = "right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=16))




#Plot of development time split by lifespan ####
data_long <- melt(data, id.vars=c("ID", "Sex", "Sire", "Dam", "Lifespan", "Inbreeding_coefficient", "Comment", "Egg_laying_day"), 
                  measure.vars=c("Egg_to_emergence","I_II", "II", "II_III", "III", "III_IV", "IV", "IV_V", "V", "V_P", "P", "P_Imago", "Imago"))
data_long$value <- as.double(data_long$value)
data_long$Stage <- as.character(data_long$variable)


data_long$Sex <- gsub("Male ", "Male", data_long$Sex)

data_long$Survival <- ifelse(data_long$Lifespan == "Imago", "Alive", "Dead")

data_long$Stage <- ifelse(data_long$Stage == "I_II", "II", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "II_III", "III",data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "III_IV", "IV", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "IV_V", "V", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "V_P", "P", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "P_Imago", "Imago", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "P", "Pupa", data_long$Stage)
data_long$Stage <- ifelse(data_long$Stage == "Egg_to_emergence", "I", data_long$Stage)



data_long$Lifespan <- ifelse(data_long$Lifespan == "Egg", "Embryo", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_I", "I",data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_II", "II", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_III", "III", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_IV", "IV", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_V", "V", data_long$Lifespan)

data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_I/II", "II", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_II/III", "III",data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_III/IV", "IV", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_IV/V", "V", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Instar_V/P", "Pupa", data_long$Lifespan)
data_long$Lifespan <- ifelse(data_long$Lifespan == "Pupa/Imago", "Imago", data_long$Lifespan)

data_long$Stage <- factor(data_long$Stage , levels=c("I", "II", "III","IV",  "V",  "Pupa", "Imago"))
data_long <- data_long[!is.na(data_long$Survival),]
data_long$Lifespan <- factor(data_long$Lifespan , levels=c("Egg_harvest", "I", "II", "III","IV",  "V",  "Pupa", "Imago"))

data_long$Survival <- ifelse(data_long$Lifespan == "Imago", "Alive", "Dead")


ggplot(data=data_long[data_long$Lifespan != "Egg_harvest",], aes(x=value/(24*60), y=Stage, col=Survival, fill=Survival, group=Lifespan)) +
  stat_summary(fun=mean, geom="line", lwd=1)+
  stat_summary(fun.args=list(conf.int=0.95), geom="ribbon", alpha=0.3, lwd=1, colour=NA)+
  xlab("Days")+
  theme_classic()+
  scale_colour_manual(values = c("#882255", "#DDCC77"))+
  scale_fill_manual(values = c("#882255", "#DDCC77"))+
  theme(aspect.ratio=1, legend.position = "right", panel.border = element_rect(colour = "black", fill=NA, size=1), plot.title = element_text(face = "bold", hjust = 0.5), axis.text=element_text(size=14, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=16))




