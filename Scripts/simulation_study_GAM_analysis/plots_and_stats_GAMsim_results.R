library(ggplot2)

GAMsim.df<-read.table(file = file.choose(), header = T)


min(GAMsim.df[GAMsim.df$Loci == 500,]$False_positives)
max(GAMsim.df[GAMsim.df$Loci == 500,]$False_positives)

mean(GAMsim.df[GAMsim.df$Rec_rate == 0,]$False_positives)/mean(GAMsim.df[GAMsim.df$Rec_rate == 3,]$False_positives)

ggplot(GAMsim.df, aes(x=as.factor(Pool_fraction), y=False_positives/10, fill=as.factor(Rec_rate)))+geom_point(size=3, shape=21, position = position_dodge2(width=0.1))+
  ylab("False positives (%)")+
  xlab("Pool fraction")+
  theme_bw()+
  scale_fill_viridis_d(name="Recombination rate (cM/Mb)", option="plasma")+
  ylim(0,12)+
  facet_wrap(~Loci, nrow=1)+
  theme( element_text(face = "bold", hjust = 0.5), strip.text = element_text(size = 14),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=11, colour="black"), axis.title=element_text(size=16), legend.text=element_text(size=14),  legend.title=element_text(size=12))
