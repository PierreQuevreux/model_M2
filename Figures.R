# Libraries ---------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(Rmisc)
library(png)
library(grid)
library(cowplot)
library(tidyr)
library(scales)
library(viridis)
theme_set(theme_gray())
path_main="results_main/"
path_chain="results_chain/"
path_time="results_time_series/"
path_rep="results_replicates/"
path_thresh="results_threshold/"
path_beta="results_beta/"
path_KLaB="results_K_L_a_B/"
path_aB="results_a_B/"
path_A="results_A/"
path_FR="results_FR/"
path_CN="results_CN/"

nPP=5
nSpecies=50
nResources=2
dim=nSpecies+nResources

# PLOT OPTIONS ####
# theme
theme<-theme_gray()+
  theme(panel.background = element_blank(),
        text = element_text(size=20, family='Times'),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank())

# colour scales
model_colour<-scale_colour_manual(values=c("chartreuse3","chocolate2","chocolate4"),
                                  name='',
                                  labels=c("SC","NC","C"),
                                  guide = guide_legend(reverse=TRUE,override.aes=list(lwd=c(1.5,1.5,2))))

model_line<-scale_linetype_manual(values=c("solid","longdash", "dashed"),
                                  name='',
                                  labels=c("SC","NC","C"),
                                  guide = guide_legend(reverse=TRUE))

model_colour_short<-scale_colour_manual(values=c("chartreuse3","chocolate4"),
                                        name='',
                                        labels=c("SC","C"),
                                        guide = guide_legend(reverse=TRUE,override.aes=list(lwd=c(1.5,2))))

model_line_short<-scale_linetype_manual(values=c("solid", "dashed"),
                                        name='',
                                        labels=c("SC","C"),
                                        guide = guide_legend(reverse=TRUE))

stab_colour<-scale_fill_manual(values=c("deeppink4","olivedrab3","deepskyblue3"),
                               name="")

stab_colour_short<-scale_colour_manual(values=c("deeppink4","deepskyblue3"),
                                       name="")

chain_colour<-scale_colour_manual(name = "",
                                  values = c("blue","chartreuse3","cadetblue4","red2","darkred"),
                                  labels = c("Mineral\nnutrients","Primary\nproducer","Herbivore","Carnivore","Top predator"),
                                  guide = guide_legend(reverse=TRUE))

chain_fill_recy<-scale_fill_manual(name = "",
                                   values = c("darkred","red2","cadetblue4","chartreuse3","orange4"),
                                   labels = c("Top predator","Carnivore","Herbivore","Primary\nproducer","Indirect\nrecycling"))

chain_fill_biomass<-scale_fill_manual(name = "",
                                      values = c("darkred","red2","cadetblue4","chartreuse3","orange4"),
                                      labels = c("Top predator","Carnivore","Herbivore","Primary\nproducer","Indirect\nrecycling"))

chain_colour_CV<-scale_colour_manual(name = "",
                                     values = c("darkred","red2","cadetblue4","chartreuse3","blue","orange4"),
                                     labels = c("Top predator","Carnivore","Herbivore","Primary\nproducer","Mineral\nnutrients","Recycling"))

chain_line_CV<-scale_linetype_manual(values=c("solid","solid","solid","solid","solid", "dashed"),
                                     name='',
                                     labels = c("Top predator","Carnivore","Herbivore","Primary\nproducer","Mineral\nnutrients","Recycling"),
                                     guide = guide_legend(override.aes=list(lwd=rep(1.5,6))))

# axis limites
lim_short<-xlim(0,300)
lim_chain<-xlim(0,200)

# axis log
x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))

# axis labels
label_nutri<-expression(paste("Mineral nutrient input ",italic("I")))

# d and delta factors
d02<-expression(paste(italic(d),"=0.2"))
d08<-expression(paste(italic(d),"=0.8"))
del02<-expression(paste(italic("\u03B4"),"=0.2"))
del08<-expression(paste(italic("\u03B4"),"=0.8"))
  
# Load data -----------------------------------------------------------------
### data ####
data<-read.table(paste(path_main,"data.txt",sep=""),sep=';',header=T)
data$persistence<-data$NbSpeciesFinal/data$NbSpeciesInit
data$d<-as.factor(data$d)
data$delta<-as.factor(data$delta)
levels(data$d)<-c(d02,d08)
levels(data$delta)<-c(del02,del08)
data$delta = factor(data$delta,levels(data$delta)[c(2,1)])
data$model = factor(data$model,levels(data$model)[c(3,2,1)])

### Species biomass ####
biomass<-read.table(paste(path_main,"biomass.txt",sep=""),sep=';',header=T)
biomass$d<-as.factor(biomass$d)
biomass$delta<-as.factor(biomass$delta)
levels(biomass$d)<-c(d02,d08)
levels(biomass$delta)<-c(del02,del08)
biomass$delta = factor(biomass$delta,levels(biomass$delta)[c(2,1)])
biomass$model = factor(biomass$model,levels(biomass$model)[c(3,2,1)])

### Species biomass CV ####
biomassCV<-read.table(paste(path_main,"biomassCV.txt",sep=""),sep=';',header=T)
biomassCV$d<-as.factor(biomassCV$d)
biomassCV$delta<-as.factor(biomassCV$delta)
levels(biomassCV$d)<-c(d02,d08)
levels(biomassCV$delta)<-c(del02,del08)
biomassCV$delta = factor(biomassCV$delta,levels(biomassCV$delta)[c(2,1)])
biomassCV$model = factor(biomassCV$model,levels(biomassCV$model)[c(3,2,1)])

### Nutrients recycled by species ####
recy<-read.table(paste(path_main,"recy.txt",sep=""),sep=';',header=T)
recy$d<-as.factor(recy$d)
recy$delta<-as.factor(recy$delta)
levels(recy$d)<-c(d02,d08)
levels(recy$delta)<-c(del02,del08)
recy$delta = factor(recy$delta,levels(recy$delta)[c(2,1)])
recy$model = factor(recy$model,levels(recy$model)[c(3,2,1)])

### Nutrients recycled CV ####
recyCV<-read.table(paste(path_main,"recyCV.txt",sep=""),sep=';',header=T)
recyCV$d<-as.factor(recyCV$d)
recyCV$delta<-as.factor(recyCV$delta)
levels(recyCV$d)<-c(d02,d08)
levels(recyCV$delta)<-c(del02,del08)
recyCV$delta = factor(recyCV$delta,levels(recyCV$delta)[c(2,1)])
recyCV$model = factor(recyCV$model,levels(recyCV$model)[c(3,2,1)])

### Detritus produced by species ####
detritus<-read.table(paste(path_main,"detritus.txt",sep=""),sep=';',header=T)
detritus$d<-as.factor(detritus$d)
detritus$delta<-as.factor(detritus$delta)
levels(detritus$d)<-c(d02,d08)
levels(detritus$delta)<-c(del02,del08)
detritus$delta = factor(detritus$delta,levels(detritus$delta)[c(2,1)])
detritus$model = factor(detritus$model,levels(detritus$model)[c(3,2,1)])

### Detritus produced by species CV ####
detritusCV<-read.table(paste(path_main,"detritusCV.txt",sep=""),sep=';',header=T)
detritusCV$d<-as.factor(detritusCV$d)
detritusCV$delta<-as.factor(detritusCV$delta)
levels(detritusCV$d)<-c(d02,d08)
levels(detritusCV$delta)<-c(del02,del08)
detritusCV$delta = factor(detritusCV$delta,levels(detritusCV$delta)[c(2,1)])
detritusCV$model = factor(detritusCV$model,levels(detritusCV$model)[c(3,2,1)])

### Trophic level ####
TL<-read.table(paste(path_main,"TL.txt",sep=""),sep=';',header=T)
TL$d<-as.factor(TL$d)
TL$delta<-as.factor(TL$delta)
levels(TL$d)<-c(d02,d08)
levels(TL$delta)<-c(del02,del08)
TL$delta = factor(TL$delta,levels(TL$delta)[c(2,1)])
TL$model = factor(TL$model,levels(TL$model)[c(3,2,1)])

### Extinction time ####
tExt<-read.table(paste(path_main,"tExtinction.txt",sep=""),sep=';',header=T)
tExt$d<-as.factor(tExt$d)
tExt$delta<-as.factor(tExt$delta)
levels(tExt$d)<-c(d02,d08)
levels(tExt$delta)<-c(del02,del08)
tExt$delta = factor(tExt$delta,levels(tExt$delta)[c(2,1)])
tExt$model = factor(tExt$model,levels(tExt$model)[c(3,2,1)])

# FIGURES FINALES - MAIN TEXT -----------------------------------------------------------------
### Figure 1 ####
# Recycled nutrients ####
recyPP <- summarySE(data[data$persistence>0 & data$model=="C",], measurevar="RecyPP", groupvars=c("I","d","delta"),na.rm=TRUE)
recyPP$recy <- recyPP$RecyPP
recyPP<-recyPP[,c(-4,-5)]
recyPP$type="Primary producers"
recySP <- summarySE(data[data$persistence>0 & data$model=="C",], measurevar="RecySP", groupvars=c("I","d","delta"),na.rm=TRUE)
recySP$recy <- recySP$RecySP
recySP<-recySP[,c(-4,-5)]
recySP$type="Consumers"
databis<-rbind(recyPP,recySP)
rm(recyPP,recySP)
RecyInd <- summarySE(data[data$persistence>0 & data$model=="C",], measurevar="RecyInd", groupvars=c("I","d","delta"),na.rm=TRUE)
RecyInd$recy <- RecyInd$RecyInd
RecyInd<-RecyInd[,c(-4,-5)]
RecyInd$type="Indirect recycling"
databis<-rbind(databis,RecyInd)
rm(RecyInd)
databis$type<-as.factor(databis$type)
databis$type = factor(databis$type,levels(databis$type)[c(1,3,2)])

p1<-ggplot(data=databis)+
      geom_area(aes(I,recy,fill=type))+
      geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("firebrick3","chartreuse3","orange4"),
                        name="")+
      theme+theme(legend.position = "bottom")+
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
      lim_short+
      xlab(label_nutri)+
      ylab("Nutrient recycled")

# Persistence ####
databis <- summarySE(data, measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)

p2<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=persistence-ci, ymax=persistence+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5, "cm"),legend.position = "bottom")+
      lim_short+
      xlab(label_nutri)+
      ylab("Species persistence")

# Final graph ####
graph<-plot_grid(p1,p2,
                 labels=c("A","B"),label_size = 20,
                 nrow=1, ncol = 2, align = "vh")
ggsave("main_figure_1.pdf",graph, width = 12, height = 7, device = cairo_pdf)

### Figure 2 ####
# average CV species ####
# biomass and total biomass
nparams=5
databis<-biomass[,-which(names(biomass)%in%c("seed","N","D"))]
databis$biomass_tot<-rowSums(databis[,(nparams+1):dim(databis)[2]])
databis[,(nparams+1):(dim(databis)[2]-1)]<-databis[,(nparams+1):(dim(databis)[2]-1)]/databis$biomass_tot
databis$biomass_tot<-NULL
# biomass CV
databis_2<-biomassCV[,-which(names(biomass)%in%c("seed","N","D"))]
databis_2$model<-droplevels(databis_2$model)
# CV weighted by the average biomass
for(i in c((nparams+1):dim(databis)[2])){
  databis[,i]<-databis[,i]*databis_2[,i]
}
databis[is.na(databis)]=0
databis$CV<-rowSums(databis[,(nparams+1):dim(databis)[2]])
# persistence
databis_2<-data
databis$persistence<-databis_2$persistence
rm(databis_2)
# Average CV 
databis<-databis[databis$persistence>0,]
databis<-summarySE(databis, measurevar="CV", groupvars=c("I","d",'delta','model'),na.rm=TRUE)

p1<-ggplot(data=databis)+
      geom_line(aes(I,CV,colour=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=CV-ci, ymax=CV+ci,colour=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5,"cm"),
                  legend.position="bottom")+
      lim_short+
      xlab(label_nutri)+
      ylab("Average weighted biomass CV")

# average CV recycling ####
databis<-data[data$model=="C" ,c("simu_ID","I","delta","d","model","NbSpeciesFinal","NbPPFinal","IrecyCV","RecyPPcv","RecySPcv","RecyIndCV")]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model","NbSpeciesFinal","NbPPFinal"),
              variable.name = "recy",
              value.name = "CV")
databis<-databis[(databis$recy=="IrecyCV" & databis$NbSpeciesFinal>0) |
                   (databis$recy=="RecyPPcv" & databis$NbPPFinal>0) |
                   (databis$recy=="RecySPcv" & (databis$NbSpeciesFinal-databis$NbPPFinal)>0) |
                   (databis$recy=="RecyIndCV" & databis$NbSpeciesFinal>0),]
databis<-summarySE(databis, measurevar="CV", groupvars=c("I","d",'delta','model','recy'),na.rm=TRUE)
levels(databis$recy)<-c("Total recycling","Primary producers","Consumers","Indirect recycling")
databis$recy = factor(databis$recy,levels(databis$recy)[c(1,4,2,3)])

p2<-ggplot(data=databis)+
      geom_line(aes(I,CV,colour=recy,linetype=recy), size=1.5)+
      geom_errorbar(aes(I,ymin=CV-ci, ymax=CV+ci,colour=recy), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_colour_manual(values=c("blue","chocolate4","chartreuse3","red"),
                          name='')+
      scale_linetype_manual(values=c("solid","longdash", "dashed","dotted"),
                            name='')+
      theme+theme(legend.key.width = unit(1.5,"cm"),
                  legend.position="bottom")+
      guides(col = guide_legend(nrow = 2))+
      lim_short+
      xlab(label_nutri)+
      ylab("Average recycling CV")

# Final graph ####
graph<-plot_grid(p1,p2,
                 labels=c("A","B"),label_size = 20,
                 nrow=1, ncol = 2, align = "vh")
ggsave("main_figure_2.pdf",graph, width = 12, height = 7, device = cairo_pdf)

### Figure 3 - stabilising effect ####
# Primary producers
# biomass
databis<-biomass[biomass$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis<-databis[,-c(which(names(databis)=="x6"):which(names(databis)=="x50"))]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis<-dcast(databis,simu_ID+I+d+delta+species~model,value.var='biomass')
colnames(databis)<-c("simu_ID","I","d","delta","species","biomass_SC","biomass_C")
# CV
databis_2<-biomassCV[biomassCV$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis_2<-databis_2[,-c(which(names(databis_2)=="x6"):which(names(databis_2)=="x50"))]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "CV")
databis_2<-dcast(databis_2,simu_ID+I+d+delta+species~model,value.var='CV')
colnames(databis_2)<-c("simu_ID","I","d","delta","species","SC","C")
# stabilising ?
databis$SC<-databis_2$SC
databis$C<-databis_2$C
rm(databis_2)
databis<-databis[databis$biomass_C>0 & databis$biomass_SC>0,]
databis$effect<-"neutral"
treshold=1e-4
databis$effect[(databis$SC-databis$C)>treshold | (databis$biomass_SC==0 & databis$biomass_C>0)]="stabilising"
databis$effect[(databis$SC-databis$C)< -treshold | (databis$biomass_SC>0 & databis$biomass_C==0)]="destabilising"
# aggregation
databis<-summarySE(databis, measurevar="C", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)
databis_2<-aggregate(databis$N,databis[,c(1:3)],FUN=sum)
names(databis_2)<-c(names(databis_2)[1:(dim(databis_2)[2]-1)],"mean")
databis<-merge(databis,databis_2,by=names(databis_2)[1:(dim(databis_2)[2]-1)])
rm(databis_2)
databis$frac<-databis$N/databis$mean

zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  effect=levels(as.factor(databis$effect)))
databis<-merge(databis,zero,by=c("I","d","delta","effect"),all = TRUE)
databis[is.na(databis)] <- 0

p1<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      geom_hline(yintercept=0.5, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("Fraction of species")+
      ggtitle("Primary producers")

# Consumers
databis<-biomass[biomass$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis<-databis[,-c(which(names(databis)=="x1"):which(names(databis)=="x5"))]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis<-dcast(databis,simu_ID+I+d+delta+species~model,value.var='biomass')
colnames(databis)<-c("simu_ID","I","d","delta","species","biomass_SC","biomass_C")
# CV
databis_2<-biomassCV[biomassCV$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis_2<-databis_2[,-c(which(names(databis_2)=="x1"):which(names(databis_2)=="x5"))]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "CV")
databis_2<-dcast(databis_2,simu_ID+I+d+delta+species~model,value.var='CV')
colnames(databis_2)<-c("simu_ID","I","d","delta","species","SC","C")
# stabilising ?
databis$SC<-databis_2$SC
databis$C<-databis_2$C
rm(databis_2)
databis<-databis[databis$biomass_C>0 & databis$biomass_SC>0,]
databis$effect<-"neutral"
treshold=1e-4
databis$effect[(databis$SC-databis$C)>treshold | (databis$biomass_SC==0 & databis$biomass_C>0)]="stabilising"
databis$effect[(databis$SC-databis$C)< -treshold | (databis$biomass_SC>0 & databis$biomass_C==0)]="destabilising"
# aggregation
databis<-summarySE(databis, measurevar="C", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)
databis_2<-aggregate(databis$N,databis[,c(1:3)],FUN=sum)
names(databis_2)<-c(names(databis_2)[1:(dim(databis_2)[2]-1)],"mean")
databis<-merge(databis,databis_2,by=names(databis_2)[1:(dim(databis_2)[2]-1)])
rm(databis_2)
databis$frac<-databis$N/databis$mean

zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  effect=levels(as.factor(databis$effect)))
databis<-merge(databis,zero,by=c("I","d","delta","effect"),all = TRUE)
databis[is.na(databis)] <- 0

p2<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      stab_colour
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      geom_hline(yintercept=0.5, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("Fraction of species")+
      ggtitle("Consumers")
      
graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.07, 0, 0.05, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 25)

ggsave("main_figure_3.pdf",graph, width = 12, height = 6, device = cairo_pdf) # 12 10

# FIGURES FINALES - APPENDIX -----------------------------------------------------------------
### Time series ####
colfunc<-colorRampPalette(c("chartreuse3","cadetblue4","red2","red3","red4","darkred"))
tmin=0
tmax=50
x_axis_TS<-scale_x_continuous(breaks=seq(tmin, tmax, 25))

TS_data<-read.table(paste(path_time,"time_series_data_0.txt",sep=""),sep=';',header=T)
TS_biomass<-read.table(paste(path_time,"time_series_biomass_0.txt",sep=""),sep=';',header=T)
TS_recy<-read.table(paste(path_time,"time_series_recy_0.txt",sep=""),sep=';',header=T)
TS_TL<-read.table(paste(path_time,"TL_0.txt",sep=""),sep=';',header=T)

TS_data$t<-TS_data$t-1000
TS_data$d<-as.factor(TS_data$d)
TS_data$delta<-as.factor(TS_data$delta)
levels(TS_data$d)<-c(d02,d08)
levels(TS_data$delta)<-c(del02,del08)
TS_data$delta = factor(TS_data$delta,levels(TS_data$delta)[c(2,1)])

TS_biomass$t<-TS_biomass$t-1000
TS_biomass$d<-as.factor(TS_biomass$d)
TS_biomass$delta<-as.factor(TS_biomass$delta)
levels(TS_biomass$d)<-c(d02,d08)
levels(TS_biomass$delta)<-c(del02,del08)
TS_biomass$delta = factor(TS_biomass$delta,levels(TS_biomass$delta)[c(2,1)])

TS_TL$t<-TS_TL$t-1000
TS_TL$d<-as.factor(TS_TL$d)
TS_TL$delta<-as.factor(TS_TL$delta)
levels(TS_TL$d)<-c(d02,d08)
levels(TS_TL$delta)<-c(del02,del08)
TS_TL$delta = factor(TS_TL$delta,levels(TS_TL$delta)[c(2,1)])

TS_recy$t<-TS_recy$t-1000
TS_recy$d<-as.factor(TS_recy$d)
TS_recy$delta<-as.factor(TS_recy$delta)
levels(TS_recy$d)<-c(d02,d08)
levels(TS_recy$delta)<-c(del02,del08)
TS_recy$delta = factor(TS_recy$delta,levels(TS_recy$delta)[c(2,1)])

# Species biomass ####
databis<-TS_biomass[TS_biomass$model=="C" & TS_biomass$t>=tmin & TS_biomass$t<=tmax,-c(2,8,9)]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model","t"),
              variable.name = "species",
              value.name = "biomass")

p1<-ggplot(data=databis)+
      geom_line(aes(t,biomass,colour=species),size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_colour_manual(values = colfunc(50))+
      theme+theme(legend.position = "none")+
      x_axis_TS+
      xlab("Time")+
      ylab("Species biomass")

# Nutrient and detritus stocks ####
databis<-TS_biomass[TS_biomass$model=="C" & TS_biomass$t>=tmin & TS_biomass$t<=tmax,c(3:9)]
databis<-melt(databis, id.vars = c("I","delta","d","model","t"),
              variable.name = "species",
              value.name = "biomass")

p2<-ggplot()+
      geom_line(data=databis,aes(t,biomass,colour=species),size=1)+
      scale_colour_manual(values = c("blue","chocolate4"),
                          name='',
                          labels=c("Mineral nutrients","Detritus"))+
      theme
legend_1<-get_legend(p2)

p2<-ggplot()+
      geom_line(data=databis,aes(t,biomass,colour=species),size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_colour_manual(values = c("blue","chocolate4"),
                          name='',
                          labels=c("Mineral nutrients","Detritus"))+
      theme+theme(legend.position = "none")+
      x_axis_TS+
      xlab("Time")+
      ylab("Nutrient stock")

# # Recycling (PP and SP)
# databis<-TS_data[TS_data$model=="C" & TS_data$t>=tmin & TS_data$t<=tmax,which(names(TS_data)%in%c("I","delta","d","model","t","RecyPP","RecySP","RecyInd"))]
# databis<-melt(databis, id.vars = c("I","delta","d","model","t"),
#               variable.name = "type",
#               value.name = "recy")
# levels(databis$type)<-c("Primary producers","Consumers","Indirect recycling")
# databis$type = factor(databis$type,levels(databis$type)[c(2,1,3)])
# 
# p3<-ggplot()+
#       geom_area(data=databis,aes(t,recy,fill=type))+
#       facet_grid(delta~d, labeller=label_parsed)+
#       scale_fill_manual(values = c("red","chartreuse3","chocolate4"),
#                         name = "")+
#       theme+
#       xlab("Time")+
#       ylab("Quantity of recycled nutrients")
# 
# graph<-plot_grid(p1, p2, p3,
#                  labels = c("A","B","C"), label_size = 25,
#                  ncol = 1, nrow=3, align = "v")
# ggsave("time_series_1.pdf",graph, width = 12, height = 30)

# Biomass per trophic level ####
databis<-TS_biomass[TS_biomass$model=="C" & TS_biomass$t>=tmin & TS_biomass$t<=tmax,-c(2,8,9)]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model","t"),
              variable.name = "species",
              value.name = "biomass")
databis_2<-TS_TL[TS_TL$model=="C" ,-c(2)]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "TL")
databis<-merge(databis,databis_2,by=c("simu_ID","I","delta","d","model","species"))
rm(databis_2)
databis<-databis[databis$TL>0,]
databis$TL[databis$TL>2 & databis$TL<=3]=3
databis$TL[databis$TL>3 & databis$TL<=4]=4
databis$TL[databis$TL>4 & databis$TL<=5]=5
databis$TL<-as.factor(databis$TL)
databis$TL = factor(databis$TL,levels(databis$TL)[seq(length(levels(databis$TL)),1,-1)])
databis$species<-NULL
databis<-aggregate(databis[,7],databis[,c(-7)],FUN=sum)

p3<-ggplot(data=databis)+
      geom_area(aes(t,x,fill=TL))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("darkred","red3","cadetblue4","chartreuse3"),
                        name='Trophic level')+
      theme+theme(legend.position = "none")+
      x_axis_TS+
      xlab("Time")+
      ylab("Average biomass")

# Recycling per trophic level ####
databis<-TS_recy[TS_recy$model=="C" & TS_recy$t>=tmin & TS_recy$t<=tmax,-c(2)]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model","t"),
              variable.name = "species",
              value.name = "recy")
databis_2<-TS_TL[TS_TL$model=="C" ,-c(2)]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "TL")
databis<-merge(databis,databis_2,by=c("simu_ID","I","delta","d","model","species"))
databis<-databis[databis$TL>0,]
databis$TL[databis$TL>2 & databis$TL<=3]=3
databis$TL[databis$TL>3 & databis$TL<=4]=4
databis$TL[databis$TL>4 & databis$TL<=5]=5
databis$TL<-as.factor(databis$TL)
databis$TL = factor(databis$TL,levels(databis$TL)[seq(length(levels(databis$TL)),1,-1)])
databis$species<-NULL
databis<-aggregate(databis[,7],databis[,c(-7)],FUN=sum)
databis_2<-TS_data[TS_data$model=="C" & TS_data$t>=tmin & TS_data$t<=tmax, which(names(TS_data)%in%c("simu_ID","I","delta","d","model","t","RecyInd"))]
names(databis_2)[dim(databis_2)[2]]="x"
databis_2$TL<-"Indirect recycling"
databis<-rbind(databis,databis_2)
rm(databis_2)

p4<-ggplot(data=databis)+
      geom_area(aes(t,x,fill=TL))+
      scale_fill_manual(values=c("darkred","red3","cadetblue4","chartreuse3","chocolate4"),
                        name='Trophic level')+
      theme
legend_2<-get_legend(p4)

p4<-ggplot(data=databis)+
      geom_area(aes(t,x,fill=TL))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("darkred","red3","cadetblue4","chartreuse3","chocolate4"),
                        name='Trophic level')+
      theme+theme(legend.position = "none")+
      x_axis_TS+
      xlab("Time")+
      ylab("Average quantity of recycled nutrients")

# Final graph ####
graph<-ggdraw(xlim = c(0, 2.35), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend_1, 2.13, 1, 0.05, 1)+
  draw_plot(legend_2, 2.13, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 25)
ggsave("supp_time_series.pdf",graph, width = 15, height = 12, device = cairo_pdf)

### Complementary results on species dynamics ####
# Nutrient CV ####
databis<-biomassCV[biomassCV$model!="NC",c(1,3:7)]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "CV")
databis<-summarySE(databis, measurevar="CV", groupvars=c("I","d",'delta','model'),na.rm=TRUE)

p1<-ggplot(data=databis)+
      geom_line(aes(I,CV,colour=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=CV-ci, ymax=CV+ci,colour=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour_short+
      model_line_short+
      theme+theme(legend.key.width = unit(1.5,"cm"),
                  legend.position="non")+
      lim_short+
      xlab(label_nutri)+
      ylab("Average mineral nutrient stock CV")

# TLmax ####
databis <- summarySE(data[data$persistence>0,], measurevar="TLmax", groupvars=c("I","d",'delta','model'),na.rm=TRUE)

p2<-ggplot(data=databis)+
      geom_line(aes(I,TLmax,color=model,linetype=model), size=2)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5, "cm"))
legend_1<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(I,TLmax,color=model,linetype=model), size=2)+
      geom_errorbar(aes(I,ymin=TLmax-ci, ymax=TLmax+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Maximum trophic level")

# TL species ####
databis<-TL[TL$model=="C",-c(2)]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "TL")
databis<-databis[databis$TL>0,-c(which(names(databis)=="species"))]
databis<-unique(databis)

p3<-ggplot(data=databis)+
      geom_point(aes(I,TL,colour=TL),alpha=1, size=2)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_color_gradientn(colours=c("chartreuse3","cadetblue4","red3","darkred","black"))+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Trophic level")

# CV and TLmax comparison ####
nparams=5
databis<-biomass[biomass$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis$biomass_tot<-rowSums(databis[,(nparams+1):dim(databis)[2]])
databis[,(nparams+1):(dim(databis)[2]-1)]<-databis[,(nparams+1):(dim(databis)[2]-1)]/databis$biomass_tot
databis$biomass_tot<-NULL
# biomass CV
databis_2<-biomassCV[biomassCV$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis_2$model<-droplevels(databis_2$model)
# CV weighted by the average biomass
for(i in c((nparams+1):dim(databis)[2])){
  databis[,i]<-databis[,i]*databis_2[,i]
}
databis[is.na(databis)]=0
databis$CV<-rowSums(databis[,(nparams+1):dim(databis)[2]])
# persistence
databis_2<-data[data$model!="NC",]
databis$persistence<-databis_2$persistence
databis$TLmax<-databis_2$TLmax
rm(databis_2)
# Average CV 
databis<-databis[databis$persistence>0,]

p4<-ggplot(data=databis)+
  geom_point(aes(I,CV,colour=TLmax),alpha=1, size=1.5)+
  scale_color_gradientn(colours=c("chartreuse3","cadetblue4","red3","darkred","black"),
                        name="Maximum\ntrophic\nlevel")+
  theme
legend_2<-get_legend(p4)

p4<-ggplot(data=databis)+
      geom_point(aes(I,CV,colour=TLmax),alpha=1, size=1.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_color_gradientn(colours=c("chartreuse3","cadetblue4","red3","darkred","black"),
                            name="Maximum\ntrophic\nlevel")+
      theme+theme(legend.position = "none")+
      lim_short+
      ylim(0,2)+
      xlab(label_nutri)+
      ylab("Average weighted biomass CV")

# Final graph ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend_1, 2.1, 1, 0.05, 1)+
  draw_plot(legend_2, 2.1, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 25)
ggsave("supp_TL_CV.png",graph, width = 14, height = 12)

### Flows and stocks ####
# biomass repartition ####
databis<-biomass[biomass$model=="C",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis$type="SP"
databis$type[which(databis$species%in%c("x1","x2","x3","x4","x5"))] = "PP"
databis<-aggregate(databis$biomass,databis[,c(1:4,dim(databis)[2])],FUN=sum)
names(databis)<-c(names(databis)[1:(dim(databis)[2]-1)],"biomass")
databis <- summarySE(databis, measurevar="biomass", groupvars=c("I","d",'delta','type'),na.rm=TRUE)

p1<-ggplot(data=databis)+
    geom_area(aes(I,biomass,fill=type))+
    scale_fill_manual(values=c("ForestGreen","red"),
                      labels=c("Primary producers","Consumers"),
                      name="")+
    theme
legend_1<-get_legend(p1)

p1<-ggplot(data=databis)+
      geom_area(aes(I,biomass,fill=type))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("ForestGreen","red"),
                        labels=c("Primary producers","Consumers"),
                        name="")+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Average total biomass")

# biomass production ####
databis<-data[data$model=="C",which(names(data)%in%c("simu_ID","I","delta","d","model","PPprod","SPprod","NbSpeciesFinal","NbPPFinal"))]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model","NbSpeciesFinal","NbPPFinal"),
              variable.name = "type",
              value.name = "production")
databis<-databis[(databis$type=="PPprod" & databis$NbPPFinal>0)
                 | (databis$type=="SPprod" & (databis$NbSpeciesFinal-databis$NbPPFinal)>0),]
databis<-summarySE(databis, measurevar="production", groupvars=c("I","d",'delta','type'),na.rm=TRUE)

p2<-ggplot(data=databis)+
      geom_area(aes(I,production,fill=type))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("ForestGreen","red"),
                        labels=c("Primary production","Secondary production"),
                        name="")+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Average biomass production")

# nutrient and detritus stocks ####
databis<-biomass[biomass$model=="C" & data$persistence>0,which(names(biomass)%in%c("simu_ID","I","delta","d","model","N","D"))]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "compartment",
              value.name = "biomass")
databis <- summarySE(databis, measurevar="biomass", groupvars=c("I","d",'delta','compartment'),na.rm=TRUE)

p3<-ggplot(data=databis)+
      geom_area(aes(I,biomass,fill=compartment))+
      scale_fill_manual(values=c("blue","chocolate4"),
                        labels=c("Mineral nutrients","Detritus"),
                        name="")+
      theme
legend_2<-get_legend(p3)

p3<-ggplot(data=databis)+
      geom_area(aes(I,biomass,fill=compartment))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("blue","chocolate4"),
                        labels=c("Mineral nutrients","Detritus"),
                        name="")+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Abiotic nutrient stock")

# density dependance mortality rate ####
databis<-data[data$model=="C",which(names(data)%in%c("simu_ID","I","delta","d","model","PPmort","SPmort","NbSpeciesFinal","NbPPFinal"))]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model","NbSpeciesFinal","NbPPFinal"),
              variable.name = "type",
              value.name = "mort")
databis<-databis[(databis$type=="PPmort" & databis$NbPPFinal>0)
                 | (databis$type=="SPmort" & (databis$NbSpeciesFinal-databis$NbPPFinal)>0),]
databis<-summarySE(databis, measurevar="mort", groupvars=c("I","d",'delta','type'),na.rm=TRUE)

p4<-ggplot(data=databis)+
      geom_area(aes(I,mort,fill=type))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("ForestGreen","red"),
                        labels=c("Primary producers","Consumers"),
                        name="")+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Average density dependant mortality")

# Final graph ####
graph<-ggdraw(xlim = c(0, 2.35), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p4, 0, 0, 1, 1)+
  draw_plot(p3, 1, 0, 1, 1)+
  draw_plot(legend_1, 2.15, 1, 0.05, 1)+
  draw_plot(legend_2, 2.15, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 25)
ggsave("supp_biomass_prod.pdf",graph, width = 14, height = 12, device = cairo_pdf)
### Recycling and detritus production per trophic level ####
# Zoom reycled nutrients ####
recyPP <- summarySE(data[data$persistence>0 & data$model=="C",], measurevar="RecyPP", groupvars=c("I","d","delta"),na.rm=TRUE)
recyPP$recy <- recyPP$RecyPP
recyPP<-recyPP[,c(-4,-5)]
recyPP$type="Primary producers"
recySP <- summarySE(data[data$persistence>0 & data$model=="C",], measurevar="RecySP", groupvars=c("I","d","delta"),na.rm=TRUE)
recySP$recy <- recySP$RecySP
recySP<-recySP[,c(-4,-5)]
recySP$type="Consumers"
databis<-rbind(recyPP,recySP)
rm(recyPP,recySP)
RecyInd <- summarySE(data[data$persistence>0 & data$model=="C",], measurevar="RecyInd", groupvars=c("I","d","delta"),na.rm=TRUE)
RecyInd$recy <- RecyInd$RecyInd
RecyInd<-RecyInd[,c(-4,-5)]
RecyInd$type="Indirect recycling"
databis<-rbind(databis,RecyInd)
rm(RecyInd)
databis$type<-as.factor(databis$type)
databis$type = factor(databis$type,levels(databis$type)[c(1,3,2)])

p1<-ggplot(data=databis)+
      geom_area(aes(I,recy,fill=type))+
      scale_fill_manual(values=c("firebrick3","chartreuse3","orange4"),
                        name="")+
      theme+theme(legend.position = "bottom")
legend_1<-get_legend(p1)

p1<-ggplot(data=databis)+
      geom_area(aes(I,recy,fill=type))+
      geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("firebrick3","chartreuse3","orange4"),
                        name="")+
      theme+theme(legend.position = "none")+
      xlim(0,100)+
      xlab(label_nutri)+
      ylab("Nutrient recycled")

# Detritus production per trophic level ####
# detritus
databis<-detritus[detritus$model=="C",-which(names(biomass)%in%c("seed"))]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "detritus")
# TL
databis_2<-TL[TL$model=="C",-which(names(biomass)%in%c("seed"))]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "TL")
# TL categories
databis$TL<-databis_2$TL
rm(databis_2)
databis$model<-NULL
databis$species<-NULL
databis<-databis[databis$TL>0,]
databis$bidon<-databis$TL
databis$TL[databis$bidon==1]="TL1"
databis$TL[databis$bidon>=2 & databis$bidon<3]="TL2+"
databis$TL[databis$bidon>=3 & databis$bidon<4]="TL3+"
databis$TL[databis$bidon>=4]="TL4+"
databis$TL<-as.factor(databis$TL)
databis$bidon<-NULL
# sum
databis<-aggregate(databis$detritus,databis[,-5],FUN=sum)
zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  TL=levels(databis$TL))
databis<-merge(databis,zero,by=c("I","d","delta","TL"),all = TRUE)
databis[is.na(databis)] <- 0
rm(zero)
databis<-summarySE(databis, measurevar="x", groupvars=c("I","d","delta","TL"),na.rm=TRUE)
# fraction of produced detritus
databis_2<-aggregate(databis$x,databis[,c(1:3)],FUN=sum)
names(databis_2)[dim(databis_2)[2]]="frac"
databis<-merge(databis,databis_2,by=c("I","d","delta"),all = TRUE)
rm(databis_2)
databis$frac<-databis$x/databis$frac
databis$TL = factor(databis$TL,levels(databis$TL)[c(4,3,2,1)])

p2<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=TL))+
      scale_fill_manual(values=c("darkred","red3","cadetblue4","chartreuse3"),
                        name="Aggregated\ntrophic levels")+
      theme
legend_2<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=TL))+
      facet_grid(delta~d, labeller=label_parsed)+
      scale_fill_manual(values=c("darkred","red3","cadetblue4","chartreuse3"),
                        name="Aggregated\ntrophic levels")+
      theme+theme(legend.position = "none")+
      xlab(label_nutri)+
      ylab("Fraction of produced detritus")

# Final graph ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 1.1)) +
  draw_plot(p1, 0, 0.1, 1, 1)+
  draw_plot(p2, 1, 0.1, 1, 1)+
  draw_plot(legend_1, 0, 0, 1, 0.1)+
  draw_plot(legend_2, 2.1, 0.1, 0.05, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1.1,1.1), size = 25)
ggsave("supp_detritus.pdf",graph, width = 14, height = 6, device = cairo_pdf)

### stabilising effect comparison data ####
# biomass
databis<-biomass[biomass$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis<-dcast(databis,simu_ID+I+d+delta+species~model,value.var='biomass')
colnames(databis)<-c("simu_ID","I","d","delta","species","biomass_SC","biomass_C")
# TL
databis_2<-TL[TL$model!="NC",-c(2)]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "TL")
databis_2<-dcast(databis_2,simu_ID+I+d+delta+species~model,value.var='TL')
colnames(databis_2)<-c("simu_ID","I","d","delta","species","SC","C")
databis_2$SC<-(databis_2$SC+databis_2$C)/2
databis_2$C<-NULL
names(databis_2)[dim(databis_2)[2]]="TL"
databis$TL<-databis_2$TL
# CV
databis_2<-biomassCV[biomassCV$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "CV")
databis_2<-dcast(databis_2,simu_ID+I+d+delta+species~model,value.var='CV')
colnames(databis_2)<-c("simu_ID","I","d","delta","species","SC","C")
databis$SC<-databis_2$SC
databis$C<-databis_2$C
rm(databis_2)
databis<-databis[databis$biomass_C>0 & databis$biomass_SC>0,]
databis$effect<-"neutral"
treshold=1e-4
databis$effect[databis$SC-databis$C>treshold | (databis$biomass_SC==0 & databis$biomass_C>0)]="stabilising"
databis$effect[databis$SC-databis$C< -treshold | (databis$biomass_SC>0 & databis$biomass_C==0)]="destabilising"
# General data
stab<-databis

# P1
databis<-stab[stab$TL==1,]
databis<-summarySE(databis, measurevar="C", groupvars=c("I","d",'delta','TL','effect'),na.rm=TRUE)
databis_2<-aggregate(databis$N,databis[,c(1:3)],FUN=sum)
names(databis_2)<-c(names(databis_2)[1:(dim(databis_2)[2]-1)],"mean")
databis<-merge(databis,databis_2,by=names(databis_2)[1:(dim(databis_2)[2]-1)])
rm(databis_2)
databis$frac<-databis$N/databis$mean

zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  effect=levels(as.factor(databis$effect)))
databis<-merge(databis,zero,by=c("I","d","delta","effect"),all = TRUE)
databis[is.na(databis)] <- 0

p1<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      geom_hline(yintercept=0.5, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("Fraction of species")+
      ggtitle("TL=1")

# P2
databis<-stab[stab$TL>=2 & stab$TL<3,]
databis<-summarySE(databis, measurevar="C", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)
databis_2<-aggregate(databis$N,databis[,c(1:3)],FUN=sum)
names(databis_2)<-c(names(databis_2)[1:(dim(databis_2)[2]-1)],"mean")
databis<-merge(databis,databis_2,by=names(databis_2)[1:(dim(databis_2)[2]-1)])
databis$frac<-databis$N/databis$mean

zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  effect=levels(as.factor(databis$effect)))
databis<-merge(databis,zero,by=c("I","d","delta","effect"),all = TRUE)
databis[is.na(databis)] <- 0

p2<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      geom_hline(yintercept=0.5, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("Fraction of species")+
      ggtitle("2⩽TL<3")

# P3
databis<-stab[stab$TL>=3 & stab$TL<4,]
databis<-summarySE(databis, measurevar="C", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)
databis_2<-aggregate(databis$N,databis[,c(1:3)],FUN=sum)
names(databis_2)<-c(names(databis_2)[1:(dim(databis_2)[2]-1)],"mean")
databis<-merge(databis,databis_2,by=names(databis_2)[1:(dim(databis_2)[2]-1)])
databis$frac<-databis$N/databis$mean

zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  effect=levels(as.factor(databis$effect)))
databis<-merge(databis,zero,by=c("I","d","delta","effect"),all = TRUE)
databis[is.na(databis)] <- 0

p3<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      geom_hline(yintercept=0.5, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("Fraction of species")+
      ggtitle("3⩽TL<4")

# P4
databis<-stab[stab$TL>=4,]
databis<-summarySE(databis, measurevar="C", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)
databis_2<-aggregate(databis$N,databis[,c(1:3)],FUN=sum)
names(databis_2)<-c(names(databis_2)[1:(dim(databis_2)[2]-1)],"mean")
databis<-merge(databis,databis_2,by=names(databis_2)[1:(dim(databis_2)[2]-1)])
rm(databis_2)
databis$frac<-databis$N/databis$mean

zero<-expand.grid(I=unique(databis$I),
                  delta=levels(databis$delta),
                  d=levels(databis$d),
                  effect=levels(as.factor(databis$effect)))
databis<-merge(databis,zero,by=c("I","d","delta","effect"),all = TRUE)
databis[is.na(databis)] <- 0

p4<-ggplot(data=databis)+
      geom_area(aes(I,frac,fill=effect))+
      stab_colour+
      theme
legend<-get_legend(p4)

p4<-ggplot(data=databis)+
  geom_area(aes(I,frac,fill=effect))+
  geom_hline(yintercept=0.5, linetype="dashed", size=1)+
  facet_grid(delta~d, labeller=label_parsed)+
  stab_colour+
  theme+theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))+
  xlab(label_nutri)+
  ylab("Fraction of species")+
  ggtitle("4⩽TL")

graph<-ggdraw(xlim = c(0, 2.35), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend, 2.13, 0.5, 0.05, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 25)
ggsave("supp_stab_complex.pdf",graph, width = 14, height = 12, device = cairo_pdf)

### stabilising effect quantitative ####
# Primary producers
# biomass
databis<-biomass[biomass$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis<-databis[,-c(which(names(databis)=="x6"):which(names(databis)=="x50"))]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis<-dcast(databis,simu_ID+I+d+delta+species~model,value.var='biomass')
colnames(databis)<-c("simu_ID","I","d","delta","species","biomass_SC","biomass_C")
# CV
databis_2<-biomassCV[biomassCV$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis_2<-databis_2[,-c(which(names(databis_2)=="x6"):which(names(databis_2)=="x50"))]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "CV")
databis_2<-dcast(databis_2,simu_ID+I+d+delta+species~model,value.var='CV')
colnames(databis_2)<-c("simu_ID","I","d","delta","species","SC","C")
# stabilising ?
databis$SC<-databis_2$SC
databis$C<-databis_2$C
rm(databis_2)
databis<-databis[databis$biomass_C>0 & databis$biomass_SC>0,]
databis$diff<-databis$SC-databis$C
databis$effect<-"neutral"
treshold=1e-4
databis$effect[(databis$SC-databis$C)>treshold | (databis$biomass_SC==0 & databis$biomass_C>0)]="stabilising"
databis$effect[(databis$SC-databis$C)< -treshold | (databis$biomass_SC>0 & databis$biomass_C==0)]="destabilising"
databis<-databis[databis$effect!="neutral",]
# aggregation
databis<-summarySE(databis, measurevar="diff", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)
p1<-ggplot(data=databis)+
      geom_line(aes(I,diff,colour=effect))+
      geom_errorbar(aes(I,ymin=diff-ci, ymax=diff+ci,colour=effect), width=10, show.legend=FALSE)+
      geom_hline(yintercept=0, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour_short+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("CV difference between SC and C models")+
      ggtitle("Primary producers")

# Consumers
databis<-biomass[biomass$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis$model<-droplevels(databis$model)
databis<-databis[,-c(which(names(databis)=="x1"):which(names(databis)=="x5"))]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis<-dcast(databis,simu_ID+I+d+delta+species~model,value.var='biomass')
colnames(databis)<-c("simu_ID","I","d","delta","species","biomass_SC","biomass_C")
# CV
databis_2<-biomassCV[biomassCV$model!="NC",-which(names(biomass)%in%c("seed","N","D"))]
databis_2<-databis_2[,-c(which(names(databis_2)=="x1"):which(names(databis_2)=="x5"))]
databis_2<-melt(databis_2, id.vars = c("simu_ID","I","delta","d","model"),
                variable.name = "species",
                value.name = "CV")
databis_2<-dcast(databis_2,simu_ID+I+d+delta+species~model,value.var='CV')
colnames(databis_2)<-c("simu_ID","I","d","delta","species","SC","C")
# stabilising ?
databis$SC<-databis_2$SC
databis$C<-databis_2$C
rm(databis_2)
databis$diff<-databis$SC-databis$C
databis<-databis[databis$biomass_C>0 & databis$biomass_SC>0,]
databis$effect<-"neutral"
treshold=1e-4
databis$effect[(databis$SC-databis$C)>treshold | (databis$biomass_SC==0 & databis$biomass_C>0)]="stabilising"
databis$effect[(databis$SC-databis$C)< -treshold | (databis$biomass_SC>0 & databis$biomass_C==0)]="destabilising"
databis<-databis[databis$effect!="neutral",]
# aggregation
databis<-summarySE(databis, measurevar="diff", groupvars=c("I","d",'delta','effect'),na.rm=TRUE)

p2<-ggplot(data=databis)+
      geom_line(aes(I,diff,colour=effect))+
      stab_colour_short+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(I,diff,colour=effect))+
      geom_errorbar(aes(I,ymin=diff-ci, ymax=diff+ci,colour=effect), width=10, show.legend=FALSE)+
      geom_hline(yintercept=0, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      stab_colour_short+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      xlab(label_nutri)+
      ylab("CV difference between SC and C models")+
      ggtitle("Consumers")

graph<-ggdraw(xlim = c(0, 2.3), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.1, 0, 0.1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 25)

ggsave("supp_stab_quantitative.pdf",graph, width = 14, height = 6, device = cairo_pdf) # 12 10

### Synchrony ####
nparams=6
# recyling
databis<-recy
databis$recy_tot<-rowSums(databis[,(nparams+1):dim(databis)[2]])
databis[,(nparams+1):(dim(databis)[2]-1)]<-databis[,(nparams+1):(dim(databis)[2]-1)]/databis$recy_tot
databis$recy_tot<-NULL
for(i in c((nparams+1):dim(databis)[2])){
  databis[,i]<-databis[,i]*recyCV[,i] # weighted CV
}
databis[is.na(databis)]=0
databis$recy_CV_mean<-rowSums(databis[,(nparams+1):dim(databis)[2]]) # weighted mean
databis$RecyDirCV<-data$RecyDirCV
databis<-databis[databis$model=="C",]
databis$phi<-databis$RecyDirCV^2/databis$recy_CV_mean^2 # synchrony
databis$phi[databis$recy_CV_mean<1e-7]=1
databis<-summarySE(databis, measurevar="phi", groupvars=c("I","d",'delta'),na.rm=TRUE)

p1<-ggplot(data=databis)+
      geom_line(aes(I,phi), size=1.5)+
      geom_errorbar(aes(I,ymin=phi-ci, ymax=phi+ci), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme(panel.grid.major = element_blank(),
            panel.background = element_blank(),
            text = element_text(size=20, family="serif"),
            axis.text = element_text(size=20),
            axis.line = element_line())+
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x)))+
      xlim(0,150)+
      xlab(expression(paste("Mineral nutrient input ",italic("I"))))+
      ylab("Synchrony of direct recycling dynamics")

# biomass
databis<-biomass
databis$PP_tot<-rowSums(databis[,(nparams+nResources+1):(nparams+nResources+nPP)])
databis[,(nparams+nResources+1):(nparams+nResources+nPP)]<-databis[,(nparams+nResources+1):(nparams+nResources+nPP)]/databis$PP_tot
databis$SP_tot<-rowSums(databis[,(nparams+nResources+nPP+1):(nparams+dim)])
databis[,(nparams+nResources+nPP+1):(nparams+dim)]<-databis[,(nparams+nResources+nPP+1):(nparams+dim)]/databis$SP_tot
databis$PP_tot<-NULL
databis$SP_tot<-NULL
for(i in c((nparams+nResources+1):(nparams+dim))){
  databis[,i]<-databis[,i]*biomassCV[,i] # weighted CV
}
databis[is.na(databis)]=0
databis$PP_CV_mean<-rowSums(databis[,(nparams+nResources+1):(nparams+nResources+nPP)]) # weighted mean
databis$SP_CV_mean<-rowSums(databis[,(nparams+nResources+nPP+1):(nparams+dim)]) # weighted mean
databis$PPbiomassCV<-data$PPbiomassCV
databis$SPbiomassCV<-data$SPbiomassCV
databis$phi_PP<-databis$PPbiomassCV^2/databis$PP_CV_mean^2 # synchrony
databis$phi_SP<-databis$SPbiomassCV^2/databis$SP_CV_mean^2 # synchrony
databis$phi_PP[databis$PP_CV_mean<1e-7]=1
databis$phi_SP[databis$SP_CV_mean<1e-7]=1
databis<-databis[,which(names(databis)%in%c("simu_ID","I","delta","d","model","phi_PP","phi_SP"))]
databis<-melt(databis, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "type",
              value.name = "phi")
databis<-summarySE(databis, measurevar="phi", groupvars=c("I","d",'delta',"model","type"),na.rm=TRUE)

p2<-ggplot(data=databis[databis$model!="NC",])+
      geom_line(aes(I,phi,color=model,linetype=type), size=1.5)+
      model_colour_short+
      scale_linetype_manual(values=c("solid", "dashed"),
                            name='',
                            labels=c("Primary\nproducers","Consumers"))+
      theme+theme(legend.key.width = unit(1.5, "cm"),
                  legend.position = "bottom")+
      guides(col = guide_legend(nrow = 2),
             linetype = guide_legend(nrow = 2))
legend<-get_legend(p2)

p2<-ggplot(data=databis[databis$model!="NC",])+
      geom_line(aes(I,phi,color=model,linetype=type), size=1.5)+
      geom_errorbar(aes(I,ymin=phi-ci, ymax=phi+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour_short+
      scale_linetype_manual(values=c("solid", "dashed"),
                            name='',
                            labels=c("Primary\nproducers","Consumers"))+
      theme+theme(legend.position = "none")+
      xlim(0,300)+
      ylim(0,1)+
      xlab(label_nutri)+
      ylab("Synchrony of biomass dynamics")

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1.15)) +
  draw_plot(p1, 0, 0.15, 1, 1)+
  draw_plot(p2, 1, 0.15, 1, 1)+
  draw_plot(legend, 1, 0, 1, 0.2)+
  draw_plot_label(c("A","B"), c(0,1), c(1.15,1.15), size = 25)
ggsave("supp_synchrony.pdf",graph, width = 14, height = 8, device = cairo_pdf)

# FIGURES FINALES - FOOD CHAIN MODEL -----------------------------------------------------------------
nParams=5
nSpecies=4
nResources=2
# CV
CV_chain<-read.table(paste(path_chain,"CV_0.txt",sep=""),sep=';',header=T)
CV_chain$TLmax=2
file<-read.table(paste(path_chain,"CV_1.txt",sep=""),sep=';',header=T)
file$TLmax=3
CV_chain<-rbind(CV_chain,file)
file<-read.table(paste(path_chain,"CV_2.txt",sep=""),sep=';',header=T)
file$TLmax=4
CV_chain<-rbind(CV_chain,file)
# mean
mean_chain<-read.table(paste(path_chain,"mean_0.txt",sep=""),sep=';',header=T)
mean_chain$TLmax=2
file<-read.table(paste(path_chain,"mean_1.txt",sep=""),sep=';',header=T)
file$TLmax=3
mean_chain<-rbind(mean_chain,file)
file<-read.table(paste(path_chain,"mean_2.txt",sep=""),sep=';',header=T)
file$TLmax=4
mean_chain<-rbind(mean_chain,file)
rm(file)

### stabilising effect ####
databis<-CV_chain[CV_chain$model!="NC",c(1:(nParams),dim(CV_chain)[2],nParams+1,(nParams+nResources+1):(nParams+nResources+nSpecies))]
databis<-melt(databis, id.vars = names(databis[,1:(nParams+1)]),
              variable.name = "species", 
              value.name = "CV")
databis<-dcast(databis,simu_ID+I+d+delta+TLmax+species~model)
databis$d<-as.factor(databis$d)
databis$delta<-as.factor(databis$delta)
levels(databis$d)<-c(d02,d08)
levels(databis$delta)<-c(del02,del08)
databis$delta = factor(databis$delta,levels(databis$delta)[c(2,1)])
databis$effect<-databis$SC-databis$C
databis$effect[databis$I>120
               & databis$TLmax==2
               & databis$delta==as.character(del08)]=NA
databis$effect[databis$I>96
               & databis$TLmax==2
               & databis$delta==as.character(del08)
               & databis$d==as.character(d08)]=NA
databis$effect[databis$I>192
               & databis$TLmax==2
               & databis$delta==as.character(del02)
               & databis$d==as.character(d08)]=NA

p1<-ggplot(data=databis[databis$TLmax==2,])+
      geom_line(aes(I,effect,color=species), size=2)+
      geom_hline(yintercept=0, linetype="dashed", size=0.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("CV difference between SC and C models")+
      ggtitle(expression(paste(TL["max"],"=2")))

databis$effect[databis$I>160
               & databis$TLmax==3
               & databis$delta==as.character(del08)]=NA

p2<-ggplot(data=databis[databis$TLmax==3,])+
      geom_line(aes(I,effect,color=species), size=2)+
      geom_hline(yintercept=0, linetype="dashed", size=0.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      ylim(-0.5,0.7)+
      xlab(label_nutri)+
      ylab("CV difference between SC and C models")+
      ggtitle(expression(paste(TL["max"],"=3")))

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_line(aes(I,effect,color=species), size=2)+
      chain_colour+
      theme
legend<-get_legend(p3)

databis$effect[databis$I>144
               & databis$TLmax==4
               & databis$delta==as.character(del08)]=NA

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_line(aes(I,effect,color=species), size=2)+
      geom_hline(yintercept=0, linetype="dashed", size=0.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_colour+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      ylim(-0.8,0.5)+
      xlab(label_nutri)+
      ylab("CV difference between SC and C models")+
      ggtitle(expression(paste(TL["max"],"=4")))

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(legend, 1.13, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,0), c(2,2,1), size = 25)
ggsave("supp_chain_stab.pdf",graph, width = 12, height = 12, device = cairo_pdf)

### recycling ####
databis<-mean_chain[mean_chain$model=="C",c(1:(nParams),dim(mean_chain)[2],24,20:23)]
databis<-melt(databis, id.vars = names(databis[,1:(nParams+1)]),
              variable.name = "species", 
              value.name = "recy")
databis$species = factor(databis$species,levels(databis$species)[seq(5,1,-1)])
databis$d<-as.factor(databis$d)
databis$delta<-as.factor(databis$delta)
levels(databis$d)<-c(d02,d08)
levels(databis$delta)<-c(del02,del08)
databis$delta = factor(databis$delta,levels(databis$delta)[c(2,1)])
databis$recy[databis$I>120
             & databis$TLmax==2
             & databis$delta==as.character(del08)]=0
databis$recy[databis$I>100
             & databis$TLmax==2
             & databis$delta==as.character(del08)
             & databis$d==as.character(d08)]=0
databis$recy[databis$I>192
             & databis$TLmax==2
             & databis$delta==as.character(del02)
             & databis$d==as.character(d08)]=0

p1<-ggplot(data=databis[databis$TLmax==2,])+
      geom_area(aes(I,recy,fill=species))+
      geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_fill_recy+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Recycled nutrients")+
      ggtitle(expression(paste(TL["max"],"=2")))

databis$recy[databis$I>160
             & databis$TLmax==3
             & databis$delta==as.character(del08)]=0

p2<-ggplot(data=databis[databis$TLmax==3,])+
      geom_area(aes(I,recy,fill=species))+
      geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_fill_recy+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Recycled nutrients")+
      ggtitle(expression(paste(TL["max"],"=3")))

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_area(aes(I,recy,fill=species))+
      chain_fill_recy+
      theme
legend<-get_legend(p3)

databis$recy[databis$I>144
             & databis$TLmax==4
             & databis$delta==as.character(del08)]=0

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_area(aes(I,recy,fill=species))+
      geom_abline(intercept=0, slope=1, linetype="dashed", size=1)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_fill_recy+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Recycled nutrients")+
      ggtitle(expression(paste(TL["max"],"=4")))

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(legend, 1.13, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,0), c(2,2,1), size = 25)
ggsave("supp_chain_recy.pdf",graph, width = 12, height = 12, device = cairo_pdf)

### biomass ####
databis<-mean_chain[mean_chain$model=="C",c(1:(nParams),dim(mean_chain)[2],8:11)]
databis<-melt(databis, id.vars = names(databis[,1:(nParams+1)]),
              variable.name = "species", 
              value.name = "biomass")
databis$species = factor(databis$species,levels(databis$species)[seq(4,1,-1)])
databis$d<-as.factor(databis$d)
databis$delta<-as.factor(databis$delta)
levels(databis$d)<-c(d02,d08)
levels(databis$delta)<-c(del02,del08)
databis$delta = factor(databis$delta,levels(databis$delta)[c(2,1)])
databis$biomass[databis$I>120
                & databis$TLmax==2
                & databis$delta==as.character(del08)]=0
databis$biomass[databis$I>100
                & databis$TLmax==2
                & databis$delta==as.character(del08)
                & databis$d==as.character(d08)]=0
databis$biomass[databis$I>192
                & databis$TLmax==2
                & databis$delta==as.character(del02)
                & databis$d==as.character(d08)]=0

p1<-ggplot(data=databis[databis$TLmax==2,])+
      geom_area(aes(I,biomass,fill=species))+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_fill_biomass+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Biomass")+
      ggtitle(expression(paste(TL["max"],"=2")))

databis$biomass[databis$I>160
                & databis$TLmax==3
                & databis$delta==as.character(del08)]=0

p2<-ggplot(data=databis[databis$TLmax==3,])+
      geom_area(aes(I,biomass,fill=species))+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_fill_biomass+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Biomass")+
      ggtitle(expression(paste(TL["max"],"=3")))

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_area(aes(I,biomass,fill=species))+
      chain_fill_biomass+
      theme+
      lim_chain
legend<-get_legend(p3)

databis$biomass[databis$I>144
                & databis$TLmax==4
                & databis$delta==as.character(del08)]=0

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_area(aes(I,biomass,fill=species))+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_fill_biomass+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Biomass")+
      ggtitle(expression(paste(TL["max"],"=4")))

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(legend, 1.13, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,0), c(2,2,1), size = 25)
ggsave("supp_chain_biomass.pdf",graph, width = 12, height = 12, device = cairo_pdf)

### CV ####
databis<-CV_chain[CV_chain$model=="C",c(1:(nParams),dim(CV_chain)[2],24,6,8:11)]
databis<-melt(databis, id.vars = names(databis[,1:(nParams+1)]),
              variable.name = "species", 
              value.name = "CV")
databis$species = factor(databis$species,levels(databis$species)[seq(6,1,-1)])
databis$d<-as.factor(databis$d)
databis$delta<-as.factor(databis$delta)
levels(databis$d)<-c(d02,d08)
levels(databis$delta)<-c(del02,del08)
databis$delta = factor(databis$delta,levels(databis$delta)[c(2,1)])
databis$CV[databis$I>120
           & databis$TLmax==2
           & databis$delta==as.character(del08)]=NA
databis$CV[databis$I>96
           & databis$TLmax==2
           & databis$delta==as.character(del08)
           & databis$d==as.character(d08)]=NA
databis$CV[databis$I>192
           & databis$TLmax==2
           & databis$delta==as.character(del02)
           & databis$d==as.character(d08)]=NA

p1<-ggplot(data=databis[databis$TLmax==2,])+
      geom_line(aes(I,CV,colour=species,linetype=species),size=1.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_colour_CV+
      chain_line_CV+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      xlab(label_nutri)+
      ylab("Biomass CV")+
      ggtitle(expression(paste(TL["max"],"=2")))

databis$CV[databis$I>155
           & databis$TLmax==3
           & databis$delta==as.character(del08)]=NA
databis$CV[databis$I>130
           & databis$TLmax==3
           & databis$delta==as.character(del08)
           & databis$d==as.character(d08)]=NA

p2<-ggplot(data=databis[databis$TLmax==3,])+
      geom_line(aes(I,CV,colour=species,linetype=species),size=1.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_colour_CV+
      chain_line_CV+   
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      ylim(0,3.5)+
      xlab(label_nutri)+
      ylab("Biomass CV")+
      ggtitle(expression(paste(TL["max"],"=3")))

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_line(aes(I,CV,colour=species,linetype=species),size=1.5)+
      chain_colour_CV+
      chain_line_CV+
      theme
legend<-get_legend(p3)

databis$CV[databis$I>138
           & databis$TLmax==4
           & databis$delta==as.character(del08)]=NA
databis$CV[databis$I>105
           & databis$TLmax==4
           & databis$delta==as.character(del08)
           & databis$d==as.character(d08)]=NA

p3<-ggplot(data=databis[databis$TLmax==4,])+
      geom_line(aes(I,CV,colour=species,linetype=species),size=1.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      chain_colour_CV+
      chain_line_CV+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_chain+
      ylim(0,3.5)+
      xlab(label_nutri)+
      ylab("Biomass CV")+
      ggtitle(expression(paste(TL["max"],"=4")))

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(legend, 1.13, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C"), c(0,1,0), c(2,2,1), size = 25)
ggsave("supp_chain_CV.pdf",graph, width = 12, height = 12, device = cairo_pdf)

# FIGURES FINALES - SENSITIVITY ANALYSIS -----------------------------------------------------------------
### Sensitivity replicates and time ####
# number of replicates ####
nParams=4
data_rep<-read.table(paste(path_rep,"data.txt",sep=""),sep=';',header=T)
data_rep$persistence<-data_rep$NbSpeciesFinal/data_rep$NbSpeciesInit
data_rep$model<-data_rep$d
data_rep$model[data_rep$model==0.8]=0.2
data_rep$model<-as.factor(data_rep$model)
levels(data_rep$model)=c("SC","NC","C")
for (i in seq(from=1,to=dim(data_rep)[1]-2,by=3)){
  data_rep$d[i] = data_rep$d[i+1]
  data_rep$d[i+2] = data_rep$d[i+1]
  data_rep$delta[i] = data_rep$delta[i+1]
  data_rep$delta[i+2] = data_rep$delta[i+1]
}
data_rep$d<-as.factor(data_rep$d)
data_rep$delta<-as.factor(data_rep$delta)
levels(data_rep$d)<-c(d02,d08)
levels(data_rep$delta)<-c(del02,del08)
data_rep$delta = factor(data_rep$delta,levels(data_rep$delta)[c(2,1)])
data_rep$model = factor(data_rep$model,levels(data_rep$model)[c(3,2,1)])

databis<-NULL
# N=100
databis_2<-summarySE(data, measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="persistence"
databis<-rbind(databis,databis_2)
databis_2<-summarySE(data, measurevar="Irecy", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="Irecy"
databis<-rbind(databis,databis_2)
databis_2<-summarySE(data, measurevar="PPprod", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="PPprod"
databis<-rbind(databis,databis_2)
databis_2<-summarySE(data, measurevar="SPprod", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="SPprod"
databis<-rbind(databis,databis_2)
databis_2<-summarySE(data, measurevar="RecyPP", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="RecyPP"
databis<-rbind(databis,databis_2)
databis_2<-summarySE(data, measurevar="RecySP", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="RecySP"
databis<-rbind(databis,databis_2)
# N=200
databis_2<-data[,which(names(data)%in%c("I","d",'delta','model',"persistence"))]
databis_2_bis<-data_rep[,which(names(data_rep)%in%c("I","d",'delta','model',"persistence"))]
databis_2<-rbind(databis_2,databis_2_bis)
databis_2<-summarySE(databis_2, measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="persistence"
databis<-rbind(databis,databis_2)
databis_2<-data[,which(names(data)%in%c("I","d",'delta','model',"Irecy"))]
databis_2_bis<-data_rep[,which(names(data_rep)%in%c("I","d",'delta','model',"Irecy"))]
databis_2<-rbind(databis_2,databis_2_bis)
databis_2<-summarySE(databis_2, measurevar="Irecy", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="Irecy"
databis<-rbind(databis,databis_2)
databis_2<-data[,which(names(data)%in%c("I","d",'delta','model',"PPprod"))]
databis_2_bis<-data_rep[,which(names(data_rep)%in%c("I","d",'delta','model',"PP"))]
names(databis_2_bis)[1]="PPprod"
databis_2<-rbind(databis_2,databis_2_bis)
databis_2<-summarySE(databis_2, measurevar="PPprod", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="PPprod"
databis<-rbind(databis,databis_2)
databis_2<-data[,which(names(data)%in%c("I","d",'delta','model',"SPprod"))]
databis_2_bis<-data_rep[,which(names(data_rep)%in%c("I","d",'delta','model',"SP"))]
names(databis_2_bis)[1]="SPprod"
databis_2<-rbind(databis_2,databis_2_bis)
databis_2<-summarySE(databis_2, measurevar="SPprod", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="SPprod"
databis<-rbind(databis,databis_2)
databis_2<-data[,which(names(data)%in%c("I","d",'delta','model',"RecyPP"))]
databis_2_bis<-data_rep[,which(names(data_rep)%in%c("I","d",'delta','model',"RecyPP"))]
databis_2<-rbind(databis_2,databis_2_bis)
databis_2<-summarySE(databis_2, measurevar="RecyPP", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="RecyPP"
databis<-rbind(databis,databis_2)
databis_2<-data[,which(names(data)%in%c("I","d",'delta','model',"RecySP"))]
databis_2_bis<-data_rep[,which(names(data_rep)%in%c("I","d",'delta','model',"RecySP"))]
databis_2<-rbind(databis_2,databis_2_bis)
databis_2<-summarySE(databis_2, measurevar="RecySP", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis_2<-databis_2[,c(1:(nParams+3))]
names(databis_2)[nParams+2]="mean"
databis_2$variable="RecySP"
databis<-rbind(databis,databis_2)
databis$N<-as.factor(databis$N)
databis$variable<-as.factor(databis$variable)
# relative mean and sd difference
databis_2<-dcast(databis,I+d+delta+model+variable~N,value.var='mean')
databis_2$mean_dif<-abs(databis_2$'100'-databis_2$'200')/databis_2$'100' # relative difference
databis_2_bis<-dcast(databis,I+d+delta+model+variable~N,value.var='sd')
databis_2_bis$sd_dif<-abs(databis_2_bis$'100'-databis_2_bis$'200')/databis_2_bis$'100' # relative difference
databis<-merge(databis_2,databis_2_bis,by=c("I","d",'delta','model','variable'))
rm(databis_2_bis)

databis_2<-summarySE(databis, measurevar="mean_dif", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
p1<-ggplot(data=databis_2)+
      geom_line(aes(I,mean_dif,colour=model,linetype=model), size=1.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position= "none",
                  panel.grid.major.y = element_line(colour='grey'))+
      xlim(0,200)+
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    limits = c(1e-3,1))+
      xlab(label_nutri)+
      ylab("Relative mean difference")

databis_2<-summarySE(databis, measurevar="sd_dif", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
p2<-ggplot(data=databis_2)+
      geom_line(aes(I,sd_dif,colour=model,linetype=model), size=1.5)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5, "cm"),
                  legend.position= "bottom")
legend<-get_legend(p2)

p2<-ggplot(data=databis_2)+
      geom_line(aes(I,sd_dif,colour=model,linetype=model), size=1.5)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  panel.grid.major.y = element_line(colour='grey'))+
      xlim(0,200)+
      scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    limits = c(1e-3,1))+
      xlab(label_nutri)+
      ylab("Relative standard deviation difference")

# number of extinction depending on time ####
databis<-melt(tExt, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "tExt")
databis<-databis[databis$tExt>0 & databis$tExt<9000,]

p3<-ggplot(data=databis)+
      geom_freqpoly(aes(tExt),binwidth = 100)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_y_continuous(breaks = c(0,1e5),
                         labels = c("0",expression(10^"5")))+
      xlab("Time")+
      ylab("Number of extinction")

# Final graph ####
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(legend, 1, 0.9, 1, 0.1)+
  draw_plot_label(c("A","B","C"), c(0,1,0), c(2,2,1), size = 25)
ggsave("supp_sensitivity_method.pdf",graph, width = 14, height = 12, device = cairo_pdf)

#### extinction threshold ####
data_thr<-read.table(paste(path_thresh,"data.txt",sep=""),sep=';',header=T)
data_thr$persistence<-data_thr$NbSpeciesFinal/data_thr$NbSpeciesInit
data_thr$d<-as.factor(data_thr$d)
data_thr$delta<-as.factor(data_thr$delta)
levels(data_thr$d)<-c(d02,d08)
levels(data_thr$delta)<-c(del02,del08)
data_thr$delta = factor(data_thr$delta,levels(data_thr$delta)[c(2,1)])
data_thr$model = factor(data_thr$model,levels(data_thr$model)[c(3,2,1)])

databis <- summarySE(data_thr, measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
p1<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=persistence-ci, ymax=persistence+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_short+
      ylim(0,0.8)+
      ggtitle(expression(paste("Extinction threshold = ",10^"-15")))+
      xlab(label_nutri)+
      ylab("Species persistence")

databis <- summarySE(data[data$d==as.character(d02) & data$delta==as.character(del02),], measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
p2<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      model_colour+
      model_line+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=persistence-ci, ymax=persistence+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_short+
      ylim(0,0.8)+
      ggtitle(expression(paste("Extinction threshold = ",10^"-30")))+
      xlab(label_nutri)+
      ylab("Species persistence")

# biomass ####
databis<-melt(biomass, id.vars = c("simu_ID","I","delta","d","model"),
              variable.name = "species",
              value.name = "biomass")
databis<-databis[databis$biomass>0,]

p3<-ggplot(data=databis[databis$model=="C" & databis$d==as.character(d02) & databis$delta==as.character(del02),])+
  geom_point(aes(I,biomass),alpha=0.2, size=1.5)+
  theme_void()+theme(legend.position = "none")+
  lim_short+
  y_axis_log10+
  xlab(label_nutri)+
  ylab("Species average biomass")
ggsave(p3,file="graphe.png")
image <- readPNG("graphe.png")

p3<-ggplot(data=databis[databis$model=="C" & databis$d==as.character(d02) & databis$delta==as.character(del02),])+
  geom_blank(aes(I,biomass))+
  annotation_raster(image, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, interpolate=TRUE)+
  facet_grid(delta~d, labeller=label_parsed)+
  theme+theme(legend.position = "none")+
  lim_short+
  y_axis_log10+
  xlab(label_nutri)+
  ylab("Species average biomass")

# Final graph ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p2, 0, 0, 1, 1)+
  #draw_plot(p3, 1, 1, 1, 1)+
  draw_plot(p1, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0, 0.05, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 25)
ggsave("supp_sensitivity_threshold.pdf",graph, width = 14, height = 6, device = cairo_pdf)

### Sensitivity beta ####
# Beta = 0
# Persistence
data_beta<-read.table(paste(path_beta,"data_0.txt",sep=""),sep=';',header=T)
data_beta$persistence<-data_beta$NbSpeciesFinal/data_beta$NbSpeciesInit
data_beta$d<-as.factor(data_beta$d)
data_beta$delta<-as.factor(data_beta$delta)
levels(data_beta$d)<-c(d02,d08)
levels(data_beta$delta)<-c(del02,del08)
data_beta$delta = factor(data_beta$delta,levels(data_beta$delta)[c(2,1)])
data_beta$model = factor(data_beta$model,levels(data_beta$model)[c(3,2,1)])

databis <- summarySE(data_beta, measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
p1<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=persistence-ci, ymax=persistence+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_short+
      ggtitle(expression(paste(italic("\u03B2"),"=0")))+
      xlab(label_nutri)+
      ylab("Species persistence")

# Biomass CV
biomass_beta<-read.table(paste(path_beta,"biomass_0.txt",sep=""),sep=';',header=T)
biomass_beta$d<-as.factor(biomass_beta$d)
biomass_beta$delta<-as.factor(biomass_beta$delta)
levels(biomass_beta$d)<-c(d02,d08)
levels(biomass_beta$delta)<-c(del02,del08)
biomass_beta$delta = factor(biomass_beta$delta,levels(biomass_beta$delta)[c(2,1)])
biomass_beta$model = factor(biomass_beta$model,levels(biomass_beta$model)[c(3,2,1)])

biomassCV_beta<-read.table(paste(path_beta,"biomassCV_0.txt",sep=""),sep=';',header=T)
biomassCV_beta$d<-as.factor(biomassCV_beta$d)
biomassCV_beta$delta<-as.factor(biomassCV_beta$delta)
levels(biomassCV_beta$d)<-c(d02,d08)
levels(biomassCV_beta$delta)<-c(del02,del08)
biomassCV_beta$delta = factor(biomassCV_beta$delta,levels(biomassCV_beta$delta)[c(2,1)])
biomassCV_beta$model = factor(biomassCV_beta$model,levels(biomassCV_beta$model)[c(3,2,1)])

nparams=5
databis<-biomass_beta[,-which(names(biomass_beta)%in%c("seed","N","D"))]
databis$biomass_tot<-rowSums(databis[,(nparams+1):dim(databis)[2]])
databis[,(nparams+1):(dim(databis)[2]-1)]<-databis[,(nparams+1):(dim(databis)[2]-1)]/databis$biomass_tot
databis$biomass_tot<-NULL
# biomass CV
databis_2<-biomassCV_beta[,-which(names(biomass)%in%c("seed","N","D"))]
databis_2$model<-droplevels(databis_2$model)
# CV weighted by the average biomass
for(i in c((nparams+1):dim(databis)[2])){
  databis[,i]<-databis[,i]*databis_2[,i]
}
databis[is.na(databis)]=0
databis$CV<-rowSums(databis[,(nparams+1):dim(databis)[2]])
# persistence
databis_2<-data_beta
databis$persistence<-databis_2$persistence
rm(databis_2)
# Average CV 
databis<-databis[databis$persistence>0,]
databis<-summarySE(databis, measurevar="CV", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
databis$ci[databis$ci>2]=0

p3<-ggplot(data=databis)+
      geom_line(aes(I,CV,colour=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=CV-ci, ymax=CV+ci,colour=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_short+
      ylim(-0.3,1.6)+
      ggtitle(expression(paste(italic("\u03B2"),"=0")))+
      xlab(label_nutri)+
      ylab("Average weighted biomass CV")

# Beta = 0.1
# Persistence
data_beta<-read.table(paste(path_beta,"data_1.txt",sep=""),sep=';',header=T)
data_beta$persistence<-data_beta$NbSpeciesFinal/data_beta$NbSpeciesInit
data_beta$d<-as.factor(data_beta$d)
data_beta$delta<-as.factor(data_beta$delta)
levels(data_beta$d)<-c(d02,d08)
levels(data_beta$delta)<-c(del02,del08)
data_beta$delta = factor(data_beta$delta,levels(data_beta$delta)[c(2,1)])
data_beta$model = factor(data_beta$model,levels(data_beta$model)[c(3,2,1)])

databis <- summarySE(data_beta, measurevar="persistence", groupvars=c("I","d",'delta','model'),na.rm=TRUE)
p2<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5,"cm"))
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(I,persistence,color=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=persistence-ci, ymax=persistence+ci,color=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_short+
      ggtitle(expression(paste(italic("\u03B2"),"=0.1")))+
      xlab(label_nutri)+
      ylab("Species persistence")

# Biomass CV
biomass_beta<-read.table(paste(path_beta,"biomass_1.txt",sep=""),sep=';',header=T)
biomass_beta$d<-as.factor(biomass_beta$d)
biomass_beta$delta<-as.factor(biomass_beta$delta)
levels(biomass_beta$d)<-c(d02,d08)
levels(biomass_beta$delta)<-c(del02,del08)
biomass_beta$delta = factor(biomass_beta$delta,levels(biomass_beta$delta)[c(2,1)])
biomass_beta$model = factor(biomass_beta$model,levels(biomass_beta$model)[c(3,2,1)])

biomassCV_beta<-read.table(paste(path_beta,"biomassCV_1.txt",sep=""),sep=';',header=T)
biomassCV_beta$d<-as.factor(biomassCV_beta$d)
biomassCV_beta$delta<-as.factor(biomassCV_beta$delta)
levels(biomassCV_beta$d)<-c(d02,d08)
levels(biomassCV_beta$delta)<-c(del02,del08)
biomassCV_beta$delta = factor(biomassCV_beta$delta,levels(biomassCV_beta$delta)[c(2,1)])
biomassCV_beta$model = factor(biomassCV_beta$model,levels(biomassCV_beta$model)[c(3,2,1)])

nparams=5
databis<-biomass_beta[,-which(names(biomass_beta)%in%c("seed","N","D"))]
databis$biomass_tot<-rowSums(databis[,(nparams+1):dim(databis)[2]])
databis[,(nparams+1):(dim(databis)[2]-1)]<-databis[,(nparams+1):(dim(databis)[2]-1)]/databis$biomass_tot
databis$biomass_tot<-NULL
# biomass CV
databis_2<-biomassCV_beta[,-which(names(biomass)%in%c("seed","N","D"))]
databis_2$model<-droplevels(databis_2$model)
# CV weighted by the average biomass
for(i in c((nparams+1):dim(databis)[2])){
  databis[,i]<-databis[,i]*databis_2[,i]
}
databis[is.na(databis)]=0
databis$CV<-rowSums(databis[,(nparams+1):dim(databis)[2]])
# persistence
databis_2<-data_beta
databis$persistence<-databis_2$persistence
rm(databis_2)
# Average CV 
databis<-databis[databis$persistence>0,]
databis<-summarySE(databis, measurevar="CV", groupvars=c("I","d",'delta','model'),na.rm=TRUE)

p4<-ggplot(data=databis)+
      geom_line(aes(I,CV,colour=model,linetype=model), size=1.5)+
      geom_errorbar(aes(I,ymin=CV-ci, ymax=CV+ci,colour=model), width=10, show.legend=FALSE)+
      facet_grid(delta~d, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5))+
      lim_short+
      ggtitle(expression(paste(italic("\u03B2"),"=0.1")))+
      xlab(label_nutri)+
      ylab("Average weighted biomass CV")

graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p3, 1, 1, 1, 1)+
  draw_plot(p2, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend, 2.07, 1, 0.05, 1)+
  draw_plot(legend, 2.07, 0, 0.05, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 25)
ggsave("supp_sensitivity_beta.pdf",graph, width = 14, height = 12, device = cairo_pdf)

### Attack rate and beta  ####
data2<-read.table(paste(path_aB,"data.txt",sep=""),sep=';',header=T)
data2[is.na(data2)] <- 0
data2$persistence<-(data2$NbSpeciesFinal/data2$NbSpeciesInit)
data2$d<-as.factor(1)
data2$delta<-as.factor(1)
levels(data2$d)<-c(d02)
levels(data2$delta)<-c(del02)

# Persistence
databis<-summarySE(data2[data2$model=="C",], measurevar="persistence", groupvars=c("a","B","I","d","delta","model"))
p1<-ggplot(data=databis)+
      geom_raster(aes(a,B,fill=persistence))+
      geom_point(aes(0.1,0.001),color="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradient(low = "light blue", high = "black","Species\npersistence")+
      x_axis_log10+
      y_axis_log10+
      xlab(expression(paste("Attack rate allometric constant ",italic("a"))))+
      ylab(expression(atop("Density dependent mortality",paste("rate allometric constant ",italic("\u03B2")))))

databis<-dcast(data2,simu_ID+a+B+I+d+delta~model,value.var="persistence")
databis$dif<-abs(databis$C-databis$SC)
databis<-summarySE(databis, measurevar="dif", groupvars=c("a","B","I","d","delta"))
p3<-ggplot(data=databis)+
      geom_raster(aes(a,B,fill=dif))+
      geom_point(aes(0.1,0.001),color="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradient(low = "white", high = "blue","Species\npersistence\ndifference\nbetween C\nand SC\nmodels")+
      x_axis_log10+
      y_axis_log10+
      xlab(expression(paste("Attack rate allometric constant ",italic("a"))))+
      ylab(expression(atop("Density dependent mortality",paste("rate allometric constant ",italic("\u03B2")))))

# Coefficient of variation
nparams=5
biomass2<-read.table(paste(path_aB,"biomass.txt",sep=""),sep=';',header=T)
biomass2<-biomass2[,-which(names(biomass2)%in%c("seed","N","D"))]
biomass2$biomass_tot<-rowSums(biomass2[,(nparams+1):dim(biomass2)[2]])
biomass2[,(nparams+1):(dim(biomass2)[2]-1)]<-biomass2[,(nparams+1):(dim(biomass2)[2]-1)]/biomass2$biomass_tot
biomass2$biomass_tot<-NULL
# biomass CV
biomassCV2<-read.table(paste(path_aB,"biomassCV.txt",sep=""),sep=';',header=T)
biomassCV2<-biomassCV2[,-which(names(biomassCV2)%in%c("seed","N","D"))]
# CV weighted by the average biomass
for(i in c((nparams+1):dim(databis)[2])){
  biomassCV2[,i]<-biomassCV2[,i]*biomass2[,i]
}
biomassCV2[is.na(biomassCV2)]=0
biomassCV2$CV<-rowSums(biomassCV2[,(nparams+1):dim(biomassCV2)[2]])
# persistence
biomassCV2$persistence<-data2$persistence
biomassCV2$d<-as.factor(1)
biomassCV2$delta<-as.factor(1)
levels(biomassCV2$d)<-c(d02)
levels(biomassCV2$delta)<-c(del02)
void<-expand.grid(a=unique(biomassCV2$a),
                  B=unique(biomassCV2$B),
                  I=unique(biomassCV2$I),
                  d=unique(biomassCV2$d),
                  delta=unique(biomassCV2$delta))
biomassCV2<-biomassCV2[biomassCV2$persistence>0,]

databis<-summarySE(biomassCV2[biomassCV2$model=="C",], measurevar="CV", groupvars=c("a","B","I","d",'delta'),na.rm=TRUE)
databis<-merge(databis,void,by=c("a","B","I","d","delta"),all.y=T)
databis$CV[is.na(databis$CV)]=0
p2<-ggplot(data=databis)+
      geom_raster(aes(a,B,fill=CV))+
      geom_point(aes(0.1,0.001),color="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradient(low = "light blue", high = "black","Average\nbiomass CV")+
      x_axis_log10+
      y_axis_log10+
      xlab(expression(paste("Attack rate allometric constant ",italic("a"))))+
      ylab(expression(atop("Density dependent mortality",paste("rate allometric constant ",italic("\u03B2")))))

databis<-biomassCV2
databis$persistence<-NULL
databis<-melt(databis,id.vars = c("simu_ID","a","B","I","d","delta","model"),
              variable.name = "species", 
              value.name = "CV")
databis<-dcast(databis,simu_ID+a+B+I+d+delta+species~model,value.var="CV")
databis$dif<-abs(databis$SC-databis$C)/databis$SC
databis<-databis[is.na(databis$dif)==F,]
threshold=1e-4
databis$dif[databis$SC<threshold]=0
databis$dif[databis$C<threshold]=0
databis<-summarySE(databis, measurevar="dif", groupvars=c("a","B","I","d","delta"))
databis<-merge(databis,void,by=c("a","B","I","d","delta"),all.y=T)
databis$dif[is.na(databis$dif)]=0
p4<-ggplot(data=databis)+
      geom_raster(aes(a,B,fill=dif))+
      geom_point(aes(0.1,0.001),color="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradient2(low = "red", mid = "white" , high = "blue","Relative\nbiomass CV\ndifference\nbetween C\nand SC\nmodels")+
      x_axis_log10+
      y_axis_log10+
      xlab(expression(paste("Attack rate allometric constant ",italic("a"))))+
      ylab(expression(atop("Density dependent mortality",paste("rate allometric constant ",italic("\u03B2")))))


# Food web regims
databis<-summarySE(biomassCV2[biomassCV2$model=="C",], measurevar="CV", groupvars=c("a","B","I","d",'delta'),na.rm=TRUE)
databis<-merge(databis,void,by=c("a","B","I","d","delta"),all.y=T)
databis$CV[is.na(databis$CV)]=0
databis2<-summarySE(data2[data2$model=="C",], measurevar="TLmax", groupvars=c("a","B","I","d","delta","model"))
databis2<-databis2[,which(names(databis2)%in%c("a","B","I","d","delta","model","TLmax"))]
databis<-merge(databis,databis2,by=c("a","B","I","d","delta"))
rm(databis2)

p5<-ggplot(data=databis)+
      geom_raster(aes(a,B,fill=TLmax))+
      geom_point(aes(0.1,0.001),colour="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradientn(name="Maximum\ntrophic\nlevel",
                           colours=c("white","chartreuse3","cadetblue4","red3","darkred"))+
      x_axis_log10+
      y_axis_log10+
      xlab(expression(paste("Attack rate allometric constant ",italic("a"))))+
      ylab(expression(atop("Density dependent mortality",paste("rate allometric constant ",italic("\u03B2")))))

databis$TLmax[databis$TLmax<1]=floor(databis$TLmax[databis$TLmax<1])
databis$regime="0"
threshold=1e-4
databis$regime[databis$TLmax==0]="collapse"
databis$regime[databis$CV<threshold & databis$TLmax>0]="fixed points"
databis$regime[databis$CV>threshold & databis$TLmax>0]="limit cycles"

p6<-ggplot(data=databis)+
      geom_raster(aes(a,B,fill=regime))+
      geom_point(aes(0.1,0.001),colour="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_manual(name = "Regime",
                        values = c("grey","gold2","darkorchid4"))+
      x_axis_log10+
      y_axis_log10+
      xlab(expression(paste("Attack rate allometric constant ",italic("a"))))+
      ylab(expression(atop("Density dependent mortality",paste("rate allometric constant ",italic("\u03B2")))))

# Final graph
graph<-plot_grid(p1, p2, p3, p4, p5, p6,
                 labels = c("A","B","C","D","E","F"), label_size = 25,
                 nrow = 3, align = "hv")
ggsave("supp_sensitivity_a_B.pdf",graph, width = 15, height = 15, device = cairo_pdf)
#ggsave("supp_sensitivity_a_B.png",graph, width = 15, height = 15)

### Leaching rate and half saturation constant ####
data1<-read.table(paste(path_KLaB,"data0.txt",sep=""),sep=';',header=T)
file<-read.table(paste(path_KLaB,"data1.txt",sep=""),sep=';',header=T)
data1<-rbind(data1,file)
rm(file)

data1[is.na(data1)] <- 0
data1$persistence<-(data1$NbSpeciesFinal/data1$NbSpeciesInit)
data1$d<-as.factor(data1$d)
data1$delta<-as.factor(data1$delta)
levels(data1$d)<-c(d02)
levels(data1$delta)<-c(del02)

# Persistence
databis <- summarySE(data1, measurevar="persistence", groupvars=c("K","L","I","d","delta"))
databis$K[databis$K==1]=0
p1<-ggplot(data=databis)+
      geom_raster(aes(K,L,fill=persistence))+
      geom_point(aes(10,0.2),color="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradient(low = "light blue", high = "black","Species\npersistence")+
      y_axis_log10+
      xlab(expression(paste("Half saturation of nutrients uptake ",italic("K"))))+
      ylab(expression(paste("Loss rate of nutrients ",italic("\u2113"))))

# Coefficient of variation
databis <- summarySE(data1[data1$persistence>0,], measurevar="CV", groupvars=c("K","L","I","d","delta"))
databis$K[databis$K==1]=0
p2<-ggplot(data=databis)+
      geom_raster(aes(K,L,fill=CV))+
      geom_point(aes(10,0.2),color="red",size=5)+
      facet_grid(delta~d, labeller=label_parsed)+
      theme+
      scale_fill_gradient(low = "light blue", high = "black","Average\nbiomass CV")+
      y_axis_log10+
      xlab(expression(paste("Half saturation of nutrients uptake ",italic("K"))))+
      ylab(expression(paste("Loss rate of nutrients ",italic("\u2113"))))

# Final graph
graph<-plot_grid(p1, p2,
                 labels = c("A","B"), label_size = 25,
                 nrow = 1, align = "hv")
ggsave("supp_sensitivity_K_l.pdf",graph, width = 14, height = 5, device = cairo_pdf)

### A ####
data1<-read.table(paste(path_A,"data.txt",sep=""),sep=';',header=T)
data1$persistence<-as.numeric(data1$persistence)
data1$d<-as.factor(data1$d)
data1$A<-as.factor(data1$A)
levels(data1$A)<-c(expression(paste(italic("A"),"=0")),
                   expression(paste(italic("A"),"=0.001")),
                   expression(paste(italic("A"),"=0.01")),
                   expression(paste(italic("A"),"=0.1")))

databis <- summarySE(data1, measurevar="persistence", groupvars=c("I","d","A"),na.rm=TRUE)
p1<-ggplot()+
      geom_line(data=databis,aes(I,persistence,color=d,linetype=d), size=2)+
      geom_errorbar(data=databis,aes(I,ymin=persistence-ci, ymax=persistence+ci,color=d), width=15, show.legend=FALSE)+
      facet_wrap(~A, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none")+
      lim_short+
      xlab(label_nutri)+
      ylab("Species persistence")

databis <- summarySE(data1[data1$persistence>0,], measurevar="CV", groupvars=c("I","d","A"),na.rm=TRUE)
p2<-ggplot()+
      geom_line(data=databis,aes(I,CV,color=d,linetype=d), size=2)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5, "cm"),
                  legend.position = "bottom")
legend<-get_legend(p2)

p2<-ggplot()+
      geom_line(data=databis,aes(I,CV,color=d,linetype=d), size=2)+
      geom_errorbar(data=databis,aes(I,ymin=CV-ci, ymax=CV+ci,color=d), width=15, show.legend=FALSE)+
      facet_wrap(~A, labeller=label_parsed)+
      model_colour+
      model_line+
      theme+theme(legend.position = "none")+
      lim_short+
      ylim(c(0,2.2))+
      xlab(label_nutri)+
      ylab("Average biomass CV")

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1.1)) +
  draw_plot(p1, 0, 0.1, 1, 1)+
  draw_plot(p2, 1, 0.1, 1, 1)+
  draw_plot(legend, 0.5, 0, 1, 0.1)+
  draw_plot_label(c("A","B"), c(0,1), c(1.1,1.1), size = 25)
ggsave("supp_A.pdf",graph, width = 12, height = 7, device = cairo_pdf)

### Functional response ####
data1<-read.table(paste(path_FR,"data.txt",sep=""),sep=';',header=T)
data1$persistence<-as.numeric(data1$persistence)
data1$model<-data1$d
data1$model<-as.factor(data1$model)
levels(data1$model)=c("SC","NC","C")
for (i in seq(from=1,to=dim(data1)[1]-2,by=3)){
  data1$d[i] = data1$d[i+1]
  data1$d[i+2] = data1$d[i+1]
  data1$delta[i] = data1$delta[i+1]
  data1$delta[i+2] = data1$delta[i+1]
}
data1$d<-as.factor(data1$d)
data1$delta<-as.factor(data1$delta)
levels(data1$d)<-c(d02)
levels(data1$delta)<-c(del02)

databis <- summarySE(data1, measurevar="persistence", groupvars=c("I","d","delta","model"),na.rm=TRUE)
p1<-ggplot()+
  geom_line(data=databis,aes(I,persistence,colour=model,linetype=model), size=2)+
  geom_errorbar(data=databis,aes(I,ymin=persistence-ci, ymax=persistence+ci,colour=model), width=15, show.legend=FALSE)+
  facet_grid(delta~d, labeller=label_parsed)+
  model_colour+
  model_line+
  theme+theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))+
  lim_short+
  ggtitle("Type III functional response")+
  xlab(label_nutri)+
  ylab("Species persistence")

databis <- summarySE(data1[data1$persistence>0,], measurevar="CV", groupvars=c("I","d","delta","model"),na.rm=TRUE)
p2<-ggplot()+
  geom_line(data=databis,aes(I,CV,colour=model,linetype=model), size=2)+
  model_colour+
  model_line+
  theme+theme(legend.key.width = unit(1.5, "cm"))
legend<-get_legend(p2)

p2<-ggplot()+
  geom_line(data=databis,aes(I,CV,colour=model,linetype=model), size=2)+
  geom_errorbar(data=databis,aes(I,ymin=CV-ci, ymax=CV+ci,colour=model), width=15, show.legend=FALSE)+
  facet_grid(delta~d, labeller=label_parsed)+
  model_colour+
  model_line+
  theme+theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5))+
  lim_short+
  ggtitle("Type III functional response")+
  xlab(label_nutri)+
  ylab("Average biomass CV")

graph<-ggdraw(xlim = c(0, 2.2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.07, 0, 0.05, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 25)
ggsave("supp_FR.pdf",graph, width = 12, height = 5, device = cairo_pdf)

### Primary producers C:N ####
data1<-read.table(paste(path_CN,"data.txt",sep=""),sep=';',header=T)
data1$persistence<-as.numeric(data1$persistence)
data1$model<-data1$d
data1$model<-as.factor(data1$model)
levels(data1$model)=c("SC","NC","C")
for (i in seq(from=1,to=dim(data1)[1]-2,by=3)){
  data1$d[i] = data1$d[i+1]
  data1$d[i+2] = data1$d[i+1]
  data1$delta[i] = data1$delta[i+1]
  data1$delta[i+2] = data1$delta[i+1]
}
data1$d<-as.factor(data1$d)
data1$delta<-as.factor(data1$delta)
levels(data1$d)<-c(d02)
levels(data1$delta)<-c(del02)
data1<-data1[,which(names(data1)%in%c("I","d","delta","model","CNPP","persistence","PP"))]

databis<-data[data$d==as.character(d02) & data$delta==as.character(del02),
                    which(names(data)%in%c("I","d","delta","model","persistence","PPprod"))]
databis$CNPP<-6.6
names(databis)[names(databis)=="PPprod"]="PP"
data1<-rbind(data1,databis)

data1$CNPP<-as.factor(data1$CNPP)
levels(data1$CNPP)<-c('C:N=6.6',"C:N=8","C:N=11")

databis <- summarySE(data1, measurevar="persistence", groupvars=c("I","d","delta","model","CNPP"),na.rm=TRUE)
p1<-ggplot()+
      geom_line(data=databis,aes(I,persistence,colour=model,linetype=model), size=2)+
      geom_errorbar(data=databis,aes(I,ymin=persistence-ci, ymax=persistence+ci,colour=model), width=15, show.legend=FALSE)+
      facet_wrap(~CNPP)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5, "cm"))+
      lim_short+
      xlab(label_nutri)+
      ylab("Species persistence")

databis <- summarySE(data1[data1$persistence>0,], measurevar="PP", groupvars=c("I","d","delta","model","CNPP"),na.rm=TRUE)
p2<-ggplot()+
      geom_line(data=databis,aes(I,PP,colour=model,linetype=model), size=2)+
      geom_errorbar(data=databis,aes(I,ymin=PP-ci, ymax=PP+ci,colour=model), width=15, show.legend=FALSE)+
      facet_wrap(~CNPP)+
      model_colour+
      model_line+
      theme+theme(legend.key.width = unit(1.5, "cm"))+
      lim_short+
      ylim(0,8500)+
      xlab(label_nutri)+
      ylab("Primary production")

graph<-plot_grid(p1, p2,
                 labels = c("A","B"), label_size = 25,
                 nrow = 2, align = "vh")
ggsave("supp_CNPP.pdf",graph, width = 12, height = 8, device = cairo_pdf)

# FIGURES FINALES - OTHER -----------------------------------------------------------------
#### Handling time ####
handlingtime<-function(b,y,Mi,Mj){
  return(Mj/Mi*b^2/(6*y*Mj^(-0.25)*(b-Mi/Mj)))
}

b=0.05
y=8*0.27
Mj=1000 # predator body mass
Mi<-seq(5,45,by=0.01) # prey body mass

htime<-handlingtime(b,y,Mi,Mj)

data1<-data.frame(Mi,htime)

p1<-ggplot(data=data1)+
  geom_line(aes(Mi,htime),size=2)+
  geom_vline(xintercept=b/2*Mj,color="green",size=2,linetype="dashed")+
  geom_vline(xintercept=b*Mj,color="red",size=2,linetype="dashed")+
  theme+
  xlab("Prey body mass (kg)")+
  ylab("Handling time (year/kg of prey)")
ggsave("supp_handling_time.pdf",p1, width = 7, height = 5, device = cairo_pdf)
