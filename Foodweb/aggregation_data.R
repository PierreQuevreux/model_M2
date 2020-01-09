path="results/"
data<-NULL
bodymass<-NULL
biomass<-NULL
biomassCV<-NULL
recy<-NULL
recyCV<-NULL
TL<-NULL
file<-NULL

nSlice=1

for (i in 0:(nSlice-1)){
  file<-read.table(paste(path,"data_",i,".txt",sep=""),sep=';',header=T)
  data<-rbind(data,file)
  file<-read.table(paste(path,"bodymass_",i,".txt",sep=""),sep=';',header=T)
  bodymass<-rbind(bodymass,file)
  file<-read.table(paste(path,"biomass_",i,".txt",sep=""),sep=';',header=T)
  biomass<-rbind(biomass,file)
  file<-read.table(paste(path,"biomassCV_",i,".txt",sep=""),sep=';',header=T)
  biomassCV<-rbind(biomassCV,file)
  file<-read.table(paste(path,"recy_",i,".txt",sep=""),sep=';',header=T)
  recy<-rbind(recy,file)
  file<-read.table(paste(path,"recyCV_",i,".txt",sep=""),sep=';',header=T)
  recyCV<-rbind(recyCV,file)
  file<-read.table(paste(path,"TL_",i,".txt",sep=""),sep=';',header=T)
  TL<-rbind(TL,file)
}
write.table(data,paste(path,"data.txt",sep=""),sep=';',row.names=F)
write.table(bodymass,paste(path,"bodymass.txt",sep=""),sep=';',row.names=F)
write.table(biomass,paste(path,"biomass.txt",sep=""),sep=';',row.names=F)
write.table(biomassCV,paste(path,"biomassCV.txt",sep=""),sep=';',row.names=F)
write.table(recy,paste(path,"recy.txt",sep=""),sep=';',row.names=F)
write.table(recyCV,paste(path,"recyCV.txt",sep=""),sep=';',row.names=F)
write.table(TL,paste(path,"TL.txt",sep=""),sep=';',row.names=F)
