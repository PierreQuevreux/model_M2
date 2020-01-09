bifurcation<-NULL
time_series<-NULL
mean<-NULL
CV<-NULL
parameters<-NULL
file<-NULL

nSlice=4

for (i in 0:(nSlice-1)){
  file<-read.table(paste("data/bifurcation_",i,".txt",sep=""),sep=';',header=T)
  bifurcation<-rbind(bifurcation,file)
  #file<-read.table(paste("data/time_series_",i,".txt",sep=""),sep=';',header=T)
  #time_series<-rbind(time_series,file)
  file<-read.table(paste("data/mean_",i,".txt",sep=""),sep=';',header=T)
  mean<-rbind(mean,file)
  file<-read.table(paste("data/CV_",i,".txt",sep=""),sep=';',header=T)
  CV<-rbind(CV,file)
  file<-read.table(paste("data/parameters_",i,".txt",sep=""),sep=';',header=T)
  parameters<-rbind(parameters,file)
}
write.table(bifurcation,"data/bifurcation.txt",sep=';',row.names=F)
#write.table(time_series,"data/time_series.txt",sep=';',row.names=F)
write.table(mean,"data/mean.txt",sep=';',row.names=F)
write.table(CV,"data/CV.txt",sep=';',row.names=F)
write.table(parameters,"data/parameters.txt",sep=';',row.names=F)