path="data/"

nSlice=1
#IDsimu=0:(nSlice-1)
IDsimu=2
#I=100
I=seq(from=1,to=400,by=3)
delta=c(0.2,0.8)
d=c(0.2,0.8)

######
# Crossed variables

params<-expand.grid(simu_ID=0,
                    I=I,
                    delta=delta,
                    d=d)
nParams=dim(params)[2]

n=dim(params)[1]
params$simu_ID<-seq(1,n)

# split the variables into nSlice sub tables
nSimu=seq(1:nSlice)
nSimu[1:nSlice]=n %/% nSlice
nSimu[nSlice]=nSimu[nSlice]+n-nSimu[nSlice]*nSlice

######
# Save the parameters

start=1
for (i in 1:nSlice){
  sub_params<-params[start:(start+nSimu[i]-1),]
  params_data<-c(nSimu[i],nParams)
  write.table(sub_params,paste(path,"parameters_",IDsimu[i],".txt",sep=""),sep=";",row.names = FALSE)
  write.table(params_data,paste(path,"parameters_data_",IDsimu[i],".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)
  start=start+nSimu[i]
}
