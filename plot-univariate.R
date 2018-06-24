##q05v18
# Make the discovered m -plot
# Make the exact and correct m -plot
# With 1dmean, mBIC and BIC
# With uniform distributed changepoint placement

lib.for.me()
#set.window(FALSE)
#m=5

filename="q05v18"

#sessionNUMBER=14
setwd("C:/Users/Kristin/Documents/v18-master/pelt-package-v01/project-imported-v01")
runs.D.now=20

## Generate a data set of correct proportions
simD <-list()
#temp=seq(from=0,to=100,by=10)
# Different to-values for mu_1
#temp_delta=seq(from=0.5,to=5.5,by=0.5)
#temp_delta=sort(c(seq(from=0.5,to=2,by=0.5),2.25,seq(from=2.5,to=6,by=0.5),0,1.25,1.6,1.7,1.8,1.90))
temp_delta=seq(from=0 , to=5,by=1)
#temp_delta=sort(c(1.5,2,seq(from=3.5,to=6,by=0.5),2.40,2.80,3.20))
#temp_delta=c(seq(from=0.35,to=1.75,0.2),2,seq(from=2.5,to=7,by=0.5))
#temp_delta=c(seq(from=0.35,to=1.75,0.2),2,seq(from=2.5,to=6.5,by=0.5),seq(from=7,to=10,by=1))
#temp_delta=seq(from=0,to=20,by=2)
#temp_delta=c(0,0.5,1,5,20)
#temp_delta=c(0,5)
#temp_delta=seq(from=0,to=4.5,by=.5)
#temp_delta=c(0,4.5,5,5.5)
#temp_delta=seq(from=0,to=1.125,by=0.25)


simD$delta = temp_delta[2:(length(temp_delta))]
simD$delta.len=length(simD$delta)
# Length of dataset
tau_1=20
m=40
n=(m+1)*tau_1
max.pen.rel(n=n,m=m)

simD$m=m
simD$length = (m+1)*tau_1
simD$length
# Placement of changepoint
close.not.even=FALSE
simD$even.not.unif=FALSE
if(close.not.even){
  simD$cpt = (1:m)
}else{
  simD$cpt = (1:m)*tau_1
}

# What is counted as an exact hit
simD$precision = 0 # 0 or
simD$hit=simD$cpt + seq(from=-simD$precision,to=simD$precision,by=1)
# Number of simulations
simD$sims=runs.D.now

# Simulates one dataset for each mean
# # datasets= length(simD$delta)*sims
# # points  = length(simD$delta)*sims*simD$length

# m simulation for each delta, simulate data
if(simD$even.not.unif){
simD$data= lapply(simD$delta,
                  function(x) simulate.mymean(mean=(0:m)*x,cpt=c(0,simD$cpt,
                                    simD$length),sims=simD$sims,seed=-1,sd=1))
}else{
simD$data= lapply(simD$delta,
                  function(x) simulate.mymean.uniform(mean=(0:m)*x, m=simD$m,
                                    n=simD$length,sims=simD$sims,seed=-1,sd=1))
}

simD$data.gpelt=lapply(simD$data,function(x) as.data.frame(t(x)))

## This is where the 'Magic' (=analysis) happens :)
simD$bic1=lapply(1:simD$delta.len,
             function(i) analysis.hitprop.w.cpts2(data=simD$data[[i]],sims=simD$sims,
                                                  pelt.yes=TRUE,hit=simD$hit,delta=simD$delta[i],
                                                  true.m=m,cpts=simD$cpt,n=simD$length))

attrb2=attrb.mycpt(p=1,n=simD$length,pen=3*log(simD$length),minseglen=1)
thecptpackage=FALSE
if(thecptpackage){
  simD$mbic=lapply(1:simD$delta.len,
                   function(i) analysis.hitprop.w.cpts2(data=simD$data[[i]],sims=simD$sims,
                                                        pelt.yes=TRUE,hit=simD$hit,delta=simD$delta[i],
                                                        true.m=m,cpts=simD$cpt,n=simD$length,
                                                        mBIC=TRUE))
}else{
  simD$mbic=lapply(1:simD$delta.len,
             function(i) analysis.hitprop.w.cpts4(data=simD$data[[i]],sims=simD$sims,
                                                  pelt.yes=TRUE,delta=simD$delta[i],
                                                  true.m=m,cpts=simD$cpt,n=simD$length,
                                                  data.gpelt=simD$data.gpelt[[i]],attrb=attrb2,type="mbic.1d.mean"))
}
beep()
##
#sessionNUMBER=25
folder=paste0("proj-plot/session",sessionNUMBER)
setwd(paste0("C:/Users/Kristin/Documents/v18-master/pelt-package-v01/",folder))


##


## Plot Frequency of number of changepoints with PELT as bar-chart
###---This is plot sim2 for PELT---
#printdeltas=c(1,2,8,14,18)

if(TRUE){
printdeltas=1:length(simD$delta)
#printdeltas=1:12

simD$print=list()
# Extract m of each
item2=cbind(lapply(printdeltas,function(i) simD$bic1[[i]]$m))
# Merge all the simD$pelt[[i]]$m into one data frame
item4=Reduce(function(...) merge(...,all=TRUE),item2)
simD$print$bic1=item4
item4$Method="BIC 1"

item5.1=cbind(lapply(printdeltas,function(i) simD$mbic[[i]]$m))
item5=Reduce(function(...) merge(...,all=TRUE),item5.1)
simD$print$mbic=item5
item5$Method="mBIC"


simD$print$both = rbind(item4,item5)

if(FALSE){
  ####
  ## Proportion with each m for different deltas, PELT
  #### (C1)

  p08.1<-barplot2.single(data=simD$print$bic1[simD$print$bic1$Freq>0.005,],
                       name = "Freq of different m found, m=5,BIC 1",
                       runs = runs.D.now)
  p08.1
  mysave.proj(p08.1,save=TRUE,filenam=filename,m=5,runs=runs.D.now,
              variety=paste0("C1-2-PELT-tau1-",simD$cpt[1],"-dmax99"),folder=paste0("proj-plot/session",sessionNUMBER))


  p08.2<-barplot2.single(data=simD$print$mbic[simD$print$mbic$Freq>0.005,],
                         name = "Freq of different m found, m=5,BS",
                         runs = runs.D.now)
  p08.2
  mysave.proj(p08.2,save=TRUE,filenam=filename,m=5,runs=runs.D.now,
              variety=paste0("C1-2-BS-tau1-",simD$cpt[1],"-dmax99"),folder=paste0("proj-plot/session",sessionNUMBER))
}
}

#load.data.both=simD$print$both[simD$print$both$Freq>0.005,]
#load2.data.both=simD$print$both[simD$print$both$Freq>0.005,]

#both=rbind(load.data.both,load2.data.both)

####
## Proportion with each m for different deltas, BS and PELT
#### (C1)
p09<-barplot2.double(data=simD$print$both[simD$print$both$Freq>0.005,],
                     name = paste0(ifelse(simD$even.not.unif,
                                          "even-","unif-"),"Freq of different m found, m=",simD$m,", >0.005"),
                     runs = runs.D.now)
p09<-p09+ theme(legend.key.height=unit(0.75,"line"))
p09+ggtitle("")

#kjort 20:21. TIl 20:22

mysave.proj(p09,save=TRUE,filenam=filename,m=simD$m,runs=runs.D.now,
            variety=paste0(ifelse(thecptpackage,"thecptpackage-",""),ifelse(simD$even.not.unif,ifelse(close.not.even,"close","even"),"unif-"),
                           "-C1-2-both-tau",paste0(simD$cpt,collapse = "-"),"-dmax3"),folder=paste0("proj-plot/session",sessionNUMBER))





















if(TRUE){
# Graph plot
## Plot proportion that found correct number of changepoints
# (1) Make df of each set of data
temp1.bic1=data.frame(
  exact=sapply(1:simD$delta.len,function(i) simD$bic1[[i]]$hit.prop),
  m.is.ok =sapply(1:simD$delta.len,function(i) simD$bic1[[i]]$ncpts.correct),
  Method=rep("BIC 1",simD$delta.len), stringsAsFactors=FALSE
)

temp1.mbic=data.frame(
  exact=sapply(1:simD$delta.len,function(i) simD$mbic[[i]]$hit.prop),
  m.is.ok =sapply(1:simD$delta.len,function(i) simD$mbic[[i]]$ncpts.correct),
  Method=rep("mBIC",simD$delta.len), stringsAsFactors=FALSE
)

# (2) Melt each, and put together

temp2.mbic= melt(temp1.mbic, id="Method", stringsAsFactors=FALSE)
temp2.bic1= melt(temp1.bic1, id="Method", stringsAsFactors=FALSE)

if(TRUE){
  i <- sapply(temp2.mbic, is.factor)
  temp2.mbic[i] <- lapply(temp2.mbic[i], as.character)

  i <- sapply(temp2.bic1, is.factor)
  temp2.bic1[i] <- lapply(temp2.bic1[i], as.character)
}


#temp2.mbic$delta=simD$delta
#temp2.bic1$delta=simD$delta

temp2= rbind(temp2.mbic,temp2.bic1)

temp2$delta=simD$delta

names(temp2)
names(temp2)=c("Method","Type", "Value","delta")
##shorten:
#delta.subset=simD$delta[1:11]
#delta.subset=simD$delta


#temp3=temp2[is.element(temp2$delta,delta.subset),]
temp3=temp2

if(FALSE){

if(FALSE){
  head(temp2)
  simD$plot.correct.ncpts1=data.frame( delta=simD$delta,
                                       BS.c=sapply(1:simD$delta.len,function(i) simD$mbic[[i]]$hit.prop),
                                       BS.m=sapply(1:simD$delta.len,function(i) simD$mbic[[i]]$ncpts.correct),
                                       PELT.c=sapply(1:simD$delta.len,function(i) simD$bic1[[i]]$hit.prop),
                                       PELT.m=sapply(1:simD$delta.len,function(i) simD$bic1[[i]]$ncpts.correct))
  simD$plot.correct.ncpts<-melt(simD$plot.correct.ncpts1,id="delta")
}


####
## Proportion with 1 changepoints (m is correct)
#### (C2)
if(FALSE){
  p03<-hitprop.plot(data=temp3,
                       name = "Proportion of simulations classifying to one changepoint",
                       runs = runs.D.now, stip=TRUE)
  p03
}
  #mysave.proj(p03,save=TRUE,filenam="q05v15",m=5,runs=runs.D.now,
#            variety="C2-both-hitprop-stipled",folder=paste0("proj-plot/session",sessionNUMBER))

if(FALSE){
  p03.2<-hitprop.plot(data=temp3,
                    name = "Proportion of simulations classifying to one changepoint",
                    runs = runs.D.now, stip=FALSE)
  p03.2
  mysave.proj(p03.2,save=TRUE,filenam="q05v15",m=5,runs=runs.D.now,
              variety=paste0("C2-both-hitprop-tau1-",simD$cpt[1],"-n",simD$length,"-dmax4"),folder=paste0("proj-plot/session",sessionNUMBER))
}


#save.my.data=temp3
#save.my.other.data=plotit

if(FALSE){
  plotit=temp3

  plotit=rbind(save.my.data[save.my.data$delta>1.25,],save.my.other.data,temp3)
  plotit$Method[plotit$Method=="BS"]="mBIC"
}
}#FALSE
}#TRUE

#load.temp3=temp3
#load2.temp3=temp3

#temp3=rbind(load.temp3,load2.temp3)

if(simD$even.not.unif){
  p03.5<-hitprop.plot.v01(data=temp3,
                      name = "Proportion of simulations classifying to one changepoint",
                      runs = runs.D.now, stip=FALSE)
  p03.5
  mysave.proj(p03.5,save=TRUE,filenam="q05v17",m=5,runs=runs.D.now,
              variety=paste0(ifelse(thecptpackage,"thecptpackage",""),ifelse(close.not.even,"-close","-even"),
                             "-C2-both-hitprop-tau",paste0(simD$cpt,collapse = "-"),"-n",simD$length,"-dmax4-lots"),
              folder=paste0("proj-plot/session",sessionNUMBER))

}else{
  temp4=temp3[temp3$Type=="m.is.ok",]
  p03.6<-hitprop.plot.v02(data=temp4,
                          name = "Proportion of simulations classifying to one changepoint",
                          runs = runs.D.now, stip=FALSE)
  p03.6
  mysave.proj(p03.6,save=TRUE,filenam="q05v17",m=5,runs=runs.D.now,
              variety=paste0("-unif-m-only-C2-both-hitprop-tau",paste0(simD$cpt,collapse = "-"),"-n",simD$length,"-dmax4-lots"),folder=paste0("proj-plot/session",sessionNUMBER))
}




