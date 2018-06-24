##plot-univariate.R
library(generalizedPELT)
lib.for.me()

runs.D.now=20

## Generate a data set of correct proportions
simD <-list()

temp_delta=seq(from=0 , to=5,by=1)



simD$delta = temp_delta[2:(length(temp_delta))]
simD$delta.len=length(simD$delta)
# Length of dataset
tau_1=20
m=4
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
}


####
## Proportion with each m for different deltas, BS and PELT
#### (C1)
p09<-barplot2.double(data=simD$print$both[simD$print$both$Freq>0.005,],
                     name = paste0(ifelse(simD$even.not.unif,
                                          "even-","unif-"),"Freq of different m found, m=",simD$m,", >0.005"),
                     runs = runs.D.now)
p09<-p09+ theme(legend.key.height=unit(0.75,"line"))
p09+ggtitle("")























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

temp3=temp2



####
## Proportion with 1 changepoints (m is correct)
#### (C2)


if(simD$even.not.unif){
  p03.5<-hitprop.plot.v01(data=temp3,
                      name = "Proportion of simulations classifying to one changepoint",
                      runs = runs.D.now, stip=FALSE)
  p03.5


}else{
  temp4=temp3[temp3$Type=="m.is.ok",]
  p03.6<-hitprop.plot.v02(data=temp4,
                          name = "Proportion of simulations classifying to one changepoint",
                          runs = runs.D.now, stip=FALSE)
  p03.6
}




