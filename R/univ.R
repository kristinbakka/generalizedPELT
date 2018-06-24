## qo1.R
# Only the functions and stuff

## Load libraries
#library('zoo') # for spcadjust
#library('spcadjust') # for cusum
lib.for.me <- function(){
  library(changepoint) # for cpt.mean(), cpt.meanvar()
  library(ggplot2)
  library(wesanderson)
  library(reshape2)
  library(stringr)
  library(beepr)
}

if(FALSE){
  ## For saving and naming of plots. Change this point only to get save
  save.yes=FALSE
  mycount=1
  version="2"
  sims.global=1000

  if(save.yes){
    dev.off(); # dev.new(width=6, height=3)
    windows(record=TRUE, width=5, height=2.5) #w=5,h=2.5 blir bra med width=0.8\textwidth i Latex (kanskje litt h?y figur)
  }
}

do.plot.save <- function(plot,version){
  name=paste(c("v", version,"plot",str_pad(mycount, 3, pad = "0")), collapse = "")
  mycount=mycount+1
  if(save.yes){
    ggsave(plot=plot+labs(title=""),  paste(c(name,".pdf"), collapse = ""), device="pdf") # Saves last plot
  }
  plot + labs(subtitle=paste(c("v", version,"plot",str_pad(mycount, 3, pad = "0")), collapse = ""))
  return(mycount)
}

# Plot the timeseries
plot.ts.different.delta <-function(simI,printdeltas=0,old.version=TRUE){
  #input: simI=list(data, delta, delta.len)
  ## Plot an example of each type of simulation
  # Multiple changepoints, differing deltas
  simI$plot=list()
  if(identical(printdeltas,0)){
    # 1-First get data such that each column is a simulation
    simI$plot$data.mat=as.data.frame(sapply(1:simI$delta.len,function(i) simI$data[[i]][1,]))
    # 2-Give the column names for each delta (str version of delta)
    names(simI$plot$data.mat)=sapply(simI$delta,toString)
  }else{
    # Else, when subset is selected
    mydeltas=simI$delta[printdeltas]
    # 1-First get data such that each column is a simulation
    simI$plot$data.mat=as.data.frame(sapply(printdeltas,function(j) simI$data[[j]][1,]))
    # 2-Give the column names for each delta (str version of delta)
    names(simI$plot$data.mat)=sapply(mydeltas,toString)
  }
  # 3-One column is t=
  simI$plot$data.mat$t=1:length(simI$plot$data.mat[,1])
  # 4-Melt to plottable long-format
  simI$plot$data.mat=melt(simI$plot$data.mat,id="t")

  #simI$data[[1]][2,]
  #plot.df<-data.frame(x=1:length(simI$data[[1]][3,]),y=simI$data[[1]][2,] )
  if(old.version){
    return(simI$plot)
  }
  # Don't need to return everything
  return(simI$plot$data.mat)
}


## Simulate a data set using this function
## sigma=1
simulate.mymean <- function(mean=c(10),cpt=c(0,100),sims=1,seed=-1,sd=1){
  # For simulating a matrix, returns only data
  m=length(mean)
  if(m!=(length(cpt)-1)){
    return(0)
  }
  # Last cpt is length of data set
  len=cpt[m+1]
  # Length of each interval
  times=cpt[2:(m+1)]-cpt[1:m]
  # Each row is a sample
  if(seed!=-1){
    set.seed(seed)
  }
  if(length(sd)==1){
    # Remark that cpt[...] is vector, so times repeats each element times times
    return(t(replicate(sims,rnorm(len,mean=rep(mean,times=times),sd=sd))))
  }else{
    # Remark that cpt[...] is vector, so times repeats each element times times
    return(t(replicate(sims,rnorm(len,mean=rep(mean,times=times),
                                  sd=rep(sd,times=times)))))
  }
}


simulate.mymean.uniform <- function(mean=c(0,10),m,n,sims=1,seed=-1,sd=1){
  # For simulating a matrix, returns only data

  if(m!=(length(mean)-1)){
    return(0)
  }

  #
  cpt = t(replicate(sims,c(0,sort(sample.int(n-1,m,replace=FALSE)),n)))
  # Length of each interval
  times=t(apply(cpt,1,function(x) x[2:(m+2)]-x[1:(m+1)]))
  # Each row is a sample
  if(seed!=-1){
    set.seed(seed)
  }
  if(length(sd)==1){
    # Remark that cpt[...] is vector, so times repeats each element times times
    return(t(sapply(1:sims,function(i) rnorm(n,mean=rep(mean,times=times[i,]),
                                             sd=rep(sd,times=n)))))
  }else{
    # Remark that cpt[...] is vector, so times repeats each element times times
    warning("I haven't tested this choice (different sd values) properly yet. <3")
    return(t(sapply(1:sims,function(i) rnorm(n,mean=rep(mean,times=times[i,]),
                                             sd=rep(sd,times=times[i,])))))
  }
}


# This is how it works
# 6 data points, point 3 is a changepoint, 10 such data sets
## sigma=1
simulate.mymean(mean=c(0,10),cpt=c(0,3,6),sims=10,seed=500)



# Simulations Old version             file
# simA  - 0 cpts, 1 sim               q01.R
# simA2 - 0 cpts, 1 sim               q01.R
# simA3 - 0 cpts  1 sim               q01.R
# simB  - 0 cpts, different n         q02.R
# simSD - 4 cpts, 1 sim               q03.R
# simC  - 1 cpt,  1 sim               q04.R

# simD  - 1 cpt, different deltas     q05.R - deleated piechart
# simE  - 1 cpt, different n          q06.R
# simF  -
# simG  -
# simH  - 5 cpt, just 1 sim           q07.R
# simI  - 5 cpts different deltas     q08.R
# simK  - 5 cpts different deltas     q09.R
          # better values maybe.


# Simulations
# file       cpts  Simulations new version
# q02v02.R    0    simB - hit proportion, ts plot




# simA functions
cpt.eval.ncpts <- function(simA){
  cat('\tMaximal number of detected internal changepoints is
      max(hat(m))=',max(simA$PELT$ncpts),' and minimum is
      min(hat(m))=',min(simA$PELT$ncpts), ' with PELT on',simA$sims,
      'datasets of length ',simA$length,'.')

  cat('\n\tMaximal number of detected internal changepoints is
      max(hat(m))=',max(simA$BS$ncpts),' and minimum is
      min(hat(m))=',min(simA$BS$ncpts), ' with BS on',simA$sims,
      'datasets of length ',simA$length,'.')

  cat('\n\tThe true changepoints are', simA$cpts,',
      so true m = ',length(simA$cpts)-2)

  cat('\n\tProportion of simulations that give estimated hat(m)=0 is
      With   BS: ', length(which(simA$BS$ncpts==0))/simA$sims,'
      With PELT: ', length(which(simA$BS$ncpts==0))/simA$sims)
}

analysis <- function(data,sims,pelt.yes=TRUE,Q=5){
  # This gives list and not vector of vector cpts
  ans = list()
  if(pelt.yes){
    run=cpt.mean(data=data,penalty="BIC",
                 param.estimates=FALSE,method="PELT")
  }else{
    run=cpt.mean(data=data,penalty="BIC", Q=Q,
                 param.estimates=FALSE,method="BinSeg")
  }
  # Save  number  of changepoints
  ans$ncpts= sapply(1:sims,function(i) ncpts(run[[i]]))
  # Save position of changepoints
  ans$cpts = lapply(1:sims,function(i) run[[i]]@cpts)
  return(ans)
}

analysis.hitprop.ARL <- function(data,sims,pelt.yes=TRUE,Q=5){
  # Proportion of sims evaluated to <<no changepoints>>
  # Only need data to return hit probabilities
  ans = list()
  if(pelt.yes){
    run=cpt.mean(data=data,penalty="BIC",
                 param.estimates=FALSE,method="PELT")
  }else{
    run=cpt.mean(data=data,penalty="BIC", Q=Q,
                 param.estimates=FALSE,method="BinSeg")
  }
  # Save  number  of changepoints
  hit= sapply(1:sims,function(i) ifelse(ncpts(run[[i]])==0,1,0))
  return(sum(hit)/sims)
}

## Measure time
measure.time <- function(times){
  a<-Sys.time()
  return(c(times,a[[1]]))
}
#how.long.time <- function(times){
#  return(time[2:length(time)]-time[1:(length(time)-1)])
#}
how.long.time <- function(times){
  ret=rep(NA,length(time))
  for (i in 2: (length(time)-1)){
    ret[i]=time[i+1][[1]]-time[i][[1]]
  }
  return(ret)
}



## simC functions
analysis.hitprop.w.cpts4 <- function(data,sims,pelt.yes=TRUE,Q=5,delta=0,
                                     true.m=1,cpts=10,n=20,data.gpelt,attrb,type){
  if(!is.element(type,c("1d.mean","1d.meanvar","1d.var",
                        "mbic.1d.mean","mbic.1d.meanvar","mbic.1d.var"))){
      warning("This only works with one parameter. :)")
  }
  # When there is in truth 1 changepoint,
  # which proportion finds the correct

  #choose whether pelt or BinSe solution

  # Q is maximum number of changepoints
  # Only works for one changepoint
  ans = list()

  # applies PELT to an entire runset (that is one configuration of parameters)
  if(pelt.yes){
    #run2=cpt.mean(data=data,penalty="BIC",
    #             class=FALSE,method="PELT")
    run.g=gpelt.both.mycpt(data.gpelt,attrb=attrb,type=type,both=FALSE)
  }else{
    warning("Only meant for pelt.yes=TRUE. :D")
    #run=cpt.mean(data=data,penalty="BIC", Q=Q,
     #            param.estimates=FALSE,method="BinSeg")
    run.g=gpelt.both.mycpt(data.gpelt,attrb=attrb,type=type,both=FALSE)
  }

  # Save  number  of changepoints
  ncpts= sapply(1:sims,function(i) length(run.g[[i]])-1)


  # Proportion with each number of changepoints
  ans$ncpts=as.data.frame(table(ncpts))
  ans$ncpts$Freq=ans$ncpts$Freq/sims

  # Also include those that are different
  different=setdiff(0:Q,ans$ncpts$ncpts)
  df=data.frame(ncpts=different,Freq=rep(0,length(different)))
  ans$m=rbind(df,ans$ncpts)
  ans$m=ans$m[order(ans$m$ncpts),]
  ans$m


  # Proportion with correct number of changepoints
  ans$ncpts.correct=ans$m$Freq[ans$m$ncpts==true.m]

  # Add column identifying this sequence
  if(delta!=0){
    ans$ncpts$delta=rep(x=delta,times=length(ans$ncpts$ncpts))
    ans$m$delta    =rep(x=delta,times=length(ans$m$ncpts))
  }

  # Save position of changepoints
  #ans$cpts = lapply(1:sims,function(i) run2[[i]])

  # Value of subjects with one changepoint, and a changepoint inside hit
  ans$hit= sapply(1:sims,function(i) ifelse(identical(as.numeric(run.g[[i]]),as.numeric(c(cpts,n))),as.integer(run.g[[i]]),NA))
  # Strip NA's
  ans$hit <- ans$hit[!is.na(ans$hit)]

  if(!length(ans$hit)){
    ans$hit=data.frame(Value=1,Proportion=0)
    ans$hit.prop=0
  }else{
    ### This is only done if there is some hits
    # Proportion of simulations with hit
    ans$hit.prop=length(ans$hit)/sims
    if(delta!=0){
      ans$correct = data.frame(correct=ans$hit.prop,
                               delta=rep(x=delta,times=length(ans$ncpts$ncpts)))
    }

    if(FALSE){
      # Proportion with each hit value
      # want each hit value represented, so add them first and subtract later
      # add:
      ans$hit<- as.data.frame(table(c(ans$hit,hit)))
      names(ans$hit)=c("Value","Proportion")
      # subtract
      ans$hit$Proportion=ans$hit$Proportion-rep(x=1,times=length(ans$hit$Proportion))
      ans$hit$Proportion=ans$hit$Proportion/sims

      # Add column identifying this sequence
      if(delta!=0){
        ans$hit$Delta=rep(x=delta,times=length(ans$hit$Proportion))
      }
    }
  }

  return(ans)
}

analysis.hitprop.w.cpts2 <- function(data,sims,hit=c(1),pelt.yes=TRUE,Q=5,delta=0,true.m=1,cpts=10,n=20,mBIC=FALSE){
  # When there is in truth 1 changepoint,
  # which proportion finds the correct

  #choose whether pelt or BinSe solution

  # Q is maximum number of changepoints
  # Only works for one changepoint
  ans = list()

  # applies PELT to an entire runset (that is one configuration of parameters)
  if(pelt.yes){
    if(mBIC){
      run2=cpt.mean(data=data,penalty="MBIC",
                    class=FALSE,method="PELT")
    }else{
      run2=cpt.mean(data=data,penalty="BIC",
                  class=FALSE,method="PELT")
    }
  }else{
    run2=cpt.mean(data=data,penalty="BIC", Q=Q,
                 class=FALSE,param.estimates=FALSE,method="BinSeg")
  }

  # Save  number  of changepoints
  ncpts= sapply(1:sims,function(i) length(run2[[i]])-1)


  # Proportion with each number of changepoints
  ans$ncpts=as.data.frame(table(ncpts))
  ans$ncpts$Freq=ans$ncpts$Freq/sims

  # Also include those that are different
  different=setdiff(0:Q,ans$ncpts$ncpts)
  df=data.frame(ncpts=different,Freq=rep(0,length(different)))
  ans$m=rbind(df,ans$ncpts)
  ans$m=ans$m[order(ans$m$ncpts),]
  ans$m


  # Proportion with correct number of changepoints
  ans$ncpts.correct=ans$m$Freq[ans$m$ncpts==true.m]

  # Add column identifying this sequence
  if(delta!=0){
    ans$ncpts$delta=rep(x=delta,times=length(ans$ncpts$ncpts))
    ans$m$delta    =rep(x=delta,times=length(ans$m$ncpts))
  }

  # Save position of changepoints
  #ans$cpts = lapply(1:sims,function(i) run2[[i]])

  # Value of subjects with one changepoint, and a changepoint inside hit
  ans$hit= sapply(1:sims,function(i) ifelse(identical(as.numeric(run2[[i]]),as.numeric(c(cpts,n))),run2[[i]],NA))
  # Strip NA's
  ans$hit <- ans$hit[!is.na(ans$hit)]

  if(!length(ans$hit)){
    ans$hit=data.frame(Value=1,Proportion=0)
    ans$hit.prop=0
  }else{
    ### This is only done if there is some hits
    # Proportion of simulations with hit
    ans$hit.prop=length(ans$hit)/sims
    if(delta!=0){
      ans$correct = data.frame(correct=ans$hit.prop,
                               delta=rep(x=delta,times=length(ans$ncpts$ncpts)))
    }

    # Proportion with each hit value
    # want each hit value represented, so add them first and subtract later
    # add:
    ans$hit<- as.data.frame(table(c(ans$hit,hit)))
    names(ans$hit)=c("Value","Proportion")
    # subtract
    ans$hit$Proportion=ans$hit$Proportion-rep(x=1,times=length(ans$hit$Proportion))
    ans$hit$Proportion=ans$hit$Proportion/sims

    # Add column identifying this sequence
    if(delta!=0){
      ans$hit$Delta=rep(x=delta,times=length(ans$hit$Proportion))
    }
  }

  return(ans)
}

# simH
simulate.multiple.means <- function(delta=2,cpt.pos=50,m=5,sims=1,seed=-1,sd=1,return.data.only=TRUE,len.max=0){
  # First changepoint at cpt.pos, each time mean increases with delta
  # There will be m internal changepoints, each segment will have cpt.pos points

  # Total number of points will be (m+1)(cpt.pos)
  len.out=(m+1)*(cpt.pos)
  means=(0:m)*delta

  # In the last part when varying lengths
  if(len.max==0){
    len.max=len.out
  }

  # Each row is a sample
  if(seed!=-1){
    set.seed(seed)
  }

  if(return.data.only){
    # Remark that cpt.pos is integer, so each repeats each element
    return(matrix(c(t(replicate(sims,rnorm(len.out,mean=rep(means,each=cpt.pos),1))),
                    matrix(NA,nrow=sims,ncol = len.max-len.out)),nrow=sims,ncol=len.max))
  }else{
    ans=list(delta=delta,m=m,cpt.pos=cpt.pos,sims=sims,data=t(replicate(sims,rnorm(len.out,mean=rep(means,each=cpt.pos),1))),Q=floor(len.out/2)+1,len.out=len.out)
    return(ans)
  }
}
# This is how it works
# 8=2*(3+1) data points, point 2,4,6 are internal changepoints, 10 such data sets
simulate.multiple.means(delta=100,cpt.pos = 2,m=3,sims=5,seed=500)
## Undress it
## THE ANALYSIS
analysis.multiple.means <- function(data,sims,pelt.yes=TRUE,Q,delta,m,cpt.pos){
  ans = list()
  if(pelt.yes){
    run=cpt.mean(data=data,penalty="BIC",
                 param.estimates=FALSE,method="PELT")
  }else{
    run=cpt.mean(data=data,penalty="BIC", Q=Q,
                 param.estimates=FALSE,method="BinSeg")
  }

  # Save  number  of changepoints
  ncpts= sapply(1:sims,function(i) ncpts(run[[i]]))
  ans$ncpts=as.data.frame(table(ncpts))
  ans$ncpts$Freq=ans$ncpts$Freq/sims

  # Add column identifying this sequence
  ans$ncpts$delta=rep(x=delta,times=length(ans$ncpts$ncpts))

  # Proportion with correct number of changepoints
  ans$correct.ncpts= sum(ifelse(ncpts==m,1,0))/sims

  ## Check if exactly correct changepoints are returned
  # True position of all changepoints
  true.cpts = (1:(m+1))*cpt.pos

  # 1 equal, 0 unequal, NA different length
  correct.pos= sapply(1:sims,function(i) ifelse(all.equal(run[[i]]@cpts,true.cpts),1,0))
  ans$correct.pos <- sum(correct.pos[!is.na(correct.pos)])/sims

  #ans$correct.pos= sum(sapply(1:sims,function(i) ifelse(all.equal(run[[i]]@cpts,true.cpts),1,0)))/sims
  ## Other for error searching
  # True position of all changepoints
  ans$true.cpts=true.cpts
  # Save position of changepoints
  ans$cpts = lapply(1:sims,function(i) run[[i]]@cpts)

  return(ans)
}




## Out of use:
## Simple plotter of pie-charts
plot.hit.pie <- function(hit.prop=0.5,miss.ncpts=0.4,title=""){
  hitdf = data.frame(Value=c(hit.prop,miss.ncpts,1-hit.prop-miss.ncpts),Group=c("Hit","Wrong ncpts","Missplaced"))

  p01<-ggplot(data=hitdf, aes(x="Proportion",y=Value,fill=Group)) + geom_col()+ coord_polar("y", start = 0)+theme_minimal() +scale_fill_manual(values=wes_palette( name="Cavalcanti"))+ggtitle(title)

  p01+ labs(subtitle=paste(c("v", version,"plot",str_pad(mycount, 3, pad = "0")), collapse = ""))
  mycount=do.plot.save(plot=p01,version=version)

}



## ----eval=FALSE----------------------------------------------------------
## simulate <- function(len=100,int.cpt=-1,cpt=-1,seed=-1){
##   # Version 1
##   if(seed !=-1){
##     set.seed(seed)
##   }else{
##     seed=NaN
##   }
##   if(cpt[1]==-1){
##     # if changepoint vector is not specified
##     if(int.cpt==-1){
##       # if number of changepoints not specified
##       # [1] Draw number of changepoints
##       int.cpt = abs(round(rnorm(1)*sqrt(len)+sqrt(len)))
##     }
##     # [2] Draw changepoint vector
##     cpt = sort(c(0,sample(1:(len-1), int.cpt, replace = FALSE),len))
##   }else{
##     # changepoint vector is specified
##     if((!is.element(0,cpt))||!(is.element(len,cpt))){
##       # if I constructed it wrong
##       cat('Element ', 0, ' is in cpt: ', is.element(0,cpt),',
##           Element ', len, 'is in cpt: ', is.element(len,cpt),'.
##           They both should be, and cpt should be sorted.')
##       return(NaN)
##     }
##     # set int. cpt accordingly
##     int.cpt=length(cpt)-2
##   }
##   # [3] Draw means
##   mean = runif(n=int.cpt+1)*20-10
##
##   # Remark that cpt[...] is vector, so times repeats each element times times
##   data=rnorm(len)+rep(mean,times=cpt[2:(int.cpt+2)]-cpt[1:(int.cpt+1)])
##   return(list(data=data,len=len, cpt=cpt,int.cpt=int.cpt, mean=mean,
##               info=list(seed=seed,version=1)))
## }
## simA = simulate(int.cpt = 10, seed=808)
## simA$beta_1= 2*log(simA$len)
##
## simA$PELT_BIC= cpt.mean(data= simA$data, penalty = "BIC",
##          method = "PELT", test.stat = "Normal", class = TRUE,
##          param.estimates = TRUE)
## simA$PELT_manual= cpt.mean(data= simA$data, penalty = "Manual", pen.value = simA$beta_1,
##          method = "PELT", test.stat = "Normal", class = TRUE,
##          param.estimates = TRUE)
##
## c('See that same answer and be happy.')
##



# analysis.hitprop.w.cpts2.02 <- function(data,sims,hit=c(1),pelt.yes=TRUE,Q=5,delta=0,true.m=1,cpts=10,n=20){
#   # When there is in truth 1 changepoint,
#   # which proportion finds the correct
#
#   #choose whether pelt or BinSe solution
#
#   # Q is maximum number of changepoints
#   # Only works for one changepoint
#   ans = list()
#
#   # applies PELT to an entire runset (that is one configuration of parameters)
#   if(pelt.yes){
#     run=cpt.mean(data=data,penalty="BIC",
#                  param.estimates=FALSE,method="PELT")
#   }else{
#     run=cpt.mean(data=data,penalty="BIC", Q=Q,
#                  param.estimates=FALSE,method="BinSeg")
#   }
#
#   # Save  number  of changepoints
#   ncpts= sapply(1:sims,function(i) ncpts(run[[i]]))
#
#
#   # Proportion with each number of changepoints
#   ans$ncpts=as.data.frame(table(ncpts))
#   ans$ncpts$Freq=ans$ncpts$Freq/sims
#
#   # Also include those that are different
#   different=setdiff(0:Q,ans$ncpts$ncpts)
#   df=data.frame(ncpts=different,Freq=rep(0,length(different)))
#   ans$m=rbind(df,ans$ncpts)
#   ans$m=ans$m[order(ans$m$ncpts),]
#   ans$m
#
#
#   # Proportion with correct number of changepoints
#   ans$ncpts.correct=ans$m$Freq[ans$m$ncpts==true.m]
#
#   # Add column identifying this sequence
#   if(delta!=0){
#     ans$ncpts$delta=rep(x=delta,times=length(ans$ncpts$ncpts))
#     ans$m$delta    =rep(x=delta,times=length(ans$m$ncpts))
#   }
#
#   # Save position of changepoints
#   #ans$cpts = lapply(1:sims,function(i) run[[i]]@cpts)
#
#   if(FALSE){
#     ## Only checks position of first changepoint
#     # Value of subjects with one changepoint, and a changepoint inside hit
#     ans$hit= sapply(1:sims,function(i) ifelse(ncpts(run[[i]])==true.m && is.element(el=run[[i]]@cpts[1],set = hit),run[[i]]@cpts[1],NA))
#     # Strip NA's
#     ans$hit <- ans$hit[!is.na(ans$hit)]
#   }else{
#     ans$hit= sapply(1:sims,function(i) ifelse(identical(as.numeric(run[[i]]@cpts),as.numeric(c(cpts,n))),run[[i]]@cpts[1],NA))
#     # Strip NA's
#     ans$hit <- ans$hit[!is.na(ans$hit)]
#   }
#
#   if(!length(ans$hit)){
#     ans$hit=data.frame(Value=1,Proportion=0)
#     ans$hit.prop=0
#   }else{
#     ### This is only done if there is some hits
#     # Proportion of simulations with hit
#     ans$hit.prop=length(ans$hit)/sims
#     if(delta!=0){
#       ans$correct = data.frame(correct=ans$hit.prop,
#                                delta=rep(x=delta,times=length(ans$ncpts$ncpts)))
#     }
#
#     # Proportion with each hit value
#     # want each hit value represented, so add them first and subtract later
#     # add:
#     ans$hit<- as.data.frame(table(c(ans$hit,hit)))
#     names(ans$hit)=c("Value","Proportion")
#     # subtract
#     ans$hit$Proportion=ans$hit$Proportion-rep(x=1,times=length(ans$hit$Proportion))
#     ans$hit$Proportion=ans$hit$Proportion/sims
#
#     # Add column identifying this sequence
#     if(delta!=0){
#       ans$hit$Delta=rep(x=delta,times=length(ans$hit$Proportion))
#     }
#   }
#
#   return(ans)
# }





# analysis.hitprop.w.cpts3 <- function(data,sims,hit=c(1),pelt.yes=TRUE,Q=5,delta=0,true.m=1,cpts=1){
#   # When there is in truth 1 changepoint,
#   # which proportion finds the correct
#
#   #choose whether pelt or BS solution
#
#   # Q is maximum number of changepoints
#   # Only works for one changepoint
#   ans = list()
#
#   # applies PELT to an entire runset (that is one configuration of parameters)
#   if(pelt.yes){
#     run=cpt.mean(data=data,penalty="BIC",
#                  param.estimates=FALSE,method="PELT")
#   }else{
#     run=cpt.mean(data=data,penalty="BIC", Q=Q,
#                  param.estimates=FALSE,method="BinSeg")
#   }
#
#   # Save  number  of changepoints
#   ncpts= sapply(1:sims,function(i) ncpts(run[[i]]))
#
#
#   # Proportion with each number of changepoints
#   ncpts2=as.data.frame(table(ncpts))
#   ncpts2$Freq=ncpts2$Freq/sims
#
#   # Also include those that are different
#   different=setdiff(0:Q,ncpts2$ncpts)
#   df=data.frame(ncpts=different,Freq=rep(0,length(different)))
#   ans$m=rbind(df,ncpts2)
#   ans$m=ans$m[order(ans$m$ncpts),]
#   ans$m
#
#
#   # Proportion with correct number of changepoints
#   ans$ncpts.correct=ans$m$Freq[ans$m$ncpts==true.m]
#
#   # Add column identifying this sequence
#   if(delta!=0){
#     ans$m$delta=rep(x=delta,times=length(ans$m$ncpts))
#     #ans$m$n    =rep(x=delta,times=length(ans$m$ncpts))
#   }
#
#   # Save position of changepoints
#   #ans$cpts = lapply(1:sims,function(i) run[[i]]@cpts)
#
#   ## Only checks position of first changepoint
#   # Value of subjects with one changepoint, and a changepoint inside hit
#   ans$hit= sapply(1:sims,function(i) ifelse(ncpts(run[[i]])==1 && is.element(el=run[[i]]@cpts[1],set = hit),run[[i]]@cpts[1],NA))
#   # Strip NA's
#   ans$hit <- ans$hit[!is.na(ans$hit)]
#
#   if(!length(ans$hit)){
#     ans$hit=data.frame(Value=1,Proportion=0)
#     ans$hit.prop=0
#   }else{
#     ### This is only done if there is some hits
#     # Proportion of simulations with hit
#     ans$hit.prop=length(ans$hit)/sims
#     if(delta!=0){
#       ans$correct = data.frame(correct=ans$hit.prop,
#                                delta=rep(x=delta,times=length(ans$ncpts$ncpts)))
#     }
#
#     # Proportion with each hit value
#     # want each hit value represented, so add them first and subtract later
#     # add:
#     ans$hit<- as.data.frame(table(c(ans$hit,hit)))
#     names(ans$hit)=c("Value","Proportion")
#     # subtract
#     ans$hit$Proportion=ans$hit$Proportion-rep(x=1,times=length(ans$hit$Proportion))
#     ans$hit$Proportion=ans$hit$Proportion/sims
#
#     # Add column identifying this sequence
#     if(delta!=0){
#       ans$hit$Delta=rep(x=delta,times=length(ans$hit$Proportion))
#       #ans$hit$n=dim(data)
#     }
#   }
#
#   return(ans)
# }





# analysis.hitprop.w.cpts2.01 <- function(data,sims,hit=c(1),pelt.yes=TRUE,Q=5,delta=0,true.m=1){
#   # When there is in truth 1 changepoint,
#   # which proportion finds the correct
#
#   #choose whether pelt or BS solution
#
#   # Q is maximum number of changepoints
#   # Only works for one changepoint
#   ans = list()
#
#   # applies PELT to an entire runset (that is one configuration of parameters)
#   if(pelt.yes){
#     run=cpt.mean(data=data,penalty="BIC",
#                  param.estimates=FALSE,method="PELT")
#   }else{
#     run=cpt.mean(data=data,penalty="BIC", Q=Q,
#                  param.estimates=FALSE,method="BinSeg")
#   }
#
#   # Save  number  of changepoints
#   ncpts= sapply(1:sims,function(i) ncpts(run[[i]]))
#
#
#   # Proportion with each number of changepoints
#   ans$ncpts=as.data.frame(table(ncpts))
#   ans$ncpts$Freq=ans$ncpts$Freq/sims
#
#   # Also include those that are different
#   different=setdiff(0:Q,ans$ncpts$ncpts)
#   df=data.frame(ncpts=different,Freq=rep(0,length(different)))
#   ans$m=rbind(df,ans$ncpts)
#   ans$m=ans$m[order(ans$m$ncpts),]
#   ans$m
#
#
#   # Proportion with correct number of changepoints
#   ans$ncpts.correct=ans$ncpts$Freq[ans$ncpts$ncpts==1]
#
#   # Add column identifying this sequence
#   if(delta!=0){
#     ans$ncpts$delta=rep(x=delta,times=length(ans$ncpts$ncpts))
#     ans$m$delta    =rep(x=delta,times=length(ans$m$ncpts))
#   }
#
#   # Save position of changepoints
#   #ans$cpts = lapply(1:sims,function(i) run[[i]]@cpts)
#
#   ## Only checks position of first changepoint
#   # Value of subjects with one changepoint, and a changepoint inside hit
#   ans$hit= sapply(1:sims,function(i) ifelse(ncpts(run[[i]])==1 && is.element(el=run[[i]]@cpts[1],set = hit),run[[i]]@cpts[1],NA))
#   # Strip NA's
#   ans$hit <- ans$hit[!is.na(ans$hit)]
#
#   if(!length(ans$hit)){
#     ans$hit=data.frame(Value=1,Proportion=0)
#     ans$hit.prop=0
#   }else{
#     ### This is only done if there is some hits
#     # Proportion of simulations with hit
#     ans$hit.prop=length(ans$hit)/sims
#     if(delta!=0){
#       ans$correct = data.frame(correct=ans$hit.prop,
#                                delta=rep(x=delta,times=length(ans$ncpts$ncpts)))
#     }
#
#     # Proportion with each hit value
#     # want each hit value represented, so add them first and subtract later
#     # add:
#     ans$hit<- as.data.frame(table(c(ans$hit,hit)))
#     names(ans$hit)=c("Value","Proportion")
#     # subtract
#     ans$hit$Proportion=ans$hit$Proportion-rep(x=1,times=length(ans$hit$Proportion))
#     ans$hit$Proportion=ans$hit$Proportion/sims
#
#     # Add column identifying this sequence
#     if(delta!=0){
#       ans$hit$Delta=rep(x=delta,times=length(ans$hit$Proportion))
#     }
#   }
#
#   return(ans)
# }


# analysis.hitprop.w.cpts <- function(data,sims,hit=c(1),pelt.yes=TRUE,Q=5,delta=0){
#   # When there is in truth 1 changepoint,
#   # which proportion finds the correct
#
#   #choose whether pelt or BS solution
#
#   # Q is maximum number of changepoints
#   # Only works for one changepoint
#   ans = list()
#
#   # applies PELT to an entire runset (that is one configuration of parameters)
#   if(pelt.yes){
#     run=cpt.mean(data=data,penalty="BIC",
#                  param.estimates=FALSE,method="PELT")
#   }else{
#     run=cpt.mean(data=data,penalty="BIC", Q=Q,
#                  param.estimates=FALSE,method="BinSeg")
#   }
#
#   # Save  number  of changepoints
#   ncpts= sapply(1:sims,function(i) ncpts(run[[i]]))
#
#   # Proportion miss with wrong number of cpts
#   ans$miss.ncpts=length(which(ncpts != 1))/sims
#
#   # Proportion with each number of changepoints
#   ans$ncpts=as.data.frame(table(ncpts))
#   ans$ncpts$Freq=ans$ncpts$Freq/sims
#
#   # Proportion with correct number of changepoints
#   ans$ncpts.correct=ans$ncpts$Freq[ans$ncpts$ncpts==1]
#
#   # Add column identifying this sequence
#   if(delta!=0){
#     ans$ncpts$delta=rep(x=delta,times=length(ans$ncpts$ncpts))
#
#   }
#
#   # Save position of changepoints
#   #ans$cpts = lapply(1:sims,function(i) run[[i]]@cpts)
#
#   ## Only checks position of first changepoint
#   # Value of subjects with one changepoint, and a changepoint inside hit
#   ans$hit= sapply(1:sims,function(i) ifelse(ncpts(run[[i]])==1 && is.element(el=run[[i]]@cpts[1],set = hit),run[[i]]@cpts[1],NA))
#   # Strip NA's
#   ans$hit <- ans$hit[!is.na(ans$hit)]
#
#   if(!length(ans$hit)){
#     ans$hit=data.frame(Value=1,Proportion=0)
#     ans$hit.prop=0
#   }else{
#     ### This is only done if there is some hits
#     # Proportion of simulations with hit
#     ans$hit.prop=length(ans$hit)/sims
#     if(delta!=0){
#       ans$correct = data.frame(correct=ans$hit.prop,
#                                delta=rep(x=delta,times=length(ans$ncpts$ncpts)))
#     }
#
#     # Proportion with each hit value
#     # want each hit value represented, so add them first and subtract later
#     # add:
#     ans$hit<- as.data.frame(table(c(ans$hit,hit)))
#     names(ans$hit)=c("Value","Proportion")
#     # subtract
#     ans$hit$Proportion=ans$hit$Proportion-rep(x=1,times=length(ans$hit$Proportion))
#     ans$hit$Proportion=ans$hit$Proportion/sims
#
#     # Add column identifying this sequence
#     if(delta!=0){
#       ans$hit$Delta=rep(x=delta,times=length(ans$hit$Proportion))
#     }
#   }
#
#   return(ans)
# }


