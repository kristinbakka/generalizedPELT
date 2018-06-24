
library(MASS)

###############################
## Objects and simulate #######

runs.mycpt <- function(runs,minseglen=75,pen=-1,sigma,mu,order,each){
  # It simulates one per run :)
  library(Matrix)
  library(stats)
  library(MASS)

  obj <- simulate.dataset.mycpt(sigma,mu,order,each,runs,
                  attrb=list(pen=pen,minseglen=minseglen))

  obj$attrb$n=length(order)*each
  obj$attrb$m=length(order)-1

  if(obj$attrb$pen==-1){obj$attrb$pen=log(obj$attrb$n)} # prepare for BIC

  return(obj)
}

simulate.1.mycpt <- function(sigma,mu,order,each){
  element=c()
  for(x in order){
    element=rbind(element,mvrnorm(n=each,mu=mu[[x]],Sigma=sigma[[x]]))
  }
  return(element)
}

simulate.dataset.mycpt <- function(sigma,mu,order,each,runs=2,return.data.only=FALSE,attrb=list()){
  if(return.data.only){
    return(replicate(n=runs,expr=simulate.1.mycpt(sigma,mu,order,each),simplify = FALSE))
  }

  attrb2 <- list(runs=runs,n=each*length(order),p=length(mu[[1]]),changepoints=(1:length(order))*each,sigma=sigma,mu=mu,order=order,each=each)

  return(list(attrb=c(attrb,attrb2),
              data=replicate(n=runs,simulate.1.mycpt(sigma,mu,order,each),simplify = FALSE)))
}

###############################
## Cost computation ###########

cost.tot.sol.mycpt <- function(data,tau.vec,type="1d.mean",pen=0){
  # Handleable form so same whether p=1 or not
  data=matrix(data)
  # Last element, also length(data) [or dim(data)[1]]
  len=tail(tau.vec,1)
  # Remove last element of tau.vec
  tau.vec=tau.vec[1:(length(tau.vec)-1)]
  # Number of internal changepoints
  m=length(tau.vec)
      # Compute interval cost for each of m+1 intervals
  int.cost = sapply(1:(m+1), function(i) cost.mycpt(intv.dat=
      # Where the i th interval is  (c(0,tau.vec)[i]+1):(c(tau.vec,len)[i])
        data[(c(0,tau.vec)[i]+1):(c(tau.vec,len)[i]),],type=type))
  return(sum(int.cost)+m*pen)
}

cost.mycpt <- function(intv.dat,type="1d.mean",n=1){
  return(
  switch(type,
         "1d.mean"=cost.1d.mean.mycpt(intv.dat=intv.dat),
         "1d.meanvar"=cost.1d.meanvar.mycpt(intv.dat=intv.dat),
         "pd.mean"=cost.pd.mean.mycpt(intv.dat=intv.dat),
         "pd.meanvar.diag"=cost.pd.meanvar.diag.mycpt(intv.dat=intv.dat),
         "pd.meanvar.full"=cost.pd.meanvar.full.mycpt(intv.dat=intv.dat),

         "mbic.1d.mean"=cost.mbic.1d.mean.mycpt(intv.dat=intv.dat,n=n),
         "mbic.1d.meanvar"=cost.mbic.1d.meanvar.mycpt(intv.dat=intv.dat,n=n),
         "mbic.pd.mean"=cost.mbic.pd.mean.mycpt(intv.dat=intv.dat,n=n),
         "mbic.pd.meanvar.diag"=cost.mbic.pd.meanvar.diag.mycpt(intv.dat=intv.dat,n=n),
         "mbic.pd.meanvar.full"=cost.mbic.pd.meanvar.full.mycpt(intv.dat=intv.dat,n=n)
  )
  )
}

cost.1d.mean.mycpt <- function(intv.dat,t=0){
  return(sum((intv.dat-mean(intv.dat))^2))
}

cost.1d.meanvar.mycpt <- function(intv.dat,t=0){
  t.n=length(intv.dat)
  #sigma.sq.hat=(t.n-1)*var(intv.dat)/t.n
  sigma.sq.hat=sum((intv.dat-mean(intv.dat))^2)/t.n
  if(sigma.sq.hat<0.0000000001){
    sigma.sq.hat=0.0000000001
  }
  return(t.n*log(sigma.sq.hat))
}

cost.pd.mean.mycpt <- function(intv.dat,t=0){
  # When Sigma is known to be I_p
  mu.hat=colMeans(intv.dat)
  return(sum((intv.dat-mu.hat)^2))
}

cost.pd.meanvar.diag.mycpt <- function(intv.dat,t=0){
  log.sigma.sq.hat = log(colSums((intv.dat-colMeans(intv.dat))^2)/dim(intv.dat)[1])
  return(dim(intv.dat)[1]*sum(log.sigma.sq.hat))
}

cost.pd.meanvar.full.mycpt <- function(intv.dat,t=0){
  ## intv.dat had one time stamp in the same row. Each column is a stream
  ## Fits a p-dim normal to data and returns cost and mean,var
  # Number of observations
  len = dim(intv.dat)[1]
  p = dim(intv.dat)[2]

  # Mean ML-estimate
  mu.hat=colMeans(intv.dat)
  # Subtract mean from data
  z=as.matrix(sweep(intv.dat,2,mu.hat))

  ## Compute sigma.hat
  # For every row compute t(x_i-mu)(x_i-mu) and sum over i
  sigma.hat=apply(z, 1, function(x) t(z)%*%(z))
  # Sum each t(x-mu)(x-mu), put into matrix, divide by normalizing
  sigma.hat=matrix(sigma.hat[,1],ncol=p,nrow=p )/len

  ## SVD
  # Is it possible to not get an eigenvalue?
  # How large does an eigenvalue need to be before it is counted in?
  # Is the req fulfilled?
  # Compute eigenvalues
  eigen = svd(x=sigma.hat,nu=0,nv=0)
  # kutte ut numerisk null
  eigen$d = eigen$d[eigen$d>10^-10]
  # Compute cost based on rank=eigen$d and det(S)=prod(eigen$d)
  cost=len*(length(eigen$d)*(log(2*pi)+1)+log(prod(eigen$d))) #cost.K is negative
  return(cost)
}



cost.mbic.1d.mean.mycpt <- function(intv.dat,t=0,n){
  return(sum((intv.dat-mean(intv.dat))^2)+log(length(intv.dat)/n))
}

cost.mbic.1d.meanvar.mycpt <- function(intv.dat,t=0,n){
  t.n=length(intv.dat)
  #sigma.sq.hat=(t.n-1)*var(intv.dat)/t.n
  sigma.sq.hat=sum((intv.dat-mean(intv.dat))^2)/t.n
  if(sigma.sq.hat<0.0000000001){
    sigma.sq.hat=0.0000000001
  }
  return(t.n*log(sigma.sq.hat)+log(length(intv.dat)/n))
}

cost.mbic.pd.mean.mycpt <- function(intv.dat,t=0,n){
  # When Sigma is known to be I_p
  mu.hat=colMeans(intv.dat)
  return(sum((intv.dat-mu.hat)^2)+log(length(intv.dat)/n))
}

cost.mbic.pd.meanvar.diag.mycpt <- function(intv.dat,t=0,n){
  log.sigma.sq.hat = log(colSums((intv.dat-colMeans(intv.dat))^2)/dim(intv.dat)[1])
  return(dim(intv.dat)[1]*sum(log.sigma.sq.hat)+log(length(intv.dat)/n))
}

cost.mbic.pd.meanvar.full.mycpt <- function(intv.dat,t=0,n){
  ## intv.dat had one time stamp in the same row. Each column is a stream
  ## Fits a p-dim normal to data and returns cost and mean,var
  # Number of observations
  len = dim(intv.dat)[1]
  p = dim(intv.dat)[2]

  # Mean ML-estimate
  mu.hat=colMeans(intv.dat)
  # Subtract mean from data
  z=as.matrix(sweep(intv.dat,2,mu.hat))

  ## Compute sigma.hat
  # For every row compute t(x_i-mu)(x_i-mu) and sum over i
  sigma.hat=apply(z, 1, function(x) t(z)%*%(z))
  # Sum each t(x-mu)(x-mu), put into matrix, divide by normalizing
  sigma.hat=matrix(sigma.hat[,1],ncol=p,nrow=p )/len

  ## SVD
  # Is it possible to not get an eigenvalue?
  # How large does an eigenvalue need to be before it is counted in?
  # Is the req fulfilled?
  # Compute eigenvalues
  eigen = svd(x=sigma.hat,nu=0,nv=0)
  # kutte ut numerisk null
  eigen$d = eigen$d[eigen$d>10^-10]
  # Compute cost based on rank=eigen$d and det(S)=prod(eigen$d)
  cost=len*(length(eigen$d)*(log(2*pi)+1)+log(prod(eigen$d))) #cost.K is negative
  return(cost+log(length(intv.dat)/n))
}


###############################
## Other important functions ##

build.solution.mycpt <-function(permanent,n){
  ## Build solution
  i=1
  tau=rep(NA,n)
  cpt=n
  while(cpt!=0){
    # assign that this is a changepoint
    tau[i]=cpt
    # Previous changepoint is r(changepoint) (=r[cpt+1])
    cpt = permanent$r[cpt+1]
    i=i+1
  }
  # Reverse vector and strip NAs
  #
  return(rev(tau[!is.na(tau)]))
}


are.cpt.vecs.identical <- function(sol1,sol2){
  return(mapply(function(x,y) identical( as.integer(x),as.integer(y)),sol1,sol2))
}

###############################
## PELT and competing #########

## PELT that allows restriction min(tau_j - tau_(j-1))>1
gpelt.both.mycpt <-function(data,attrb,type="1d.mean",both=TRUE){
  # Calculate PELT
  #pelt.mat=lapply(data,function(x) pelt4.mycpt(attrb=attrb,dat=x,type=type,mBIC.style=FALSE))
  pelt.mat=lapply(data,function(x) gpelt.mycpt(attrb=attrb,dat=x,type=type))


  #Should be same
  #pelt.mat=lapply(data,pelt4.mycpt,attrb=attrb,type=type)


  pelt.cpts = lapply(pelt.mat,function(x) build.solution.mycpt(permanent = x, n=attrb$n))
  if(both){
    return(pelt=list(cpts=pelt.cpts,permanent=pelt.mat))
  }
  return(pelt.cpts)
}

gpelt.mycpt <- function(attrb,dat,type="1d.mean"){
  # Not manually debugging
  my.debug=FALSE

  # This is an overly complex way to do it, but gives a table
  # "permantent"
  # that is easier to interpret to understand the algorithm


  # Is type among the selection of cost functions
  if(attrb$p==1){
    if(!is.element(type,c("1d.mean","1d.meanvar","pd.meanvar.diag",
                          "mbic.1d.mean","mbic.1d.meanvar","mbic.pd.meanvar.diag"))){
      return("Type is not valid.")
    }
    dat=matrix(dat,ncol=1)
  }else{
    if(!is.element(type,c("pd.mean","pd.meanvar.diag","pd.meanvar.full",
                          "mbic.pd.mean","mbic.pd.meanvar.diag","mbic.pd.meanvar.full"))){
      return("Type is not valid.")
    }
  }


  ### Initialize first step such that
  # inherit = 0, F(0) = -\pen, s.set={0}, r(0)=0
  # Outer data frame of t,F,r
  permanent <- data.frame(t=seq(0,attrb$n),F.val=rep(NA,attrb$n+1),r=rep(NA,attrb$n+1))
  permanent[1,2:3]=c(-attrb$pen,0)


  ### Initialize first step such that
  for(t in attrb$minseglen:min(2*attrb$minseglen-1,attrb$n)){
    # predecessor is 0th data point
    permanent[permanent$t==t,2:3]=
      c(cost.mycpt(intv.dat=dat[(1):t,],type=type,n=attrb$n),0)
  }
  # Return if finished
  if(attrb$n<2*attrb$minseglen){
    return(permanent)
  }
  # Else construct Inherit such that
  # When we inherit fromm time t, we get the s.set at Inherit[[t+1]]
  #inherit$q is the data point we inherit from,
  # inherit$s is the pruned s.set at the time we inherit from
  Inherit=as.list(c(rep(0,2*attrb$minseglen),rep(NA,attrb$n-3*attrb$minseglen+1)))

  if(my.debug){t=2*attrb$minseglen-1}


  ####
  ## Compute for the rest of the data points
  for(t in (2*attrb$minseglen):(attrb$n)){
    if(my.debug&&(t%%25==0)){cat("t=",t,".\n")}
    if(my.debug){t=t+1}
    ###  Combine inherited and earned data points to get s.set
    s.set=c(Inherit[[t-attrb$minseglen+1]], #inherited
            max(attrb$minseglen,t-2*attrb$minseglen+1):(t-attrb$minseglen)) #earned


    ###  For a changepoint at t find best most recent changepoint s
    # Use cost function to compute int.cost C(s+1,t) for all s in s.set
    temp<-data.frame(
      s=s.set,
      int.cost = sapply(s.set, function(x) cost.mycpt(intv.dat=dat[(x+1):t,],type=type,n=attrb$n))
    )
    ## Compute full cost and pruning cost
    temp$full.cost <- permanent[s.set+1,2] + temp$int.cost + attrb$pen
    temp$prune.cost<- permanent[s.set+1,2] + temp$int.cost

    ## Determine smallest (optimal) full cost
    # Save smallest (optimal) full cost
    permanent$F.val[t+1]=min(temp$full.cost)

    # Save previous changepoint, the s with smallest full cost
    # That is the last s for which F.val is minimal
    permanent$r[t+1]=tail(temp$s[temp$full.cost==permanent$F.val[t+1]],1)


    ### Remove non-optimal predecessors
    ### Remember which data points to inherit
    # s with smaller pre-beta cost, the ones to keep
    A=temp$prune.cost<=permanent$F.val[t+1]

    ## Only add element to next s.set if it has a valid predecessor
    if(length(A==TRUE)==0){
      Inherit[[t+1]]=NULL
    }else{
      Inherit[[t+1]]=temp$s[A]
    }

    # Debug
    if(my.debug){t}
    if(my.debug){s.set} #current s.set to go through, out to be 0 until 2*minseglen

    if(my.debug){t}
    if(my.debug){temp}
    if(my.debug){Inherit[[t+1]]}
  #  if(my.debug){inherit$s[inherit$q==t-(attrb$minseglen)]} #Inherited, first part of s.set
    if(my.debug){t-(attrb$minseglen)} # Inherited from
    if(my.debug){inherit$s[inherit$q==t]} # legacy (inheritance passed on from this node (ouht to be 0 until 2*attrb$minseglen)


    if(my.debug){temp}
    if(my.debug){permanent}
    if(my.debug){View(permanent)}
  }
  if(FALSE){cat('\n1 run of gPELT performed.\n')}
  return(permanent)
}


## The OP method that allows restriction min(tau_j - tau_(j-1))>1
op.both.mycpt <-function(data,attrb,type="1d.mean",both=TRUE){
  # Calculate PELT
  pelt.mat=lapply(data,function(x) op.mycpt(attrb=attrb,dat=x,type=type))
  pelt.cpts = lapply(pelt.mat,function(x) build.solution.mycpt(permanent = x, n=attrb$n))
  if(both){
    return(pelt=list(cpts=pelt.cpts,permanent=pelt.mat))
  }
  return(pelt.cpts)
}

op.mycpt <- function(attrb,dat,type="1d.mean"){
  # exact same as OP, nothing is ever pruned
  # Not manually debugging
  my.debug=FALSE
  # Is type among the selection of cost functions
  if(attrb$p==1){
    if(!is.element(type,c("1d.mean","1d.meanvar","pd.meanvar.diag",
                          "mbic.1d.mean","mbic.1d.meanvar","mbic.pd.meanvar.diag"))){
      return("Type is not valid.")
    }
    dat=matrix(dat,ncol=1)
  }else{
    if(!is.element(type,c("pd.mean","pd.meanvar.diag","pd.meanvar.full",
                          "mbic.pd.mean","mbic.pd.meanvar.diag","mbic.pd.meanvar.full"))){
      return("Type is not valid.")
    }
  }

  ## Initialize first step such that
  # s = 0, F(0) = -\pen, s.set={0}, r(0)=0
  # Outer data frame of t,F,r
  permanent <- data.frame(t=seq(0,attrb$n),F.val=rep(NA,attrb$n+1),r=rep(NA,attrb$n+1))
  permanent[1,2:3]=c(-attrb$pen,0)
  s.set=c(0)
  if(my.debug){t=attrb$minseglen-1}

  ## Compute for all data sets lengths shorter than attrb$n+1
  # Work in delay by starting at minseglen
  for(t in (attrb$minseglen):attrb$n){
    if(my.debug&&(t%%25==0)){cat("t=",t,".\n")}
    if(my.debug){t=t+1}

    ## Use cost function to compute int.cost C(s+1,t) for all s in s.set
    # This is the only place the cost function is evaluated
    temp<-data.frame(
      s=s.set,
      int.cost = sapply(s.set, function(x) cost.mycpt(intv.dat=dat[(x+1):t,],type=type,n=attrb$n))
    )

    # This is an overly complex way to do it, but gives a table
    # "permantent" that is easier to interpret to unerstand the algorithm
    ## Compute full cost and pruning cost
    temp$full.cost <- permanent[s.set+1,2] + temp$int.cost + attrb$pen
    temp$prune.cost<- permanent[s.set+1,2] + temp$int.cost

    ## Determine smallest (optimal) full cost
    # Save smallest (optimal) full cost
    #¤Edit [t+1] to [permanent$s==t], maybe if not slower ;)
    permanent$F.val[t+1]=min(temp$full.cost)

    # Save previous changepoint, the s with smallest full cost
    # That is the last s for which F.val is minimal
    permanent$r[t+1]=tail(temp$s[temp$full.cost==permanent$F.val[t+1]],1)

    ## Prune - prepare next s.set
    # s with smaller pre-beta cost
    A=temp$prune.cost<=permanent$F.val[t+1]
    # or superceding t
    B=temp$s>permanent$r[t+1] #####
    #    B=rep(FALSE,length(A))
    #    if(B&!A){
    #      warning(paste("B&!A for t=",t,".\n"))
    #    }

    ## Only add element to set if it has a valid predecessor
    if(t>=(2*attrb$minseglen-1)){s.set = c(s.set,t+1-(attrb$minseglen))}

    # Debug
    if(my.debug){cat("t=",t,".\n")}
    if(my.debug){temp}
    if(my.debug){permanent}
    if(my.debug){View(permanent)}
  }
  return(permanent)
}


## Working title of op
pelt5.both.mycpt<- function(data,attrb,type="1d.mean",both=TRUE){
  return(op.both.mycpt(data,attrb,type,both))
}
pelt5.mycpt<- function(attrb,dat,type="1d.mean"){
  return(op.mycpt(attrb,dat,type))
}

## Gives same result as PELT in changepoint package, but with my implementation of mBIC
pelt2.both.mycpt <-function(data,attrb,type="1d.mean",both=TRUE){
  # Calculate PELT
  pelt.mat=lapply(data,function(x) pelt2.mycpt(attrb=attrb,dat=x,type=type))
  pelt.cpts = lapply(pelt.mat,function(x) build.solution.mycpt(permanent = x, n=attrb$n))
  if(both){
    return(pelt=list(cpts=pelt.cpts,permanent=pelt.mat))
  }
  return(pelt.cpts)
}

pelt2.mycpt <- function(attrb,dat,type="1d.mean"){
  # Not manually debugging
  my.debug=FALSE
if(attrb$p==1){
  if(!is.element(type,c("1d.mean","1d.meanvar","pd.meanvar.diag",
                        "mbic.1d.mean","mbic.1d.meanvar","mbic.pd.meanvar.diag"))){
    return("Type is not valid.")
  }
  dat=matrix(dat,ncol=1)
}else{
  if(!is.element(type,c("pd.mean","pd.meanvar.diag","pd.meanvar.full",
                        "mbic.pd.mean","mbic.pd.meanvar.diag","mbic.pd.meanvar.full"))){
    return("Type is not valid.")
  }
}
  ## Initialize first step such that
  # s = 0, F(0) = -\pen, s.set={0}, r(0)=0
  # Outer data frame of t,F,r
  permanent <- data.frame(t=seq(0,attrb$n),F.val=rep(NA,attrb$n+1),r=rep(NA,attrb$n+1))
  permanent[1,2:3]=c(-attrb$pen,0)
  s.set=c(0)
  if(my.debug){t=attrb$minseglen-1}

  ## Compute for all data sets lengths shorter than attrb$n+1
  # Work in delay by starting at minseglen
  for(t in (attrb$minseglen):attrb$n){
    if(my.debug&&(t%%25==0)){cat("t=",t,".\n")}
    if(my.debug){t=t+1}

    ## Use cost function to compute int.cost C(s+1,t) for all s in s.set
    # This is the only place the cost function is evaluated
    temp<-data.frame(
      s=s.set,
      int.cost = sapply(s.set, function(x) cost.mycpt(intv.dat=dat[(x+1):t,],type=type,n=attrb$n))
    )

    # This is an overly complex way to do it, but gives a table
    # "permantent" that is easier to interpret to unerstand the algorithm
    ## Compute full cost and pruning cost
    temp$full.cost <- permanent[s.set+1,2] + temp$int.cost + attrb$pen
    temp$prune.cost<- permanent[s.set+1,2] + temp$int.cost

    ## Determine smallest (optimal) full cost
    # Save smallest (optimal) full cost
    #¤Edit [t+1] to [permanent$s==t], maybe if not slower ;)
    permanent$F.val[t+1]=min(temp$full.cost)

    # Save previous changepoint, the s with smallest full cost
    # That is the last s for which F.val is minimal
    permanent$r[t+1]=tail(temp$s[temp$full.cost==permanent$F.val[t+1]],1)

    ## Prune - prepare next s.set
    # s with smaller pre-beta cost
    A=temp$prune.cost<=permanent$F.val[t+1]
    # or superceding t
    #skip this

    ## Only add element to next s.set if it has a valid predecessor
    if(t>=(2*attrb$minseglen-1)){s.set = c(temp$s[A],t-(attrb$minseglen-1))}

    # Debug
    if(my.debug){temp}
    if(my.debug){permanent}
    if(my.debug){View(permanent)}
  }
  return(permanent)
}


###############################
## Other functions ############

est.param.mycpt<-function(obj,tau.vec,type){
  est=list(type="1d.mean",tau.vec=tau.vec)
  tvt=c(0,tau.vec)

  # Set names
  if(is.element(type,c("1d.mean","pd.meanvar.diag"))){
    # Estimate the mean
    param <- lapply(1:length(tau.vec), function(i) {
      list(cpt =  tau.vec[i], mean= colMeans(data[(tvt[i]+1):tvt[i+1],]))})
  }else if(is.element(type,c("1d.meanvar","pd.meanvar.diag"))){
    # Estimate the mean and the variance of each
    param <- lapply(1:length(tau.vec), function(i) {
      list(cpt =  tau.vec[i], mean= colMeans(data[(tvt[i]+1):tvt[i+1],]),
           sigma= colSums((data[(tvt[i]+1):tvt[i+1],]-
                             colMeans(data[(tvt[i]+1):tvt[i+1],]))^2)/(tvt[i+1]-tvt[i])
      )})
  }else if(type=="pd.meanvar.full"){
    # Estimate the mean and the variance of each
    param <- lapply(1:length(tau.vec), function(i) {
      list(cpt =  tau.vec[i], mean= colMeans(data[(tvt[i]+1):tvt[i+1],]),
           sigma= est.param.full.mycpt(data[(tvt[i]+1):tvt[i+1],])
      )})

  }else{
    warning(paste("Not valid type =",type,".\nChoose among types:\n","1d.mean, ","1d.meanvar, ",
                  "pd.mean, \n","pd.meanvar.diag, ","pd.meanvar.full."))
  }

  setNames(param,sapply(1:length(tau.vec),function(i) paste0("c",i)))
  est$param=param
  return(est)
}

est.list.mean.mycpt <- function(tau){
  return(list(cpt=cpt,mean=colMeans(intv.dat)))
}


simulate.mycpt <-function(n=10,p=5,fixed=FALSE){
  warning("Can I delete this, or is it in use? [simulate.mycpt]")
  ## No changes in the data set
  library('mvtnorm')
  # Simulates a data set that I may use later
  # Every column
  if(fixed){set.seed(0)}
  mu=round(rnorm(p,mean=10,sd=8))
  q=matrix(rnorm(p*p,mean=2,sd=5),nrow=p,ncol=p)
  sigma=q%*%t(q)
  return(list(dat=rmvnorm(n=n, mean=mu,sigma=sigma),mean=mu,sigma=sigma))
}

#types=c("1d.mean","1d.meanvar","pd.mean","pd.meanvar.diag","pd.meanvar.full")
mycpt <- function(dat=0,sim=FALSE,n=100,p=1,minseglen=75,pen=-1){
  # It simulates one single
  library(Matrix)
  library(stats)

  warning("Check if this is still in use, I don't think so. :)")

  ## Simulate dataset
  if(sim){
    sim1 <- simulate.mycpt(n=n,p=p)
    dat  <- sim1$dat
    sigma<- sim1$sigma
    mean <- sim1$mean
  }



  ## Save parameters

  attrb <- list(minseglen=minseglen)
  if(sim){attrb$genparams <- list(mean=mean,sigma=sigma)}
  attrb$n=ifelse(is.null(dim(dat)),length(dat),dim(dat)[1])
  attrb$p=ifelse(is.null(dim(dat)),1,dim(dat)[2])

  # Set penalty term \beta
  attrb$pen <- ifelse(pen==-1,0,pen)

  obj=list(attrb=attrb,dat=as.matrix(dat))
  class(obj) <- 'mycpt'
  return(obj)
}



