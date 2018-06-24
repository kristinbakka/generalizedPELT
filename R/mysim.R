#-----------------------------------------------
set.window<- function(off=TRUE){
  if(off){dev.off();} # dev.new(width=6, height=3)
  windows(record=TRUE, width=5, height=2.5) #w=5,h=2.5 blir bra med width=0.8\textwidth i Latex (kanskje litt h?y figur)
}
set.window.2x<- function(off=TRUE){
  if(off){dev.off();} # dev.new(width=6, height=3)
  windows(record=TRUE, width=5, height=4.5) #w=5,h=2.5 blir bra med width=0.8\textwidth i Latex (kanskje litt h?y figur)
}
#-----------------------------------------------
## Functions needed for simulation

## Make data set
sigma.set.1 <- function(sig1){
  # Datasett med sig1 i sigma, p=1
  sigma=list(s1=c(1),s2=c(sig1))
  return(runs.mycpt(runs=5,minseglen=2,sigma=sigma,mu=list(m1=c(0),m2=c(0)),
                    order=c(1,2),each=10))
}

sigma.set.2 <- function(sig1){
  # Datasett med diag(sig1,sig1) i sigma, p=2
  sigma=list(s1=matrix(c(1,0,0,1),nrow=2),s2=matrix(c(sig1,0,0,sig1),nrow=2))
  return(runs.mycpt(runs=5,minseglen=2,pen=0,sigma=sigma,mu=list(mu1=c(0,0),mu2=c(0,0)),
                    order=c(1,2),each=10))
}

## Make data set from other data set
sep.1st.stream <- function(obj){
  # From p streams, separate out 1st stream
  # and adjust penalty
  obj$attrb$p=1
  obj$attrb$pen=3*log(obj$attrb$n)
  obj$data=lapply(obj$data,function(y) return(y=y[,1]))
  return(obj)
}

#---------- No changepoints, BIC penalty
# Generate data and set attributes

sim.0cpts.n <- function(vector.each,runs,minseglen,sigma,mu,order=c(1),multiplier=1,m=0){
  return(sim.cpts.n(vector.each,runs,minseglen,sigma,mu,order=c(1),multiplier=1,m=0))
}
sim.cpts.n <- function(vector.each,runs,minseglen,sigma,mu,order=c(1),multiplier=1,m=0){
  # Create one set of "runs" simulations for each element in "vector.each"
  # whith penalty log(n), so that it may easily be changed to the BIC-penalty
  # m comes in through order (m=(length(order)-1))
  temp1=lapply(vector.each,sim.cpts.n.call,runs=runs,minseglen=minseglen,
               sigma=sigma,mu=mu,order=order,multiplier=multiplier)
  return(temp1)
}

sim.cpts.n.call <- function(each,runs,minseglen,sigma,mu,order,multiplier){
  # Create a set of "runs" simulations for these parameters
  # with penalty prepared for BIC-penalty.
  pen=log(each*length(order))*multiplier
  return(runs.mycpt(runs=runs,minseglen=minseglen,pen=pen,sigma=sigma,mu=mu,
                    order=order,each=each))
  # Actually, it looks like this might be re-used as-is when not 0cpts.
}

regulate.attrb.pen.vec<-function(pen.mult,pen.const,type.names,sim){
  ## Multiply each penalty by pen.multiplier, add pen.const
  if(FALSE){
    #Debug
    sim=sim.b$dat1d

    #length skulle vært 7
    class(all.sims)
    length(all.sims)

    class(all.sims[[1]])
    length(all.sims[[1]])
    class(sim)
    length(sim)

    class(all.sims[[1]][[1]])
    length(all.sims[[1]][[1]])

    all.sims[[1]][[1]]$attrb
    all.sims[[1]][[1]]$attrb

  }
  if(length(pen.const)!=length(pen.const)){
    warning("\nThere is a problem. You see, \nlength(pen.const)!=length(pen.const), deary <3.\n")
    return(NA)
  }
  ## Multiply each penalty by pen.multiplier, add pen.const
  pen.len = length(pen.const)
  itrs=length(sim)
  # all.sims2<-rep(sim,pen.len) fungerer av uante grunner ikke (gir pen.len*itrs elementer, ikke pen.len)
  all.sims=list()
  for(j in 1:pen.len){
    all.sims[[j]]<-sim
  }
  for(k in 1:itrs){
    for(j in 1:pen.len){
    all.sims[[j]][[k]]$attrb$pen=sim[[k]]$attrb$pen*pen.mult[j]+pen.const[j]
    }
  }
  if(length(type.names)==pen.len){names(all.sims)=type.names}
  return(all.sims)
}


regulate.attrb.penalty<-function(pen.multiplier,pen.const=0,sim){
  ## Multiply each penalty by pen.multiplier
  itrs=length(sim)
  for(k in 1:itrs){
    sim[[k]]$attrb$pen=sim[[k]]$attrb$pen*pen.multiplier+pen.const
  }
  return(sim)
}

extract.atts <-function(sim){
  atts=lapply(sim,function(x) x$attrb)
  names(atts)=rep("attrb",length(atts))
  return(atts)
}
extract.dats <-function(sim){
  dats=lapply(sim,function(x) x$data)
  return(dats)
}

regulate.atts.penalty<-function(atts,pen.multiplier,pen.const=0){
  ## Multiply each penalty by pen.multiplier
  itrs=length(atts)
  for(k in 1:itrs){
    atts[[k]]$attrb$pen=atts[[k]]$attrb$pen*pen.multiplier+pen.const
  }
  return(atts)
}

#---------- No changepoints, BIC penalty--------
# Calculate statistics
sim.hitprop.0cpts <-function(sim.sol){
  return(sapply(sim.sol,hitprop.0cpts,runs=length(sim.sol[[1]])))
}
hitprop.0cpts <- function(runset.sol,runs){
  return(sum(sapply(runset.sol,function(x) ifelse(length(x)==1,1,0)))/runs)
}

sim.hitprop.ncpts <-function(sim.sol,sims){
  # Proportion of simulations with correct m
  # Assumes m in constant through observations
  m=length(sims[[1]]$attrb$changepoints)-1
  return(sapply(sim.sol,hitprop.ncpts,runs=length(sim.sol[[1]]),m=m))
}

hitprop.ncpts <- function(runset.sol,runs,m){
  return(sum(sapply(runset.sol,function(x) ifelse((length(x)-1)==m,1,0)))/runs)
}

#/// Exact
sim.hitprop.exact <-function(sim.sol,sims){
  #k=mapply(function(x,y) cat(x,y,"\n"),c("a","b","c"),c(1,2,3))
  temp=mapply(hitprop.exact, sim.sol,sims,runs=length(sim.sol[[1]]))
  return(temp)
}

hitprop.exact <- function(runset.sol,runset.sim,runs){
  # Proportion with exact same changepoint vector
  cpt.vec=runset.sim$attrb$changepoints
  return(sum(sapply(runset.sol,
    function(x) ifelse(identical(x,cpt.vec),1,0)))/runs)
}

## Only proportion of changepoints
sim.prop <-function(sim.sol,n.vec){
  #k=mapply(function(x,y) cat(x,y,"-"),c("a","b","c"),c(1,2,3))
  #m=length(sims[[1]]$attrb$changepoints)-1
  if(FALSE){
    sim.sol=pelt.a$meanvar$c0
    sims=sim.a$d1$meanvar$c0

    runset.sol=sim.sol[[1]]

    runset.sim$attrb$changepoints
    i=1
    i=1+1
    runset.sol[[i]]
    runset.sol[[10]]
    identical(runset.sol[[i]],runset.sim$attrb$changepoints)

    sims[[1]]$attrb$changepoints
    runset.sim$attrb$changepoints
  }
  library(dplyr)
  temp=bind_rows(mapply(prop,sim.sol,n.vec,runs=length(sim.sol[[1]]),SIMPLIFY = FALSE))
  colnames(temp)[1:2]=c("m","value")
  return(temp)
}

prop <- function(runset.sol,n,runs){
  # Proportion with exact same changepoint vector
  if(FALSE){
    # Debug
    x=runset.sol[[1]][[1]]
    runset.sol[[2]]
    runset.sol=sim.sol[[10]]
  }
  ms=as.data.frame(table(sapply(runset.sol,
               function(x) length(x)-1))/runs)
  ms[,1]=as.numeric(as.character(ms[,1]))
  ms[,2]=as.numeric((ms[,2]))
  ms$n=rep(n,length(ms[[1]]))
  return(ms)
}


#-------------


## Run PELT
attrb.mycpt <- function(p=1,n=100,pen=2*log(100),minseglen=1){
  return(list(p=p,n=n,pen=pen,minseglen=minseglen))
}

run.pelt.sims<-function(sims,type,count=TRUE){
  if(count)(cat('\nFinished all sims of this variety :D (tot=length(b.n))\n'))
  #templapply(sims,run.pelt.on.sim,type=type)
  return(lapply(sims,run.pelt.on.sim,type=type))
}

run.pelt.on.sim <-function(obj,type,count=TRUE){
  if(count)(cat('PELT has run 1 time (tot=runs)\n'))
  temp7<-pelt4.both.mycpt(obj$data,obj$attrb,type=type,both=FALSE)
  return(temp7)
}

run.pelt.on.attrb.data <-function(data,attrb,type,count=TRUE){
  #types=c("1d.mean","1d.meanvar","pd.mean","pd.meanvar.diag","pd.meanvar.full")
  if(count)(cat('PELT has run 1 time\n'))
  return(pelt4.both.mycpt(data,attrb,type=type,both=FALSE))
}

#------------
max.pen<-function(n,m){
  return(3*m*log(n)-(m+1)*log(m+1))
}

max.pen.rel<-function(n,m){
  # Ideally a number between 1 and 2
  return(2-(m+1)*log(m+1)/(m*log(n)))
}

min.pen.rel<-function(n,m){
  return(1+log(1-m/n)/(m*log(n)))
}


min.pen<-function(n,m){
  return(2*m*log(n)+log(1-m/n))
}
