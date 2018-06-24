### Test gPELT and gOP with multivariate Gaussian

types1=c("1d.mean","1d.meanvar","pd.meanvar.diag",
         "mbic.1d.mean","mbic.1d.meanvar","mbic.pd.meanvar.diag")
types2=c("pd.mean","pd.meanvar.diag","pd.meanvar.full",
         "mbic.pd.mean","mbic.pd.meanvar.diag","mbic.pd.meanvar.full")

q=list()

# Set Sigma
a1=matrix(rnorm(9),nrow=3)
a2=matrix(rnorm(9),nrow=3)
q$sigma=list(s1=a1%*%t(a1),s2=a2%*%t(a2))
# Set mean
q$mu= list(mu1=rnorm(3)*1,mu2=c(0,0,0))
q$each=5
q$order=c(1,2,1,2)

## One sigma
# Runs is number of simulations
obj1<-runs.mycpt(runs=10,minseglen=3,pen=-1,q$sigma,q$mu,q$order,q$each)
true.changepoints=obj1$attrb$changepoints

## In this detection method off-diagonal elements of Sigma are zero
attrb1=obj1$attrb
attrb1$pen=2*obj1$attrb$p*log(obj1$attrb$n) # The value of the penalty is an open research question. This is number of estimated parameters like in BIC.
changepoints.detected.diag=gpelt.both.mycpt(obj1$data,attrb1,
                                       type="pd.meanvar.diag",both=FALSE)

## In this detection method there are no restrictions on Sigma
attrb2=obj1$attrb
attrb2$pen=(obj1$attrb$p+obj1$attrb$p^2)*log(obj1$attrb$n) # The value of the penalty is an open research question. This is number of estimated parameters like in BIC.
changepoints.detected.full=gpelt.both.mycpt(obj1$data,attrb2,
                                            type="pd.meanvar.full",both=FALSE)



## Same with mBIC inspired cost
## In this detection method off-diagonal elements of Sigma are zero
attrb3=obj1$attrb
attrb3$pen=(1+2*obj1$attrb$p)*log(obj1$attrb$n) # The value of the penalty is an open research question. This is number of estimated parameters like in BIC.
changepoints.detected.diag.mbic=gpelt.both.mycpt(obj1$data,attrb3,
                                            type="mbic.pd.meanvar.diag",both=FALSE)

## In this detection method there are no restrictions on Sigma
attrb4=obj1$attrb
attrb4$pen=(1+obj1$attrb$p+obj1$attrb$p^2)*log(obj1$attrb$n) # The value of the penalty is an open research question. This is number of estimated parameters like in BIC.
changepoints.detected.full.mbic=gpelt.both.mycpt(obj1$data,attrb4,
                                            type="mbic.pd.meanvar.full",both=FALSE)

## The detected and the true changepoints
true.changepoints
changepoints.detected.diag
changepoints.detected.full

changepoints.detected.diag.mbic
changepoints.detected.full.mbic

