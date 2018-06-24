## Here we see that gOP and gPELT give the same solutions,
# but that they are sometimes different from the solutions from the package changepoint

library(changepoint)
library(MASS)

types1=c("1d.mean","1d.meanvar","pd.meanvar.diag",
         "mbic.1d.mean","mbic.1d.meanvar","mbic.pd.meanvar.diag")
types2=c("pd.mean","pd.meanvar.diag","pd.meanvar.full",
         "mbic.pd.mean","mbic.pd.meanvar.diag","mbic.pd.meanvar.full")
##########
## Test the different algorithms on one short data set
## Dataset 1 - mbic
attrb.1=attrb.mycpt(p=1,n=7,pen=3*log(7),minseglen=1)
dat.1=c(-4.19 , -3.35 , -6.17 , 2.84 , -0.197 , 1.75 , 1.36)
permKi.cpt = cpt.mean(data=dat.1,penalty="MBIC",method="PELT",minseglen = 1, class=FALSE)
permOP = op.mycpt(dat=dat.1,
                  attrb=attrb.1,
                  type="mbic.1d.mean")
permGp = gpelt.mycpt(dat=dat.1,
                     attrb=attrb.1,
                     type="mbic.1d.mean")
permOP.cpt=build.solution.mycpt(permOP,n=7)
permKi.cpt
permOP.cpt

## Dataset 2 - bic
attrb.2=attrb.mycpt(p=1,n=8,pen=2*log(8),minseglen=1)
dat.2=c(2.15, -0.021, -1.35, -1.14,  2.22,  0.56,  1.03,  1.79)
permKi.cpt = cpt.mean(data=dat.2,penalty="BIC",method="PELT",minseglen = 1, class=FALSE)
permOP = op.mycpt(dat=dat.2,
                  attrb=attrb.2,
                  type="1d.mean")
permGp = gpelt.mycpt(dat=dat.2,
                     attrb=attrb.2,
                     type="1d.mean")
permOP.cpt=build.solution.mycpt(permOP,n=8)
permKi.cpt
permOP.cpt



##########
## Compare the solutions from of the algorithms on several data sets
## Only mean varies
data.set.1 <-runs.mycpt(runs=50,minseglen=1,pen=0,sigma=list(s1=1,s2=1),mu=list(mu1=0,mu2=1),
                      order=c(1,1,2,1,2,1,1,2),each=1)
data.set.1$data = lapply(data.set.1$data,function(x) as.vector(x))

####
## With BIC - only mean varies
attrb=data.set.1$attrb
attrb$pen=2*log(data.set.1$attrb$n)
# This is OP solution - gOP, and should always be correct
op.sol<- op.both.mycpt(data=data.set.1$data,attrb=attrb,type="1d.mean",both=FALSE)
# This is new gPELT
gpelt.sol <- gpelt.both.mycpt(data=data.set.1$data,attrb=attrb,type="1d.mean",both=FALSE)
# This is erroneus straight forward pelt
wrongpelt.sol <- pelt2.both.mycpt(data=data.set.1$data,attrb=attrb,type="1d.mean",both=FALSE)
# This solution is from package changepoint
firstpelt.sol<-lapply(data.set.1$data,cpt.mean,penalty="BIC",method="PELT",minseglen = 1, class=FALSE )

# This is the case handeles in the OP
are.cpt.vecs.identical(op.sol,gpelt.sol)
are.cpt.vecs.identical(wrongpelt.sol,firstpelt.sol)
are.cpt.vecs.identical(op.sol,firstpelt.sol)


####
## With mBIC - only mean varies
attrb.mbic=data.set.1$attrb
attrb.mbic$pen=3*log(data.set.1$attrb$n)
# This is OP solution - gOP, and should always be correct
op.sol.mbic <- op.both.mycpt(data=data.set.1$data,attrb=attrb.mbic,type="mbic.1d.mean",both=FALSE)
# This is new gPELT
gpelt.sol.mbic <- gpelt.both.mycpt(data=data.set.1$data,attrb=attrb.mbic,type="mbic.1d.mean",both=FALSE)
# This is erroneus straight forward pelt
wrongpelt.sol.mbic <- pelt2.both.mycpt(data=data.set.1$data,attrb=attrb.mbic,type="mbic.1d.mean",both=FALSE)
# This solution is from package changepoint
firstpelt.sol.mbic<-lapply(data.set.1$data,cpt.mean,penalty="MBIC",method="PELT",minseglen = 1, class=FALSE )

## Here they are all identical as it is the case handled in Killick (2012).
are.cpt.vecs.identical(op.sol.mbic,gpelt.sol.mbic)
are.cpt.vecs.identical(op.sol.mbic,wrongpelt.sol.mbic)
are.cpt.vecs.identical(op.sol.mbic,firstpelt.sol.mbic)

##########
## Compare the solutions from of the algorithms on several data sets
## Only mean and variance varies, minseglen is introduced

## With BIC - mean and variance
data.set.2 <-runs.mycpt(runs=50,minseglen=2,pen=0,sigma=list(s1=1,s2=1),mu=list(mu1=0,mu2=1),
                      order=c(1,1,2,2,1,2,2,1,1,2,2,1,1,1),each=1)
data.set.2$data = lapply(data.set.2$data,function(x) as.vector(x))

attrb.mv=data.set.2$attrb
attrb.mv$pen=3*log(data.set.2$attrb$n)
# This is OP solution - gOP, and should always be correct
op.sol.mv <- op.both.mycpt(data=data.set.2$data,attrb=attrb.mv,type="1d.meanvar",both=FALSE)
# This is new gPELT
gpelt.sol.mv <- gpelt.both.mycpt(data=data.set.2$data,attrb=attrb.mv,type="1d.meanvar",both=FALSE)
# This is erroneus straight forward pelt
wrongpelt.sol.mv <- pelt2.both.mycpt(data=data.set.2$data,attrb=attrb.mv,type="1d.meanvar",both=FALSE)
# This solution is from package changepoint
firstpelt.sol.mv<-lapply(data.set.2$data,cpt.meanvar,penalty="Manual",
                         pen.value=attrb.mv$pen,method="PELT",minseglen = 3, class=FALSE )


are.cpt.vecs.identical(op.sol.mv,gpelt.sol.mv)
are.cpt.vecs.identical(op.sol.mv,wrongpelt.sol.mv)
are.cpt.vecs.identical(op.sol.mv,firstpelt.sol.mv)





