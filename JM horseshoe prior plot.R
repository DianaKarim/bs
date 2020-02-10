
library(MCMCpack)
library(logspline)
library(rmutil)
b <- 1
set.seed(123)
gammavec <- rgamma(1e5,shape=.7,rate=1/(2*b))
tau2vec <- rinvgamma(1e5,shape=.5,scale=gammavec/2)
betavec <- rnorm(1e5,mean=0,sd=sqrt(tau2vec))
# only plot the right side of 0 (density is symmetric around 0 anyway)
select1 <- which((betavec>0)*(betavec<10)==1)
plot(logspline(betavec[select1],lbound=0),ylim=c(0,3.5),xlim=c(0,10))
#compare with standard normal prior
thetaseq <- seq(0,10,length=1000)
lines(thetaseq,2*dnorm(thetaseq,mean=0,sd=1),col=2)
#compare with Cauchy prior
lines(thetaseq,2*dt(thetaseq,df=1),col=4)
lines(thetaseq,2*dlaplace(thetaseq),col=3)
legend(2,1.5,lty=rep(1,3),col=c(1,4,2,3),legend=c("horseshoe","Cauchy","normal","LaPlace"))



