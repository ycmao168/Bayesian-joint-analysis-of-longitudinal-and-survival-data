# intial for ita is a vector of (splineknotslongitudinal+splinedegree-1) components, for r is a vector of (splineknotssurvival+splinedegree-1) components, and the number of components in the theta and beta vectors is the number of columns of the covariate matrix
jointanalysis <- function(longitudinal,survival,splineknotslongitudinal,splineknotssurvival,splinedegree,varu,varxi0,varxi1,covxi,vare,c,ita,theta,r0,r,beta, burnin,n){
y <- longitudinal[[1]]
time <- longitudinal[[2]]
subject <- longitudinal[[3]]
delta1 <- survival[[1]]
delta2 <- survival[[2]]
delta3 <- survival[[3]]
L <- survival[[4]]
R <- survival[[5]]
x <- matrix(cbind(survival[[6]],survival[[7]]),ncol=2)
nobs <- survival[[9]]
# sample size
B <- length(nobs)
# number of covariates
p <- ncol(x)
# number of spline functions
jd1 <- splineknotslongitudinal+splinedegree-1
jd2 <- splineknotssurvival+splinedegree-1
# generate spline basis functions
union <- c(L[which(!is.na(L))],R[which(!is.na(R))])
s <- seq(1,B,1)
union <- c(L,R)
unionid <- c(s,s)
unionnewid <- unionid[which(!is.na(union))]
newunion <- union[which(!is.na(union))]
knots1 <- quantile(time,prob=seq(0,1,1/(splineknotslongitudinal-1)))[-1][-(splineknotslongitudinal-1)]
knots2 <- quantile(newunion,prob=seq(0,1,1/(splineknotssurvival-1)))[-1][-(splineknotssurvival-1)]
library(splines2)
msMat <- mSpline(time,knots=knots1,degree=splinedegree,intercept=TRUE)
isMat <- iSpline(newunion,knots=knots2,degree=splinedegree,intercept=TRUE)
isMatl <- matrix(NA,B,jd2)
isMatr <- matrix(NA,B,jd2)
for(i in 1:length(unionnewid)){
if(unionnewid[i+1] < unionnewid[i])
{k=i+1;
break}
}
for(pp in 1:k-1){
isMatl[unionnewid[pp],] <- isMat[pp,]
}
for(q in k:length(unionnewid)){
isMatr[unionnewid[q],] <- isMat[q,]
}
talpha <- matrix(0,B,jd2)
for(l in 1:B){
if(delta1[l]==1){talpha[l,]<-isMatr[l,]}else{talpha[l,]<-isMatl[l,]}}
# function to generate random effects
generateu <- function(xiexp,zexp,varu,vare,c,theta,beta,ita,r0,r){
mu <- msMat%*%ita
xtheta <- x%*%theta
sum1 <- y-mu-rep(xtheta,nobs)-rep(xiexp[,1],nobs)-rep(xiexp[,2],nobs)*time
dat1 <- data.frame(cbind(subject,sum1))
totl <- aggregate(sum1~subject,dat1,sum)[,2]
tots <- zexp-talpha%*%r-r0-x%*%beta
bb <- rep(0,B)
aa <- rep(0,B)
for(l in 1:B){
bb[l] <- totl[l]/vare+tots[l]*c
aa[l] <- (nobs[l]*varu+c^2*vare*varu+vare)/(2*vare*varu)
}
expu <- rep(0,B)
for(i in 1:B){
expu[i] <- rnorm(1,bb[i]/(2*aa[i]),sqrt(1/(2*aa[i])))
}
return(expu)
}
# function to generate random effects over time
generatexi <- function(uexp,varxi0,varxi1,covxi,vare,theta,ita){
dat1 <- data.frame(cbind(subject,time,time2=time^2))
tot1 <- aggregate(time~subject,dat1,sum)[,2]
tot2 <- aggregate(time2~subject,dat1,sum)[,2]
ximatrix <- array(0,dim=c(2,2,B))
for(l in 1:B){
ximatrix[,,l] <- matrix(c(nobs[l],tot1[l],tot1[l],tot2[l]),ncol=2)
}
sumxi <- y-msMat%*%ita-rep(uexp,nobs)-rep(x%*%theta,nobs)
dat2 <- data.frame(cbind(subject,sumxi,sumxi*time))
tot3 <- matrix(cbind(aggregate(dat2[,2]~subject,dat2,sum)[,2], aggregate(dat2[,3]~subject,dat2,sum)[,2]),ncol=2)
sigmaxi0 <- matrix(c(varxi0,covxi,covxi,varxi1),ncol=2)
sigmaxi <- array(0,dim=c(2,2,B))
xitilda <- matrix(0,B,2)
for(l in 1:B){
sigmaxi[,,l] <- solve(solve(sigmaxi0)+1/vare*ximatrix[,,l])
xitilda[l,] <- solve(vare*solve(sigmaxi0)+ximatrix[,,l])%*%tot3[l,]
}
library(mvtnorm)
xiexp <- matrix(0,B,2)
for(l in 1:B){
xiexp[l,] <- rmvnorm(1,mean=xitilda[l,],sigma=sigmaxi[,,l])
}
return(xiexp)
}
# function to generate latent variable
generatez <- function(uexp,beta,r0,r,c){
zmeanr <- isMatr%*%r+r0+x%*%beta
zmeanl <- isMatl%*%r+r0+x%*%beta
zexp <- rep(0,B)
for(i in 1:B){
if(delta1[i]==1){zexp[i] <- qnorm(runif(1)*pnorm(zmeanr[i]+c*uexp[i])+1-pnorm(zmeanr[i]+ c*uexp[i]))+zmeanr[i]+c*uexp[i]}
if(delta2[i]==1){zexp[i] <- qnorm(runif(1)*(pnorm(zmeanr[i]+c*uexp[i])-pnorm(zmeanl[i]+c*uexp[i]))+1-pnorm(zmeanr[i]+c*uexp[i]))+zmeanl[i]+c*uexp[i]}
if(delta3[i]==1){zexp[i] <- qnorm(runif(1)*(1-pnorm(zmeanl[i]+c*uexp[i])))+zmeanl[i]+c*uexp[i]}
}
return(zexp)
}
# function to generate parameters
generateparam <- function(uexp,xiexp,zexp,varu,varxi0,varxi1,covxi,vare,c,theta,beta,ita,r0,r){
au <- 0.1
bu <- 0.1
newvaru <- 1/rgamma(1,B/2+au,sum(uexp^2)/2+bu)
# introduce library for drawing sample from inverse Wishart distribution
library(LaplacesDemon)
S <- t(xiexp)%*%xiexp/B
M0 <- 3
V0 <- matrix(c(0.5,0,0,0.5),2,2)
newvarmatrix <- rinvwishart(B+M0,B*S+V0)
newvarxi0 <- newvarmatrix[1,1]
newcovxi <- newvarmatrix[1,2]
newvarxi1 <- newvarmatrix[2,2]
sumc <- zexp-talpha%*%r-r0-x%*%beta
vc <- 1
newc <- rnorm(1, sum(sumc*uexp)/(sum(uexp^2)+vc),sqrt(1/(sum(uexp^2)+vc)))
sumr0 <- zexp-talpha%*%r-x%*%beta-newc*uexp
m0 <- -6
v0 <- 1
newr0 <- rnorm(1,(sum(sumr0)+m0*v0)/(B+v0),sqrt(1/(B+v0)))
newr <- rep(0,jd2)
for(i in 1:jd2){
wl <- sum(talpha[,i]^2)
alambda <- 1
blambda <- 1
lambda <- rgamma(1,alambda+jd2,blambda+sum(r))
# introduce library for drawing sample from truncated normal distribution
library(truncnorm)
if(wl==0){newr[i] <- rexp(lambda)}
if(wl > 0){el <- 1/wl*(sum(talpha[,i]*(zexp-newr0-talpha[,-i]%*%r[-i]-x%*%beta-newc*uexp))-lambda); aaa=-zexp[which(delta2==1)]-isMatr[which(delta2==1),][,-i]%*%r[-i]+isMatl[which(delta2==1),][,-i]%*%r[-i]; bbb <- isMatr[which(delta2==1),][,i]-isMatl[which(delta2==1),][,i]; cl <- max(aaa[which(bbb>1e-8)]/bbb[which(bbb>1e-8)]); dl=max(cl,0); newr[i] <- rtruncnorm(1,dl,Inf,el,sqrt(1/wl))}
r[i] <- newr[i]
}
betavar <- 100*diag(p)
sigmabeta0 <- betavar
sigmabeta <- solve(solve(sigmabeta0)+t(x)%*%x)
betatilda <- sigmabeta%*%t(x)%*%(zexp-(talpha%*%newr+newr0)-newc*uexp)
# introduce library to draw sample from multivariate normal distribution
library(mvtnorm)
newbeta <- c(rmvnorm(1,mean=betatilda,sigma=sigmabeta))
itavar <- 100*diag(jd1)
sigmaita0 <- itavar
sigmaita <- solve(solve(sigmaita0)+1/vare*t(msMat)%*%msMat)
newsumita <- (y-rep(x%*%theta,nobs)-rep(uexp,nobs)-rep(xiexp[,1],nobs)-rep(xiexp[,2],nobs)*time)*msMat
itatilda <- sigmaita%*%colSums(newsumita)*1/vare
newita <- c(rmvnorm(1,mean=itatilda,sigma=sigmaita))
thetavar <- 100*diag(p)
sigmatheta0 <- thetavar
sigmatheta <- solve(solve(sigmatheta0)+1/vare*t(x)%*%(nobs*x))
longx <- matrix(0,length(y),p)
for(l in 1:p){
longx[,l] <- rep(x[,l],nobs) 
}
thetatilda <- sigmatheta%*%t(longx)%*%(y-msMat%*%newita-rep(uexp,nobs)-rep(xiexp[,1],nobs)-rep(xiexp[,2],nobs)*time)*1/vare
newtheta <- c(rmvnorm(1,mean=thetatilda,sigma=sigmatheta))
sume <- y-msMat%*%newita-rep(x%*%newtheta,nobs)-rep(uexp,nobs)-rep(xiexp[,1],nobs)-rep(xiexp[,2],nobs)*time
ae <- 0.1
be <- 0.1
newvare <- 1/rgamma(1,sum(nobs)/2+ae,sum(sume^2)/2+be)
return(list(newvaru,newvarxi0,newvarxi1,newcovxi,newvare,newc,newtheta,newbeta,newita,newr0,newr))
}
# function for Gibbs sampler
gibbs <- function(n,varu,varxi0,varxi1,covxi,vare,c,theta,beta,ita,r0,r){
expu <- matrix(0,B,n)
zexp <- matrix(0,B,n)
expxi <- array(0,dim=c(B,2,n))
newvaru <- rep(0,n)
newvarxi0 <- rep(0,n)
newvarxi1 <- rep(0,n)
newcovxi <- rep(0,n)
newvare <- rep(0,n)
newc <- rep(0,n)
newtheta <- matrix(0,n,p)
newbeta <- matrix(0,n,p)
newita <- matrix(0,n,jd1)
newr0 <- rep(0,n)
newr <- matrix(0,n,jd2)
uexp <- rnorm(B,0,sqrt(varu))
library(mvtnorm)
xiexp <- rmvnorm(B,mean=c(0,0),sigma=matrix(c(varxi0,covxi,covxi,varxi1),ncol=2))
for(l in 1:n){
zexp[,l] <- generatez(uexp,beta,r0,r,c)
param <- generateparam(uexp,xiexp,zexp[,l],varu,varxi0,varxi1,covxi,vare,c,theta,beta,ita,r0,r)
newvaru[l] <- param[[1]]
newvarxi0[l] <- param[[2]]
newvarxi1[l] <- param[[3]]
newcovxi[l] <- param[[4]]
newvare[l] <- param[[5]]
newc[l] <- param[[6]]
newtheta[l,] <- param[[7]]
newbeta[l,] <- param[[8]]
newita[l,] <- param[[9]]
newr0[l] <- param[[10]]
newr[l,] <- param[[11]]
varu <- newvaru[l]
varxi0 <- newvarxi0[l]
varxi1 <- newvarxi1[l]
covxi <- newcovxi[l]
vare <- newvare[l]
c <- newc[l]
theta <- newtheta[l,]
beta <- newbeta[l,]
ita <- newita[l,]
r0 <- newr0[l]
r <- newr[l,]
expxi[,,l] <- generatexi(uexp,varxi0,varxi1,covxi,vare,theta,ita)
expu[,l] <- generateu(expxi[,,l],zexp[,l],varu,vare,c,theta,beta,ita,r0,r)
xiexp <- expxi[,,l]
uexp <- expu[,l]
}
return(list(newvaru,newvarxi0,newvarxi1,newcovxi,newvare,newc,newtheta,newbeta,newita,newr0,newr,expu,expxi,zexp))
}
# run Gibbs sampler
parameter <- gibbs(n,varu,varxi0,varxi1,covxi,vare,c,theta,beta,ita,r0,r)
# results
varuvalue <- parameter[[1]][(burnin+1):n]
varxi0value <- parameter[[2]][(burnin+1):n]
varxi1value <- parameter[[3]][(burnin+1):n]
covxivalue <- parameter[[4]][(burnin+1):n]
varevalue <- parameter[[5]][(burnin+1):n]
cvalue <- parameter[[6]][(burnin+1):n]
thetavalue <- parameter[[7]][((burnin+1):n),]
betavalue <- parameter[[8]][((burnin+1):n),]
itavalue <- parameter[[9]][((burnin+1):n),]
r0value <- parameter[[10]][(burnin+1):n]
rvalue <- parameter[[11]][((burnin+1):n),]
postvaru <- mean(varuvalue)
lowervaru <- quantile(varuvalue, probs = c(0.025, 0.975))[[1]]
uppervaru <- quantile(varuvalue, probs = c(0.025, 0.975))[[2]]
sdvaru <- sd(varuvalue)
postvarxi0 <- mean(varxi0value)
lowervarxi0 <- quantile(varxi0value, probs = c(0.025, 0.975))[[1]]
uppervarxi0 <- quantile(varxi0value, probs = c(0.025, 0.975))[[2]]
sdvarxi0 <- sd(varxi0value)
postvarxi1 <- mean(varxi1value)
lowervarxi1 <- quantile(varxi1value, probs = c(0.025, 0.975))[[1]]
uppervarxi1 <- quantile(varxi1value, probs = c(0.025, 0.975))[[2]]
sdvarxi1 <- sd(varxi1value)
postcovxi <- mean(covxivalue)
lowercovxi <- quantile(covxivalue, probs = c(0.025, 0.975))[[1]]
uppercovxi <- quantile(covxivalue, probs = c(0.025, 0.975))[[2]]
sdcovxi <- sd(covxivalue)
postvare <- mean(varevalue)
lowervare <- quantile(varevalue, probs = c(0.025, 0.975))[[1]]
uppervare <- quantile(varevalue, probs = c(0.025, 0.975))[[2]]
sdvare <- sd(varevalue)
postc <- mean(cvalue)
lowerc <- quantile(cvalue, probs = c(0.025, 0.975))[[1]]
upperc <- quantile(cvalue, probs = c(0.025, 0.975))[[2]]
sdc <- sd(cvalue)
posttheta <- rep(0,p)
lowertheta <- rep(0,p)
uppertheta <- rep(0,p)
sdtheta <- rep(0,p)
postbeta <- rep(0,p)
lowerbeta <- rep(0,p)
upperbeta <- rep(0,p)
sdbeta <- rep(0,p)
postita <- rep(0,jd1)
lowerita <- rep(0,jd1)
upperita <- rep(0,jd1)
sdita <- rep(0,jd1)
postr0 <- mean(r0value)
lowerr0 <- quantile(r0value, probs = c(0.025, 0.975))[[1]]
upperr0 <- quantile(r0value, probs = c(0.025, 0.975))[[2]]
sdr0 <- sd(r0value)
postr <- rep(0,jd2)
lowerr <- rep(0,jd2)
upperr <- rep(0,jd2)
sdr <- rep(0,jd2)
for(i in 1:p){
posttheta[i] <- mean(thetavalue[,i])
lowertheta[i] <- quantile(thetavalue[,i], probs = c(0.025, 0.975))[[1]]
uppertheta[i] <- quantile(thetavalue[,i], probs = c(0.025, 0.975))[[2]]
sdtheta[i] <- sd(thetavalue[,i])
postbeta[i] <- mean(betavalue[,i])
lowerbeta[i] <- quantile(betavalue[,i], probs = c(0.025, 0.975))[[1]]
upperbeta[i] <- quantile(betavalue[,i], probs = c(0.025, 0.975))[[2]]
sdbeta[i] <- sd(betavalue[,i])
}
for(i in 1:jd1){
postita[i] <- mean(itavalue[,i])
lowerita[i] <- quantile(itavalue[,i], probs = c(0.025, 0.975))[[1]]
upperita[i] <- quantile(itavalue[,i], probs = c(0.025, 0.975))[[2]]
sdita[i] <- sd(itavalue[,i])
}
for(i in 1:jd2){
postr[i] <- mean(rvalue[,i])
lowerr[i] <- quantile(rvalue[,i], probs = c(0.025, 0.975))[[1]]
upperr[i] <- quantile(rvalue[,i], probs = c(0.025, 0.975))[[2]]
sdr[i] <- sd(rvalue[,i])
}
group1long <- msMat%*%postita+rep(posttheta%*%c(mean(x[,1]),1),length(y))
group0long <- msMat%*%postita+rep(posttheta%*%c(mean(x[,1]),0),length(y))
usediteration <- n-burnin
cdf1matrix <- matrix(0,B,usediteration)
for(m in 1:B){
for(l in 1:usediteration){
cdf1matrix[m,l] <- pnorm((r0value[l]+talpha[m,]%*%rvalue[l,]+betavalue[l,]%*%c(mean(x[,1]),1))/sqrt(1+cvalue[l]^2*varuvalue[l]^2))
}
}
surgroup1 <- rep(0,B)
for(m in 1:B){
surgroup1[m] <- 1-mean(cdf1matrix[m,])
}
cdf0matrix <- matrix(0,B,usediteration)
for(m in 1:B){
for(l in 1:usediteration){
cdf0matrix[m,l] <- pnorm((r0value[l]+talpha[m,]%*%rvalue[l,]+betavalue[l,]%*%c(mean(x[,1]),0))/sqrt(1+cvalue[l]^2*varuvalue[l]^2))
}
}
surgroup0 <- rep(0,B)
for(m in 1:B){
surgroup0[m] <- 1-mean(cdf0matrix[m,])
}
return(list(postvaru,postvarxi0,postvarxi1,postcovxi,postvare,postc,posttheta,postbeta,sdvaru,sdvarxi0,sdvarxi1,sdcovxi,sdvare,sdc,sdtheta,sdbeta,lowervaru,uppervaru,lowervarxi0,uppervarxi0,lowervarxi1,uppervarxi1,lowercovxi,uppercovxi,lowervare,uppervare,lowerc,upperc,lowertheta,uppertheta,lowerbeta,upperbeta,postita,sdita,lowerita,upperita,postr0,sdr0,lowerr0,upperr0,postr,sdr,lowerr,upperr,group1long,group0long,surgroup1,surgroup0))
}
# Create table function
createtable <- function(survival){
x <- matrix(cbind(survival[[6]],survival[[7]]),ncol=2)
p <- ncol(x)
parameters <- rep(0,(p+p+6))
for(l in 1:p){
parameters[l] <- paste0("theta",l)
}
for(l in (p+1):(p+p)){
parameters[l] <- paste0("beta",l)
}
parameters[(p+p+1):(p+p+6)] <- c("varu","varxi0","varxi1","covxi01","vare","c")
table <- data.frame(parameters,estimate=round(c(results[[7]],results[[8]],results[[1]],results[[2]],results[[3]],results[[4]],results[[5]],results[[6]]),3))
table$lowerestimate <- round(c(results[[29]],results[[31]],results[[17]],results[[19]],results[[21]],results[[23]],results[[25]],results[[27]]),3)
table$upperestimate <- round(c(results[[30]],results[[32]],results[[18]],results[[20]],results[[22]],results[[24]],results[[26]],results[[28]]),3)
# display table
return(table)
}



