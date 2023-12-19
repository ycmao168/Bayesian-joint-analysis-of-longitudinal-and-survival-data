jointanalysis <- function(id_long,time_long,dep_long,x_long,z_long=NULL,L_surv,R_surv,x_surv,censoring=NULL,spline.long.specification,spline.surv.specification,varu=NULL,xisigma=NULL,vare=NULL,c=NULL,ita=NULL,theta=NULL,r0=NULL,r=NULL,beta=NULL,n,burnin,grid=NULL){
# id_long is the ID corresponding to longitudinal response
# time_long is the time point when each longitudinal response is collected
# dep_long is longitudinal response
# x_long is the covariate matrix in the longitudinal submodel
# z_long is the design matrix associated with subject-specific random effects
# L_surv and R_surv are the lower and upper bound of an interval during which the event occurs
# x_surv is the covariate matrix in the survival submodel
# censoring includes left,interval, and right censoring
# spline.long.specification=c(lnumber, ldegree) includes both the number of knots, and the degree for M-spline in modeling longitudinal baseline mean function
# spline.surv.specification=c(snumber,sdegree) is the entry including both the number of knots and the degree for I-spline in modeling survival baseline mean function
# varu is the variance of u, where u is the shared frailty
# xisigma is the variance-covariance matrix of the subject-specific random effects
# vare is the variance of the measurement error
# c is a constant, used in survival submodel as the coefficient of u to balance the scale difference and control the sign of the statistical association between the two types of responses
# ita is the vector of M spline coefficients; r0 is the I-spline intercept, and r is the vector of I-spline coeffcients
# theta and beta are covariate coefficients in longitudinal and survival submodel, respectively.
# B below is the sample size
# n is the number of iterations in the Gibbs sampling, and burnin is the number of iterations not used in the summary.
# nobs is the number of longitudinal responses for each subject 
# model is fully described in Mao, et al. (2023+): Bayesian joint analysis of longitudinal and interval-censored failure time data



nobs <- as.data.frame(table(id_long))[,2]
L <- L_surv
R <- R_surv
if(is.null(z_long)){z_long <- cbind(1,time_long)}
if(is.null(censoring)){delta1 <- ifelse(L==0,1,0);delta2 <- ifelse(L!=0 & R!=Inf,1,0);delta3 <- ifelse(R==Inf,1,0)}
q <- ncol(z_long)
if(is.null(varu)){varu=0.5}
if(is.null(xisigma)){xisigma=0.5*diag(q)}
if(is.null(vare)){vare=1}
if(is.null(c)){c=1}
if(is.null(ita)){ita=rep(0.2,14)}
if(is.null(theta)){theta=rep(0.2,2)}
if(is.null(r0)){r0=-1}
if(is.null(r)){r=rep(0.2,14)}
if(is.null(beta)){beta=rep(0.2,2)}
B <- length(nobs)
p_long <- ncol(x_long)
p_surv <- ncol(x_surv)
# jd1 and jd2 are the number of spline functions in the longitudinal and survival submodel, respectively
jd1 <- sum(spline.long.specification)-1
jd2 <- sum(spline.surv.specification)-1
# generate spline basis functions
union <- c(L,R)
idwithbound <- rep(1:B,2)[which(union!=0 & union!=Inf)]
time_surv <- union[which(union!=0 & union!=Inf)]
knots_long <- quantile(time_long,prob=seq(0,1,1/(spline.long.specification[1]-1)))[-1][-(spline.long.specification[1]-1)]
knots_surv <- quantile(time_surv,prob=seq(0,1,1/(spline.surv.specification[1]-1)))[-1][-(spline.surv.specification[1]-1)]
library(splines2)
msMat <- mSpline(time_long,knots=knots_long,degree=spline.long.specification[2],intercept=TRUE)
isMat <- iSpline(time_surv,knots=knots_surv,degree=spline.surv.specification[2],intercept=TRUE)
if(is.null(grid)){grid <- seq(quantile(time_long,0.05),quantile(time_long,0.95),length=50)}
lg <- length(grid)
longspline <- mSpline(grid,knots=knots_long,degree=spline.long.specification[2],intercept=TRUE)
survspline <- iSpline(grid,knots=knots_surv,degree=spline.surv.specification[2],intercept=TRUE)
isMatl <- matrix(NA,B,jd2)
isMatr <- matrix(NA,B,jd2)
for(i in 1:length(idwithbound)){
if(idwithbound[i+1] < idwithbound[i])
{k=i+1;
break}
}
for(l in 1:k-1){
isMatl[idwithbound[l],] <- isMat[l,]
}
for(m in k:length(idwithbound)){
isMatr[idwithbound[m],] <- isMat[m,]
}
talpha <- matrix(0,B,jd2)
for(l in 1:B){
if(delta1[l]==1){talpha[l,]<-isMatr[l,]}else{talpha[l,]<-isMatl[l,]}}
# function to generate random effects
generateu <- function(xiexp,zexp,varu,vare,c,theta,beta,ita,r0,r){
sumu <- dep_long-msMat%*%ita-x_long%*%theta-rowSums(apply(xiexp,2,function(c) rep(c,nobs))*z_long)
dat1 <- data.frame(cbind(id_long,sumu))
totl <- aggregate(sumu~id_long,dat1,sum)[,2]
tots <- zexp-talpha%*%r-r0-x_surv%*%beta
bb <- totl/vare+tots*c
aa <- (nobs*varu+c^2*vare*varu+vare)/(2*vare*varu)
expu <- rnorm(B,bb/(2*aa),sqrt(1/(2*aa)))
return(expu)
}
# function to generate random effects over time
generatexi <- function(uexp,xisigma,vare,theta,ita){
subz <- lapply(split(z_long,rep(c(1:B),nobs)),matrix,ncol=q)
sumxi <- dep_long-msMat%*%ita-rep(uexp,nobs)-x_long%*%theta
subsumxi <- lapply(split(sumxi,rep(c(1:B),nobs)),matrix,ncol=1)
sigmaxi <- array(0,dim=c(q,q,B))
xitilda <- matrix(0,B,q)
xiexp <- matrix(0,B,q)
for(l in 1:B){
sigmaxi <- solve(solve(xisigma)+1/vare*t(subz[[l]])%*%subz[[l]])
xitilda <- solve(vare*solve(xisigma)+t(subz[[l]])%*%subz[[l]])%*%(t(subz[[l]])%*%subsumxi[[l]])
library(mvtnorm)
xiexp[l,] <- rmvnorm(1,mean=xitilda,sigma=sigmaxi)
}
return(xiexp)
}
# function to generate latent variable
generatez <- function(uexp,beta,r0,r,c){
zmeanr <- isMatr%*%r+r0+x_surv%*%beta
zmeanl <- isMatl%*%r+r0+x_surv%*%beta
zexp <- rep(0,B)
index1=(delta1==1)
zexp[index1]<- qnorm(runif(sum(index1))*pnorm(zmeanr[index1]+c*uexp[index1])+1-pnorm(zmeanr[index1]+ c*uexp[index1]))+zmeanr[index1]+c*uexp[index1]
index2=(delta2==1)
zexp[index2] <- qnorm(runif(sum(index2))*(pnorm(zmeanr[index2]+c*uexp[index2])-pnorm(zmeanl[index2]+c*uexp[index2]))+1-pnorm(zmeanr[index2]+c*uexp[index2]))+zmeanl[index2]+c*uexp[index2]
index3=(delta3==1)
zexp[index3] <- qnorm(runif(sum(index3))*(1-pnorm(zmeanl[index3]+c*uexp[index3])))+zmeanl[index3]+c*uexp[index3]
return(zexp)
}
# function to generate parameters
generateparam <- function(uexp,xiexp,zexp,varu,xisigma,vare,c,theta,beta,ita,r0,r){
au <- 0.1
bu <- 0.1
newvaru <- 1/rgamma(1,B/2+au,sum(uexp^2)/2+bu)
# introduce library for drawing sample from inverse Wishart distribution
library(LaplacesDemon)
S <- t(xiexp)%*%xiexp/B
M0 <- q+1
V0 <- 0.5*diag(q)
newxisigma <- rinvwishart(B+M0,B*S+V0)
sumc <- zexp-talpha%*%r-r0-x_surv%*%beta
vc <- 1
newc <- rnorm(1, sum(sumc*uexp)/(sum(uexp^2)+vc),sqrt(1/(sum(uexp^2)+vc)))
sumr0 <- zexp-talpha%*%r-x_surv%*%beta-newc*uexp
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
if(wl==0){newr[i] <- rexp(1,lambda)}
if(wl > 0){el <- 1/wl*(sum(talpha[,i]*(zexp-newr0-talpha[,-i]%*%r[-i]-x_surv%*%beta-newc*uexp))-lambda); aaa=-zexp[which(delta2==1)]-isMatr[which(delta2==1),][,-i]%*%r[-i]+isMatl[which(delta2==1),][,-i]%*%r[-i]; bbb <- isMatr[which(delta2==1),][,i]-isMatl[which(delta2==1),][,i]; cl <- max(aaa[which(bbb>1e-8)]/bbb[which(bbb>1e-8)]); dl=max(cl,0); newr[i] <- rtruncnorm(1,dl,Inf,el,sqrt(1/wl))}
r[i] <- newr[i]
}
betavar <- 100*diag(p_surv)
sigmabeta0 <- betavar
sigmabeta <- solve(solve(sigmabeta0)+t(x_surv)%*%x_surv)
betatilda <- sigmabeta%*%t(x_surv)%*%(zexp-(talpha%*%newr+newr0)-newc*uexp)
# introduce library to draw sample from multivariate normal distribution
library(mvtnorm)
newbeta <- c(rmvnorm(1,mean=betatilda,sigma=sigmabeta))
itavar <- 100*diag(jd1)
sigmaita0 <- itavar
sigmaita <- solve(solve(sigmaita0)+1/vare*t(msMat)%*%msMat)
newsumita <- (dep_long-as.vector(x_long%*%theta)-rep(uexp,nobs)-rowSums(apply(xiexp,2,function(c) rep(c,nobs))*z_long))*msMat
itatilda <- sigmaita%*%colSums(newsumita)*1/vare
newita <- c(rmvnorm(1,mean=itatilda,sigma=sigmaita))
thetavar <- 100*diag(p_long)
sigmatheta0 <- thetavar
sigmatheta <- solve(solve(sigmatheta0)+1/vare*t(x_long)%*%(x_long))
thetatilda <- sigmatheta%*%t(x_long)%*%(dep_long-msMat%*%newita-rep(uexp,nobs)-rowSums(apply(xiexp,2,function(c) rep(c,nobs))*z_long))*1/vare
newtheta <- c(rmvnorm(1,mean=thetatilda,sigma=sigmatheta))
sume <- dep_long-msMat%*%newita-x_long%*%newtheta-rep(uexp,nobs)-rowSums(apply(xiexp,2,function(c) rep(c,nobs))*z_long)
ae <- 0.1
be <- 0.1
newvare <- 1/rgamma(1,sum(nobs)/2+ae,sum(sume^2)/2+be)
return(list(newvaru,newxisigma,newvare,newc,newtheta,newbeta,newita,newr0,newr))
}
gibbs <- function(n,varu,xisigma,vare,c,theta,beta,ita,r0,r){
expu <- matrix(0,B,n)
zexp <- matrix(0,B,n)
expxi <- array(0,dim=c(B,q,n))
newvaru <- rep(0,n)
newxisigma <- array(0,dim=c(q,q,n))
newvare <- rep(0,n)
newc <- rep(0,n)
newtheta <- matrix(0,n,p_long)
newbeta <- matrix(0,n,p_surv)
newita <- matrix(0,n,jd1)
newr0 <- rep(0,n)
newr <- matrix(0,n,jd2)
uexp <- rnorm(B,0,sqrt(varu))
library(mvtnorm)
xiexp <- rmvnorm(B,mean=rep(0,q),sigma=xisigma)
for(l in 1:n){
zexp[,l] <- generatez(uexp,beta,r0,r,c)
param <- generateparam(uexp,xiexp,zexp[,l],varu,xisigma,vare,c,theta,beta,ita,r0,r)
newvaru[l] <- param[[1]]
newxisigma[,,l] <- param[[2]]
newvare[l] <- param[[3]]
newc[l] <- param[[4]]
newtheta[l,] <- param[[5]]
newbeta[l,] <- param[[6]]
newita[l,] <- param[[7]]
newr0[l] <- param[[8]]
newr[l,] <- param[[9]]
varu <- newvaru[l]
xisigma <- newxisigma[,,l]
vare <- newvare[l]
c <- newc[l]
theta <- newtheta[l,]
beta <- newbeta[l,]
ita <- newita[l,]
r0 <- newr0[l]
r <- newr[l,]
expxi[,,l] <- generatexi(uexp,xisigma,vare,theta,ita)
expu[,l] <- generateu(expxi[,,l],zexp[,l],varu,vare,c,theta,beta,ita,r0,r)
xiexp <- expxi[,,l]
uexp <- expu[,l]
}
return(list(newvaru,newxisigma,newvare,newc,newtheta,newbeta,newita,newr0,newr,expu,expxi,zexp))
}
# run Gibbs sampler
parameter <- gibbs(n,varu,xisigma,vare,c,theta,beta,ita,r0,r)
varuvalue <- parameter[[1]][(burnin+1):n]
xisigmavalue <- parameter[[2]][,,(burnin+1):n]
varevalue <- parameter[[3]][(burnin+1):n]
cvalue <- parameter[[4]][(burnin+1):n]
thetavalue <- parameter[[5]][((burnin+1):n),]
betavalue <- parameter[[6]][((burnin+1):n),]
itavalue <- parameter[[7]][((burnin+1):n),]
r0value <- parameter[[8]][(burnin+1):n]
rvalue <- parameter[[9]][((burnin+1):n),]
estimate <- function(value){
post <- mean(value)
lower <- quantile(value, probs = c(0.025, 0.975))[[1]]
upper <- quantile(value, probs = c(0.025, 0.975))[[2]]
return(c(post,lower,upper))
}
estimatecoefficient <- function(value){
estimated <- matrix(0,p_surv,3)
for(l in 1:p_surv){
estimated[l,] <- estimate(value[,l])
}
return(estimated)
}
estimatevaru <- estimate(varuvalue)
estimatevare <- estimate(varevalue)
estimatec <- estimate(cvalue)
estimatetheta <- estimatecoefficient(thetavalue)
estimatebeta <- estimatecoefficient(betavalue)
meanxivar <- matrix(0,q,q)
xivar025 <- matrix(0,q,q)
xivar975 <- matrix(0,q,q)
varname <- matrix(0,q,q)
for(l in 1:q){
for(m in 1:q){
meanxivar[l,m] <- mean(xisigmavalue[l,m,])
xivar025[l,m] <- quantile(xisigmavalue[l,m,],0.025)
xivar975[l,m] <- quantile(xisigmavalue[l,m,],0.975)
varname[l,m] <- paste0("covxi",m-1,l-1)
}
}
for(l in 1:q){
varname[l,l] <- paste0("varxi",l-1,l-1)
}
meanxivarvector <- meanxivar[lower.tri(meanxivar,diag=TRUE)]
xivar025vector <- xivar025[lower.tri(xivar025,diag=TRUE)]
xivar975vector <- xivar975[lower.tri(xivar975,diag=TRUE)]
varnamevector <- varname[lower.tri(varname,diag=TRUE)]
estimatexisigma <- matrix(cbind(meanxivarvector,xivar025vector,xivar975vector),ncol=3)
parameters <- rep(0,(p_long+p_surv+3+q*(q+1)/2))
for(l in 1:p_long){
parameters[l] <- paste0("theta",l)
}
for(l in (p_long+1):(p_long+p_surv)){
parameters[l] <- paste0("beta",l-p_long)
}
parameters[(p_long+p_surv+1):(p_long+p_surv+3)] <- c("varu","vare","c")
parameters[(p_long+p_surv+4):(p_long+p_surv+3+q*(q+1)/2)] <- varnamevector
table <- rbind(estimatetheta,estimatebeta,estimatevaru,estimatevare,estimatec,estimatexisigma)
newtable <- data.frame(parameters,estimate=round(table[,1],3))
newtable$percentile_2.5th <- round(table[,2],3)
newtable$percentile_97.5th <- round(table[,3],3)
usediteration <- n-burnin
ygrid1 <- longspline%*%t(itavalue)+t(matrix(rep(thetavalue%*%c(mean(x_surv[,1]),1),lg),ncol=lg))
ygrid0 <- longspline%*%t(itavalue)+t(matrix(rep(thetavalue%*%c(mean(x_surv[,1]),0),lg),ncol=lg))
meanygrid1 <- rowMeans(ygrid1)
meanygrid0 <- rowMeans(ygrid0)
margicdf1matrix <- pnorm((t(matrix(rep(r0value,lg),ncol=lg))+survspline%*%t(rvalue)+t(matrix(rep(betavalue%*%c(mean(x_surv[,1]),1),lg),ncol=lg)))/sqrt(t(matrix(rep(rep(1,usediteration)+cvalue^2*varuvalue^2,lg),ncol=lg))))
margicdf0matrix <- pnorm((t(matrix(rep(r0value,lg),ncol=lg))+survspline%*%t(rvalue)+t(matrix(rep(betavalue%*%c(mean(x_surv[,1]),0),lg),ncol=lg)))/sqrt(t(matrix(rep(rep(1,usediteration)+cvalue^2*varuvalue^2,lg),ncol=lg))))
condicdf1matrix <- pnorm((t(matrix(rep(r0value,lg),ncol=lg))+survspline%*%t(rvalue)+t(matrix(rep(betavalue%*%c(mean(x_surv[,1]),1),lg),ncol=lg))))
condicdf0matrix <- pnorm((t(matrix(rep(r0value,lg),ncol=lg))+survspline%*%t(rvalue)+t(matrix(rep(betavalue%*%c(mean(x_surv[,1]),0),lg),ncol=lg))))
survgroup1margi <- 1-rowMeans(margicdf1matrix)
survgroup0margi <- 1-rowMeans(margicdf0matrix)
survgroup1condi <- 1-rowMeans(condicdf1matrix)
survgroup0condi <- 1-rowMeans(condicdf0matrix)
longfigure <- data.frame(grid,meanygrid1,meanygrid0)
survfigure <- data.frame(grid,survgroup1margi,survgroup0margi,survgroup1condi,survgroup0condi)
return(list(newtable,longfigure,survfigure))
}
