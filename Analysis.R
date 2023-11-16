source("Supportfunctions.R")

# input data
load("sampledata.RData")

results <- jointanalysis(longitudinal,survival,12,12,3,0.5,0.5,0.5,0,1,1,rep(0.2,14),rep(0.2,2),-1,rep(0.2,14),rep(0.2,2),1000,5000)


# Create tables
table <- createtable(survival)
table

# The results is a table containing estimates, and 95% credible intervals for the parameters 


# Create Figure 1: the estimated mean response over time by two groups
y <- longitudinal[[1]]
time <- longitudinal[[2]]
group1long <- results[[45]]
group0long <- results[[46]]
newdataframe <- data.frame(time,group1long,group0long)
library(dplyr)
newdataframe <- newdataframe %>% arrange(time)
plot(newdataframe$time,newdataframe$group1long,xlab="Time",ylab="Estimated Mean Response", ylim=c((-2+min(y)),(max(y)+2)),type="l",lty=1)
lines(newdataframe$time,newdataframe$group0long,type="l",lty=2)
legend("topright",legend=c("Group 1","Group 0"),lty=c(1,2))



# The figure 1 shows estimated mean response over time for subjects who had mean of the standard normal variable in group 1 and group 0, respectively.

# Create Figure 2: the estimated survival time by two groups
delta1 <- survival[[1]]
delta2 <- survival[[2]]
delta3 <- survival[[3]]
L <- survival[[4]]
R <- survival[[5]]
nobs <- survival[[9]]
B <- length(nobs)
survtime <- rep(0,B)
for(l in 1:B){
if(delta1[l]==1){survtime[l] <- R[l]}
if(delta2[l]==1){survtime[l] <- L[l]}
if(delta3[l]==1){survtime[l] <- L[l]}
}
surgroup1 <- results[[47]]
surgroup0 <- results[[48]]
plot(sort(survtime),sort(surgroup1,decreasing=TRUE),xlab="Time",ylab="Estimated Survival Function", type="l",lty=1)
lines(sort(survtime),sort(surgroup0,decreasing=TRUE),type="l",lty=2)
legend("topright",legend=c("Group 1","Group 0"),lty=c(1,2))


