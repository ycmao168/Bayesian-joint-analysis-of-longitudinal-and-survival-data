source("Supportfunctions.R")

# input data
load("sampledata.RData")

results <- jointanalysis(id_long=longitudinal$subject,time_long=longitudinal$time,dep_long=longitudinal$y,x_long=cbind(longitudinal$longx1,longitudinal$longx2),z_long=cbind(1,longitudinal$time),L_surv=survival$L,R_surv=survival$R,x_surv=cbind(survival$x1,survival$x2),spline.long.specification=c(12,3),spline.surv.specification=c(12,3),n=500,burnin=100)

# Present summary results for all parameters including two sets of regression parameters and all variance parameters
# Summarized results include the posterior means and 95% credible intervals
summary <- results[[1]]
summary

# A summary example is as follows:
#    parameters  estimate    percentile_2.5th  percentile_97.5th
#1      theta1   -0.980           -1.075            -0.884
#2      theta2    1.014            0.837             1.190
#3       beta1    1.028            0.803             1.265
#4       beta2   -1.096           -1.449            -0.764
#5        varu    0.276            0.128             0.474
#6        vare    1.064            0.986             1.151
#7           c    1.131            0.432             2.069
#8     varxi00    0.194            0.063             0.350
#9     covxi01    0.041           -0.040             0.107
#10    varxi11    0.203            0.143             0.279


# Create a figure to plot the estimated mean longitudinal response over time for two different subgroups. 
longfigure <- results[[2]]
plot(longfigure[,1],longfigure[,2],xlab="Time",ylab="Estimated Mean Response",ylim=c(-5,4),type="l",lty=1)
lines(longfigure[,1],longfigure[,3],type="l",lty=2)
legend("topright",legend=c("Group 1","Group 0"),lty=c(1,2))

# Create a figure to plot the estimated marginal survival functions of the failure time for two different subgroups.
survfigure <- results[[3]]
plot(survfigure[,1],sort(survfigure[,2],decreasing=TRUE),xlab="Time",ylab="Estimated Survival Function",xlim=c(0,3),type="l",lty=1)
lines(survfigure[,1],sort(survfigure[,3],decreasing=TRUE),type="l",lty=2)
legend("topright",legend=c("Group 1","Group 0"),lty=c(1,2))

# Create a figure to plot the estimated conditional survival functions of the failure time given frailty=0 for two different subgroups.
plot(survfigure[,1],sort(survfigure[,4],decreasing=TRUE),xlab="Time",ylab="Estimated Survival Function",xlim=c(0,3),type="l",lty=1)
lines(survfigure[,1],sort(survfigure[,5],decreasing=TRUE),type="l",lty=2)
legend("topright",legend=c("Group 1","Group 0"),lty=c(1,2))

