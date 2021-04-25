
library(MASS)
library(ggplot2)

download.file("folk.ntnu.no/joeid/MA8702/sensorA.txt", destfile="/Users/haakogry/OneDrive - NTNU/Doktorgrad/Videregående beregningskrevende statistiske metoder (MA8702)/Prosjekter/MA8702/project-3/sensorA.txt")
download.file("folk.ntnu.no/joeid/MA8702/sensorB.txt", destfile="/Users/haakogry/OneDrive - NTNU/Doktorgrad/Videregående beregningskrevende statistiske metoder (MA8702)/Prosjekter/MA8702/project-3/sensorB.txt")


# Creating vectors, matrices and constants

mu_1 = c(10, 30, 10, -10) # mean for initial state vector
Sigma_1 = diag(c(10^2, 10^2, 5^2,5^2)) # covariance matrix for initial state vector
delta = 1/60 
B = 1000 # ensemble members
T_steps = 50 # number of time steps

A = diag(4) # forward function
A[1,3] = A[2,4] = delta

epsilon_mean = rep(0,4) # E(epsilon_(t+1))
epsilon_var = diag(c(0.1^2, 0.1^2, 0.5^2, 0.5^2)) # Var(epsilon_(t+1))

sensorA = read.delim("sensorA.txt",header=FALSE)
sensorB = read.delim("sensorB.txt",header=FALSE)
y_vec = cbind(sensorA, sensorB)

y_mean = rep(0,2)
y_var =  diag(c(0.1^2, 0.1^2))

# Ensemble Kalman filter:


EnKF = function(B, T_steps){
# Initializing state space vector and ensemble
  mean_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps) # row t contains the ensemble average of E_t and N_t
  lower_bound_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps)
  upper_bound_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps)
  x_t = (mvrnorm(n=1, mu = mu_1, Sigma = Sigma_1))
  ensemble = t(mvrnorm(n = B, mu = mu_1, Sigma = Sigma_1))
  
  
  for (t in 1:T_steps){
    
    #y_t_sampled = cbind(rep(atan(x_t[1]/x_t[2]),B), rep(atan((40-x_t[2])/(40-x_t[1])),B))+mvrnorm(n=B, mu = y_mean, Sigma = y_var)
    y_t_sampled = matrix(0, nrow=2, ncol=B)
    for (b in 1:B){
      y_t_sampled[,b] = c(atan(ensemble[1,b]/ensemble[2,b]), atan((40-ensemble[2,b])/(40-ensemble[1,b])))+mvrnorm(n=1, mu = y_mean, Sigma = y_var)
    }
    cov_mat_y = cov(t(y_t_sampled)) # empirical covariance matrix
    cov_mat_xy = cov(t(ensemble), t(y_t_sampled))
    
    y_t = y_vec[t,]
    
    for (b in 1:B){
      ensemble[,b] = ensemble[,b] + t(cov_mat_xy%*%solve(cov_mat_y)%*%t(y_t-y_t_sampled[,b]))
    }
    
    # moving ensemble and state space vector one step forward in time
    
    for (b in 1:B){
      ensemble[,b] = A%*%ensemble[,b]+mvrnorm(n=1, mu = epsilon_mean, Sigma = epsilon_var)
    }
   
    x_t = A%*%x_t+mvrnorm(n=1, mu = epsilon_mean, Sigma = epsilon_var)
   
    mean_ensemble_t[,t] = c(mean(ensemble[1,]), mean(ensemble[2,]))
    #mean_ensemble_t[2,t] = mean(ensemble[2,])
    lower_bound_ensemble_t[,t] = c(sort(ensemble[1,])[as.integer(B*0.05)], sort(ensemble[2,])[as.integer(B*0.05)])
    upper_bound_ensemble_t[,t] = c(sort(ensemble[1,])[as.integer(B*0.95)], sort(ensemble[2,])[as.integer(B*0.95)])
    
    
    cat("t = ",t, "\n")
  }
  print(mean_ensemble_t)
  print(lower_bound_ensemble_t)
  print(upper_bound_ensemble_t)
  #df_mean_lower_upper = data.frame(mean_ensemble_t, lower_bound_ensemble_t, upper_bound_ensemble_t)
  #colnames(df_mean_lower_upper) =c("mean", "lower","upper")
  return (list("mean"=mean_ensemble_t, "lower"=lower_bound_ensemble_t, "upper"=upper_bound_ensemble_t))
}


EnKF_output = EnKF(1000, 50)

mean_ensemble_t = EnKF_output$mean
lower_bound_ensemble_t =  EnKF_output$lower
upper_bound_ensemble_t =  EnKF_output$upper


#df_lower_bound = data.frame(lower_bound_ensemble_t[1,], lower_bound_ensemble_t[2,])
#df_upper_bound = data.frame(upper_bound_ensemble_t[1,], upper_bound_ensemble_t[2,])
df_positions = data.frame("x_mean"=mean_ensemble_t[1,],"y_mean"= mean_ensemble_t[2,], "x_lower"=lower_bound_ensemble_t[1,],"y_lower"=lower_bound_ensemble_t[2,],"x_upper"=upper_bound_ensemble_t[1,], "y_upper"=upper_bound_ensemble_t[2,])
#df_positions_1 = data.frame("x"=c(mean_ensemble_t[1,],lower_bound_ensemble_t[1,],upper_bound_ensemble_t[1,]), "y"= c(mean_ensemble_t[2,],lower_bound_ensemble_t[2,],upper_bound_ensemble_t[2,]))


ggplot(df_positions)+geom_point(aes(x=x_mean, y=y_mean, color="mean"))+geom_point(aes(x=x_lower, y=y_lower, color="lower"))+geom_point(aes(x=x_upper, y=y_upper, color="upper"))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+xlim(c(0,40))+ylim(c(0,40))




ggplot(df_positions_1)+geom_point(aes(x=x, y=y))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+xlim(c(0,40))+ylim(c(0,40))
rlang::last_error()


ggplot(df_lower_bound)+geom_point(aes(x=lower_bound_ensemble_t[1,], y=lower_bound_ensemble_t[2,]))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+xlim(c(0,40))+ylim(c(0,40))
ggplot(df_upper_bound)+geom_point(aes(x=upper_bound_ensemble_t[1,], y=upper_bound_ensemble_t[2,]))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+xlim(c(0,40))+ylim(c(0,40))


plot(t(mean_ensemble_t))

