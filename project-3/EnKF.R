
library(MASS)
library(ggplot2)
library(ggpubr)


set.seed(10)
# Creating vectors, matrices and constants

mu_1 = c(10, 30, 10, -10) # mean for initial state vector
Sigma_1 = diag(c(10^2, 10^2, 5^2,5^2)) # covariance matrix for initial state vector
delta = 1/60 

T_steps = 50 # number of time steps

A = diag(4) # forward function
A[1,3] = A[2,4] = delta

epsilon_mean = rep(0,4) # E(epsilon_(t+1))
epsilon_var = diag(c(0.1^2, 0.1^2, 0.5^2, 0.5^2)) # Var(epsilon_(t+1))

# reading sensor data from file
sensorA = read.delim("sensorA.txt",header=FALSE)
sensorB = read.delim("sensorB.txt",header=FALSE)
y_vec = cbind(sensorA, sensorB)

y_mean = rep(0,2)
y_var =  diag(c(0.1^2, 0.1^2))

# Ensemble Kalman filter:


EnKF = function(B, T_steps){
  # creating matrices to containing the ensemble averages in each iteration, in addition to confidence bounds 
  mean_pos_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps) # row number t contains the ensemble average of E_t and N_t
  lower_bound_pos_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps) # ma
  upper_bound_pos_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps)
  
  mean_vel_t = matrix(data=0, nrow=2, ncol=T_steps) # # row t contains the ensemble average of v_t and u_t
  lower_bound_vel_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps)
  upper_bound_vel_ensemble_t = matrix(data=0, nrow=2, ncol=T_steps)
  
  # Initializing  ensemble
  ensemble = t(mvrnorm(n = B, mu = mu_1, Sigma = Sigma_1))
  
  
  for (t in 1:T_steps){
    
    # creating forecast data
    y_t_sampled = matrix(0, nrow=2, ncol=B)
    for (b in 1:B){
      y_t_sampled[,b] = c(atan(ensemble[1,b]/ensemble[2,b]), atan((40-ensemble[2,b])/(40-ensemble[1,b])))+mvrnorm(n=1, mu = y_mean, Sigma = y_var)
    }
    cov_mat_y = cov(t(y_t_sampled)) # empirical covariance matrix for y
    cov_mat_xy = cov(t(ensemble), t(y_t_sampled)) # empirical covariance matrix for x and y
    
    y_t = y_vec[t,] # y_t is the sensor data at time t
    for (b in 1:B){
      ensemble[,b] = ensemble[,b] + t(cov_mat_xy%*%solve(cov_mat_y)%*%t(y_t-y_t_sampled[,b])) # conditioning each ensemble member on the data
    }
    
    # moving the ensemble one step forward in time

    for (b in 1:B){
      ensemble[,b] = A%*%ensemble[,b]+mvrnorm(n=1, mu = epsilon_mean, Sigma = epsilon_var)
    }
   
   # storing the mean and confidence bounds in matrices
    mean_pos_ensemble_t[,t] = c(mean(ensemble[1,]), mean(ensemble[2,]))
    mean_vel_t[,t] = c(mean(ensemble[3,]), mean(ensemble[4,]))
    
    lower_bound_pos_ensemble_t[,t] = c(sort(ensemble[1,])[as.integer(B*0.05)], sort(ensemble[2,])[as.integer(B*0.05)])
    upper_bound_pos_ensemble_t[,t] = c(sort(ensemble[1,])[as.integer(B*0.95)], sort(ensemble[2,])[as.integer(B*0.95)])
    lower_bound_vel_ensemble_t[,t] = c(sort(ensemble[3,])[as.integer(B*0.05)], sort(ensemble[4,])[as.integer(B*0.05)])
    upper_bound_vel_ensemble_t[,t] = c(sort(ensemble[3,])[as.integer(B*0.95)], sort(ensemble[4,])[as.integer(B*0.95)])
    
    cat("t = ",t, "\n")
   
  }
  
  return (list("mean_pos"=mean_pos_ensemble_t, "lower_pos"=lower_bound_pos_ensemble_t, "upper_pos"=upper_bound_pos_ensemble_t, "mean_vel"=mean_vel_t, "lower_vel"=lower_bound_vel_ensemble_t, "upper_vel"=upper_bound_vel_ensemble_t))
}

B = 1000
EnKF_output = EnKF(B, T_steps) # running EnKF

# storing the output
mean_pos_ensemble_t = EnKF_output$mean_pos
lower_bound_pos_ensemble_t =  EnKF_output$lower_pos
upper_bound_pos_ensemble_t =  EnKF_output$upper_pos

mean_vel_ensemble_t = EnKF_output$mean_vel
lower_bound_vel_ensemble_t =  EnKF_output$lower_vel
upper_bound_vel_ensemble_t =  EnKF_output$upper_vel

# creating data frames for plotting
df_positions = data.frame("x_mean"=mean_pos_ensemble_t[1,],"y_mean"= mean_pos_ensemble_t[2,], "x_lower"=lower_bound_pos_ensemble_t[1,],"y_lower"=lower_bound_pos_ensemble_t[2,],"x_upper"=upper_bound_pos_ensemble_t[1,], "y_upper"=upper_bound_pos_ensemble_t[2,])
df_x_velocities = data.frame("x_vel_mean"=mean_vel_ensemble_t[1,], "x_lower"=lower_bound_vel_ensemble_t[1,],"x_upper"=upper_bound_vel_ensemble_t[1,])
df_y_velocities = data.frame("y_vel_mean"= mean_vel_ensemble_t[2,],"y_lower"=lower_bound_vel_ensemble_t[2,], "y_upper"=upper_bound_vel_ensemble_t[2,])


# plotting trajectory
ggplot(df_positions)+geom_point(aes(x=x_mean, y=y_mean, color="mean"))+geom_point(aes(x=x_lower, y=y_lower, color="95 % confidence bounds"))+geom_point(aes(x=x_upper, y=y_upper))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+xlim(c(0,40))+ylim(c(0,40))+scale_color_manual(values=c("black","red"))+theme(legend.position = "bottom")


# plotting velocities
x_vel = ggplot(df_x_velocities)+geom_point(aes(x=1:T_steps, y=x_vel_mean, color="mean"))+geom_point(aes(x=1:T_steps, y=x_lower, color="95 % confidence bounds"))+geom_point(aes(x=1:T_steps, y=x_upper))+coord_fixed(ratio = 1)+xlab("time")+ylab("velocity")+ylim(c(0,25))+xlim(c(0,50))+scale_color_manual(values=c("black","red"))
y_vel = ggplot(df_y_velocities)+geom_point(aes(x=1:T_steps, y=y_vel_mean, color="mean"))+geom_point(aes(x=1:T_steps, y=y_lower, color="95 % confidence bounds"))+geom_point(aes(x=1:T_steps, y=y_upper))+coord_fixed(ratio = 0.7)+xlab("time")+ylab("velocity")+xlim(c(0,50))+ylim(c(-35,0))+scale_color_manual(values=c("black","red"))

ggarrange(x_vel, y_vel, nrow=2)

