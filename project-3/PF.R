
library(MASS)
library(ggplot2)
library(ggpubr)
library(tictoc)


set.seed(10)
# Creating vectors, matrices and constants

mu_1 = c(10, 30, 10, -10) # mean for initial state vector
Sigma_1 = diag(c(10^2, 10^2, 5^2,5^2)) # covariance matrix for initial state vector
delta = 1/60

T_steps = 50 # number of time steps

A = diag(4) # forward function
A[1,3] = A[2,4] = delta

epsilon_mean = rep(0,4) # E(epsilon_(t+1))
epsilon_cov = diag(c(0.1^2, 0.1^2, 0.5^2, 0.5^2)) # Var(epsilon_(t+1))

# reading sensor data from file
sensorA = read.delim("sensorA.txt",header=FALSE)
sensorB = read.delim("sensorB.txt",header=FALSE)
y_vec = t(as.matrix(cbind(sensorA, sensorB))) # combining observations in a matrix

y_mean = rep(0,2) # assumed mean of gaussian noise of observations
y_cov =  diag(c(0.1^2, 0.1^2)) # assumed coviarance of  gaussian noise

#  computing inverse of y_cov used in get_l() function
iR =diag(1/diag(y_cov))

# h(x)
h = function(mu){
  return(c(atan(mu[1]/mu[2]),atan((40-mu[2])/(40-mu[1]))))
}

# function calculating unnormalized particle weights
get_l = function(particles,y){
  return(sapply(seq(ncol(particles)),function(x){exp(-1/2*t(y - h(particles[,x]))%*%iR%*%(y - h(particles[,x])))}))
}

# particle filter function
PF = function(B, T_steps){
  tic() # time start
  
  # pre assigning matricies
  mu_t = matrix(data=0, nrow=4, ncol=T_steps)
  cov_t = array(data = 0, dim= c(4,4,T_steps))
  particles = t(mvrnorm(n=B, mu = mu_1, Sigma = Sigma_1))
  lower_t = matrix(data = 0, nrow=4, ncol = T_steps)
  upper_t = matrix(data = 0, nrow=4, ncol = T_steps)
  
  for (t in 1:T_steps){
    l_t = get_l(particles,y_vec[,t]) # calc unnormalized particle weights
    w_t = l_t/sum(l_t) # normalizing weights
    
    # weighted resampling 
    particles = particles[,sample(seq(ncol(particles)),size = B,replace = TRUE, prob = w_t)]
    
    # emperical estimates
    mu_t[,t] = rowMeans(particles)
    cov_t[,,t] = rowSums(particles- mu_t[,t])%*%t(rowSums(particles - mu_t[,t]))/B

    # empirical upper and lower bounds
    lower_t[,t] = c(sort(particles[1,])[as.integer(B*0.05)],
                    sort(particles[2,])[as.integer(B*0.05)],
                    sort(particles[3,])[as.integer(B*0.05)],
                    sort(particles[4,])[as.integer(B*0.05)])
    upper_t[,t] = c(sort(particles[1,])[as.integer(B*0.95)],
                    sort(particles[2,])[as.integer(B*0.95)],
                    sort(particles[3,])[as.integer(B*0.95)],
                    sort(particles[4,])[as.integer(B*0.95)])

    # forcasting step
    particles = A%*%particles + t(mvrnorm(n=B, mu = epsilon_mean, Sigma = epsilon_cov))

    cat("t = ",t, "\n")

  }
  toc() # time end
  return(list(mu_t,cov_t,lower_t,upper_t))
}

B = 100 # number of particles
res = PF(B, T_steps) # running particle filter


# creating data frames for plotting
pos_df = data.frame("x_mean"=res[[1]][1,],
                    "y_mean"= res[[1]][2,],
                    "x_lower"=res[[3]][1,],
                    "y_lower"=res[[3]][2,],
                    "x_upper"=res[[4]][1,],
                    "y_upper"=res[[4]][2,])
x_vel_df = data.frame("x_vel_mean"=res[[1]][3,],
                      "x_lower"=res[[3]][3,],
                      "x_upper"=res[[4]][3,])
y_vel_df = data.frame("y_vel_mean"= res[[1]][4,],
                      "y_lower"=res[[3]][4,],
                      "y_upper"=res[[4]][4,])


# plotting trajectory
ggplot(pos_df)+
  geom_point(aes(x=x_mean, y=y_mean, color="mean"))+
  geom_point(aes(x=x_lower, y=y_lower, color="90 % confidence bounds"))+
  geom_point(aes(x=x_upper, y=y_upper))+
  coord_fixed(ratio = 1)+
  labs(x = "easting", y = "northing", color = "")+
  xlim(c(0,40))+
  ylim(c(0,40))+
  scale_color_manual(values=c("black","red"))+
  theme_bw() +
  theme(legend.position = "bottom")


# plotting velocities
x_vel = ggplot(x_vel_df)+
  geom_point(aes(x=1:T_steps, y=x_vel_mean, color="mean"))+
  geom_point(aes(x=1:T_steps, y=x_lower, color="90 % confidence bounds"))+
  geom_point(aes(x=1:T_steps, y=x_upper))+
  coord_fixed(ratio = 0.7)+
  labs(x = "time", y = "velocity", color = "",title = "a")+
  ylim(c(-5,25))+xlim(c(0,50))+
  theme_bw() +
  scale_color_manual(values=c("black","red")) + 
  theme(plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = 12), panel.background = element_blank())
x_vel

y_vel = ggplot(y_vel_df)+
  geom_point(aes(x=1:T_steps, y=y_vel_mean, color="mean"))+
  geom_point(aes(x=1:T_steps, y=y_lower, color="90 % confidence bounds"))+
  geom_point(aes(x=1:T_steps, y=y_upper))+
  coord_fixed(ratio = 0.6)+
  labs(x = "time", y = "velocity", color = "", title="b")+
  xlim(c(0,50))+ylim(c(-35,0))+
  theme_bw() +
  scale_color_manual(values=c("black","red")) + 
  theme(plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = 12), panel.background = element_blank())
y_vel

ggarrange(x_vel, y_vel, nrow=2,common.legend = T, legend="bottom" )

