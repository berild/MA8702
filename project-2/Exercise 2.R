
# Exercise 2

# Task 2.1: Simulation


library(fields)
library(MASS)
library(ggplot2)
library(mvtnorm)

# defining covariance function and mean function

matern_cov <- function(sigma_sq, phi, dist_matrix){ 
  # returns the covariance given parameters sigma^2 and phi, and distance matrix, assuming Matérn covariance function
  rows = dim(dist_matrix)[1] # number of rows in distance matrix
  cols = dim(dist_matrix)[2] # number of columns in distance matrix
  Sigma = matrix(NA, nrow=rows, ncol=cols) # covariance matrix should have same dimension as distance matrix
  for (i in seq(1,rows)){ 
    for (j in seq(1,cols)){
      Sigma[i,j] = sigma_sq*(1+phi*dist_matrix[i,j])*exp(-phi*dist_matrix[i,j]) # Cov(x(s_i), x(s_j)) = sigma^2*(1+phi*h)*exp(-phi*h)
    }
  }
  return (Sigma)
}

# defining a function returning mean value at site (x_coord, y_coord) for a specific value alpha
mean_value <- function(alpha, x_coord, y_coord){ # returns mu_j, where coord is the coordinate for j
  return (alpha*((x_coord-0.5)+(y_coord-0.5))) # mu_j = alpha*((s_j1-0.5)+(s_j2-0.5))
}

# parameters
set.seed(10)
N = 200 # number of random observation sites
sigma_sq = 1 # sigma^2
phi = 10 # phi
tau_sq = 0.05^2 # phi^2
alpha = 1 # alpha


x_coords = runif(n=N) # sampling x-coordinates for the random sites
y_coords = runif(n=N) # sampling y-coordinates for the random sites
coords = cbind(x_coords, y_coords) # coordinates of the N random sites
plot(x=coords[,1], y=coords[,2]) # plot of the random sites


# Creating covariance matrix:

dist_matrix = rdist(coords) # matrix containing distance between every pair of points
Sigma = matern_cov(sigma_sq = sigma_sq, phi = phi, dist_matrix = dist_matrix) # defining Sigma through the Matérn covariance function

# sampling x

L = t(chol(Sigma)) # L is the Cholesky decomposition of the covariance matrix Sigma
z = rnorm(N) # z ~ N(0,1)
x = L%*%z # x ~ N(0, Sigma)

mu_vec = rep(0, N)
for (i in seq(1,N)){
  mu_vec[i] = mean_value(alpha = alpha, x_coord = coords[i,1], y_coord = coords[i,2])
}

x = x+mu_vec # x ~ N(mu, Sigma)

y = x+rnorm(n=N,sd=sqrt(tau_sq)) # y ~ N(mu, Sigma+tau*I_N)  
y = as.vector(y)

# plotting the random sites of the Gaussian field x
df_x = data.frame(coords, x)
ggplot(df_x)+geom_point(aes(x=coords[,1], y=coords[,2],color=x))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")
dev.print(pdf, "x_21.pdf")

# plotting the observations y
df_y = data.frame(coords, y)
ggplot(df_y)+geom_point(aes(x=coords[,1], y=coords[,2],color=y))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")
dev.print(pdf, "y_21.pdf")



# Task 2.2: Parameter estimation

# defining matrices in regression model. Using same notation as in lecture notes
C = Sigma+tau_sq*diag(N) # C = Sigma+tau^2* I_N
Q = solve(C) # Q = C^(-1)

X = coords[,1]+coords[,2]-rep(1,N) # X is vector where row j of X is s_j1+sj2-1


# creating the partial derivates of C

# returns dC/d(sigma^2) between the elements in x
C_deriv_sigma_sq <- function(phi, dist_matrix){ 
  rows = dim(dist_matrix)[1] # number of rows in distance matrix
  cols = dim(dist_matrix)[2] # number of columns in distance matrix
  C_sigma_sq = matrix(1, ncol=rows, nrow=cols) # dC/d(sigma^2) has same dimensions as distance matrix
  for (i in seq(1,rows)){
    for (j in seq(1,cols)){
      C_sigma_sq[i,j] = (1+phi*dist_matrix[i,j])*exp(-phi*dist_matrix[i,j]) # dC/ d(sigma^2) = (1+phi*h)*exp(-phi*h)
    }
  }
  return (C_sigma_sq)
}

# returns dC/d(phi)
C_deriv_phi <- function(sigma_sq, phi, dist_matrix){  
  rows = dim(dist_matrix)[1] # number of rows in distance matrix
  cols = dim(dist_matrix)[2] # number of columns in distance matrix
  C_phi= matrix(1, ncol=rows, nrow=cols) # dC/d(phi) has same dimensions as distance matrix
  for (i in seq(1,rows)){
    for (j in seq(1,cols)){
      C_phi[i,j] = -sigma_sq*phi*dist_matrix[i,j]^2*exp(-phi*dist_matrix[i,j]) # dC/d(phi) = -sigma^2*phi*h^2*exp(-phi*h)
    }
  }
  return (C_phi)
}

# returns a vector containing the partial derivates of C, i.e. dC/d(sigma^2), dC/d(tau^2) and dC/d(phi)
compute_C_partial_derivatives = function(theta_hat, dist_matrix){ # Computes dC/d(theta)
  C_partial_sigma_sq = C_deriv_sigma_sq(theta_hat[3], dist_matrix) # dC/ d(sigma^2)
  C_partial_tau_sq = diag(N)  # dC/d(tau^2) = I_N
  C_partial_phi = C_deriv_phi(theta_hat[1], theta_hat[3], dist_matrix) # dC/d(phi)
  C_partial = list(C_partial_sigma_sq, C_partial_tau_sq, C_partial_phi) # gathering all the partial derivatives in a vector, for practical purposes
  return (C_partial)
}

# returns the score vector (dl/d(sigma^2), dl/d(tau^2), dl/d(phi))
compute_score_vec = function(C_partial, Q, Z){ 
  score_vec = rep(0,3) # vector of length 3 containing the three partial derivatives of the log-likelihood
  for (i in c(1,2,3)){ # dl/d(theta_i) =(-1/2)*trace(Q*dC/d(theta_i))+(1/2)*Z^T*Q*dC/d(theta_i)*Q*Z     theta_i is either sigma^2, tau^2 or phi
    score_vec[i] = (-1/2)*sum(diag(Q%*%C_partial[[i]]))+(1/2)*t(Z)%*%Q%*%C_partial[[i]]%*%Q%*%Z
  }
  return (score_vec)
}

# returns the expected hessian E(d^2l/d(theta)^2)
compute_hessian = function(C_partial, Q){ # computing hessian d^2 l/d(theta^2)
  hessian = matrix(data=NA, ncol=3, nrow=3) # 3x3 matrix containing the hessian of the log-likelihood
  for (i in c(1,2,3)){
    for (j in c(1,2,3)){# d^2 l/(d(theta_i) d(theta_j)) = (-1/2)*trace(Q* (dC/d(theta_i))* Q* (dC/d(theta_j)))
      hessian[i,j] = -(1/2)*sum(diag(Q%*%C_partial[[i]]%*%Q%*%C_partial[[j]]))
    }
  }
  return (hessian)
}


# proposing start values of Fisher scoring iterations
sigma_sq_prop = 5
tau_sq_prop = 0.2^2 
phi_prop = 20
alpha_prop = 5


beta_hat = alpha_prop # beta_hat = alpha
theta_hat = c(sigma_sq_prop, tau_sq_prop, phi_prop) # theta_hat = (sigma^2, tau^2, phi)

prev_loglikelihood = -100
rho = 100 # difference between likelihood in one iteration and the previous
tol = 10^{-10} # stopping criterion for Fisher scoring iterations
iter = 1 # counter for the iterations in the Fisher scoring algorithm


# Fisher scoring algorithm:


while (rho > tol){ # stopping the iterations when rho is smaller than tol
  
  # Computing Q = Q(theta_p)
  Sigma_prop = matern_cov(theta_hat[1], theta_hat[3], dist_matrix) # computing proposed Sigma-matrix with proposed parameter values
  C_prop = Sigma_prop+theta_hat[2]*diag(N) # C = Sigma + tau^2 I_N
  Q_prop = solve(C_prop) # Q = C^(-1)
  
  
  # computing beta_hat
  beta_hat = solve(t(X)%*%Q_prop%*%X)%*%t(X)%*%Q_prop%*%y # beta_hat = (X^T Q X)^{-1} X^T Q y
  
  
  # computing partial derivatives of C and Z = y-X*beta_hat
  C_partial = compute_C_partial_derivatives(theta_hat, dist_matrix)
  
  Z_prop = y-X*as.vector(beta_hat) # Z = y-X*beta_hat
  
  
  # computing score vector dl/d(theta)
  score_vec = compute_score_vec(C_partial, Q_prop, Z_prop)
  
  
  # computing hessian d^2 l/d(theta)^2
  hessian = compute_hessian(C_partial, Q_prop)
  
  # theta_hat_(p+1) = theta_hat(p+1)-hessian^(-1)*score_vec
  theta_hat = theta_hat - solve(hessian)%*%score_vec
  
  
  loglikelihood = dmvnorm(x=y, mean=X*as.vector(beta_hat), sigma=C_prop, log=TRUE) # computing loglikelihood 
  rho = abs(loglikelihood-prev_loglikelihood) # rho = | likelihood_(p+1)-likelihood_p|
  prev_loglikelihood = loglikelihood # the new likelihood is the old likelihood in next iteration
  
  
  # Printing number of iterations, likelihood, theta_hat and beta_hat
  
  cat("Iterations:", iter,"\n")
  
  cat("loglikelihood:", loglikelihood, "\n")
  
  cat("theta_hat:", theta_hat, "\n")
  
  cat("beta_hat:", beta_hat, "\n", "\n")
  
  iter = iter+1 # incrementing iter
}



# Task 2.3: Kriging

# creating grid for kriging prediction

n = 25 # making n times n grid
grid = expand.grid(x=seq(from=0, to = 1, length.out = n), y=seq(from = 0, to = 1, length.out = n)) # making grid sites
grid = cbind(grid[,1],grid[,2]) # removing first column in grid matrix
X_0 = grid%*%c(1,1)-rep(1,n^2) # X_0 is a vector of length n^2 = 625. jth entry of X_0 is s_j1+s_j2-1


# creating distance matrix and covariance matrix for prediction sites

grid_distances = rdist(grid) # matrix containing distances between all pairs of grid nodes
C_0 = matern_cov(theta_hat[1],theta_hat[3], grid_distances) # covariance matrix for prediction sites


# creating distance matrix and covariance matrix between prediction sites and measurement sites

pred_measure_dist = rdist(grid, coords) # n^2 times N-matrix (625 x 200) containing distances between prediction and measurement sites

C_pred_measure = matern_cov(theta_hat[1], theta_hat[3],pred_measure_dist) # C_{0,.} Covariances between prediction and measurement sites (625 x 200)


# computing kriging mean vector

kriging_mean = X_0*as.vector(beta_hat)+C_pred_measure%*%solve(C)%*%(y-X*as.vector(beta_hat))   # E(Y_0|Y) = X_0*beta_hat+C_{0,.}*C^(-1)(Y-X*beta-hat) 


# plotting kriging mean

df_kriging_mean = data.frame(grid, kriging_mean)
ggplot(df_kriging_mean)+geom_point(aes(x=grid[,1], y=grid[,2], color=kriging_mean))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+scale_color_continuous("")
dev.print(pdf, "krig_mean_23.pdf")
ggplot(df_kriging_mean, aes(x=grid[,1], y=grid[,2], fill = kriging_mean)) + geom_tile()+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+scale_fill_continuous("") # plotting map
dev.print(pdf, "krig_mean_map_23.pdf")

# comparing against x and y

df_x = data.frame(coords, x)
ggplot(df_x)+geom_point(aes(x=coords[,1], y=coords[,2],color=x)) # plotting x

df_y = data.frame(coords, y)
ggplot(df_y)+geom_point(aes(x=coords[,1], y=coords[,2],color=y)) # plotting observations


# computing variance matrix for prediction sites

kriging_var = C_0-C_pred_measure%*%solve(C)%*%t(C_pred_measure) # Var(Y_0|Y) = C_0-C_{0,.}C^{-1}C_{0,.}^T
diag_kriging_var = diag(kriging_var) # making a matrix containing diagonal elements of Var(Y_0|Y)

# plotting the kriging variance
df_kriging_var = data.frame(grid, diag_kriging_var)
ggplot(df_kriging_var)+geom_point(aes(x=grid[,1], y=grid[,2], color=diag_kriging_var))+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+scale_color_continuous("") # plotting each kriging variance in each grid node
dev.print(pdf, "krig_var_23.pdf")
ggplot(df_kriging_var, aes(x=grid[,1], y=grid[,2], fill = diag_kriging_var)) + geom_tile()+coord_fixed(ratio = 1)+xlab("easting")+ylab("northing")+scale_fill_continuous("") # plotting the kriging variance as a map
dev.print(pdf, "krig_var_map_23.pdf")

