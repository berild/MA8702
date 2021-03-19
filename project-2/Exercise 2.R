
# Exercise 2

# Task 2.1


library(fields)

matern_cov <- function(sigma_sq, phi, h){ # returns the covariance between the two elements of x, assuming Matérn covariance function
  return (sigma_sq*(1+phi*h)*exp(-phi*h))
  
}

C_deriv_phi <- function(sigma_sq, phi, h){
  return (sigma_sq*(2*h+h^2)*exp(-phi*h))
}


mean_value <- function(alpha, x_coord, y_coord){ # returns mu_j, where coord is the coordinate for j
  return (alpha*((x_coord-0.5)+(y_coord-0.5)))
}

N = 200 # number of random sites
sigma_sq = 1
phi = 10
tau_sq = 0.05^2
alpha = 1

x_coords = runif(n=N) # sampling x-coordinates for the random sites
y_coords = runif(n=N) # sampling y-coordinates for the random sites
coords = cbind(x_coords, y_coords) # coordinates of the N random sites
plot(x=coords[,1], y=coords[,2]) # plot of the random sites


# making covariance matrix:

dist_matrix = rdist(coords) # matrix containing distance between every pair of points
Sigma = matrix(data=NA, nrow = N, ncol = N) # covariance matrix of x
Sigma = matern_cov(sigma_sq = sigma_sq, phi = phi, h = dist_matrix) # defining Sigma through the Matérn covariance function


L = chol(Sigma) # L is the cholesky decomposition of the covariance matrix Sigma
z = mvrnorm(n=1, mu=rep(0,N), Sigma=diag(N))
x = L%*%z

mu_vec = mapply(FUN=mean_value, x_coord=coords[,1], y_coord=coords[,2], alpha=1) # 
x = x+mu_vec

y = x+t(rmvnorm(n=1, mean=rep(0,N), sigma=tau_sq*diag(N))) # observations y equals Gaussian field x plus gaussian noise with variance tau^2
plot(x=coords[,1], y=coords[,2]) # plot of the random sites
plot(mu_vec)


df = data.frame(coords, x)
ggplot(df)+geom_point(aes(x=coords[,1], y=coords[,2],color=x))


# Task 2.2

C = Sigma+tau_sq*diag(N)
Q = solve(C)
X = matrix(cbind(rep(1,N), rep(1,N),rep(-1,N)),ncol=3,nrow=N)

sigma_sq_prop = 10
tau_sq_prop = 0.5
phi_prop = 2
alpha_prop = 2
iter = 15

for (i in seq(1, iter)){
  
  beta_prop = alpha_prop
    
  Z_prop = y-X%*%beta_prop
  Sigma_prop = matern_cov(sigma_sq_prop, tau_sq_prop, dist_matrix)
  C_prop = Sigma_prop+tau_sq_prop*diag(N)
  Q_prop = solve(C_prop)
  beta_hat = solve(t(X)%*%Q_prop%*%X)%*%t(X)%*%Q_prop%*%y
  score_vec = rep(0,3)
  hessian = matrix(data=NA, ncol=3, nrow=3)
  
  partial_sigma_sq = C/sigma_sq
  partial_tau_sq = diag(N)
  partial_deriv_phi = C_deriv_phi(sigma_sq, phi, dist_matrix)
  C_partial = c(partial_sigma_sq, partial_tau_sq, partial_deriv_phi)
  
  #score_vec[1] = (-1/(2*sigma_sq))*sum(diag(Q%*%C))+(1/(2*sigma_sq))*t(Z_prop)%*%Q%*%C%*%Q%*%Z_prop
  #score_vec[2] = (-1/2)*sum(diag(Q))+(1/2)*t(Z_prop)%*%Q%*%Q%*%Z_prop
  #score_vec[3] = (1/2)*sum(diag(Q%*%partial_deriv_phi))+(1/2)*t(C)%*%Q%*%partial_deriv_phi%*%Q%*%Z_prop
  
  for (i in c(1,2,3)){
    score_vec[i] = (-1/2)*sum(diag(Q%*%C_partial[i]))+(1/2)*t(Z_prop)%*%Q%*%C_partial[i]%*%Q%*%Z_prop
  }
  
  for (i in c(1,2,3)){
    for (j in c(1,2,3)){
      hessian[i,j] = -(1/2)*sum(diag(Q%*%C_partial[i]%*%Q%*%C_partial[j]))
    }
  }

  
  
  theta_hat = theta_hat + solve(hessian)%*%score_vec

}
X%*%beta

# Task 2.3

