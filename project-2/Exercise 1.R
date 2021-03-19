
library(mvtnorm)
library(ggplot2)
library(gtools)

# a)

n = 100 # dimension of x
mu = rep(0,n) # mean E(x_i)
sigma = 1 # variance Var(x_i)



corr_func <- function(i, j){ # returns correlation between x_i and x_j
  return (exp(-0.1*abs(i-j)))
}


Sigma = matrix(data=NA, nrow=n, ncol=n) # covariance matrix of x
for (i in seq(1,n)){
  for (j in seq(1,n)){
    Sigma[i,j] = corr_func(i,j)
  }
}


x = mvrnorm(n=1, mu=mu, Sigma = Sigma) # vector from multivariate normal distribution
plot(x)

image(Sigma)

# b)

L = chol(Sigma) # Cholesky decomposition of Sigma. Sigma = L L^T 
image(L)


# c)

z = mvrnorm(n=1, mu=rep(0,n), Sigma=diag(n)) # vector of multivariate standard normal
x_chol = L%*%z # x = Lz
plot(x_chol)



# d)

Q = solve(Sigma) # Q = Sigma^{-1}
image(Q)
L_Q = chol(Q) # Sigma^{-1} = L_Q L_Q^T
image(L_Q)



# e)

z = mvrnorm(n=1, mu=rep(0,n), Sigma=diag(n))
x_chol_Q = solve(L_Q, z)
plot(x_chol_Q)



# f)

perm_pos = permute(seq(1,n)) # permuted positions 

Sigma_perm = matrix(data=NA, nrow=n, ncol=n)
for (i in seq(1,n)){
  for (j in seq(1,n)){
    Sigma_perm[i,j] = corr_func(perm_pos[i], perm_pos[j])
  }
}

image(Sigma_perm)
x_perm = mvrnorm(n=1, mu=mu, Sigma = Sigma) # vector from multivariate normal distribution
plot(seq(1,n),x_perm)




L_perm = chol(Sigma_perm) # Cholesky decomposition of Sigma. Sigma = L L^T 
image(L_perm)


z = mvrnorm(n=1, mu=rep(0,n), Sigma=diag(n)) # vector of multivariate standard normal
x_chol_perm = L_perm%*%z # x = Lz
plot(x_chol_perm)

Q_perm = solve(Sigma_perm) # Q = Sigma^{-1}
image(Q_perm)
L_Q_perm = chol(Q_perm) # Sigma^{-1} = L_Q L_Q^T
image(L_Q_perm)


z = mvrnorm(n=1, mu=rep(0,n), Sigma=diag(n))
x_chol_Q_perm = solve(L_Q_perm, z)
plot(x_chol_Q_perm)









