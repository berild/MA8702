
library(mvtnorm)
library(ggplot2)
library(purrr)

# a)

mu = 0
sigma = 1
n = 100
x = mvrnorm(n=n, mu=mu, Sigma = sigma) # vector from multivariate normal


corr_func <- function(i, j){ # returns correlation between x_i and x_j
  return (exp(-0.1*abs(i-j)))
}


Sigma = matrix(data=NA, nrow=n, ncol=n)
for (i in seq(1,n)){
  for (j in seq(1,n)){
    Sigma[i,j] = corr_func(i,j)
  }
}

#ggplot(data.frame(x = seq(1,n),y=seq(1,n),z=Sigma)) + geom_contour_filled(aes(x=x,y=y,z=z))
image(Sigma)

# b)

L = chol(Sigma)
image(L)


# c)

z = mvrnorm(n=n, mu=rep(0,n), Sigma=diag(n))
x_chol = L%*%z
plot(x_chol)

# d)

Q = solve(Sigma)
image(Q)
L_Q = chol(Q)
image(L_Q)

# e)

z = mvrnorm(n=n, mu=rep(0,n), Sigma=diag(n))
x_chol_Q = solve(L_Q, z)
plot(x_chol_Q)

# f)

















