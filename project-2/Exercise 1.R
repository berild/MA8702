
library(mvtnorm)
library(ggplot2)
library(gtools)
library(MASS)

set.seed(10)
# a)

n = 100 # dimension of x
mu = rep(0,n) # mean E(x_i) = 0
sigma = 1 # variance Var(x_i) = 1



corr_func <- function(i, j){ # returns correlation between x_i and x_j
  return (exp(-0.1*abs(i-j))) # Corr(x_i, x_j) = exp(-0.1*|i-j|)
}


Sigma = matrix(data=NA, nrow=n, ncol=n) # covariance matrix of x
for (i in seq(1,n)){
  for (j in seq(1,n)){
    Sigma[i,j] = corr_func(i,j)
  }
}

# plotting Sigma
grid = expand.grid(x=seq(from=0, to = n, length.out = n), y=-seq(from = 0, to = n, length.out = n))
df_Sigma = data.frame(grid, Sigma)
ggplot(df_Sigma, asp=1)+geom_point(aes(x=grid[,1], y=grid[,2], color=Sigma))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.print(pdf, "Sigma_1a.pdf")

# b)

L = t(chol(Sigma)) # Cholesky decomposition of Sigma. Sigma = L L^T 

# plotting L
df_L = data.frame(grid, L)
ggplot(df_L)+geom_point(aes(x=grid[,1], y=grid[,2], color=t(L)))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+scale_color_continuous("L")
dev.print(pdf, "L_1b.pdf")


# c)

# simulating x = L*z
z = mvrnorm(n=1, mu=rep(0,n), Sigma=diag(n)) # z ~ N(0,I)
x_chol = L%*%z # x = L*z ~ N(0, LL^T) = N(0, Sigma)

# plotting x
plot(x_chol, ylab="")
dev.print(pdf, "x_1c.pdf")


# d)

# computing Q
Q = solve(Sigma) # Q = Sigma^{-1}

# plotting Q
df_Q = data.frame(grid, Q)
ggplot(df_Q)+geom_point(aes(x=grid[,1], y=grid[,2], color=Q))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.print(pdf, "Q_1d.pdf")

# computing L_Q
L_Q = t(chol(Q)) # Sigma^{-1} = L_Q L_Q^T


# plotting L_Q
df_L_Q = data.frame(grid, L_Q)
ggplot(df_L_Q)+geom_point(aes(x=grid[,1], y=grid[,2], color=t(L_Q)))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+scale_color_continuous("L_Q")
dev.print(pdf, "L_Q_1d.pdf")


# e)

# sampling x through L_Q

x_chol_Q = solve((L_Q), z) # L_Q^T x = z. x ~ N(0, (L_Q*L_Q^T)^(-1)) = N(0, Q^(-1)) = N(0, Sigma)

# plotting x
plot(x_chol_Q, ylab="")
dev.print(pdf, "x_1e.pdf")


# f)


perm_pos = permute(seq(1,n)) # permuted positions 


# constructing permutation matrix "perm_matrix" such that perm_matrix*Sigma_perm*perm_matrix = Sigma
perm_matrix = matrix(0, nrow=n, ncol=n)
for (i in seq(1:n)){
  perm_matrix[perm_pos[i],i] = 1
}


# computing the covariance matrix Sigma where the positions are permuted, denoted Sigma_perm
Sigma_perm = matrix(data=NA, nrow=n, ncol=n) 
for (i in seq(1,n)){
  for (j in seq(1,n)){
    Sigma_perm[i,j] = corr_func(perm_pos[i], perm_pos[j])
  }
}


# plotting Sigma_perm
df_Sigma_perm = data.frame(grid, Sigma)
ggplot(df_Sigma_perm)+geom_point(aes(x=grid[,1], y=grid[,2], color=Sigma_perm))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+scale_color_continuous("Sigma")
dev.print(pdf, "Sigma_1f.pdf")


# computing L_perm satisfying Sigma_perm = L_perm*L_perm^T
L_perm = t(chol(Sigma_perm)) # Cholesky decomposition of Sigma. Sigma = L L^T 
df_L_perm = data.frame(grid, L_perm)
ggplot(df_L_perm)+geom_point(aes(x=grid[,1], y=grid[,2], color=t(L_perm)))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+scale_color_continuous("L")
dev.print(pdf, "L_1f.pdf")

# sampling x through Sigma_perm
x_chol_perm = perm_matrix%*%L_perm%*%z # x = L_perm z. x ~ N(0, Sigma_perm)

# plotting x
plot(x_chol_perm, ylab = "", xlab = "Index", ylim=c(-3,2))
dev.print(pdf, "x_Sigma_perm_1f.pdf")



# repeating the same procedure. This time sampling x using precision matrix Q_perm

Q_perm = solve(Sigma_perm) # Q_perm = Sigma_perm^{-1}

# plotting Q_perm 
df_Q_perm = data.frame(grid, Q_perm)
ggplot(df_L_perm)+geom_point(aes(x=grid[,1], y=grid[,2], color=Q_perm))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+scale_color_continuous("")
dev.print(pdf, "Q_perm_1f.pdf")

# computing and plotting L_Q_perm
L_Q_perm = t(chol(Q_perm)) # Sigma^{-1} = L_Q L_Q^T
df_L_Q_perm = data.frame(grid, L_Q_perm)
ggplot(df_L_Q_perm)+geom_point(aes(x=grid[,1], y=grid[,2], color=t(L_Q_perm)))+theme(axis.title.x = element_blank(), axis.title.y = element_blank())+scale_color_continuous("")
dev.print(pdf, "L_Q_perm_1f.pdf")

# sampling x through Q_perm

x_chol_Q_perm = perm_matrix%*%solve(t(L_Q_perm), z) # x ~ N(0, (L_Q*L_Q^T)^(-1)) = N(0, Q^(-1)) = N(0, Sigma)

# plotting x
plot(x_chol_Q_perm, ylab="")
dev.print(pdf, "x_Q_perm_1f.pdf")






