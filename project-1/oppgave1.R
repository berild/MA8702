
library(ggplot2)

dStNorm <- function(x,log=TRUE){
  mvtnorm::dmvnorm(x,mean = rep(0,2),sigma = matrix(c(1,0.9,0.9,1),nrow=2,ncol=2),log = log)
}

gradient.dStNorm <-function(x){
  solve(matrix(c(1,0.9,0.9,1),nrow=2,ncol=2),x)
}

dMMod <- function(x,log=TRUE){
  means = matrix(c(-1.5,-1.5,1.5,1.5,-2,2),nrow = 2, ncol = 3)
  vars = c(1,1,0.8)
  weights = rep(1/3,3)
  res = 0 
  for (i  in seq(3)){
    res = res + weights[i]*mvtnorm::dmvnorm(x,mean = means[,i], sigma = diag(vars[i],nrow=2))
  }
  if (log){
    return(log(res))
  }else{
    return(res)
  }
}

gradient.dMMod <-function(x){
  means = matrix(c(-1.5,-1.5,1.5,1.5,-2,2),nrow = 2, ncol = 3)
  vars = c(1,1,0.8)
  weights = rep(1/3,3)
  res = 0 
  like = dMMod(x,log=FALSE)
  for (i in seq(3)){
    res = res + weights[i]*solve(diag(vars[i],nrow=2),means[,i]-x)*mvtnorm::dmvnorm(x,mean = means[,i], sigma = diag(vars[i],nrow=2))
  }
  return(res/like)
}

dVolcano <- function(x,log =TRUE){
  x = as.vector(x)
  if (log){
    log(1/(2*pi)*exp(-1/2*t(x)%*%x)*(t(x)%*%x + 0.25))
  }else
    1/(2*pi)*exp(-1/2*t(x)%*%x)*(t(x)%*%x + 0.25)
}

gradient.dVolcano <- function(x){
  x = as.vector(x)
  return(as.vector(-x + 2*x/(t(x)%*%x+0.25)))
}

plotdensity <- function(f){
  x = as.matrix(expand.grid(seq(-5,5,0.1),seq(-5,5,0.1)))
  y = numeric(nrow(x))
  for (i in seq(nrow(x))){
    y[i] = f(x[i,])
  }
  ggplot(data.frame(x = x[,1],y=x[,2],z=y))  + 
    geom_contour_filled(aes(x=x,y=y,z=z))
}


plotdensity(dStNorm)

plotdensity(dMMod)

plotdensity(dVolcano)


dProp <- function(x,y,sigma = rep(0.5,2),log=TRUE){
  dnorm(x,mean = y, sd = sigma,log=log)
}

rProp <- function(y,sigma = rep(0.5,2),log=TRUE){
  rnorm(2, mean=y,sd = sigma)
}

MH <- function(init,dTar,dProp,rProp,nsamples=10000){
  eta = matrix(data = NA,nrow = nsamples, ncol = 2)
  acc.vec = numeric(nsamples)
  acc.prob = numeric(nsamples)
  eta[1,] = init$mu
  starttime = Sys.time()
  res = list()
  pb <- txtProgressBar(min = 0, max = nsamples, style = 3)
  for (i in seq(2,nsamples)){
    setTxtProgressBar(pb, i)
    eta.new = rProp(y = eta[i-1,],init$sigma)
    acc.prob[i] = min(1,exp(dTar(eta.new)  - dTar(eta[i-1,])))
    if (runif(1) < acc.prob[i]){
      eta[i,] = eta.new
      acc.vec[i] = T
    }else{
      eta[i,] = eta[i-1,]
      acc.vec[i] = F
    }
  }
  res$eta = eta
  res$acc.vec = acc.vec
  res$acc.prob = acc.prob
  return(res)
}

langevinMH <- function(init,dTar,dProp,rProp,grad,nsamples=10000){
  eta = matrix(data = NA,nrow = nsamples, ncol = 2)
  acc.vec = numeric(nsamples)
  acc.prob = numeric(nsamples)
  eta[1,] = init$mu
  starttime = Sys.time()
  res = list()
  pb <- txtProgressBar(min = 0, max = nsamples, style = 3)
  for (i in seq(2,nsamples)){
    setTxtProgressBar(pb, i)
    eta.new = rProp(y = eta[i-1,] + init$sigma^2/2*grad(eta[i-1,]),init$sigma)
    acc.prob[i] = min(1,exp(dTar(eta.new) + dProp(x = eta[i-1,], y = eta.new+init$sigma^2/2*grad(eta.new),sigma = init$sigma)  - dTar(eta[i-1,]) - dProp(x = eta.new, y = eta[i-1,]+init$sigma^2/2*grad(eta[i-1,]), sigma = init$sigma)))
    if (runif(1) < acc.prob[i]){
      eta[i,] = eta.new
      acc.vec[i] = T
    }else{
      eta[i,] = eta[i-1,]
      acc.vec[i] = F
    }
  }
  res$eta = eta
  res$acc.vec = acc.vec
  res$acc.prob = acc.prob
  return(res)
}


hamiltonianMH <- function(init,dTar,dProp,rProp,grad,epsilon,steps,nsamples=10000){
  eta = matrix(data = NA,nrow = nsamples, ncol = 2)
  acc.vec = numeric(nsamples)
  acc.prob = numeric(nsamples)
  eta[1,] = init$mu
  starttime = Sys.time()
  res = list()
  pb <- txtProgressBar(min = 0, max = nsamples, style = 3)
  for (i in seq(2,nsamples)){
    setTxtProgressBar(pb, i)
    eta.new = rProp(y = eta[i-1,],init$sigma)
    p.new = rnorm(n = length(eta.new),mean=0,sd = 1)
    p.init = p.new
    for (j in seq(steps)){
      p.new = p.new - epsilon/2*grad(eta.new)
      eta.new = eta.new + epsilon*p.new
      p.new = p.new - epsilon/2*grad(eta.new)
    }
    p.new = -p.new
    acc.prob[i] = min(1,exp(dTar(eta.new) - dTar(eta[i-1,]) + sum(p.init^2)/2 - sum(p.new^2)/2))
    if (is.nan(acc.prob[i])){
      acc.prob[i]=0
    }
    if (runif(1) < acc.prob[i]){
      eta[i,] = eta.new
      acc.vec[i] = T
    }else{
      eta[i,] = eta[i-1,]
      acc.vec[i] = F
    }
  }
  res$eta = eta
  res$acc.vec = acc.vec
  res$acc.prob = acc.prob
  return(res)
}

plottrace <- function(eta){
  t2 <- ggplot(data.frame(x=seq(nrow(eta)),y=eta[,1]))+
    geom_path(aes(x=x,y=y)) + 
    labs(x="Index",y="")
  t1 <- ggplot(data.frame(x = seq(nrow(eta)),y=eta[,2]))+
    geom_path(aes(x=x,y=y)) + 
    labs(x="Index",y="")
  gridExtra::grid.arrange(t1,t2)
}

plothist <- function(eta){
  ggplot(data.frame(x = eta[,1],y = eta[,2]),aes(x=x,y=y))+
    geom_point(alpha =0.4) + 
    geom_density_2d_filled(alpha = 0.9) 
}



### Metropolis-Hastings MH

# Standard Normal MH
init = list(mu = c(0,0),sigma = rep(0.7,2))
res_StNorm = MH(init, dStNorm, dProp, rProp, nsamples=10000)
plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)


# Multimodal Normal MH
init = list(mu = c(0,0),sigma = rep(1.5,2))
res_MMod = MH(init, dMMod, dProp, rProp, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)


# Volcano MH
init = list(mu =c(0,0),sigma = rep(1,2))
res_Volcano = MH(init, dVolcano, dProp, rProp, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)


### Langevin MH

# Standard Normal Langevin MH
init = list(mu = c(0,0),sigma = rep(0.3,2))
res_StNorm = langevinMH(init, dStNorm, dProp, rProp,gradient.dStNorm, nsamples=10000)
plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)


# Multimodal Normal Langevin MH
init = list(mu = c(0,0),sigma = rep(0.9,2))
res_MMod = langevinMH(init, dMMod, dProp, rProp,gradient.dMMod, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)


# Volcano Langevin MH
init = list(mu =c(0,0),sigma = rep(0.8,2))
res_Volcano = langevinMH(init, dVolcano, dProp, rProp,gradient.dVolcano, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)



### Hamiltonian MH

# Standard Normal Hamiltonian MH
init = list(mu = c(0,0),sigma = rep(0.3,2))
res_StNorm = hamiltonianMH(init, dStNorm, dProp, rProp,gradient.dStNorm,epsilon=0.1,steps=10,nsamples=10000)
plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)


# Multimodal Normal Hamiltonian MH
init = list(mu = c(0,0),sigma = rep(0.9,2))
res_MMod = hamiltonianMH(init, dMMod, dProp, rProp,gradient.dMMod,epsilon=0.2,steps=15, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)


# Volcano Hamiltonian MH
init = list(mu =c(0,0),sigma = rep(0.8,2))
res_Volcano = hamiltonianMH(init, dVolcano, dProp, rProp,gradient.dVolcano,epsilon=0.08,steps=17, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)
