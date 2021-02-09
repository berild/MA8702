
library(ggplot2)

dStNorm <- function(x,log=TRUE){
  mvtnorm::dmvnorm(x,mean = rep(0,2),sigma = matrix(c(1,0.9,0.9,1),nrow=2,ncol=2),log = log)
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

dVolcano <- function(x,log =TRUE){
  x = c(x[1],x[2])
  if (log){
    log(1/(2*pi)*exp(-1/2*t(x)%*%x)*(t(x)%*%x + 0.25))
  }else
    1/(2*pi)*exp(-1/2*t(x)%*%x)*(t(x)%*%x + 0.25)
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


prior <- function(x, y = rep(0,2), sigma = diag(1,2),log = TRUE) {
  mvtnorm::dmvnorm(x,mean = y,sigma = sigma,log=log)
}

dProp <- function(x,y,sigma = diag(0.5,2),log=TRUE){
  mvtnorm::dmvnorm(x,mean = y,sigma = sigma,log=log)
}

rProp <- function(y,sigma = diag(0.5,2),log=TRUE){
  mvtnorm::rmvnorm(1,mean = y,sigma = sigma,log=log)
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
    eta.new = rProp(y = eta[i-1,], sigma = init$cov)
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

LangevinMH <- function(init,dTar,dProp,rProp,nsamples=10000){
  
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
init = list(mu = c(0,0),cov = diag(0.6,2))
res_StNorm = MH(init, dStNorm, dProp, rProp, prior, nsamples=10000)
plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)


# Multimodal Normal MH
init = list(mu = c(0,0),cov = diag(3,2))
res_MMod = MH(init, dMMod, dProp, rProp, prior, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)


# Volcano MH
init = list(mu =c(0,0),cov = diag(1.5,2))
res_Volcano = MH(init, dVolcano, dProp, rProp, prior, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)


### Langevin MH

