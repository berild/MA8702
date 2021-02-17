


# Exercise 1.1:


plotdensity(dStNorm)

plotdensity(dMMod)

plotdensity(dVolcano)

# Exercise 1.2:



### Metropolis-Hastings MH

# Standard Normal MH
init = list(mu = c(0,0),sigma = rep(5,2))
res_StNorm = MH(init, dStNorm, dProp, rProp, nsamples=10000)
plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)

mean(res_StNorm$acc.prob)
acf(res_StNorm$eta[,2],main="Autocorrelation for x_2")


res_StNorm_sigma_005 = MH(list(mu = c(0,0),sigma = rep(0.05,2)), dStNorm, dProp, rProp, nsamples=10000)
res_StNorm_sigma_05 = MH(list(mu = c(0,0),sigma = rep(0.5,2)), dStNorm, dProp, rProp, nsamples=10000)
res_StNorm_sigma_5 = MH(list(mu = c(0,0),sigma = rep(5,2)), dStNorm, dProp, rProp, nsamples=10000)

res = data.frame("sigma: 0.05" = res_StNorm_sigma_005$eta[,1], "sigma: 0.05" = res_StNorm_sigma_005$eta[,2],  "sigma: 0.5" = res_StNorm_sigma_05$eta[,1],
                 "sigma: 0.5" = res_StNorm_sigma_05$eta[,2], "sigma: 5" = res_StNorm_sigma_5$eta[,1],  "sigma: 5" = res_StNorm_sigma_5$eta[,2])

res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()


bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()





# Multimodal Normal MH
init = list(mu = c(0,0),sigma = rep(5,2))
res_MMod = MH(init, dMMod, dProp, rProp, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)

mean(res_MMod$acc.prob)
acf(res_MMod$eta[,2],main="Autocorrelation for x_2")


res_MMod_sigma_005 = MH(list(mu = c(0,0),sigma = rep(0.05,2)), dMMod, dProp, rProp, nsamples=10000)
res_MMod_sigma_05 = MH(list(mu = c(0,0),sigma = rep(0.5,2)), dMMod, dProp, rProp, nsamples=10000)
res_MMod_sigma_5 = MH(list(mu = c(0,0),sigma = rep(5,2)), dMMod, dProp, rProp, nsamples=10000)

res = data.frame("sigma: 0.05" = res_MMod_sigma_005$eta[,1], "sigma: 0.05" = res_MMod_sigma_005$eta[,2],  "sigma: 0.5" = res_MMod_sigma_05$eta[,1],
                 "sigma: 0.5" = res_MMod_sigma_05$eta[,2], "sigma: 5" = res_MMod_sigma_5$eta[,1],  "sigma: 5" = res_MMod_sigma_5$eta[,2])

res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()



bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()





# Volcano MH
init = list(mu =c(0,0),sigma = rep(5,2))
res_Volcano = MH(init, dVolcano, dProp, rProp, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)

mean(res_Volcano$acc.prob)
acf(res_Volcano$eta[,2],main="Autocorrelation for x_2")


res_Volcano_sigma_005 = MH(list(mu = c(0,0),sigma = rep(0.05,2)), dVolcano, dProp, rProp, nsamples=10000)
res_Volcano_sigma_05 = MH(list(mu = c(0,0),sigma = rep(0.5,2)), dVolcano, dProp, rProp, nsamples=10000)
res_Volcano_sigma_5 = MH(list(mu = c(0,0),sigma = rep(5,2)), dVolcano, dProp, rProp, nsamples=10000)

res = data.frame("sigma: 0.05" = res_Volcano_sigma_005$eta[,1], "sigma: 0.05" = res_Volcano_sigma_005$eta[,2],  "sigma: 0.5" = res_Volcano_sigma_05$eta[,1],
                 "sigma: 0.5" = res_Volcano_sigma_05$eta[,2], "sigma: 5" = res_Volcano_sigma_5$eta[,1],  "sigma: 5" = res_Volcano_sigma_5$eta[,2])


res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()



bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()



# Exercise 1.3

### Langevin MH

# Standard Normal Langevin MH
init = list(mu = c(0,0),sigma = rep(0.5,2))
res_StNorm = langevinMH(init, dStNorm, dProp, rProp,gradient.dStNorm, nsamples=10000)
plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)

mean(res_StNorm$acc.prob)
acf(res_StNorm$eta[,1],main="Autocorrelation for x_1")
acf(res_StNorm$eta[,2],main="Autocorrelation for x_2")

res_StNorm_sigma_02 = langevinMH(list(mu = c(0,0),sigma = rep(0.2,2)), dStNorm, dProp, rProp,gradient.dStNorm, nsamples=10000)
res_StNorm_sigma_05 = langevinMH(list(mu = c(0,0),sigma = rep(0.5,2)), dStNorm, dProp, rProp,gradient.dStNorm, nsamples=10000)
res_StNorm_sigma_1 = langevinMH(list(mu = c(0,0),sigma = rep(1,2)), dStNorm, dProp, rProp,gradient.dStNorm, nsamples=10000)

res = data.frame("sigma: 0.2" = res_StNorm_sigma_02$eta[,1], "sigma: 0.2" = res_StNorm_sigma_02$eta[,2], "sigma: 0.5" = res_StNorm_sigma_05$eta[,1],  
                 "sigma: 0.5" = res_StNorm_sigma_05$eta[,2], "sigma: 1" = res_StNorm_sigma_1$eta[,1], "sigma: 1" = res_StNorm_sigma_1$eta[,2])



res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()




bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()





# Multimodal Normal Langevin MH
# 0.5, 1, 0.2
init = list(mu = c(0,0),sigma = rep(0.5,2))
res_MMod = langevinMH(init, dMMod, dProp, rProp,gradient.dMMod, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)

mean(res_MMod$acc.prob)
acf(res_MMod$eta[,1],main="Autocorrelation for x_1")
acf(res_MMod$eta[,2],main="Autocorrelation for x_2")


res_MMod_sigma_02 = langevinMH(list(mu = c(0,0),sigma = rep(0.2,2)), dMMod, dProp, rProp,gradient.dMMod, nsamples=10000)
res_MMod_sigma_05 = langevinMH(list(mu = c(0,0),sigma = rep(0.5,2)), dMMod, dProp, rProp,gradient.dMMod, nsamples=10000)
res_MMod_sigma_1 = langevinMH(list(mu = c(0,0),sigma = rep(1,2)), dMMod, dProp, rProp,gradient.dMMod, nsamples=10000)

res = data.frame("sigma: 0.2" = res_MMod_sigma_02$eta[,1], "sigma: 0.2" = res_MMod_sigma_02$eta[,2], "sigma: 0.5" = res_MMod_sigma_05$eta[,1],  
                  "sigma: 0.5" = res_MMod_sigma_05$eta[,2], "sigma: 1" = res_MMod_sigma_1$eta[,1], "sigma: 1" = res_MMod_sigma_1$eta[,2])



res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()



bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)

bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()




# Volcano Langevin MH
init = list(mu =c(0,0),sigma = rep(0.5,2))
res_Volcano = langevinMH(init, dVolcano, dProp, rProp,gradient.dVolcano, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)

mean(res_Volcano$acc.prob)
acf(res_Volcano$eta[,1],main="Autocorrelation for x_1")
acf(res_Volcano$eta[,2],main="Autocorrelation for x_2")


res_Volcano_sigma_02 = langevinMH(list(mu = c(0,0),sigma = rep(0.2,2)), dVolcano, dProp, rProp,gradient.dVolcano, nsamples=10000)
res_Volcano_sigma_05 = langevinMH(list(mu = c(0,0),sigma = rep(0.5,2)), dVolcano, dProp, rProp,gradient.dVolcano, nsamples=10000)
res_Volcano_sigma_1 = langevinMH(list(mu = c(0,0),sigma = rep(1,2)), dVolcano, dProp, rProp,gradient.dVolcano, nsamples=10000)

res = data.frame("sigma: 0.2" = res_Volcano_sigma_02$eta[,1], "sigma: 0.2" = res_Volcano_sigma_02$eta[,2], "sigma: 0.5" = res_Volcano_sigma_05$eta[,1], 
                  "sigma: 0.5" = res_Volcano_sigma_05$eta[,2], "sigma: 1" = res_Volcano_sigma_1$eta[,1],  "sigma: 1" = res_Volcano_sigma_1$eta[,2])



res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()




bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)

bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()



# Exercise 1.4
### Hamiltonian MH

# alt 1: epsilon = 0.05, n_steps = 10,
# alt 2: epsilon = 0.1, n_steps = 10
# alt 3: epsilon = 0.1, n_steps = 15

# Standard Normal Hamiltonian MH
init = list(mu = c(0,0),sigma = rep(0.5,2))

res_StNorm = hamiltonianMH(init, dStNorm, dProp, rProp,gradient.dStNorm,epsilon=epsilon_1,steps=n_steps,nsamples=10000)

plottrace(res_StNorm$eta)
plothist(res_StNorm$eta)

mean(res_StNorm$acc.prob)
acf(res_StNorm$eta[,1],main="Autocorrelation for x_1")
acf(res_StNorm$eta[,2],main="Autocorrelation for x_2")



epsilon_1 = 0.1
epsilon_2 = 0.5
epsilon_3 = 0.64
n_steps_1 = 10


res_StNorm_11  = hamiltonianMH(init, dStNorm, dProp, rProp,gradient.dStNorm,epsilon=epsilon_1,steps=n_steps_1,nsamples=10000)
res_StNorm_21 = hamiltonianMH(init, dStNorm, dProp, rProp,gradient.dStNorm,epsilon=epsilon_2,steps=n_steps_1,nsamples=10000)
res_StNorm_31 = hamiltonianMH(init, dStNorm, dProp, rProp,gradient.dStNorm,epsilon=epsilon_3,steps=n_steps_1,nsamples=10000)

plottrace(res_StNorm_11$eta)
mean(res_StNorm_11$acc.prob)
plothist(res_StNorm_11$eta)

plottrace(res_StNorm_21$eta)
mean(res_StNorm_21$acc.prob)
plothist(res_StNorm_21$eta)

plottrace(res_StNorm_31$eta)
mean(res_StNorm_31$acc.prob)
plothist(res_StNorm_31$eta)


res = data.frame("epsilon: 0.1" = res_StNorm_11$eta[,1], "epsilon: 0.1" = res_StNorm_11$eta[,2], "epsilon: 0.5" = res_StNorm_21$eta[,1], 
                  "epsilon: 0.5" = res_StNorm_21$eta[,2], "epsilon: 0.64" = res_StNorm_31$eta[,1], "epsilon: 0.64" = res_StNorm_31$eta[,2])




res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()



bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()




# Multimodal Normal Hamiltonian MH
init = list(mu = c(0,0),sigma = rep(0.5,2))
epsilon_1 = 0.1
step_size_1 = 15

res_MMod = hamiltonianMH(init, dMMod, dProp, rProp,gradient.dMMod,epsilon=epsilon_1,steps=step_size_1, nsamples=10000)
plottrace(res_MMod$eta)
plothist(res_MMod$eta)

mean(res_MMod$acc.prob)
acf(res_MMod$eta[,1],main="Autocorrelation for x_1")
acf(res_MMod$eta[,2],main="Autocorrelation for x_2")




epsilon_1 = 0.05
epsilon_2 = 0.1
n_steps_1 = 5
n_steps_2 = 15

res_MMod_11  = hamiltonianMH(init, dMMod, dProp, rProp,gradient.dMMod,epsilon=epsilon_1,steps=n_steps_1,nsamples=10000)
res_MMod_21  = hamiltonianMH(init, dMMod, dProp, rProp,gradient.dMMod,epsilon=epsilon_2,steps=n_steps_1,nsamples=10000)
res_MMod_22  = hamiltonianMH(init, dMMod, dProp, rProp,gradient.dMMod,epsilon=epsilon_2,steps=n_steps_2,nsamples=10000)



plottrace(res_MMod_11$eta)
mean(res_MMod_11$acc.prob)
plothist(res_MMod_11$eta)

plottrace(res_MMod_21$eta)
mean(res_MMod_21$acc.prob)
plothist(res_MMod_21$eta)

plottrace(res_MMod_22$eta)
mean(res_MMod_22$acc.prob)
plothist(res_MMod_22$eta)


res = data.frame("Combination 1" = res_MMod_11$eta[,1], "Combination 1" = res_MMod_11$eta[,2],"Combination 2" = res_MMod_21$eta[,1], 
                  "Combination 2" = res_MMod_21$eta[,2],"Combination 2" = res_MMod_22$eta[,1], "Combination 2" = res_MMod_22$eta[,2])




res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()


bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()






# Volcano Hamiltonian MH
init = list(mu =c(0,0),sigma = rep(0.5,2))
epsilon = 0.1
step_size_1 = 15
res_Volcano = hamiltonianMH(init, dVolcano, dProp, rProp,gradient.dVolcano,epsilon=epsilon_1,steps=step_size_1, nsamples=10000)
plottrace(res_Volcano$eta)
plothist(res_Volcano$eta)

mean(res_Volcano$acc.prob)
acf(res_Volcano$eta[,1],main="Autocorrelation for x_1")
acf(res_Volcano$eta[,2],main="Autocorrelation for x_2")




epsilon_1 = 0.01
epsilon_2 = 0.1
n_steps_1 = 10
n_steps_2 = 30

res_Volcano_11  = hamiltonianMH(init, dVolcano, dProp, rProp,gradient.dVolcano,epsilon=epsilon_1,steps=n_steps_1,nsamples=10000)
res_Volcano_21  = hamiltonianMH(init, dVolcano, dProp, rProp,gradient.dVolcano,epsilon=epsilon_2,steps=n_steps_1,nsamples=10000)
res_Volcano_22  = hamiltonianMH(init, dVolcano, dProp, rProp,gradient.dVolcano,epsilon=epsilon_2,steps=n_steps_2,nsamples=10000)




plottrace(res_Volcano_11$eta)
mean(res_Volcano_11$acc.prob)
plothist(res_Volcano_11$eta)

plottrace(res_Volcano_21$eta)
mean(res_Volcano_21$acc.prob)
plothist(res_Volcano_21$eta)

plottrace(res_Volcano_22$eta)
mean(res_Volcano_22$acc.prob)
plothist(res_Volcano_22$eta)


res = data.frame("Combination 1" = res_Volcano_11$eta[,1],"Combination 1" = res_Volcano_11$eta[,2], "Combination 2" = res_Volcano_21$eta[,1], 
                  "Combination 2" = res_Volcano_21$eta[,2],  "Combination 3" = res_Volcano_22$eta[,1], "Combination 3" = res_Volcano_22$eta[,2])



res_2 = data.frame(steps=seq(nrow(res)), res)
res_2 = reshape2::melt(res_2, id.vars='steps',variable.name='vars')

ggplot(res_2, aes(x=steps, y=value))+geom_line()+facet_wrap(vars~ .,scales = 'free', ncol=2)+ theme_bw()



bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 40,plot = FALSE)$acf))
bacdf$lag <- 0:(nrow(bacdf) - 1)
significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))
bacdf <- bacdf %>% select(lag, everything())
bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')
ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point()  + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()



