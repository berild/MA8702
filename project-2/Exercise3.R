# load required libraries
library(INLA)
library(ggplot2)

# load data
data = read.delim(file = "skijump.txt",header=T)

plm <- ggplot() +
  geom_point(data = as.data.frame(data),aes(x=Year,y=Length)) +
  geom_line(data = data.frame(Year = data$Year, Length = predict(lm(Length~Year, data= data),data=data$Year)),aes(x=Year,y=Length),color = "firebrick") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12), axis.title = element_text(size=14))
plm

# fit simple model
res = inla(Length ~ Year, data = data)

# inspecting summary of INLA object
summary(res)

# plotting marginal density of the fixed effect (Year)
# Intercept
p1 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.fixed$`(Intercept)`),aes(x=x,y=y)) +
  labs(y="",x = expression(beta[0]),title = "a") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p1
# Fixed effect of Year
p2 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.fixed$Year),aes(x=x,y=y)) +
  labs(y="",x = expression(beta["Year"]),title = "b") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p2

# plotting marginal density of the hyperparameter tau (precision)
p3 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.hyperpar$`Precision for the Gaussian observations`),aes(x=x,y=y)) +
  labs(y="",x = expression(tau),title = "c") +
  theme_bw() +
  coord_cartesian(xlim=c(0,1)) +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p3

# computing the transformation of the precision to the variance
sigma_df <- as.data.frame(inla.tmarginal(function(x){1/sqrt(x)},res$marginals.hyperpar$`Precision for the Gaussian observations`))

# plotting marginal density of the hyperparameter sigma (variance)
p4 <- ggplot() +
  geom_line(data = sigma_df,aes(x=x,y=y)) +
  labs(y="",x = expression(sigma),title = "d") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p4

# computing posterior estimates of sigma
inla.zmarginal(sigma_df)

# loading Seed dataset
data(Seeds)

# specifying the formula
formula = r ~ x1 + x2 + f(plate, model = "iid")

# fitting binomial model
# default link: logit
res = inla(formula,data=Seeds,family="binomial",Ntrials=n)

# inspecting summary of INLA object
summary(res)

# plotting marginal densities of model parameters
# intercept
p1 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.fixed$`(Intercept)`),aes(x=x,y=y)) +
  labs(y="",x = expression(a[0]),title = "a") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.01,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p1

# plotting marginal density fixed effect of x1
p2 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.fixed$x1),aes(x=x,y=y)) +
  labs(y="",x = expression(a[1]),title = "b") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.01,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p2

# plotting marginal density fixed effect of x2
p3 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.fixed$x2),aes(x=x,y=y)) +
  labs(y="",x = expression(a[2]),title = "c") +
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.01,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p3

# plotting marginal density precision of iid random effect
p4 <- ggplot() +
  geom_line(data = as.data.frame(res$marginals.hyperpar$`Precision for plate`),aes(x=x,y=y)) +
  labs(y="",x = expression(tau),title = "d") +
  theme_bw() +
  coord_cartesian(xlim = c(0,0.4e+05))+
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.01,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p4

# Checking if variance looks better than marginal of precision
sigma_df <- as.data.frame(inla.tmarginal(function(x){1/sqrt(x)},res$marginals.hyperpar$`Precision for plate`))
# plotting marginal density of the hyperparameter sigma (variance)
p5 <- ggplot() +
  geom_line(data = sigma_df,aes(x=x,y=y)) +
  labs(y="",x = expression(sigma),title = "d") +
  theme_bw() +
  coord_cartesian(xlim = c(0,0.2)) +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.01,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p5
