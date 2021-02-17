library(rstan) # observe startup messages
library(dplyr)
library(purrr)
library(coda)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
pump_data <- list(N = 10,
                    y = c(5,  1, 5,  14, 3,  19, 1, 1, 4, 22),
                    t = c(94.3, 15.7, 62.9, 126.0,  5.24, 31.4, 1.05, 1.05, 2.1, 10.5))


fit <- stan(
  file = "pump.stan",  # Stan program
  data = pump_data,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 3000,            # total number of iterations per chain
  refresh = 1000,         # show progress every 1000 iterations
  seed = 123
)

res <- as.data.frame(fit)
res = res[,-ncol(res)]

stan_trace(fit,pars = colnames(res),ncol = 3,inc_warmup = TRUE)
stan_ess(fit,pars = colnames(res))
stan_mcse(fit,pars = colnames(res))
stan_plot(fit,point_est = "median",pars = colnames(res),show_density = TRUE,fill_color="lightskyblue2",outline_color="dodgerblue3",est_color = "firebrick")

effectiveSize(res)
sapply(seq(ncol(res)), function(x){sd(res[,x])})
sapply(seq(ncol(res)), function(x){mean(res[,x])})
sapply(seq(ncol(res)), function(x){median(res[,x])})

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

sapply(seq(ncol(res)), function(x){getmode(res[,x])})

df = reshape2::melt(res ,  id.vars = 'steps', variable.name = 'vars')
ggplot(df, aes(steps,value)) + geom_line() + facet_grid(vars ~ .,scales = "free")
ggplot(df, aes(x = value, y=..density..)) + geom_histogram(color="dodgerblue4",fill="lightskyblue",bins = 30) + geom_density(color = "firebrick",size = 1) + facet_wrap(vars ~ .,scales = "free",ncol=3) +theme_bw()

bacdf <- map_df(res, function(ts) as.vector(acf(ts, lag.max = 10,plot = FALSE)$acf))

bacdf$lag <- 0:(nrow(bacdf) - 1)

significance_level <- qnorm((1 + 0.95)/2)/sqrt(nrow(res))

bacdf <- bacdf %>% select(lag, everything())

bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')

ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point() + geom_hline(yintercept = c(significance_level, -significance_level), lty = 3, color = "blue") + labs(x="Lag",y = "ACF") + facet_wrap(vars ~ .,scales = "fixed",ncol = 3) + theme_bw()

