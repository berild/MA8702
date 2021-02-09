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
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  refresh = 1000,         # show progress every 1000 iterations
  seed = 123
)

res <- as.data.frame(fit)
res = res[,-ncol(res)]

res = cbind(steps = seq(nrow(res)),res)


effectiveSize(res)
res
df = reshape2::melt(res ,  id.vars = 'steps', variable.name = 'vars')
ggplot(df, aes(steps,value)) + geom_line() + facet_grid(vars ~ .,scales = "free")
ggplot(df, aes(x = value, y=..density..)) + geom_histogram() + facet_wrap(vars ~ .,scales = "free")

bacdf <- map_df(res, function(ts) as.vector(acf(ts, plot = FALSE)$acf))

# The lags are all the same just 0 through the number of rows minus 1

bacdf$lag <- 0:(nrow(bacdf) - 1)

# reorder things and eliminate `Date` and unclean AAPL which is actually identical to AAPL

bacdf <- bacdf %>% select(lag, everything(), -steps)

bacdf = reshape2::melt(bacdf ,  id.vars = 'lag', variable.name = 'vars')

ggplot(bacdf, aes(lag,value)) + geom_segment(mapping = aes(xend = lag, yend = 0)) + geom_point() +facet_grid(vars ~ .,scales = "free") + theme_minimal()

