library(rstan) # observe startup messages
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

library(coda)

effectiveSize(res)
