library(rstan)

# load data
source("mungeData.R")

# compile model
multi.modFL <- stan_model(file = "multilevel_multistate.stan")

# sampling call to reproduce results in manuscript
multi.fitFL.test <- sampling(multi.modFL, data = mod_data, iter = 2000, seed = 1234,
                                 init = 0,
                       control = list(adapt_delta = 0.999,
                                      stepsize = 0.01))

# example results printout: estimated cumulative survival for Yolo 2012 fish:
print(multi.fitFL.test, pars = "pred_survYolo12")
