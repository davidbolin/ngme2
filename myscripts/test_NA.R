# # Example Y = Y1 NA Y2 NA Y3
# library(devtools); load_all()

# set.seed(7)
# set.seed(Sys.time())

# n_obs <- 5
# ar_mu <- 0
# ar_sigma <- 1.3
# ar_eta <- 0.8

# ar1 <- ngme.simulate(
#   f(1:n_obs,
#     model = "ar1",
#     theta_K = 0.6,
#     noise = ngme.noise.nig(
#       theta_mu = ar_mu,
#       theta_sigma = ar_sigma,
#       theta_V = ar_eta
#     )
#   ),
#   seed = 123
# )

# # normal noise
# Y <- ar1$realization + rnorm(n_obs)

# # set NA for prediction
# pred_index <- c(2, 4)
# Y[pred_index] <- NA

# data_index <- which(!is.na(Y))
# est_A <- ngme.ts.make.A(data_index, rep(1, length(data_index)))
# pred_A <- ngme.ts.make.A(pred_index, rep(1, length(pred_index)))

# est_A
# pred_A

# ngme.ts.make.A(1:3, NULL)

# ngme_control <- ngme.control(
#   estimation = TRUE,
#   burnin = 2,
#   iterations = 2,
#   gibbs_sample = 5,
#   stepsize = 1,
#   kill_var = FALSE,
#   threshold = 1e-4,
#   opt_beta = T
# )
# ngme_out <- ngme(
#   Y ~ 0 +
#   f(1:n_obs,
#     replicates = NULL,
#     model = "ar1",
#     theta_K = 0.6,
#     W = ar1$realization,
#     noise = ngme.noise(
#       theta_mu = ar_mu,
#       theta_sigma = ar_sigma,
#       theta_V = ar_eta,
#       V = ar1$noise$V
#     ),
#     control = ngme.control.f(
#       numer_grad       = FALSE,
#       use_precond      = TRUE,

#       fix_operator     = FALSE,
#       fix_W            = FALSE
#     ),
#     debug = TRUE
#   ),
#   data = data.frame(),
#   control = ngme_control,
#   noise = ngme.noise.normal(),
#   debug = ngme.debug(
#     debug = TRUE,
#     not_run = F
#   ),
#   seed = 1
# )

# str(ngme_out$est_output)
# # plot_out(ngme_out$opt_trajectory, start = 5)
# plot_out(ngme_out$opt_trajectory, start = 2)


# ngme.noise.nig(
#   theta_mu = 0,
#   theta_sigma = 0,
#   theta_V = 1,
#   V = NULL,
#   B_mu = matrix(1),
#   B_sigma = matrix(1),
#   fix_mu = TRUE
# )

library(devtools); load_all()
Y1 <- 1:5
ngme.ts.make.A(c(1,3,5))

Y <- c(0.5, NA, 0.2, NA, 0.3)

Y_data <- na.omit(Y); Y_data
index_data <- which(!is.na(Y))        # 1 3 5
index_pred <- which(is.na(Y))        # 2 4

# estimation
# out <- ngme(
#   # Y_data ~ f(c(1,3,5), model = "ar1", noise = ngme.noise.nig()),
#   Y1 ~ f(1:5, model = "ar1", noise = ngme.noise.nig()),
#   noise = ngme.noise.normal(),
#   data = data.frame(Y1 = Y1),
#   control = ngme.control(burnin=5, iteration=1)
# )

out <- ngme(
  Y_data ~ f(index_data, model = "ar1", noise = ngme.noise.nig()),
  noise = ngme.noise.normal(),
  data = data.frame(Y_data = Y_data),
  control = ngme.control(burnin=5, iteration=1)
)

# prediction
# A_pred <- ngme.ts.make.A(index_pred, start=1, end=5) # whole mesh is 1:5
# prediction <- A_pred * out$W

load_all()
ngme.ar1(index = c(1, 3, 5), replicates = NULL)
