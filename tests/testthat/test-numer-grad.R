test_that("numerical gradient works", {
  # 1. test OU process
  n_obs <- 1000
  kernel_num <- 5

  mean_BK<-round(seq(1, n_obs, length=kernel_num))
  dens_k <- list()
  mesh = seq(1, 4*kernel_num, length = n_obs)
  dens_k[[1]]<-rep(1, n_obs)
  for (i in 2:(kernel_num+1)) {
    dens_k[[i]] <- dnorm(x=mesh,mean=mesh[mean_BK[[i-1]]], sd=1)
  }
  B_K <-as.matrix(as.data.frame(dens_k))

  theta_k <- c(-1, 1, 1.5, -1.5, 2, -1)
  ou <- f(1:n_obs, model="ou",
    B_K = B_K,
    theta_K = theta_k,
    noise = noise_nig()
  )
  W <- simulate(ou)[[1]]
  Y <- W + rnorm(n_obs, sd=0.01)

  out <- ngme(
    Y ~ 0 + f(1:n_obs, model="ou", name="ou",
      # fix_W = T, W = W,
      # fix_theta_K = T,
      B_K = B_K,
      theta_K = rep(0, 6),
      noise=noise_nig(
        # fix_theta_mu = T,
        # fix_V = T, V = attr(eps, "noise")$V
      )
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      burnin = 100,
      iterations = 100,
      n_parallel_chain = 4,
      print_check_info = FALSE,
      verbose = F
    ),
    control_ngme =control_ngme(n_gibbs_samples = 5)
  )
  out
  traceplot(out, "ou")
  out$replicates[[1]]$latent[[1]]$theta_K

  plot(exp(B_K %*% out$replicates[[1]]$models$ou$theta_K),main="estimation of K")
  lines(exp(B_K %*% theta_k),main="true K", col="red")

  expect_true(TRUE)
})
