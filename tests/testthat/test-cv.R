# test_that("test merged replicate", {
#   # creating ar1 process of 2 replicates
#   n <- 5
#   df <- data.frame(
#     Y = 1:10,
#     idx = 1:10,
#     group = rep(1:2, each = n)
#   )

#   ar_idx <- c(1,2,3,4,5, 3,2,1,5,4)
#   rw_idx <- c(5,4,3,2,1, 1,4,3,5,2)
#   mod.replicates <- ngme(
#     formula = Y ~ 1 + f(
#       ar_idx,
#       model = "ar1",
#       noise = noise_nig()
#     ) + f(
#       rw_idx,
#       model = "rw1",
#       noise = noise_nig()
#     ),
#     replicate = df$group,
#     data = df,
#     control_opt = control_opt(
#       print_check_info = FALSE
#     )
#   )

#   ret <- merge_replicates(
#     mod.replicates, 
#     train_idx = c(1:4, 6:9),
#     test_idx = c(5, 10)
#   )

#   expect_equal(ret$test_Y, c(5,10))

#   m <- ret$merged_rep
#   expect_equal(m$Y, c(1,2,3,4,6,7,8,9))
#   print(m$models[[1]]$A)
#   print(m$models[[2]]$A)

#   print(ret$test_A_block)
#   # [1,] . . . 0 1 1 0 . . .
#   # [2,] . . . 1 0 . 1 0 . .
#   # Correct


#   ########## Compute scores ##########
#   s = compute_scores(
#     m, 10, 10, 1,
#     A_pred_block = ret$test_A_block,
#     test_noise = ret$test_noise,
#     y_data = ret$test_Y,
#     group_data = ret$test_group,
#     X_pred = ret$test_X,
#     transform = identity,
#     thining_gap = 0
#   )
#   s
# })

# # This test ensures the posterior sampling procedure in both C++ and R are the same!
# test_that("test cross validation (NIG and gaussian)", {
#   # load fitted models
#   ret_gauss <- readRDS("examples/models/ret_gauss.rds")
#   ret_nig <- readRDS("examples/models/ret_nig.rds")

#   seed <- 500
#   cv <- cross_validation(
#     list(gauss=ret_gauss, nig=ret_nig),
#     k=20,
#     parallel=TRUE,
#     n_gibbs_samples=100,
#     print=TRUE,
#     N_sim=4,
#     seed=seed+50
#   )
#   cv
#   gauss_mse <- cv$mean.scores$MSE[1]; gauss_mae <- cv$mean.scores$MAE[1]
#   nig_mse <- cv$mean.scores$MSE[2]; nig_mae <- cv$mean.scores$MAE[2]
#   expect_true(gauss_mse > nig_mse)
#   expect_true(gauss_mae > nig_mae)
# })


# test_that("test R and C version of prediction are the same", {
#   # Set Gaussian to FALSE to test NIG case
#   # gaussian <- FALSE
#   library(Matrix)
#   for (gaussian in c(TRUE, FALSE)){
#     seed <- 500
#     set.seed(seed)
    
#     n_obs <- 100
#     rho <- 0.5
#     nu <- 0.2
#     mu <- 2
#     h <- rep(1, n_obs)
#     delta <- -mu * h
#     sigma <- 1
#     sigma_eps <- 1
#     K <- ngme2::ar1(1:n_obs, rho = rho)$K
#     A_obs <- Matrix::Diagonal(n_obs)

#     # Generate true data
#     trueV <- ngme2::rig(n_obs, nu, nu * h^2, seed = seed)
#     mynoise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)

#     trueW <- solve(K, mynoise)
#     Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)

#     # Cholesky-based sampler for N(Q^{-1} m, Q^{-1})
#     sample_mvn_cholesky <- function(Q, m) {
#       R <- chol(Q)
#       mu <- solve(Q, m)
#       n <- length(m)
#       z <- rnorm(n)
#       sample <- mu + solve(R, z)
#       return(list(sample = sample, mu = mu))
#     }

#     burnin <- 100
#     iter <- 100

#     i <- 5
#     A <- A_obs
#     A <- A[-i,]
#     Y_new <- Y[-i]
#     post_W <- matrix(0, nrow=n_obs, ncol=iter)
#     post_m <- matrix(0, nrow=n_obs, ncol=iter)
    

#     if (gaussian){
#         V <- rep(1, n_obs)
#     } else{
#         V <- ngme2::rig(n_obs, nu, nu * h^2)
#     }
    
#     for (j in 1:iter) {
#       Q = sigma^(-2) * t(K) %*% diag(1/V) %*% K + sigma_eps^(-2) * t(A) %*% A
#       m = sigma^(-2) * t(K) %*% diag(1/V) %*% (mu*(V-h)) + sigma_eps^(-2) * t(A) %*% Y_new
#       m = as.numeric(m)
      
#       # sample N(Q^-1 m, Q^-1)
#       aux_W <- sample_mvn_cholesky(Q, m)
#       post_W[,j] <- aux_W$sample
#       post_m[,j] <- aux_W$mu
      
#       if (gaussian){
#         V <- rep(1, n_obs)
#       } else{
#         p <- rep(-1, n_obs)
#         a <- rep(nu, n_obs)
#         a <- a + (mu / sigma)^2
#         M <- K %*% post_W[,j] + mu * h
#         M <- as.vector(M)
#         b <- nu * h^2
#         b <- b + (M / sigma)^2

#         V <- ngme2::rgig(n_obs, p, a, b)         
#       }
#     }
#     last_W <- post_W[,iter]
#     last_m <- post_m[,iter]
#     last_W
#     last_m

#     idx <- 1:n_obs
#     idx <- idx[-i]
#     Y_new <- Y[-i]
#     ngme_model <- ngme2::ngme(
#       formula = Y ~ 0 + f(
#         idx,
#         model = "ar1",
#         mesh = fm_mesh_1d(1:n_obs),
#         rho = rho,
#         noise = if (gaussian) 
#             noise_normal(sigma=sigma) 
#           else noise_nig(
#             mu = mu,
#             sigma = sigma,
#             nu = nu
#           )
#       ),
#       data = data.frame(Y = Y_new),
#       family = noise_normal(sigma=sigma_eps),
#       control_opt = control_opt(
#         estimation = FALSE
#       )
#     )

#     ngme_repl <- ngme_model$replicates[[1]]
#     # Check that the A matrix is the same as the one used in the C++ code
#     expect_true(all(ngme_model$replicates[[1]]$models[[1]]$A == A_obs[-i, ]))
    
#     time_start <- Sys.time()
#     ret = sampling_cpp(
#       ngme_repl,
#       n = 100,
#       n_burnin = 100,
#       seed = seed,
#       posterior = TRUE
#     )
#     time_end <- Sys.time()
#     print(paste0("time taken: ", time_end - time_start))

#     if (gaussian){
#       last_m_cpp <- ret$cond_W[[100]]
#       plot(last_m_cpp, ylim = c(min(last_m, last_m_cpp), max(last_m, last_m_cpp)))
#       lines(last_m, col = "red")

#       # Looks good
#       print(paste0("the difference in gaussian case is ", mean(abs(last_m - last_m_cpp))))

#       first_m_cpp <- ret$cond_W[[1]]
#       first_m <- post_m[,1] 
#       plot(first_m_cpp, type = "l")
#       lines(first_m, col = "red")

#       # fixed
#       expect_true(mean(abs(last_m - last_m_cpp)) < 1e-10)
#       expect_true(mean(abs(first_m - first_m_cpp)) < 1e-10)
#       expect_true(all(first_m_cpp == last_m_cpp))
#       expect_true(all(first_m == last_m))
#     } else {
#       first_m <- post_m[,1]
#       first_m_cpp <- ret$cond_W[[1]]
#       plot(first_m_cpp, type="l", ylim = c(min(first_m, first_m_cpp), max(first_m, first_m_cpp)))
#       lines(first_m, col = "red")
#       # Looks poor
#       print(paste0("the difference in NIG case is (first m) ", mean(abs(first_m - first_m_cpp))))

#       last_m <- post_m[,iter]
#       last_m_cpp <- ret$cond_W[[100]]
#       plot(last_m_cpp, type="l", ylim = c(min(last_m, last_m_cpp), max(last_m, last_m_cpp)))
#       lines(last_m, col = "red")
#       # Looks close
#       print(paste0("the difference in NIG case is (last m) ", mean(abs(last_m - last_m_cpp))))

#       mean_m <- apply(post_m, 1, mean)
#       mean_m_cpp <- mean_list(ret$cond_W)
#       plot(mean_m, type = "l")
#       lines(mean_m_cpp, col = "red")
#       # Looks good
#       print(paste0("the difference in NIG case is (mean m) ", mean(abs(mean_m - mean_m_cpp))))
#       expect_true(mean(abs(mean_m - mean_m_cpp)) < 0.2)
#     }
#   }
# })

# # This test ensures CRPS and sCRPS are computed correctly
# test_that("test iid CRPS, sCRPS", {
#   set.seed(345)
#   n_obs <- 100
#   ys <- rnorm(n_obs, 0, 1)

#   n_sim <- 1000
#   y_sim1 <- rnorm(n_sim, 0, 1)
#   y_sim2 <- rnorm(n_sim, 0, 1)
#   CRPS <- numeric(n_obs)
#   sCRPS <- numeric(n_obs)
#   for (i in 1:n_obs) {
#     y <- ys[i]
#     E_sim_data <- mean(abs(y_sim1 - y))
#     E_sim_sim <- mean(abs(y_sim1 - y_sim2))
#     CRPS[i] <- 0.5 * E_sim_sim - E_sim_data
#     sCRPS[i] <- -E_sim_data / E_sim_sim - 0.5 * log(E_sim_sim)
#   }
#   mean(CRPS)
#   mean(sCRPS)

#   ng <- ngme(
#     y ~ 0,
#     data = data.frame(y = ys),
#     family = noise_normal(
#       sigma = 1
#     ),
#     control_opt = control_opt(estimation = FALSE)
#   )
#   ng

#   # Cross-validation of the ngme model
#   cv_iid <- ngme2::cross_validation(
#     list(
#       ng = ng
#     ),
#     type = "loo",
#     n_gibbs_samples = 1000,
#     print = TRUE,
#     thining_gap = 0
#   )
#   cv_iid

#   expect_true(abs(cv_iid$mean.scores$neg.CRPS + mean(CRPS)) < 0.1)
#   expect_true(abs(cv_iid$mean.scores$neg.sCRPS + mean(sCRPS)) < 0.1)
# })