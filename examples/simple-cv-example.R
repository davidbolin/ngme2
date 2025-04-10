library(ngme2)
library(Matrix)

# Define seeds before using it
seeds = 1:10

result = data.frame(
  mae_gauss = double(length(seeds)), 
  mae_gauss_rao = double(length(seeds)), 
  mae_nig = double(length(seeds)), 
  mae_nig_rao = double(length(seeds)),
  mse_gauss = double(length(seeds)), 
  mse_gauss_rao = double(length(seeds)), 
  mse_nig = double(length(seeds)),
  mse_nig_rao = double(length(seeds))
)
# Use parallel processing for the seeds
library(parallel)

# Function to process a single seed
process_seed <- function(seed) {
  set.seed(seed)
  n_obs <- 500
  rho <- 0.5
  K <- ngme2::ar1(1:n_obs, rho = rho)$K
  h <- rep(1, n_obs)
  A_obs <- Matrix::Diagonal(n_obs)

  # library(fmesher)
  # mesh <- fmesher::fm_mesh_1d(1:n_obs)
  # A_obs <- fm_basis(mesh, 1:n_obs)
  # fem <- fm_fem(mesh)
  # C = fem$c0
  # G = fem$g1
  # kappa <- 1
  # K <- kappa^2 * C + G
  # h <- Matrix::diag(fem$c0) 

  # Common parameters across the entire code
  nu <- 0.2
  mu <- 2
  delta <- -mu * h
  sigma <- 1
  sigma_eps <- 1

  # Generate true data
  trueV <- ngme2::rig(n_obs, nu, nu * h^2, seed = seed)
  mynoise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
  # trueW2 <- Reduce(function(x,y){y + rho*x}, mynoise, accumulate = TRUE)
  trueW <- solve(K, mynoise)
  Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)

  # Cholesky-based sampler for N(Q^{-1} m, Q^{-1})
  sample_mvn_cholesky <- function(Q, m) {
    R <- chol(Q)
    mu <- solve(Q, m)
    n <- length(m)
    z <- rnorm(n)
    sample <- mu + solve(R, z)
    return(list(sample = sample, mu = mu))
  }

  burnin <- 10
  iter <- 100

  loo_predict <- function(n_obs, seed, iter=iter, burnin=burnin, gaussian=FALSE, indices=1:n_obs, K, Y, sigma, sigma_eps, mu, h, nu, A_obs) {
    set.seed(seed)
    if(gaussian){
      cat(paste("Seed", seed, "- Gaussian\n"))
    } else{
      cat(paste("Seed", seed, "- NIG\n"))
    }

    pred <- double(length(indices))
    pred_rao <- double(length(indices))
    pred_rao_avg <- double(length(indices))
    for (i in indices) {
      A <- A_obs
      A <- A[-i,]
      Y_new <- Y[-i]
      post_W <- matrix(0, nrow=n_obs, ncol=iter)
      post_m <- matrix(0, nrow=n_obs, ncol=iter)
      post_m_avg <- matrix(0, nrow=n_obs, ncol=iter)

      if (gaussian){
          V <- rep(1, n_obs)
      } else{
          V <- ngme2::rig(n_obs, nu, nu * h^2)
      }
      cat(paste0("Seed ", seed, " - i = ", i, "\n"))
      
      for (j in 1:iter) {
        Q = sigma^(-2) * t(K) %*% diag(1/V) %*% K + sigma_eps^(-2) * t(A) %*% A
        m = sigma^(-2) * t(K) %*% diag(1/V) %*% (mu*(V-h)) + sigma_eps^(-2) * t(A) %*% Y_new
        m = as.numeric(m)
        
        # sample N(Q^-1 m, Q^-1)
        aux_W <- sample_mvn_cholesky(Q, m)
        post_W[,j] <- aux_W$sample
        post_m[,j] <- aux_W$mu

        if (gaussian){
          V <- rep(1, n_obs)
        } else{
          p <- rep(-1, n_obs)
          a <- rep(nu, n_obs)
          a <- a + (mu / sigma)^2
          M <- K %*% post_W[,j] + mu * h
          M <- as.vector(M)
          b <- nu * h^2
          b <- b + (M / sigma)^2

          V <- ngme2::rgig(n_obs, p, a, b)         
        }

      }
      
      pred[which(indices==i)] <- mean(post_W[i,(burnin+1):iter])
      pred_rao[which(indices==i)] <- mean(post_m[i,(burnin+1):iter])
      
      cat(paste0("Seed ", seed, " - Y[",i,"] = ", Y[i], "; pred = ",pred[which(indices==i)], "; pred_rao = ",pred_rao[which(indices==i)], "\n"))      
    }
    
    list(
      pred = pred, 
      pred_rao = pred_rao
    )
  }

  # Run Gaussian and Non-Gaussian versions in parallel
  indices = c(300)
  # Create a list of parameters for both models
  model_params <- list(
    list(n_obs=n_obs, seed=seed, iter=iter, burnin=burnin, gaussian=TRUE, indices=indices, 
        K=K, Y=Y, sigma=sigma, sigma_eps=sigma_eps, mu=mu, h=h, nu=nu, A_obs=A_obs),
    list(n_obs=n_obs, seed=seed, iter=iter, burnin=burnin, gaussian=FALSE, indices=indices, 
        K=K, Y=Y, sigma=sigma, sigma_eps=sigma_eps, mu=mu, h=h, nu=nu, A_obs=A_obs)
  )

  # Run both models in parallel
  results <- mclapply(model_params, function(params) {
    do.call(loo_predict, params)
  }, mc.cores = 2)

  # print(results)

  # Extract results
  post_W_gauss <- results[[1]]
  post_W_nig <- results[[2]]

  # Process Gaussian results
  pred_gauss <- post_W_gauss[["pred"]]
  pred_gauss_rao <- post_W_gauss[["pred_rao"]]
  mae_gauss <- mean(abs(pred_gauss - Y[indices]))
  mse_gauss <- mean((pred_gauss - Y[indices])^2)
  mae_gauss_rao <- mean(abs(pred_gauss_rao - Y[indices]))
  mse_gauss_rao <- mean((pred_gauss_rao - Y[indices])^2)

  # Process Non-Gaussian results
  pred_nig <- post_W_nig[["pred"]]
  pred_nig_rao <- post_W_nig[["pred_rao"]]
  mae_nig <- mean(abs(pred_nig - Y[indices]))
  mse_nig <- mean((pred_nig - Y[indices])^2)
  mae_nig_rao <- mean(abs(pred_nig_rao - Y[indices]))
  mse_nig_rao <- mean((pred_nig_rao - Y[indices])^2)

  return(c(
    mae_gauss, 
    mae_gauss_rao, 
    mae_nig,
    mae_nig_rao, 
    mse_gauss, 
    mse_gauss_rao, 
    mse_nig, 
    mse_nig_rao
  ))
}

# Determine number of cores to use (leave one free for system)
num_cores <- min(detectCores() - 1, length(seeds))
cat(paste("Running with", num_cores, "cores\n"))

# Process all seeds in parallel
results_list <- mclapply(seeds, process_seed, mc.cores = num_cores)

# Convert results list to data frame
for (i in 1:length(seeds)) {
  result[i, ] <- results_list[[i]]
}

# res_ind_400_seed_1000 <- colMeans(result)
# res_ind_400_seed_1000

result_ind_300_seed_300 <- colMeans(result)
result_ind_300_seed_300


