#' Simulates an ARMA(p,q) time series.
#'
#' Extends the ARMA(1,1) function to support vectorized 'ar' and 'ma' parameters.
#'
#' @param n integer, The length of the series to be generated.
#' @param ar numeric vector, The AR(p) coefficients (phi_1, phi_2, ... phi_p). 
#'           Set to NULL or numeric(0) for p=0.
#' @param ma numeric vector, The MA(q) coefficients (theta_1, theta_2, ... theta_q).
#'           Set to NULL or numeric(0) for q=0.
#' @param sigma numeric, Standard deviation of the noise (used for 'normal' noise).
#' @param mu numeric, Parameter 'mu' for 'nig' noise.
#' @param nu numeric, Parameter 'nu' for 'nig' noise.
#' @param noise string, Type of innovation: "normal" or "nig".
#' @param seed integer, Random seed for reproducibility.
#'
#' @return A numeric vector of length n, representing the simulated ARMA(p,q) series.
#' @note 
#' The 'nig' noise implementation relies on an 'rnig' function being defined 
#' in the environment, which accepts 'delta', 'mu', 'sigma', 'nu' as arguments.
#'
simulate_armapq <- function(
  n = 100, 
  ar = NULL,
  ma = NULL,
  sigma = 1.0, 
  mu = 0, nu = 1, 
  noise = "normal", 
  seed = 123
) {
  set.seed(seed)
  
  # --- 1. Generate Innovations ---
  if (noise == "normal") {
    e <- rnorm(n, mean = 0, sd = sigma)
  } else if (noise == "nig") {
    # Using the exact call provided in the original code
    e <- rnig(n, delta = -mu, mu = mu, sigma = sigma, nu = nu)
  } else {
    stop("noise must be either 'normal' or 'nig'")
  }
  
  # --- 2. Get p and q ---
  # length(NULL) returns 0, which is correct for ARMA(0,q) or ARMA(p,0)
  p <- length(ar)
  q <- length(ma)
  
  # --- 3. Initialize series ---
  x <- numeric(n)
  
  # --- 4. Iterate ARMA(p,q) ---
  # This loop correctly handles the "warm-up" phase where t <= p or t <= q
  
  for (t in 1:n) {
    ar_part <- 0
    ma_part <- 0
    
    # --- AR(p) Component ---
    if (p > 0 && t > 1) {
      # Determine the number of available AR lags
      # e.g., at t=2, even if p=3, we only use 1 lag (x[t-1])
      n_ar_lags <- min(p, t - 1)
      
      # Extract available lagged x values
      # Note: Lags are x[t-1], x[t-2], ...
      past_x <- x[ (t - 1) : (t - n_ar_lags) ]
      
      # Extract corresponding AR coefficients (phi_1, phi_2, ...)
      ar_coeffs <- ar[1:n_ar_lags]
      
      # Calculate the contribution of the AR component
      # AR part = phi_1*x[t-1] + phi_2*x[t-2] + ...
      ar_part <- sum(ar_coeffs * past_x)
    }
    
    # --- MA(q) Component ---
    if (q > 0 && t > 1) {
      # Determine the number of available MA lags
      n_ma_lags <- min(q, t - 1)
      
      # Extract available lagged e values
      # Note: Lags are e[t-1], e[t-2], ...
      past_e <- e[ (t - 1) : (t - n_ma_lags) ]
      
      # Extract corresponding MA coefficients (theta_1, theta_2, ...)
      ma_coeffs <- ma[1:n_ma_lags]
      
      # Calculate the contribution of the MA component
      # MA part = theta_1*e[t-1] + theta_2*e[t-2] + ...
      ma_part <- sum(ma_coeffs * past_e)
    }
    
    # --- Combination ---
    # x[t] = AR_part + e[t] + MA_part
    x[t] <- ar_part + e[t] + ma_part
  }
  
  # --- 5. Return result ---
  # ts.plot(x, main = "Simulated ARMA(p,q) Series")
  
  invisible(x)
}


test_that("test fit ARMA(2,2) with normal noise (precond)", {
  n <- 1000
  rho <- c(-0.5, 0.3)
  phi <- c(0.8, 0.2)
  # phi <- c(0.8, -0.5)
  # phi <- c(0.0, 0.0)
  sigma <- 4
  sigma_e <- 1
  mesh <- fmesher::fm_mesh_1d(1:n)
  w <- simulate_armapq(n = n, ar = rho, ma = phi, sigma = sigma, seed = 123)
  sd(w)
  mean(w)
  y <- w + rnorm(n, mean = 0, sd = sigma_e)

  fit <- ngme(
    formula = y ~ 0 + f(
      x, model = "arma", 
      ar_order = 2, 
      ma_order = 2,
      noise = noise_normal()
    ),
    data = data.frame(x = 1:n, y = y),
    control_ngme = control_ngme(n_gibbs_samples=4),
    control_opt = control_opt(
      seed = 123,
      burnin = 100,
      iterations = 2000,
      stop_points = 1,
      # verbose=TRUE,
      # rao_blackwellization = TRUE,
      # optimizer = precond_sgd(),
      n_parallel_chain = 4,
      start_sd = 0.01
    )
  )
  fit
  traceplot(fit, hline=sigma_e)
  traceplot(fit, "field1", hline=c(rho, phi, sigma))

  arma_result <- as.numeric(ngme_result(fit)$field1)
  score <- abs((arma_result -  c(rho, phi, sigma)) / c(rho, phi, sigma))
  expect_true(all(score < 0.4))
})


test_that("test fit ARMA(1,1) with nig noise (RB, preconditioner)", {
  n <- 1000
  rho <- 0.5
  phi <- -0.8
  mu <- 2
  sigma <- 4
  nu <- 0.3
  sigma_e <- 2
  mesh <- fmesher::fm_mesh_1d(1:n)
  w <- simulate_armapq(n = n, ar = rho, ma = phi, 
    mu = mu, sigma = sigma, nu = nu, seed = 123, noise="nig")
  sd(w)
  mean(w)
  y <- w + rnorm(n, mean = 0, sd = sigma_e)

  fit <- ngme(
    formula = y ~ 0 + f(
      map=x, model = "arma", 
      ar_order = 1, 
      ma_order = 1,
      ar = 0.5,
      ma = 0.0,
      noise = noise_nig()
    ),
    data = data.frame(x = 1:n, y = y),
    control_ngme = control_ngme(n_gibbs_samples=4),
    control_opt = control_opt(
      seed = 123,
      burnin = 100,
      iterations = 500,
      stop_points = 1,
      # rao_blackwellization = TRUE,
      optimizer = precond_sgd(),
      n_parallel_chain = 4,
      start_sd = 0
    )
  )
  fit
  traceplot(fit, hline=sigma_e)
  traceplot(fit, "field1", hline=c(rho, phi, mu, sigma, nu))

  arma_result <- as.numeric(ngme_result(fit)$field1)
  score <- abs((arma_result -  c(rho, phi, mu, sigma, nu)) / c(rho, phi, mu, sigma, nu))
  score
  expect_true(all(score < c(0.1, 0.1, 0.3, 0.3, 0.7)))

  est_sigma_e <- ngme_result(fit)$data$sigma
  expect_true(all(abs((est_sigma_e - sigma_e) / sigma_e) < 0.2))
})
