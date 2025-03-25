extract_rspde_fit <- function(rspde_fit, rspde_model) {
  result_fit <- rspde.result(rspde_fit, "field", rspde_model)
  prec_seaDist <- rspde_fit$summary.hyperpar["Precision for distSea", "mean"]
  sigma_seaDist <- sqrt(1/prec_seaDist); sigma_seaDist
  # (ngme_result <- ngme_result(m_gauss_gauss, "rw1"))

  # spde model
  tau <- result_fit$summary.tau$mean
  kappa <- result_fit$summary.kappa$mean
  sigma <- 1/tau
  kappa; sigma
  # (ngme_result <- ngme_result(m_gauss_gauss, "spde"))

  # measurement error model
  prec_mn <- rspde_fit$summary.hyperpar["Precision for the Gaussian observations", "mean"]
  sigma_mn <- sqrt(1/prec_mn); sigma_mn
  # ngme_result(m_gauss_gauss)$noise

  # fixed effect
  feff <- rspde_fit$summary.fixed$mean
  # ngme_result(m_gauss_gauss)$feff

  list(
    kappa = kappa,
    sigma = sigma,
    sigma_seaDist = sigma_seaDist,
    sigma_mn = sigma_mn,
    feff = feff
  )
}

# This tests the cross-validation of the rSPDE package against the ngme2 package
test_that("test CV (rSPDE))", {
  library(ggplot2)
  library(INLA)
  library(inlabru)
  library(splancs)
  library(viridis)
  library(sf)
  library(ngme2)
  library(rSPDE)

  data(PRprec)
  data(PRborder)
  set.seed(123)

  Y <- rowMeans(PRprec[, 3 + 1:31])
  ind <- !is.na(Y)

  # sample a subset of points
  ind <- sample(which(ind), 100)
  Y <- Y[ind]
  coords <- as.matrix(PRprec[ind, 1:2])
  alt <- PRprec$Altitude[ind]
  seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
    longlat = TRUE
  ), 1, min)

  prdomain <- fm_nonconvex_hull(coords, -0.03, -0.05, resolution = c(100, 100))
  prmesh <- fm_mesh_2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.3)

  # plot(prmesh, asp = 1, main = "")
  # lines(PRborder, col = 3)
  # points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")

  prdata_df <- data.frame(long = coords[,1], lat = coords[,2], 
                          seaDist = inla.group(seaDist), y = Y)
  prdata <- st_as_sf(prdata_df, coords = c("long", "lat"), crs = 4326)

  rspde_model <- rspde.matern(
    mesh = prmesh,
    nu = 1
  )

  cmp <- y ~ Intercept(1) + 
    distSea(seaDist, model="rw1") +
    field(geometry, model = rspde_model)

  rspde_fit <- bru(cmp, data = prdata,
    family = "gaussian",
    options = list(
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE,
      num.threads = "1:1")
  )
  rspde_fit

  rspde_result <- extract_rspde_fit(rspde_fit, rspde_model)
  rspde_result

  # Construct the ngme model (using rSPDE parameters)
  m_gauss_gauss <- ngme(
    formula = y ~ 1 +
      f(seaDist, name="rw1", model = "rw1", noise = noise_normal(
        sigma = rspde_result$sigma_seaDist
      )) +
      f(~long+lat, model = "matern", 
        mesh = prmesh, 
        name="spde", 
        theta_K = log(rspde_result$kappa),
        noise = noise_normal(
          sigma = rspde_result$sigma
        )),
    data = prdata_df,
    family = noise_normal(
      sigma = rspde_result$sigma_mn
    ),
    control_opt = control_opt(estimation = FALSE)
  )
  m_gauss_gauss
  rspde_result

  # Cross-validation of the rSPDE model
  rspde_cv <- rSPDE::cross_validation(
    list(
      rspde = rspde_fit
    ),
    cv_type = "k-fold",
    k = 5,
    print = TRUE,
    return_train_test = TRUE
  )
  
  # Given trains and tests for ngme CV
  trains <- sapply(rspde_cv$train_test, function(x) x$train)
  tests <- sapply(rspde_cv$train_test, function(x) x$test)

  # Cross-validation of the ngme model
  ngme_cv <- ngme2::cross_validation(
    list(
      gauss_gauss = m_gauss_gauss
    ),
    N = 10,
    type = "custom",
    train_idx = trains,
    test_idx = tests,
    n_gibbs_samples = 2000,
    print = TRUE,
    thining_gap = 10
  )

  # read the results
  ngme_cv$mean.scores
  rspde_cv$scores_df

  n_batches <- length(trains)

  # Create train_test_folds for rSPDE::cross_validation from ngme CV
  # This converts
  train_test_folds <- lapply(1:n_batches, function(i) list(train = list(trains[[i]]), test = list(tests[[i]])))

  # Cross-validation of the rSPDE model
  rspde_cv_2 <- rSPDE::cross_validation(
    list(
      rspde = rspde_fit
    ),
    print = TRUE,
    train_test_indexes = train_test_folds
  )

  # Compare the results
  mse_diff <- as.numeric(rspde_cv_2$mse[[1]]) - ngme_cv$mean.scores$MSE
  mae_diff <- as.numeric(rspde_cv_2$mae[[1]]) - ngme_cv$mean.scores$MAE
  crps_diff <- as.numeric(rspde_cv_2$crps[[1]]) - ngme_cv$mean.scores$neg.CRPS
  sCRPS_diff <- as.numeric(rspde_cv_2$scrps[[1]]) - ngme_cv$mean.scores$neg.sCRPS

  expect_true(abs(mse_diff) < 0.05)
  expect_true(abs(mae_diff) < 0.05)
  expect_true(abs(crps_diff) < 0.05)
  expect_true(abs(sCRPS_diff) < 0.05)
})


test_that("test iid CRPS, sCRPS", {
  set.seed(345)
  n_obs <- 100
  ys <- rnorm(n_obs, 0, 1)

  n_sim <- 1000
  y_sim1 <- rnorm(n_sim, 0, 1)
  y_sim2 <- rnorm(n_sim, 0, 1)
  CRPS <- numeric(n_obs)
  sCRPS <- numeric(n_obs)
  for (i in 1:n_obs) {
    y <- ys[i]
    E_sim_data <- mean(abs(y_sim1 - y))
    E_sim_sim <- mean(abs(y_sim1 - y_sim2))
    CRPS[i] <- 0.5 * E_sim_sim - E_sim_data
    sCRPS[i] <- -E_sim_data / E_sim_sim - 0.5 * log(E_sim_sim)
  }
  mean(CRPS)
  mean(sCRPS)

  ng <- ngme(
    y ~ 0,
    data = data.frame(y = ys),
    family = noise_normal(
      sigma = 1
    ),
    control_opt = control_opt(estimation = FALSE)
  )
  ng

  # Cross-validation of the ngme model
  cv_iid <- ngme2::cross_validation(
    list(
      ng = ng
    ),
    type = "loo",
    n_gibbs_samples = 1000,
    print = TRUE,
    thining_gap = 0
  )
  cv_iid

  expect_true(abs(cv_iid$mean.scores$neg.CRPS + mean(CRPS)) < 0.1)
  expect_true(abs(cv_iid$mean.scores$neg.sCRPS + mean(sCRPS)) < 0.1)
})
