# Benchmark model suite for ngme2.
#
# A small spread of shapes that between them exercise the paths where
# performance work tends to land: Gaussian against non-Gaussian latent noise
# (the latter runs n_gibbs passes instead of one Rao-Blackwell pass), several
# operator families, more than one latent component, non-Gaussian measurement
# noise, and replicates.
#
# ngme() resolves its formula in the CALLING frame, so everything the formula
# names has to be visible there. Hence prefixed globals and plain top-level fit
# functions rather than closures -- the same shape the package's own tests use.

suppressMessages({library(ngme2); library(fmesher)})

# Posterior sampling runs after estimation and is irrelevant to fit cost, so it
# is cut to the minimum: these timings should measure the fit. n_gibbs_samples
# lives in control_ngme, not control_opt.
BM_GIBBS <- 5L
# Spacetime stabilization: TRUE is the package default. The analytic first
# derivative covers it; the second does not and falls back to differencing, so
# this is exposed to keep that cost measurable.
BM_STAB <- as.logical(Sys.getenv('NGME2_BM_STAB', 'TRUE'))
bm_ngme <- function() control_ngme(n_post_samples = 1L, post_burnin = 0L,
                                   n_gibbs_samples = BM_GIBBS)

## ---- AR(1), Gaussian ------------------------------------------------------
setup_ar1 <- function(nt = 2000, rho = 0.7, sig = 1, meas = 0.4, seed = 1) {
  set.seed(seed); idx <- 1:nt
  W <- as.numeric(simulate(f(map = idx, model = ar1(mesh = idx, rho = rho),
                             noise = noise_normal(sigma = sig)), seed = seed)[[1]])
  ar_dat <<- data.frame(Y = W + rnorm(nt, 0, meas), idx = idx); ar_idx <<- idx
  c(rho = rho, sigma = sig, meas_sigma = meas)
}
fit_ar1 <- function(...)
  ngme(Y ~ 0 + f(map = idx, model = ar1(mesh = ar_idx), noise = noise_normal(), name = "ar"),
       data = ar_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- AR(1), NIG latent ----------------------------------------------------
setup_ar1_nig <- function(nt = 2000, rho = 0.7, sig = 1, mu = 1, nu = 1, meas = 0.4, seed = 1) {
  set.seed(seed); idx <- 1:nt
  W <- as.numeric(simulate(f(map = idx, model = ar1(mesh = idx, rho = rho),
             noise = noise_nig(mu = mu, sigma = sig, nu = nu)), seed = seed)[[1]])
  nar_dat <<- data.frame(Y = W + rnorm(nt, 0, meas), idx = idx); nar_idx <<- idx
  c(rho = rho, mu = mu, sigma = sig, nu = nu, meas_sigma = meas)
}
fit_ar1_nig <- function(...)
  ngme(Y ~ 0 + f(map = idx, model = ar1(mesh = nar_idx), noise = noise_nig(), name = "ar"),
       data = nar_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- RW(1), Gaussian: a non-square operator -------------------------------
setup_rw1 <- function(nt = 1500, sig = 1, meas = 0.4, seed = 1) {
  set.seed(seed); idx <- 1:nt
  W <- cumsum(rnorm(nt, 0, sig))
  rw_dat <<- data.frame(Y = W + rnorm(nt, 0, meas), idx = idx); rw_idx <<- idx
  c(sigma = sig, meas_sigma = meas)
}
fit_rw1 <- function(...)
  ngme(Y ~ 0 + f(map = idx, model = rw1(rw_idx), noise = noise_normal(), name = "rw"),
       data = rw_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Matern 2-D, Gaussian and NIG -----------------------------------------
setup_matern <- function(n = 900, kappa = 1.5, sig = 1, meas = 0.4, seed = 1,
                         edge = c(1, 2), cut = 0.3) {
  set.seed(seed); loc <- cbind(runif(n, 0, 10), runif(n, 0, 10))
  mt_mesh <<- fm_mesh_2d(loc = loc, max.edge = edge, cutoff = cut); mt_loc <<- loc
  W <- as.numeric(simulate(f(map = loc, model = matern(mesh = mt_mesh, kappa = kappa),
             noise = noise_normal(sigma = sig)), seed = seed)[[1]])
  mt_dat <<- data.frame(Y = W + rnorm(n, 0, meas))
  c(kappa = kappa, sigma = sig, meas_sigma = meas)
}
fit_matern <- function(...)
  ngme(Y ~ 0 + f(map = mt_loc, model = matern(mesh = mt_mesh), noise = noise_normal(), name = "sp"),
       data = mt_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

setup_matern_nig <- function(n = 900, kappa = 1.5, sig = 1, mu = 1, nu = 1, meas = 0.4,
                             seed = 1, edge = c(1, 2), cut = 0.3) {
  set.seed(seed); loc <- cbind(runif(n, 0, 10), runif(n, 0, 10))
  mn_mesh <<- fm_mesh_2d(loc = loc, max.edge = edge, cutoff = cut); mn_loc <<- loc
  W <- as.numeric(simulate(f(map = loc, model = matern(mesh = mn_mesh, kappa = kappa),
             noise = noise_nig(mu = mu, sigma = sig, nu = nu)), seed = seed)[[1]])
  mn_dat <<- data.frame(Y = W + rnorm(n, 0, meas))
  c(kappa = kappa, mu = mu, sigma = sig, nu = nu, meas_sigma = meas)
}
fit_matern_nig <- function(...)
  ngme(Y ~ 0 + f(map = mn_loc, model = matern(mesh = mn_mesh), noise = noise_nig(), name = "sp"),
       data = mn_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Tensor product AR(1) x Matern, Gaussian and NIG ----------------------
setup_tp <- function(nsp = 60, ntime = 12, rho = 0.6, kappa = 2, sig = 1, meas = 0.4, seed = 1) {
  set.seed(seed); loc <- cbind(runif(nsp, 0, 10), runif(nsp, 0, 10))
  tp_mesh_s <<- fm_mesh_2d(loc = loc, max.edge = c(1.5, 4), cutoff = 0.8); tp_nt <<- ntime
  df <- data.frame(time = rep(1:ntime, each = nsp),
                   long = rep(loc[, 1], ntime), lat = rep(loc[, 2], ntime))
  spec <- f(map = list(df$time, ~ long + lat),
            model = tp(first = ar1(mesh = 1:ntime, rho = rho),
                       second = matern(mesh = tp_mesh_s, kappa = kappa)),
            noise = noise_normal(sigma = sig), data = df)
  df$Y <- as.numeric(simulate(spec, seed = seed)[[1]]) + rnorm(nrow(df), 0, meas)
  tp_dat <<- df
  c(rho = rho, kappa = kappa, sigma = sig, meas_sigma = meas)
}
fit_tp <- function(...)
  ngme(Y ~ 0 + f(map = list(time, ~ long + lat),
                 model = tp(first = ar1(mesh = 1:tp_nt), second = matern(mesh = tp_mesh_s)),
                 noise = noise_normal(), name = "tp"),
       data = tp_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

setup_tp_nig <- function(nsp = 60, ntime = 12, rho = 0.6, kappa = 2, sig = 1,
                         mu = 1, nu = 1, meas = 0.4, seed = 1,
                         edge = c(1.5, 4), cut = 0.8) {
  set.seed(seed); loc <- cbind(runif(nsp, 0, 10), runif(nsp, 0, 10))
  tn_mesh_s <<- fm_mesh_2d(loc = loc, max.edge = edge, cutoff = cut); tn_nt <<- ntime
  df <- data.frame(time = rep(1:ntime, each = nsp),
                   long = rep(loc[, 1], ntime), lat = rep(loc[, 2], ntime))
  spec <- f(map = list(df$time, ~ long + lat),
            model = tp(first = ar1(mesh = 1:ntime, rho = rho),
                       second = matern(mesh = tn_mesh_s, kappa = kappa)),
            noise = noise_nig(mu = mu, sigma = sig, nu = nu), data = df)
  df$Y <- as.numeric(simulate(spec, seed = seed)[[1]]) + rnorm(nrow(df), 0, meas)
  tn_dat <<- df
  c(rho = rho, kappa = kappa, mu = mu, sigma = sig, nu = nu, meas_sigma = meas)
}
fit_tp_nig <- function(...)
  ngme(Y ~ 0 + f(map = list(time, ~ long + lat),
                 model = tp(first = ar1(mesh = 1:tn_nt), second = matern(mesh = tn_mesh_s)),
                 noise = noise_nig(), name = "tp"),
       data = tn_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Non-separable spacetime (advection-diffusion) ------------------------
setup_st <- function(nsp = 250, ntime = 6, kappa = 1.5, sig = 1, meas = 0.4, seed = 1,
                     edge = c(1, 2.5), cut = 0.5) {
  set.seed(seed); loc <- cbind(runif(nsp, 0, 5), runif(nsp, 0, 5))
  st_mesh <<- list(fm_mesh_1d(1:ntime),
                   fm_mesh_2d(loc = loc, max.edge = edge, cutoff = cut))
  df <- data.frame(time = rep(1:ntime, each = nsp),
                   c1 = rep(loc[, 1], ntime), c2 = rep(loc[, 2], ntime))
  spec <- f(list(df$time, cbind(df$c1, df$c2)),
            model = spacetime(mesh = st_mesh, alpha = 2, kappa = kappa,
                              stabilization = BM_STAB,
                              theta_gamma_x = 0.3, theta_gamma_y = -0.2),
            noise = noise_normal(sigma = sig), data = df)
  df$Y <- as.numeric(simulate(spec, seed = seed)[[1]]) + rnorm(nrow(df), 0, meas)
  st_dat <<- df
  c(kappa = kappa, gamma_x = 0.3, gamma_y = -0.2, sigma = sig, meas_sigma = meas)
}
fit_st <- function(...)
  ngme(Y ~ 0 + f(list(time, cbind(c1, c2)),
                 model = spacetime(mesh = st_mesh, alpha = 2, stabilization = BM_STAB),
                 noise = noise_normal(), name = "st"),
       data = st_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Two latent components in one formula ---------------------------------
setup_two <- function(nt = 1200, nsp = 400, rho = 0.6, kappa = 1.5,
                      s1 = 1, s2 = 0.8, meas = 0.4, seed = 1) {
  set.seed(seed); idx <- 1:nt; loc <- cbind(runif(nt, 0, 10), runif(nt, 0, 10))
  tw_mesh <<- fm_mesh_2d(loc = loc[seq_len(nsp), , drop = FALSE],
                         max.edge = c(1.5, 4), cutoff = 0.8)
  tw_idx <<- idx; tw_loc <<- loc
  W1 <- as.numeric(simulate(f(map = idx, model = ar1(mesh = 1:nt, rho = rho),
             noise = noise_normal(sigma = s1)), seed = seed)[[1]])
  W2 <- as.numeric(simulate(f(map = loc, model = matern(mesh = tw_mesh, kappa = kappa),
             noise = noise_normal(sigma = s2)), seed = seed + 1)[[1]])
  tw_dat <<- data.frame(Y = W1 + W2 + rnorm(nt, 0, meas))
  c(rho = rho, sigma1 = s1, kappa = kappa, sigma2 = s2, meas_sigma = meas)
}
fit_two <- function(...)
  ngme(Y ~ 0 + f(map = tw_idx, model = ar1(mesh = 1:length(tw_idx)),
                 noise = noise_normal(), name = "t") +
                f(map = tw_loc, model = matern(mesh = tw_mesh),
                 noise = noise_normal(), name = "s"),
       data = tw_dat, family = "normal", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Gaussian latent, NIG MEASUREMENT noise -------------------------------
setup_nigmeas <- function(n = 900, kappa = 1.5, sig = 1, meas = 0.4,
                          mu_e = 1, nu_e = 1, seed = 1) {
  set.seed(seed); loc <- cbind(runif(n, 0, 10), runif(n, 0, 10))
  nm_mesh <<- fm_mesh_2d(loc = loc, max.edge = c(1, 2), cutoff = 0.3); nm_loc <<- loc
  W <- as.numeric(simulate(f(map = loc, model = matern(mesh = nm_mesh, kappa = kappa),
             noise = noise_normal(sigma = sig)), seed = seed)[[1]])
  V <- ngme2:::rGIG_cpp(rep(-0.5, n), rep(nu_e, n), rep(nu_e, n), seed)
  e <- mu_e * (V - 1) + meas * sqrt(V) * rnorm(n)
  nm_dat <<- data.frame(Y = W + e)
  c(kappa = kappa, sigma = sig, mu = mu_e, meas_sigma = meas, nu = nu_e)
}
fit_nigmeas <- function(...)
  ngme(Y ~ 0 + f(map = nm_loc, model = matern(mesh = nm_mesh),
                 noise = noise_normal(), name = "sp"),
       data = nm_dat, family = "nig", control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Replicates, with fixed effects ---------------------------------------
setup_repl <- function(nt = 600, nrep = 3, rho = 0.7, sig = 1, meas = 0.4, seed = 1) {
  set.seed(seed); idx <- 1:nt
  Ws <- lapply(seq_len(nrep), function(r)
    as.numeric(simulate(f(map = idx, model = ar1(mesh = idx, rho = rho),
               noise = noise_normal(sigma = sig)), seed = seed + r)[[1]]))
  x1 <- runif(nt * nrep)
  rp_dat <<- data.frame(Y = unlist(Ws) + 1.5 + 2 * x1 + rnorm(nt * nrep, 0, meas),
                        x1 = x1, idx = rep(idx, nrep), rep = rep(seq_len(nrep), each = nt))
  rp_idx <<- idx
  c(rho = rho, sigma = sig, meas_sigma = meas, beta0 = 1.5, beta1 = 2)
}
fit_repl <- function(...)
  ngme(Y ~ x1 + f(map = idx, model = ar1(mesh = rp_idx), noise = noise_normal(), name = "ar"),
       data = rp_dat, replicate = rp_dat$rep, family = "normal",
       control_opt = control_opt(...), control_ngme = bm_ngme())

## ---- Larger variants, run for only a few iterations ------------------------
# The small models above answer "did this change the numbers" and "is a pass
# cheaper". They cannot answer "does the change still help when the system is
# ten times bigger", which is where sparse work usually stops behaving. These
# are sized so one pass is expensive and the run is short.

setup_matern_big <- function(n = 6000, ...)
  setup_matern(n = n, edge = c(0.22, 0.6), cut = 0.09, ...)
fit_matern_big <- function(...) fit_matern(...)

setup_matern_nig_big <- function(n = 6000, ...)
  setup_matern_nig(n = n, edge = c(0.22, 0.6), cut = 0.09, ...)
fit_matern_nig_big <- function(...) fit_matern_nig(...)

setup_tp_nig_big <- function(nsp = 500, ntime = 16, ...)
  setup_tp_nig(nsp = nsp, ntime = ntime, edge = c(0.7, 2), cut = 0.35, ...)
fit_tp_nig_big <- function(...) fit_tp_nig(...)

setup_st_big <- function(nsp = 900, ntime = 14, ...)
  setup_st(nsp = nsp, ntime = ntime, edge = c(0.28, 0.8), cut = 0.12, ...)
fit_st_big <- function(...) fit_st(...)

# `iters` overrides the harness default for that model; `dim` is filled in by
# the harness from the fitted object and reported so the scaling is visible.
BM_MODELS <- list(
  ar1        = list(setup = setup_ar1,        fit = fit_ar1),
  ar1_nig    = list(setup = setup_ar1_nig,    fit = fit_ar1_nig),
  rw1        = list(setup = setup_rw1,        fit = fit_rw1),
  matern     = list(setup = setup_matern,     fit = fit_matern),
  matern_nig = list(setup = setup_matern_nig, fit = fit_matern_nig),
  tp         = list(setup = setup_tp,         fit = fit_tp),
  tp_nig     = list(setup = setup_tp_nig,     fit = fit_tp_nig),
  spacetime  = list(setup = setup_st,         fit = fit_st),
  two_field  = list(setup = setup_two,        fit = fit_two),
  nig_meas   = list(setup = setup_nigmeas,    fit = fit_nigmeas),
  replicates = list(setup = setup_repl,       fit = fit_repl),
  # larger, few iterations
  matern_big     = list(setup = setup_matern_big,     fit = fit_matern_big,     iters = 20L, large = TRUE),
  matern_nig_big = list(setup = setup_matern_nig_big, fit = fit_matern_nig_big, iters = 20L, large = TRUE),
  tp_nig_big     = list(setup = setup_tp_nig_big,     fit = fit_tp_nig_big,     iters = 20L, large = TRUE),
  spacetime_big  = list(setup = setup_st_big,         fit = fit_st_big,         iters = 20L, large = TRUE))
