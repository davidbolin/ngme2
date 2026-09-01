#' Predict function of ngme2
#' predict using ngme after estimation
#'
#' @param object a ngme object
#' @param map a named list (or dataframe) of the locations to make the prediction
#' @param data a data.frame or matrix of covariates (used for fixed effects)
#'  names(loc) corresponding to the name each latent model
#'  vector or matrix (n * 2) for spatial coords
#' @param type what type of prediction, c("fe", "lp", <model_name>)
#' "fe" is fixed effect prediction
#' <model_name> is prediction of a specific model
#' "lp" is linear predictor (including fixed effect and all sub-models)
#' "response" is the linear predictor plus a fresh measurement-noise draw
#' @param group which filed to predict
#'   (used for bivariate model, should be of same length as map)
#' @param estimator what type of estimator. Options include:
#'   - "mean", "median", "mode", "sd": standard estimators
#'   - "0.XXXq": any quantile specified as probability (e.g., "0.025q", "0.5q", "0.975q")
#' @param sampling_size size of posterior sampling
#' @param burnin_size size of posterior burnin
#' @param seed random seed. Defaults to a seed drawn from the current R
#'   random number stream, so \code{set.seed()} makes the result reproducible.
#' @param replicate which replicate to predict for. A replicate LEVEL (the name
#'   used in the fit, e.g. `2007` or `"site_a"`) or a positive integer position
#'   in `object$replicates`. Defaults to the first replicate.
#' @param train_idx optional vector of training indices to use for posterior sampling.
#'   If provided, only these indices from the original data will be used for training,
#'   similar to cross-validation. If NULL, uses all original training data.
#' @param chain_combine how to combine multiple optimization chains:
#'   \itemize{
#'     \item \code{"param_mean"}: default behavior using the fitted object parameters.
#'     \item \code{"predictive_average"}: run prediction for each optimization chain and
#'       average predictions across chains.
#'   }
#' @param return_samples logical; when `TRUE`, attach sample draws for the
#'   requested output in `attr(ret, "samples")`. For `type = "response"`, the
#'   attached samples are response predictive draws.
#' @param ... additional arguments (currently unused)
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
predict.ngme <- function(
    object,
    map,
    data = NULL,
    type = "lp",
    group = NULL,
    estimator = c("mean", "sd", "0.05q", "0.95q", "median", "mode"),
    sampling_size = 500,
    burnin_size = 100,
    seed = ngme_random_seed(),
    replicate = NULL,
    train_idx = NULL,
    chain_combine = c("param_mean", "predictive_average"),
    return_samples = FALSE,
    ...) {
  chain_combine <- match.arg(chain_combine)
  stopifnot(
    sampling_size > 0,
    "Make sure the object is of class 'ngme'." = inherits(object, "ngme")
  )

  if (chain_combine == "predictive_average") {
    chain_fits <- resolve_ngme_chain_fits(object)
    if (length(chain_fits) > 1) {
      chain_preds <- lapply(seq_along(chain_fits), function(i) {
        predict_ngme_param_mean(
          object = chain_fits[[i]],
          map = map,
          data = data,
          type = type,
          group = group,
          estimator = estimator,
          sampling_size = sampling_size,
          burnin_size = burnin_size,
          seed = seed + i - 1L,
          replicate = replicate,
          train_idx = train_idx,
          return_samples = return_samples
        )
      })

      ret <- list()
      for (est in names(chain_preds[[1]])) {
        ret[[est]] <- rowMeans(do.call(cbind, lapply(chain_preds, `[[`, est)))
      }

      samples_list <- lapply(chain_preds, function(x) attr(x, "samples"))
      if (all(vapply(samples_list, is.matrix, logical(1)))) {
        attr(ret, "samples") <- do.call(cbind, samples_list)
      }

      return(ret)
    }
  }

  predict_ngme_param_mean(
    object = object,
    map = map,
    data = data,
    type = type,
    group = group,
    estimator = estimator,
    sampling_size = sampling_size,
    burnin_size = burnin_size,
    seed = seed,
    replicate = replicate,
    train_idx = train_idx,
    return_samples = return_samples
  )
}


# Resolve a user-supplied `replicate` to an index into object$replicates.
resolve_predict_replicate <- function(object, replicate) {
  levels_avail <- names(object$replicates)
  n <- length(object$replicates)
  if (is.null(replicate)) return(1L)
  if (length(replicate) != 1L) {
    stop("`replicate` must be a single replicate level or position.", call. = FALSE)
  }
  # a whole number is a position; anything else is matched as a level
  if (is.numeric(replicate) && !is.na(replicate) &&
      isTRUE(all.equal(replicate, round(replicate)))) {
    idx <- as.integer(round(replicate))
    if (idx < 1L || idx > n) {
      stop("`replicate` position ", idx, " is out of range: the fit has ", n,
           " replicate(s)",
           if (!is.null(levels_avail)) paste0(" (", paste(levels_avail, collapse = ", "), ")"),
           ".", call. = FALSE)
    }
    return(idx)
  }
  if (is.null(levels_avail)) {
    stop("This fit's replicates are unnamed; pass `replicate` as a position.",
         call. = FALSE)
  }
  idx <- match(as.character(replicate), levels_avail)
  if (is.na(idx)) {
    stop("`replicate` = ", sQuote(as.character(replicate)),
         " is not a replicate of this fit. Available: ",
         paste(levels_avail, collapse = ", "), ".", call. = FALSE)
  }
  idx
}

predict_ngme_param_mean <- function(
    object,
    map,
    data = NULL,
    type = "lp",
    group = NULL,
    estimator = c("mean", "sd", "0.05q", "0.95q", "median", "mode"),
    sampling_size = 500,
    burnin_size = 100,
    seed = ngme_random_seed(),
    replicate = NULL,
    train_idx = NULL,
    return_samples = FALSE) {
  fm <- attr(object, "fit")$formula
  ngme <- object$replicates[[resolve_predict_replicate(object, replicate)]]
  seed_int <- normalize_prediction_seed(seed)
  type_names <- prediction_type_names(type, ngme)
  requested_model_names <- prediction_requested_model_names(type_names, ngme)

  # If train_idx is provided, subset the data for posterior sampling
  if (!is.null(train_idx)) {
    # Validate train_idx
    n_data <- attr(object, "fit")$n_data
    if (any(train_idx > n_data) || any(train_idx < 1)) {
      stop("train_idx out of bounds. Data has ", n_data, " observations.")
    }

    # Subset training data (similar to cross-validation logic)
    ngme$X <- ngme$X[train_idx, , drop = FALSE]
    ngme$Y <- ngme$Y[train_idx]
    ngme$noise <- subset_noise(
      ngme$noise,
      sub_idx = train_idx, compute_corr = TRUE
    )

    # Subset A matrices for training data
    for (i in seq_along(ngme$models)) {
      ngme$models[[i]]$A <- ngme$models[[i]]$A[train_idx, , drop = FALSE]
    }
  }

  samples_W <- sampling_cpp(ngme,
    n = sampling_size,
    n_burnin = burnin_size,
    posterior = TRUE,
    seed = seed_int
  )[["W"]]

  # Convert samples_W from list to matrix
  if (is.list(samples_W)) {
    samples_W <- do.call(cbind, samples_W)
  }

  use_output_samples <- identical(type, "response") ||
    isTRUE(return_samples) ||
    prediction_needs_output_samples_for_composite_noise(
      type_names = type_names,
      estimator = estimator,
      ngme = ngme
    )

  if (use_output_samples) {
    output_samples <- build_prediction_output_samples(
      ngme = ngme,
      fm = fm,
      map = map,
      data = data,
      type = type,
      group = group,
      samples_W = samples_W
    )

    if (identical(type, "response")) {
      if (isTRUE(ngme$noise$corr_measurement)) {
        stop("type = 'response' does not currently support correlated measurement noise.")
      }

      noise_draws <- simulate_predictive_noise_draws(
        noise = ngme$noise,
        n_pred = nrow(output_samples),
        n_samp = ncol(output_samples),
        seed = seed_int + 100000L
      )
      output_samples <- ensure_prediction_sample_matrix(output_samples + noise_draws)
    }

    ret <- summarize_prediction_samples(output_samples, estimator)
    attr(ret, "samples") <- output_samples
    attr(ret, "latent_samples") <- samples_W
    return(ret)
  }

  ret <- NULL
  for (estimator in estimator) {
    # Check if estimator matches pattern like "0.025q", "0.975q", etc.
    if (grepl("^0\\.[0-9]*[1-9][0-9]*q$", estimator)) {
      # Extract the probability from the estimator string
      quantile_prob <- as.numeric(gsub("q$", "", estimator))
      # Validate the quantile probability (should be between 0 and 1, exclusive)
      if (quantile_prob <= 0 || quantile_prob >= 1) {
        stop("Quantile probability should be between 0 and 1 (exclusive). Got: ", quantile_prob)
      }
      post_W <- apply(samples_W, 1, function(x) {
        quantile(x, quantile_prob)
      })
    } else {
      post_W <- switch(estimator,
        "mean"      = apply(samples_W, 1, mean),
        "median"    = apply(samples_W, 1, median),
        "sd"        = apply(samples_W, 1, sd),
        "mode"      = apply(samples_W, 1, emprical_mode),
        stop("No such estimator available: ", estimator)
      )
    }

    # update post W (notice here W is concated)
    j <- 1
    for (i in seq_along(ngme$models)) {
      sz <- ngme$models[[i]]$W_size
      ngme$models[[i]][[estimator]] <- post_W[j:(j + sz - 1)]
      j <- j + sz
    }

    AW <- list()
    if (length(requested_model_names) > 0) {
      validate_prediction_map(map, requested_model_names)

      for (model_name in requested_model_names) {
        model <- ngme$models[[model_name]]
        loc <- map[[model_name]]
        group_use <- group
        if (inherits(model$operator$mesh, "metric_graph")) {
          loc <- as.data.frame(loc)
        } else if (model$model != "tp") {
          loc <- as.matrix(loc)
        } else {
          stopifnot(
            length(loc) == 2, # map 1 and map 2
            length_map(loc[[1]]) == length_map(loc[[2]])
          )
        }

        if (model$model %in% c("bv", "bv2", "bv_matern") && !is.null(group_use)) {
          group_use <- rep(group_use, length.out = length_map(loc))
        }

        A <- ngme_build_A(
          model$model,
          model$operator$mesh,
          loc,
          model$operator,
          group_use,
          group_levels = levels(ngme$group)
        )
        A <- expand_prediction_A_for_composite_noise(A, model)

        AW[[model_name]] <- as.numeric(A %*% model[[estimator]])
      }
    }

    preds <- 0
    for (i in seq_along(type_names)) {
      name <- type_names[[i]]
      if (name == "fe" && length(ngme$feff) > 0) {
        X_pred <- if (is.null(data) && attr(terms(fm), "intercept")) {
          matrix(1, nrow = length(AW[[1]]), ncol = 1)
        } else {
          stopifnot("Please provide covariates for predictions" = !is.null(data))
          # build plain_fm
          tf <- terms.formula(fm, specials = c("f"))
          terms <- attr(tf, "term.labels")
          intercept <- attr(tf, "intercept")
          spec_order <- attr(tf, "specials")$f - 1
          fixf <- if (length(spec_order) == 0) terms else terms[-spec_order]
          plain_fm_str <- paste("~", intercept, paste(c("", fixf), collapse = " + "))
          plain_fm <- formula(plain_fm_str)
          model.matrix(plain_fm, data = data)
        }
        if (!is.null(data) && nrow(X_pred) < nrow(data)) {
          stop("The data has NA values in the covariates.")
        }
        preds <- preds + as.numeric(X_pred %*% ngme$feff)
      } else if (name %in% names(ngme$models)) {
        preds <- preds + AW[[name]]
      }
    }
    ret[[estimator]] <- preds
  }
  attr(ret, "samples") <- samples_W
  ret
}

normalize_prediction_seed <- function(seed) {
  if (inherits(seed, "POSIXt")) {
    seed <- as.numeric(seed)
  }

  seed <- as.numeric(seed)[1]
  if (!is.finite(seed)) {
    stop("seed must be coercible to a finite numeric value.")
  }

  as.integer(abs(seed) %% 2147483647)
}


build_prediction_output_samples <- function(
    ngme,
    fm,
    map,
    data,
    type,
    group,
    samples_W) {
  n_samp <- ncol(samples_W)

  type_names <- prediction_type_names(type, ngme)
  requested_model_names <- prediction_requested_model_names(type_names, ngme)
  needs_map <- length(requested_model_names) > 0

  if (needs_map) {
    validate_prediction_map(map, requested_model_names)
  }

  fixed_effect <- NULL
  n_pred <- NULL
  if ("fe" %in% type_names && length(ngme$feff) > 0) {
    X_pred <- build_predict_model_matrix(
      fm = fm,
      ngme = ngme,
      data = data,
      n_pred = if (!is.null(data)) {
        nrow(data)
      } else {
        prediction_map_size(map, requested_model_names)
      }
    )
    fixed_effect <- as.numeric(X_pred %*% ngme$feff)
    n_pred <- length(fixed_effect)
  }

  model_draws <- list()
  if (needs_map) {
    j <- 1
    for (i in seq_along(ngme$models)) {
      model <- ngme$models[[i]]
      sz <- model$W_size
      if (model$name %in% requested_model_names) {
        loc <- map[[model$name]]
        if (inherits(model$operator$mesh, "metric_graph")) {
          loc <- as.data.frame(loc)
        } else if (model$model != "tp") {
          loc <- as.matrix(loc)
        } else {
          stopifnot(
            length(loc) == 2,
            length_map(loc[[1]]) == length_map(loc[[2]])
          )
        }

        A <- ngme_build_A(
          model$model,
          model$operator$mesh,
          loc,
          model$operator,
          group,
          group_levels = levels(ngme$group)
        )
        A <- expand_prediction_A_for_composite_noise(A, model)

        model_draws[[model$name]] <- A %*% samples_W[j:(j + sz - 1), , drop = FALSE]
        if (is.null(n_pred)) n_pred <- nrow(model_draws[[model$name]])
      }
      j <- j + sz
    }
  }

  if (is.null(n_pred)) {
    stop("Could not determine prediction size. Please provide either `map` or `data`.")
  }

  output_samples <- matrix(0, nrow = n_pred, ncol = n_samp)

  for (name in type_names) {
    if (name == "fe" && !is.null(fixed_effect)) {
      output_samples <- output_samples +
        matrix(rep(fixed_effect, n_samp), nrow = n_pred, ncol = n_samp)
    } else if (name %in% names(model_draws)) {
      output_samples <- output_samples + model_draws[[name]]
    }
  }

  ensure_prediction_sample_matrix(output_samples)
}


expand_prediction_A_for_composite_noise <- function(A, model) {
  composite_noise <- all(model$noise$noise_type %in% c(
    "normal_nig", "normal_gal", "nig_gal"
  ))

  if (composite_noise && ncol(A) * 2 == model$W_size) {
    return(cbind(A, A))
  }

  A
}


prediction_needs_output_samples_for_composite_noise <- function(
    type_names,
    estimator,
    ngme) {
  if (all(estimator == "mean")) {
    return(FALSE)
  }

  requested_model_names <- prediction_requested_model_names(type_names, ngme)
  any(vapply(requested_model_names, function(model_name) {
    model <- ngme$models[[model_name]]
    all(model$noise$noise_type %in% c("normal_nig", "normal_gal", "nig_gal"))
  }, logical(1)))
}


build_predict_model_matrix <- function(fm, ngme, data, n_pred = NULL) {
  if (is.null(data) && attr(terms(fm), "intercept")) {
    if (is.null(n_pred)) stop("Please provide covariates for predictions")
    return(matrix(1, nrow = n_pred, ncol = 1))
  }

  stopifnot("Please provide covariates for predictions" = !is.null(data))
  tf <- terms.formula(fm, specials = c("f"))
  terms <- attr(tf, "term.labels")
  intercept <- attr(tf, "intercept")
  spec_order <- attr(tf, "specials")$f - 1
  fixf <- if (length(spec_order) == 0) terms else terms[-spec_order]
  plain_fm_str <- paste("~", intercept, paste(c("", fixf), collapse = " + "))
  plain_fm <- formula(plain_fm_str)
  X_pred <- model.matrix(plain_fm, data = data)
  if (!is.null(data) && nrow(X_pred) < nrow(data)) {
    stop("The data has NA values in the covariates.")
  }
  X_pred
}


summarize_prediction_samples <- function(samples, estimator) {
  samples <- ensure_prediction_sample_matrix(samples)
  ret <- list()
  for (est in estimator) {
    if (grepl("^0\\.[0-9]*[1-9][0-9]*q$", est)) {
      quantile_prob <- as.numeric(gsub("q$", "", est))
      if (quantile_prob <= 0 || quantile_prob >= 1) {
        stop("Quantile probability should be between 0 and 1 (exclusive). Got: ", quantile_prob)
      }
      ret[[est]] <- apply(samples, 1, quantile, quantile_prob)
    } else {
      ret[[est]] <- switch(est,
        "mean" = rowMeans(samples),
        "median" = apply(samples, 1, median),
        "sd" = apply(samples, 1, sd),
        "mode" = apply(samples, 1, emprical_mode),
        stop("No such estimator available: ", est)
      )
    }
  }
  ret
}


ensure_prediction_sample_matrix <- function(samples) {
  if (is.matrix(samples)) {
    return(samples)
  }

  if (inherits(samples, "Matrix")) {
    return(as.matrix(samples))
  }

  if (is.null(dim(samples))) {
    return(matrix(samples, ncol = 1))
  }

  as.matrix(samples)
}


prediction_type_names <- function(type, ngme) {
  if (type %in% c("lp", "response")) c("fe", names(ngme$models)) else type
}


prediction_requested_model_names <- function(type_names, ngme) {
  intersect(type_names, names(ngme$models))
}


validate_prediction_map <- function(map, required_model_names) {
  if (length(required_model_names) == 0) {
    return(invisible(NULL))
  }

  stopifnot(
    "map should be a named list (name for each model)" = is.list(map) && !is.null(names(map))
  )

  missing_names <- setdiff(required_model_names, names(map))
  if (length(missing_names) > 0) {
    stop(
      "The loction for model ",
      paste(missing_names, collapse = ", "),
      " is not provided"
    )
  }

  invisible(NULL)
}


prediction_map_size <- function(map, model_names = names(map)) {
  if (is.null(map) || length(model_names) == 0) {
    return(NULL)
  }

  first_name <- model_names[[1]]
  loc <- map[[first_name]]
  if (is.null(loc)) {
    return(NULL)
  }

  if (is.list(loc) && !inherits(loc, c("data.frame", "matrix"))) {
    stopifnot(
      length(loc) == 2,
      length_map(loc[[1]]) == length_map(loc[[2]])
    )
    return(length_map(loc[[1]]))
  }

  length_map(loc)
}


simulate_predictive_noise_draws <- function(noise, n_pred, n_samp, seed) {
  if (length(noise$noise_type) == 2) {
    stop("type = 'response' does not currently support bivariate measurement noise.")
  }

  h_vec <- rep(1, n_pred)
  mu_vec <- recycle_predictive_noise_param(
    as.numeric(noise$B_mu %*% noise$theta_mu),
    n_pred = n_pred,
    param_name = "mu"
  )
  sigma_vec <- recycle_predictive_noise_param(
    as.numeric(exp(noise$B_sigma %*% noise$theta_sigma)),
    n_pred = n_pred,
    param_name = "sigma"
  )
  nu_vec <- recycle_predictive_noise_param(
    as.numeric(noise$nu_lower_bound + exp(noise$B_nu %*% noise$theta_nu)),
    n_pred = n_pred,
    param_name = "nu"
  )

  noise_draws <- vapply(seq_len(n_samp), function(i) {
    simulate_noise(
      noise_type = noise$noise_type,
      h_vec = h_vec,
      mu_vec = mu_vec,
      sigma_vec = sigma_vec,
      nu_vec = nu_vec,
      seed = seed + i,
      single_V = noise$single_V
    )
  }, numeric(n_pred))

  matrix(noise_draws, nrow = n_pred, ncol = n_samp)
}


recycle_predictive_noise_param <- function(param_vec, n_pred, param_name) {
  if (length(param_vec) == n_pred) {
    return(param_vec)
  }

  if (length(param_vec) == 1) {
    return(rep(param_vec, n_pred))
  }

  if (all(isTRUE(all.equal(param_vec, rep(param_vec[[1]], length(param_vec)))))) {
    return(rep(param_vec[[1]], n_pred))
  }

  stop(
    "type = 'response' requires measurement-noise ", param_name,
    " to be scalar or already match prediction length. Got length ",
    length(param_vec), " for ", n_pred, " prediction locations."
  )
}
