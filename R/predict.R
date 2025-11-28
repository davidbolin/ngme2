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
#' @param group which filed to predict
#'   (used for bivariate model, should be of same length as map)
#' @param estimator what type of estimator. Options include:
#'   - "mean", "median", "mode", "sd": standard estimators
#'   - "0.XXXq": any quantile specified as probability (e.g., "0.025q", "0.5q", "0.975q")
#' @param sampling_size size of posterior sampling
#' @param burnin_size size of posterior burnin
#' @param seed random seed
#' @param train_idx optional vector of training indices to use for posterior sampling.
#'   If provided, only these indices from the original data will be used for training,
#'   similar to cross-validation. If NULL, uses all original training data.
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
    seed = Sys.time(),
    train_idx = NULL,
    ...) {
  fm <- attr(object, "fit")$formula
  ngme <- object$replicate[[1]]
  stopifnot(
    sampling_size > 0,
    "Make sure the object is of class 'ngme'." = inherits(object, "ngme")
  )

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
    seed = seed
  )[["W"]]

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
      post_W <- apply(as.data.frame(samples_W), 1, function(x) {
        quantile(x, quantile_prob)
      })
    } else {
      post_W <- switch(estimator,
        "mean"      = mean_list(samples_W),
        "median"    = apply(as.data.frame(samples_W), 1, median),
        "sd"        = apply(as.data.frame(samples_W), 1, sd),
        "mode"      = apply(as.data.frame(samples_W), 1, emprical_mode),
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

    if (!is.null(map)) {
      stopifnot(
        "map should be a named list (name for each model)" = is.list(map) && !is.null(names(map))
      )
      names <- names(map)
      stopifnot(length(names) == length(ngme$models))

      AW <- list()
      for (i in seq_along(ngme$models)) {
        loc <- map[[ngme$models[[i]]$name]]
        if (is.null(loc)) stop("The loction for model ", ngme$models[[i]]$name, " is not provided")
        if (inherits(ngme$models[[i]]$operator$mesh, "metric_graph")) {
          loc <- as.data.frame(loc)
        } else if (ngme$models[[i]]$model != "tp") {
          loc <- as.matrix(loc)
        } else {
          stopifnot(
            length(loc) == 2, # map 1 and map 2
            length_map(loc[[1]]) == length_map(loc[[2]])
          )
        }

        mesh <- ngme$models[[i]]$operator$mesh
        W <- ngme$models[[i]][[estimator]]

        A <- ngme_build_A(
          ngme$models[[i]]$model,
          mesh,
          loc,
          ngme$models[[i]]$operator,
          group,
          group_levels = levels(ngme$group)
        )

        AW[[ngme$models[[i]]$name]] <- as.numeric(A %*% W)
      }
    }
    # lp case is just fe + A1 * W1 + A2 * W2
    # e.g. names <- c("fe", "field1", "field2")
    type_names <- if (type == "lp") c("fe", names(ngme$models)) else type

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
