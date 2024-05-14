# This file contains function related to model predict function

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
#' @param estimator what type of estimator, c("mean", "median", "mode", "quantile")
#' @param sampling_size size of posterior sampling
#' @param burnin_size size of posterior burnin
#' @param seed random seed
#' @param q quantile if using "quantile"
#' @param ... extra argument from 0 to 1 if using "quantile"
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
predict.ngme <- function(
  object,
  map,
  data = NULL,
  type = "lp",
  estimator = c("mean", "sd", "5quantile", "95quantile", "median", "mode"),
  sampling_size = 100,
  burnin_size = 100,
  q = NULL,
  seed = Sys.time(),
  ...
) {
  fm <- attr(object, "fit")$formula
  ngme <- object$replicate[[1]]
  stopifnot(
    sampling_size > 0,
    "Make sure the object is of class 'ngme'." = inherits(object, "ngme")
  )
  samples_W <- sampling_cpp(ngme, 
    n=sampling_size, 
    n_burnin = sampling_size, 
    posterior=TRUE, 
    seed=seed
  )[["W"]]

  ret <- NULL
  for (estimator in estimator) {
    post_W <- switch(estimator,
      "mean"      = mean_list(samples_W),
      "median"    = apply(as.data.frame(samples_W), 1, median),
      "sd"        = apply(as.data.frame(samples_W), 1, sd),
      "mode"      = apply(as.data.frame(samples_W), 1, emprical_mode),
      "5quantile" = apply(as.data.frame(samples_W), 1, function(x) {quantile(x, 0.05)}),
      "95quantile" = apply(as.data.frame(samples_W), 1, function(x) {quantile(x, 0.95)}),
      "quantile"  = {
        stopifnot("please provide quantile argument q between 0 to 1"
          = !is.null(q) && length(q) == 1 && q > 0 && q < 1)
        apply(as.data.frame(samples_W), 1, function(x) {quantile(x, q)})
      },
      stop("No such estimator available")
    )

    # update post W (notice here W is concated)
    j <- 1
    for (i in seq_along(ngme$models)) {
      sz <- ngme$models[[i]]$W_size
      ngme$models[[i]][[estimator]] <- post_W[j:(j + sz - 1)]
      j <- j + sz
    }

    if (!is.null(map)) {
      stopifnot(
        "map should be a named list (name for each model)"
          = is.list(map) && !is.null(names(map))
      )
      names <- names(map)
      stopifnot(length(names) == length(ngme$models))

      AW <- list()
      for (i in seq_along(ngme$models)) {
        loc <- map[[ngme$models[[i]]$name]]
    if (is.null(loc)) stop("The loction for model ", ngme$models[[i]]$name, " is not provided")
        if (ngme$models[[i]]$model != "tp")
          loc <- as.matrix(loc)
        else {
          stopifnot(
            length(loc) == 2, # map 1 and map 2
            length_map(loc[[1]]) == length_map(loc[[2]])
          )
        }

        AW[[ngme$models[[i]]$name]] <- with(ngme$models[[i]], {
          mesh <- operator$mesh
          W <- ngme$models[[i]][[estimator]]
          A <- if (inherits(mesh, "metric_graph"))
              mesh$fem_basis(loc)
            else if (model=="tp") {
              A1 <- fmesher::fm_basis(mesh[[1]], loc = loc[[1]])
              A2 <- fmesher::fm_basis(mesh[[2]], loc = as.matrix(loc[[2]]))
              fmesher::fm_row_kron(A1, A2)
            }
            else fmesher::fm_basis(mesh, loc = loc)

          if (model == "bv") A <- Matrix::bdiag(A, A)
          stopifnot(ncol(A) == length(W))
          as.numeric(A %*% W)
        })
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
          preds <- preds + as.numeric(X_pred %*% ngme$feff)
      } else if (name %in% names(ngme$models)) {
        preds <- preds + AW[[name]]
      }
    }
    ret[[estimator]] <- preds
  }
  ret
}
