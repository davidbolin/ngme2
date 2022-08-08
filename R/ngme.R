#' Fit a non-guassian model
#'
#' @param formula formula
#' @param family  distribution family
#' @param data    a dataframe contains data
#' @param controls control variables
#' @param debug  debug option
#' @param start  1. last fitting object 2. ngme.start() for block model
#'
#' @return a list of outputs
#' @export
#'
#' @examples
#' ngme(formula = (y1 | y2 ~ x1 + x2 + f(x1, model="SPDE", var="nig") + f(W) | x3 + f(X|I, model="ar1", var="nig")),
#'
#' ngme(y ~ f(x, model="AR1", var="NIG"), family="normal", data=df)
#'
ngme <- function(formula,
                 data,
                 family   = "normal",
                 controls = ngme.control(),
                 debug    = ngme.debug(),
                 start    = ngme.start(),
                 notrun = FALSE)
{
  time.start <- Sys.time()

  # 1. check up
  if (is.null(formula)) {
    stop("Usage: ngme(formula, family, data, ...); see ?ngme\n")
  }

  if (is.null(data)) {
    stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")
  }

  if (!is.data.frame(data) && !is.list(data)) {
    stop("\n\tArgument `data' must be a data.frame or a list.")
  }

  # # store data into temp. file
  # path <- file.path(tempdir(), "ngme_df.rda")
  # save(data, file = path)

  # generate debug option
  if (!is.null(debug$trueW)) {
    debug$fixW = TRUE
  }

  # 2. parse the formula
  fm = Formula::Formula(formula)
# Y1|Y2|Y3 ~ ..|..|..
  if (all(length(fm)==c(2,2))) { ######################### bivariate model
    lfm = formula(fm, lhs=1, rhs=1)
    rfm = formula(fm, lhs=2, rhs=2)
    ########## to-do
  }
  else if (all(length(fm)==c(1,1))) {  ########################## univariate case
    fm = formula(fm)

    # 1. extract f and eval  2. get the formula without f function
    res = ngme.interpret.formula(fm, data)
    latents_in = res$latents_in
    plain.fm = res$plain.fm

    # eval part without f (fm = y ~ x1 + x2)
    Y = model.frame(plain.fm, data)[[1]]
    X = model.matrix(plain.fm, data) # design matrix

    ############### n_meshs is the dim of the block matrix
    n_meshs     = sum(unlist(lapply(latents_in, function(x) x["n_mesh"] )))
    n_la_params = sum(unlist(lapply(latents_in, function(x) x["n_la_params"] )))
    model.types = unlist(lapply(latents_in, function(x) x["model_type"] ))
    var.types   = unlist(lapply(latents_in, function(x) x["var.type"] ))

    n_feff = ncol(X);
    n_merr = 1
    n_params = n_la_params + n_feff + n_merr

    # 3. prepare in_list for estimate
    lm.model = lm.fit(X, Y)
    general_in <- list(
      Y                = Y,
      X                = X,
      family           = "normal",
      n_meshs          = n_meshs,
      n_la_params      = n_la_params,
      n_params         = n_params # how many param to opt. in total
    )

####### 4. Set starting point / initial values (beta, sigma_eps, latents)
    # use prevous ngme object
    if (inherits(start, "ngme")) {
      start <- start$output

      # put the last estimates into latents_in
      for (i in seq_along(latents_in)) {
        estimates = start$latent.model[[i]]$estimates

        # If not specified, use previous fitting
        # if (is.null(latents_in[[i]]$start$theta_K))
          latents_in[[i]]$start$theta_K = estimates[[1]]
        # if (is.null(latents_in[[i]]$start$theta_mu))
          latents_in[[i]]$start$theta_mu = estimates[["theta.mu"]]
        # if (is.null(latents_in[[i]]$start$theta_sigma))
          latents_in[[i]]$start$theta_sigma = estimates[["theta.sigma"]]
        # if (is.null(latents_in[[i]]$start$theta_noise))
          latents_in[[i]]$start$theta_noise = estimates[["theta.noise"]]

        latents_in[[i]]$start$V = start$latent.model[[i]][["V"]]
      }
    }

    if (is.null(start$block.W) && isTRUE(controls$fixW)) stop("if fixing W, the initial W should be provided.")
    if (is.null(start$fixed.effects)) start$fixed.effects = lm.model$coeff
    if (is.null(start$mesurement.noise)) start$mesurement.noise = sd(lm.model$residuals)

    in_list = list(general_in = general_in,
                  latents_in  = latents_in,
                  start_in    = start,
                  control_in  = controls,
                  debug       = debug)

  } else {
    stop("unknown structure of formula")
  }

# print
if (debug$debug) print(str(in_list))

  ################# Run CPP ####################

  if (notrun) {
    print(str(in_list))
    stop()
  }
  out = estimate_cpp(in_list)

  # construct output
    out$input = in_list

    # fix_eff
    out$n_fe     = ncol(X)

    # m_err
    out$family = general_in$family

    # operator
    out$n_la_params = unlist(lapply(latents_in, function(x) x["n_la_params"] ))

    # process
    out$n_latent = length(latents_in)
    out$model.types = model.types
    out$var.types   = var.types

    # generate output
    out$result <- out$output
    out$result$block.W <- NULL
    for (i in seq_along(out$result$latent.model)) {
      out$result$latent.model[[i]]$W <- NULL
      out$result$latent.model[[i]]$V <- NULL
    }

  class(out)   = "ngme"

  print(paste("total time is", Sys.time()-time.start))
  return (out)
}

