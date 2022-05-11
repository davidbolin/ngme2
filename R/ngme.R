#' Fit a non-guassian model
#'
#' @param formula formula
#' @param family  distribution family
#' @param data    a dataframe contains data
#' @param init    initial values 
#' @param controls control variables
#'
#' @return a list of outputs
#' @export
#'
#' @examples
#' ngme(formula = (y1 | y2 ~ x1 + x2 + f(x1, model="SPDE", var="nig") + f(W) | x3 + f(X|I, model="ar1", var="nig")),
#'
#' ngme(y ~ f(x, model="AR1", var="NIG"), family="normal", data=df)


ngme <- function(formula,
                 data,
                 family  = "normal",
                 init    = list(),
                 controls = list(burnin           = 100,
                                iterations        = 100,
                                gibbs_sample      = 5,
                                stepsize          = 0.5,

                                opt_fix_effect    = TRUE,
                                fix_trueVW        = FALSE,
                                trueSV            = NULL,
                                trueW             = NULL),
                  debug = FALSE)
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

  # 2. parse the formula
  fm = Formula::Formula(formula)

  if (all(length(fm)==c(2,2))) {
    # bivariate model
    ########## todo
    lfm = formula(fm, lhs=1, rhs=1)
    rfm = formula(fm, lhs=2, rhs=2)

    # todo
  }
  else if (all(length(fm)==c(1,1))) {
    ####### univariate case  
    fm = formula(fm)
    
    # 1. extract f and eval  2. get the formula without f function
    res = ngme.interpret.formula(fm, data)
    latents_in = res$latents_in
    plain.fm = res$plain.fm

    # eval part without f (fm = y ~ x1 + x2)
    Y = model.frame(plain.fm, data)[[1]]
    X = model.matrix(plain.fm, data) # design matrix
    n =length(Y)

    ############### n_regs is the dim of the block matrix
    n_regs = sum(unlist(lapply(latents_in, function(x) x["n_reg"] )))

    # 3. prepare in_list for estimate
    general_in <- list( n                = n,
                      Y                = Y,
                      X                = X,
                      family           = "normal",
                      n_regs           = n,
                      init             = list(beta=lm.fit(X, Y)$coeff,
                                              sigma_eps = 1),
                      debug=debug)

    in_list = list(general_in = general_in,
                  latents_in = latents_in,
                  config_in  = controls)

  } else {
    stop("unknown structure of formula")
  }

  # debug
  if (debug) print(in_list)

  # estimate
  out = estimate_cpp(in_list)
  
  # construct output
  out$n_fe     = ncol(X)
  out$n_latent = length(latents_in)
  class(out)   = "ngme"

  print(paste("total time is", Sys.time()-time.start))
  return (out)
}

