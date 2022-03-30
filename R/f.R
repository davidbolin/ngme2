#' f function
#'
#' @param ... a symbol like x
#' @param data  provide data.frame or
#' @param model specify a model
#' @param noise noise type
#' @param response response type
#'
#'
#' @return a list of objects
#' @export
f <- function(...,
              data   = NULL,
              model  = NULL,
              debug  = FALSE,
              var    = "NIG",
              init   = list(kappa     = 0.5,
                            mu        = 0,
                            sigma     = 1,
                            nu        = 0.5,),
              config = list(opt_kappa     = TRUE,
                            opt_mu        = TRUE,
                            opt_sigma     = TRUE,
                            opt_var       = TRUE,)
              ) {
  ## in ... is the name of the covariate  and possibly the location of the weights
  ## like f(covariate, weights)

  vars <- as.list(substitute(list(...)))[-1]
  d <- length(vars)
  if (d == 0L) {
    stop(paste("Missing covariate in f(...) for model=", model))
  }
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  if (debug) {
    print(vars)
  }

  ## the second term in ... is the (possible) weights for the selected covariate!
  if (d == 1) {
    weights <- NULL
  } else if (d == 2) {
    weights <- deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)
  } else if (d > 2) {
    stop(paste("To many variables included in f():", inla.paste(vars)))
  } else if (d == 0) {
    stop("At least one variable in f() needs to be defined")
  }

  ## get the weights
  term <- attr(terms(reformulate(term)), "term.labels")
  if (d > 1) {
    weigths <- attr(terms(reformulate(weights)), "weights.labels")
  }

  # use model.frame to extract the covariates

  x = model.frame(paste0("~", term), data=data)[[term]]

  ## get the model
  if (model=="ar1") {
    latent_in <- ngme.model.ar1(x, var=var, init=init, config=config)
  }

  return (latent_in)

  # expr_formula <- rlang::enexpr(formula)
  # index <- NULL
  # formula_string <- NULL
  #
  # if (class(expr_formula)=="call" && expr_formula[[1]]=="|") {  # 2 arguments
  #   first  <- eval(expr_formula[[2]], envir = data)
  #   index <- eval(expr_formula[[3]], envir = data)
  #   formula_string <- "response ~ first | index"
  #
  # } else { # 1 argument
  #   first  <- eval(expr_formula, envir = data)
  #   formula_string <- "response ~ first"
  # }
  #


  # # random effect case
  #   if (is.null(model) && noise=="gaussian") {
  #     n <- length(first)
  #
  #     rtrn <- list()
  #     lf <- lme4::lFormula(formula_string)
  #     idx <- levels(lf$reTrms$flist$index)
  #     diags <- length(idx)
  #     nr <- dim(lf$reTrms$Zt)[1] / diags
  #     nc <- dim(lf$reTrms$Zt)[2] / diags
  #
  #     for (i in 1:length(idx)) {
  #       # K matrix
  #       K <- function(params) {
  #         n <- nc
  #         tmp <- matrix(0, nrow=n, ncol=n)
  #         diag(tmp) <- exp(params[1:n])
  #         count <- n; j <- n-1
  #         for (i in 1:(n-2)) {
  #           diag(tmp[(1+i):n, 1:(n-i)]) <- params[(count+1):(count+j)]
  #           count <- count+j; j <- j-1
  #         }
  #         tmp[n, 1] <- params[length(params)]
  #         tmp
  #       }
  #
  #       # mean
  #       m <- function(a) {
  #         rep(0, nc)
  #       }
  #
  #       # A matrix
  #       tmp <- lf$reTrms$Zt[((i-1)*nr+1):(i*nr), ((i-1)*nc+1):(i*nc)]
  #       A <- t(as.matrix(tmp))
  #
  #       dK <- NULL
  #       R <- NULL
  #
  #       tmp <- list(A, m, K, dK, R, "gaussian")
  #       names(tmp) <- c("A", "m", "K", "dK", "R", "noise")
  #       rtrn[[i]] <- tmp
  #     }
  #   }
  #
  #     class(rtrn) <- "f"
  #     return (rtrn)
  # }
}


