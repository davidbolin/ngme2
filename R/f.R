#' f function
#'
#' @param x numerical vector, which is the covariate
#' @param model specify a model
#' @param var variance component type
#' @param init initial value
#' @param config control variables
#' @param debug debug
#'
#' @return a list of objects
#' @export
f <- function(x = NULL,
              model  = "ar1",
              var    = "NIG",
              init   = list(kappa     = 0.5,
                            mu        = 0,
                            sigma     = 1,
                            nu        = 0.5),
              control = control.f(),
              debug  = FALSE
              ) {
  ## in ... is the name of the covariate  and possibly the location of the weights
  ## like f(x, weights)

  n = length(x)

  # construct operator
  if (model=="ar1") {
    # G
      G <- Matrix::Matrix(diag(n));
      G <- as(G, "dgCMatrix");

    # C
      C <- Matrix::Matrix(0, n, n)
      C[seq(2, n*n, by=n+1)] <- -1
      C <- as(C, "dgCMatrix");

    # A
      A <- as(Matrix::Matrix(diag(n)), "dgCMatrix");

    # h
      h <- rep(1.0, n)

    operator_in   = list(C=C, G=G, numerical_dK=FALSE)
  }

  # construct variance component
  if (var=="nig" || var=="NIG") {
    var_in = list(type = "ind_IG")
  }

  if (var=="normal") {
    var_in = list(type = "normal")
  }

  # construct latent_in
  la_in <- list(type          = model,
                var.type      = var,
                n_reg         = n,        # !: make sure this is the second place
                A             = A,
                h             = h,
                opt_kappa     = control$opt_kappa,
                opt_mu        = control$opt_mu,
                opt_sigma     = control$opt_sigma,
                opt_var       = control$opt_var,
                numer_grad    = control$numer_grad,
                use_precond   = control$use_precond,
                eps           = control$eps,
                debug         = debug,
                operator_in   = operator_in,
                var_in        = var_in,
                init_value    = init)

  return (la_in)
}

  # vars <- as.list(substitute(list(...)))[-1]
  # d <- length(vars)
  # if (d == 0L) {
  #   stop(paste("Missing covariate in f(...) for model=", model))
  # }
  # term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  # if (debug) {
  #   print(vars)
  # }

  # ## get the weights
  # term <- attr(terms(reformulate(term)), "term.labels")
  # if (d > 1) {
  #   weigths <- attr(terms(reformulate(weights)), "weights.labels")
  # }

  # # use model.frame to extract the covariates

  # x = model.frame(paste0("~", term), data=data)[[term]]

  # ## get the model
  # if (model=="ar1") {
  #   latent_in <- ngme.model.ar1(x, var=var, init=init, config=config)
  # }

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
