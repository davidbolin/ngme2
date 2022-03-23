#' f function
#'
#' @param formula provide a formula
#' @param data  provide data
#' @param model specify a model
#' @param noise noise type
#' @param response response type
#'
#'
#' @return a list of objects
#' @export

f <- function(...,
              model=NULL,
              noise="gaussian",
              init=NULL,
              debug=FALSE,
              config=NULL)
  {
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

  ## check the model
  if (model=="AR1") {
    ar_in <- ngme.model.ar1(term)
  }

  ar_in


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


#
# library(lme4)
# str(sleepstudy)
# lmod <- lFormula(Reaction ~ (0+Days|Subject), sleepstudy)
# lmod$X
# lmod$reTrms
#
# library(rlang)
# # test
# X <- c(1,2,3,4)
# f_rtrn <- f(X, list(X=X), model="AR(1)")
# f_rtrn[[1]]$K(0.5)
# f_rtrn[[1]]$m(0.5)
#
# class(expr(X|I))
#
# I <- 1:5
# f_rtrn <- f(X|I, list(X=X, I=I), model="AR(1)")
# f_rtrn
#
#
#
# X <- diag(x=1, nrow=5, ncol=5)
# X
# X[seq(2, 5*5, by=6)] <- 0.5
# X
#
# expression(X|I)
#
# class(expr(X))
# class(expr(X|I))
# class(expr(X|I)[[2]])
# eval(expr(X), envir = data.frame(X=1))
#
#
#
#
# for (i in 1:5) {
#   l[[length(l) + 1]] = 1
# }
# l
#
#
# lmod <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy)
# names(lmod)
# lmod$reTrms
