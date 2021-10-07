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
#'
#'

f <- function(formula, data, model=NULL, noise="gaussian", response=NULL) {

  expr_formula <- rlang::enexpr(formula)
  index <- NULL
  formula_string <- NULL

  if (class(expr_formula)=="call" && expr_formula[[1]]=="|") {  # 2 arguments
    first  <- eval(expr_formula[[2]], envir = data)
    index <- eval(expr_formula[[3]], envir = data)
    formula_string <- "response ~ first | index"

  } else { # 1 argument
    first  <- eval(expr_formula, envir = data)
    formula_string <- "response ~ first"
  }

  # one argument case
  if (is.null(index)) {
    # AR-1
    if (model == "AR(1)") {
      n <- length(first)
      # construct A
      A <- diag(n)

      # return a function of alpha
      # construct K(a)
      K <- function(a){
        K_temp <- diag(n)
        K_temp[seq(2, n*n, by=n+1)] <- -a
        K_temp
      }

      # construct R - placeholder for rational approximation
      R <- NULL

      # m(theta, v)
      m <- function(a) {
        rep(0, n)
      }

      # construct the derivative of K
      dK <- function(a) {
        K_temp <- matrix(data=0, nrow=n, ncol=n)
        K_temp[seq(2, n*n, by=n+1)] <- -1
        K_temp
      }

      tmp <- list(A, m, K, dK, R, noise)
      names(tmp) <- c("A", "m", "K", "dK", "R", "noise")
      rtrn <- list(tmp)
      class(rtrn) <- "f"
      return (rtrn)
    } else {
      stop("Model support yet")
    }
  } else {
  # random effect case
    if (is.null(model) && noise=="gaussian") {
      n <- length(first)

      rtrn <- list()
      lf <- lme4::lFormula(formula_string)
      idx <- levels(lf$reTrms$flist$index)
      diags <- length(idx)
      nr <- dim(lf$reTrms$Zt)[1] / diags
      nc <- dim(lf$reTrms$Zt)[2] / diags

      for (i in 1:length(idx)) {
        # K matrix
        K <- function(params) {
          n <- nc
          tmp <- matrix(0, nrow=n, ncol=n)
          diag(tmp) <- exp(params[1:n])
          count <- n; j <- n-1
          for (i in 1:(n-2)) {
            diag(tmp[(1+i):n, 1:(n-i)]) <- params[(count+1):(count+j)]
            count <- count+j; j <- j-1
          }
          tmp[n, 1] <- params[length(params)]
          tmp
        }

        # mean
        m <- function(a) {
          rep(0, nc)
        }

        # A matrix
        tmp <- lf$reTrms$Zt[((i-1)*nr+1):(i*nr), ((i-1)*nc+1):(i*nc)]
        A <- t(as.matrix(tmp))

        dK <- NULL
        R <- NULL

        tmp <- list(A, m, K, dK, R, "gaussian")
        names(tmp) <- c("A", "m", "K", "dK", "R", "noise")
        rtrn[[i]] <- tmp
      }
    }

      class(rtrn) <- "f"
      return (rtrn)
  }
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
