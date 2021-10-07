#' Fit a non-guassian model
#'
#' @param formula some formula
#' @param data a list with data input
#' @param ... additional arguments
#'
#' @return a list of outputs
#' @export
#'
#' @examples
#' ngme(formula = Formula::Formula(Y | Z ~ X1 + X2 + f(S, model="SPDE", type="..") + f(W) | X3 + X4 + f(X|I, model="SPDE")),
#' data=list(Z=1:3, Y=2:4, X1=3:5, X2=4:6, X3=5:7, X4=6:8))
#'

# X1 <- 1:3; X2 <- 2:4; X3 <- 7:9; X4 <- 5:7;  Z <- 3:5; Y <- 10:12
# ff <- Y | Z ~ X1 + X2 + f(S, model="SPDE", type="..") + f(W) | X3 + X4 + f(X|I, model="SPDE")
# ff <- Formula::Formula(ff)
# ngme(ff)

ngme <- function(formula, data, ...) {
  ngme_call <- match.call()
  MF <- match.call(expand.dots = FALSE)

  # if data is not provided, verify the current R workspace
  if (missing(data)) {
    data <- environment(formula)
  }

  formula <- Formula::Formula(formula)

  rtrn <- list()
  rtrn$call <- ngme_call
  rtrn$first <- formula(formula, lhs=1, rhs=1)
  rtrn$second <- formula(formula, lhs=2, rhs=2)

  # index for functional and non-functional terms
  first_drop <- grep("^f[(]", attr(terms(formula(rtrn$first)), "term.labels"))
  second_drop <- grep("^f[(]", attr(terms(formula(rtrn$second)), "term.labels"))

  first_keep <- 1:length(terms(rtrn$first))
  second_keep <- 1:length(terms(rtrn$second))

  first_keep <- first_keep[! first_keep %in% first_drop]
  second_keep <- second_keep[! second_keep %in% second_drop]

  # extract non-functional terms
  if (length(first_drop)) {
    rtrn$first_v <- formula(drop.terms(terms(rtrn$first), first_drop, keep.response = T))
  } else {
    rtrn$first_v <- rtrn$first
  }
  if (length(second_drop)) {
    rtrn$second_v <- formula(drop.terms(terms(rtrn$second), second_drop, keep.response = T))
  } else {
    rtrn$second_v <- rtrn$second
  }

  # extract functional terms
  if (length(first_keep)) {
    rtrn$first_f <- formula(drop.terms(terms(rtrn$first), first_keep, keep.response = T))
  } else {
    rtrn$first_f <- rtrn$first
  }
  if (length(second_keep)) {
    rtrn$second_f <- formula(drop.terms(terms(rtrn$second), second_keep, keep.response = T))
  } else {
    rtrn$second_f <- rtrn$second
  }
  rtrn$first_f <- attr(terms(rtrn$first_f), "term.labels")
  rtrn$second_f <- attr(terms(rtrn$second_f), "term.labels")

  # split formula into 2 parts
  formula <- Formula::Formula(formula)
  terms_formula <- terms(formula)
  term_labels <- attr(terms_formula, "term.labels")

  # model matrix
  rtrn$first_mf <- model.frame(rtrn$first_v, data)
  rtrn$second_mf <- model.frame(rtrn$second_v, data)

  class(rtrn) <- "ngme"
  return (rtrn)
}




