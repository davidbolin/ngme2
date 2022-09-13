#' Interpret the formula for ngme function
#'
#' @param gf formula
#' @param data data.frame
#' @param debug
#'
#' @return
#'  1. plain formula without f function
#'  2. latents_in - from each f function
#' @export
#'
#' @examples
ngme.interpret.formula <- function(
  gf,
  data,
  debug=FALSE
) {
  # eval the response variable to see NA
  # Y <- eval(gf[[2]], envir = data)
  # index_prd <- which(is.na(Y))
  # index_est <- which(!is.na(Y))

  # adding special mark
  tf <- terms.formula(gf, specials = c("f"))

  terms <- attr(tf, "term.labels")
  intercept <- attr(tf, "intercept")

  latents_in <- list()
  # order of f terms in labels
  spec_order <- attr(tf, "specials")$f - 1
  for (i in spec_order) {
    if (!grepl(re2 <- "data *=", terms[i])) {
      # adding data=data if not specified
      str <- gsub("^f\\(", "ngme2::f(data=data,", terms[i])
    } else if (grepl("data *= *NULL", terms[i])) {
      # change data=NULL to data=data
      str <- gsub("^f\\(", "ngme2::f(", terms[i])
      str <- gsub("data *= *NULL", "data=data", str)
    } else {
      # keep data=sth.
      str <- gsub("^f\\(", "ngme2::f(", terms[i])
    }

    # adding 1 term for furthur use in f
    # data$ngme_response <- Y
    res <- eval(parse(text = str), envir = data)
    latents_in[[length(latents_in) + 1]] <- res
  }
  fixf <- terms[-spec_order]

  # construct plain formula without f
  fm <- as.character(attr(tf, "variables")[[2]])
  fm <- paste(fm, "~", intercept, paste(c("", fixf), collapse = " + "))

  list(
    latents_in = latents_in,
    plain.fm = formula(fm)
    # ,
    # index_prd = index_prd,
    # index_est = index_est
  )
}