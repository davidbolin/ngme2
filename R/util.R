# toc
# ngme.ts.make.A
# ngme.parse.formula

#' Make observation matrix for time series
#'
#' @param loc   integers (after sorting, no gaps > 1)
#' @param replicates indicating replicate measure at same location
#' @param range range for the mesh
#'  by default range=(min(loc), max(loc))
#'
#' @return A matrix (length(loc) * length(unique(loc)))
#' @export
#'
#' @examples ngme.ts.make.A(c(1, 2, 2), end = 5, replicates = c(1, 1, 2))
#'
ngme.ts.make.A <- function(loc, replicates = NULL, range = c(1, max(loc))) {
  if (is.null(loc)) return (NULL)

  n_loc <- length(loc)
  nrep <- 1

  start <- range[1]; end <- range[2]
  n_range <- end - start + 1

  if (is.null(replicates)) {
    replicates <- rep(1, n_range)
  }

  unique_rep <- unique(replicates)
  nrep <- length(unique_rep)

  A <- matrix(0, nrow = n_loc, ncol = n_range * nrep)

  for (i in 1:n_loc) {
    ncol_rep <- which(unique_rep == replicates[i])
    A[i, (ncol_rep - 1) * n_range + loc[i] - start + 1] <- 1
  }
  as(A, "dgCMatrix")
}


# ngme.ts.make.A <- function(loc) {
#   n_loc = length(loc)
#   n_range = max(loc)-min(loc)+1
#   if (any((diff(sort(loc))) > 1)) stop("no gaps allowed")
#
#   A = matrix(0, nrow=n_loc, ncol=n_range)
#   for (i in 1:n_loc) {
#     A[i, loc[i]-min(loc)+1] = 1
#   }
# #Do an na.rm
#   as(A, "dgCMatrix")
# }


#' Convert sparse matrix into sparse dgCMatrix
#'
#' @return sparse dgCMatrix
#' @export
#'
#' @examples
ngme.as.sparse <- function(G) {
  tryCatch(
    expr={
      G <- as(G, "dgCMatrix")
    },
    error = function(e) {
      G <- as(G, "dgTMatrix")
      idx <- which(G@i <= G@j)
      G = Matrix::sparseMatrix(i=G@i[idx], j = G@j[idx], x= G@x[idx],
                               symmetric=FALSE, index1 = FALSE)
      G = as(G, "dgCMatrix")
    },
    finally = {
      G
    }
  )
}

#' Make index for the matern model
#'
#' @param name
#' @param n.spde
#' @param n.repl
#' @param mesh
#' @param dim
#'
#' @return
#' @export
#'
#' @examples
ngme.matern.make.index <- function(n.spde=NULL,
                                   n.repl = 1,
                                   mesh = NULL,
                                   dim = NULL){

  if(is.null(n.spde)&&is.null(mesh)){
    stop("You should provide either n.spde or mesh!")
  }

  if(!is.null(mesh)){
    n_mesh = mesh$n

    if(mesh$manifold == "R1"){
      dim = 1
    } else if(mesh$manifold == "R2"){
      dim = 2
    } else{
      stop("The domain must be flat manifolds of dimension 1 or 2, that is,
         the domain must be a line or a plane.")
    }

  } else{
    n_mesh <- n.spde
    if(is.null(dim)) {
      stop("You should provide the dimension d!")
    }
  }

  out <- list()
  out$index <- rep(1:n_mesh, times = n.repl)

  out$replicates <- rep(1:n.repl, each = n_mesh)
  return(out)
}


#' Split the data according to NA in response variable
#'
#' @param formula
#' @param data
#'
#' @return
#' @export
#'
#' @examples
parse_formula_NA <- function(formula, data) {
  stopifnot("Please provide a valid formula" = inherits(formula, "formula"))

  Y_full <- model.frame(formula, data, na.action = NULL)[[1]]

  mf <- model.frame(formula, data)
  Y_data <- mf[[1]]
  index_data  <- as.numeric(rownames(mf))
  X_data      <- model.matrix((terms(formula)), data)

  # watch out! X_full can be 0*0
  X_full      <- model.matrix(delete.response(terms(formula)), data)
  contain_NA  <- !is.null(attr(mf, "na.action"))

  # if there is NA in response variable
  if (!is.null(attr(mf, "na.action"))) {
    index_NA     <- as.numeric(attr(mf, "na.action"))

    # X_NA         <- X_full[index_NA, , drop = FALSE]
    if (all(dim(X_full) == c(0, 0)))
      X_NA      <- X_full
    else
      X_NA      <- X_full[index_NA, , drop = FALSE]
  } else {
    index_NA <- X_NA <- NULL
  }

  list(
    length      = length(Y_data),
    Y_data      = Y_data,
    X_data      = X_data,
    index_data  = index_data,
    contain_NA  = contain_NA,
    index_NA    = index_NA,
    X_NA        = X_NA

    # Y_full      = Y_full # only for f to use
  )
}


# post.sampleW(ngme) {
#   W
# }


#' @name get_inla_mesh_dimension
#' @title Get the dimension of an INLA mesh
#' @description Get the dimension of an INLA mesh
#' @param inla_mesh An INLA mesh
#' @return The dimension of an INLA mesh.
#' @noRd
#'
get_inla_mesh_dimension <- function(inla_mesh) {
  cond1 <- inherits(inla_mesh, "inla.mesh.1d")
  cond2 <- inherits(inla_mesh, "inla.mesh")
  stopifnot(cond1 || cond2)
  if (inla_mesh$manifold == "R1") {
    d <- 1
  } else if (inla_mesh$manifold == "R2") {
    d <- 2
  } else {
    stop("The mesh should be from a flat manifold.")
  }
  return(d)
}

#' Parse the formula for ngme function
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
ngme.parse.formula <- function(
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

ngme.format <- function(x) {
  format(x, digits = 3)
}