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
#' @examples
#' ngme_ts_make_A(c(1, 2, 2), replicates = c(1, 1, 2))
#' ngme_ts_make_A(c(1, 2, 2), range = c(1, 5))
ngme_ts_make_A <- function(
  loc,
  replicates = NULL,
  range = c(1, max(loc))
) {
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
ngme_as_sparse <- function(G) {
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
ngme_parse_formula <- function(
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
    res <- eval(parse(text = str), envir = data, enclos = parent.frame())
    latents_in[[length(latents_in) + 1]] <- res
  }
  # watch out! terms[-double(0)] -> character(0)
  fixf <- if (length(spec_order) == 0) terms else terms[-spec_order]

  # construct plain formula without f
  fm <- as.character(attr(tf, "variables")[[2]])
  fm <- paste(fm, "~", intercept, paste(c("", fixf), collapse = " + "))

  list(
    latents_in = latents_in,
    plain_fm = formula(fm)
    # ,
    # index_prd = index_prd,
    # index_est = index_est
  )
}

ngme_format <- function(param, val, model = NULL) {
  stationary <- (length(val) == 1)
  dne <- (length(val) == 0)

  if (is.null(model)) { # noise
    if (stationary)
      val <- if (param == "sigma") format(exp(val), digits = 3) else format(val, digits = 3)
    else
      val <-  paste0(format(val, digits = 3), collapse = ", ")

    switch(param,
      "sigma" = if (stationary) paste0("sigma = ", val)
                else paste0("theta_sigma = ", val),
      "mu"    = if (stationary) paste0("mu = ", val)
                else paste0("theta_mu = ", val),
      "nu"    = paste0("nu = ", val),
      "beta"  = if (dne) "No fixed effects" else paste0("beta = ", val)
    )
  } else { # model
    switch(model,
      "ar1"     = paste0("alpha = ", format(ar1_th2a(val), digits = 3)),
      "matern"  = if (stationary)
          paste0("kappa = ", format(exp(val), digits = 3))
        else
          paste0("theta_kappa = ", paste0(format(val, digits = 3), collapse = ", "))
    )
  }
}

#' taking mean over a list of nested lists
#'
#' @param lls a list
#'
#' @return a list of nested lists
#' @export
#'
#' @examples
#' ls <- list(
#'   list(a=1, b=2, t="nig", ll=list(a=1,b=2, w="ab")),
#'   list(a=3, b=5, t="nig", ll=list(a=1,b=6, w="ab")),
#'   list(a=5, b=5, t="nig", ll=list(a=4,b=2, w="ab"))
#' )
#' mean_list(ls)
mean_list <- function(lls) {
  # helpers
  nest_list_add <- function(l1, l2) {
    for (i in seq_along(l2)) {
      if (is.numeric(l2[[i]])) {
        l1[[i]] <- l1[[i]] + l2[[i]]
      }
      if (is.list(l2[[i]])) {
        l1[[i]] <- nest_list_add(l1[[i]], l2[[i]])
      }
    }
    l1
  }

  nest_list_divide <- function(l, n) {
    for (i in seq_along(l)) {
      if (is.numeric(l[[i]])) l[[i]] <- l[[i]] / n
      if (is.list(l[[i]]))    l[[i]] <- nest_list_divide(l[[i]], n)
    }
    l
  }

  l <- Reduce(nest_list_add, lls[-1], init = lls[[1]])
  nest_list_divide(l, length(lls))
}

# helper functions
# ar1 alpha (0~1) to theta_K
ar1_a2th <- function(a) {
  log((-1 - a) / (-1 + a))
}

# theta_K to ar1 alpha (0~1)
ar1_th2a <- function(th) {
 -1 + (2 * exp(th)) / (1 + exp(th))
}


# #' Make index for the matern model
# #'
# #' @param name
# #' @param n.spde
# #' @param n.repl
# #' @param mesh
# #' @param dim
# #'
# #' @return
# #' @export
# #'
# #' @examples
# ngme.matern.make.index <- function(
#   n.spde=NULL,
#   n.repl = 1,
#   mesh = NULL,
#   dim = NULL
# ){
#   if(is.null(n.spde)&&is.null(mesh)){
#     stop("You should provide either n.spde or mesh!")
#   }

#   if(!is.null(mesh)){
#     n_mesh = mesh$n

#     if(mesh$manifold == "R1"){
#       dim = 1
#     } else if(mesh$manifold == "R2"){
#       dim = 2
#     } else{
#       stop("The domain must be flat manifolds of dimension 1 or 2, that is,
#          the domain must be a line or a plane.")
#     }

#   } else{
#     n_mesh <- n.spde
#     if(is.null(dim)) {
#       stop("You should provide the dimension d!")
#     }
#   }

#   out <- list()
#   out$index <- rep(1:n_mesh, times = n.repl)

#   out$replicates <- rep(1:n.repl, each = n_mesh)
#   return(out)
# }