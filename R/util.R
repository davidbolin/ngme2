
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
#' @param G matrix
#'
#' @return sparse dgCMatrix
#' @export
ngme_as_sparse <- function(G) {
  tryCatch(
    expr={
      G = as(as(G, "CsparseMatrix"), "generalMatrix")
      # G <- as(as(G, "dMatrix"), "generalMatrix")
    },
    error = function(e) {
      G <- as(G, "dgTMatrix")
      idx <- which(G@i <= G@j)
      G = Matrix::sparseMatrix(i=G@i[idx], j = G@j[idx], x= G@x[idx],
                               symmetric=FALSE, index1 = FALSE)
      G = as(as(G, "CsparseMatrix"), "generalMatrix")
      # G <- as(as(G, "dMatrix"), "generalMatrix")
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
#' @param index_NA index of unknown
#'
#' @return
#'  1. plain formula without f function
#'  2. latents_in - from each f function
ngme_parse_formula <- function(
  gf,
  data,
  index_NA
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
    if (!grepl("data *=", terms[i])) {
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

    # add information of index_NA
    str <- gsub("ngme2::f\\(", "ngme2::f(index_NA=index_NA,", str)
    # adding 1 term for furthur use in f
    # data$ngme_response <- Y
    res <- eval(parse(text = str), envir = data, enclos = parent.frame())

    # default name
    if (res$name == "field") res$name <- paste0("field", length(latents_in) + 1)
    # latents_in[[length(latents_in) + 1]] <- res
    latents_in[[res$name]] <- res
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

# format output
ngme_format <- function(param, val, model = NULL) {
  stationary <- (length(val) == 1)
  dne <- (length(val) == 0)

  if (is.null(model)) { # noise
    if (stationary)
      val <- if (grepl("sigma", param, fixed=TRUE))
        format(exp(val), digits = 3) else format(val, digits = 3)
    else
      val <-  paste0(format(val, digits = 3), collapse = ", ")

    switch(param,
      "sigma" = if (stationary) paste0("sigma = ", val)
                else paste0("theta_sigma = ", val),
      "sigma_nig" = if (stationary) paste0("sigma_nig = ", val)
                else paste0("theta_sigma_nig = ", val),
      "sigma_normal" = if (stationary) paste0("sigma_normal = ", val)
                else paste0("theta_sigma_normal = ", val),
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

#   out$replicate <- rep(1:n.repl, each = n_mesh)
#   return(out)
# }

# helper function to modify ngme_model's A and A_pred function
# using idx_NA
ngme_make_A <- function(
  mesh,
  map,
  n_map,
  idx_NA,
  replicate
) {
  # make it logical vector
  if (is.null(idx_NA))    idx_NA <- rep(FALSE, n_map)
  if (is.numeric(idx_NA)) idx_NA <- 1:n_map %in% idx_NA

  # make A and A_pred according to mesh
  A_pred <- NULL
  if (is.numeric(mesh) && is.null(dim(mesh)))
    mesh <- INLA::inla.mesh.1d(mesh)
  if (inherits(mesh, "inla.mesh.1d")) {
    A <- INLA::inla.spde.make.A(mesh = mesh, loc = map[!idx_NA], repl=replicate[!idx_NA])
    if (any(idx_NA)) A_pred <- INLA::inla.spde.make.A(mesh = mesh, loc = map[idx_NA])
  } else if (inherits(mesh, "inla.mesh")) {
    # 2d location
    A <- INLA::inla.spde.make.A(mesh = mesh, loc = map[!idx_NA, ,drop = FALSE], repl=replicate[!idx_NA])
    if (any(idx_NA)) A_pred <- INLA::inla.spde.make.A(mesh = mesh, loc = map[idx_NA, ,drop = FALSE])
  } else {
    stop("this mesh not implement yet!!!")
  }

  list(A = A, A_pred = A_pred)
}

#' Make observation matrix for time series
#'
#' @param loc   integers (after sorting, no gaps > 1)
#' @param replicate indicating replicate measure at same location
#' @param range range for the mesh
#'  by default range=(min(loc), max(loc))
#'
#' @return A matrix (length(loc) * length(unique(loc)))
#' @export
#'
#' @examples
#' ngme_ts_make_A(c(1, 2, 2), replicate = c(1, 1, 2))
#' ngme_ts_make_A(c(1, 2, 2), range = c(1, 5))
ngme_ts_make_A <- function(
  loc,
  replicate = NULL,
  range = c(1, max(loc))
) {
  if (is.null(loc) || length(loc) == 0) return (NULL)

  n_loc <- length(loc)
  nrep <- 1

  start <- range[1]; end <- range[2]
  n_range <- end - start + 1

  if (is.null(replicate)) {
    replicate <- rep(1, n_range)
  }

  unique_rep <- unique(replicate)
  nrep <- length(unique_rep)

  A <- matrix(0, nrow = n_loc, ncol = n_range * nrep)

  for (i in 1:n_loc) {
    ncol_rep <- which(unique_rep == replicate[i])
    A[i, (ncol_rep - 1) * n_range + loc[i] - start + 1] <- 1
  }
  # as(A, "dgCMatrix")
  as(as(A, "dMatrix"), "generalMatrix")
}

# compute the mode of data
emprical_mode <- function(x, breaks = max(20, length(x) / 20)) {
  h <- hist(x, breaks=breaks, plot=FALSE)
  idx <- which.max(h$counts)
  h$mids[idx]
}

# build <- function(A1, A2) {
# }

# t, (x, y)

# f(mesh = list(mesh1, mesh2), replicate=112233) + f(ar1 )


# given a list of replicate
merge_repls <- function(repls) {
  # merge the list of data frames
  # input: list of numerics
  # output: merged list
  # assert of equal length
  stopifnot(length(unique(as.numeric(lapply(repls, length)))) == 1)
  # helper function of merge 2 repls
  merge_repl <- function(repl, group) {
    halas <- FALSE
    while (!halas) {
      unique_group <- unique(group)
      halas <- TRUE
      if (length(unique_group) == 1) return (group)
      for (i in 1:(length(unique_group)-1)) {
        for (j in (i+1):length(unique_group)) {
          if (length(intersect(repl[group == unique_group[i]], repl[group == unique_group[j]])) > 0) {
            group[group == unique_group[j]] <- unique_group[i]
            halas <- FALSE
          }
          if (!halas) break
        }
        if (!halas) break
      }
    }
    group
  }

  Reduce(function(x, y) merge_repl(x, y), repls)
}
