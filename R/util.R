
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
  } else if (inla_mesh$manifold %in% c("R2", "S2")) {
    d <- 2
  } else {
    stop("The mesh should be from a flat manifold.")
  }
  return(d)
}

# format output
ngme_format <- function(param, val, model = NULL, ...) {
  stationary <- (length(val) == 1)
  dne <- (length(val) == 0)

  if (is.null(model)) { # noise
    if (stationary)
      val <- if (grepl("sigma", param, fixed=TRUE) || grepl("nu", param, fixed=TRUE))
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
      "nu"    = if (stationary) paste0("nu = ", val)
                else paste0("theta_nu = ", val),
      "feff"  = if (dne) "No fixed effects" else paste0("feff = ", val)
    )
  } else { # model
    switch(model,
      "ar1"     = paste0("alpha = ", format(ar1_th2a(val), digits = 3)),
      "matern"  = paste0("theta_kappa = ", paste0(format(val, digits = 3), collapse = ", ")),
      "ou"      = paste0("theta_K = ", paste0(format(val, digits = 3), collapse = ", ")),
      "re"      = {
        invisible(print(vecK_to_Sigma(val, list(...)[[1]])))
      }
    )
  }
}

#' taking mean over a list of nested lists
#'
#' @param lls a list
#' @param weights weights of each list
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
mean_list <- function(lls, weights=NULL) {
  n <- length(lls)
  weights <- if (is.null(weights)) rep(1 / n, n)
    else weights / sum(weights)

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
  nest_list_mult <- function(l, n) {
    for (i in seq_along(l)) {
      if (is.numeric(l[[i]])) l[[i]] <- l[[i]] * n
      if (is.list(l[[i]]))    l[[i]] <- nest_list_mult(l[[i]], n)
    }
    l
  }

  ret <- nest_list_mult(lls[[1]], 0)
  for (i in seq_along(lls)) {
    tmp <- nest_list_mult(lls[[i]], weights[[i]])
    ret <- nest_list_add(ret, tmp)
  }
  ret
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
  range = c(min(loc), max(loc))
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
# merge_repls <- function(repls) {
#   # merge the list of data frames
#   # input: list of numerics
#   # output: merged list
#   # assert of equal length
#   stopifnot(length(unique(as.numeric(lapply(repls, length)))) == 1)
#   # helper function of merge 2 repls
#   merge_repl <- function(repl, group) {
#     halas <- FALSE
#     while (!halas) {
#       unique_group <- unique(group)
#       halas <- TRUE
#       if (length(unique_group) == 1) return (group)
#       for (i in 1:(length(unique_group)-1)) {
#         for (j in (i+1):length(unique_group)) {
#           if (length(intersect(repl[group == unique_group[i]], repl[group == unique_group[j]])) > 0) {
#             group[group == unique_group[j]] <- unique_group[i]
#             halas <- FALSE
#           }
#           if (!halas) break
#         }
#         if (!halas) break
#       }
#     }
#     group
#   }

#   Reduce(function(x, y) merge_repl(x, y), repls)
# }

split_matrix <- function(mat, repl) {
  split_mat <- lapply(split(mat, repl, drop = FALSE),
    matrix, ncol = ncol(mat))
  split_mat
}

# helper function to unfiy way of accessing 1d and 2d index
as_map <- function(locs) {
  if (inherits(locs, c("data.frame", "matrix"))) {
    as.matrix(locs)
  } else {
    as.numeric(locs)
  }
}

sub_map <- function(locs, idx) {
  if (inherits(locs, c("data.frame", "matrix"))) {
    locs[idx, , drop = FALSE]
  } else {
    locs[idx]
  }
}

dim_map <- function(map) {
  if (inherits(map, c("data.frame", "matrix"))) ncol(map) else 1
}

rep_map <- function(map, times) {
  if (inherits(map, c("data.frame", "matrix"))) {
    do.call(rbind, replicate(times, map, simplify=FALSE))
  } else {
    rep(map, times = times)
  }
}

length_map <- function(map) {
  if (inherits(map, c("data.frame", "matrix"))) {
    nrow(map)
  } else {
    length(map)
  }
}

# vectorize a matrix
veci <- function(v, n, m) {
  if (length(v) != n * m) {
    cat("Wrong dimensions in reshape:", length(v), "(", n, ",", m, ")\n")
  }
  M <- matrix(0, nrow = n, ncol = m)
  count <- 1
  for (i in 1:m) {
    M[, i] <- v[count:(count+n-1)]
    count <- count + n
  }
  return(M)
}

# M <- matrix(c(1,2,3,2,4,5,3,5,6), 3, 3); M
vech <- function (M) {
  stopifnot(nrow(M) == ncol(M))
  n = nrow(M)
  V <- double(n * (n - 1) / 2 + n)
	V[1:n] = diag(M)

  k = n+1
  for (i in seq_len(n-1)) {
    V[k : (k + n-i-1)] <- tail(M[, i], n - i)
    k = k + n - i
  }
  return(V)
}

vech_to_mat <- function(vech, n) {
  stopifnot(n*(n-1)/2 + n == length(vech))
  mat <- matrix(0, nrow = n, ncol = n)
  mat[lower.tri(mat)] <- tail(vech, length(vech) - n)
  mat <- mat + t(mat)
  diag(mat) <- vech[1:n]
  mat
}

vecK_to_Sigma <- function(theta_K, n) {
  tmp <- theta_K; tmp[1:n] <- exp(theta_K[1:n])
  K <- vech_to_mat(tmp, n)
  Kinv = solve(K)
  Sigma = Kinv %*% t(Kinv)
  Sigma
}

build_D <- function(theta, rho) {
  d11 <- cos(theta) + rho * sin(theta)
  d12 <- -sin(theta) * (sqrt(1+rho^2))
  d21 <- sin(theta) - rho * cos(theta)
  d22 <- cos(theta) * (sqrt(1+rho^2))
  matrix(c(d11, d21, d12, d22), nrow = 2, ncol = 2)
}

build_effect_K <- function(n_reff, theta_K) {
  n_theta_K <- sum(1:n_reff)
  K <- diag(n_reff); diag(K) <- exp(theta_K[1:n_reff])
  if (n_reff > 1)
    K[lower.tri(K)] <- theta_K[(n_reff+1):n_theta_K]
  K
}

#' @title ngme prior specification
#' @description
#' Prior specification for internal representation of the parameters.
#' We will list all the available priors here.
#' Their PDFs can be found in vignette("ngme2-prior").
#'
#' Flat prior with no parameter: \eqn{p(\theta) \propto 1}
#'
#' Normal prior with parameter mean \eqn{\mu} and precision \eqn{\tau}: \eqn{\theta \sim N(\mu, 1/ \tau)}
#'
#' Penalized complexity prior with parameter \eqn{\lambda}
#'
#' @param type type of prior, say ngme_prior_types()
#' @param param parameters of the prior
#'
#' @return a list of prior specification
#' @export
ngme_prior <- function(type, param=NULL) {
  stopifnot("Please check if the prior name is in ngme_prior_types()" =
    type %in% ngme2::ngme_prior_types())

  # check if num. of parameter is correct
  switch(type,
    "flat" = stopifnot(length(param) == 0),
    "normal" = stopifnot(length(param) == 2),
    # internally log precision
    "pc.sd" = {
      # lambda
      # param = - log(alpha) / u
      stopifnot(length(param) == 1)
    },
    "pc.cor0" = stopifnot(length(param) == 1)
  )

  structure(
    list(type = type, param = param),
    class = "ngme_prior"
  )
}

# given a basis matrix, find if it is stationary
is_stationary <- function(B) {
  stopifnot(is.matrix(B))
  ncol(B) == 1 && all(B == 1)
}


#' @title ngme make mesh for different replicates
#' @description
#' Make different mesh for different replicates
#'
#' @param data provide the data.frame
#' @param map provide the map to make mesh, i.g. ~x+y, x and y will be extracted from data to make a 2d mesh
#' @param replicate provide the replicate information, i.g. ~id
#' @param mesh_type type of mesh, "regular" means use all the point from same replicate to make a mesh
#'
#' @return a list of mesh of length of different replicates
#' @export
ngme_make_mesh_repls <- function(
  data,
  map,
  replicate,
  mesh_type = "regular"
) {

  if (inherits(map, "formula")) {
    map <- model.matrix(map, data)[, -1]
  }

  if (inherits(replicate, "formula")) {
    replicate <- model.matrix(replicate, data)[, -1]
  }

  stopifnot(length_map(map) == length(replicate))

  mesh_repls <- NULL
  for (repl in replicate) {
    map_repl <- subset(map, replicate == repl)
    if (dim_map(map) == 1)
      mesh_repls[[as.character(repl)]] <- tryCatch(
        fmesher::fm_mesh_1d(map_repl),
        error = function(e) {
          stop("The nodes for making mesh is not valid for replicate id=" , repl)
        }
      )
    else if (dim_map(map) == 2) {
      stop("Not implemented yet.")
      mesh_repls[[repl]] <- fmesher::fm_mesh_2d(map_repl)
    } else {
      stop("The dimension of the mesh should be 1 or 2.")
    }
  }

  mesh_repls
}


# check if group is valid
# return as integer
validate_rep_or_group <- function(replicate, data) {
  if (is.null(replicate))
    return (rep(1, nrow(data)))

  if (inherits(replicate, "formula")) {
    # input as: replicate = ~id
    stopifnot("Allow 1 variable (column in data) as replicate. i.g. replicate=~id"
      = length(replicate) == 2 && length(replicate[[2]]) == 1)

    replicate <- eval(replicate[[2]], envir = data, enclos = parent.frame())
  }

  if (inherits(replicate, "character") && (length(replicate) == 1)) {
    # input as: replicate = "id"
    replicate <- data[[replicate]]
  }

  stopifnot(
    "Please make sure the length of replicate is equal to the number of rows of data"
     = nrow(data) == length(replicate)
  )

  return (as.factor(replicate))
}

#' @title Use more interpretable parameterization of the matern model
#' @description
#' From SPDE parameter (kappa, sigma, alpha) to (theta_kappa, theta_sigma, theta_alpha)
#'
#' @param ope provide the operator
#'
#' @return a list of mesh of length of different replicates
#' @export
matern_result <- function(ope) {
  stop("to-do")
}