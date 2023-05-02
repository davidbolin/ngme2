
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
  } else if (inla_mesh$manifold == "R2") {
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

# split block by repl
split_block <- function(block, repl) {
  # split Y, X
  Ys <- split(block$Y, repl)
  Xs <- split_matrix(block$X, repl)

  # split latents A
  As <- lapply(block$latents, function(latent) {
    split_matrix(latent$A, repl)
  })

  latentss <- list()
  for (i in unique(repl)) {
    latentss[[i]] <- block$latents
    for (lat in seq_along(block$latents)) {
      latentss[[i]][[lat]]$A <- As[[lat]][[i]]
    }
  }

  blocks <- list()
  # build new blocks
  for (i in unique(repl)) {
    blocks[[i]] <- ngme_replicate(
      Y                 = Ys[[i]],
      X                 = Xs[[i]],
      latents           = latentss[[i]],
      beta              = block$beta,
      W_sizes           = block$W_sizes,
      V_sizes           = block$V_sizes,
      n_la_params       = block$n_la_params,
      n_params          = block$n_params,
      noise             = block$noise,
      seed              = block$seed,
      debug             = block$debug,
      control           = block$control
    )
  }
  blocks
}

split_matrix <- function(mat, repl) {
  split_mat <- lapply(split(mat, repl, drop = FALSE),
    matrix, ncol = ncol(mat))
  split_mat
}

# help to build a list of mesh for different replicates
build_mesh <- function(
  model,
  map,
  ...
) {
  if (model == "ar1") {
    mesh <- INLA::inla.mesh.1d(min(map):max(map))
  } else if (model == "tp") {
    # check later
    mesh = NULL
  } else {
    stop("Please provide mesh")
  }

  mesh
}

# helper function to unfiy way of accessing 1d and 2d index
sub_locs <- function(locs, idx) {
  if (inherits(locs, c("data.frame", "matrix"))) {
    locs[idx, , drop = FALSE]
  } else {
    locs[idx]
  }
}

dim_map <- function(map) {
  if (inherits(map, c("data.frame", "matrix"))) {
    2
  } else {
    1
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