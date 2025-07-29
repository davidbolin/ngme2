#' Create Time Series Cross-Validation Indices
#' 
#' @description Creates indices for time series cross-validation with options for expanding window
#' or sliding window approaches. Supports both single-step and multi-step forecasting, and can
#' handle replicated observations.
#' 
#' @param time_idx A numeric vector of time indices in ascending order
#' @param train_length An integer specifying the fixed length of training sets. 
#'        If NULL (default), an expanding window approach is used.
#' @param test_length An integer specifying the number of observations to include in each test set.
#'        Default is 1 (single-step forecasting).
#' @param replicate An optional vector of the same length as time_idx, indicating which observations
#'        belong to the same replicate group. When provided, ensures that all observations with the
#'        same replicate value are either entirely in the training set or entirely in the test set.
#' @param gap An integer specifying the gap between the training set and test set.
#'        Default is 0 (no gap, test set starts immediately after training set).
#'        For example, gap=1 means skip one time point between training and test (useful for 2-step ahead forecasting).
#' 
#' @return A list with two components:
#' \describe{
#'   \item{train}{A list of numeric vectors, where each vector contains the time indices for training in that fold}
#'   \item{test}{A list of numeric vectors, where each vector contains the time indices for testing in that fold}
#' }
#' 
#' @details
#' Time series cross-validation requires respecting the temporal order of observations.
#' This function implements two common approaches:
#' 
#' 1. Expanding window (when train_length = NULL): The training set grows with each fold, 
#'    starting with a minimal set and expanding to include all but the test data.
#' 
#' 2. Sliding window (when train_length is specified): Uses a fixed-length window 
#'    that slides through the time series, maintaining the same training size across folds.
#' 
#' The test_length parameter allows for multi-step forecasting evaluation.
#' 
#' When replicate is provided, the function ensures that all observations with the same 
#' replicate value are kept together, either all in the training set or all in the test set.
#' This is useful for scenarios where multiple observations at the same time point should
#' be treated as a group.
#' 
#' The gap parameter creates a separation between training and test sets, which is useful
#' for multi-step ahead forecasting validation.
#' 
#' @examples
#' # Expanding window approach with single-step forecasting
#' cv_expanding <- make_time_series_cv_index(1:10)
#' 
#' # Sliding window approach with window size 3 and single-step forecasting
#' cv_sliding <- make_time_series_cv_index(1:10, train_length = 3)
#' 
#' # Sliding window with multi-step forecasting (predict 2 steps ahead)
#' cv_multistep <- make_time_series_cv_index(1:10, train_length = 3, test_length = 2)
#' 
#' # Working with replicates
#' time_idx <- c(1, 1, 1, 2, 2, 3, 3)
#' replicates <- c(1, 1, 1, 2, 2, 3, 3)
#' cv_with_replicates <- make_time_series_cv_index(time_idx, replicate = replicates)
#' 
#' # 2-step ahead forecasting with a gap
#' cv_with_gap <- make_time_series_cv_index(1:10, gap = 1)
#' 
#' @export
make_time_series_cv_index <- function(
  time_idx, 
  train_length = NULL,
  test_length = 1,
  replicate = time_idx,
  gap = 0
) {
  # Validate inputs
  if (!is.numeric(time_idx)) {
    stop("time_idx must be numeric")
  }
  
  n <- length(time_idx)
  if (n < 2) {
    stop("time_idx must have at least 2 elements")
  }
  
  # Validate test_length
  if (test_length < 1) {
    stop("test_length must be at least 1")
  }
  
  # Validate gap
  if (gap < 0) {
    stop("gap must be non-negative")
  }
  
  # Initialize lists for train and test indices
  train_indices <- list()
  test_indices <- list()
  
  # Validate and process replicate information if provided
  if (!is.null(replicate)) {
    if (length(replicate) != length(time_idx)) {
      stop("replicate must have the same length as time_idx")
    }
    # Create a data frame to keep track of time points and their replicates
    data <- data.frame(
      idx = 1:length(time_idx),
      time = time_idx,
      replicate = replicate
    )
    
    # Sort by time and then by replicate to ensure proper ordering
    data <- data[order(data$time, data$replicate), ]
    
    # Get unique replicates in order of their first appearance
    unique_replicates <- unique(data$replicate)
    replicate_first_times <- tapply(data$time, data$replicate, min)
    replicate_order <- order(replicate_first_times)
    sorted_unique_replicates <- unique_replicates[replicate_order]
    
    # We'll work with replicate groups instead of individual time points
    time_groups <- sort(unique(replicate_first_times))
    n_groups <- length(time_groups)
    
    # Time series CV with replicates - we'll work with groups
    if (is.null(train_length)) {
      # Expanding window mode with replicates
      # We'll create folds by using all earlier replicates for training and later ones for testing
      
      # Calculate the maximum number of folds considering the gap
      max_folds <- n_groups - gap - test_length
      
      # Ensure we have at least one fold
      if (max_folds < 1) {
        stop("Not enough replicate groups for the specified gap and test_length")
      }
      
      # For each possible train/test split of replicate groups
      for (i in 1:max_folds) {
        # Identify which replicates are in training set
        train_replicates <- sorted_unique_replicates[1:i]
        
        # Get all rows where the replicate is in the training set
        train_idx <- data$idx[data$replicate %in% train_replicates]
        train_indices[[i]] <- time_idx[train_idx]
        
        # Get the test replicates after the gap
        test_replicates <- sorted_unique_replicates[
          (i+1+gap):min(i+gap+test_length, n_groups)
        ]
        
        # Get all rows where the replicate is in the test set
        test_idx <- data$idx[data$replicate %in% test_replicates]
        test_indices[[i]] <- time_idx[test_idx]
      }
    } else {
      # Fixed-length sliding window with replicates
      
      # Validate train_length - now in terms of number of replicate groups
      if (train_length <= 0 || train_length >= n_groups) {
        stop("train_length must be between 1 and the number of replicate groups - 1")
      }
      
      # The maximum number of folds we can create, accounting for the gap
      max_folds <- n_groups - train_length - gap - test_length + 1
      
      # Ensure we have at least one fold
      if (max_folds < 1) {
        stop("Not enough replicate groups for the specified train_length, gap, and test_length")
      }
      
      # Generate sliding window CV folds based on replicate groups
      for (i in 1:max_folds) {
        # Training set: replicates from position i to i+train_length-1
        train_replicates <- sorted_unique_replicates[i:(i+train_length-1)]
        
        # Get all rows where the replicate is in the training set
        train_idx <- data$idx[data$replicate %in% train_replicates]
        train_indices[[i]] <- time_idx[train_idx]
        
        # Test set: the next test_length replicates after the gap
        test_replicates <- sorted_unique_replicates[
          (i+train_length+gap):(min(i+train_length+gap+test_length-1, n_groups))
        ]
        
        # Get all rows where the replicate is in the test set
        test_idx <- data$idx[data$replicate %in% test_replicates]
        test_indices[[i]] <- time_idx[test_idx]
      }
    }
  } else {
    # Standard time series CV without replicates
    # Sort time indices to ensure chronological order if they aren't already
    if (!all(diff(time_idx) >= 0)) {
      time_idx <- sort(unique(time_idx))
      n <- length(time_idx)
    }
    
    if (is.null(train_length)) {
      # Expanding window mode (default behavior)
      # We'll create as many folds as possible while ensuring enough test points
      # and accounting for the gap
      max_folds <- n - gap - test_length
      
      if (max_folds < 1) {
        stop("Not enough data points for the specified gap and test_length")
      }
      
      for (i in 1:max_folds) {
        # Training set: indices from 1 to i
        train_indices[[i]] <- time_idx[1:i]
        
        # Test set: observations after training and gap
        test_start <- i + 1 + gap
        test_end <- min(test_start + test_length - 1, n)
        test_indices[[i]] <- time_idx[test_start:test_end]
      }
    } else {
      # Fixed-length sliding window mode
      
      # Validate train_length
      if (train_length <= 0 || train_length >= n) {
        stop("train_length must be between 1 and length(time_idx)-1")
      }
      
      # The maximum number of folds we can create, accounting for the gap
      max_folds <- n - train_length - gap - test_length + 1
      
      # Ensure we have at least one fold
      if (max_folds < 1) {
        stop("Not enough data points for the specified train_length, gap, and test_length")
      }
      
      # Generate sliding window CV folds
      for (i in 1:max_folds) {
        # For sliding window: start with indices i:(i+train_length-1)
        start_idx <- i
        end_idx <- i + train_length - 1
        
        # Training set: window of length train_length
        train_indices[[i]] <- time_idx[start_idx:end_idx]
        
        # Test set: observations after training and gap
        test_start <- end_idx + 1 + gap
        test_end <- min(test_start + test_length - 1, n)
        test_indices[[i]] <- time_idx[test_start:test_end]
      }
    }
  }
  
  # Return the result as a list containing train and test indices
  return(list(
    train = train_indices,
    test = test_indices
  ))
}



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
# ar1 alpha (-1~1) to theta_K
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
  if (is.null(data)) {
    return (as.factor(replicate))
  }

  if (is.null(replicate))
    replicate <- rep(1, nrow(data))

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
    "Please make sure the length of replicate/group is equal to the number of rows of data"
     = nrow(data) == length(replicate)
  )

  return (as.factor(replicate))
}

# #' @title Use more interpretable parameterization of the matern model
# #' @description
# #' From SPDE parameter (kappa, sigma, alpha) to (theta_kappa, theta_sigma, theta_alpha)
# #'
# #' @param ope provide the operator
# #'
# #' @return a list of mesh of length of different replicates
# #' @export
# matern_result <- function(ope) {
#   stop("to-do")
# }


# return the relative idx among all latent models
# using W_size assume W_size == V_size
idx_range <- function(ngme_rep, name_or_idx) {
  if (is.character(name_or_idx)) {
    idx <- which(sapply(ngme_rep$models, function(model) model$name) == name_or_idx)
  } else {
    idx <- name_or_idx
  }

  models <- ngme_rep$models
  tmp <- sapply(models, function(model) model$W_size)
  sizes = Reduce(`+`, tmp, accumulate = TRUE)

  if (idx == 1) {
    return (1:sizes[1])
  } else {
    return ((sizes[idx-1]+1):sizes[idx])
  }
}

#' @title posterior samples of different latent models
#' @description
#' Extract the posterior samples of different latent models
#'
#' @param ngme_object ngme object
#' @param model_name name of the model, or index of the model
#' @param type type of samples, "W" or "V"
#' @param replicate which replicate
#'
#' @return a data.frame of posterior samples (mesh_size * n_post_samples)
#' @export
ngme_post_samples <- function(
  ngme_object,
  model_name=1,
  type = "W",
  replicate = 1
) {
  ngme_rep <- ngme_object$replicates[[replicate]]
  stopifnot(
    inherits(ngme_rep, "ngme_replicate"),
    type %in% c("W", "V")
  )
  idx <- idx_range(ngme_rep, model_name)

  if (type == "W") {
    return (ngme_rep$post_W[idx, , ])
  } else {
    return (ngme_rep$post_V[idx, , ])
  }
}

#' @title variance of the data or the latent field
#' @description
#' Compute the variance of the data or the latent field
#'
#' @param ngme ngme_model
#' @param model_name
#'   if model_name = "data", then return the covariance matrix of the data (without measurement noise)
#'   if the model_name is the name or index of the latent, then return the covariance matrix of the latent field
#' @param replicate which replicate (default = 1)
#'
#' @return a data.frame of posterior samples (mesh_size * n_post_samples)
#' @export
ngme_cov_matrix <- function(
  ngme_object,
  model_name = "data",
  replicate = 1
) {
  ngme_rep <- ngme_object$replicates[[replicate]]

  stopifnot(
    "Please provide the correct ngme object." = inherits(ngme_object, "ngme"),
    "Please provide the correct model name or index." =
      model_name %in% c(
        "data",
        sapply(ngme_rep$models, function(model) model$name),
        1:length(ngme_rep$models)
      )
  )

  if (model_name == "data") {
    V_mean = apply(ngme_rep$post_V, 1, mean)

    diag_K = Matrix::bdiag(
      sapply(ngme_rep$models, function(model) model$operator$K))

    Q <- diag_K %*% diag(1 / V_mean) %*% Matrix::t(diag_K)
    block_A <- Reduce(cbind, sapply(ngme_rep$models, function(model) model$A))

    # return A Q^(-1) A^t,
    return (block_A %*% Matrix::solve(Q, Matrix::t(block_A)))
  } else {
    K <- ngme_rep$models[[model_name]]$operator$K

    idx <- idx_range(ngme_rep, model_name)
    V_mean = apply(ngme_rep$post_V[idx, , ], 1, mean)

    Q <- K %*% diag(1 / V_mean) %*% Matrix::t(K)
    return (solve(Q))
  }
}


ngme_build_A <- function(model, mesh, map, operator, group, group_levels=NULL) {
  group <- validate_rep_or_group(group, NULL)
  
  if (is.null(group_levels)) group_levels <- levels(group)

  if (inherits(mesh, "metric_graph")) {
    A <- if (is.null(map)) NULL else mesh$fem_basis(map)
    return (A)
  }

  if (model %in% c("tp", "spacetime")) {
    mesh_t <- if (model=="tp") operator$first$mesh else operator$mesh[[1]]
    mesh_s <- if (model=="tp") operator$second$mesh else operator$mesh[[2]]
    stopifnot("Now only support first to be 1d model"
      = inherits(mesh_t, "inla.mesh.1d"))

    # watch-out! Should use as.factor to start at 1, not min(map[[1]])
    # blk_group <- as.integer(as.factor(map[[1]]))
    
    blk_group <- as.integer(map[[1]])
    # important to start with 1
    min_loc <- mesh_t$interval[[1]]
    blk_group <- blk_group - min_loc + 1
    stopifnot("The map[[1]] of the tp model is not valid, out of range." = all(blk_group <= (mesh_t$range[2] + 1)))

    blk <- fmesher::fm_block(blk_group, n_block = mesh_t$n)
    basis <- fmesher::fm_basis(mesh_s, loc=map[[2]])
    A0 = fmesher::fm_row_kron(Matrix::t(blk), basis)

    if (model=="tp" && operator$second$model %in% c("bv", "bv_matern_normal", "bv_normal", "bv_matern_nig")) {
      # for tp-bv model
      # check if group is valid
      if (length(group) == 0) stop("Please provide the `group` argument.")
      all(group %in% group_levels) || stop("The group is not valid.")
      length(group) == length_map(map[[1]]) || stop("The length of group should be equal to the length of map.")
      # tp bv model
      # 1. expand A <- cbind(A, 0), double the column
      # 2. move 2nd field to the 2nd half
      # 3. re-order the 1st and 2nd
      # e.g., (each field is of size 2)
      # 1 2 3 4 5 6 7 8 to
      # 1 2 5 6 3 4 7 8
      A_expand <- cbind(
        A0,
        matrix(0, nrow=nrow(A0), ncol=ncol(A0))
      )
      
      # select row (of the 2nd field)
        row_2nd_field <- group == group_levels[[2]]
        half_1st <- 1:ncol(A0)
        half_2nd <- 1:ncol(A0) + ncol(A0)
      
      # Move the 2nd field to the right
        A_expand[row_2nd_field, half_2nd] <-
          A_expand[row_2nd_field, half_1st]

        if (any(row_2nd_field)) 
          A_expand[row_2nd_field, half_1st] <- 0

      # Re-order (f1 f2 f1 f2 ...)
        n <- operator$first$mesh$n
        bv_mesh_size <-  operator$second$mesh$n

        idx <- 1:ncol(A_expand)
        field_1st_idx <- ceiling(idx / bv_mesh_size) %% bv_mesh_size == 1
        field_2nd_idx <- !field_1st_idx
        reorder_idx = c(which(field_1st_idx), which(field_2nd_idx))
      
      # return after re-order
        return (ngme_as_sparse(A_expand[, reorder_idx]))
    } else {
      return (ngme_as_sparse(A0))
    }
  }
  
  # bivariate model
  if (model %in% c("bv", "bv_normal", "bv_matern_normal", "bv_matern_nig")) {
    # check if group is valid
    if (length(group) == 0) stop("Please provide the `group` argument.")
    all(group %in% group_levels) || stop("The group is not valid.")
    length(group) == length_map(map) || stop("The length of group should be equal to the length of map.")

    blk_group <- as.integer(as.factor(group))
    blk <- fmesher::fm_block(blk_group)
    basis <- fmesher::fm_basis(mesh, loc=map)
    A <- fmesher::fm_row_kron(Matrix::t(blk), basis)
    return (ngme_as_sparse(A))
  }

  if (model == "re") {
    return (ngme_as_sparse(operator$B_theta_K))
  }

  if (is.factor(map)) map <- as.integer(map) # for iid model
  return (fmesher::fm_basis(mesh, loc = map))
}




contain_bv_model <- function(ngme){
  any(sapply(ngme$replicates[[1]]$models, function(model) {
    model$model %in% c("bv", "bv_normal", "bv_matern_normal", "bv_matern_nig")
  })) || 
  any(sapply(ngme$replicates[[1]]$models, function(model) {
    model$model %in% c("tp") && model$operator$second$model %in% c("bv", "bv_normal", "bv_matern_normal", "bv_matern_nig")
  }))
}


#' @title get the trajectories of parameters of the model
#' @description
#' Get the trajectories of the parameters of the model
#'
#' @param ngme_object ngme object
#' @param model_name name of the model, if NULL, return the trajectories of the parameters of the measurement noise and fixed effects
#' @return a list of trajectories, each element is a matrix of trajectories (n_samples * n_trajectories)
#' @export
get_trajectories <- function(
  ngme_object,
  model_name
) {
  stopifnot(inherits(ngme_object, "ngme"))

  ngme <- ngme_object$replicates[[1]]
  if (!is.null(model_name) && model_name %in% names(ngme$models)) {
    traj <- attr(ngme$models[[model_name]], "lat_traj")
    ts <- get_latent_info(ngme$models[[model_name]])
  } else {
    traj <- attr(ngme, "block_traj")
    ts <- get_noise_info(ngme$noise)
  }
  #  str(traj)
  #  List of 4
  #  $ :'data.frame':       4 obs. of  1000 variables:
  
  convert_to_3d_array <- function(traj) {
    # Extract the dimensions
    num_replicates <- length(traj) # Number of replicates
    num_variables <- nrow(traj[[1]]) # Number of variables (4)
    num_iterations <- ncol(traj[[1]]) # Number of iterations (1000)
    
    # Initialize an empty 3D array
    array_3d <- array(NA, dim = c(num_variables, num_iterations, num_replicates))
    
    # Fill the array
    for (replicate_index in seq_len(num_replicates)) {
      array_3d[,,replicate_index] <- as.matrix(traj[[replicate_index]])
    }
    
    return(array_3d)
  }
  traj_3d <- convert_to_3d_array(traj)

  ret <- list()
  for (i in seq_along(ts$name)) {
    name <- ts$name[i]
    data <- traj_3d[i, , ]
    ret[[name]] <- ts$trans[[i]](data)
  }

  return (ret)
}




#' @title get the parameters of the noise
#' @description
#' Get the parameters of the noise
#'
#' @param noise noise object
#' @param prefix prefix of the parameter
#' @return a list with parameter name as key and parameter value as value
get_noise_param <- function(noise, prefix) {
  B      <- noise[[paste0("B_", prefix)]]
  theta  <- noise[[paste0("theta_", prefix)]]
  if (length(theta) == 0) return (NULL)

  if (is_stationary(B)) {
    value <- if (prefix == "mu") theta else exp(theta)
  } else {
    value <- theta
  }

  name <- if (is_stationary(B)) prefix else paste0("theta_", prefix)
  result <- list()
  result[[name]] <- value
  result
}



# extract parameters from all models in ngme object
extract_parameters <- function(ngme_object) {
  stopifnot(inherits(ngme_object, "ngme"))
  ngme_1rep <- ngme_object$replicates[[1]]
  result_transformed <- list() # transformed parameters (user friendly)
  result_raw <- list() # raw parameters (for optimizer)

  # LATENT MODELS
  models <- ngme_1rep$models
  model_names <- names(models)

  for (j in seq_along(models)) {
    names <- models[[j]]$operator$param_name
    trans <- models[[j]]$operator$param_trans
    result_model_transformed <- list()
    result_model_raw <- list()
    for (i in seq_along(names)) {
      result_model_transformed[[names[i]]] <- trans[[i]](models[[j]]$operator$theta_K[[i]])
      if (is.function(trans[[i]]) && trans[[i]](3) == 3) {
        raw_name <- names[i]
      } else {
        raw_name <- paste0("theta_", names[i])
      }
      result_model_raw[[raw_name]] <- models[[j]]$operator$theta_K[[i]]
    }

    # Use get_noise_param for noise parameters
    noise <- models[[j]]$noise
    noise_params_transformed <- list()
    noise_params_raw <- list()
    for (p in c("mu", "sigma", "nu")) {
      noise_param <- get_noise_param(noise, p)
      
      # special case for rw1 and rw2
      if (models[[j]]$model %in% c("rw1", "rw2") && p == "sigma") {
        noise_param$theta_sigma <- noise_param$theta_sigma[-1]
        if (length(noise_param$theta_sigma) == 1) {
          noise_param <- list(sigma = exp(noise_param$theta_sigma))
        }
      }

      if (!is.null(noise_param)) {
        noise_params_transformed <- c(noise_params_transformed, noise_param)
        noise_params_raw <- c(noise_params_raw, noise_param)
      }
    }

    result_model_transformed <- c(result_model_transformed, noise_params_transformed)
    result_model_raw <- c(result_model_raw, noise_params_raw)

    result_transformed[[model_names[j]]] <- result_model_transformed
    result_raw[[model_names[j]]] <- result_model_raw
  }

  # FIXED EFFECTS and MEASUREMENT NOISE
  fixed_effects <- ngme_1rep$feff
  result_transformed$data$fixed_effects <- fixed_effects
  result_raw$data$fixed_effects <- fixed_effects

  m_noise <- ngme_1rep$noise
  m_noise_params_transformed <- list()
  m_noise_params_raw <- list()
  for (p in c("mu", "sigma", "nu")) {
    noise_param <- get_noise_param(m_noise, p)
    if (!is.null(noise_param)) {
      m_noise_params_transformed <- c(m_noise_params_transformed, noise_param)
      m_noise_params_raw <- c(m_noise_params_raw, noise_param)
    }
  }

  result_transformed$data <- c(result_transformed$data, m_noise_params_transformed)
  result_raw$data <- c(result_raw$data, m_noise_params_raw)

  list(transformed = result_transformed, raw = result_raw)
}


