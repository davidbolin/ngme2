#' Specifying a latent process model (wrapper function for each model)
#'
#' Function used for defining of smooth and spatial terms
#' within ngme model formulae.
#' The function is a wrapper function for specific submodels.
#' (see ngme_models_types() for available models).
#'
#' @param map  symbol or numerical value: index or covariates to build index
#' @param model  string, model type, see ngme_model_types()
#' @param noise  ngme_noise object, noise_nig() or noise_gal()
#' @param name   name of the field, for later use, if not provided, will be "field1" etc.
#' @param data      specifed or inherit from ngme() function
#' @param group   group factor indicate resposne variable, can be inherited from ngme() function, (used for bivariate model)
#' @param which_group  belong to which group
#' @param replicate factor indicating replicate structure for INLA-style replicates.
#'   When provided, operators and A matrices will be block-diagonalized across replicates.
#' @param A  observation matrix, automatically computed given map and model
#' @param W      starting value of the process
#' @param fix_K  fix the estimation of parameter of K
#' @param prior prior specification created by \code{prior_*()} or
#'   \code{priors(...)} for operator parameters (\code{theta_K}).
#' @param fix_W  stop sampling for W
#' @param debug     debug mode
#' @param subset    subset of the model
#'
#' @details When using different meshes for different replicates, provide the mesh parameter
#' as a list of mesh objects. The number of meshes should match the number of replicates.
#' For example: \code{mesh = list(mesh1, mesh2, mesh3)} for 3 replicates.
#'
#' When using the `replicate` argument, the operator matrices will be block-diagonalized.
#' For a generic operator K = param1 * Matrix1 + param2 * Matrix2, with 2 replicates,
#' the resulting operator becomes:
#' K = param1 * bdiag(Matrix1_rep1, Matrix1_rep2) + param2 * bdiag(Matrix2_rep1, Matrix2_rep2)
#' Similarly, the A matrix will be block-diagonalized as bdiag(A_rep1, A_rep2, ...).
#'
#' @return a list for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. Information about K matrix
#' 2. Information about noise
#' 3. Control variables
#'
#' @export
f <- function(
    map,
    model,
    noise = noise_normal(),
    name = "field",
    data = NULL,
    group = NULL,
    which_group = NULL,
    replicate = NULL,
    A = NULL,
    W = NULL,
    fix_W = FALSE,
    fix_K = FALSE,
    prior = NULL,
    subset = rep(TRUE, length_map(map)),
    debug = FALSE) {
  # examine the noise
  if (is.list(noise) && !inherits(noise, "ngme_noise")) {
    noise <- convert_noise_list_to_normal_nig(noise)
  }

  # helper: check whether all provided noises are normal
  is_noise_all_normal <- function(noise_obj) {
    if (inherits(noise_obj, "ngme_noise")) {
      return(noise_obj$noise_type == "normal")
    }
    if (is.list(noise_obj) && length(noise_obj) > 0) {
      return(all(vapply(noise_obj, function(x) inherits(x, "ngme_noise") && x$noise_type == "normal", logical(1))))
    }
    FALSE
  }

  maybe_apply_default_nu_prior <- function(noise_obj, h_vec) {
    if (isTRUE(noise_obj$fix_theta_nu) || noise_obj$n_theta_nu == 0) {
      return(noise_obj)
    }
    if (isTRUE(noise_obj$prior_nu_user)) {
      return(noise_obj)
    }
    if (!any(noise_obj$noise_type %in% c("nig", "normal_nig"))) {
      return(noise_obj)
    }

    h_star <- stats::median(h_vec, na.rm = TRUE)
    if (!is.finite(h_star) || h_star <= 0) {
      return(noise_obj)
    }

    lambda <- log(2) / h_star
    if (!is.finite(lambda) || lambda <= 0) {
      return(noise_obj)
    }

    lower <- if (!is.null(noise_obj$nu_lower_bound)) noise_obj$nu_lower_bound else 0
    noise_obj$prior_nu <- as_internal_prior(
      prior_inv_exp(lambda = lambda, lower = lower, target = "coef")
    )
    noise_obj
  }

  noise_all_normal <- is_noise_all_normal(noise)
  # If the user builds a bv/bv2/bv_matern model inline and all noises are normal,
  # ensure fix_theta = TRUE in that model call.
  model_expr <- substitute(model)
  # Check if model is a call to bv/bv2/bv_matern
  if (noise_all_normal && is.call(model_expr)) {
    op_name <- tryCatch(as.character(model_expr[[1]]), error = function(e) "")
    if (op_name %in% c("bv", "bv2", "bv_matern")) {
      model_list <- as.list(model_expr)
      name_vec <- names(model_list)
      idx <- which(name_vec %in% "fix_theta")

      if (length(idx) > 0 && isFALSE(model_list[[idx]])) {
        stop("For bv/bv2/bv_matern models with all normal noises, fix_theta must be TRUE. Please set fix_theta = TRUE.")
      } else if (length(idx) == 0) {
        model_list <- c(model_list, list(fix_theta = TRUE))
        model_expr <- as.call(model_list)
      } else { # fix_theta is present and not FALSE
        model_list[[idx]] <- TRUE
        model_expr <- as.call(model_list)
      }
    }
  }

  # Evaluate model expression (promises are forced here so we can inspect it consistently)
  model <- eval(model_expr, envir = parent.frame())
  model_type <- model$model
  map <- eval(substitute(map), envir = data, enclos = parent.frame())

  # Evaluate map expression (promises are forced here so we can inspect it consistently)
  if (inherits(map, "formula")) {
    map <- get_data_from_formula(map, data)
  } else if (is.list(map)) {
    map <- lapply(map, function(m) {
      if (inherits(m, "formula")) {
        get_data_from_formula(m, data)
      } else {
        m
      }
    })
  }

  if (!is.null(data) && is.null(A)) {
    # check if model is a string or an object

    if (!model_type %in% c("tp", "spacetime")) {
      stopifnot(
        "Please make sure length of map is same as nrow(data) in f()" =
          length_map(map) == nrow(data)
      )
    } else {
      stopifnot(
        "Please make sure length of map is same as nrow(data) in f()" =
          length_map(map[[1]]) == nrow(data) &&
            length_map(map[[2]]) == nrow(data)
      )
    }
  }

  group <- validate_rep_or_group(group, data)

  # Validate and process replicate argument
  if (!is.null(replicate)) {
    replicate <- validate_rep_or_group(replicate, data)
    if (!is.null(data)) {
      stopifnot(
        "Length of replicate must match number of observations" =
          length(replicate) == nrow(data)
      )
    } else {
      stopifnot(
        "Length of replicate must match length of map" =
          length(replicate) == length_map(map)
      )
    }
  }

  # set the subset if provide group and which_group
  if (!is.null(which_group)) {
    stopifnot(
      "Please provide group factor" = !is.null(group),
      "Please check if which_group is in group" = which_group %in% levels(group)
    )
    subset <- group %in% which_group
  }

  stopifnot(
    "Please provide model from ngme_model_types():" = !is.null(model),
    "Please specify model as ngme_operator or ngme_operator_def" = inherits(model, "ngme_operator") || inherits(model, "ngme_operator_def")
  )

  # 1. Extract or build mesh
  mesh <- NULL
  mesh_list <- NULL

  if (inherits(model, "ngme_operator")) {
    mesh <- model$mesh
  } else if (inherits(model, "ngme_operator_def")) {
    if (model$model == "re") {
      model$args$map <- map
      mesh <- NULL
    } else if (!is.null(model$args$mesh)) {
      mesh <- model$args$mesh
    } else {
      # Build mesh from map if not provided
      mesh <- ngme_build_mesh(sub_map(map, subset), model$model)
      model$args$mesh <- mesh
    }
  } else {
    stop("Model must be ngme_operator or ngme_operator_def")
  }

  # 2. Check if mesh is a list (for replicates)
  is_replicate_mesh <- is.list(mesh) && !inherits(mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))

  if (is_replicate_mesh && model_type %in% c("tp", "spacetime")) {
    # Check if it is a list of meshes for the model itself (not replicates)
    if (inherits(mesh[[1]], c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) {
      is_replicate_mesh <- FALSE
    }
  }

  if (is_replicate_mesh) {
    mesh_list <- mesh
    # Use the first mesh as a template for building the operator
    mesh <- mesh_list[[1]]

    # Update model definition to use the template mesh
    if (inherits(model, "ngme_operator_def")) {
      model$args$mesh <- mesh
    }
  }

  # 3. Convert ngme_operator_def to ngme_operator
  if (inherits(model, "ngme_operator_def")) {
    operator <- do.call(model$model, model$args)
  } else {
    operator <- model
  }
  model_name <- operator$model

  # At this point, 'operator' is an ngme_operator object (possibly a template)
  # and 'model' is the model name string.

  # Build A matrix
  A <- if (is.null(A)) {
    # If mesh_list is present, we defer A matrix construction for replicates
    if (!is.null(mesh_list)) {
      NULL
    } else {
      ngme_build_A(model_name, mesh, map, operator, group)
    }
  } else {
    ngme_as_sparse(A)
  }

  # subset the A matrix only if it was built
  if (!is.null(A) && !all(subset)) {
    A <- A[subset, , drop = FALSE]
  }

  # Sanity warning for ARMA with fixed W and free MA/AR params
  if (model_name == "arma" && isTRUE(fix_W)) {
    ar_mask <- if (!is.null(operator$fix_rho)) operator$fix_rho else logical(0)
    ma_mask <- if (!is.null(operator$fix_phi)) operator$fix_phi else logical(0)
    if ((length(ar_mask) == 0 || !all(ar_mask)) || (length(ma_mask) == 0 || !all(ma_mask))) {
      warning("For ARMA with fix_W=TRUE, ar/ma determine Z and W=P^{-1}x; fixing W while estimating ar/ma leads to mismatch. Either (i) set fix_ar/fix_ma=TRUE and pass W=P^{-1}x, or (ii) allow W to be estimated.")
    }
  }

  # Handle replicates: block-diagonalize operator matrices and A matrix
  if (!is.null(replicate)) {
    n_repl <- length(levels(replicate))

    # Build separate operators for each replicate using their respective meshes
    if (!is.null(mesh_list)) {
      operator_list <- list()
      A_list <- list()
      h_total <- NULL

      for (i in 1:n_repl) {
        # Get indices for this replicate
        repl_idx <- which(replicate == levels(replicate)[i])

        # Get map for this replicate
        if (model_name %in% c("tp", "spacetime")) {
          map_repl <- list(map[[1]][repl_idx], map[[2]][repl_idx])
        } else {
          map_repl <- sub_map(map, repl_idx)
        }

        # Get mesh for this replicate
        mesh_repl <- if (length(mesh_list) >= i) mesh_list[[i]] else mesh_list[[1]]

        # Build operator for this replicate with its specific mesh
        # If model is ngme_operator, we might need to clone it and update mesh?
        # Or if it's ngme_operator_def, we build it.
        # The original code called build_operator.

        if (inherits(model, "ngme_operator")) {
          # If it's already an operator, we can't easily change the mesh if it's baked in?
          # Usually ngme_operator has a mesh.
          # If the user passed a single operator but multiple meshes in f(),
          # we need to re-instantiate the operator with the new mesh.
          # But we don't have the constructor arguments easily if it's already an operator.
          # This suggests that for replicates with different meshes, passing ngme_operator_def is better.

          # If it is an operator, we can try to update the mesh if the operator supports it.
          operator_repl <- model
          operator_repl$mesh <- mesh_repl
          # We might need to recompute K or other things?
          # This is tricky. If the user provides a fully built operator, they might expect it to be used as is.
          # But if mesh changes, the operator matrices change.

          warning("Using a pre-built ngme_operator with replicate and different meshes. The mesh in the operator will be replaced, but internal matrices might not update if they depend on the mesh during construction. Consider using ngme_operator_def (e.g. ar1(..., mesh=NULL)) instead.")
        } else if (inherits(model, "ngme_operator_def")) {
          args <- model$args
          args$mesh <- NULL
          operator_repl <- do.call(model$model, c(list(mesh = mesh_repl), args))
        } else {
          stop("Model must be ngme_operator or ngme_operator_def.")
        }

        operator_list[[i]] <- operator_repl

        # Build A matrix for this replicate
        if (is.null(A)) {
          A_repl <- ngme_build_A(
            model_name, mesh_repl, map_repl, operator_repl,
            if (!is.null(group)) group[repl_idx] else NULL
          )
          A_list[[i]] <- A_repl
        }

        # Accumulate h vectors
        h_total <- c(h_total, operator_repl$h)
      }

      # Create block-diagonal operator matrices
      if (operator_list[[1]]$generic_type %in% c("generic", "generic_ns")) {
        # For generic operators, block-diagonalize each matrix type
        n_matrices <- length(operator_list[[1]]$matrices)
        combined_matrices <- list()

        for (j in 1:n_matrices) {
          matrix_list <- lapply(operator_list, function(op) op$matrices[[j]])
          combined_matrices[[j]] <- do.call(Matrix::bdiag, matrix_list)
        }

        operator$matrices <- combined_matrices

        # Update K matrix using the new block-diagonal matrices
        if (!is.null(operator$update_K)) {
          operator$K <- operator$update_K(operator$theta_K)
        } else {
          # Block-diagonalize K matrices
          K_list <- lapply(operator_list, function(op) op$K)
          operator$K <- do.call(Matrix::bdiag, K_list)
        }
      } else {
        # For non-generic operators, block-diagonalize K directly
        K_list <- lapply(operator_list, function(op) op$K)
        operator$K <- do.call(Matrix::bdiag, K_list)

        # If the operator has component matrices, block-diagonalize them too
        if (!is.null(operator_list[[1]]$C)) {
          operator$C <- do.call(Matrix::bdiag, lapply(operator_list, function(op) op$C))
        }
        if (!is.null(operator_list[[1]]$G)) {
          operator$G <- do.call(Matrix::bdiag, lapply(operator_list, function(op) op$G))
        }
        if (!is.null(operator_list[[1]]$Ci)) {
          operator$Ci <- do.call(Matrix::bdiag, lapply(operator_list, function(op) op$Ci))
        }
        if (!is.null(operator_list[[1]]$B_K)) {
          operator$B_K <- do.call(rbind, lapply(operator_list, function(op) op$B_K))
        }
        if (!is.null(operator_list[[1]]$Z)) {
          operator$Z <- do.call(Matrix::bdiag, lapply(operator_list, function(op) op$Z))
        }

        # Update update_K to return block-diagonal K for shared parameters
        if (!is.null(operator_list[[1]]$update_K)) {
          operator$update_K <- function(theta_K) {
            K_list <- lapply(operator_list, function(op) op$update_K(theta_K))
            do.call(Matrix::bdiag, K_list)
          }
        }
      }

      # Update h vector
      operator$h <- h_total

      # Block-diagonalize A matrices if they were built
      if (!is.null(A_list) && length(A_list) > 0) {
        A <- do.call(Matrix::bdiag, A_list)
      }
    } else {
      # Original logic for same mesh across replicates
      # Block-diagonalize operator matrices if it's a generic operator
      if (operator$generic_type %in% c("generic", "generic_ns")) {
        # Create block-diagonal versions of all matrices
        operator$matrices <- lapply(operator$matrices, function(mat) {
          # Create list of matrices for each replicate
          mat_list <- replicate(n_repl, mat, simplify = FALSE)
          # Block diagonalize
          do.call(Matrix::bdiag, mat_list)
        })

        # Upda simple cases, just block-diagonalize K directly
        K_list <- replicate(n_repl, operator$K, simplify = FALSE)
        operator$K <- do.call(Matrix::bdiag, K_list)
        # Update h vector (concatenate for each replicate)
        operator$h <- rep(operator$h, n_repl)
      } else {
        stop("Not implemented replicate for non-generic operators yet.")
        # For non-generic operators, block-diagonalize K directly
        K_list <- replicate(n_repl, operator$K, simplify = FALSE)
        operator$K <- do.call(Matrix::bdiag, K_list)
        operator$h <- rep(operator$h, n_repl)
      }

      # Block-diagonalize A matrix if it exists
      if (!is.null(A)) {
        # Build A matrix for each replicate
        A_list <- list()
        for (i in 1:n_repl) {
          # Get indices for this replicate
          repl_idx <- which(replicate == levels(replicate)[i])

          # Get map for this replicate
          if (model_name %in% c("tp", "spacetime")) {
            map_repl <- list(map[[1]][repl_idx], map[[2]][repl_idx])
          } else {
            map_repl <- sub_map(map, repl_idx)
          }

          # Build A matrix for this replicate using the same mesh
          A_repl <- ngme_build_A(
            model_name, mesh, map_repl, operator,
            if (!is.null(group)) group[repl_idx] else NULL
          )
          A_list[[i]] <- A_repl
        }

        # Block-diagonalize A matrices
        A <- do.call(Matrix::bdiag, A_list)
      }
    }

    # Apply subset if needed
    if (!is.null(A) && !all(subset)) {
      A <- A[subset, , drop = FALSE]
    }

    # Update W if provided (replicate for each replicate)
    if (!is.null(W)) {
      W <- rep(W, n_repl)
    }
  }

  # 2. build noise given operator
  # bivariate noise
  # 2. build noise given operator
  # bivariate noise
  if (model_name %in% c("bv", "bv_matern_normal", "bv2", "bv_matern_nig")) {
    stopifnot(
      "Please specify noise for each field" = length(noise) >= 2,
      "Input: noise=list(a=<noise>,b=<noise>)" = inherits(noise[[1]], "ngme_noise"),
      "Input: noise=list(a=<noise>,b=<noise>)" = inherits(noise[[2]], "ngme_noise"),
      "Please specify noise with same name as in the sub_models argument!" = all(names(noise[1:2]) %in% operator$model_names),
      "Keep the noise same if you want to specify single V for each noise" = noise[[1]]$single_V == noise[[2]]$single_V,
      "Two noise should be the same type" = noise[[1]]$noise_type == noise[[2]]$noise_type
    )

    noise1 <- update_noise(noise[[operator$model_names[[1]]]],
      n = length(operator$h) / 2
    )
    noise2 <- update_noise(noise[[operator$model_names[[2]]]],
      n = length(operator$h) / 2
    )
    bv_noises <- list(noise1, noise2)
    names(bv_noises) <- operator$model_names
    share_V <- !is.null(noise$share_V) && noise$share_V
    if (share_V) {
      stopifnot(
        "share_V option is only supported for 2 NIG noise." = noise1$noise_type == noise2$noise_type,
        "share_V option requires nu from both noise are same." = noise1$theta_nu == noise2$theta_nu
      )
    }

    noise <- ngme_noise(
      noise_type = c(noise1$noise_type, noise2$noise_type),
      B_mu = as.matrix(Matrix::bdiag(noise1$B_mu, noise2$B_mu)),
      B_sigma = as.matrix(Matrix::bdiag(noise1$B_sigma, noise2$B_sigma)),
      B_nu = as.matrix(Matrix::bdiag(noise1$B_nu, noise2$B_nu)),
      theta_mu = c(noise1$theta_mu, noise2$theta_mu),
      theta_sigma = c(noise1$theta_sigma, noise2$theta_sigma),
      theta_nu = c(noise1$theta_nu, noise2$theta_nu),
      fix_theta_nu = noise1$fix_theta_nu && noise2$fix_theta_nu,
      fix_theta_mu = noise1$fix_theta_mu && noise2$fix_theta_mu,
      fix_theta_sigma = c(noise1$fix_theta_sigma, noise2$fix_theta_sigma),
      share_V = !is.null(noise$share_V) && noise$share_V,
      single_V = noise1$single_V,
      bv_noises = bv_noises,
      fix_V = !is.null(noise$fix_V) && noise$fix_V,
      V = if (!is.null(noise$V)) noise$V else NULL
    )
  } else {
    n <- length(operator$h)
    if (model_name %in% c("rw1", "rw2") && operator$cyclic) {
      # Compensate the constraint for cyclic rw1/rw2 (we introduce lagrange in definition)
      if (nrow(noise$B_sigma) > 1) { # user-defined B_sigma
        noise$B_sigma <- rbind(0, noise$B_sigma)
      }
      if (nrow(noise$B_mu) > 1) { # user-defined B_mu
        noise$B_mu <- rbind(0, noise$B_mu)
      }
      if (nrow(noise$B_nu) > 1) { # user-defined B_nu
        noise$B_nu <- rbind(0, noise$B_nu)
      }
    }
    noise <- update_noise(noise, n = n)
  }

  if (noise$share_V &&
    !(model_name %in% c("bv", "bv2", "bv_matern_normal", "bv_matern_nig"))) {
    stop("Not allow for share_V for univariate model")
  }

  if (model_name %in% c("rw1", "rw2")) {
    # force the constrain (e.g. fixed the 1st position to be 0)
    theta_sigma <- c(-10, noise$theta_sigma)
    fix_theta_sigma <- c(TRUE, noise$fix_theta_sigma)
    B_sigma <- cbind(0, noise$B_sigma)
    B_sigma[1, 1] <- 1
    noise$B_sigma[1, -1] <- 0
    if (any(abs(eigen(t(B_sigma) %*% B_sigma)$values) < 1e-10)) {
      # cannot apply transformation
      stop("Do not fix the 1st position of sigma for non-cyclic rw1")
    } else {
      noise$B_sigma <- B_sigma
      noise$theta_sigma <- theta_sigma
      noise$fix_theta_sigma <- fix_theta_sigma
    }

    if (model_name == "rw2" && !operator$cyclic) {
      noise$B_sigma[2, 1] <- 1
      noise$B_sigma[2, -1] <- 0
      noise$theta_sigma[1] <- -6 # set larger for numerical stability
    }

    if (operator$cyclic) A <- cbind(0, A) # ignore the 1st auxiliary element
  }

  # Reshape the noise structure for normal_nig
  if (all(noise$noise_type %in% c(
    "normal_nig", "normal_gal", "nig_gal"
  ))) {
    # update operator matrices M to be diag(M, M)
    if (operator$generic_type %in% c("generic", "generic_ns")) {
      operator$matrices <- lapply(operator$matrices, function(x) {
        Matrix::bdiag(x, x)
      })
      operator$K <- operator$K <- Matrix::bdiag(operator$K, operator$K)
      operator$h <- operator$h <- c(operator$h, operator$h)
    }
    # A <- [A A]
    A <- cbind(A, A)

    W <- c(W, W)
  }

  noise <- maybe_apply_default_nu_prior(noise, operator$h)

  operator_prior_names <- operator$param_name
  if (is.null(operator_prior_names) ||
      length(operator_prior_names) != length(operator$theta_K)) {
    operator_prior_names <- paste0("theta", seq_along(operator$theta_K))
  }
  prior_theta_K <- compile_operator_priors(prior, operator_prior_names)

  ngme_model(
    model = model_name,
    operator = operator,
    noise = noise,
    W_size = ncol(operator$K),
    V_size = nrow(operator$K),
    A = A,
    map = map,
    mesh = mesh,
    mesh_list = mesh_list,
    n_map = length_map(map),
    W = W,
    fix_W = fix_W,
    name = name,
    fix_K = fix_K,
    prior_theta_K = prior_theta_K,
    debug = debug,
    replicate = replicate
  )
}

# build operator


# Convert a list containing both nig and normal noise to a normal_nig noise object
convert_noise_list_to_normal_nig <- function(noise_list) {
  # Check if we have a list with both nig and normal noise types
  if (!(length(noise_list) == 2 &&
    ((inherits(noise_list[[1]], "ngme_noise") && noise_list[[1]]$noise_type == "nig" &&
      inherits(noise_list[[2]], "ngme_noise") && noise_list[[2]]$noise_type == "normal") ||
      (inherits(noise_list[[1]], "ngme_noise") && noise_list[[1]]$noise_type == "normal" &&
        inherits(noise_list[[2]], "ngme_noise") && noise_list[[2]]$noise_type == "nig")))) {
    return(noise_list)
  }

  # Extract the nig and normal components
  nig_noise <- if (noise_list[[1]]$noise_type == "nig") noise_list[[1]] else noise_list[[2]]
  normal_noise <- if (noise_list[[1]]$noise_type == "normal") noise_list[[1]] else noise_list[[2]]

  # Create a normal_nig noise object
  noise_normal_nig(
    mu = nig_noise$theta_mu,
    sigma_nig = exp(nig_noise$theta_sigma),
    sigma_normal = exp(normal_noise$theta_sigma),
    nu = nig_noise$nu_lower_bound + exp(nig_noise$theta_nu),
    nu_lower_bound = nig_noise$nu_lower_bound,
    B_mu = nig_noise$B_mu,
    B_sigma_nig = nig_noise$B_sigma,
    B_sigma_normal = normal_noise$B_sigma,
    B_nu = nig_noise$B_nu,
    fix_theta_mu = nig_noise$fix_theta_mu,
    fix_theta_sigma = c(
      nig_noise$fix_theta_sigma,
      normal_noise$fix_theta_sigma
    ),
    fix_theta_nu = nig_noise$fix_theta_nu,
    corr_measurement = nig_noise$corr_measurement || normal_noise$corr_measurement,
    V = nig_noise$V,
    fix_V = nig_noise$fix_V
  )
}

# help to build a list of mesh for different replicates
ngme_build_mesh <- function(
    loc,
    model = NULL,
    ...) {
  if (inherits(loc, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) {
    return(loc)
  }

  # Check if loc is a list of meshes
  if (is.list(loc) && length(loc) > 0 && all(vapply(loc, function(x) inherits(x, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph")), logical(1)))) {
    return(loc)
  }

  if (!is.null(model)) {
    model_name <- if (inherits(model, "ngme_operator")) model$model else model
    if (model_name %in% c("re", "tp")) {
      return(NULL)
    }
    if (model_name == "iid") loc <- as.integer(as.factor(loc))
    if (model_name %in% c("ar", "ar1", "arma")) {
      stopifnot(
        "The map should be integers." = is.numeric(loc) && all(loc == round(loc))
      )
      return(fmesher::fm_mesh_1d(loc = min(loc):max(loc)))
      # return (fmesher::fm_mesh_1d(as.integer(as.factor(loc))))
    }
  }

  if (is.matrix(loc) && ncol(loc) == 2) {
    stop("Please build and provide the mesh for spatial data using fmesher::fm_mesh_2d()")
  } else if (is.numeric(loc)) {
    mesh <- fmesher::fm_mesh_1d(loc = loc)
  } else {
    stop("The mesh provided is invalid.")
  }

  mesh
}
