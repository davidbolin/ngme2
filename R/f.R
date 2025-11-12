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
#' @param mesh   mesh for the model, if not provided, will be built from map. 
#'   Can be a single mesh object or a list of mesh objects for different replicates.
#'   When using replicates, provide mesh as a list where mesh[[i]] corresponds to replicate i.
#' @param name   name of the field, for later use, if not provided, will be "field1" etc.
#' @param data      specifed or inherit from ngme() function
#' @param group   group factor indicate resposne variable, can be inherited from ngme() function, (used for bivariate model)
#' @param which_group  belong to which group
#' @param replicate factor indicating replicate structure for INLA-style replicates. 
#'   When provided, operators and A matrices will be block-diagonalized across replicates.
#' @param A  observation matrix, automatically computed given map and model except for Bayesian regression
#' @param W      starting value of the process
#' @param fix_W  stop sampling for W
#' @param fix_theta_K fix the estimation for theta_K.
#' @param prior_theta_K prior for theta_K
#' @param debug     debug mode
#' @param subset    subset of the model
#' @param ...       additional arguments (e.g. parameters for model)
#'  inherit the data from ngme function
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
  noise       = noise_normal(),
  mesh        = NULL,
  name        = "field",
  data        = NULL,
  group       = NULL,
  which_group = NULL,
  replicate   = NULL,
  A           = NULL,
  W           = NULL,
  fix_W       = FALSE,
  fix_theta_K = FALSE,
  prior_theta_K = ngme_prior("normal", param=c(0, 0.001)),
  subset      = rep(TRUE, length_map(map)),
  debug       = FALSE,
  ...
) {
  # examine the noise 
  if (is.list(noise) && !inherits(noise, "ngme_noise")) {
    noise <- convert_noise_list_to_normal_nig(noise)
  }
  
  stopifnot(
    "Please provide noise as ngme_noise object" = 
      model %in% c("bv", "bv_normal", "bv_matern_normal", "bv_matern_nig") ||
      inherits(noise, "ngme_noise")
  )

  if ((missing(map) || (is.null(map))) && inherits(mesh, "metric_graph")) {
    stopifnot("To use metric graph model, please install MetricGraph package"
      = rlang::is_installed("MetricGraph"))

    # extract the map
    graph_data <- (tryCatch(mesh$get_data(), error=function(e) NULL))
    map <- if (is.null(graph_data)) NULL
      else with(graph_data, cbind(`.edge_number`, `.distance_on_edge`))
  }

  map <- eval(substitute(map), envir = data, enclos = parent.frame())

  if (inherits(map, "formula")) {
    if (model == "re") {
      map <- model.matrix(map, data)
    } else {
      map <- model.matrix(map, data)[, -1]
    }
  }

  if (!is.null(data) && is.null(A)) {
    if (! model %in% c("tp", "spacetime")) stopifnot(
      "Please make sure length of map is same as nrow(data) in f()" =
      length_map(map) == nrow(data)
    ) else {
      "Please make sure length of map is same as nrow(data) in f()" =
      length_map(map[[1]]) == nrow(data) &&
      length_map(map[[2]]) == nrow(data)
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
      "Please check if which_group is in group"
        = which_group %in% levels(group))
    subset <- group %in% which_group
  }

  stopifnot(
    "Please provide model from ngme_model_types():" = !is.null(model),
    "Please specify model as character" = is.character(model),
    "prior_theta_K is not specified properly, please use ngme_prior(..)"
      = class(prior_theta_K) == "ngme_prior"
  )

  # Validate mesh input - can be single mesh, list of meshes, or NULL
  mesh_list <- NULL
  if (!is.null(mesh)) {
    if (
      model != "spacetime" &&
      is.list(mesh) &&
      !inherits(mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))
    ) {
      # mesh is a list of meshes for different replicates
      mesh_list <- mesh
      mesh <- NULL  # Set to NULL to defer mesh selection to parsing stage
      stopifnot(
        "All elements in mesh list must be valid mesh objects" =
        all(sapply(mesh_list, function(m) inherits(m, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))))
      )
    }
  }

  # 0. build mesh if not specified
  if (is.null(mesh) && is.null(A) && is.null(mesh_list)) {
    mesh <- ngme_build_mesh(sub_map(map, subset), model)
  }

  # remove NULL in arguments
  f_args <- Filter(Negate(is.null),  as.list(environment()))
  # add arguments in ...
  f_args <- c(f_args, list(...))

  # Remove mesh_list from f_args since it's not needed for operator building
  f_args$mesh_list <- NULL
  # Remove replicate from f_args for now since individual operators don't handle it
  f_args$replicate <- NULL
  
  # If we have mesh_list, remove mesh from f_args so operator building doesn't fail
  if (!is.null(mesh_list)) {
    f_args$mesh <- NULL
  }

  # if (model %in% c("tp", "spacetime")) {
  if (model %in% c("tp")) {
    stopifnot(
      "Please specify map as a list of length two (for 2 sub-models)" =
      is.list(map) && length(map) == 2
    )
    
    map <- lapply(map, function(x) {
      if (inherits(x, "formula")) {
        model.matrix(x, data)[, -1]
      } else {
        x
      }
    })

    # examine argument for tp model
    first <- list(...)$first; second <- list(...)$second
    stopifnot(
      "The length of map of 2 sub_models should be same"
        = length_map(map[[1]]) == length_map(map[[2]]),
      "Please provide f(first = ..., second = ...), see ?tp."
        = !is.null(first) && !is.null(second),
      "Please make sure first is a list and second is a list"
        = is.list(first) && is.list(second),
      "Please provide the model argument in first and second"
        = !is.null(first$model) && !is.null(second$model),
      "Please provide mesh individually for first and second"
        = is.null(mesh)
    )

    if (is.null(first$mesh)) first$mesh <- ngme_build_mesh(map[[1]], model=first$model)
    if (is.null(second$mesh)) second$mesh <- ngme_build_mesh(map[[2]],
    model=second$model)
    f_args$first <- if (inherits(first, "ngme_operator")) first else
      build_operator(first$model, first)
    f_args$second <- if (inherits(second, "ngme_operator")) second else
      build_operator(second$model, second)
  }

  # build the operator
  if (!is.null(mesh_list)) {
    # Use the first mesh as a template for building the operator
    # The correct mesh will be set during parsing stage
    temp_mesh <- mesh_list[[1]]
    f_args$mesh <- temp_mesh
    operator <- build_operator(model, f_args)
    
    # Defer A matrix construction to parsing stage since we used wrong mesh
    A <- NULL
  } else {
    operator <- build_operator(model, f_args)
    
    # Build A matrix 
    A <- if (is.null(A)) ngme_build_A(model, mesh, map, operator, group)
          else ngme_as_sparse(A)
    
    # subset the A matrix only if it was built
    if (!is.null(A) && !all(subset)) {
      A <- A[subset, , drop = FALSE]
    }
  }

  # Sanity warning for ARMA with fixed W and free MA/AR params
  if (model == "arma" && isTRUE(fix_W)) {
    ar_mask <- if (!is.null(operator$fix_rho)) operator$fix_rho else logical(0)
    ma_mask <- if (!is.null(operator$fix_phi)) operator$fix_phi else logical(0)
    if ((length(ar_mask)==0 || !all(ar_mask)) || (length(ma_mask)==0 || !all(ma_mask))) {
      warning("For ARMA with fix_W=TRUE, ar/ma determine Z and W=P^{-1}x; fixing W while estimating ar/ma leads to mismatch. Either (i) set fix_ar/fix_ma=TRUE and pass W=P^{-1}x, or (ii) allow W to be estimated.")
    }
  }
  
  # Handle replicates: block-diagonalize operator matrices and A matrix
  if (!is.null(replicate)) {
    n_repl <- length(levels(replicate))
    
    # If we have mesh_list, we need to build operators and A matrices for each replicate separately
    if (!is.null(mesh_list)) {
      # Build separate operators for each replicate using their respective meshes
      operator_list <- list()
      A_list <- list()
      h_total <- NULL
      
      for (i in 1:n_repl) {
        # Get indices for this replicate
        repl_idx <- which(replicate == levels(replicate)[i])
        
        # Get map for this replicate
        if (model %in% c("tp", "spacetime")) {
          map_repl <- list(map[[1]][repl_idx], map[[2]][repl_idx])
        } else {
          map_repl <- sub_map(map, repl_idx)
        }
        
        # Get mesh for this replicate
        mesh_repl <- if (length(mesh_list) >= i) mesh_list[[i]] else mesh_list[[1]]
        
        # Build operator for this replicate with its specific mesh
        f_args_repl <- f_args
        f_args_repl$mesh <- mesh_repl
        f_args_repl$map <- map_repl
        operator_repl <- build_operator(model, f_args_repl)
        operator_list[[i]] <- operator_repl
        
        # Build A matrix for this replicate
        if (is.null(A)) {
          A_repl <- ngme_build_A(model, mesh_repl, map_repl, operator_repl, 
                                if (!is.null(group)) group[repl_idx] else NULL)
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
          if (model %in% c("tp", "spacetime")) {
            map_repl <- list(map[[1]][repl_idx], map[[2]][repl_idx])
          } else {
            map_repl <- sub_map(map, repl_idx)
          }
          
          # Build A matrix for this replicate using the same mesh
          A_repl <- ngme_build_A(model, mesh, map_repl, operator, 
                                if (!is.null(group)) group[repl_idx] else NULL)
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
  if (model %in% c("bv", "bv_matern_normal", "bv_normal", "bv_matern_nig")) {
    stopifnot(
      "Please specify noise for each field" = length(noise) >= 2,
      "Input: noise=list(a=<noise>,b=<noise>)" = inherits(noise[[1]], "ngme_noise"),
      "Input: noise=list(a=<noise>,b=<noise>)" = inherits(noise[[2]], "ngme_noise"),
      "Please specify noise with same name as in the sub_models argument!"
        = all(names(noise[1:2]) %in% operator$model_names),
      "Keep the noise same if you want to specify single V for each noise"
        = noise[[1]]$single_V == noise[[2]]$single_V,
      "Two noise should be the same type"
        = noise[[1]]$noise_type == noise[[2]]$noise_type
    )

# make sure the noise is in the same order as the model_names
if (noise[[1]]$noise_type == "normal") {
  stopifnot(
    "Please use model=bv_normal/bv_matern_normal for Gaussian noise (then rotation is fixed)" = model %in% c("bv_normal", "bv_matern_normal")
  )
}

if (noise[[1]]$noise_type != "normal") {
  stopifnot(
    "Please use model=bv for non-Gaussian noise" = 
      model %in% c("bv", "bv_matern_nig")
  )
}
    noise1 <- update_noise(noise[[operator$model_names[[1]]]],
      n=length(operator$h)/2)
    noise2 <- update_noise(noise[[operator$model_names[[2]]]],
      n=length(operator$h)/2)
    bv_noises <- list(noise1, noise2); names(bv_noises) <- operator$model_names
    share_V <- !is.null(noise$share_V) && noise$share_V
    if (share_V) {
      stopifnot(
        "share_V option is only supported for 2 NIG noise."
          = noise1$noise_type == noise2$noise_type,
        "share_V option requires nu from both noise are same."
          = noise1$theta_nu == noise2$theta_nu
      )
    }

    noise <- ngme_noise(
      noise_type = c(noise1$noise_type, noise2$noise_type),
      B_mu    = as.matrix(Matrix::bdiag(noise1$B_mu, noise2$B_mu)),
      B_sigma = as.matrix(Matrix::bdiag(noise1$B_sigma, noise2$B_sigma)),
      B_nu    = as.matrix(Matrix::bdiag(noise1$B_nu, noise2$B_nu)),
      theta_mu    = c(noise1$theta_mu, noise2$theta_mu),
      theta_sigma = c(noise1$theta_sigma, noise2$theta_sigma),
      theta_nu    = c(noise1$theta_nu, noise2$theta_nu),
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
    if (model %in% c("rw1", "rw2") && operator$cyclic) {
      # Compensate the constraint for cyclic rw1/rw2 (we introduce lagrange in definition)
      if (nrow(noise$B_sigma) > 1)   # user-defined B_sigma
        noise$B_sigma <- rbind(0, noise$B_sigma)
      if (nrow(noise$B_mu) > 1)   # user-defined B_mu
        noise$B_mu <- rbind(0, noise$B_mu)
      if (nrow(noise$B_nu) > 1)   # user-defined B_nu
        noise$B_nu <- rbind(0, noise$B_nu)
    }
    noise <- update_noise(noise, n = n)
  }

  if (noise$share_V && 
    !(model %in% c("bv", "bv_normal", "bv_matern_normal", "bv_matern_nig")))
    stop("Not allow for share_V for univariate model")

  if (model %in% c("re", "bv_normal", "bv_matern_normal", "bv_matern_nig")) {
    noise$fix_theta_sigma <- TRUE
    noise$n_params  <- noise$n_params - noise$n_theta_sigma
    noise$n_theta_sigma  <- 0

    # for printing to ignore sigma
    if (model != "re") {
      noise$bv_noises[[1]]$fix_theta_sigma <- TRUE
      noise$bv_noises[[2]]$fix_theta_sigma <- TRUE
    }
  }

  if (model %in% c("rw1", "rw2")) {
    # force the constrain (e.g. fixed the 1st position to be 0)
    theta_sigma <- c(-10, noise$theta_sigma)
    fix_theta_sigma <- c(TRUE, noise$fix_theta_sigma)
    B_sigma <- cbind(0, noise$B_sigma); 
    B_sigma[1, 1] <- 1; noise$B_sigma[1, -1] <- 0
    if (any(abs(eigen(t(B_sigma) %*% B_sigma)$values) < 1e-10)) {
      # cannot apply transformation
      stop("Do not fix the 1st position of sigma for non-cyclic rw1")
    } else {
      noise$B_sigma <- B_sigma
      noise$theta_sigma <- theta_sigma
      noise$fix_theta_sigma <- fix_theta_sigma
    }

    if (model == "rw2" && !operator$cyclic) {
      noise$B_sigma[2, 1] <- 1; noise$B_sigma[2, -1] <- 0
      noise$theta_sigma[1] <- -6  # set larger for numerical stability
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

  ngme_model(
    model     = model,
    operator  = operator,
    noise     = noise,
    W_size    = ncol(operator$K),
    V_size    = nrow(operator$K),
    A         = A,
    map       = map,
    mesh      = mesh,
    mesh_list = mesh_list,
    n_map     = length_map(map),
    W         = W,
    fix_W      = fix_W,
    name      = name,
    debug     = debug,
    fix_theta_K = fix_theta_K,
    prior_theta_K = prior_theta_K,
    replicate = replicate
  )
}

# build operator
build_operator <- function(model_name, args_list) {
  stopifnot(
    is.character(model_name),
    is.list(args_list)
  )

  switch(model_name,
    tp  = do.call(tp, args_list),
    bv  = do.call(bv, args_list),
    bv_normal = do.call(bv_normal, args_list),
    bv_matern_normal = do.call(bv_matern_normal, args_list),
    bv_matern_nig = do.call(bv_matern_nig, args_list),
    ar  = do.call(ar, args_list),
    arma = do.call(arma, args_list),
    ar1 = do.call(ar1, args_list),
    rw1 = do.call(rw1, args_list),
    rw2 = do.call(rw2, args_list),
    ou  = do.call(ou, args_list),
    matern = do.call(matern, args_list),
    re  = do.call(re, args_list),
    spacetime = do.call(spacetime, args_list),
    generic = do.call(generic, args_list),
    generic_ns = do.call(generic_ns, args_list),
    iid = do.call(iid, args_list),
    stop("Unknown models, please check if model name is in ngme_model_types()")
  )
}

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
    nu = exp(nig_noise$theta_nu),
    B_mu = nig_noise$B_mu,
    B_sigma_nig = nig_noise$B_sigma,
    B_sigma_normal = normal_noise$B_sigma,
    B_nu = nig_noise$B_nu,
    fix_theta_mu = nig_noise$fix_theta_mu,
    fix_theta_sigma = c(nig_noise$fix_theta_sigma,
                        normal_noise$fix_theta_sigma),
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
  ...
) {
  if (inherits(loc, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) return(loc)

  if (!is.null(model)) {
    if (model %in% c("re", "tp")) return(NULL)
    if (model == "iid") loc <- as.integer(as.factor(loc))
    if (model %in% c("ar", "ar1", "arma")) {
      stopifnot("The map should be integers."
        = is.numeric(loc) && all(loc == round(loc)))
      return (fmesher::fm_mesh_1d(loc = min(loc):max(loc)))
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
