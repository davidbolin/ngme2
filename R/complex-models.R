

#' Ngme bivariate model specification
#'
#' Giving 2 sub_models, build a correlated bivaraite operator based on K = D(theta, eta) %*% diag(K_1, K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (total 2 sub_models)
#' @param mesh mesh for build the model
#' @param group group vector, can be inherited from ngme() function
#' @param theta the rotation parameter (starting point should
#'   from -pi/4 to pi/4)
#' @param rho the parameter related to correlation
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
bv <- function(
  mesh,
  sub_models,
  theta = 0,
  rho = 0,
  group = NULL,
  share_param = FALSE,
  fix_bv_theta = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  model_names <- sort(names(sub_models))

  # read group argument from parent frame
  if (is.null(group) &&
      exists("group", parent.frame())){
    group = get("group", parent.frame())
  }

  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields"
      = !is.null(group),
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group"
    = all(model_names %in% levels(as.factor(group)))
  )
  group <- factor(group, levels = model_names)

  # build 2 sub_models (pass environment to sub_models)
  # delete theta and rho in env_args (to avoid sub_models use them)
  env_args <- as.list(environment())
  env_args$theta <- NULL; env_args$rho <- NULL
  arg1 <- sub_models[[model_names[1]]]
  if (!is.list(arg1)) {
    first <- build_operator(arg1, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg1$model))
    first <- build_operator(arg1$model, modifyList(env_args, arg1))
  }
  arg2 <- sub_models[[model_names[2]]]
  if (!is.list(arg2)) {
    second <- build_operator(arg2, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg2$model))
    second <- build_operator(arg2$model, modifyList(env_args, arg2))
  }

  stopifnot(
    "the theta should be in (-pi/2, pi/2)"
      = theta >= -pi/2 & theta <= pi/2
  )
  eta <- tan(theta)
  theta_K <- c(eta, rho, first$theta_K, second$theta_K)

  # pass the theta_K to first and second
  first$theta_K <- theta_K[3:(2 + first$n_theta_K)]
  second$theta_K <- theta_K[(3 + first$n_theta_K):length(theta_K)]

  D <- build_D(theta, rho)
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

  # the range of theta is (-pi/4, pi/4), we expand it to (-pi/2, pi/2)
  # so that the maximum can be reached in either direction
  # in the end, we will NOT transform it back to (-pi/4, pi/4)
  eta_to_theta <- function(eta) atan(eta)

  update_K <- function(theta_K) {
    eta <- theta_K[1]
    rho <- theta_K[2]
    theta <- atan(eta)
    D <- build_D(theta, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))
    theta_K1 <- theta_K[3:(2 + first$n_theta_K)]
    theta_K2 <- theta_K[(3 + first$n_theta_K):length(theta_K)]
    first$K <- first$update_K(theta_K1)
    second$K <- second$update_K(theta_K2)
    bigD %*% Matrix::bdiag(first$K, second$K)
  }

  ngme_operator(
    mesh        = mesh,
    model       = "bv",
    first        = first,
    second      = second,
    theta_K     = theta_K,
    update_K    = update_K,
    K           = bigD %*% Matrix::bdiag(first$K, second$K),
    h           = c(first$h, second$h),
    symmetric   = FALSE,
    zero_trace  = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_bv_theta = fix_bv_theta,
    param_name  = c("theta", "rho",
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(list(eta_to_theta, identity), first$param_trans, second$param_trans)
  )
}


#' ngme tensor-product model specification
#'
#' Given 2 operator (first and second), build a tensor-product operator based on K = K_first x K_second (here x is Kronecker product)
#'
#' @param first ngme_operator, left side of kronecker model (usually a temporal or iid model)
#' @param second ngme_operator, right side of kronecker model (ususally a temporal or spatial model)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
tp <- function(
  first,
  second,
  ...
) {
  stopifnot(
    inherits(first, "ngme_operator"),
    inherits(second, "ngme_operator")
  )

  theta_K <- c(first$theta_K, second$theta_K)
  stopifnot(length(theta_K) == first$n_theta_K + second$n_theta_K)

  update_K <- function(theta_K) {
    first$K %x% second$K
  }
  ngme_operator(
    mesh = list(first$mesh, second$mesh),
    model = "tp",
    first = first,
    second = second,
    theta_K = theta_K,
    update_K = update_K,
    K = first$K %x% second$K,
    h = first$h %x% second$h,
    symmetric = first$symmetric & second$symmetric,
    zero_trace = first$zero_trace & second$zero_trace,
    param_name = c(first$param_name, second$param_name),
    param_trans = c(first$param_trans, second$param_trans)
  )
}


#' Ngme bivariate model specification on operator (theta=0)
#'
#' Giving 2 sub_models, build a correlated bivaraite operator based on K = D(theta, eta) %*% diag(K_1, K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#' Exact same implementation as in paper.
#' Fix noise sd=1.
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (total 2 sub_models)
#' @param mesh mesh for build the model
#' @param group group vector, can be inherited from ngme() function
#' @param rho the parameter related to correlation
#' @param c1 the noise sd for 1st sub_model
#' @param c2 the noise sd for 2nd sub_model
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
bv_normal <- function(
  mesh,
  sub_models,
  rho = 0,
  c1 = 1, 
  c2 = 1,
  group = NULL,
  share_param = FALSE,
  fix_bv_theta = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  model_names <- sort(names(sub_models))

  # read group argument from parent frame
  if (is.null(group) &&
      exists("group", parent.frame())){
    group = get("group", parent.frame())
  }

  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields"
      = !is.null(group),
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group"
    = all(model_names %in% levels(as.factor(group)))
  )
  group <- factor(group, levels = model_names)

  # build 2 sub_models (pass environment to sub_models)
  # delete theta and rho in env_args (to avoid sub_models use them)
  env_args <- as.list(environment())
  env_args$theta <- NULL;
  env_args$rho <- NULL
  arg1 <- sub_models[[model_names[1]]]
  if (!is.list(arg1)) {
    first <- build_operator(arg1, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg1$model))
    first <- build_operator(arg1$model, modifyList(env_args, arg1))
  }
  arg2 <- sub_models[[model_names[2]]]
  if (!is.list(arg2)) {
    second <- build_operator(arg2, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg2$model))
    second <- build_operator(arg2$model, modifyList(env_args, arg2))
  }

  stopifnot(c1 > 0, c2 > 0)
  theta_K <- c(rho, log(c1), log(c2), first$theta_K, second$theta_K)

  update_K <- function(theta_K) {
    rho <- theta_K[1]
    c1 <- exp(theta_K[2])
    c2 <- exp(theta_K[3])

    D <- build_D(0, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

    # pass the theta_K to first and second
    first$theta_K <- theta_K[4:(3 + first$n_theta_K)]
    second$theta_K <- theta_K[(4 + first$n_theta_K):length(theta_K)]

    first$K <- first$update_K(first$theta_K)
    second$K <- second$update_K(second$theta_K)

    bigD %*% Matrix::bdiag(c1 * first$K, c2 * second$K)
  }

  ngme_operator(
    mesh        = mesh,
    model       = "bv_normal",
    first        = first,
    second      = second,
    theta_K     = theta_K,
    update_K    = update_K,
    K           = update_K(theta_K),
    h           = c(first$h, second$h),
    symmetric   = FALSE,
    zero_trace  = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_bv_theta = fix_bv_theta,
    param_name  = c("rho",
      "c1",
      "c2",
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(identity, exp, exp, first$param_trans, second$param_trans)
  )
}


#' Ngme bivariate model specification (theta=0)
#'
#' Does not fit noise sd.
#' 
#' Giving 2 sub_models, build a correlated bivaraite operator based on K = D(theta, eta) %*% diag(K_1, K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (should be two matern models)
#' @param mesh mesh for build the model
#' @param group group vector, can be inherited from ngme() function
#' @param rho the parameter related to correlation
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
bv_matern_normal <- function(
  mesh,
  sub_models,
  rho = 0,
  sd1=1, sd2=1,
  group = NULL,
  share_param = FALSE,
  fix_bv_theta = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  d <- switch(mesh$manifold, 
    "R1" = 1,
    "R2" = 2,
    stop("manifold not supported")
  )

  # sort the sub_models
  model_names <- sort(names(sub_models))

  # read group argument from parent frame
  if (is.null(group) &&
      exists("group", parent.frame())){
    group = get("group", parent.frame())
  }

  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields"
      = !is.null(group),
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group"
    = all(model_names %in% levels(as.factor(group))),
    "sd1 and sd2 should be positive" = sd1 > 0 & sd2 > 0
  )
  group <- factor(group, levels = model_names)

  # build 2 sub_models (pass environment to sub_models)
  # delete theta and rho in env_args (to avoid sub_models use them)
  env_args <- as.list(environment())
  env_args$theta <- NULL;
  env_args$rho <- NULL
  arg1 <- sub_models[[model_names[1]]]
  if (!is.list(arg1)) {
    first <- build_operator(arg1, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg1$model))
    first <- build_operator(arg1$model, modifyList(env_args, arg1))
  }
  arg2 <- sub_models[[model_names[2]]]
  if (!is.list(arg2)) {
    second <- build_operator(arg2, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg2$model))
    second <- build_operator(arg2$model, modifyList(env_args, arg2))
  }

  theta_K <- c(
    rho, log(sd1), log(sd2),
    first$theta_K, second$theta_K
  )

  update_K <- function(theta_K) {
    theta <- 0
    rho <- theta_K[1]
    sd1 <- exp(theta_K[2])
    sd2 <- exp(theta_K[3])
    theta_K1 <- theta_K[4:(3 + first$n_theta_K)]
    theta_K2 <- theta_K[(4 + first$n_theta_K):length(theta_K)]
    
    alpha1 <- first$alpha; alpha2 <- second$alpha
    kappa1 <- exp(first$theta_K); kappa2 <- exp(second$theta_K)

    nu1 = alpha1 - d/2; nu2 = alpha2 - d/2
    c1 = sqrt(gamma(nu1) / (gamma(alpha1))) / (kappa1^(nu1) * (4*pi)^(d/4) * sd1)
    c2 = sqrt(gamma(nu2) / (gamma(alpha2))) / (kappa2^(nu2) * (4*pi)^(d/4) * sd2)

    D <- build_D(theta, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))
    first$K <- first$update_K(theta_K1)
    second$K <- second$update_K(theta_K2)
    bigD %*% Matrix::bdiag(c1 * first$K, c2 * second$K)
  }

  ngme_operator(
    mesh        = mesh,
    model       = "bv_matern_normal",
    first        = first,
    second      = second,
    theta_K     = theta_K,
    update_K    = update_K,
    K           = update_K(theta_K),
    h           = c(first$h, second$h),
    dim         = d,
    symmetric   = FALSE,
    zero_trace  = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_bv_theta = fix_bv_theta,
    param_name  = c("rho", "sd1", "sd2",
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(identity, exp, exp, first$param_trans, second$param_trans)
  )
}



#' ngme space-time non-separable model specification
#'
#' Given a spatial and temporal model, build a non-separable space-time model
#'
#'   can use ~ map_s + map_t to specify the model
#' @param mesh a list of 2 objects,
#' .   list(mesh_t, mesh_s), mesh_t is temporal mesh, mesh_s is spatial mesh
# ' @param alpha 2 or 4, SPDE smoothness parameter
#' @param c parameter
#' @param kappa kappa parameter from matern SPDE
#' @param lambda the spatial damping parameter
#' @param gamma  2d vector, direction of the transport term
#' @param method discretization method,
#' @param alpha 2 or 4, SPDE smoothness parameter
#' choose "galerkin" or "euler" for implicit euler
#' @param normalize_gamma TRUE if normalize gamma to unit length
#' @param stabilization TRUE if use stabilization term (for implicit euler)
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
spacetime <- function(
  mesh,
  gamma = c(1, 1), # gamma (s): function of mesh node
  lambda = 1, # fixed
  alpha = 2, # alpha = 2, 4, fixed
  method = "galerkin", # galerkin, implicit euler
  # parameters
  c = 1,
  kappa = 1,
  normalize_gamma = TRUE,
  stabilization = TRUE,
  ...
) {

  stopifnot(
    "Please provide mesh as a list of length 2" = length(mesh) == 2,
    "alpha should be 2 or 4" = alpha == 2 || alpha == 4,
    "gamma should be of length 2" = length(gamma) == 2,
    "method should be galerkin or euler" = method %in% c("galerkin", "euler"),
    "First mesh should be 1d" = fmesher::fm_manifold_dim(mesh[[1]]) == 1,
    "Second mesh should be 2d" =
      fmesher::fm_manifold_dim(mesh[[2]]) == 2
  )

  if (normalize_gamma) gamma <- gamma / sqrt(sum(gamma^2))

  mesh_t <- mesh[[1]]
  mesh_s <- mesh[[2]]

##### temporal FEM #####
  fem_t <- fmesher::fm_fem(mesh_t, order = 2)
  nt <- mesh_t$n
  Ct <- fem_t$c0
  Gt <- fem_t$g1
  
  # Build Bt
  Bt <- Matrix::bandSparse(n = nt, m = nt, k = c(-1, 0, 1),
                   diagonals = cbind(rep(0.5,nt), rep(0,nt), rep(-0.5,nt)))
  Bt[1,1] = -0.5
  Bt[nt,nt] = 0.5

##### space FEM #####
  fem_s <- fmesher::fm_fem(mesh_s, order = alpha)
  ns <- mesh_s$n
  Cs <- fem_s$c0
  Gs <- fem_s$g1
  
  # Build Bs
  ret = ngme_make_Bs(mesh_s, gamma)
  Bs <- ret$H
  S <- ret$S  #stability matrix
  
  # compute h
  if (method == "galerkin"){
    h <- Matrix::diag(Ct) %x% Matrix::diag(Cs)
  } else {
    dt =c(1, diff(mesh_t$loc))
    h <- dt %x% Matrix::diag(Cs) / c
  }

  # FV = mesh$graph$tv
  # P <- sf::st_coordinates(fmesher::fm_vertices(mesh))[,1:2]
  # fem2d <- rSPDE.fem2d(FV = FV, P = P)

  theta_K <- c(log(c), log(kappa))
  update_K <- function(theta_K) {
    c <- exp(theta_K[1])
    kappa <- exp(theta_K[2])

    # compute L_s
    L = kappa^2 * Cs + lambda * Gs + Bs
  
    Cs_inv <- Matrix::Diagonal(n = ns, x = 1/Matrix::diag(Cs))
    
  # watch-out! make sure which t
    if (alpha == 4) L = L %*% Cs_inv %*% Matrix::t(L)
    
    if (method == "galerkin") {
      K <- Bt %x% Cs + 1/c * Ct %x% L
    } else if (method == "euler") {
      # implicit euler
      dt =c(1, diff(mesh_t$loc))

      null_matrix <- Matrix::Diagonal(n = ns, x = 0)
      
      if (stabilization) {
        L = L + S
      }
      
      diag_L <- Matrix::bdiag(
        lapply(1:(nt-1), function(i) L)
      )
      
      K <- rw1(1:nt)$K %x% Cs + 1/c *
        Matrix::bdiag(null_matrix, diag_L)
    }
    return (K)
  }
  K <- update_K(theta_K)

  BtCs = if (method == "galerkin")
    ngme_as_sparse(Bt %x% Cs) 
    else ngme_as_sparse(rw1(1:nt)$K %x% Cs )

  ngme_operator(
    model = "spacetime",
    mesh = mesh,
    alpha = alpha,
    gamma = gamma,
    lambda = lambda,
    method = method,
    theta_K = theta_K,
    K = ngme_as_sparse(K),
    update_K = update_K,
    h = h,
    Ct_diag = Matrix::diag(Ct),
    Cs_diag = Matrix::diag(Cs),
    BtCs = BtCs,
    Ct = ngme_as_sparse(Ct),
    Cs = ngme_as_sparse(Cs),
    Gs = ngme_as_sparse(Gs),
    Bs = ngme_as_sparse(Bs),
    S = ngme_as_sparse(S),
    stabilization = stabilization,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name =
      c("c", "kappa", "lambda"),
    param_trans =
      c(identity, exp, identity)
  )
}



