# From David's Nested SPDE code
make.nested.spde <- function(mesh, spde, B.vec, prior.vec.mu, prior.vec.Q) {
  n.vec <- dim(B.vec[[1]])[1]

  if (missing(prior.vec.mu)) prior.vec.mu <- rep(0, n.vec)

  if (missing(prior.vec.Q)) {
    prior.vec.Q <- Matrix::sparseMatrix(
      i = 1:n.vec, j = 1:n.vec, x = rep(0.01, n.vec),
      dims = c(n.vec, n.vec)
    )
  }
  # ?inla.fmesher.smorg
  # fem2 = inla.fmesher.smorg(mesh$loc, mesh$graph$tv, fem = 2, gradients=TRUE)

  fem <- fmesher::fm_basis(mesh, loc = mesh$loc, derivatives = TRUE)

  n <- dim(B.vec[[1]])[2]
  
  Bvec <- list()
  for (i in 1:n.vec) {
    Bvec[[i]] <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = B.vec[[1]][i, ], dims = c(n, n)) %*% fem$dx +
      Matrix::sparseMatrix(i = 1:n, j = 1:n, x = B.vec[[2]][i, ], dims = c(n, n)) %*% fem$dy +
      Matrix::sparseMatrix(i = 1:n, j = 1:n, x = B.vec[[3]][i, ], dims = c(n, n)) %*% fem$dz
  }

  prior.vec <- list(mu = prior.vec.mu, Q = prior.vec.Q)
  nested.spde <- list(
    mesh = mesh,
    spde = spde,
    model = "nested spde",
    B.vec = Bvec,
    n.vec = n.vec,
    prior.vec = prior.vec
  )
  return(nested.spde)
}


nested.spde.H <- function(nested.spde, theta.vec) {
  
  n <- nested.spde$mesh$n 
  # browser()
  # n <- dim(nested.spde$spde$param.inla$M0)[1]
  
  Hl <- 0
  if (!is.na(theta.vec[1])) {
    for (i in 1:length(theta.vec)) {
      Hl <- Hl + theta.vec[i] * nested.spde$B.vec[[i]]
    }
  }
  H <- Hl + Matrix::sparseMatrix(i = 1:n, j = 1:n, x = rep(1, n), dims = c(n, n))
  return(H)
}

nested.spde.S <- function(nested.spde, theta.vec) {
  n <- nested.spde$mesh$n 
  # n <- dim(nested.spde$spde$param.inla$M0)[1]
  
  Hl <- 0
  if (!is.na(theta.vec[1])) {
    for (i in 1:length(theta.vec)) {
      Hl <- Hl + theta.vec[i] * nested.spde$B.vec[[i]]
    }
  }
  H <- Hl + Matrix::sparseMatrix(i = 1:n, j = 1:n, x = rep(1, n), dims = c(n, n))
  
  return(H %*% Matrix::t(H))

  # gamma1^2 Hxx + gamma2^2 Hyy + gamma1 * gamma2 (Hxy + Hyx)
}


vector.basis.2d <- function(mesh, FVc = FALSE) {
  if (FVc) {
    n <- dim(mesh$graph$tv)[1]
    FV <- mesh$graph$tv
    nF <- dim(FV)[1]
    Ptri <- sapply(1:nF, function(f) rowMeans(P[, FV[f, ]]))
    P <- sapply(1:nF, function(f) Ptri[, f] / sqrt(t(Ptri[, f]) %*% Ptri[, f]))
  } else {
    n <- dim(mesh$loc)[1]
    P <- mesh$loc
  }
  B <- list()
  for (i in 1:3) {
    B[[i]] <- Matrix::Matrix(0, 3, n)
    B[[i]][i, ] <- 1
  }
  return(list(B = B, P = P))
}

# make Bs given a 2d mesh on R2
ngme_make_Bs <- function(mesh, gamma.vec) {
  stopifnot(
    "2d mesh is required." = fmesher::fm_manifold_dim(mesh) == 2,
    length(gamma.vec) == 2
  )
  gamma.vec <- c(gamma.vec)

  v.b = vector.basis.2d(mesh=mesh)
  
  # spde <- inla.spde2.matern(mesh, alpha=2)
  spde <- NULL # not needed
  
  n.spde = make.nested.spde(mesh, spde, v.b$B)
  H = nested.spde.H(n.spde, gamma.vec)
  S = nested.spde.S(n.spde, gamma.vec)
  
  list(
    H = H,
    S = S
  )
}