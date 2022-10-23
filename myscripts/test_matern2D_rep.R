# Assuming that two variables are observed at same locations with same number of replicates
library(INLA)
library(MASS)
set.seed(41)


# Define locations
# n = 10 points on a [0,1]*[0,1] for each process
nrep <- 3 # 5 years observed at the same location(inside box)
n <- round(runif(nrep, 7, 10)) # length of each replicate varies
n_total <- sum(n) # total number of observations
loc_2d_mesh <- list()
for (i in 1:nrep) {
    loc_2d_mesh[[i]] <- 5 * matrix(runif(n[i] * 2), n[i], 2)
    plot(loc_2d_mesh[[i]])
}
# choose index of a replicate (year in context of Argo) with max number of observations to create the mesh
imax <- which.max(n)

# Create mesh
mesh_2d <- inla.mesh.2d(
    loc = loc_2d_mesh[[imax]],
    cutoff = 0.1,
    # the inner edge and outer edge
    max.edge = c(1, 5),
    # offset extension distance inner and outer extension
    offset = c(0.5, 1)
)
plot(mesh_2d, draw.vertices = TRUE, main = "Common mesh for 2 processes")

# # For process_1 & process_2 in case of equal number of observations per replicates
# fem_1d.2 <- fem_1d.1 <- inla.mesh.fem(mesh_2d)
# C2 <- C1 <- fem_1d.1$c1 # Creates the mass matrix
# G2 <- G1 <- fem_1d.1$g1 # Creates the stiffness matrix
# h2 <- h1 <- diag(fem_1d.1$c0) # Creates the "volume" vector
# m <- length(h1) # number of basis functions
# A = inla.spde.make.A(mesh_2d, loc_2d_mesh)
# A = bdiag(A, A)


# Parameters
kappa1 <- 2
kappa2 <- 3
sigma1 <- 2
sigma2 <- 3
sigma_eps1 <- 0.2
sigma_eps2 <- 0.5
beta1 <- 10
beta2 <- 20
alpha1 <- 2
alpha2 <- 2
d <- 2
sigma <- c(sigma1, sigma2)
alpha <- c(alpha1, alpha2)
kappa <- c(kappa1, kappa2)
rho <- 0.1
theta <- atan(rho)
rho_eps <- 0.8
# rescaling constant c_i
tau <- sqrt(gamma(alpha - d / 2) / (gamma(alpha) * (4 * pi)^(d / 2) * kappa^(2 * (alpha - d / 2)) * sigma^2))
c1 <- tau[1]
c2 <- tau[2]
# additional parameters for NIG
nu1 <- 3 # eta
nu2 <- 7
mu1 <- 1
mu2 <- 2


# parametrizations to separate control of variances & cross-correlations
D <- matrix(c(
    cos(theta) + rho * sin(theta),
    sin(theta) - rho * cos(theta),
    -sin(theta) * sqrt(1 + rho^2),
    cos(theta) * sqrt(1 + rho^2)
), 2)
# covaraince of measurement noise
cov_eps <- matrix(c(
    sigma_eps1^2, sigma_eps1 * sigma_eps2 * rho_eps,
    sigma_eps1 * sigma_eps2 * rho_eps, sigma_eps2^2
), 2)

# Projection matrix
A <- inla.spde.make.A(
    mesh = mesh_2d,
    loc = loc_2d_mesh[[1]],
    index = 1:n[1],
    repl = rep(1, n[1])
)

for (j in 2:nrep) {
    Atmp <- inla.spde.make.A(
        mesh = mesh_2d,
        loc = loc_2d_mesh[[j]],
        index = 1:n[j],
        repl = rep(1, n[j])
    )
    A <- bdiag(A, Atmp)
}

image(A)
A <- bdiag(A, A) #for bivariate field
fem_mesh <- inla.mesh.fem(mesh_2d)
C <- kronecker(Diagonal(nrep, 1), fem_mesh$c0)
G <- kronecker(Diagonal(nrep, 1), fem_mesh$g1)
h <- diag(fem_mesh$c0)
m <- length(h)
image(C)
image(G)
K1 <- kappa1^2 * C + G
K2 <- kappa2^2 * C + G
K <- bdiag(c1 * K1, c2 * K2)
K_D <- (kronecker(D, diag(nrep * m))) %*% K
image(K_D)
image(K_D)

Beta <- c(rep(beta1, n_total), rep(beta2, n_total))
# measurement noise
noise <- mvrnorm(n_total, mu = rep(0, 2), Sigma = cov_eps)
noise <- c(noise[, 1], noise[, 2])
V1 <- ngme2::rig(m, a = rep(nu1, m), b = h^2 * nu1, sample.int(10^6, 1))
V2 <- ngme2::rig(m, a = rep(nu2, m), b = h^2 * nu2, sample.int(10^6, 1))
f.V1 <- mvrnorm(1, mu = mu1 * (-h + V1), Sigma = diag(V1))
f.V2 <- mvrnorm(1, mu = mu2 * (-h + V2), Sigma = diag(V2))
temp <- c(f.V1, f.V2) # joining the vector of M_dot - NIG noise process
temp <- rep(temp, nrep)
X <- solve(K_D, temp) # sample of Gaussian Process
X1 <- X[1:n_total, ]
X2 <- X[(n_total + 1):(2 * n_total), ]

Y <- Beta + A %*% X + noise
Y1 = Y[1:n_total, ]
Y2 = Y[(n_total + 1):(2 * n_total), ]



