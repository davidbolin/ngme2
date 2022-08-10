library(INLA)
dat <- read.table('myscripts/Argo/ocean_data.txt',header=TRUE)
str(dat)

max.edge    = 1
bound.outer = 5
locY = unique(cbind(dat$lon, dat$lat))
# nrow(locY) == nrow(dat) no replicates

mesh = inla.mesh.2d(loc=locY,
                    # the inner edge and outer edge
                    max.edge = c(1,5) * max.edge,
                    # offset extension distance inner and outer extenstion
                    offset = c(max.edge, bound.outer)
)
plot(mesh)
points(locY, col = "red")
fem_mesh = inla.mesh.fem(mesh)
Ce <- fem_mesh$c1 #<phi_i, phi_j>
C <- fem_mesh$c0 #approximation of Ce
G <- fem_mesh$g1
A <- inla.spde.make.A(mesh, locY) #dim(A) = data loc * vertices

#There are some columns in the projector matrix all of whose elements equal zero:
table(colSums(A) > 0)

# #These columns correspond to triangles with no point location inside.
# #These columns can be dropped.
# P <-mesh$loc
# FV<-mesh$graph$tv


# ------- ngme fit --------
library(devtools)
load_all()

formula1 <- Y ~ 0 + f(
  1:mesh$n,
  model = matern <- ngme.matern(
    alpha = 2,
    mesh = mesh,
    kappa = 1
  ),
  A = A,
  debug = TRUE,
  theta.mu = 0,
  theta.sigma = 0,
  noise = "nig",
  control = ngme.control.f(
    numer_grad       = FALSE,
    use_precond      = TRUE,

    fix_operator     = FALSE,
    fix_mu           = FALSE,
    fix_sigma        = FALSE,
    fix_noise        = FALSE,
    use_iter_solver  = FALSE
  )
)

ngme_out1 <- ngme(
  formula = formula1,
  data = data.frame(
    Y = dat$temp
  ),
  family = "normal",
  control=ngme.control(
    burnin=100,
    iterations=100,
    gibbs_sample = 5
  ),
  debug=ngme.debug(
    debug = TRUE,
    not_run = F
  )
)

ngme_out1$result

# --------- using non-stationary mu ------------

formula2 <- Y ~ 0 + f(
  1:mesh$n,
  model = matern <- ngme.matern(
    alpha = 2,
    mesh = mesh,
    kappa = 1
  ),
  A = A,
  debug = TRUE,
  theta.mu = c(0.5, 0.4),
  theta.sigma = 0,
  noise = "nig",
  control = ngme.control.f(
    numer_grad       = FALSE,
    use_precond      = TRUE,

    fix_operator     = FALSE,
    fix_mu           = FALSE,
    fix_sigma        = FALSE,
    fix_noise        = FALSE,
    use_iter_solver  = FALSE
  )
)

ngme_out2 <- ngme(
  formula = formula2,
  data = data.frame(
    Y = dat$temp
  ),
  family = "normal",
  control=ngme.control(
    burnin=100,
    iterations=300,
    gibbs_sample = 5
  ),
  debug=ngme.debug(
    debug = TRUE,
    not_run = F
  )
#   start = ngme_out1
)

ngme_out2$result
