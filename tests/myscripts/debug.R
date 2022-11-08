library(devtools)
usethis::use_pkgdown()

pkgdown::build_site(devel = TRUE)

# build website


library(INLA)
mesh1 <- inla.mesh.1d(rnorm(100))
mesh1

inla.spde.make.A(mesh=mesh1, loc = 4)
