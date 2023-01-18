# test inla functionality

library(INLA)

data(SPDEtoy)
toy_sp <- SPDEtoy
coordinates(toy_sp) <- ~ s1 + s2

m <- inla(
  y ~  s1 + f(s2, model="rw1", scale.model=TRUE),
  data=SPDEtoy
)
summary(m)


# ngme
load_all()
ngme(y ~ s1 + f(s2, model="rw1"), data=SPDEtoy)

