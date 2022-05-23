# inla
library(INLA)
data("SPDEtoy")

SPDEtoy.sp <- SPDEtoy
summary(SPDEtoy)
coordinates(SPDEtoy.sp) <- ~ s1 + s2

?bubble

# showMethods(coordinates)
bubble(SPDEtoy.sp, "y", key.entries = c(5, 7.5, 10, 12.5, 15),
       maxsize = 2, xlab = "s1", ylab = "s2")

m0 <- inla(y ~ s1 + s2, data = SPDEtoy)
summary(m0)

# random walk 1
f.rw1 <- y ~ f(s1, model = "rw1", scale.model = TRUE) +
  f(s2, model = "rw1", scale.model = TRUE)

m1 <- inla(f.rw1, data = SPDEtoy)
summary(m1)

####
f(s1, model = "rw1", scale.model = TRUE)

#

# y ~ f(1:n, model="AR1", noise="NIG")

library(INLA)
vignette('SPDEhowto')
vignette('SPDE1d')

#
attributes(a_init)
typeof(SPDEtoy)
class(SPDEtoy)


