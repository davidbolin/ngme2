library(devtools)
load_all()
usethis::use_pkgdown()

pkgdown::build_site(devel = TRUE)

# build website


library(INLA)
mesh1 <- inla.mesh.1d(rnorm(100))
mesh1

inla.spde.make.A(mesh=mesh1, loc = 4)

t.data.frame
predict.glm

?scale

?print.lm
?lm

y <- 1:100
x <- rnorm(100)
m1 <- (lm(y~x))
??print.lm
print.summary.lm

??print.ngme
?print.lm
??ngme

