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

?use_import_from()
library(devtools)
use_import_from("rlang", ".data")

use_vignette("Argo_flat_data")

library(grid)
grid.arrange(rectGrob(), rectGrob())
## Not run:
library(ggplot2)
pl <- lapply(1:3, function(.x) qplot(1:10, rnorm(10), main=paste("plot", .x)))
ml <- marrangeGrob(pl, nrow=2, ncol=2)
length(pl)
## non-interactive use, multipage pdf
ggsave("multipage.pdf", ml)
## interactive use; open new devices
ml

str(out$noise)
?switch
pl
c(list(a=1, b=2), list(a=2))

str(out$latents[[1]][["theta_mu"]])
