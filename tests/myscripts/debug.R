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

library(devtools)
use_vignette("model-estimation")

devtools::load_all()
m = model_rw(c(1.1, 2.2, 3.5, 5.6), order=1); m$C + m$G

m = model_ar1(c(1, 2, 3, 10), order=1, noise=noise_nig()); m$C + m$G

str(m)
check()


# sim rw1 model
devtools::load_all()
x <- sort(rexp(10))
model_rw(x)$noise$n_noise

solve(matrix(rexp(6), nrow=2), b=c(1,2))

ff <- function(a=1, ...) {
  print(list(...))
}
ff(b=3)

# double dlambda_V(const double loglambda,
#                  const Eigen::VectorXd &V,
#                  const Eigen::VectorXd &h,
#                  const int GAL)
# {
#   double dlambda = 0;
#   if(GAL){
#   for(int i=0; i < h.size(); i++){
#     double h_lambda = exp(loglambda) * h[i];
#     //digamma(0.1) = digamma(1.1) - 1/0.1;
#     if(h_lambda > 1){
#         dlambda -=  h_lambda * R::digamma(h_lambda);
#       }else
#       {
#         dlambda -=  h_lambda * R::digamma(h_lambda + 1) - 1.;
#       }
#     dlambda += h_lambda *  log(V(i) ) ;
#   }
#   }else{
#     double srqt_two = pow(2, 0.5);
#     for(int i=0; i < h.size(); i++){
#       dlambda +=  1 -  ( pow(h(i), 2) / V(i) ) * exp( 2 * loglambda);
#       dlambda += srqt_two * h(i)  * exp(loglambda);
#     }

#   }

#   return(dlambda);
# }


# to-do
# replicates
m3 <- model_rw(c(1.1, 2.2, 2.2, 3.3, 4.4), replicates = c(1, 1, 2, 1, 1))
m3$K
m3$A
#  [1,] 1 . . . . . . .
#  [2,] . 1 . . . . . .
#  [3,] . . . . . 1 . .
#  [4,] . . 1 . . . . .
#  [5,] . . . 1 . . . .

library(INLA)
mesh = inla.mesh.1d(loc=c(1.1, 2.2, 2.2, 3.3, 4.4))
inla.mesh.1d.A(mesh=mesh, loc=c(1,2,1))
..$ :List of 21
  .. ..$ model      : chr "rw1"
  .. ..$ noise_type : chr "normal"
  .. ..$ W_size     : int 25
  .. ..$ theta_K    : num 1
  .. ..$ n_theta_K  : int 1
  .. ..$ A          :Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
  .. .. .. ..@ i       : int [1:1194] 291 293 295 298 302 304 306 307 308 560 ...
  .. .. .. ..@ j       : int [1:1194] 0 0 0 0 0 0 0 0 0 0 ...
  .. .. .. ..@ Dim     : int [1:2] 597 25
  .. .. .. ..@ Dimnames:List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : NULL
  .. .. .. ..@ x       : num [1:1194] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..@ factors : list()
  .. ..$ A_pred     : NULL
  .. ..$ noise      :List of 17
  .. .. ..$ n_noise        : int 24
  .. .. ..$ h              : num [1:24] 24 28.5 23.2 26.2 15.8 ...
  .. .. ..$ noise_type     : chr "normal"
  .. .. ..$ theta_V        : num 1
  .. .. ..$ V              : NULL
  .. .. ..$ theta_mu       : num 0
  .. .. ..$ theta_sigma    : num 0
  .. .. ..$ B_mu           : num [1:24, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..$ B_sigma        : num [1:24, 1] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..$ n_theta_mu     : int 1
  .. .. ..$ n_theta_sigma  : int 1
  .. .. ..$ n_theta_V      : num 1
  .. .. ..$ fix_theta_mu   : logi FALSE
  .. .. ..$ fix_theta_sigma: logi FALSE
  .. .. ..$ fix_theta_V    : logi FALSE
  .. .. ..$ fix_V          : logi FALSE
  .. .. ..$ n_params       : int 1
  .. .. ..- attr(*, "class")= chr "ngme_noise"
  .. ..$ W          : NULL
  .. ..$ fix_W      : logi FALSE
  .. ..$ fix_theta_K: logi TRUE
  .. ..$ V_size     : num 24
  .. ..$ control    :List of 5
  .. .. ..$ numer_grad     : logi FALSE
  .. .. ..$ use_precond    : logi FALSE
  .. .. ..$ use_num_hess   : logi TRUE
  .. .. ..$ eps            : num 0.01
  .. .. ..$ use_iter_solver: logi FALSE
  .. .. ..- attr(*, "class")= chr "ngme_control_f"
  .. ..$ n_params   : num 4
  .. ..$ debug      : logi FALSE
  .. ..$ par_string : chr " ignored    mu_1 sigma_1    nu_1"
  .. ..$ h          : num [1:24] 24 28.5 23.2 26.2 15.8 ...
  .. ..$ index      : num [1:597] 383 383 358 358 332 ...
  .. ..$ K          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  .. .. .. ..@ i       : int [1:48] 0 0 1 1 2 2 3 3 4 4 ...
  .. .. .. ..@ p       : int [1:26] 0 1 3 5 7 9 11 13 15 17 ...
  .. .. .. ..@ Dim     : int [1:2] 24 25
  .. .. .. ..@ Dimnames:List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : NULL
  .. .. .. ..@ x       : num [1:48] 1 -1 1 -1 1 -1 1 -1 1 -1 ...
  .. .. .. ..@ factors : list()
  .. ..$ G          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  .. .. .. ..@ i       : int [1:24] 0 1 2 3 4 5 6 7 8 9 ...
  .. .. .. ..@ p       : int [1:26] 0 1 2 3 4 5 6 7 8 9 ...
  .. .. .. ..@ Dim     : int [1:2] 24 25
  .. .. .. ..@ Dimnames:List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : NULL
  .. .. .. ..@ x       : num [1:24] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..@ factors : list()
  .. ..$ C          :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  .. .. .. ..@ i       : int [1:24] 0 1 2 3 4 5 6 7 8 9 ...
  .. .. .. ..@ p       : int [1:26] 0 0 1 2 3 4 5 6 7 8 ...
  .. .. .. ..@ Dim     : int [1:2] 24 25
  .. .. .. ..@ Dimnames:List of 2
  .. .. .. .. ..$ : NULL
  .. .. .. .. ..$ : NULL
  .. .. .. ..@ x       : num [1:24] -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
  .. .. .. ..@ factors : list()
  .. ..- attr(*, "class")= chr "ngme_model"