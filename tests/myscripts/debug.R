# library(devtools)
# load_all()
# usethis::use_pkgdown()

# pkgdown::build_site(devel = TRUE)

# # build website


# library(INLA)
# mesh1 <- inla.mesh.1d(rnorm(100))
# mesh1

# inla.spde.make.A(mesh=mesh1, loc = 4)

# t.data.frame
# predict.glm

# ?scale

# ?print.lm
# ?lm

# y <- 1:100
# x <- rnorm(100)
# m1 <- (lm(y~x))
# ??print.lm
# print.summary.lm

# ??print.ngme
# ?print.lm
# ??ngme

# ?use_import_from()
# library(devtools)
# use_import_from("rlang", ".data")


# library(grid)
# grid.arrange(rectGrob(), rectGrob())
# ## Not run:
# library(ggplot2)
# pl <- lapply(1:3, function(.x) qplot(1:10, rnorm(10), main=paste("plot", .x)))
# ml <- marrangeGrob(pl, nrow=2, ncol=2)
# length(pl)
# ## non-interactive use, multipage pdf
# ggsave("multipage.pdf", ml)
# ## interactive use; open new devices
# ml

# str(out$noise)
# ?switch
# pl
# c(list(a=1, b=2), list(a=2))

# str(out$latents[[1]][["theta_mu"]])

# library(devtools)
# use_vignette("model-estimation")

# devtools::load_all()
# m = model_rw(c(1.1, 2.2, 3.5, 5.6), order=1); m$C + m$G

# m = model_ar1(c(1, 2, 3, 10), order=1, noise=noise_nig()); m$C + m$G

# str(m)
# check()


# # sim rw1 model
# devtools::load_all()
# x <- sort(rexp(10))
# model_rw(x)$noise$n_noise

# solve(matrix(rexp(6), nrow=2), b=c(1,2))

# ff <- function(a=1, ...) {
#   print(list(...))
# }
# ff(b=3)

# # double dlambda_V(const double loglambda,
# #                  const Eigen::VectorXd &V,
# #                  const Eigen::VectorXd &h,
# #                  const int GAL)
# # {
# #   double dlambda = 0;
# #   if(GAL){
# #   for(int i=0; i < h.size(); i++){
# #     double h_lambda = exp(loglambda) * h[i];
# #     //digamma(0.1) = digamma(1.1) - 1/0.1;
# #     if(h_lambda > 1){
# #         dlambda -=  h_lambda * R::digamma(h_lambda);
# #       }else
# #       {
# #         dlambda -=  h_lambda * R::digamma(h_lambda + 1) - 1.;
# #       }
# #     dlambda += h_lambda *  log(V(i) ) ;
# #   }
# #   }else{
# #     double srqt_two = pow(2, 0.5);
# #     for(int i=0; i < h.size(); i++){
# #       dlambda +=  1 -  ( pow(h(i), 2) / V(i) ) * exp( 2 * loglambda);
# #       dlambda += srqt_two * h(i)  * exp(loglambda);
# #     }

# #   }

# #   return(dlambda);
# # }


# # to-do
# # replicates
# m3 <- model_rw(c(1.1, 2.2, 2.2, 3.3, 4.4), replicates = c(1, 1, 2, 1, 1))
# m3$K
# m3$A
# #  [1,] 1 . . . . . . .
# #  [2,] . 1 . . . . . .
# #  [3,] . . . . . 1 . .
# #  [4,] . . 1 . . . . .
# #  [5,] . . . 1 . . . .

######## test gal case
# rgal(100, delta = 0, mu = 5, sigma = 1, nu = 1)
library(devtools)
load_all()
n_obs <- 500
gal <- noise_gal(mu=0, sigma=1, nu=1, n=n_obs)
es <- simulate(gal)
# plot(es, type="l")

######## test ar with gal noise
ww <- simulate(model_ar1(1:n_obs, alpha=0.7, noise=noise_gal(mu=1, sigma=2, nu=1.5)))
yy <- ww + rnorm(n_obs)

# make sure length(nosie$V) == W
res <- ngme(
  y ~ 0 + f(1:n_obs,
    model = "ar1",
    noise = noise_gal(
      V = attr(es, "noise")$V,
      # fix_V = TRUE
    ),
    debug=TRUE),
  family = "normal",
  data = list(y = yy),
  control = ngme_control(
    gibbs_sample = 5,
    iterations = 500,
    n_parallel_chain = 1
  ),
  debug = TRUE
)
res
?rgamma




xx <- seq(-5,5,0.01)
plot(xx, 2*dnorm(xx, 0, 2), type="l")
lines(density(rnorm(100000, 0, 1)), col="red")

# rejection sampling
niter <- 100000
accept_x <- 0
for (iter in 1:niter) {
  x <- rnorm(1, 0, 2)
  alpha <- dnorm(x, 0, 1) /  (2*dnorm(x, 0, 2))
  if (runif(1) < alpha) {
    accept_x <- c(accept_x, x)
  }
}
length(accept_x) / niter
plot(density(accept_x))

# test distribution function
# Gamma(x; shape=h*nu, rate=nu)
# h, x known, maximize nu
h <- 1

grad_nu <- function(nu, x=data, h=1) {
  g <- h - x - h*log(1/nu) + h*log(x) - h*digamma(h*nu)
  mean(g)
}

data <- rgamma(100, shape=1.98, rate=1.98)
log_post <- function(nu) {
  sum(dgamma(data, shape=nu, rate=nu, log=TRUE))
}

optim(par=1, fn=function(nu) {
  -log_post(nu)
}, method="Brent", lower=0.01, upper=100)

grad_nu(x=data, h=1, nu=2.62)



#### test ig distribution

grad_nu_ig <- function(x, nu, h=1) {
  mean(0.5*(2*h + 1/nu - h^2/x - x))
}

h = 1; nu=1.2
igs <- rig(10000, nu, nu*h^2)

log_post_ig<- function(nu) {
  sum(dig(igs, a=nu, b=nu*h^2, log=TRUE))
}

optim(par=1, fn=function(nu) {
  -log_post_ig(nu)
}, method="Brent", lower=0.01, upper=100)

grad_nu_ig(igs, 1.2)

library(numDeriv)

grad(log_post, 1.98)

grad_nu(1)


# testing noise_normal_nig in R
load_all()
eps <- noise_normal_nig(
  sigma_normal = 2,
  mu = 1,
  nu = 0.5,
  sigma_nig = 1.3
)
eps
str(eps)

str(f(1:10, model="ar1", noise=noise_gal()))
str(f(1:10, model="ar1", noise=eps))

# testing noise_normal_nig in C
