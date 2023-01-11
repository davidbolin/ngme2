load_all()
mu <- 2; nu <- 1; sigma <- 1

xx <- seq(-5, 5, length.out = 1000)
dd <- dnig(xx, delta = -mu, mu, nu, sigma)
plot(xx, dd, type = "l")


?dnig
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)
nig_df_fun <-  function(x, sigma, mu){
  data.frame(x = x, dnig = ngme2::dnig(x=x,
                                      delta=0,
                                      mu=mu, nu=1, sigma))
}
params <- expand.grid(sigma = c(1, 2, 4, 8), mu = c(-5,0,5))
nig_par <- mdply(params, nig_df_fun, x = seq(-20, 20, length=400))
nig_par <- nig_par %>% mutate(label = paste0("mu = ",mu))

ggplot(nig_par, mapping = aes(x = x, y=dnig, colour=factor(sigma)))+
  facet_wrap(~label,scales="free", nrow=3) +
  geom_line() +
  ylab("NIG density") + labs(colour = "sigma")


plot.noise(ngme.noise.nig(
    theta_mu = 0.5,
    theta_sigma = 1,
    nu = 1
))

plot(dnig(xx, delta = -1, 1, nu, sigma))
lines(dnig(xx, delta = -4, 4, nu, sigma))
