
# convergence check code

n_param <- 4
N <- N_chain <- 2

m = matrix(rnorm(8), nrow=N_chain, ncol=n_param)
sigma2 <- matrix(abs(rnorm(8)), nrow=N_chain, ncol=n_param)

n.points <- 3
sigma2[(N-n.points+1):N,i]

for(i in 1:n.test){
  std.satisfied <- sqrt(sigma2[N,i])/abs(m[N,1]) < std.lim
  Sigma <- diag(sigma2[(N-n.points+1):N,i])
  Q <- solve(t(B)%*%solve(Sigma,B))
  beta <- Q%*%(t(B)%*%solve(Sigma,m[(N-n.points+1):N,i]))
  slope.satisfied <- abs(beta[2])-2*sqrt(Q[2,2])<trend.lim*abs(beta[1]) #no significant trend
  output[i] = std.satisfied&slope.satisfied
}

