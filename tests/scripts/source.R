set.seed(16)
library(INLA)
load_all()

##############################  simulation
mesh2d <- inla.mesh.2d(
  loc.domain = cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5),
  max.edge = c(0.5, 2)
)
mesh2d$n
loc <- cbind(runif(300, 0, 10), runif(300, 0, 5))
# mesh2d$n
n_ar <- 10
arr <<- model_ar1(1:n_ar, alpha=0.7, debug=T)

plot(mesh2d)
points(loc, col="red", pch=19)

matern <<- model_matern(
  loc=loc,
  mesh=mesh2d, kappa = 2,
  debug=T)
K <- arr$K %x% matern$K
n <- nrow(K)

eps <- simulate(noise_nig(mu=-2, sigma=1, nu=1, n = n))
W <- drop(solve(K, eps))
A <- arr$A %x% matern$A
dim(A)

AW <- drop(A %*% W)
n_obs <- length(AW)
Y <- AW + rnorm(n_obs, sd=0.5)

# f(model=matern(mesh), group=ar1(1:3))

tp <- f(model=matern, group=arr); tp
expect_true(all(tp$A == A))
str(tp$model_right$noise)
##############################  estimation
load_all()
# str(f(model=matern, group=ar, noise=noise_nig())$model_right$noise)
matern$theta_K = -2
arr$theta_K = 1
out <- ngme(
  Y ~ 0 + f(
    model=matern, group=arr, noise=noise_nig(),
    debug=TRUE, control=ngme_control_f(
      numer_grad = T, use_precond = F, eps=0.01),
    fix_W = TRUE, W = W,
    fix_V = TRUE, V = attr(eps, "noise")$V
  ),
  data = list(Y=Y),
  family = "normal",
  control = ngme_control(
    iterations = 100,
    n_parallel_chain = 4,
    gibbs_sample = 5,
    estimation = T,
    max_relative_step = 1,
    max_absolute_step = 2
  ),
  debug = TRUE
)

out
traceplot(out, 1)

plot(out$latents[[1]]$noise,
 noise_nig(mu=-2, sigma=1, nu=1))



library(dplyr)

# Create N lists of numbers
list1 <- c(1, 2, 3, 4, 5)
list2 <- c(2, 4, 6, 8, 10)
list3 <- c(1, 3, 5, 7, 9)
lists <- list(list1, list2, list3)

# Define a function to split a single list into non-overlapping groups
split_list <- function(l, n) {
  # Split the list into n groups of equal length
  groups <- split(l, rep(1:n, length(l))[1:length(l)])
  # Return the list of groups
  return(groups)
}

# Apply the split_list() function to each list
split_lists <- lapply(lists, split_list, n = 2)

# Combine the non-overlapping groups from all lists into a single data frame
df <- data.frame()
for (i in 1:2) {
  for (j in 1:3) {
    df <- bind_rows(df, data.frame(group = i, list = j, values = split_lists[[j]][[i]]))
  }
}

# Split the combined data frame by the group variable
grouped_df <- df %>% group_split(group)

# Print the non-overlapping groups from all lists
for (i in 1:2) {
  cat(paste0("Group ", i, ":\n"))
  for (j in 1:3) {
    cat(paste0("  List ", j, ": ", paste0(grouped_df[[i]]$values[grouped_df[[i]]$list == j], collapse = " "), "\n"))
  }
}


# Define the matrix
mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                2, 4, 6, 8, 10, 1, 3, 5, 7, 9), ncol = 5)

# Define a function to split a matrix into non-overlapping submatrices
split_matrix <- function(m, n) {
  # Calculate the number of submatrices to create
  num_submatrices <- ceiling(ncol(m) / n)

  # Initialize a list to store the submatrices
  submatrices <- list()

  # Split the matrix into non-overlapping submatrices
  for (i in 1:num_submatrices) {
    submatrix_cols <- (i-1)*n+1:min(i*n, ncol(m))
    submatrix <- m[, submatrix_cols, drop = FALSE]
    submatrices[[i]] <- submatrix
  }

  # Return the list of submatrices
  return(submatrices)
}

# Apply the split_matrix() function to the matrix
submatrices <- split_matrix(mat, 2)

# Print the non-overlapping submatrices
for (i in 1:length(submatrices)) {
  cat(paste0("Submatrix ", i, ":\n"))
  print(submatrices[[i]])
}
