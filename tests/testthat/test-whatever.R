# test_that("modify_ngme_with_idx_NA", {
#   # load_all()
#   X2 = c(3,4,5,9,2)
#   spde1 <<- model_matern(mesh = 1:10, map=X2)
#   mm <- ngme(
#     YY ~ 1 + X1 + f(X2, model="rw1") + f(model="matern", mesh=spde1$mesh, map=X2),
#     data = data.frame(
#       YY=c(0,2,3,4,5),
#       X1=c(3,4,5,6,7),
#       X2=X2
#     ),
#     control_opt=control_opt(iterations = 10, print_check_info = FALSE)
#   )
#   mm
#   new_model <- modify_ngme_with_idx_NA(mm, idx_NA = 3)
#   str(new_model$noise)
#   # new_model$models[[2]]$A_pred

#   # new_model$models[[1]]$A_pred
#   expect_true(length(new_model$Y) == 4)
#   # new_model$models[[1]]
# })

# test_that("merge_reps works", {
#   repls <- c("1 1 2 2 2 3 4 5 5",
#              "1 2 3 1 5 6 4 1 2",
#              "1 1 1 1 1 2 2 3 3",
#              "1 1 1 1 1 2 2 1 1")
#   repls <- lapply(repls, function(x) as.numeric(strsplit(x, " ")[[1]]))
#   repls

#   # equal to 1 1 1 1 1 2 2 1 1
#   expect_true(all(merge_repls(repls)
#     == as.numeric(strsplit("1 1 1 1 1 2 2 1 1", " ")[[1]])))
# })


test_that("compute 3d precision matrix", {
  cor_mat <- matrix(0, nrow=3, ncol=3)
  cor_mat[1,3] = 0.1
  cor_mat[1,2] = 0.2
  cor_mat[2,3] = 0.3
})

test_that("add priors", {
  # ar1
  out <- ngme(
    YY ~ 1 + f(X1,
      model="ar1",
      noise=noise_nig()
    ),
    data = data.frame(
      YY=c(0,2,3,4,5),
      X1=c(3,4,5,6,7)
    ),
    control_opt=control_opt(
      iterations = 10,
      print_check_info = FALSE,
      estimation = T
    )
  )

  out$replicates[[1]]$noise$prior_mu

  # parameters, mu, sigma, nu
  # priors on mu, sigma, nu
})



