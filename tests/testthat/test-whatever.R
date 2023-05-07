
test_that("test names", {
  expect_equal(f(1:10, model="ar1", name="ar")$name, "ar")

  bm <- ngme(
    Y ~ f(1:10, name="f1", model="ar1") +
      f(1:10, name="f2", model="ar1"),
    data = data.frame(Y = rnorm(10)),
    control_opt = control_opt(estimation = FALSE)
  )

  expect_equal(names(bm[[1]]$models), c("f1", "f2"))
})

# tests on posterior sampling


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

test_that("merge_reps works", {
  repls <- c("1 1 2 2 2 3 4 5 5",
             "1 2 3 1 5 6 4 1 2",
             "1 1 1 1 1 2 2 3 3",
             "1 1 1 1 1 2 2 1 1")
  repls <- lapply(repls, function(x) as.numeric(strsplit(x, " ")[[1]]))
  repls

  # equal to 1 1 1 1 1 2 2 1 1
  expect_true(all(merge_repls(repls)
    == as.numeric(strsplit("1 1 1 1 1 2 2 1 1", " ")[[1]])))
})


test_that("new operator", {
  build_operator(
    model = "ar1",
    mesh = inla.mesh.1d(1:5),
    map = c(1:5, 1:3),
    replicate = c(rep(1, 5), rep(2, 3)),
    alpha = -0.7
  )

  load_all()
  a = f(model="ar1", map=c(1:5, 1:3), replicate=c(rep(1, 5), rep(2,3)), alpha=-0.7, eval=T)
str(a)

bv(first=ar1(1:3), second=ar1(1:3))

ar1(1:3)
  kronecker(diag(2), C)
  inla.spde.make.A(loc=c(4,5,6), mesh = inla.mesh.1d(4:6))

  # R interface: a few examples
  f(model="bv", first = ar1(1:3), second = ar1(1:3),
    noise = c("normal", "nig"),
    eval=T)


})