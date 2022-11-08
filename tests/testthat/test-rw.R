test_that("the order of W same as order of index?", {
  devtools::load_all()
  XX    <- c(1.1, 3.1, 2.2, 4.5)
  myloc <- c(3.2, 1.2)

  mesh1 <- inla.mesh.1d(XX); mesh1
  rw1 <- f(XX, model = "rw1", index_NA = c(F, F, F, T))
  rw1$A
  rw1$A_pred

  # use INLA make A? the order is bad
  A <- inla.spde.make.A(mesh=mesh1, loc = myloc); A
# 2 x 4 sparse Matrix of class "dgCMatrix"
# [1,] .         .         0.7142857 0.2857143
# [2,] 0.4545455 0.5454545 .         .

  # the order is wrong, recover by permutation matrix
  A %*% as(mesh1$idx$loc, "pMatrix")

  rw1$C + rw1$G
# 3 x 4 sparse Matrix of class "dgCMatrix"
# [1,] 1 -1  .  .
# [2,] .  1 -1  .
# [3,] .  .  1 -1

  A %*% mesh1$loc

  sorted_index <- sort(XX, index.return = TRUE)
  A2 <- Matrix::sparseMatrix(seq_along(XX), sorted_index$ix, x = 1); A2

  # rw2
  rw2 <- model_rw(1:4, order = 2, circular = F)
  rw2$C + rw2$G

  sparseMatrix(1:10, 1:10, x = 1)
})
