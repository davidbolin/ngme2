test_that("f() accepts mesh in model", {
  mesh <- fmesher::fm_mesh_1d(1:10)
  # f_args$mesh <- temp_mesh # mesh is not in f_args
  #
  # print(paste("DEBUG: model class:", paste(class(model), collapse = ", ")))
  out <- f(map = 1:10, model = ar1(mesh = mesh), noise = noise_normal())
  expect_true(inherits(out, "ngme_model"))
  expect_equal(out$operator$model, "ar1")
})

test_that("f() accepts list of meshes in model for replicates", {
  mesh1 <- fmesher::fm_mesh_1d(1:10)
  mesh2 <- fmesher::fm_mesh_1d(1:10)
  mesh_list <- list(mesh1, mesh2)

  # Replicate with different meshes
  out <- f(
    map = 1:20, replicate = rep(1:2, each = 10),
    model = ar1(mesh = mesh_list),
    noise = noise_normal()
  )

  expect_true(inherits(out, "ngme_model"))
  # Check if operator list is created internally (not directly exposed in ngme_model usually, but we can check structure)
  # Actually f() returns ngme_model which has operator.
  # If replicates, operator might be block diagonalized?
  # Yes, f() block diagonalizes.
  expect_equal(length(out$operator$h), mesh1$n + mesh2$n)
})

test_that("f() builds mesh from map if not provided", {
  out <- f(map = 1:10, model = ar1(), noise = noise_normal())
  expect_true(inherits(out, "ngme_model"))
  # 1:10 creates mesh of size 10 (min to max)
  expect_equal(out$mesh$n, 10)
})
