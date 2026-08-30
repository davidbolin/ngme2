# One mesh per replicate: f(map, model = <op>(mesh = <list of meshes>)) together
# with ngme(replicate = ...) must build the operator and the A matrix of each
# replicate against that replicate's own mesh.

make_1d_repl_data <- function(n_per, seed = 23) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(n_per), function(i) {
    data.frame(
      y = rnorm(n_per[i]),
      t = sort(runif(n_per[i], 0, 10)),
      rep = i
    )
  }))
}

mesh_sizes <- function(mesh_list) {
  vapply(mesh_list, function(m) as.numeric(m$n), numeric(1), USE.NAMES = FALSE)
}

test_that("ngme_make_mesh_repls builds one mesh per replicate level", {
  n_per <- c(20, 30, 40)
  dm <- make_1d_repl_data(n_per)

  mr <- ngme_make_mesh_repls(dm, dm$t, dm$rep)

  expect_length(mr, length(n_per))
  expect_equal(names(mr), c("1", "2", "3"))

  expected <- vapply(
    seq_along(n_per),
    function(i) as.numeric(fmesher::fm_mesh_1d(sort(dm$t[dm$rep == i]))$n),
    numeric(1)
  )
  expect_equal(mesh_sizes(mr), expected)

  # the mesh of a replicate must not depend on the row order of the data
  shuffled <- dm[sample(nrow(dm)), ]
  mr_shuffled <- ngme_make_mesh_repls(shuffled, shuffled$t, shuffled$rep)
  expect_equal(names(mr_shuffled), names(mr))
  expect_equal(mesh_sizes(mr_shuffled), expected)
})

test_that("per-replicate 1D mesh list gives each replicate its own operator", {
  n_per <- c(60, 90, 120)
  dm <- make_1d_repl_data(n_per)
  mr <- ngme_make_mesh_repls(dm, dm$t, dm$rep)

  n_nodes <- mesh_sizes(mr)
  expect_length(mr, 3)
  expect_false(any(duplicated(n_nodes))) # meshes of unequal size

  fit <- ngme(
    y ~ f(t, model = rw1(mesh = mr), noise = noise_normal()),
    data = dm,
    replicate = dm$rep,
    family = "normal",
    control_opt = control_opt(
      seed = 1, burnin = 5, iterations = 10, n_batch = 1
    )
  )

  expect_length(fit$replicates, 3)
  for (i in seq_along(n_per)) {
    model <- fit$replicates[[i]]$models[[1]]
    expect_equal(as.numeric(model$mesh$n), n_nodes[i])
    expect_equal(as.numeric(model$W_size), n_nodes[i])
    expect_equal(ncol(model$A), n_nodes[i])
    expect_equal(nrow(model$A), n_per[i])
    expect_equal(ncol(model$operator$K), n_nodes[i])
  }
})

test_that("per-replicate 2D mesh list gives each replicate its own operator", {
  set.seed(7)
  boundary <- cbind(c(0, 5, 5, 0, 0), c(0, 0, 5, 5, 0))
  meshes <- list(
    fmesher::fm_mesh_2d(loc.domain = boundary, max.edge = c(1.5, 3), cutoff = 0.6),
    fmesher::fm_mesh_2d(loc.domain = boundary, max.edge = c(0.8, 3), cutoff = 0.3)
  )
  n_nodes <- mesh_sizes(meshes)
  expect_false(any(duplicated(n_nodes))) # meshes of unequal size

  n_per <- c(80, 120)
  dm <- do.call(rbind, lapply(seq_along(n_per), function(i) {
    data.frame(
      y = rnorm(n_per[i]),
      x1 = runif(n_per[i], 0.5, 4.5),
      x2 = runif(n_per[i], 0.5, 4.5),
      rep = i
    )
  }))

  fit <- ngme(
    y ~ f(cbind(x1, x2), model = matern(mesh = meshes), noise = noise_normal()),
    data = dm,
    replicate = dm$rep,
    family = "normal",
    control_opt = control_opt(
      seed = 1, burnin = 5, iterations = 10, n_batch = 1
    )
  )

  expect_length(fit$replicates, 2)
  for (i in seq_along(n_per)) {
    model <- fit$replicates[[i]]$models[[1]]
    expect_equal(as.numeric(model$mesh$n), n_nodes[i])
    expect_equal(as.numeric(model$W_size), n_nodes[i])
    expect_equal(ncol(model$A), n_nodes[i])
    expect_equal(nrow(model$A), n_per[i])
  }
})

test_that("per-replicate mesh list is picked up from a pre-built operator def", {
  n_per <- c(20, 30)
  dm <- make_1d_repl_data(n_per)
  mr <- ngme_make_mesh_repls(dm, dm$t, dm$rep)
  op <- rw1(mesh = mr)
  expect_s3_class(op, "ngme_operator_def")

  fit <- ngme(
    y ~ f(t, model = op, noise = noise_normal()),
    data = dm,
    replicate = dm$rep,
    family = "normal",
    control_opt = control_opt(estimation = FALSE)
  )

  for (i in seq_along(n_per)) {
    model <- fit$replicates[[i]]$models[[1]]
    expect_equal(as.numeric(model$W_size), mesh_sizes(mr)[i])
    expect_equal(ncol(model$A), mesh_sizes(mr)[i])
  }
})

test_that("too few meshes for the number of replicates is reported", {
  n_per <- c(20, 30, 40)
  dm <- make_1d_repl_data(n_per)
  mr <- ngme_make_mesh_repls(dm, dm$t, dm$rep)[1:2]

  expect_error(
    ngme(
      y ~ f(t, model = rw1(mesh = mr), noise = noise_normal()),
      data = dm,
      replicate = dm$rep,
      family = "normal",
      control_opt = control_opt(estimation = FALSE)
    ),
    "Insufficient meshes"
  )
})

test_that("f() rejects a mesh list without replicate information", {
  meshes <- list(fmesher::fm_mesh_1d(1:10), fmesher::fm_mesh_1d(1:15))
  expect_error(
    f(map = 1:20, model = rw1(mesh = meshes), noise = noise_normal()),
    "replicate"
  )
})
