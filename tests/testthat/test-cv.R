test_that("test merged replicate", {
  # creating ar1 process of 2 replicates
  n <- 5
  df <- data.frame(
    Y = 1:10,
    idx = 1:10,
    group = rep(1:2, each = n)
  )

  ar_idx <- c(1,2,3,4,5, 3,2,1,5,4)
  rw_idx <- c(5,4,3,2,1, 1,4,3,5,2)
  mod.replicates <- ngme(
    formula = Y ~ 1 + f(
      ar_idx,
      model = "ar1",
      noise = noise_nig()
    ) + f(
      rw_idx,
      model = "rw1",
      noise = noise_nig()
    ),
    replicate = df$group,
    data = df,
    control_opt = control_opt(
      print_check_info = FALSE
    )
  )

  ret <- merge_replicates(
    mod.replicates, 
    train_idx = c(1:4, 6:9),
    test_idx = c(5, 10)
  )

  expect_equal(ret$test_Y, c(5,10))

  m <- ret$merged_rep
  expect_equal(m$Y, c(1,2,3,4,6,7,8,9))
  print(m$models[[1]]$A)
  print(m$models[[2]]$A)

  print(ret$test_A_block)
  # [1,] . . . 0 1 1 0 . . .
  # [2,] . . . 1 0 . 1 0 . .
  # Correct


  ########## Compute scores ##########
  s = compute_scores(
    m, 10, 10, 1,
    A_pred_block = ret$test_A_block,
    test_noise = ret$test_noise,
    y_data = ret$test_Y,
    group_data = ret$test_group,
    X_pred = ret$test_X,
    transform = identity,
    thining_gap = 0
  )
  s
})


test_that("test cross validation (NIG and gaussian)", {
  # load fitted models
  ret_gauss <- readRDS("examples/models/ret_gauss.rds")
  ret_nig <- readRDS("examples/models/ret_nig.rds")

  seed <- 500
  cv <- cross_validation(
    list(gauss=ret_gauss, nig=ret_nig),
    k=20,
    parallel=FALSE,
    n_gibbs_samples=100,
    print=TRUE,
    N_sim=4,
    seed=seed+50
  )
  cv
  gauss_mse <- cv$mean.scores$MSE[1]; gauss_mae <- cv$mean.scores$MAE[1]
  nig_mse <- cv$mean.scores$MSE[2]; nig_mae <- cv$mean.scores$MAE[2]
  expect_true(gauss_mse > nig_mse)
  expect_true(gauss_mae > nig_mae)
})