test_that("Matern2d formula works", {
  fm <- y1 | y2 ~ f(model = "matern2d",
    loc = list(x1, x2), mesh = list(mesh1, mesh2))
  Fm <- Formula::Formula(fm)
  length(Formula(y1|y2 ~ x1))
  Formula(y1|y2 ~ x1)

  library(Formula)
  load_all()
  data <- data.frame(y1 = 1:10, y2 = 2:11, x1 = 1:10, x2 = 1:10)
  res <- ngme_parse_formula(fm, data, control_ngme(), noise_nig())

  y1 <- formula(Fm, lhs=1, rhs=0)[[2]]
  eval(y1, data)
## create a simple Formula with one response and two regressor parts
})
