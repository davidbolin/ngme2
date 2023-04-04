test_that("Matern2d formula works", {


  Y ~ f(model="matern2d", map=list(x1, x2), names=list("matern1", "matern2"), rho=0.5, theta=1, alpha=2)
})
