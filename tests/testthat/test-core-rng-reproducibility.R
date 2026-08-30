# The GIG-based samplers draw from their own C++ stream rather than R's. They
# take a `seed`, and rGIG_cpp treats seed == 0 as "seed from the system clock".
# Their default must therefore come from R's stream (ngme_random_seed()), or
# set.seed() would silently fail to make them reproducible.

test_that("rnig(), rig() and rgig() are reproducible under set.seed()", {
  draw <- function(f) {
    set.seed(20240101)
    f()
  }
  samplers <- list(
    rnig = function() rnig(10, delta = -2, mu = 2, sigma = 4, nu = 0.3),
    rig  = function() rig(10, a = 1, b = 1),
    rgig = function() rgig(10, p = -0.5, a = 1, b = 1)
  )
  for (nm in names(samplers)) {
    expect_equal(draw(samplers[[nm]]), draw(samplers[[nm]]),
                 info = paste(nm, "must respect set.seed()"))
  }
})

test_that("different seeds give different draws", {
  set.seed(1); a <- rnig(10, delta = -2, mu = 2, sigma = 4, nu = 0.3)
  set.seed(2); b <- rnig(10, delta = -2, mu = 2, sigma = 4, nu = 0.3)
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("an explicit seed still drives the GIG stream", {
  # rnig() mixes two streams: the mixing variable V comes from the C++ GIG
  # sampler (driven by `seed`), the Gaussian part from R's own generator. Both
  # have to be pinned for the draw to repeat.
  draw <- function() {
    set.seed(11)
    rnig(10, delta = -2, mu = 2, sigma = 4, nu = 0.3, seed = 7)
  }
  expect_equal(draw(), draw())

  set.seed(11); s7 <- rnig(10, delta = -2, mu = 2, sigma = 4, nu = 0.3, seed = 7)
  set.seed(11); s8 <- rnig(10, delta = -2, mu = 2, sigma = 4, nu = 0.3, seed = 8)
  expect_false(isTRUE(all.equal(s7, s8)))
})

test_that("simulate() of a latent model is reproducible under set.seed()", {
  m <- f(1:50, model = ar1(rho = 0.4), noise = noise_nig())
  set.seed(42); s1 <- simulate(m, nsim = 1)
  set.seed(42); s2 <- simulate(m, nsim = 1)
  expect_equal(s1, s2)
})
