f <- function(a) g(a)
g <- function(b) h(b)
h <- function(c) i(c)
i <- function(d) {
  if (!is.numeric(d)) {
    stop("`d` must be numeric", call. = FALSE)
  }
  d + 10
}
f(b=1)

g <- function(b) {
  browser()
  h(b)
}
f(10)

use_test("f")
