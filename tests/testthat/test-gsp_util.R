test_that("Polyshift works", {
  expect_equal(
    {
      coef <- 1:4
      x <- 4
      coef*(x^(0:3))
      sum((PolyShift(3,4) %*% coef)*(x-3)^(0:3) ) -
      sum(coef*(x^(0:3)))
    }, 0)
})


