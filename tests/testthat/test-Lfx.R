test_that("Lfx works", {
  expect_equal(
    {
      fit <- lm(mpg ~ cyl + disp, mtcars)
      L <- Lfx(fit, list(0,0*cyl,1))
      round(c(wald(fit, L)[[1]]$anova$`F-value`),4)
    }, 4.0268)
})
