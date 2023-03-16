test_that("Testing the number of objects in my output", {
  a0 <- c(-2, 2, 2, 0, 1)
  b0 <- c(-1, 0.5, 0.5, -1, .5, -0.5, 1, 0, .5)
  source("~/Modana/Modana/R/moderation.R")
  datt <- sim_data(n = 100, b0, a0 = NULL, binary.Xs = FALSE,
                   sigma = 1, uniform = FALSE, c0 = 1,
                   link.function = "logistic", rho = 0.2,
                   observational = FALSE, trt.p = 0.5,
                   interaction = 1:2, details = FALSE)
  getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3,
                       detail = TRUE, y = "y", trt = "trt",
                       data = datt, effmod = c("x1", "x2"),
                       corstr = "independence")
  expect_equal(length(getres), 5)
})

