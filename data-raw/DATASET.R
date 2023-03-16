## code to prepare `DATASET` dataset goes here
a0 <- c(-2, 2, 2, 0, 1)
b0 <- c(-1, 0.5, 0.5, -1, .5, -0.5, 1, 0, .5)
source("~/Modana/Modana/R/moderation.R")
DATASET <- read.csv("C:/Users/Admin/OneDrive - University of Utah/Documents/wart.csv")
SimData <- sim_data(n = 100, b0, a0 = NULL, binary.Xs = FALSE,
                    sigma = 1, uniform = FALSE, c0 = 1,
                    link.function = "logistic", rho = 0.2,
                    observational = FALSE, trt.p = 0.5,
                    interaction = 1:2, details = FALSE)

usethis::use_data(DATASET, SimData, overwrite = TRUE)
