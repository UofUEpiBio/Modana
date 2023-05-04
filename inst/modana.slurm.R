library(dplyr)
library(tibble)
library(parallel)
library(ggplot2)
library(slurmR) # This also loads the parallel package
library(Modana)

#set.seed(111001)
set.seed(100)
#b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, 2, 1.5)
b0 <- c(-1, 1, 1, -1, 1, -1, 2, 0) # Multiple Interaction Terms 
b0 <- c(-1, 0.5, 0.5, -1, 1, 0, 0, 0) # Verify property one (Main Effect Estimation)
#b0 <- c(-1, 0.5, 0.5, -0.5, 1, 2) # Moderating Effect Estimation
N = c(100, 180, 500, 1000)
vartype <- c("continuous", "binary")
corr <- c(0.2, 0.5, 0.8)
nrun <- 1000

simul <- function(nrun){
  OUT <- NULL
  result <- NULL
  for(j in 1:length(N)){
    n <- N[j]
    for(k in 1:length(vartype)){
      predictors <- vartype[k]
      if(predictors=="continuous")
        uniform <- binary.Xs <- FALSE
      else uniform <- binary.Xs <- TRUE
      for(m in 1:length(corr)){
        rho <- corr[m]
        datt <- sim_data(n = n, b0 = b0, a0 = NULL, binary.Xs = binary.Xs,
                         sigma = 1, uniform = uniform, c0 = 1,
                         link.function = "logistic", rho = rho,
                         observational = FALSE, trt.p = 0.5,
                         interaction = 1:3, details = FALSE)
        getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3, # y ~ trt + x1 + x2 + x3 + x4,
                             detail = FALSE, y = "y",
                             trt = "trt", data = datt,
                             effmod = c("x1", "x2", "x3"),
                             corstr = "independence")
        
        # Verify property one (Main Effect Estimation)
        geecoefs <- summary(getres)[c(2,10:12),  1:2]
        inversecoefs <-summary(getres, model = 2)[c(2,6:8),  1:2]
        directcoefs <- summary(getres, model = 1)[c(2,6:8),  1:2]
        geeconfint <- confint(getres, parm = c(2,10:12), level = .95, model = 3)[, -1]
        inverseconfint<- confint(getres, parm = c(2,6:8), level = .95, model = 2)[, -1]
        directconfint <- confint(getres, parm = c(2,6:8), level = .95, model = 1)[, -1]

        # Moderating Effect Estimation
        
        #geecoefs <- summary(getres)[c(2, 10),  1:2]
        #inversecoefs <-summary(getres, model = 2)[c(2, 6),  1:2]
        #directcoefs <- summary(getres, model = 1)[c(2, 6),  1:2]
        #geeconfint <- confint(getres, parm = c(2, 10), level = .95, model = 3)[, -1]
        #inverseconfint<- confint(getres, parm = c(2, 6), level = .95, model = 2)[, -1]
        #directconfint <- confint(getres, parm = c(2, 6), level = .95, model = 1)[, -1]
        
        out1 <- rbind(geecoefs, directcoefs, inversecoefs)
        out2 <-  rbind(cbind(geeconfint, method = "gee"), cbind(directconfint, method = "direct"),
                       cbind(inverseconfint, method= "inverse"))
        out3 <- cbind(out1, out2, N = n, vartype = predictors, corcoef = rho,
                      Est = c(row.names(geecoefs), row.names(directcoefs), 
                              row.names(inversecoefs)))
        OUT <- rbind(OUT, out3)
      }
    }
  }
rbind(result, OUT)
}


ncores <- 2
njobs <- 20
niter <- 1000


job <- do.call(rbind.data.frame, Slurm_lapply(
  1:niter, simul,
  njobs    = njobs,
  mc.cores = ncores,
  plan     = "collect",
  tmp_path = "/scratch/general/nfs1/u1373115", # This is where all temp files will be exported
  sbatch_opt = list(
  account = "greene", #"notchpeak-shared-short",
   partition =  "lonepeak"), #"notchpeak-shared-short"), 
  export = c("b0", "N", "vartype", "corr")
  ))
res <- "~/Modana/jobarrays"
saveRDS(job, file = file.path(res, "maineff.rds"))