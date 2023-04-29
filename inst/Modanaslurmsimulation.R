library(Modana)

#set.seed(111001)
set.seed(1001)
b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, 2, 1.5)
a0 <- c(-1, 1.5, -1, 1.5)
N = c(100, 180, 500, 1000)
vartype <- c("continuous", "binary")
corr <- c(0.2, 0.5, 0.8)
nrun <- 10
grid <- expand.grid(N = N, vartype = vartype, corr = corr)




simf <- function(kk){
    results <- NULL
    OUT <- NULL
  for(j in 1:dim(grid)[1]){
    n <- grid[j, "N"]
    predictors <- grid[j, "vartype"]
    rho <- grid[j, "corr"]
      if(predictors=="continuous")
        uniform <- binary.Xs <- FALSE
      else uniform <- binary.Xs <- TRUE
        resu <- do.call(rbind.data.frame, lapply(1:kk, \(dat){
          datt <- sim_data(n = n, b0 = b0,
                          a0 = NULL, binary.Xs = binary.Xs,
                         sigma = 1, uniform = uniform, c0 = 1,
                         link.function = "logistic", rho = rho,
                         observational = FALSE, trt.p = 0.5,
                         interaction = 1:2, details = FALSE)
        getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3 + x4,
                             detail = FALSE, y = "y",
                             trt = "trt", data = datt,
                             effmod = c("x1", "x2"),
                             corstr = "independence")
        geecoefs <- summary(getres)[c(2,12,13),  1:2]
        inversecoefs <-summary(getres, model = 2)[c(2,7, 8),  1:2]
        directcoefs <- summary(getres, model = 1)[c(2,7,8),  1:2]
        geeconfint <- confint(getres, parm = c(2,12,13), level = .95, model = 3)[, -1]
        inverseconfint<- confint(getres, parm = c(2,7,8), level = .95, model = 2)[, -1]
        directconfint <- confint(getres, parm = c(2,7,8), level = .95, model = 1)[, -1]
        out1 <- rbind(geecoefs, directcoefs, inversecoefs)
        out2 <-  rbind(cbind(geeconfint, method = "gee"), cbind(directconfint, method = "direct"),
                       cbind(inverseconfint, method= "inverse"))
        out3 <- cbind(out1, out2, N = n, vartype = predictors, corcoef = rho)
        #OUT <- rbind(OUT, out3)
        return(out3)
      }))
        OUT <- rbind(OUT, resu) 
  }
    return(OUT)
  }
   

getR <- simf(10)












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
                         interaction = 1:2, details = FALSE)
        getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3 + x4,
                             detail = FALSE, y = "y",
                             trt = "trt", data = datt,
                             effmod = c("x1", "x2"),
                             corstr = "independence")
        geecoefs <- summary(getres)[c(2,12,13),  1:2]
        inversecoefs <-summary(getres, model = 2)[c(2,7, 8),  1:2]
        directcoefs <- summary(getres, model = 1)[c(2,7,8),  1:2]
        geeconfint <- confint(getres, parm = c(2,12,13), level = .95, model = 3)[, -1]
        inverseconfint<- confint(getres, parm = c(2,7,8), level = .95, model = 2)[, -1]
        directconfint <- confint(getres, parm = c(2,7,8), level = .95, model = 1)[, -1]
        out1 <- rbind(geecoefs, directcoefs, inversecoefs)
        out2 <-  rbind(cbind(geeconfint, method = "gee"), cbind(directconfint, method = "direct"),
                       cbind(inverseconfint, method= "inverse"))
        out3 <- cbind(out1, out2, N = n, vartype = predictors, corcoef = rho)
        OUT <- rbind(OUT, out3)
      }
    }
  }
rbind(result, OUT)
}




cl <- parallel::makeCluster(10)
parallel::clusterExport(cl, list("b0", "N", "vartype", "corr"))
parallel::clusterEvalQ(cl, library(Modana))
fullres <- do.call(rbind.data.frame, parallel::parLapply(cl, 1:10, simul))
parallel::stopCluster(cl)



#set.seed(111001)
set.seed(1001)
b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, 2, 1.5)
a0 <- c(-1, 1.5, -1, 1.5)
N = c(100, 180, 500, 1000)
vartype <- c("continuous", "binary")
corr <- c(0.2, 0.5, 0.8)
nrun <- 10
OUT <- NULL
results <- NULL
for(i in 1:nrun){
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
                         interaction = 1:2, details = FALSE)
        getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3 + x4,
                             detail = FALSE, y = "y",
                             trt = "trt", data = datt,
                             effmod = c("x1", "x2"),
                             corstr = "independence")
        geecoefs <- summary(getres)[c(2,12,13),  1:2]
        inversecoefs <-summary(getres, model = 2)[c(2,7, 8),  1:2]
        directcoefs <- summary(getres, model = 1)[c(2,7,8),  1:2]
        geeconfint <- confint(getres, parm = c(2,12,13), level = .95, model = 3)[, -1]
        inverseconfint<- confint(getres, parm = c(2,7,8), level = .95, model = 2)[, -1]
        directconfint <- confint(getres, parm = c(2,7,8), level = .95, model = 1)[, -1]
        out1 <- rbind(geecoefs, directcoefs, inversecoefs)
        out2 <-  rbind(cbind(geeconfint, method = "gee"), cbind(directconfint, method = "direct"),
                       cbind(inverseconfint, method= "inverse"))
        out3 <- cbind(out1, out2, N = n, vartype = predictors, corcoef = rho)
        OUT <- rbind(OUT, out3)
      }
    }
  }
  #results <- rbind(results, OUT) 
}









