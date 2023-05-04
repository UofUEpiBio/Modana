library(Modana)

set.seed(111001)
#set.seed(100)
#b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, 2, 1.5)
b0 <- c(-1, 1, 1, -1, 1, -1, 2, 0) # Multiple Interaction Terms 
b0 <- c(-1, 0.5, 0.5, -1, 1, 0, 0, 0) # Verify property one (Main Effect Estimation)
b0 <- c(-1, 0.5, 1, -1.5, 1, 1.5) # Moderating Effect Estimation
#a0 <- c(-1, 1.5, -1, 1.5)
N = c(100, 180, 500, 1000)
vartype <- c("continuous", "binary")
corr <- c(0.2, 0.5, 0.8)
nrun <- 500
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
                         interaction = 1, details = FALSE)
        getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3, # y ~ trt + x1 + x2 + x3 + x4,
                             detail = FALSE, y = "y",
                             trt = "trt", data = datt,
                             effmod = c("x1"), #c("x1")
                             corstr = "independence")
        # geecoefs <- summary(getres)[c(2,12,13),  1:2]
        # inversecoefs <-summary(getres, model = 2)[c(2,7, 8),  1:2]
        # directcoefs <- summary(getres, model = 1)[c(2,7,8),  1:2]
        # geeconfint <- confint(getres, parm = c(2,12,13), level = .95, model = 3)[, -1]
        # inverseconfint<- confint(getres, parm = c(2,7,8), level = .95, model = 2)[, -1]
        # directconfint <- confint(getres, parm = c(2,7,8), level = .95, model = 1)[, -1]
        
      # Verify property one (Main Effect Estimation)
        # geecoefs <- summary(getres)[c(2,10:12),  1:2]
        # inversecoefs <-summary(getres, model = 2)[c(2,6:8),  1:2]
        # directcoefs <- summary(getres, model = 1)[c(2,6:8),  1:2]
        # geeconfint <- confint(getres, parm = c(2,10:12), level = .95, model = 3)[, -1]
        # inverseconfint<- confint(getres, parm = c(2,6:8), level = .95, model = 2)[, -1]
        # directconfint <- confint(getres, parm = c(2,6:8), level = .95, model = 1)[, -1]

        # Moderating Effect Estimation
        
        geecoefs <- summary(getres)[c(2, 10),  1:2]
        inversecoefs <-summary(getres, model = 2)[c(2, 6),  1:2]
        directcoefs <- summary(getres, model = 1)[c(2, 6),  1:2]
        geeconfint <- confint(getres, parm = c(2, 10), level = .95, model = 3)[, -1]
        inverseconfint<- confint(getres, parm = c(2, 6), level = .95, model = 2)[, -1]
        directconfint <- confint(getres, parm = c(2, 6), level = .95, model = 1)[, -1]
        
        out1 <- rbind(geecoefs, directcoefs, inversecoefs)
        out2 <-  rbind(cbind(geeconfint, method = "gee"), cbind(directconfint, method = "direct"),
                       cbind(inverseconfint, method= "inverse"))
        out3 <- cbind(out1, out2, N = n, vartype = predictors, corcoef = rho,
                      Est = c(row.names(geecoefs), row.names(directcoefs), 
                              row.names(inversecoefs)))
        #OUT <- rbind(OUT, out3)
        return(out3)
      }))
        OUT <- rbind(OUT, resu) 
  }
    return(OUT)
}




   

getR <- simf(kk=500)
 
singleeffmod <- getR 
fn <- system.file(
  file.path("maineff.rds"),
  package = "Modana"
)

maineff <- file.path(.libPaths(), "Modana", "inst", "maineff.rds")
maineff <- readRDS("maineff.rds")
effmod <- readRDS("effmod.rds")
library(dplyr)
library(ggplot2)
func <- \(x) as.numeric(as.character(x))
plotdata <- maineff %>% dplyr::mutate_at(c("lwr", "upr"), func) %>% 
  tibble::remove_rownames() %>% as.data.frame()

#Summary table for Main effect
plotdata %>% 
  filter(Est %in% c("trt", "y"), corcoef == .5,
         N %in% c(100, 1000),
         vartype %in% "continuous") %>% 
  group_by(method, N) %>% 
  summarise(Means = mean(Estimate), SD = sd(Estimate),
            "Average SE" = mean(`Std. Error`))


#Plot of main effect
plotdata %>% filter(Est %in% c("y", "trt"), 
              corcoef == .5, N %in% c(100, 500),
              vartype %in% "continuous") %>% 
  group_by(method) %>%
  ggplot(aes(Estimate)) +
  geom_density(col = "red", fill = "cadetblue") +
 # facet_grid(rows = vars(method), scales = "free") +
  facet_wrap(vars(N, method))


#Summary table for moderating effect
plotdata <- getR %>% dplyr::mutate_at(c("lwr", "upr"), func) %>% 
  tibble::remove_rownames() %>% as.data.frame() %>%
  mutate(Z = pchisq((Estimate/`Std. Error`)^2, df=1, lower.tail=FALSE),
         Pind = ifelse(Z<=0.05, 1, 0))
#Plot of main effect
plotdata %>% filter(Est %in% c("y:x1", "trt:x1"), 
                    corcoef == .5, N %in% c(100, 500),
                    vartype %in% "continuous") %>% 
  #group_by(method) %>%
  ggplot(aes(Estimate)) +
  geom_density(col = "red", fill = "cadetblue") +
  # facet_grid(rows = vars(method), scales = "free") +
  facet_wrap(vars(N, method))


plotdata %>% dplyr::mutate_at(c("lwr", "upr"), func) %>% 
  tibble::remove_rownames() %>% as.data.frame() %>%
  mutate(Z = pchisq((Estimate/`Std. Error`)^2, df=1, lower.tail=FALSE),
         Pind = ifelse(Z <= 0.05, 1, 0)) %>%
filter(Est %in% c("y:x1", "trt:x1"), corcoef == .5,
       N %in% c(100, 500),
       vartype %in% "continuous") %>%
  group_by(method, N) %>% summarise(
    Means = mean(Estimate), SD = sd(Estimate),
    "Average SE" = mean(`Std. Error`),
    POWER = mean(Pind))

# Z <- pchisq((BETA0/SE0)^2, df=1, lower.tail=FALSE)
# POWER <- matrix(apply(Z<=0.05, 2, FUN=mean), nrow=3, byrow=FALSE)





#=============================================================================
plotdata %>% filter(Est %in% c("trt", "y")) %>%  group_by(method) %>%
  ggplot(aes(Estimate, vartype)) +
  geom_boxplot() +
  facet_wrap(vars(method, N))

plotdata %>% filter(Est %in% c("trt", "y"), corcoef == .2) %>%  group_by(method) %>%
  ggplot(aes(`Std. Error`)) +
  geom_density() +
  facet_wrap(vars(method, N))


plotdata %>% filter(Est %in% c("trt:x1", "trt:x2", "y:x1", "y:x2"), corcoef == .5) %>%
  group_by(method) %>%
  ggplot(aes(`Std. Error`)) +
  geom_boxplot() +
  facet_wrap(vars(method, N))

plotdata %>% filter(Est %in% c("trt", "y")) %>%  group_by(method) %>%
ggplot(aes(`Std. Error`, Estimate)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  facet_wrap(vars(method, N))


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
        # geecoefs <- summary(getres)[c(2,12,13),  1:2]
        # inversecoefs <-summary(getres, model = 2)[c(2,7, 8),  1:2]
        # directcoefs <- summary(getres, model = 1)[c(2,7,8),  1:2]
        # geeconfint <- confint(getres, parm = c(2,12,13), level = .95, model = 3)[, -1]
        # inverseconfint<- confint(getres, parm = c(2,7,8), level = .95, model = 2)[, -1]
        # directconfint <- confint(getres, parm = c(2,7,8), level = .95, model = 1)[, -1]
        
        geecoefs <- summary(getres)[c(2,10:12),  1:2]
        inversecoefs <-summary(getres, model = 2)[c(2,6:8),  1:2]
        directcoefs <- summary(getres, model = 1)[c(2,6:8),  1:2]
        geeconfint <- confint(getres, parm = c(2,10:12), 
                              level = .95, model = 3)[, -1]
        inverseconfint<- confint(getres, parm = c(2,6:8), 
                                 level = .95, model = 2)[, -1]
        directconfint <- confint(getres, parm = c(2,6:8), 
                                 level = .95, model = 1)[, -1]
        out1 <- rbind(geecoefs, directcoefs, inversecoefs)
        out2 <-  rbind(cbind(geeconfint, method = "gee"),
                       cbind(directconfint, method = "direct"),
                       cbind(inverseconfint, method= "inverse"))
        out3 <- cbind(out1, out2, N = n, vartype = predictors,
                      corcoef = rho, Est = c(row.names(geecoefs),
                      row.names(directcoefs),row.names(inversecoefs)))
        OUT <- rbind(OUT, out3)
      }
    }
  }
rbind(result, OUT)
}




cl <- parallel::makeCluster(10)
parallel::clusterExport(cl, list("b0", "N", "vartype", "corr"))
parallel::clusterEvalQ(cl, library(Modana))
fullres <- do.call(rbind.data.frame, parallel::parLapply(cl, 1:500, simul))
parallel::stopCluster(cl)



#set.seed(111001)
set.seed(1001)
b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, 2, 1.5)
a0 <- c(-1, 1.5, -1, 1.5)
N = c(100, 180, 500, 1000)
vartype <- c("continuous", "binary")
corr <- c(0.2, 0.5, 0.8)
nrun <- 10
results <- NULL
for(i in 1:nrun){
  OUT <- NULL 
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
  results <- rbind(results, OUT) 
}









