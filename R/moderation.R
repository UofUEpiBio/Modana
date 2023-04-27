#' @title Refined moderation Analysis with binary outcomes
#' @description This package employs the symmetry property concerning the ratio of odds ratios, which suggests that heterogeneous treatment effects could
#' be equivalently estimated via a role exchange between the outcome and treatment variable in logistic regression
#' models. Refined inference on moderating effects is obtained by rearranging data and combing two models into
#' one via a generalized estimating equation (GEE) approach.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param y  name of a binary response column in data.
#' @param trt name of a binary treatment column in data.
#' @param effmod a vector of potential effect modifiers used in the data.
#' @param detail a logical argument specified to print summary results of the direct and inverse models
#' @param ...  Other arguments, which are passed to the functions in \code{\link[geepack:geeglm]{geeglm}}.
#' @details Three models of interest are fitted. The generalized linear models \code{\link[stats:glm]{glm}} (for the direct and inverse models) are fitted to verify the symmetry property of OR/RR in moderation analysis for the main treatment effect as well as the moderating effects.
#' The main idea of this computing trick is to rearrange data into a format that micmicks a clustered or longitudinal data. A clustered data of size `2` is obtained simply by swapping the roles of the response variable \code{y} and treatment variable \code{trt}
#' columns in the original data \code{data} and then stacking the resultant data on top of the original \code{data}. The dummy variable \code{z}  is introduced to distinguish the two glm models.
#' @param data original data frame containing the variables in the model.
#' @return An object of type \code{glm} and \code{geeglm}.
#' @references Anto E, Su X. Refined moderation analysis with binary outcomes in precision medicine research. Statistical Methods in Medical Research. 2023;32(4):732-747. doi:10.1177/09622802231151206.
#' @examples
#' set.seed(1001)
#' b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, 2, 1.5)
#' a0 <- c(-1, 1.5, -1, 1.5)
#' datt <- sim_data(n = 100, b0 = b0, a0 = NULL, binary.Xs = FALSE,
#'                   sigma = 1, uniform = FALSE, c0 = 1,
#'                    link.function = "logistic", rho = 0.2,
#'                    observational = FALSE, trt.p = 0.5,
#'                    interaction = 1:2, details = FALSE)
#'
#' getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3 + x4,
#'                      detail = FALSE, y = "y",
#'                     trt = "trt", data = datt,
#'                      effmod = c("x1", "x2"),
#'                       corstr = "independence")
#' names(getres)
#' print(getres)
#' dsu <- summary(getres, 1); dsu
#' summary(getres, model = NULL)
#' plot(getres)
#' confint(getres, model = 3)
#' @export
refinedmod <- function(formula, y = "y", trt = "trt",
                       effmod = NULL, data, detail = FALSE, ...){
  dat0 <- swaparr(data, y, trt)
  dat0 <- data.frame(dat0[order(dat0$ID), ], row.names = NULL)
  ID = factor(dat0$ID)
  data0 <- data.table::setDT(dat0)
  data1 <- data.table::setDT(data)
  form <- formula
  vars <- all.vars(form)
  covar <- vars[!vars %in% c(y, trt)]
  covars <- c(trt, covar)
  covars.inv <- c(y, covar)
  if(is.null(effmod)){
    form.gee <- as.formula(paste(y, "~",  trt, "+(",
                                 paste0(covar, sep = "", collapse = "+"),")*z"))
    form.inverse <- as.formula(paste(trt, "~",
                                     paste0(covars.inv, sep = "", collapse = "+")))
    form.direct <- as.formula(paste(y, "~",
                                    paste0(covars, sep = "", collapse = "+")))
  } else{
    form.gee <- as.formula(paste(y, "~",  trt, "+(",
                                 paste0(covar, sep = "", collapse = "+"),")*z +(",
                                 paste0(effmod, sep = "", collapse = "+"),"):",trt))

    form.inverse <- as.formula(paste(trt, "~",
                                     paste0(covars.inv, sep = "", collapse = "+"), "+(",
                                     paste0(effmod, sep = "", collapse = "+"),"):", y))

    form.direct <- as.formula(paste(y, "~",
                                    paste0(covars, sep = "", collapse = "+"), "+(",
                                    paste0(effmod, sep = "", collapse = "+"),"):", trt))
  }
  #Fitting the models
  fit.gee <- data0[, geepack::geeglm(form.gee, family = binomial(link = "logit"),
                          data = .SD, id = ID, ...)]
  fit.direct <- data1[, glm(form.direct, family =  binomial(link = "logit"), data = .SD)]
  fit.inverse <- data1[, glm(form.inverse, family = binomial(link = "logit"), data = .SD)]
  # fit.gee <- geepack::geeglm(formula = formula0, data = dat0, id = ID,
  #                            family = binomial, ...)
  # fit.direct <- stats::glm(form.direct, data = data,
  #                          family = binomial(link = "logit"))
  # fit.inverse <- stats::glm(form.inverse, data = data,
  #                    family = binomial(link = "logit"))
  if (detail){
    cat("The direct modeling fitting: \n")
    cat("=================================================================\n")
    print(summary(fit.direct))
  }

  if (detail){
    cat("The inverse modeling fitting: \n")
    cat("=================================================================\n")
    print(summary(fit.inverse))
  }
  out.direct <- summary(fit.direct)$"coefficients"
  out.inverse <- summary(fit.inverse)$"coefficients"
  out.gee <- summary(fit.gee)$coefficients
  #names(result.gee)
  out.gee$Wald <- (out.gee$Estimate)/(out.gee$Std.err )
  out.gee$'Pr(>|W|)' <- 2*pnorm(-abs(out.gee$Wald ))

  names(out.gee)[names(out.gee) == "Wald"] <- "z value"
  names(out.gee)[names(out.gee) == "Pr(>|W|)"] <- "Pr(>|z|)"
  names(out.gee)[names(out.gee) == "Std.err"] <- "Std. Error"
  if (detail){
    cat("The GEE modeling fitting: \n")
    cat("=================================================================\n")
    print(out.gee)
  }
  # alres <- list(Direct.model = summary(fit.direct),
  #               Inverse.model = summary(fit.inverse),
  #               Combinedgee.model = out.gee)
  
  pred.direct<- as.vector(stats::predict(fit.direct, type = "response"))
  pred.inverse <- as.vector(stats::predict(fit.inverse, type = "response"))
  pred.gee <- as.vector(stats::predict(fit.gee, type = "response"))
  resid.direct <- as.vector(fit.direct$residuals)
  resid.inverse <- as.vector(fit.inverse$residuals)
  resid.gee <- as.vector(fit.gee$residuals)
  pltdata <- rbind.data.frame(cbind(pred = pred.gee,
                              residual = resid.gee,
                              model = "gee"),
                              cbind(pred = pred.direct,
                              residual = resid.direct,
                              model = "direct"),
                              cbind(pred = pred.inverse,
                              residual = resid.inverse,
                              model = "inverse")
                              )
  pltdata$pred<-as.numeric(as.character(pltdata$pred))
  pltdata$residual<-as.numeric(as.character(pltdata$residual))
  structure(list(
    allsummary.coef = 
    list(out.direct = out.direct,
          out.inverse = out.inverse,
          out.gee   = out.gee),
    Directmodel  = fit.direct,
    Inversemodel = fit.inverse,
    Refinedmodel = fit.gee,
    pltdata = pltdata),
    class = "refinedmod")
}


#' @rdname refinedmod
#' @param data original data frame containing the variables in the model.
#' @param y  name of a binary response column in data.
#' @param trt name of a binary treatment column in data.
#' @details swaps the columns specified by \code{y} and \code{trt}, appends the swapped data frame to the original data frame \code{data}, and adds two new columns \code{(group and ID)} to the merged data frame
#' @return  New data frame created by swapping selected columns \code{y} and \code{trt} from \code{data} using the sequence of indices constructed mimicking a longitudinal with 2 clusters.
#' @export

swaparr <- function(data, y = "y", trt = "trt")
{
  y <- if (class(y)=="character" & is.na(suppressWarnings(as.integer(y))))
    which(colnames(data)==y) else as.integer(y)
  trt <- if (class(trt)=="character" & is.na(suppressWarnings(as.integer(trt))))
    which(colnames(data)==trt) else as.integer(trt)

  if (!(1<=y & y<=length(data))) stop( "`y` represents invalid index!" )
  if (!(1<=trt & trt<=length(data))) stop("`trt` represents invalid index!")
  DF <- data[if(y==trt) 1:length(data)
             else c((if (min(y, trt)==1) c()
                     else 1:(min(y, trt)-1) ),
                    (if (min(y,trt)+1 == max(y,trt))
                      (min(y,trt)+1):(max(y,trt)-1)
                     else c( max(y,trt), (min(y,trt)+1):(max(y,trt)-1),
                             min(y,trt))), (if (max(y,trt)==length(data))
                               c() else (max(y,trt)+1):length(data)))]
  n <- nrow(data)
  names(DF) <- names(data)
  Dat <- data.frame(rbind(data, DF))
  Dat$z <- rep(c(0, 1), c(n, n))
  Dat$ID <- rep(1:n, 2)
  return(Dat)
}


#is.even <- function(x) x %% 2 == 0

#expit <- function(x) (tanh(x/2)+1)/2  	# CORRECT & VERIFIED
# expit0 <- function(x) exp(x)/(1+exp(x))  # LOGISTIC FUNCTION, INVERSE OF LOGIT FUNCTION
#' @rdname refinedmod
#' @details The function \code{sim_data()} generates simulated data for a randomized controlled trial or an observational study
#' @param n the number of observations to generate.
#' @param b0 a vector of coefficients for the linear predictor in the outcome model. The length of \code{b0} should be equal to the number of covariates plus one (for the intercept term). The first element of \code{b0} corresponds to the intercept, and the remaining elements correspond to the coefficients for the covariates which may include include interaction coefficients.
#' @param a0 a vector of coefficients for the linear predictor in the treatment assignment model. The length of \code{a0} should be equal to the number of covariates plus one (for the intercept term). The first element of a0 corresponds to the intercept, and the remaining elements correspond to the coefficients for the covariates. This argument is only used if `observational = TRUE`.
#' @param a0 a vector of coefficients for the linear predictor in the treatment assignment model. The length of \code{a0} should be equal to the number of covariates plus one (for the intercept term). The first element of a0 corresponds to the intercept, and the remaining elements correspond to the coefficients for the covariates. This argument is only used if `observational = TRUE`.
#' @param link.function a character string specifying the link function to use in the outcome model. The default is `logistic`, which corresponds to the logistic link function. The other option is `exp`, which corresponds to the log-linear link function.
#' @param sigma the standard deviation of the covariate distribution. The default value is `1`.
#' @param uniform a logical value indicating whether to transform the covariates to a uniform distribution on the interval `[0, c0]`. The default value is `FALSE`.
#' @param c0 the upper limit of the uniform distribution used to transform the covariates. The default value is 1.
#' @param binary.Xs a logical value indicating whether to change some columns of continuous covariates into binary variables. The default value is `FALSE`.
#' @param rho the correlation parameter used to generate the covariance matrix. The default value is `0.2`.
#' @param observational a logical value indicating whether to generate data for an observational study. The default value is FALSE.
#' @param trt.p the probability of receiving the treatment. This argument is only used if `observational = FALSE`. The default value is 0.5.
#' @param interaction a vector of indices indicating which covariates should be used to create interaction terms with the treatment variable. The length of interaction should be less than or equal to the number of covariates. The default value is `NULL`, which indicates no interaction terms should be used.
#' @param details a logical value indicating whether to print details about the generated data. The default value is FALSE.
#' @export
sim_data <- function(n = 100, b0, a0 = NULL, link.function = "logistic",
                     sigma=1, uniform=TRUE, c0 = 1, binary.Xs=FALSE,
                     rho = 0.2, observational = FALSE, trt.p = 0.5,
                     interaction = NULL, details = F) {

  # Check b0 length is even
  #if (!is.even(length(b0)) || length(b0) < 2) stop("Length of b0 should be an even number!")
  #p <- (length(b0) - 2) / 2
  # Number of covariates
  if (length(b0) < 2 || !is.numeric(b0)) stop("Invalid b0 argument")
  p <- (length(b0) - 2)
  if(!is.null(interaction)){
    p <- length(b0) - (2 + length(interaction))
  }
  if (!is.null(a0) && (length(a0) != (p + 1) || !is.numeric(a0))) stop("Invalid a0 argument")
  if (!is.null(interaction) &&
      (length(interaction) > p || !is.numeric(interaction))) stop("Invalid interactions argument")
  cat("There are a total of", p, "covariates \n")
  # Generate covariance matrix
  S <- sapply(1:p, function(i){
    sapply(1:p, function(j){
      sigma^2*rho^(abs(i-j))})})

  # Generate X
  #X <- draw.d.variate.uniform(no.row = n, d = p, cov.mat = S)
  mu <- rep(0, p)
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = S)
  if (uniform) X <- pnorm(X)*c0 	# PROB INTEGRAL TRANSFORM TO UNIFORM [0,c0]
  # CHANGE SOME COLUMNS OF X INTO BINARY
  if(binary.Xs){
    threshold <- ifelse(uniform, c0/2, 0);
    X <- (X > threshold)
  }
  # Treatment assignment
  if (observational) {
    if (is.null(a0)) stop("Observational study requires specifying a0.")
    if (length(a0) != (p + 1)) stop("Length of a0 must be consistent with b0 (i.e. length of a0 should be equal to ", paste(sum(p+1)))
    eta.trt <- as.vector(cbind(1, X)%*% a0)
    pi.trt <- (tanh(eta.trt/2)+1)/2
    trt <- stats::rbinom(n = n, size = 1, prob = pi.trt)
    if (details) print(cbind(p.treatment = mean(pi.trt), p.trt = mean(trt)))
  } else {
    trt <- stats::rbinom(n = n, size = 1, prob = trt.p)
  }

  # Model matrix
  X0 <- cbind(1, trt, X)
  if (!is.null(interaction)) {
    if (length(interaction) != length(unique(interaction))) stop("Duplicate entries in interactions argument")
    if (any(interaction < 1) || any(interaction > p)) stop("Invalid entries in interactions argument")
    if (length(interaction) > p) stop("Number of interaction terms should be less than or equal to the number of covariates.", "In your case,", paste("the number of covariates is ", p, "which is less than the number of interaction terms specified which is", length(interaction), ". Consider setting the length of the vector for the `interaction` argument to at most", p))
    X_int <- X[, c(interaction), drop = FALSE]
    X0 <- cbind(X0, X_int*trt)
  }

  # Model response
  eta <- as.vector(X0 %*% b0)
  if (details) print(cbind(eta.mean = mean(eta), eta.min = max(eta),
                           prop.positive = sum(eta > 0)/length(eta)))
  if (link.function == "exp") {
    pi0 <- exp(eta * (eta <= 0))  # LOGLINEAR
  } else {
    pi0 <-  (tanh(eta/2)+1)/2 #expit(eta)  # LOGISTIC
  }
  y <- stats::rbinom(n, size = 1, prob = as.vector(pi0))

  # Return data
  dat <- data.frame(cbind(y, trt, X))
  colnames(dat) <- c("y", "trt", paste0("x", 1:p))
  if (!is.null(interaction)) {
    datint <- cbind(dat, X_int*trt)
    colnames(datint)[(p+3):(p+(2+length(interaction)))] <- paste0("z", 1:length(interaction))
    if(details) print(head(datint))
  }
  dat
}

#' @rdname refinedmod
#' @export
#' @param object an object class of \code{refinedmod}
#' @param model index indicating the summary call to print
#' @param ... other argument not in use at the moment 
summary.refinedmod <- function(object, model = NULL, ...){
                if(!is.null(model)){
                if(model==1){
                out <- capture.output(object[["allsummary.coef"]][["out.direct"]])
                v <- c("Model summary of the direct estimation:", out, "\n") # char vec
                s <- paste(v, collapse = "\n") # single string
                }else if(model == 2){
                  out <- capture.output(object[["allsummary.coef"]][["out.inverse"]])
                  v <- c("Model summary of the inverse estimation:", out, "\n") # char vec
                  s <- paste(v, collapse = "\n") # single string
                }else{
                  out <- capture.output(object[["allsummary.coef"]])
                  v <- c("Model summary of the direct, inverse and inverse estimation:", out, "\n") # char vec
                  s <- paste(v, collapse = "\n")}
                }else{
                  out <- capture.output(object[["allsummary.coef"]][["out.gee"]])
                  v <- c("Model summary of the refined GEE estimation:", out, "\n") # char vec
                  s <- paste(v, collapse = "\n") # single string
                }
             # cat(s)
            structure(s, class = "summary.refinedmod")
}

#' @export 
 print.summary.refinedmod <- function(x,...){
   cat(x)
  invisible(x)
 }

#' @rdname refinedmod
#' @export
#' @param x an object class of \code{refinedmod}
#' @param ... other argument not in use at the moment 
print.refinedmod <- function(x, model = 1, ...){
  if(model==1){
    out <- coef(summary(x[["Directmodel"]]))
  }
  else if(model==2){
    out <- coef(summary(x[["Inversemodel"]]))
  }else{
    out <- summary(x[["Refinedmodel"]])$coefficients
    out$Wald <- (out$Estimate)/(out$Std.err)
    out$'Pr(>|W|)' <- 2*pnorm(-abs(out$Wald))
    names(out)[names(out) == "Wald"] <- "z value"
    names(out)[names(out) == "Pr(>|W|)"] <- "Pr(>|z|)"
    names(out)[names(out) == "Std.err"]  <- "Std. Error"
    
}
    return(out)
  #invisible(x)
}

#' @rdname refinedmod
#' @param x an object class of \code{refinedmod}
#' @param ... Other arguments, which are passed to the functions in \code{ggplot2}.
##@param y not in use at the moment
#' @export
plot.refinedmod <- function(x, ...){
  plotdata <- x[["pltdata"]]
 pl <- ggplot2::ggplot(plotdata, ggplot2::aes(pred, residual)) + 
   ggplot2::labs(x = "Fitted values", y = "Residuals") +
    ggplot2::geom_point() + ggplot2::facet_wrap(~model)
 pl
}

#coef.refinedmod <- function(object, model = 1, ...){
  
#}

#' @rdname refinedmod
#' @export
#' @param object an object fitted model of class of \code{refinedmod}
#' @param model index indicating the summary call to print
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level the confidence level required.
#' @param ... other argument not in use at the moment 
confint.refinedmod <- function(object, parm, level = 0.90, model = 1, ...) {
  
  #if(any(colnames(summary(x)$coefficients)=="Wald"))
  if(model==1){
  out <- coef(summary(object[["Directmodel"]]))
  }
  else if(model==2){
    out <- coef(summary(object[["Inversemodel"]]))
  }else{
    out <- summary(object[["Refinedmodel"]])$coefficients
    out$Wald <- (out$Estimate)/(out$Std.err)
    out$'Pr(>|W|)' <- 2*pnorm(-abs(out$Wald))
    names(out)[names(out) == "Wald"] <- "z value"
    names(out)[names(out) == "Pr(>|W|)"] <- "Pr(>|z|)"
    names(out)[names(out) == "Std.err"]  <- "Std. Error"
  }
  mult <- qnorm((1+level)/2)
  citab <- with(as.data.frame(out),
                cbind(Coef = Estimate, lwr=Estimate-mult*`Std. Error`,
                      upr=Estimate+mult*`Std. Error`))
  rownames(citab) <- rownames(out)
  citab[parm, ]
}
