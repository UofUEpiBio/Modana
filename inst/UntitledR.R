#' #' @rdname refinedmod
#' #' @export
#' #' @param object an object class of \code{refinedmod}
#' #' @param model index indicating the summary call to print
#' #' @param ... other argument not in use at the moment 
#' summary.refinedmod <- function(object, model = NULL, ...){
#'   if(!is.null(model)){
#'     if(model==1){
#'       out <- capture.output(object[["allsummary.coef"]][["out.direct"]])
#'       v <- c("Model summary of the direct estimation:", out, "\n") # char vec
#'       s <- paste(v, collapse = "\n") # single string
#'     }else if(model == 2){
#'       out <- capture.output(object[["allsummary.coef"]][["out.inverse"]])
#'       v <- c("Model summary of the inverse estimation:", out, "\n") # char vec
#'       s <- paste(v, collapse = "\n") # single string
#'     }else{
#'       out <- capture.output(object[["allsummary.coef"]])
#'       v <- c("Model summary of the direct, inverse and inverse estimation:", out, "\n") # char vec
#'       s <- paste(v, collapse = "\n")}
#'   }else{
#'     out <- capture.output(object[["allsummary.coef"]][["out.gee"]])
#'     v <- c("Model summary of the refined GEE estimation:", out, "\n") # char vec
#'     s <- paste(v, collapse = "\n") # single string
#'   }
#'   # cat(s)
#'   structure(s, class = "summary.refinedmod")
#' }
#' 
#' 
#' #' @export 
#' print.summary.refinedmod <- function(x,...){
#'   cat(x)
#'   invisible(x)
#' }
#' 
#' #' #' @export 
#' #'  print.summary2.refinedmod <- function(x,...){
#' #'    cat(attr(x, "message"))
#' #'    print.data.frame(x)
#' #'    invisible(x)
#' #'  }
#' 



#```{r}
library(ggplot2)
library(dplyr)
# fn <- system.file(file.path("~/Modana/inst",
#   "array.slurm", "modana.slurm.R"),
#   package = "Modana"
#   )
maineff <- readRDS("~/Modana/inst/maineff.rds")
func <- \(x) as.numeric(as.character(x))
plotdata <- maineff %>% dplyr::mutate_at(c("lwr", "upr"), func) %>% 
  tibble::remove_rownames() %>% as.data.frame()

#Plot of main effect
plotdata %>% filter(Est %in% c("y", "trt"), 
                    corcoef == .5, N %in% c(100, 500),
                    vartype %in% "continuous") %>% 
  group_by(method) %>%
  ggplot(aes(Estimate)) +
  geom_density(col = "red", fill = "cadetblue") +
  # facet_grid(rows = vars(method), scales = "free") +
  facet_wrap(vars(N, method))

#```


#```{r}
#Summary table for Main effect
summtab <- plotdata %>% 
  filter(Est %in% c("trt", "y"), corcoef == .5,
         N %in% c(100, 500),
         vartype %in% "continuous") %>% 
  group_by(method, N) %>% 
  summarise(Means = mean(Estimate), SD = sd(Estimate),
            "Average SE" = mean(`Std. Error`))
knitr::kable(summtab)

#```




## Simulation Results (cont'd)


#```{r}
#| fig-cap: Empirical densities from simulation results for assessing the moderating effect  
#| warning: false

effmod <- readRDS("~/Modana/inst/effmod.rds")
func <- \(x) as.numeric(as.character(x))
plotdata <- effmod %>% dplyr::mutate_at(c("lwr", "upr"), func) %>% 
  tibble::remove_rownames() %>% as.data.frame()

#Plot of moderating effect
plotdata %>% filter(Est %in% c("y:x1", "trt:x1"), 
                    corcoef == .5, N %in% c(100, 500),
                    vartype %in% "continuous") %>% 
  #group_by(method) %>%
  ggplot(aes(Estimate)) +
  geom_density(col = "red", fill = "cadetblue") +
  # facet_grid(rows = vars(method), scales = "free") +
  facet_wrap(vars(N, method))

#```


#```{r}

#Summary table for Moderating effect
sumtab <- plotdata %>% dplyr::mutate_at(c("lwr", "upr"), func) %>% 
  tibble::remove_rownames() %>% as.data.frame() %>%
  mutate(Z = pchisq((Estimate/`Std. Error`)^2, df = 1, lower.tail=FALSE),
         Pind = ifelse(Z <= 0.05, 1, 0)) %>%
  filter(Est %in% c("y:x1", "trt:x1"), corcoef == .5,
         N %in% c(100, 1000),
         vartype %in% "continuous") %>%
  group_by(method, N) %>% summarise(
    Means = mean(Estimate), SD = sd(Estimate),
    "Average SE" = mean(`Std. Error`),
    # POWER = mean(Pind),
    lower = Means-qnorm((1+.95)/2)*`Average SE`,
    upper = Means+qnorm((1+.95)/2)*`Average SE`,
    conf_width = upper - lower)
knitr::kable(sumtab)
#```


## Results

#```{r}
#| fig-cap: Empirical confidence interval from simulation results for assessing the moderating effect 
plotdata %>% filter(Est %in% c("y:x1", "trt:x1"), corcoef == .5,
                    N %in% c(100, 500),
                    vartype %in% "continuous") %>%
  #group_by(method) %>%
  ggplot(aes(Estimate)) +        # ggplot2 plot with confidence intervals
  #geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr, col = "cadetblue")) +
  facet_wrap(vars(N, method))

#```

