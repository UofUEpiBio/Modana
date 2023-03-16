
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modana

<!-- badges: start -->
<!-- badges: end -->

The goal `Modana` package is to implement a refinement of moderation
analysis with binary outcomes, as proposed by [Anto and Su
(2023)](https://journals.sagepub.com/doi/abs/10.1177/09622802231151206?journalCode=smma).
The function fits three models of interest: a direct model, an inverse
model, and a generalized estimating equation (GEE) model. The direct and
inverse models are fitted using the `glm` function to verify the
symmetry property of odds ratio/relative risk in moderation analysis for
the main treatment effect as well as the moderating effects. The GEE
model is fitted using the `geeglm` function and is used to estimate the
treatment effect accounting for within-cluster correlation. The `Modana`
package has three different built functions.

## Installation

You can install the development version of Modana like so:

``` r
install.packages("Modana_0.1.0.tar.gz", repos = NULL)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(Modana)
## basic example code
```

You can simulate some data to run moderation analysis

``` r
set.seed(1099)
b0 <- c(1, 1.2, 1, 1.5, -2, -.5, 1.2, 1.5)
a0 <- c(-1, 1.5, -1, 1.5)
#a0 <- c(-2, 2, 2, 0, 1)
b0 <- c(-1.5, 1.5, 1.5, -2, 2, -.5, -1, 1.5)
 datt <- sim_data(n = 100, b0, a0 = NULL, binary.Xs = FALSE,
                      sigma = 1, uniform = FALSE, c0 = 1,
                      link.function = "logistic", rho = 0.2,
                      observational = FALSE, trt.p = 0.5,
                      interaction = 1:2, details = FALSE)
#> There are a total of 4 covariates
head(datt)
#>   y trt         x1         x2         x3         x4
#> 1 0   0 -1.0272485  1.5804971 -0.1492128  0.2379859
#> 2 1   0 -2.3314105 -1.6163845  1.9313538  1.8053311
#> 3 0   1  0.4191473 -0.7286691 -0.5131350  1.8232865
#> 4 0   0 -1.8775643  1.0827282  1.3265782 -1.1509396
#> 5 1   0  0.5414926 -0.1900738  0.4512489  0.4825833
#> 6 1   0 -0.6896955  0.6951233  1.9000859 -0.3240436
```

The `refinedmod()` function implements refined moderation analysis to
estimate both treatment effect and moderating effects associated with
two logistic regression models with binary outcomes/treatment via
generalized estimating equations (GEE). The function first performs a
role swap algorithm on the original data to obtain a clustered data of
size `2` (i.e., swapping the roles of the response variable and
treatment variable columns in the original data data and then stacking
the resultant data on top of the original data).

``` r
getres <- refinedmod(formula = y ~ trt + x1 + x2 + x3 + x4,
                     detail = TRUE, y = "y",
                    trt = "trt", data = datt,
                     effmod = c("x1", "x2"),
                      corstr = "independence")
#> The direct modeling fitting: 
#> =================================================================
#> 
#> Call:
#> glm(formula = form.direct, family = binomial(link = "logit"), 
#>     data = .SD)
#> 
#> Deviance Residuals: 
#>      Min        1Q    Median        3Q       Max  
#> -2.12526  -0.54951  -0.08492   0.53698   2.68620  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -2.1411     0.6469  -3.310 0.000934 ***
#> trt           1.9063     0.7249   2.630 0.008545 ** 
#> x1            1.6592     0.6746   2.459 0.013921 *  
#> x2           -2.3898     0.6830  -3.499 0.000467 ***
#> x3            2.4399     0.5453   4.474 7.67e-06 ***
#> x4           -0.3810     0.2994  -1.273 0.203182    
#> trt:x1       -1.6094     0.7477  -2.153 0.031357 *  
#> trt:x2        1.4938     0.7740   1.930 0.053598 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 136.058  on 99  degrees of freedom
#> Residual deviance:  71.649  on 92  degrees of freedom
#> AIC: 87.649
#> 
#> Number of Fisher Scoring iterations: 6
#> 
#> The inverse modeling fitting: 
#> =================================================================
#> 
#> Call:
#> glm(formula = form.inverse, family = binomial(link = "logit"), 
#>     data = .SD)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -1.8781  -0.9950   0.4503   1.0532   1.8161  
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)  
#> (Intercept)  -0.3216     0.3075  -1.046   0.2957  
#> y             1.6005     0.6344   2.523   0.0116 *
#> x1            0.1812     0.2571   0.705   0.4810  
#> x2           -0.2929     0.3610  -0.811   0.4172  
#> x3           -0.2798     0.3029  -0.924   0.3556  
#> x4            0.4367     0.2297   1.901   0.0573 .
#> y:x1         -0.6432     0.4373  -1.471   0.1414  
#> y:x2          1.2298     0.5382   2.285   0.0223 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 138.27  on 99  degrees of freedom
#> Residual deviance: 120.78  on 92  degrees of freedom
#> AIC: 136.78
#> 
#> Number of Fisher Scoring iterations: 4
#> 
#> The GEE modeling fitting: 
#> =================================================================
#>               Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept) -1.9882558  0.6313724 -3.149102 1.637732e-03
#> trt          1.7420214  0.6742747  2.583548 9.778975e-03
#> x1           1.1263252  0.3902621  2.886073 3.900812e-03
#> x2          -2.1327750  0.4347750 -4.905468 9.320477e-07
#> x3           2.2584618  0.4912020  4.597827 4.269210e-06
#> x4          -0.4006063  0.2991221 -1.339273 1.804817e-01
#> z            1.6338187  0.5822168  2.806203 5.012902e-03
#> x1:z        -0.8611383  0.3730072 -2.308637 2.096371e-02
#> x2:z         1.8597407  0.4994187  3.723811 1.962380e-04
#> x3:z        -2.5671182  0.6832798 -3.757053 1.719264e-04
#> x4:z         0.8363693  0.3954295  2.115091 3.442221e-02
#> trt:x1      -0.9294065  0.4424454 -2.100613 3.567497e-02
#> trt:x2       1.2861794  0.5361274  2.399018 1.643911e-02
names(getres)
#> [1] "Directmodel"       "Inversemodel"      "Geerefinedmodel"  
#> [4] "out.combined"      "comparisonresults"
```

You can call each elements of listed models

Mostly importantly, we can print the results of the refined estimation.

``` r
print(getres$Geerefinedmodel)
#>               Estimate Std. Error   z value     Pr(>|z|)
#> (Intercept) -1.9882558  0.6313724 -3.149102 1.637732e-03
#> trt          1.7420214  0.6742747  2.583548 9.778975e-03
#> x1           1.1263252  0.3902621  2.886073 3.900812e-03
#> x2          -2.1327750  0.4347750 -4.905468 9.320477e-07
#> x3           2.2584618  0.4912020  4.597827 4.269210e-06
#> x4          -0.4006063  0.2991221 -1.339273 1.804817e-01
#> z            1.6338187  0.5822168  2.806203 5.012902e-03
#> x1:z        -0.8611383  0.3730072 -2.308637 2.096371e-02
#> x2:z         1.8597407  0.4994187  3.723811 1.962380e-04
#> x3:z        -2.5671182  0.6832798 -3.757053 1.719264e-04
#> x4:z         0.8363693  0.3954295  2.115091 3.442221e-02
#> trt:x1      -0.9294065  0.4424454 -2.100613 3.567497e-02
#> trt:x2       1.2861794  0.5361274  2.399018 1.643911e-02
```

Example usage on a real-world data

``` r
#install.packages("CERFIT")
library(CERFIT)
data(warts, package =  "CERFIT")
dat <- warts
#data wrangling
dat$Type <- ifelse(dat$Type==3, 1, 0)
dat$response <- dat$Result_of_Treatment
```

``` r
res <- refinedmod(formula = response ~ treatment + age + Time + Type,
                     detail = FALSE, y = "response",
                    trt = "treatment", data = dat,
                     effmod = c("age", "Type"),
                      corstr = "independence")
print(res$comparisonresults)
#> $Direct.model
#> 
#> Call:
#> glm(formula = form.direct, family = binomial(link = "logit"), 
#>     data = .SD)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -3.4483  -0.2552   0.2177   0.5278   1.8712  
#> 
#> Coefficients:
#>                Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     7.43389    1.44174   5.156 2.52e-07 ***
#> treatment       0.59584    1.31994   0.451   0.6517    
#> age            -0.03101    0.02545  -1.219   0.2230    
#> Time           -0.58627    0.11036  -5.312 1.08e-07 ***
#> Type           -0.54809    0.87790  -0.624   0.5324    
#> treatment:age  -0.05274    0.04053  -1.301   0.1931    
#> treatment:Type -2.84916    1.43163  -1.990   0.0466 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 230.51  on 179  degrees of freedom
#> Residual deviance: 119.85  on 173  degrees of freedom
#> AIC: 133.85
#> 
#> Number of Fisher Scoring iterations: 6
#> 
#> 
#> $Inverse.model
#> 
#> Call:
#> glm(formula = form.inverse, family = binomial(link = "logit"), 
#>     data = .SD)
#> 
#> Deviance Residuals: 
#>     Min       1Q   Median       3Q      Max  
#> -2.2753  -1.0578   0.1139   1.0293   1.9538  
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    1.33013    1.12689   1.180 0.237861    
#> response       1.12250    1.10380   1.017 0.309181    
#> age           -0.00712    0.02312  -0.308 0.758123    
#> Time          -0.09123    0.06679  -1.366 0.171959    
#> Type           1.87854    0.70810   2.653 0.007980 ** 
#> response:age  -0.07198    0.03264  -2.205 0.027433 *  
#> response:Type -3.63225    0.95771  -3.793 0.000149 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 249.53  on 179  degrees of freedom
#> Residual deviance: 202.03  on 173  degrees of freedom
#> AIC: 216.03
#> 
#> Number of Fisher Scoring iterations: 4
#> 
#> 
#> $Combinedgee.model
#>                   Estimate Std. Error    z value     Pr(>|z|)
#> (Intercept)     7.42036779 2.12457045  3.4926438 4.782640e-04
#> treatment       0.94854255 1.09332337  0.8675773 3.856258e-01
#> age            -0.02664027 0.02350883 -1.1332029 2.571291e-01
#> Time           -0.60056451 0.16830032 -3.5684098 3.591544e-04
#> Type           -0.35357228 0.91818191 -0.3850787 7.001791e-01
#> z              -6.07463142 1.39817721 -4.3446792 1.394794e-05
#> age:z           0.01717663 0.02031301  0.8455974 3.977774e-01
#> Time:z          0.51570766 0.13299353  3.8776897 1.054531e-04
#> Type:z          2.13458124 0.65420080  3.2628839 1.102847e-03
#> treatment:age  -0.06485767 0.03032095 -2.1390378 3.243260e-02
#> treatment:Type -3.40781554 0.91201807 -3.7365658 1.865506e-04
```
