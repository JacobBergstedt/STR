
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STR

<!-- badges: start -->
<!-- badges: end -->

The goal of STR is to provide a tool-set for fitting and evaluating
classical twin models.

## Installation

You can install the development version of STR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("JacobBergstedt/STR")
```

## Example

Here I first simulate a bivariate binary classical twin model with two
covariates, age and sex. Then I fit a model on the simulated data and
extract the parameters.

``` r
library(STR)
library(tidyverse)


# Sim data ----------------------------------------------------------------

### Choose parameter values, as "Cholesky" (axx^2+cxx^2 must be less than 1)
a11 <- .6; a12 <- .2; a22 <-.5;
c11 <- .5; c12 <- 0.2; c22 <- 0.6;
e11 <- sqrt(1-a11^2-c11^2); e12 <- .35; e22 <- sqrt(1-a22^2-c22^2-a12^2-c12^2-e12^2)

### The model matrices, using the path coefficients, which we later are going to fit
aa <- matrix(c(a11, a12, 0, a22), 2, 2)
cc <- matrix(c(c11, c12, 0, c22), 2, 2)
ee <- matrix(c(e11, e12, 0, e22), 2, 2)

### Modeled correlation within individual
aa%*%t(aa)+cc%*%t(cc)+ee%*%t(ee)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.4385749
#> [2,] 0.4385749 1.0000000

### Modeled correlation between MZ twins
aa%*%t(aa)+cc%*%t(cc)
#>      [,1] [,2]
#> [1,] 0.61 0.22
#> [2,] 0.22 0.69

### Modeled correlation between DZ twins
0.5*aa%*%t(aa)+cc%*%t(cc)
#>      [,1]  [,2]
#> [1,] 0.43 0.160
#> [2,] 0.16 0.545

### Number of twin pairs
nMZ <- 5000
nDZ <- 6000

### Seed
set.seed(9876)

### Generate the data
datMZ <- MASS::mvrnorm(nMZ,c(0,0,0,0),
                 rbind( cbind(aa%*%t(aa)+cc%*%t(cc)+ee%*%t(ee) , aa%*%t(aa)+cc%*%t(cc)  ),
                        cbind( aa%*%t(aa)+cc%*%t(cc)        , aa%*%t(aa)+cc%*%t(cc)+ee%*%t(ee))) )
datDZ <- MASS::mvrnorm(nDZ,c(0,0,0,0),
                 rbind(cbind(aa%*%t(aa)+cc%*%t(cc)+ee%*%t(ee) , .5*aa%*%t(aa)+cc%*%t(cc)),
                       cbind( .5*aa%*%t(aa)+cc%*%t(cc)       , aa%*%t(aa)+cc%*%t(cc)+ee%*%t(ee))) )

datMZ[ 1: 10 , ]
#>              [,1]       [,2]       [,3]        [,4]
#>  [1,]  0.50956089  1.0292550  1.0914929  0.48137454
#>  [2,] -1.73719313 -0.7429699 -1.0481252 -0.04576988
#>  [3,]  0.38929756 -0.1656650 -0.7266638 -0.07911859
#>  [4,] -0.74208084  0.8600272 -0.7278043  0.24944315
#>  [5,] -0.01933607 -0.5000609  0.4757190  0.12617020
#>  [6,] -0.28361857 -1.0325298 -1.0187149 -1.13326146
#>  [7,]  0.87166941 -0.2268021  0.8816732  0.33442157
#>  [8,] -0.39185237 -1.3733630 -1.3790755 -2.31824631
#>  [9,]  0.48187154  0.3586654  0.5668263  0.54068732
#> [10,] -0.36129294 -0.4010857  0.4050917  0.28739268

### Generate "sex"
sexMZ <- data.frame(S1 = rbinom(nMZ, 1, .49),
                    S2 = 0 )
sexMZ$S2<-sexMZ$S1
sexDZ <- data.frame(S1 = rbinom(nDZ, 1, .49),
                    S2 = rbinom(nDZ, 1, .49))
sexMZ[1:10,]
#>    S1 S2
#> 1   0  0
#> 2   1  1
#> 3   0  0
#> 4   1  1
#> 5   1  1
#> 6   1  1
#> 7   1  1
#> 8   0  0
#> 9   1  1
#> 10  1  1

### Generate "birth year"
byMZ <- data.frame(BY1 = round(runif(nMZ,0,1),1), BY2 = 0 )
byMZ$BY2 <- byMZ$BY1
byDZ <- data.frame(BY1 = round(runif(nDZ,0,1),1), BY2 = 0 )
byDZ$BY2 <- byDZ$BY1
byMZ[1:10,]
#>    BY1 BY2
#> 1  0.4 0.4
#> 2  0.6 0.6
#> 3  0.8 0.8
#> 4  0.2 0.2
#> 5  0.0 0.0
#> 6  0.4 0.4
#> 7  0.6 0.6
#> 8  0.7 0.7
#> 9  0.5 0.5
#> 10 0.8 0.8

### Note: Birth year is the samer across twins, but coded so it's allowed to be
###       different below, to allow for other variables.

### Chose threshold for variable 1 and variable 2
thr1 <- 1
thr2 <- 1.6

### Chose regression coefficients for S on X and Y
betaX <- 0
betaY <- 0

### Chose regression coefficients for BY on X and Y
betaX2 <- 0
betaY2 <- 0

### Binarize the data
datMZb <- cbind(
  1*( datMZ[,1]+betaX*sexMZ$S1+betaX2*byMZ$BY1 > thr1 ) ,
  1*( datMZ[,2]+betaY*sexMZ$S1+betaY2*byMZ$BY1 > thr2 ) ,
  1*( datMZ[,3]+betaX*sexMZ$S2+betaX2*byMZ$BY2 > thr1 ) ,
  1*( datMZ[,4]+betaY*sexMZ$S2+betaY2*byMZ$BY2 > thr2 ) )

datDZb <- cbind(
  1*( datDZ[,1]+betaX*sexDZ$S1+betaX2*byDZ$BY1 > thr1 ) ,
  1*( datDZ[,2]+betaY*sexDZ$S1+betaY2*byDZ$BY1 > thr2 ) ,
  1*( datDZ[,3]+betaX*sexDZ$S2+betaX2*byDZ$BY2 > thr1 ) ,
  1*( datDZ[,4]+betaY*sexDZ$S2+betaY2*byDZ$BY2 > thr2 ) )



datMZ1 <- tibble(X = 1*( datMZ[,1]+betaX*sexMZ$S1+betaX2*byMZ$BY1 > thr1 ),
                 Y = 1*( datMZ[,2]+betaY*sexMZ$S1+betaY2*byMZ$BY1 > thr2 ),
                 Female = sexMZ$S1,
                 b_year = byMZ$BY1,
                 age = 2 * b_year,
                 pairnnr = 1:nMZ,
                 twinnr = paste0(pairnnr, 1),
                 Zyg = "MZ")

datMZ2 <- tibble(X = 1*( datMZ[,3]+betaX*sexMZ$S2+betaX2*byMZ$BY2 > thr1 ),
                 Y = 1*( datMZ[,4]+betaY*sexMZ$S2+betaY2*byMZ$BY2 > thr2 ),
                 Female = sexMZ$S2,
                 b_year = byMZ$BY2,
                 age = 2 * b_year,
                 pairnnr = 1:nMZ,
                 twinnr = paste0(pairnnr, 2),
                 Zyg = "MZ")

datDZ1 <- tibble(X = 1*( datDZ[,1]+betaX*sexDZ$S1+betaX2*byDZ$BY1 > thr1 ),
                 Y = 1*( datDZ[,2]+betaY*sexDZ$S1+betaY2*byDZ$BY1 > thr2 ),
                 Female = sexDZ$S1,
                 b_year = byDZ$BY1,
                 pairnnr = (nMZ + 1):(nMZ + nDZ),
                 twinnr = paste0(pairnnr, 1),
                 Zyg = "DZ")

datDZ2 <- tibble(X = 1*( datDZ[,3]+betaX*sexDZ$S2+betaX2 * byDZ$BY2 > thr1 ),
                 Y = 1*( datDZ[,4]+betaY*sexDZ$S2+betaY2 * byDZ$BY2 > thr2 ),
                 Female = sexDZ$S2,
                 b_year = byDZ$BY2,
                 pairnnr = (nMZ + 1):(nMZ + nDZ),
                 twinnr = paste0(pairnnr, 2),
                 Zyg = "DZ")

datDZ <- bind_rows(datDZ1, datDZ2) |>
  group_by(pairnnr) |>
  mutate(Zyg = if_else(Female[1] == Female[2], "DZ_same_sex", "DZ_diff_sex")) |>
  mutate(age = 2 * b_year) |>
  ungroup()

twinData <- bind_rows(datMZ1, datMZ2, datDZ) |>
  mutate(X = factor(X, levels = c(0, 1)),
         Y = factor(Y, levels = c(0, 1)))


twinData
#> # A tibble: 22,000 x 8
#>    X     Y     Female b_year   age pairnnr twinnr Zyg  
#>    <fct> <fct>  <int>  <dbl> <dbl>   <int> <chr>  <chr>
#>  1 0     0          0    0.4   0.8       1 11     MZ   
#>  2 0     0          1    0.6   1.2       2 21     MZ   
#>  3 0     0          0    0.8   1.6       3 31     MZ   
#>  4 0     0          1    0.2   0.4       4 41     MZ   
#>  5 0     0          1    0     0         5 51     MZ   
#>  6 0     0          1    0.4   0.8       6 61     MZ   
#>  7 0     0          1    0.6   1.2       7 71     MZ   
#>  8 0     0          0    0.7   1.4       8 81     MZ   
#>  9 0     0          1    0.5   1         9 91     MZ   
#> 10 0     0          1    0.8   1.6      10 101    MZ   
#> # i 21,990 more rows

# -------------------------------------------------------------------------

# Select Variables for Analysis

prep <- prep_bivariate_data_non_expand(db = twinData,
                                       traitX = "X",
                                       traitY = "Y",
                                       covs = c("Female", "age"),
                                       response_typeX = "binary",
                                       response_typeY = "binary",
                                       same_sex = FALSE)

m <- fit_ACE(prep, covs = "Female", extra_tries = 1)
```

``` r
get_estimates(m)
#> # A tibble: 147 x 4
#>    Param Value     SD model  
#>    <chr> <dbl>  <dbl> <chr>  
#>  1 cov   1.25  1.20   ACE_biv
#>  2 covA  0.346 0.370  ACE_biv
#>  3 covC  0.342 0.360  ACE_biv
#>  4 covE  0.562 0.540  ACE_biv
#>  5 bivA  0.277 0.139  ACE_biv
#>  6 bivC  0.274 0.110  ACE_biv
#>  7 bivE  0.450 0.0460 ACE_biv
#>  8 rg    0.499 0.225  ACE_biv
#>  9 rc    0.350 0.132  ACE_biv
#> 10 re    0.596 0.0544 ACE_biv
#> # i 137 more rows
```
