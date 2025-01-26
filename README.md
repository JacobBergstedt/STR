
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
a11 <- .8; a12 <- .2; a22 <-.6;
c11 <- .1; c12 <- 0; c22 <- .05;
e11 <- sqrt(1-a11^2-c11^2); e12 <- .35; e22 <- sqrt(1-a22^2-c22^2-a12^2-c12^2-e12^2)

### The model matrices, using the path coefficients, which we later are going to fit
aa <- matrix(c(a11, a12, 0, a22), 2, 2)
cc <- matrix(c(c11, c12, 0, c22), 2, 2)
ee <- matrix(c(e11, e12, 0, e22), 2, 2)

### Modeled correlation within individual
aa%*%t(aa)+cc%*%t(cc)+ee%*%t(ee)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.3670628
#> [2,] 0.3670628 1.0000000

### Modeled correlation between MZ twins
aa%*%t(aa)+cc%*%t(cc)
#>      [,1]   [,2]
#> [1,] 0.65 0.1600
#> [2,] 0.16 0.4025

### Modeled correlation between DZ twins
0.5*aa%*%t(aa)+cc%*%t(cc)
#>      [,1]   [,2]
#> [1,] 0.33 0.0800
#> [2,] 0.08 0.2025

### Number of twin pairs
nMZ <- 10000
nDZ <- 12000

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
#>             [,1]       [,2]       [,3]       [,4]
#>  [1,] -1.0185179 -0.6690920 -0.5325441 -0.7209690
#>  [2,]  2.1771622  1.0055201  0.4592814 -0.5496150
#>  [3,] -0.8333055 -0.2668675  0.4836718  1.3318920
#>  [4,] -0.5079851  0.1206256  0.3462524  0.3872831
#>  [5,]  0.3057839  1.1489179  0.3007260 -1.9760223
#>  [6,]  1.0891023  1.9501608  0.2296970  0.1479644
#>  [7,] -1.2310997 -1.3357307  0.5809768  0.2061921
#>  [8,]  0.8496573  1.9012625  0.9330501  1.7933927
#>  [9,] -0.2545374  0.2710359 -0.6159503 -1.2651324
#> [10,] -0.0320923 -0.3720217  0.9100579 -0.6555562

### Generate "sex"
sexMZ <- data.frame(S1 = rbinom(nMZ, 1, .49),
                    S2 = 0 )
sexMZ$S2<-sexMZ$S1
sexDZ <- data.frame(S1 = rbinom(nDZ, 1, .49),
                    S2 = rbinom(nDZ, 1, .49))
sexMZ[1:10,]
#>    S1 S2
#> 1   1  1
#> 2   0  0
#> 3   0  0
#> 4   1  1
#> 5   0  0
#> 6   1  1
#> 7   1  1
#> 8   1  1
#> 9   0  0
#> 10  1  1

### Generate "birth year"
byMZ <- data.frame(BY1 = round(runif(nMZ,0,1),1), BY2 = 0 )
byMZ$BY2 <- byMZ$BY1
byDZ <- data.frame(BY1 = round(runif(nDZ,0,1),1), BY2 = 0 )
byDZ$BY2 <- byDZ$BY1
byMZ[1:10,]
#>    BY1 BY2
#> 1  0.4 0.4
#> 2  0.7 0.7
#> 3  0.0 0.0
#> 4  0.7 0.7
#> 5  0.1 0.1
#> 6  0.5 0.5
#> 7  0.8 0.8
#> 8  0.9 0.9
#> 9  0.5 0.5
#> 10 0.3 0.3

### Note: Birth year is the same across twins, but coded so it's allowed to be
###       different below, to allow for other variables.

### Chose threshold for variable 1 and variable 2
thr1 <- 0.8
thr2 <- 1.2

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



# -------------------------------------------------------------------------

# Select Variables for Analysis

prep <- prep_bivariate_data_non_expand(db = twinData,
                                       traitX = "X",
                                       traitY = "Y",
                                       response_typeX = "binary",
                                       response_typeY = "binary",
                                       same_sex = FALSE)

# m <- umxACE(selDVs = c("X", "Y"), selCovs = c("Birth_year_first", "Birth_year_second"), sep = "", dzData = prep$DZ, mzData = prep$MZ, tryHard = "yes")
# mv1 <- umxACEv(selDVs = c("X", "Y"), selCovs = c("age"), sep = "", dzData = prep$DZ, mzData = prep$MZ, tryHard = "yes")


m <- fit_ACE(prep, extra_tries = 1)
```

``` r
p <- get_estimates(m)
p
#> # A tibble: 147 x 4
#>    Param   Value       SD model  
#>    <chr>   <dbl>    <dbl> <chr>  
#>  1 cov   0.332   NaN      ACE_biv
#>  2 covA  0.115     0.0337 ACE_biv
#>  3 covC  0.00991   0.0310 ACE_biv
#>  4 covE  0.206     0.0111 ACE_biv
#>  5 bivA  0.348     0.116  ACE_biv
#>  6 bivC  0.0299    0.0925 ACE_biv
#>  7 bivE  0.622     0.0387 ACE_biv
#>  8 rg    0.258     0.0815 ACE_biv
#>  9 rc    1.00      0.0214 ACE_biv
#> 10 re    0.485     0.0294 ACE_biv
#> # i 137 more rows
```
