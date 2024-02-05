
#' @export
fit_saturated <- function(x, ...) {

  UseMethod("fit_saturated")

}


#' @export
fit_saturated.prep.uni.bin <- function(x, ...) {

  MZ <- x$MZ
  DZ <- x$DZ

  nv <- 1 # number of variables
  ntv <- nv * 2 # number of total variables
  selVars <- paste("X", c(rep(1, nv), rep(2, nv)), sep = "")

  # Set Starting Values
  svBsex <- 1 # start value for regressions
  svB_birth_year_first <- 0.1 # start value for regressions
  svB_birth_year_second <- 0.1

  svTh <- .8 # start value for thresholds
  svCor <- .5 # start value for correlations
  lbCor <- -0.99 # lower bounds for correlations
  ubCor <- 0.99 # upper bounds for correlations

  # Create Matrices for Covariates and linear Regression Coefficients
  def_female_MZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Female1"), name = "def_female_MZ1")
  def_female_MZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Female2"), name = "def_female_MZ2")
  def_female_DZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Female1"), name = "def_female_DZ1")
  def_female_DZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Female2"), name = "def_female_DZ2")

  def_birth_year_first_MZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_first1"), name = "def_birth_year_first_MZ1")
  def_birth_year_first_MZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_first2"), name = "def_birth_year_first_MZ2")
  def_birth_year_first_DZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_first1"), name = "def_birth_year_first_DZ1")
  def_birth_year_first_DZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_first2"), name = "def_birth_year_first_DZ2")

  def_birth_year_second_MZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_second1"), name = "def_birth_year_second_MZ1")
  def_birth_year_second_MZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_second2"), name = "def_birth_year_second_MZ2")
  def_birth_year_second_DZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_second1"), name = "def_birth_year_second_DZ1")
  def_birth_year_second_DZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_second2"), name = "def_birth_year_second_DZ2")

  path_b_sex <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svBsex, label = "beta_fem", name = "bfem" )
  path_b_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_birth_year_first", name = "b_birth_year_first")
  path_b_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_birth_year_second", name = "b_birth_year_second")

  # Create Algebra for expected Mean & Threshold Matrices
  meanG <- mxMatrix(type = "Zero", nrow = 1, ncol = ntv, name = "meanG" )
  expMeanMZ <- mxAlgebra(expression = meanG +
                           cbind(def_female_MZ1 %*% bfem + def_birth_year_first_MZ1 %*% b_birth_year_first + def_birth_year_second_MZ1 %*% b_birth_year_second,
                                 def_female_MZ2 %*% bfem + def_birth_year_first_MZ2 %*% b_birth_year_first + def_birth_year_second_MZ2 %*% b_birth_year_second), name = "expMeanMZ")

  expMeanDZ <- mxAlgebra(expression = meanG + cbind(def_female_DZ1 %*% bfem + def_birth_year_first_DZ1 %*% b_birth_year_first + def_birth_year_second_DZ1 %*% b_birth_year_second,
                                                    def_female_DZ2 %*% bfem + def_birth_year_first_DZ2 %*% b_birth_year_first + def_birth_year_second_DZ2 %*% b_birth_year_second), name = "expMeanDZ")

  threMZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svTh, labels=c("tMZ1","tMZ2"), name="threMZ" )
  threDZ <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svTh, labels=c("tDZ1","tDZ2"), name="threDZ" )

  # Create Algebra for expected Correlation Matrices
  corMZ <-
    mxMatrix(
      type = "Stand",
      nrow = ntv,
      ncol = ntv,
      free = TRUE,
      values = svCor,
      lbound = lbCor,
      ubound = ubCor,
      labels = "rMZ",
      name = "corMZ"
    )

  corDZ <-
    mxMatrix(
      type = "Stand",
      nrow = ntv,
      ncol = ntv,
      free = TRUE,
      values = svCor,
      lbound = lbCor,
      ubound = ubCor,
      labels = "rDZ",
      name = "corDZ"
    )

  # Create Data Objects for Multiple Groups
  dataMZ <- mxData(observed = MZ, type = "raw")
  dataDZ <- mxData(observed = DZ, type = "raw")

  # Create Expectation Objects for Multiple Groups
  expMZ <-
    mxExpectationNormal(
      covariance = "corMZ",
      means = "expMeanMZ",
      dimnames = selVars,
      thresholds = "threMZ"
    )
  expDZ <-
    mxExpectationNormal(
      covariance = "corDZ",
      means = "expMeanDZ",
      dimnames = selVars,
      thresholds = "threDZ"
    )

  funML <- mxFitFunctionML()

  # Create Model Objects for Multiple Groups
  pars <- list(meanG, path_b_sex, path_b_birthyear_first, path_b_birthyear_second)
  defs <- list(def_female_MZ1, def_female_MZ2, def_female_DZ1, def_female_DZ2,
               def_birth_year_first_MZ1, def_birth_year_first_MZ2, def_birth_year_first_DZ1, def_birth_year_first_DZ2,
               def_birth_year_second_MZ1, def_birth_year_second_MZ2, def_birth_year_second_DZ1, def_birth_year_second_DZ2)

  modelMZ <-
    mxModel("MZ",
            pars,
            defs,
            expMeanMZ,
            corMZ,
            threMZ,
            dataMZ,
            expMZ,
            funML,
            name = "MZ")

  modelDZ <-
    mxModel("DZ",
            pars,
            defs,
            expMeanDZ,
            corDZ,
            threDZ,
            dataDZ,
            expDZ,
            funML,
            name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))

  # Create Confidence Interval Objects
  ciCor <- mxCI(c('MZ.corMZ', 'DZ.corDZ'))
  ciThre <- mxCI(c('MZ.threMZ', 'DZ.threDZ'))

  # Build Saturated Model with Confidence Intervals
  model_saturated <- mxModel("saturated", pars, modelMZ, modelDZ, multi, ciCor, ciThre)
  fit_saturated <- mxTryHard(model_saturated, intervals = TRUE)

  # RUN SUBMODELS
  model_no_cov  <- mxModel(fit_saturated, name = "no_cov" )
  model_no_cov  <- omxSetParameters(model_no_cov, label = c("beta_fem", "beta_birth_year_first", "beta_birth_year_second"), free = FALSE, values = 0)
  fit_no_cov    <- mxTryHard(model_no_cov)


  # Constrain expected Thresholds to be equal across Twin Order

  init_MZ <- mean(mxEval(c(tMZ1, tMZ2), fit_saturated))
  init_DZ <- mean(mxEval(c(tDZ1, tDZ2), fit_saturated))

  modelETO <- mxModel(fit_saturated, name = "ETO") %>%
    omxSetParameters(
      label = c("tMZ1", "tMZ2"),
      free = TRUE,
      values = init_MZ,
      newlabels = 'tMZ'
    ) %>% omxSetParameters(
      label = c("tDZ1", "tDZ2"),
      free = TRUE,
      values = init_DZ,
      newlabels = 'tDZ'
    )

  fitETO <- mxTryHard(modelETO)

  # Constrain expected Thresholds to be equal across Twin Order and Zygosity


  init_thresh <- mean(mxEval(c(tMZ, tDZ), fitETO))
  modelETOZ <- mxModel(fitETO, name = "ETOZ") %>%
    omxSetParameters(
      label = c("tMZ", "tDZ"),
      free = TRUE,
      values = init_thresh,
      newlabels = 'tZ'
    )

  fitETOZ <- mxTryHard(modelETOZ)

  out <- list(sat = fit_saturated,
              no_cov = fit_no_cov,
              ETO = fitETO,
              ETOZ = fitETOZ,
              trait = x$trait,
              response_type = x$response_type)

  class(out) <- "saturated.binary"

  out


}


#' @export
fit_saturated.prep.uni.num <- function(x, ...) {

  MZ <- x$MZ
  DZ <- x$DZ

  nv <- 1 # number of variables
  ntv <- nv * 2 # number of total variables
  selVars <- paste("X", c(rep(1, nv), rep(2, nv)), sep = "")


  svBsex <- 1 # start value for regressions
  svB_birth_year_first <- 0.1 # start value for regressions
  svB_birth_year_second <- 0.1
  svMe <- 1 # start value for means
  svVa <- .8 # start value for variance
  lbVa <- .0001

  def_female_MZ1      <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Female1"), name = "def_female_MZ1")
  def_female_MZ2      <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Female2"), name = "def_female_MZ2")
  def_female_DZ1      <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Female1"), name = "def_female_DZ1")
  def_female_DZ2      <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Female2"), name = "def_female_DZ2")

  def_birth_year_first_MZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_first1"), name = "def_birth_year_first_MZ1")
  def_birth_year_first_MZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_first2"), name = "def_birth_year_first_MZ2")
  def_birth_year_first_DZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_first1"), name = "def_birth_year_first_DZ1")
  def_birth_year_first_DZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_first2"), name = "def_birth_year_first_DZ2")

  def_birth_year_second_MZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_second1"), name = "def_birth_year_second_MZ1")
  def_birth_year_second_MZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("MZ.data.Birth_year_second2"), name = "def_birth_year_second_MZ2")
  def_birth_year_second_DZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_second1"), name = "def_birth_year_second_DZ1")
  def_birth_year_second_DZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("DZ.data.Birth_year_second2"), name = "def_birth_year_second_DZ2")

  path_b_sex <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svBsex, label = "beta_fem", name = "bfem" )
  path_b_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_birth_year_first", name = "b_birth_year_first")
  path_b_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_birth_year_second", name = "b_birth_year_second")

  # Create Algebra for expected Mean Matrices
  meanMZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mMZ1", name = "meanMZ1")
  meanMZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mMZ2", name = "meanMZ2")
  meanDZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mDZ1", name = "meanDZ1")
  meanDZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mDZ2", name = "meanDZ2")

  # Create Algebra for expected Mean & Threshold Matrices
  expMeanMZ   <- mxAlgebra(expression = cbind(meanMZ1 + def_female_MZ1 %*% bfem + def_birth_year_first_MZ1 %*% b_birth_year_first + def_birth_year_second_MZ1 %*% b_birth_year_second,
                                              meanMZ2 + def_female_MZ2 %*% bfem + def_birth_year_first_MZ2 %*% b_birth_year_first + def_birth_year_second_MZ2 %*% b_birth_year_second), name = "expMeanMZ")
  expMeanDZ   <- mxAlgebra(expression = cbind(meanDZ1 + def_female_DZ1 %*% bfem + def_birth_year_first_DZ1 %*% b_birth_year_first + def_birth_year_second_DZ1 %*% b_birth_year_second,
                                              meanDZ2 + def_female_DZ2 %*% bfem + def_birth_year_first_DZ2 %*% b_birth_year_first + def_birth_year_second_DZ2 %*% b_birth_year_second), name = "expMeanDZ")

  # Create Algebra for expected Variance/Covariance Matrices
  covMZ <- mxMatrix(
    type = "Symm",
    nrow = ntv,
    ncol = ntv,
    free = TRUE,
    values = diag(svVa, ntv, ntv),
    lbound = diag(lbVa, ntv, ntv),
    labels = c("vMZ1", "cMZ21", "vMZ2"),
    name = "covMZ"
  )

  covDZ <- mxMatrix(
    type = "Symm",
    nrow = ntv,
    ncol = ntv,
    free = TRUE,
    values = diag(svVa, ntv, ntv),
    lbound = diag(lbVa, ntv, ntv),
    labels = c("vDZ1", "cDZ21", "vDZ2"),
    name = "covDZ"
  )

  # Create Data Objects for Multiple Groups
  dataMZ <- mxData(observed = MZ, type = "raw")
  dataDZ <- mxData(observed = DZ, type = "raw")

  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal(covariance = "covMZ",
                               means = "expMeanMZ",
                               dimnames = selVars)

  expDZ <- mxExpectationNormal(covariance = "covDZ",
                               means = "expMeanDZ",
                               dimnames = selVars)

  funML <- mxFitFunctionML()

  pars <- list(path_b_sex, path_b_birthyear_first, path_b_birthyear_second, meanMZ1, meanMZ2, meanDZ1, meanDZ2)
  defs <- list(def_female_MZ1, def_female_MZ2, def_female_DZ1, def_female_DZ2,
               def_birth_year_first_MZ1, def_birth_year_first_MZ2, def_birth_year_first_DZ1, def_birth_year_first_DZ2,
               def_birth_year_second_MZ1, def_birth_year_second_MZ2, def_birth_year_second_DZ1, def_birth_year_second_DZ2)

  # Create Model Objects for Multiple Groups
  modelMZ <- mxModel(pars, defs, covMZ, dataMZ, expMeanMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, defs, covDZ, dataDZ, expMeanDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))

  # Build Saturated Model with Confidence Intervals
  model_saturated <- mxModel("saturated", modelMZ, modelDZ, multi)
  fit_saturated <- mxTryHard(model_saturated)

  # model without covariates
  model_no_cov  <- mxModel(fit_saturated, name = "no_cov" )
  model_no_cov  <- omxSetParameters(model_no_cov, label = c("beta_fem", "beta_birth_year_first", "beta_birth_year_second"), free = FALSE, values = 0)
  fit_no_cov    <- mxTryHard(model_no_cov)


  init_MZ <- mean(mxEval(c(mMZ1, mMZ2), fit_saturated))
  init_DZ <- mean(mxEval(c(mDZ1, mDZ2), fit_saturated))


  modelEMO <- mxModel(fit_saturated, name = "EMO") %>%
    omxSetParameters(label = c("mMZ1", "mMZ2"), free = TRUE, values = svMe, newlabels = 'mMZ') %>%
    omxSetParameters(label = c("mDZ1", "mDZ2"), free = TRUE, values = svMe, newlabels = 'mDZ')

  fitEMO <- mxTryHard(modelEMO)


  init_MZ <- mean(mxEval(c(vMZ1, vMZ2), fitEMO))
  init_DZ <- mean(mxEval(c(vDZ1, vDZ2), fitEMO))

  modelEMVO <- mxModel(fitEMO, name = "EMVO") %>%
    omxSetParameters(label = c("vMZ1", "vMZ2"), free = TRUE, values = init_MZ, newlabels = 'vMZ') %>%
    omxSetParameters(label = c("vDZ1", "vDZ2"), free = TRUE, values = init_DZ, newlabels = 'vDZ')

  fitEMVO <- mxTryHard(modelEMVO)


  # Constrain expected Means and Variances to be equal across Twin Order and Zygosity

  init_m <- mean(mxEval(c(mMZ, mDZ), fitEMVO))
  init_V <- mean(mxEval(c(vMZ, vDZ), fitEMVO))


  modelEMVOZ <- mxModel(fitEMVO, name = "EMVOZ") %>%
    omxSetParameters(label = c("mMZ", "mDZ"), free = TRUE, values = init_m, newlabels = 'mZ') %>%
    omxSetParameters(label = c("vMZ", "vDZ"), free = TRUE, values = init_V, newlabels = 'vZ')

  fitEMVOZ <- mxTryHard(modelEMVOZ)

  out <- list(
    sat = fit_saturated,
    no_cov = fit_no_cov,
    EMO = fitEMO,
    EMVO = fitEMVO,
    EMVOZ = fitEMVOZ,
    trait = x$trait,
    response_type = x$response_type
  )

  class(out) <- "saturated.num"

  out

}


#' @export
fit_saturated.prep.uni.5group.binary <- function(x, ...) {

  # Select Variables for Analysis
  vars      <- "X"
  nv        <- 1
  ntv       <- nv * 2
  selVars   <- paste(vars, c(rep(1, nv), rep(2, nv)), sep="")

  svB_birth_year_first <- 0.05 # start value for regressions
  svB_birth_year_second <- 0.05
  svTh      <- 0.8                       # start value for thresholds
  svCor     <- 0.5                       # start value for correlations
  lbCor     <- -0.99                     # lower bounds for correlations
  ubCor     <- 0.99                      # upper bounds for correlations



  # Definition variables

  def_birth_year_first <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("data.Birth_year_first1"), name = "def_birth_year_first")
  def_birth_year_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("data.Birth_year_second1"), name = "def_birth_year_second")


  # Regression parameters

  path_bm_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_m_birth_year_first", name = "bm_birth_year_first")
  path_bm_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_m_birth_year_second", name = "bm_birth_year_second")

  path_bf_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_f_birth_year_first", name = "bf_birth_year_first")
  path_bf_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_f_birth_year_second", name = "bf_birth_year_second")


  # Create Algebra for expected Mean & Threshold Matrices
  meanG <- mxMatrix(type = "Zero", nrow = 1, ncol = ntv, name = "meanG" )

  expMean_zf <- mxAlgebra(expression = meanG +
                            cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                  def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second),
                          name = "expMean_zf")


  expMean_zm <- mxAlgebra(expression = meanG +
                            cbind(def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second,
                                  def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                          name = "expMean_zm")

  expMean_zo <- mxAlgebra(expression = meanG +
                            cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                  def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                          name = "expMean_zo")

  funML     <- mxFitFunctionML()

  # Create algebra for thresholds
  thre_mzf   <- mxMatrix(type="Full", nrow = 1, ncol = 2, free = TRUE, values = svTh, labels = c("tmzf1", "tmzf2"), name = "thre_mzf")
  thre_mzm   <- mxMatrix(type="Full", nrow = 1, ncol = 2, free = TRUE, values = svTh, labels = c("tmzm1", "tmzm2"), name = "thre_mzm")
  thre_dzf   <- mxMatrix(type="Full", nrow = 1, ncol = 2, free = TRUE, values = svTh, labels = c("tdzf1", "tdzf2"), name = "thre_dzf")
  thre_dzm   <- mxMatrix(type="Full", nrow = 1, ncol = 2, free = TRUE, values = svTh, labels = c("tdzm1", "tdzm2"), name = "thre_dzm")
  thre_dzo   <- mxMatrix(type="Full", nrow = 1, ncol = 2, free = TRUE, values = svTh, labels = c("tdzo1", "tdzo2"), name = "thre_dzo")


  # Create Algebra for expected Correlation Matrices
  cor_mzf    <- mxMatrix(type = "Stand", nrow = 2, ncol = 2, free = TRUE, values = svCor, lbound = lbCor, ubound = ubCor, labels = "r_mzf_21", name = "cor_mzf")
  cor_mzm    <- mxMatrix(type = "Stand", nrow = 2, ncol = 2, free = TRUE, values = svCor, lbound = lbCor, ubound = ubCor, labels = "r_mzm_21", name = "cor_mzm")
  cor_dzf    <- mxMatrix(type = "Stand", nrow = 2, ncol = 2, free = TRUE, values = svCor, lbound = lbCor, ubound = ubCor, labels = "r_dzf_21", name = "cor_dzf")
  cor_dzm    <- mxMatrix(type = "Stand", nrow = 2, ncol = 2, free = TRUE, values = svCor, lbound = lbCor, ubound = ubCor, labels = "r_dzm_21", name = "cor_dzm")
  cor_dzo    <- mxMatrix(type = "Stand", nrow = 2, ncol = 2, free = TRUE, values = svCor, lbound = lbCor, ubound = ubCor, labels = "r_dzo_21", name = "cor_dzo")

  # Create data objects
  data_mzf <- mxData(observed = x$mzf, type = "raw")
  data_mzm <- mxData(observed = x$mzm, type = "raw")
  data_dzf <- mxData(observed = x$dzf, type = "raw")
  data_dzm <- mxData(observed = x$dzm, type = "raw")
  data_dzo <- mxData(observed = x$dzo, type = "raw")

  # Create expectation objects
  exp_mzf    <- mxExpectationNormal(covariance = "cor_mzf", means = "expMean_zf", dimnames = selVars, thresholds = "thre_mzf" )
  exp_mzm    <- mxExpectationNormal(covariance = "cor_mzm", means = "expMean_zm", dimnames = selVars, thresholds = "thre_mzm" )
  exp_dzf    <- mxExpectationNormal(covariance = "cor_dzf", means = "expMean_zf", dimnames = selVars, thresholds = "thre_dzf" )
  exp_dzm    <- mxExpectationNormal(covariance = "cor_dzm", means = "expMean_zm", dimnames = selVars, thresholds = "thre_dzm" )
  exp_dzo    <- mxExpectationNormal(covariance = "cor_dzo", means = "expMean_zo", dimnames = selVars, thresholds = "thre_dzo" )

  # Create model objects for multiple groups
  pars      <- list(path_bm_birthyear_first, path_bm_birthyear_second, path_bf_birthyear_first, path_bf_birthyear_second, meanG)
  defs      <- list(def_birth_year_first, def_birth_year_second)

  model_mzf  <- mxModel(pars, defs, expMean_zf, cor_mzf, thre_mzf, data_mzf, exp_mzf, funML, name="mzf")
  model_mzm  <- mxModel(pars, defs, expMean_zm, cor_mzm, thre_mzm, data_mzm, exp_mzm, funML, name="mzm")
  model_dzf  <- mxModel(pars, defs, expMean_zf, cor_dzf, thre_dzf, data_dzf, exp_dzf, funML, name="dzf")
  model_dzm  <- mxModel(pars, defs, expMean_zm, cor_dzm, thre_dzm, data_dzm, exp_dzm, funML, name="dzm")
  model_dzo  <- mxModel(pars, defs, expMean_zo, cor_dzo, thre_dzo, data_dzo, exp_dzo, funML, name="dzo")

  multi     <- mxFitFunctionMultigroup(c("mzf","mzm", "dzf","dzm","dzo"))

  model_sat  <- mxModel( "sat_5group_binary",
                         pars,
                         model_mzf,
                         model_mzm,
                         model_dzf,
                         model_dzm,
                         model_dzo,
                         multi)

  fit_saturated <- mxTryHard(model_sat)

  # RUN SUBMODELS

  # Test Significance of Covariate
  model_no_cov  <- mxModel(fit_saturated, name = "no_cov")
  model_no_cov  <- omxSetParameters(model_no_cov, label = c("beta_m_birth_year_first", "beta_m_birth_year_second", "beta_f_birth_year_first", "beta_f_birth_year_second"), free = FALSE, values = 0)
  fit_no_cov    <- mxTryHard(model_no_cov)

  # Test Sex Difference in Covariate
  model_no_sex_diff_cov <- mxModel(fit_saturated, name="no_sex_diff_cov" )
  model_no_sex_diff_cov <- omxSetParameters(model_no_sex_diff_cov, label=c("beta_m_birth_year_first","beta_f_birth_year_first"), free=TRUE, values = svB_birth_year_first, newlabels = "beta_birth_year_first")
  model_no_sex_diff_cov <- omxSetParameters(model_no_sex_diff_cov, label=c("beta_m_birth_year_second","beta_f_birth_year_second"), free=TRUE, values = svB_birth_year_second, newlabels = "beta_birth_year_second")
  fit_no_sex_diff_cov   <- mxTryHard(model_no_sex_diff_cov)


  # Constrain expected Thresholds to be equal across twin order
  modelETO  <- mxModel(fit_saturated, name = "ETO")
  modelETO  <- omxSetParameters(modelETO, label=c("tmzf1","tmzf2"), free=TRUE, values=svTh, newlabels='tmzf')
  modelETO  <- omxSetParameters(modelETO, label=c("tdzf1","tdzf2"), free=TRUE, values=svTh, newlabels='tdzf')
  modelETO  <- omxSetParameters(modelETO, label=c("tmzm1","tmzm2"), free=TRUE, values=svTh, newlabels='tmzm')
  modelETO  <- omxSetParameters(modelETO, label=c("tdzm1","tdzm2"), free=TRUE, values=svTh, newlabels='tdzm')
  fitETO    <- mxTryHard(modelETO)

  # Constrain expected Thresholds to be equal across twin order and zygosity
  modelETOZ  <- mxModel(fitETO, name = "ETOZ" )
  modelETOZ  <- omxSetParameters(modelETOZ, label = c("tmzf", "tdzf"), free = TRUE, values = svTh, newlabels = "tzf")
  modelETOZ  <- omxSetParameters(modelETOZ, label = c("tmzm", "tdzm"), free = TRUE, values = svTh, newlabels = "tzm")
  fitETOZ    <- mxTryHard(modelETOZ)

  # Constrain expected Thresholds to be equal across twin order and zygosity and SS/OS
  modelETP  <- mxModel(fitETOZ, name="ETP")
  modelETP  <- omxSetParameters(modelETP, label=c("tzf","tdzo1"), free=TRUE, values=svTh, newlabels="tzf")
  modelETP  <- omxSetParameters(modelETP, label=c("tzm","tdzo2"), free=TRUE, values=svTh, newlabels="tzm")
  fitETP    <- mxTryHard(modelETP)

  # Constrain expected Thresholds to be equal across twin order and zygosity and SS/OS and sex
  modelETS  <- mxModel(fitETP, name="ETS" )
  modelETS  <- omxSetParameters(modelETS, label=c("tzf","tzm"), free = TRUE, values = svTh, newlabels = "tZ")
  fitETS    <- mxTryHard(modelETS)

  out <- list(sat = fit_saturated,
              no_cov = fit_no_cov,
              no_sex_diff_cov = fit_no_sex_diff_cov,
              ETO = fitETO,
              ETOZ = fitETOZ,
              ETP = fitETP,
              ETS =fitETS,
              response_type = x$response_type,
              trait = x$trait)

  class(out) <- "saturated.5group.binary"
  out


}



#' @export
fit_saturated.prep.uni.5group.num <- function(x, ...) {

  # Select Variables for Analysis
  vars      <- "X"
  nv        <- 1
  ntv       <- nv * 2
  selVars   <- paste(vars, c(rep(1, nv), rep(2, nv)), sep="")


  # Set Starting Values
  svB_birth_year_first <- 0.05 # start value for regressions
  svB_birth_year_second <- 0.05
  svMe <- 20
  svVa      <- 0.8                       # start value for variance
  lbVa      <- 0.0001                    # start value for lower bounds

  # Definition variables

  def_birth_year_first <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("data.Birth_year_first1"), name = "def_birth_year_first")
  def_birth_year_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = FALSE, labels = c("data.Birth_year_second1"), name = "def_birth_year_second")


  # Regression parameters

  path_bm_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_m_birth_year_first", name = "bm_birth_year_first")
  path_bm_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_m_birth_year_second", name = "bm_birth_year_second")

  path_bf_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_f_birth_year_first", name = "bf_birth_year_first")
  path_bf_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_f_birth_year_second", name = "bf_birth_year_second")


  # Create Algebra for expected Mean Matrices
  mean_mzf   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZf1","mMZf2"), name="mean_mzf")
  mean_dzf   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZf1","mDZf2"), name="mean_dzf")
  mean_mzm   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZm1","mMZm2"), name="mean_mzm")
  mean_dzm   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZm1","mDZm2"), name="mean_dzm")
  mean_dzo   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZo1","mDZo2"), name="mean_dzo")

  expMean_mzf <- mxAlgebra(expression = mean_mzf +
                             cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                   def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second),
                           name = "expMean_mzf")

  expMean_dzf <- mxAlgebra(expression = mean_dzf +
                             cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                   def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second),
                           name = "expMean_dzf")

  expMean_mzm <- mxAlgebra(expression = mean_mzm +
                             cbind(def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second,
                                   def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                           name = "expMean_mzm")

  expMean_dzm <- mxAlgebra(expression = mean_dzm +
                             cbind(def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second,
                                   def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                           name = "expMean_dzm")

  expMean_dzo <- mxAlgebra(expression = mean_dzo +
                             cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                   def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                           name = "expMean_dzo")

  funML     <- mxFitFunctionML()

  # Create Algebra for expected Variance/Covariance Matrices
  cov_mzf <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values = diag(svVa, nrow = 2, ncol = 2), lbound = diag(lbVa, nrow = 2, ncol = 2), labels=c("vMZf1","cMZf21","vMZf2"), name="cov_mzf")
  cov_dzf <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values = diag(svVa, nrow = 2, ncol = 2), lbound = diag(lbVa, nrow = 2, ncol = 2), labels=c("vDZf1","cDZf21","vDZf2"), name="cov_dzf")
  cov_mzm <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values = diag(svVa, nrow = 2, ncol = 2), lbound = diag(lbVa, nrow = 2, ncol = 2), labels=c("vMZm1","cMZm21","vMZm2"), name="cov_mzm")
  cov_dzm <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values = diag(svVa, nrow = 2, ncol = 2), lbound = diag(lbVa, nrow = 2, ncol = 2), labels=c("vDZm1","cDZm21","vDZm2"), name="cov_dzm")
  cov_dzo <- mxMatrix(type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values = diag(svVa, nrow = 2, ncol = 2), lbound = diag(lbVa, nrow = 2, ncol = 2), labels=c("vDZo1","cDZo21","vDZo2"), name="cov_dzo")

  # Create Algebra for Maximum Likelihood Estimates of Twin Correlations
  matI      <- mxMatrix(type = "Iden", nrow = ntv, ncol = ntv, name = "I")
  cor_mzf    <- mxAlgebra(solve(sqrt(I * cov_mzf)) %&% cov_mzf, name = "cor_mzf")
  cor_dzf    <- mxAlgebra(solve(sqrt(I * cov_dzf)) %&% cov_dzf, name = "cor_dzf")
  cor_mzm    <- mxAlgebra(solve(sqrt(I * cov_mzm)) %&% cov_mzm, name = "cor_mzm")
  cor_dzm    <- mxAlgebra(solve(sqrt(I * cov_dzm)) %&% cov_dzm, name = "cor_dzm")
  cor_dzo    <- mxAlgebra(solve(sqrt(I * cov_dzo)) %&% cov_dzo, name = "cor_dzo")

  # Create data objects
  data_mzf <- mxData(observed = x$mzf, type = "raw")
  data_dzf <- mxData(observed = x$dzf, type = "raw")
  data_mzm <- mxData(observed = x$mzm, type = "raw")
  data_dzm <- mxData(observed = x$dzm, type = "raw")
  data_dzo <- mxData(observed = x$dzo, type = "raw")

  # Create expectation objects
  exp_mzf    <- mxExpectationNormal(covariance = "cov_mzf", means = "expMean_mzf", dimnames = selVars)
  exp_dzf   <- mxExpectationNormal(covariance = "cov_dzf", means = "expMean_dzf", dimnames = selVars)
  exp_mzm   <- mxExpectationNormal(covariance = "cov_mzm", means = "expMean_mzm", dimnames = selVars)
  exp_dzm    <- mxExpectationNormal(covariance = "cov_dzm", means = "expMean_dzm", dimnames = selVars)
  exp_dzo    <- mxExpectationNormal(covariance = "cov_dzo", means = "expMean_dzo", dimnames = selVars)
  funML     <- mxFitFunctionML()

  # Create model objects for multiple groups
  pars      <- list(path_bm_birthyear_first, path_bm_birthyear_second, path_bf_birthyear_first, path_bf_birthyear_second)
  defs      <- list(def_birth_year_first, def_birth_year_second)

  model_mzf  <- mxModel(pars, defs, mean_mzf, expMean_mzf, cov_mzf, cor_mzf, matI, data_mzf, exp_mzf, funML, name="mzf")
  model_dzf  <- mxModel(pars, defs, mean_dzf, expMean_dzf, cov_dzf, cor_dzf, matI, data_dzf, exp_dzf, funML, name="dzf")
  model_mzm  <- mxModel(pars, defs, mean_mzm, expMean_mzm, cov_mzm, cor_mzm, matI, data_mzm, exp_mzm, funML, name="mzm")
  model_dzm  <- mxModel(pars, defs, mean_dzm, expMean_dzm, cov_dzm, cor_dzm, matI, data_dzm, exp_dzm, funML, name="dzm")
  model_dzo  <- mxModel(pars, defs, mean_dzo, expMean_dzo, cov_dzo, cor_dzo, matI, data_dzo, exp_dzo, funML, name="dzo")

  multi     <- mxFitFunctionMultigroup(c("mzf", "dzf", "mzm", "dzm","dzo"))

  model_sat  <- mxModel( "sat_5group_num",
                         pars,
                         model_mzf,
                         model_dzf,
                         model_mzm,
                         model_dzm,
                         model_dzo,
                         multi)

  fit_sat <- mxTryHard(model_sat)

  # Test covariates
  model_no_cov  <- mxModel(fit_sat, name = "no_cov")
  model_no_cov  <- omxSetParameters(model_no_cov, label = c("beta_m_birth_year_first", "beta_m_birth_year_second", "beta_f_birth_year_first", "beta_f_birth_year_second"), free = FALSE, values = 0)
  fit_no_cov    <- mxTryHard(model_no_cov)

  # Test Sex Difference in Covariate
  model_no_sex_diff_cov <- mxModel(fit_sat, name="no_sex_diff_cov" )
  model_no_sex_diff_cov <- omxSetParameters(model_no_sex_diff_cov, label=c("beta_m_birth_year_first","beta_f_birth_year_first"), free = TRUE, values = svB_birth_year_first, newlabels = "beta_birth_year_first")
  model_no_sex_diff_cov <- omxSetParameters(model_no_sex_diff_cov, label=c("beta_m_birth_year_second","beta_f_birth_year_second"), free = TRUE, values = svB_birth_year_second, newlabels = "beta_birth_year_second")
  fit_no_sex_diff_cov   <- mxTryHard(model_no_sex_diff_cov)


  # Constrain expected Means to be equal across twin order
  modelEMO  <- mxModel(fit_sat, name="EMO" )
  modelEMO  <- omxSetParameters(modelEMO, label=c("mMZf1","mMZf2"), free=TRUE, values=svMe, newlabels='mMZf' )
  modelEMO  <- omxSetParameters(modelEMO, label=c("mDZf1","mDZf2"), free=TRUE, values=svMe, newlabels='mDZf' )
  modelEMO  <- omxSetParameters(modelEMO, label=c("mMZm1","mMZm2"), free=TRUE, values=svMe, newlabels='mMZm' )
  modelEMO  <- omxSetParameters(modelEMO, label=c("mDZm1","mDZm2"), free=TRUE, values=svMe, newlabels='mDZm' )
  fitEMO    <- mxTryHard(modelEMO)

  # Constrain expected Means and Variances to be equal across twin order
  modelEMVO <- mxModel(fitEMO, name="EMVO" )
  modelEMVO <- omxSetParameters(modelEMVO, label=c("vMZf1","vMZf2"), free=TRUE, values=svVa, newlabels='vMZf')
  modelEMVO <- omxSetParameters(modelEMVO, label=c("vDZf1","vDZf2"), free=TRUE, values=svVa, newlabels='vDZf')
  modelEMVO <- omxSetParameters(modelEMVO, label=c("vMZm1","vMZm2"), free=TRUE, values=svVa, newlabels='vMZm')
  modelEMVO <- omxSetParameters(modelEMVO, label=c("vDZm1","vDZm2"), free=TRUE, values=svVa, newlabels='vDZm')
  fitEMVO   <- mxTryHard(modelEMVO)

  # Constrain expected Means and Variances to be equal across twin order and zygosity
  modelEMVZ <- mxModel( fitEMVO, name="EMVZ" )
  modelEMVZ <- omxSetParameters(modelEMVZ, label=c("mMZf","mDZf"), free=TRUE, values=svMe, newlabels='mZf' )
  modelEMVZ <- omxSetParameters(modelEMVZ, label=c("vMZf","vDZf"), free=TRUE, values=svVa, newlabels='vZf' )
  modelEMVZ <- omxSetParameters(modelEMVZ, label=c("mMZm","mDZm"), free=TRUE, values=svMe, newlabels='mZm' )
  modelEMVZ <- omxSetParameters(modelEMVZ, label=c("vMZm","vDZm"), free=TRUE, values=svVa, newlabels='vZm' )
  fitEMVZ   <- mxTryHard(modelEMVZ)

  # Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS
  modelEMVP <- mxModel(fitEMVZ, name="EMVP" )
  modelEMVP <- omxSetParameters(modelEMVP, label=c("mZf","mDZo1"), free=TRUE, values=svMe, newlabels='mZf')
  modelEMVP <- omxSetParameters(modelEMVP, label=c("vZf","vDZo1"), free=TRUE, values=svVa, newlabels='vZf')
  modelEMVP <- omxSetParameters(modelEMVP, label=c("mZm","mDZo2"), free=TRUE, values=svMe, newlabels='mZm')
  modelEMVP <- omxSetParameters(modelEMVP, label=c("vZm","vDZo2"), free=TRUE, values=svVa, newlabels='vZm')
  fitEMVP   <- mxTryHard(modelEMVP)

  # Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS and sex
  modelEMVS <- mxModel(fitEMVP, name = "EMVS")
  modelEMVS <- omxSetParameters(modelEMVS, label=c("mZf","mZm"), free=TRUE, values=svMe, newlabels='mZ')
  modelEMVS <- omxSetParameters(modelEMVS, label=c("vZf","vZm"), free=TRUE, values=svVa, newlabels='vZ')
  fitEMVS   <- mxRun(modelEMVS)


  out <- list(sat = fit_sat,
              no_cov = fit_no_cov,
              no_sex_diff_cov = fit_no_sex_diff_cov,
              EMO = fitEMO,
              EMVO = fitEMVO,
              EMVZ = fitEMVZ,
              EMVP = fitEMVP,
              EMVS = fitEMVS,
              response_type = x$response_type,
              trait = x$trait)

  class(out) <- "saturated.5group.num"
  out


}



