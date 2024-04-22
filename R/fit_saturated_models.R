
#' @export
fit_saturated <- function(x, ...) {

  UseMethod("fit_saturated")

}



#' @export
fit_saturated.prep.uni.bin <- function(x, covs, extra_tries = 10, ...) {

  nv <- 1 # number of variables
  ntv <- nv * 2 # number of total variables
  selVars <- paste("X", c(rep(1, nv), rep(2, nv)), sep = "")

  # Set Starting Values
  svB <- 1

  svTh <- .8 # start value for thresholds
  svCor <- .5 # start value for correlations
  lbCor <- -0.99 # lower bounds for correlations
  ubCor <- 0.99 # upper bounds for correlations

  # Create Matrices for Covariates and linear Regression Coefficients
  defs1 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "1"), name = "defs1")
  defs2 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "2"), name = "defs2")

  path_b  <- mxMatrix(type = "Full", nrow = length(covs), ncol = 1, free = TRUE, values = rep(0, length(covs)), label = paste0("beta_", covs), name = "b")

  # Create Algebra for expected Mean & Threshold Matrices
  meanG <- mxMatrix(type = "Zero", nrow = 1, ncol = ntv, name = "meanG" )
  expMeanMZ <- mxAlgebra(expression = meanG + cbind(defs1 %*% b, defs2 %*% b), name = "expMeanMZ")
  expMeanDZ <- mxAlgebra(expression = meanG + cbind(defs1 %*% b, defs2 %*% b), name = "expMeanDZ")

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
  dataMZ <- mxData(observed = as.data.frame(x$MZ), type = "raw")
  dataDZ <- mxData(observed = as.data.frame(x$DZ), type = "raw")

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
  pars <- list(meanG, path_b)
  defs <- list(defs1, defs2)

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

  # Build Saturated Model with Confidence Intervals
  model_saturated <- mxModel("sat", pars, modelMZ, modelDZ, multi)
  fit_saturated <- mxTryHardOrdinal(model_saturated, extraTries = extra_tries)

  # RUN SUBMODELS
  if (!is_null(covs)) {

    # Test covariates
    model_no_cov  <- mxModel(fit_saturated, name = "sat_no_cov")
    model_no_cov  <- omxSetParameters(model_no_cov, label = c(paste0("beta_", covs)), free = FALSE, values = 0)
    fit_no_cov    <- mxTryHardOrdinal(model_no_cov, extraTries = extra_tries)

  }


  # Constrain expected Thresholds to be equal across Twin Order

  init_MZ <- mean(mxEval(c(tMZ1, tMZ2), fit_saturated))
  init_DZ <- mean(mxEval(c(tDZ1, tDZ2), fit_saturated))

  modelETO <- mxModel(fit_saturated, name = "sat_equal_order") %>%
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

  fitETO <- mxTryHardOrdinal(modelETO, extraTries = extra_tries)

  # Constrain expected Thresholds to be equal across Twin Order and Zygosity


  init_thresh <- mean(mxEval(c(tMZ, tDZ), fitETO))
  modelETOZ <- mxModel(fitETO, name = "sat_equal_order_zyg") %>%
    omxSetParameters(
      label = c("tMZ", "tDZ"),
      free = TRUE,
      values = init_thresh,
      newlabels = 'tZ'
    )

  fitETOZ <- mxTryHardOrdinal(modelETOZ, extraTries = extra_tries)

  out <- list(sat = fit_saturated,
              ETO = fitETO,
              ETOZ = fitETOZ,
              trait = x$trait,
              response_type = x$response_type)

  if (!is_null(covs)) out$no_cov <- fit_no_cov

  class(out) <- "saturated.binary"

  out


}





#' @export
fit_saturated.prep.uni.num <- function(x, covs, extra_tries = 10,...) {

  MZ <- x$MZ
  DZ <- x$DZ

  nv <- 1 # number of variables
  ntv <- nv * 2 # number of total variables
  selVars <- paste("X", c(rep(1, nv), rep(2, nv)), sep = "")


  svB <- 1 # start value for regressions
  svMe <- 1 # start value for means
  svVa <- .8 # start value for variance
  lbVa <- .0001

  defs1 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "1"), name = "defs1")
  defs2 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "2"), name = "defs2")

  path_b  <- mxMatrix(type = "Full", nrow = length(covs), ncol = 1, free = TRUE, values = rep(0, length(covs)), label = paste0("beta_", covs), name = "b")

  # Create Algebra for expected Mean Matrices
  meanMZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mMZ1", name = "meanMZ1")
  meanMZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mMZ2", name = "meanMZ2")
  meanDZ1 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mDZ1", name = "meanDZ1")
  meanDZ2 <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svMe, labels = "mDZ2", name = "meanDZ2")

  # Create Algebra for expected Mean & Threshold Matrices
  expMeanMZ   <- mxAlgebra(expression = cbind(meanMZ1 + defs1 %*% b, meanMZ2 + defs2 %*% b), name = "expMeanMZ")
  expMeanDZ   <- mxAlgebra(expression = cbind(meanDZ1 + defs1 %*% b, meanDZ2 + defs2 %*% b), name = "expMeanDZ")

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
  dataMZ <- mxData(observed = as.data.frame(MZ), type = "raw")
  dataDZ <- mxData(observed = as.data.frame(DZ), type = "raw")

  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal(covariance = "covMZ",
                               means = "expMeanMZ",
                               dimnames = selVars)

  expDZ <- mxExpectationNormal(covariance = "covDZ",
                               means = "expMeanDZ",
                               dimnames = selVars)

  funML <- mxFitFunctionML()

  pars <- list(path_b, meanMZ1, meanMZ2, meanDZ1, meanDZ2)
  defs <- list(defs1, defs2)

  # Create Model Objects for Multiple Groups
  modelMZ <- mxModel(pars, defs, covMZ, dataMZ, expMeanMZ, expMZ, funML, name = "MZ")
  modelDZ <- mxModel(pars, defs, covDZ, dataDZ, expMeanDZ, expDZ, funML, name = "DZ")
  multi <- mxFitFunctionMultigroup(c("MZ", "DZ"))

  # Build Saturated Model with Confidence Intervals
  model_saturated <- mxModel("sat", modelMZ, modelDZ, multi)
  fit_saturated <- mxTryHard(model_saturated, exhaustive = TRUE, extraTries = extra_tries)

  # RUN SUBMODELS
  if (!is_null(covs)) {

    # Test covariates
    model_no_cov  <- mxModel(fit_saturated, name = "sat_no_cov")
    model_no_cov  <- omxSetParameters(model_no_cov, label = c(paste0("beta_", covs)), free = FALSE, values = 0)
    fit_no_cov    <- mxTryHard(model_no_cov, exhaustive = TRUE, extraTries = extra_tries)

  }


  init_MZ <- mean(mxEval(c(mMZ1, mMZ2), fit_saturated))
  init_DZ <- mean(mxEval(c(mDZ1, mDZ2), fit_saturated))


  modelEMO <- mxModel(fit_saturated, name = "sat_equal_mean_order") %>%
    omxSetParameters(label = c("mMZ1", "mMZ2"), free = TRUE, values = svMe, newlabels = 'mMZ') %>%
    omxSetParameters(label = c("mDZ1", "mDZ2"), free = TRUE, values = svMe, newlabels = 'mDZ')

  fitEMO <- mxTryHard(modelEMO, exhaustive = TRUE, extraTries = extra_tries)


  init_MZ <- mean(mxEval(c(vMZ1, vMZ2), fitEMO))
  init_DZ <- mean(mxEval(c(vDZ1, vDZ2), fitEMO))

  modelEMVO <- mxModel(fitEMO, name = "sat_equal_mean_var_order") %>%
    omxSetParameters(label = c("vMZ1", "vMZ2"), free = TRUE, values = init_MZ, newlabels = 'vMZ') %>%
    omxSetParameters(label = c("vDZ1", "vDZ2"), free = TRUE, values = init_DZ, newlabels = 'vDZ')

  fitEMVO <- mxTryHard(modelEMVO, exhaustive = TRUE, extraTries = extra_tries)


  # Constrain expected Means and Variances to be equal across Twin Order and Zygosity

  init_m <- mean(mxEval(c(mMZ, mDZ), fitEMVO))
  init_V <- mean(mxEval(c(vMZ, vDZ), fitEMVO))


  modelEMVOZ <- mxModel(fitEMVO, name = "sat_equal_mean_var_order_zyg") %>%
    omxSetParameters(label = c("mMZ", "mDZ"), free = TRUE, values = init_m, newlabels = 'mZ') %>%
    omxSetParameters(label = c("vMZ", "vDZ"), free = TRUE, values = init_V, newlabels = 'vZ')

  fitEMVOZ <- mxTryHard(modelEMVOZ, exhaustive = TRUE, extraTries = extra_tries)

  out <- list(
    sat = fit_saturated,
    EMO = fitEMO,
    EMVO = fitEMVO,
    EMVOZ = fitEMVOZ,
    trait = x$trait,
    response_type = x$response_type
  )

  if (!is_null(covs)) out$no_cov <- fit_no_cov

  class(out) <- "saturated.num"

  out

}


#' @export
fit_saturated.prep.uni.5group.binary <- function(x, covs, extra_tries = 10, ...) {

  # Select Variables for Analysis
  vars      <- "X"
  nv        <- 1
  ntv       <- nv * 2
  selVars   <- paste(vars, c(rep(1, nv), rep(2, nv)), sep="")

  svB <- 20
  svTh      <- 0.8                       # start value for thresholds
  svCor     <- 0.5                       # start value for correlations
  lbCor     <- -0.99                     # lower bounds for correlations
  ubCor     <- 0.99                      # upper bounds for correlations



  # Definition variables
  defs1 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "1"), name = "defs1")
  defs2 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "2"), name = "defs2")

  # Regression parameters
  path_bm  <- mxMatrix(type = "Full", nrow = length(covs), ncol = 1, free = TRUE, values = rep(svB, length(covs)), label = paste0("beta_m_", covs), name = "bm")
  path_bf  <- mxMatrix(type = "Full", nrow = length(covs), ncol = 1, free = TRUE, values = rep(svB, length(covs)), label = paste0("beta_f_", covs), name = "bf")


  # Create Algebra for expected Mean & Threshold Matrices
  meanG <- mxMatrix(type = "Zero", nrow = 1, ncol = ntv, name = "meanG" )

  expMean_zf <- mxAlgebra(expression = meanG + cbind(defs1 %*% bf, defs2 %*% bf),
                          name = "expMean_zf")


  expMean_zm <- mxAlgebra(expression = meanG + cbind(defs1 %*% bm, defs2 %*% bm),
                          name = "expMean_zm")

  expMean_zo <- mxAlgebra(expression = meanG + cbind(defs1 %*% bf, defs2 %*% bm),
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
  data_mzf <- mxData(observed = as.data.frame(x$mzf), type = "raw")
  data_mzm <- mxData(observed = as.data.frame(x$mzm), type = "raw")
  data_dzf <- mxData(observed = as.data.frame(x$dzf), type = "raw")
  data_dzm <- mxData(observed = as.data.frame(x$dzm), type = "raw")
  data_dzo <- mxData(observed = as.data.frame(x$dzo), type = "raw")

  # Create expectation objects
  exp_mzf    <- mxExpectationNormal(covariance = "cor_mzf", means = "expMean_zf", dimnames = selVars, thresholds = "thre_mzf" )
  exp_mzm    <- mxExpectationNormal(covariance = "cor_mzm", means = "expMean_zm", dimnames = selVars, thresholds = "thre_mzm" )
  exp_dzf    <- mxExpectationNormal(covariance = "cor_dzf", means = "expMean_zf", dimnames = selVars, thresholds = "thre_dzf" )
  exp_dzm    <- mxExpectationNormal(covariance = "cor_dzm", means = "expMean_zm", dimnames = selVars, thresholds = "thre_dzm" )
  exp_dzo    <- mxExpectationNormal(covariance = "cor_dzo", means = "expMean_zo", dimnames = selVars, thresholds = "thre_dzo" )

  # Create model objects for multiple groups
  pars      <- list(meanG, path_bm, path_bf)
  defs      <- list(defs1, defs2)

  model_mzf  <- mxModel(pars, defs, expMean_zf, cor_mzf, thre_mzf, data_mzf, exp_mzf, funML, name="mzf")
  model_mzm  <- mxModel(pars, defs, expMean_zm, cor_mzm, thre_mzm, data_mzm, exp_mzm, funML, name="mzm")
  model_dzf  <- mxModel(pars, defs, expMean_zf, cor_dzf, thre_dzf, data_dzf, exp_dzf, funML, name="dzf")
  model_dzm  <- mxModel(pars, defs, expMean_zm, cor_dzm, thre_dzm, data_dzm, exp_dzm, funML, name="dzm")
  model_dzo  <- mxModel(pars, defs, expMean_zo, cor_dzo, thre_dzo, data_dzo, exp_dzo, funML, name="dzo")

  multi     <- mxFitFunctionMultigroup(c("mzf","mzm", "dzf","dzm","dzo"))

  model_sat  <- mxModel("sat_5group_binary",
                        pars,
                        model_mzf,
                        model_mzm,
                        model_dzf,
                        model_dzm,
                        model_dzo,
                        multi)

  fit_saturated <- mxTryHardOrdinal(model_sat, extraTries = extra_tries)

  # RUN SUBMODELS
  if (!is_null(covs)) {

    # Test covariates
    model_no_cov  <- mxModel(fit_saturated, name = "sat_5group_no_cov")
    model_no_cov  <- omxSetParameters(model_no_cov, label = c(paste0("beta_m_", covs), paste0("beta_f_", covs)), free = FALSE, values = 0)
    fit_no_cov    <- mxTryHardOrdinal(model_no_cov, extraTries = extra_tries)

    # Test Sex Difference in Covariate
    model_no_sex_diff_cov <- mxModel(fit_saturated, name = "sat_5group_no_sex_diff_cov")
    for (cov in covs) {

      model_no_sex_diff_cov <- omxSetParameters(model_no_sex_diff_cov, label = c(paste0("beta_m_", cov), paste0("beta_f_", cov)), free = TRUE, values = 0, newlabels = paste0("beta_", cov))

    }

    fit_no_sex_diff_cov   <- mxTryHardOrdinal(model_no_sex_diff_cov, extraTries = extra_tries)

  }


  # Constrain expected Thresholds to be equal across twin order
  model_equal_order <- mxModel(fit_saturated, name = "sat_5group_equal_order")
  model_equal_order  <- omxSetParameters(model_equal_order, label=c("tmzf1","tmzf2"), free=TRUE, values=svTh, newlabels='tmzf')
  model_equal_order  <- omxSetParameters(model_equal_order, label=c("tdzf1","tdzf2"), free=TRUE, values=svTh, newlabels='tdzf')
  model_equal_order  <- omxSetParameters(model_equal_order, label=c("tmzm1","tmzm2"), free=TRUE, values=svTh, newlabels='tmzm')
  model_equal_order  <- omxSetParameters(model_equal_order, label=c("tdzm1","tdzm2"), free=TRUE, values=svTh, newlabels='tdzm')
  fit_equal_order   <- mxTryHardOrdinal(model_equal_order, extraTries = extra_tries)

  # Constrain expected Thresholds to be equal across twin order and zygosity
  model_equal_order_zyg  <- mxModel(fit_equal_order, name = "sat_5group_equal_order_zyg" )
  model_equal_order_zyg  <- omxSetParameters(model_equal_order_zyg, label = c("tmzf", "tdzf"), free = TRUE, values = svTh, newlabels = "tzf")
  model_equal_order_zyg  <- omxSetParameters(model_equal_order_zyg, label = c("tmzm", "tdzm"), free = TRUE, values = svTh, newlabels = "tzm")
  fit_equal_order_zyg    <- mxTryHardOrdinal(model_equal_order_zyg, extraTries = extra_tries)

  # Constrain expected Thresholds to be equal across twin order and zygosity and SS/OS
  model_equal_order_zyg_ss_os  <- mxModel(fit_equal_order_zyg, name="sat_5group_equal_order_zyg_ss_os")
  model_equal_order_zyg_ss_os  <- omxSetParameters(model_equal_order_zyg_ss_os, label=c("tzf","tdzo1"), free=TRUE, values=svTh, newlabels="tzf")
  model_equal_order_zyg_ss_os  <- omxSetParameters(model_equal_order_zyg_ss_os, label=c("tzm","tdzo2"), free=TRUE, values=svTh, newlabels="tzm")
  fit_equal_order_zyg_ss_os    <- mxTryHardOrdinal(model_equal_order_zyg_ss_os, extraTries = extra_tries)

  # Constrain expected Thresholds to be equal across twin order and zygosity and SS/OS and sex
  model_equal_order_zyg_ss_os_sex  <- mxModel(fit_equal_order_zyg_ss_os, name="sat_5group_equal_order_zyg_ss_os_sex" )
  model_equal_order_zyg_ss_os_sex  <- omxSetParameters(model_equal_order_zyg_ss_os_sex, label=c("tzf","tzm"), free = TRUE, values = svTh, newlabels = "tZ")
  fit_equal_order_zyg_ss_os_sex    <- mxTryHardOrdinal(model_equal_order_zyg_ss_os_sex, extraTries = extra_tries)

  out <- list(sat = fit_saturated,
              equal_order = fit_equal_order,
              equal_order_zyg = fit_equal_order_zyg,
              equal_order_zyg_ss_os = fit_equal_order_zyg_ss_os,
              equal_order_zyg_ss_os_sex = fit_equal_order_zyg_ss_os_sex,
              response_type = x$response_type,
              trait = x$trait)

  if (!is_null(covs)) {

    out$no_cov = fit_no_cov
    out$no_sex_diff_cov = fit_no_sex_diff_cov

  }

  class(out) <- "saturated.5group.binary"
  out


}



#' @export
fit_saturated.prep.uni.5group.num <- function(x, covs, extra_tries = 10, ...) {


  # Select Variables for Analysis
  vars      <- "X"
  nv        <- 1
  ntv       <- nv * 2
  selVars   <- paste(vars, c(rep(1, nv), rep(2, nv)), sep="")


  # Set Starting Values
  svB <- 20 # start value for regressions
  svMe <- 5
  svVa      <- 0.8                       # start value for variance
  lbVa      <- 0.0001                    # start value for lower bounds

  # Definition variables
  defs1 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "1"), name = "defs1")
  defs2 <- mxMatrix(type = "Full", nrow = 1, ncol = length(covs), free = FALSE, labels = paste0("data.", covs, "2"), name = "defs2")

  # Regression parameters
  path_bm  <- mxMatrix(type = "Full", nrow = length(covs), ncol = 1, free = TRUE, values = rep(svB, length(covs)), label = paste0("beta_m_", covs), name = "bm")
  path_bf  <- mxMatrix(type = "Full", nrow = length(covs), ncol = 1, free = TRUE, values = rep(svB, length(covs)), label = paste0("beta_f_", covs), name = "bf")

  # Create Algebra for expected Mean Matrices
  mean_mzf   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZf1","mMZf2"), name="mean_mzf")
  mean_dzf   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZf1","mDZf2"), name="mean_dzf")
  mean_mzm   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mMZm1","mMZm2"), name="mean_mzm")
  mean_dzm   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZm1","mDZm2"), name="mean_dzm")
  mean_dzo   <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mDZo1","mDZo2"), name="mean_dzo")

  expMean_mzf <- mxAlgebra(expression = mean_mzf + cbind(defs1 %*% bf, defs2 %*% bf),
                           name = "expMean_mzf")

  expMean_dzf <- mxAlgebra(expression = mean_dzf + cbind(defs1 %*% bf, defs2 %*% bf),
                           name = "expMean_dzf")

  expMean_mzm <- mxAlgebra(expression = mean_mzm + cbind(defs1 %*% bm, defs2 %*% bm),
                           name = "expMean_mzm")

  expMean_dzm <- mxAlgebra(expression = mean_dzm + cbind(defs1 %*% bm, defs2 %*% bm),
                           name = "expMean_dzm")

  expMean_dzo <- mxAlgebra(expression = mean_dzo + cbind(defs1 %*% bf, defs2 %*% bm),
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
  data_mzf <- mxData(observed = as.data.frame(x$mzf), type = "raw")
  data_dzf <- mxData(observed = as.data.frame(x$dzf), type = "raw")
  data_mzm <- mxData(observed = as.data.frame(x$mzm), type = "raw")
  data_dzm <- mxData(observed = as.data.frame(x$dzm), type = "raw")
  data_dzo <- mxData(observed = as.data.frame(x$dzo), type = "raw")

  # Create expectation objects
  exp_mzf    <- mxExpectationNormal(covariance = "cov_mzf", means = "expMean_mzf", dimnames = selVars)
  exp_dzf   <- mxExpectationNormal(covariance = "cov_dzf", means = "expMean_dzf", dimnames = selVars)
  exp_mzm   <- mxExpectationNormal(covariance = "cov_mzm", means = "expMean_mzm", dimnames = selVars)
  exp_dzm    <- mxExpectationNormal(covariance = "cov_dzm", means = "expMean_dzm", dimnames = selVars)
  exp_dzo    <- mxExpectationNormal(covariance = "cov_dzo", means = "expMean_dzo", dimnames = selVars)
  funML     <- mxFitFunctionML()

  # Create model objects for multiple groups
  pars      <- list(path_bm, path_bf)
  defs      <- list(defs1, defs2)

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

  fit_sat <- mxTryHard(model_sat, exhaustive = TRUE, extraTries = extra_tries)

  if (!is_null(covs)) {


    # Test covariates
    model_no_cov  <- mxModel(fit_sat, name = "sat_5group_no_cov")
    model_no_cov  <- omxSetParameters(model_no_cov, label = c(paste0("beta_m_", covs), paste0("beta_f_", covs)), free = FALSE, values = 0)
    fit_no_cov    <- mxTryHard(model_no_cov, exhaustive = TRUE, extraTries = extra_tries)

    # Test Sex Difference in Covariate
    model_no_sex_diff_cov <- mxModel(fit_sat, name="sat_5group_no_sex_diff_cov")

    for (cov in covs) {

      model_no_sex_diff_cov <- omxSetParameters(model_no_sex_diff_cov, label = c(paste0("beta_m_", cov), paste0("beta_f_", cov)), free = TRUE, values = 0, newlabels = paste0("beta_", cov))

    }

    fit_no_sex_diff_cov   <- mxTryHard(model_no_sex_diff_cov, exhaustive = TRUE, extraTries = extra_tries)

  }

  # Constrain expected Means to be equal across twin order
  model_equal_mean_order  <- mxModel(fit_sat, name="sat_5group_equal_mean_order" )
  model_equal_mean_order  <- omxSetParameters(model_equal_mean_order, label=c("mMZf1","mMZf2"), free=TRUE, values=svMe, newlabels='mMZf' )
  model_equal_mean_order  <- omxSetParameters(model_equal_mean_order, label=c("mDZf1","mDZf2"), free=TRUE, values=svMe, newlabels='mDZf' )
  model_equal_mean_order  <- omxSetParameters(model_equal_mean_order, label=c("mMZm1","mMZm2"), free=TRUE, values=svMe, newlabels='mMZm' )
  model_equal_mean_order  <- omxSetParameters(model_equal_mean_order, label=c("mDZm1","mDZm2"), free=TRUE, values=svMe, newlabels='mDZm' )
  fit_equal_mean_order    <- mxTryHard(model_equal_mean_order, exhaustive = TRUE, extraTries = extra_tries)

  # Constrain expected Means and Variances to be equal across twin order
  model_equal_mean_var_order <- mxModel(fit_equal_mean_order, name="sat_5group_equal_mean_var_order" )
  model_equal_mean_var_order <- omxSetParameters(model_equal_mean_var_order, label=c("vMZf1","vMZf2"), free=TRUE, values=svVa, newlabels='vMZf')
  model_equal_mean_var_order <- omxSetParameters(model_equal_mean_var_order, label=c("vDZf1","vDZf2"), free=TRUE, values=svVa, newlabels='vDZf')
  model_equal_mean_var_order <- omxSetParameters(model_equal_mean_var_order, label=c("vMZm1","vMZm2"), free=TRUE, values=svVa, newlabels='vMZm')
  model_equal_mean_var_order <- omxSetParameters(model_equal_mean_var_order, label=c("vDZm1","vDZm2"), free=TRUE, values=svVa, newlabels='vDZm')
  fit_equal_mean_var_order  <- mxTryHard(model_equal_mean_var_order, exhaustive = TRUE, extraTries = extra_tries)

  # Constrain expected Means and Variances to be equal across twin order and zygosity
  model_equal_mean_var_order_zyg <- mxModel(fit_equal_mean_var_order, name="sat_5group_equal_mean_var_order_zyg" )
  model_equal_mean_var_order_zyg <- omxSetParameters(model_equal_mean_var_order_zyg, label=c("mMZf","mDZf"), free=TRUE, values = svMe, newlabels='mZf')
  model_equal_mean_var_order_zyg <- omxSetParameters(model_equal_mean_var_order_zyg, label=c("vMZf","vDZf"), free=TRUE, values = svVa, newlabels='vZf')
  model_equal_mean_var_order_zyg <- omxSetParameters(model_equal_mean_var_order_zyg, label=c("mMZm","mDZm"), free=TRUE, values = svMe, newlabels='mZm')
  model_equal_mean_var_order_zyg <- omxSetParameters(model_equal_mean_var_order_zyg, label=c("vMZm","vDZm"), free=TRUE, values = svVa, newlabels='vZm')
  fit_equal_mean_var_order_zyg  <- mxTryHard(model_equal_mean_var_order_zyg, exhaustive = TRUE, extraTries = extra_tries)

  # Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS
  model_equal_mean_var_order_zyg_ss_os <- mxModel(fit_equal_mean_var_order_zyg, name="sat_5group_equal_mean_var_order_zyg_ss_os" )
  model_equal_mean_var_order_zyg_ss_os <- omxSetParameters(model_equal_mean_var_order_zyg_ss_os, label=c("mZf","mDZo1"), free=TRUE, values=svMe, newlabels='mZf')
  model_equal_mean_var_order_zyg_ss_os <- omxSetParameters(model_equal_mean_var_order_zyg_ss_os, label=c("vZf","vDZo1"), free=TRUE, values=svVa, newlabels='vZf')
  model_equal_mean_var_order_zyg_ss_os <- omxSetParameters(model_equal_mean_var_order_zyg_ss_os, label=c("mZm","mDZo2"), free=TRUE, values=svMe, newlabels='mZm')
  model_equal_mean_var_order_zyg_ss_os <- omxSetParameters(model_equal_mean_var_order_zyg_ss_os, label=c("vZm","vDZo2"), free=TRUE, values=svVa, newlabels='vZm')
  fit_equal_mean_var_order_zyg_ss_os   <- mxTryHard(model_equal_mean_var_order_zyg_ss_os, exhaustive = TRUE, extraTries = extra_tries)
  # Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS and sex
  model_equal_mean_var_order_zyg_ss_os_sex <- mxModel(fit_equal_mean_var_order_zyg_ss_os, name = "sat_5group_equal_mean_var_order_zyg_ss_os_sex")
  model_equal_mean_var_order_zyg_ss_os_sex <- omxSetParameters(model_equal_mean_var_order_zyg_ss_os_sex, label=c("mZf","mZm"), free=TRUE, values=svMe, newlabels='mZ')
  model_equal_mean_var_order_zyg_ss_os_sex <- omxSetParameters(model_equal_mean_var_order_zyg_ss_os_sex, label=c("vZf","vZm"), free=TRUE, values=svVa, newlabels='vZ')
  fit_equal_mean_var_order_zyg_ss_os_sex   <- mxTryHard(model_equal_mean_var_order_zyg_ss_os_sex, exhaustive = TRUE, extraTries = extra_tries)


  out <- list(sat = fit_sat,
              equal_mean_order = fit_equal_mean_order,
              equal_mean_var_order = fit_equal_mean_var_order,
              equal_mean_var_order_zyg = fit_equal_mean_var_order_zyg,
              equal_mean_var_order_zyg_ss_os = fit_equal_mean_var_order_zyg_ss_os,
              equal_mean_var_order_zyg_ss_os_sex = fit_equal_mean_var_order_zyg_ss_os_sex,
              response_type = x$response_type,
              trait = x$trait)

  if (!is_null(covs)) {

    out$no_cov <- fit_no_cov
    out$no_sex_diff_cov <- fit_no_sex_diff_cov


  }

  class(out) <- "saturated.5group.num"
  out


}


