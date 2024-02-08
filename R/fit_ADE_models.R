


# Internal functions ------------------------------------------------------

fit_ADE_5group_sub_models <- function(fitADEra, svPe, svPa) {

  # Run ADErd Model - Qualitative (Rd) & Quantative Sex Differences ADE model
  modelADErd <- mxModel(fitADEra, name = "ADE_5group_rd" )
  modelADErd <- omxSetParameters(modelADErd, labels = c("VAms11"), free=FALSE, values=0)
  modelADErd <- omxSetParameters(modelADErd, labels = c("VDms11"), free=TRUE, values=0.01)
  modelADErd <- omxSetParameters(modelADErd, labels = c("VEf11","VEm11"), free=TRUE, values=svPe)
  fitADErd   <- mxTryHard(modelADErd)

  # Run ADEq model - Quantitative non-scalar Sex Differences ADE model
  modelADEq <- mxModel(fitADEra, name="ADE_5group_qa")
  modelADEq <- omxSetParameters(modelADEq, labels=c("VAms11"), free=FALSE, values=0)
  modelADEq <- omxSetParameters(modelADEq, labels=c("VEf11","VEm11"), free=TRUE, values=svPe)
  fitADEq   <- mxTryHard(modelADEq)

  # Run ADE model - No Sex differences ADE model
  modelADE  <- mxModel(fitADEq, name = "ADE_5group")
  modelADE  <- omxSetParameters(modelADE, labels = c("VAf11", "VAm11"), free = TRUE, values = svPa, newlabels = "VA11")
  modelADE  <- omxSetParameters(modelADE, labels = c("VDf11", "VDm11"), free = TRUE, values = svPa, newlabels = "VD11")
  modelADE  <- omxSetParameters(modelADE, labels = c("VEf11", "VEm11"), free = TRUE, values = svPe, newlabels = "VE11")
  fitADE   <- mxTryHard(modelADE)

  # Test Significance of Sources of Variance of ADEra/rd model with Qualitative and Quantitative Sex differences
  # Run AEra model
  modelAEra <- mxModel(fitADEra, name = "AE_5group_ra")
  modelAEra <- omxSetParameters( modelAEra, labels=c("VDf11","VDm11"), free=FALSE, values=0 )
  fitAEra   <- mxTryHard(modelAEra)

  # Run AErd model
  modelAErd <- mxModel(fitADErd, name = "AE_5group_rd" )
  modelAErd <- omxSetParameters(modelAErd, labels = c("VDf11","VDm11"), free = FALSE, values = 0)
  fitAErd   <- mxTryHard(modelAErd)

  # Test Significance of Sources of Variance of ADEq model with Quantitative Sex differences
  # Run AEq model
  modelAEq  <- mxModel(fitADEq, name="AE_5group_qa" )
  modelAEq  <- omxSetParameters(modelAEq, labels = c("VDf11","VDm11"), free = FALSE, values = 0)
  fitAEq    <- mxTryHard(modelAEq)

  # Run Eq model
  modelEq   <- mxModel(fitAEq, name = "E_5group_qa")
  modelEq   <- omxSetParameters(modelEq, labels = c("VAf11","VAm11"), free = FALSE, values = 0)
  modelEq   <- omxSetParameters(modelEq, labels = c("VEf11","VEm11"), free = TRUE, values = svPe)
  fitEq     <- mxTryHard(modelEq)

  # Test Significance of Sources of Variance of ADE model without Sex differences
  # Run AE model
  modelAE   <- mxModel(fitADE, name = "AE_5group" )
  modelAE   <- omxSetParameters(modelAE, labels=c("VD11"), free=FALSE, values=0)
  fitAE     <- mxTryHard(modelAE)

  # Run E model
  modelE    <- mxModel( fitAE, name="E_5group_qa")
  modelE    <- omxSetParameters(modelE, labels = c("VA11"), free=FALSE, values = 0)
  fitE      <- mxTryHard(modelE)



  list(
    ADEra = fitADEra,
    ADErc = fitADErd,
    ADEq = fitADEq,
    ADE = fitADE,
    AEra = fitAEra,
    AErd = fitAErd,
    AEq = fitAEq,
    Eq = fitEq,
    AE = fitAE,
    E = fitE
  )



}




# Outward facing functions ------------------------------------------------



#' @export
fit_ADE <- function(x, ...) {

  UseMethod("fit_ADE")

}



#' @export
fit_ADE.prep.uni <- function(x, covs = c("Female", "Birth_year_first", "Birth_year_second"), constrained = TRUE) {



  if (constrained) {

    m <- umxACE(
      name = "ACE_ICD",
      selDVs = "X",
      selCovs = covs,
      sep = "",
      dzData = x$DZ,
      mzData = x$MZ,
      opt = "NPSOL",
      dzCr=.25,
      autoRun = FALSE,
      addCI = FALSE
    )

    m <- mxTryHard(m)

    # AE Model
    m_AE <- umxModify(m, update = "c_r1c1", name = "AE", tryHard = "yes")

    # CE Model
    m_CE <- umxModify(m, update = "a_r1c1", name = "CE", tryHard = "yes")


  } else {

    m <- umxACEv(
      name = "ACE_ICD",
      selDVs = "X",
      selCovs = covs,
      sep = "",
      dzCr =.25,
      dzData = x$DZ,
      mzData = x$MZ,
      autoRun = FALSE,
      opt = "NPSOL",
      addCI = FALSE)

    m <- mxTryHard(m)



    # AE Model
    m_AE    <- umxModify(m, update = "C_r1c1", name = "AE", tryHard = "yes")

    # CE Model
    m_DE    <- umxModify(m, update = "A_r1c1", name = "CE", tryHard = "yes")


  }

  out <- list(ACE = m, AE = m_AE, DE = m_DE, response_type = x$response_type, trait = x$trait, constrained = constrained)
  class(out) <- c("ADE.uni")
  out


}



# 5 group ----------------------------------------------------------------



#' @export
fit_ADE.prep.uni.5group.binary <- function(x) {

  nv        <- 1                         # number of variables
  ntv       <- nv*2                      # number of total variables
  selVars   <- paste("X", c(rep(1, nv), rep(2, nv)),sep="")


  svB_birth_year_first <- 0.05 # start value for regressions
  svB_birth_year_second <- 0.05

  svTh      <- 0.8                       # start value for thresholds
  svCor     <- 0.5                       # start value for correlations
  lbCor     <- -0.99                     # lower bounds for correlations
  ubCor     <- 0.99                      # upper bounds for correlations

  svPa      <- 0.4                       # start value for path coefficient
  svPe      <- 0.8                       # start value for path coefficient for e


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


  threGf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svTh, labels=c("tZf","tZf"), name="threGf" )
  threGm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svTh, labels=c("tZm","tZm"), name="threGm" )
  threGo    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svTh, labels=c("tZf","tZm"), name="threGo" )

  # Create Matrices for expected Direct Symmetric Covariances
  covAf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VAf11", name="VAf")
  covDf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VDf11", name="VDf")
  covEf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VEf11", name="VEf")
  covAm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VAm11", name="VAm")
  covDm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VDm11", name="VDm")
  covEm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VEm11", name="VEm")
  covAms    <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="VAms11", name="VAms")
  covDms    <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="VDms11", name="VDms")

  # Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
  signA     <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA")
  signD     <- mxAlgebra( ((-1)^omxLessThan(VDf,0))*((-1)^omxLessThan(VDm,0)), name="signD")

  # Calculate absolute covariation between Males and Females and then un-absoluting the product
  covAos    <- mxAlgebra(signA * (sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name = "VAos")
  covDos    <- mxAlgebra(signD * (sqrt(abs(VDf))*t(sqrt(abs(VDm)))), name = "VDos")

  # Calculate rg/rc from reparameterized model
  pathRg    <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
  pathRd    <- mxAlgebra( signD*(sqrt(abs(VDf))*t(sqrt(abs(VDm))))/sqrt(VDf*(VDm+VDms)), name="rd")

  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covPf     <- mxAlgebra(expression = VAf + VDf + VEf, name = "Vf")
  covPm     <- mxAlgebra(expression = VAm + VDm + VEm, name = "Vm")
  covMZf    <- mxAlgebra(expression = VAf + VDf, name = "cMZf")
  covDZf    <- mxAlgebra(expression = 0.5 %x% VAf + 0.25 %x% VDf, name ="cDZf")
  covMZm    <- mxAlgebra(expression = VAm + VDm, name = "cMZm")
  covDZm    <- mxAlgebra(expression = 0.5 %x% VAm + 0.25%x%VDm, name="cDZm")
  covDZo    <- mxAlgebra(expression = 0.5 %x% VAos + VDos, name = "cDZo")

  expCovMZf <- mxAlgebra(expression = rbind(cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name = "expCovMZf")
  expCovDZf <- mxAlgebra(expression = rbind(cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name = "expCovDZf")
  expCovMZm <- mxAlgebra(expression = rbind(cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name = "expCovMZm")
  expCovDZm <- mxAlgebra(expression = rbind(cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name = "expCovDZm")
  expCovDZo <- mxAlgebra(expression = rbind(cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name = "expCovDZo")

  # Constrain Variance of Binary Variables
  var1f     <- mxConstraint(expression = diag2vec(Vf) == 1, name = "Var1f")
  var1m     <- mxConstraint(expression = diag2vec(Vm) == 1, name = "Var1m")


  # Create data objects
  data_mzf <- mxData(observed = as.data.frame(x$mzf), type = "raw")
  data_dzf <- mxData(observed = as.data.frame(x$dzf), type = "raw")
  data_mzm <- mxData(observed = as.data.frame(x$mzm), type = "raw")
  data_dzm <- mxData(observed = as.data.frame(x$dzm), type = "raw")
  data_dzo <- mxData(observed = as.data.frame(x$dzo), type = "raw")

  # Create Expectation Objects for Multiple Groups
  expMZf    <- mxExpectationNormal(covariance = "expCovMZf", means = "expMean_zf", dimnames = selVars, thresholds = "threGf" )
  expDZf    <- mxExpectationNormal(covariance = "expCovDZf", means = "expMean_zf", dimnames = selVars, thresholds = "threGf" )
  expMZm    <- mxExpectationNormal(covariance = "expCovMZm", means = "expMean_zm", dimnames=selVars, thresholds = "threGm" )
  expDZm    <- mxExpectationNormal(covariance = "expCovDZm", means = "expMean_zm", dimnames=selVars, thresholds = "threGm" )
  expDZo    <- mxExpectationNormal(covariance = "expCovDZo", means = "expMean_zo", dimnames=selVars, thresholds = "threGo" )
  funML     <- mxFitFunctionML()

  # Create Model Objects for Multiple Groups

  parsZf    <- list(path_bf_birthyear_first, path_bf_birthyear_second, meanG, threGf, covAf, covDf, covEf, covPf )
  parsZm    <- list(path_bm_birthyear_first,  path_bm_birthyear_second, meanG, threGm, covAm, covDm, covEm, covPm, covAms, covDms)
  parsZo    <- list(parsZm, parsZf, meanG, threGo, signA, signD, covAos, covDos, pathRg, pathRd)
  defs      <- list(def_birth_year_first, def_birth_year_second)

  modelMZf  <- mxModel(parsZf, defs, expMean_zf, covMZf, expCovMZf, data_mzf, expMZf, funML, name = "MZf")
  modelDZf  <- mxModel(parsZf, defs, expMean_zf, covDZf, expCovDZf, data_dzf, expDZf, funML, name = "DZf")
  modelMZm  <- mxModel(parsZm, defs, expMean_zm, covMZm, expCovMZm, data_mzm, expMZm, funML, name = "MZm")
  modelDZm  <- mxModel(parsZm, defs, expMean_zm, covDZm, expCovDZm, data_dzm, expDZm, funML, name = "DZm")
  modelDZo  <- mxModel(parsZo, defs, expMean_zo, covDZo, expCovDZo, data_dzo, expDZo, funML, name = "DZo")
  multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

  # Create Algebra for Variance Components
  rowUS     <- rep('US',nv)
  colUS     <- rep(c('VAf','VDf','VEf','SAf','SDf','SEf','VAm','VDm','VEm','SAm','SDm','SEm','rg','rc'), each = nv)
  estUS     <- mxAlgebra(expression = cbind(VAf,VDf,VEf,VAf/Vf,VDf/Vf,VEf/Vf,VAm+VAms,VDm+VDms,VEm,(VAm+VAms)/Vm,(VDm+VDms)/Vm,VEm/Vm,rg,rd), name="US", dimnames=list(rowUS,colUS))

  # Build Model with Confidence Intervals
  modelADEra <- mxModel("ADE_5group_ra", parsZf, parsZm, parsZo, var1f, var1m, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, estUS)

  # RUN MODEL

  # Run ADEra Model - Qualitative (Ra) & Quantative Sex Differences ACE model
  fitADEra  <- mxTryHard(modelADEra)
  out <- fit_ADE_5group_sub_models(fitADEra, svPe = svPe, svPa = svPa)
  out$response_type <- x$response_type
  out$trait <- x$trait
  class(out) <- "ADE.5group"
  out

}


