


# Internal functions ------------------------------------------------------

fit_ACE_5group_sub_models <- function(fitACEra, svPe, svPa) {


  modelACErc <- mxModel(fitACEra, name = "ACE_5group_rc" )
  modelACErc <- omxSetParameters(modelACErc, labels=c("VAms11"), free=FALSE, values=0 )
  modelACErc <- omxSetParameters(modelACErc, labels=c("VCms11"), free=TRUE, values=0.01 )
  modelACErc <- omxSetParameters(modelACErc, labels=c("VEf11","VEm11"), free=TRUE, values=svPe)
  fitACErc   <- mxTryHard(modelACErc)

  # Run ACEq model - Quantitative non-scalar Sex Differences ACE model
  modelACEq <- mxModel( fitACEra, name="ACE_5group_qa" )
  modelACEq <- omxSetParameters(modelACEq, labels = c("VAms11"), free=FALSE, values=0)
  modelACEq <- omxSetParameters(modelACEq, labels = c("VEf11", "VEm11"), free =TRUE, values = svPe)
  fitACEq   <- mxTryHard(modelACEq)

  # Run ACE model - No Sex differences ACE model
  modelACE  <- mxModel(fitACEq, name = "ACE_5group")
  modelACE  <- omxSetParameters(modelACE, labels = c("VAf11", "VAm11"), free = TRUE, values = svPa, newlabels = 'VA11' )
  modelACE  <- omxSetParameters(modelACE, labels = c("VCf11", "VCm11"), free = TRUE, values = svPa, newlabels = 'VC11' )
  modelACE  <- omxSetParameters(modelACE, labels = c("VEf11", "VEm11"), free = TRUE, values = svPe, newlabels = 'VE11' )
  fitACE    <- mxTryHard(modelACE)

  # Test Significance of Sources of Variance of ACEra/rc model with Qualitative and Quantitative Sex differences
  # Run AEra model
  modelAEra <- mxModel(fitACEra, name="AE_5group_ra")
  modelAEra <- omxSetParameters(modelAEra, labels=c("VCf11", "VCm11"), free = FALSE, values = 0)
  fitAEra   <- mxTryHard(modelAEra)

  # Run CErc model
  modelCErc <- mxModel(fitACErc, name="CE_5group_rc" )
  modelCErc <- omxSetParameters( modelCErc, labels=c("VAf11","VAm11"), free=FALSE, values=0 )
  modelCErc <- omxSetParameters( modelCErc, labels=c("VCf11","VCm11"), free=TRUE, values=svPa )
  modelCErc <- omxSetParameters( modelCErc, labels=c("VEf11","VEm11"), free=TRUE, values=svPe )
  fitCErc   <- mxTryHard(modelCErc)


  # Test Significance of Sources of Variance of ACEq model with Quantitative Sex differences
  # Run AEq model
  modelAEq  <- mxModel(fitACEq, name="AE_5group_qa" )
  modelAEq  <- omxSetParameters(modelAEq, labels=c("VCf11","VCm11"), free=FALSE, values=0 )
  fitAEq    <- mxTryHard(modelAEq)

  # Run CEq model
  modelCEq  <- mxModel( fitACEq, name="CE_5group_qa" )
  modelCEq  <- omxSetParameters( modelCEq, labels=c("VAf11","VAm11"), free=FALSE, values=0 )
  modelCEq  <- omxSetParameters( modelCEq, labels=c("VCf11","VCm11"), free=TRUE, values=svPa )
  modelCEq  <- omxSetParameters( modelCEq, labels=c("VEf11","VEm11"), free=TRUE, values=svPe )
  fitCEq    <- mxTryHard( modelCEq)

  # Test Significance of Sources of Variance of ACE model without Sex differences
  # Run AE model
  modelAE   <- mxModel( fitACE, name="AE_5group" )
  modelAE   <- omxSetParameters( modelAE, labels=c("VC11"), free=FALSE, values=0 )
  fitAE     <- mxTryHard( modelAE)

  # Run CE model
  modelCE   <- mxModel(fitACE, name="CE_5group" )
  modelCE   <- omxSetParameters(modelCE, labels=c("VA11"), free=FALSE, values=0)
  fitCE     <- mxTryHard( modelCE)

  list(
    ACEra = fitACEra,
    ACErc = fitACErc,
    ACEq = fitACEq,
    ACE = fitACE,
    AEra = fitAEra,
    CErc = fitCErc,
    AEq = fitAEq,
    CEq = fitCEq,
    AE = fitAE,
    CE = fitCE
  )




}


# Outward facing functions ------------------------------------------------



#' @export
fit_ACE <- function(x, ...) {

  UseMethod("fit_ACE")

}



# Univariate --------------------------------------------------------------



#' @export
fit_ACE.prep.uni <- function(x, covs = c("Female", "Birth_year_first", "Birth_year_second"), constrained = TRUE) {



  if (constrained) {

    m <- umxACE(
      name = "ACE_ICD",
      selDVs = "X",
      selCovs = covs,
      sep = "",
      dzData = x$DZ,
      mzData = x$MZ,
      opt = "NPSOL",
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
      dzData = x$DZ,
      mzData = x$MZ,
      autoRun = FALSE,
      opt = "NPSOL",
      addCI = FALSE)

    m <- mxTryHard(m)




    m_AE    <- umxModify(m, update = "C_r1c1", name = "AE", tryHard = "yes")
    m_CE    <- umxModify(m, update = "A_r1c1", name = "CE", tryHard = "yes")


  }

  out <- list(ACE = m, AE = m_AE, CE = m_CE, response_type = x$response_type, trait = x$trait, constrained = constrained)
  class(out) <- c("ACE.uni")
  out


}



# 5 group ----------------------------------------------------------------



#' @export
fit_ACE.prep.uni.5group.binary <- function(x) {

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
  covAf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VAf11", name="VAf" )
  covCf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VCf11", name="VCf" )
  covEf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VEf11", name="VEf" )
  covAm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VAm11", name="VAm" )
  covCm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VCm11", name="VCm" )
  covEm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VEm11", name="VEm" )
  covAms    <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="VAms11", name="VAms" )
  covCms    <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="VCms11", name="VCms" )

  # Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
  signA     <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA")
  signC     <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")

  # Calculate absolute covariation between Males and Females and then un-absoluting the product
  covAos    <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
  covCos    <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")

  # Calculate rg/rc from reparameterized model
  pathRg    <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
  pathRc    <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covPf     <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
  covPm     <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
  covMZf    <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
  covDZf    <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
  covMZm    <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
  covDZm    <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
  covDZo    <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )

  expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
  expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
  expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
  expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
  expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

  # Constrain Variance of Binary Variables
  var1f     <- mxConstraint( expression=diag2vec(Vf)==1, name="Var1f" )
  var1m     <- mxConstraint( expression=diag2vec(Vm)==1, name="Var1m" )


  # Create data objects
  data_mzf <- mxData(observed = x$mzf, type = "raw")
  data_mzm <- mxData(observed = x$mzm, type = "raw")
  data_dzf <- mxData(observed = x$dzf, type = "raw")
  data_dzm <- mxData(observed = x$dzm, type = "raw")
  data_dzo <- mxData(observed = x$dzo, type = "raw")

  # Create Expectation Objects for Multiple Groups
  expMZf    <- mxExpectationNormal( covariance="expCovMZf", means="expMean_zf", dimnames=selVars, thresholds="threGf" )
  expDZf    <- mxExpectationNormal( covariance="expCovDZf", means="expMean_zf", dimnames=selVars, thresholds="threGf" )
  expMZm    <- mxExpectationNormal( covariance="expCovMZm", means="expMean_zm", dimnames=selVars, thresholds="threGm" )
  expDZm    <- mxExpectationNormal( covariance="expCovDZm", means="expMean_zm", dimnames=selVars, thresholds="threGm" )
  expDZo    <- mxExpectationNormal( covariance="expCovDZo", means="expMean_zo", dimnames=selVars, thresholds="threGo" )
  funML     <- mxFitFunctionML()

  # Create Model Objects for Multiple Groups

  parsZf    <- list(path_bf_birthyear_first, path_bf_birthyear_second, meanG, threGf, covAf, covCf, covEf, covPf )
  parsZm    <- list(path_bm_birthyear_first,  path_bm_birthyear_second, meanG, threGm, covAm, covCm, covEm, covPm, covAms, covCms)
  parsZo    <- list(parsZm, parsZf, meanG, threGo, signA, signC, covAos, covCos, pathRg, pathRc)
  defs      <- list(def_birth_year_first, def_birth_year_second)

  modelMZf  <- mxModel( parsZf, defs, expMean_zf, covMZf, expCovMZf, data_mzf, expMZf, funML, name="MZf" )
  modelMZm  <- mxModel( parsZm, defs, expMean_zm, covMZm, expCovMZm, data_mzm, expMZm, funML, name="MZm" )
  modelDZf  <- mxModel( parsZf, defs, expMean_zf, covDZf, expCovDZf, data_dzf, expDZf, funML, name="DZf" )
  modelDZm  <- mxModel( parsZm, defs, expMean_zm, covDZm, expCovDZm, data_dzm, expDZm, funML, name="DZm" )
  modelDZo  <- mxModel( parsZo, defs, expMean_zo, covDZo, expCovDZo, data_dzo, expDZo, funML, name="DZo" )
  multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

  # Create Algebra for Variance Components
  rowUS     <- rep('US', nv)
  colUS     <- rep(c('VAf','VCf','VEf','SAf','SCf','SEf','VAm','VCm','VEm','SAm','SCm','SEm','rg','rc'), each=nv)
  estUS     <- mxAlgebra(expression= cbind(VAf,VCf,VEf,VAf/Vf,VCf/Vf,VEf/Vf,VAm+VAms,VCm+VCms,VEm,(VAm+VAms)/Vm,(VCm+VCms)/Vm,VEm/Vm,rg,rc), name="US", dimnames=list(rowUS, colUS))

  # Build Model with Confidence Intervals
  modelACEra <- mxModel("ACE_5group_ra", parsZf, parsZm, parsZo, var1f, var1m, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, estUS)

  # RUN MODEL

  # Run ACEra Model - Qualitative (Ra) & Quantative Sex Differences ACE model
  fitACEra  <- mxTryHard(modelACEra)
  out <- fit_ACE_5group_sub_models(fitACEra, svPe = svPe, svPa = svPa)
  out$response_type <- x$response_type
  out$trait <- x$trait
  class(out) <- "ACE.5group"
  out

}





#' @export
fit_ACE.prep.uni.5group.num <- function(x) {

  nv        <- 1                         # number of variables
  ntv       <- nv*2                      # number of total variables
  selVars   <- paste("X", c(rep(1, nv), rep(2, nv)),sep="")


  svB_birth_year_first <- 0.05 # start value for regressions
  svB_birth_year_second <- 0.05
  svMe <- 1
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
  meanGf    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mZf","mZf"), name="meanGf")
  meanGm    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mZm","mZm"), name="meanGm")
  meanGo    <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=c("mZf","mZm"), name="meanGo")

  meanG <- mxMatrix(type = "Zero", nrow = 1, ncol = ntv, name = "meanG" )

  expMean_zf <- mxAlgebra(expression = meanGf +
                            cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                  def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second),
                          name = "expMean_zf")


  expMean_zm <- mxAlgebra(expression = meanGm +
                            cbind(def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second,
                                  def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                          name = "expMean_zm")

  expMean_zo <- mxAlgebra(expression = meanGo +
                            cbind(def_birth_year_first %*% bf_birth_year_first + def_birth_year_second %*% bf_birth_year_second,
                                  def_birth_year_first %*% bm_birth_year_first + def_birth_year_second %*% bm_birth_year_second),
                          name = "expMean_zo")

  # Create Matrices for expected Direct Symmetric Covariances
  covAf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VAf11", name="VAf" )
  covCf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VCf11", name="VCf" )
  covEf     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VEf11", name="VEf" )
  covAm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VAm11", name="VAm" )
  covCm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VCm11", name="VCm" )
  covEm     <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPe, label="VEm11", name="VEm" )
  covAms    <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=0, label="VAms11", name="VAms" )
  covCms    <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=0, label="VCms11", name="VCms" )

  # Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
  signA     <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA")
  signC     <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")

  # Calculate absolute covariation between Males and Females and then un-absoluting the product
  covAos    <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
  covCos    <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")

  # Calculate rg/rc from reparameterized model
  pathRg    <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
  pathRc    <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covPf     <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
  covPm     <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
  covMZf    <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
  covDZf    <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
  covMZm    <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
  covDZm    <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
  covDZo    <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )
  expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
  expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
  expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
  expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
  expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

  # Create Data Objects for Multiple Groups
  dataMZf   <- mxData(observed = x$mzf, type = "raw")
  dataDZf   <- mxData(observed = x$dzf, type = "raw")
  dataMZm   <- mxData(observed = x$mzm, type = "raw")
  dataDZm   <- mxData(observed = x$dzm, type = "raw")
  dataDZo   <- mxData(observed = x$dzo, type = "raw")

  # Create Expectation Objects for Multiple Groups
  expMZf    <- mxExpectationNormal(covariance="expCovMZf", means="expMean_zf", dimnames = selVars)
  expDZf    <- mxExpectationNormal(covariance="expCovDZf", means="expMean_zf", dimnames = selVars)
  expMZm    <- mxExpectationNormal(covariance="expCovMZm", means="expMean_zm", dimnames = selVars)
  expDZm    <- mxExpectationNormal(covariance="expCovDZm", means="expMean_zm", dimnames = selVars)
  expDZo    <- mxExpectationNormal(covariance="expCovDZo", means="expMean_zo", dimnames = selVars)
  funML     <- mxFitFunctionML()

  # Regression parameters

  path_bm_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_m_birth_year_first", name = "bm_birth_year_first")
  path_bm_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_m_birth_year_second", name = "bm_birth_year_second")

  path_bf_birthyear_first  <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_first, label = "beta_f_birth_year_first", name = "bf_birth_year_first")
  path_bf_birthyear_second <- mxMatrix(type = "Full", nrow = 1, ncol = 1, free = TRUE, values = svB_birth_year_second, label = "beta_f_birth_year_second", name = "bf_birth_year_second")


  # Create Model Objects for Multiple Groups
  parsZf    <- list( path_bf_birthyear_first, path_bf_birthyear_second, covAf, covCf, covEf, covPf )
  parsZm    <- list( path_bm_birthyear_first, path_bm_birthyear_second, covAm, covCm, covEm, covPm, covAms, covCms )
  parsZo    <- list( parsZm, parsZf, signA, signC, covAos, covCos, pathRg, pathRc )
  defs      <- list( def_birth_year_first, def_birth_year_second)
  modelMZf  <- mxModel(parsZf, defs, meanGf, expMean_zf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf")
  modelDZf  <- mxModel(parsZf, defs, meanGf, expMean_zf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf")
  modelMZm  <- mxModel(parsZm, defs, meanGm, expMean_zm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm")
  modelDZm  <- mxModel(parsZm, defs, meanGm, expMean_zm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm")
  modelDZo  <- mxModel(parsZo, defs, meanGo, expMean_zo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo")
  multi     <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

  # Create Algebra for Variance Components
  rowUS     <- rep('US',nv)
  colUS     <- rep(c('VAf','VCf','VEf','SAf','SCf','SEf','VAm','VCm','VEm','SAm','SCm','SEm','rg','rc'),each=nv)
  estUS     <- mxAlgebra( expression=cbind(VAf,VCf,VEf,VAf/Vf,VCf/Vf,VEf/Vf,VAm+VAms,VCm+VCms,VEm,(VAm+VAms)/Vm,(VCm+VCms)/Vm,VEm/Vm,rg,rc), name="US", dimnames=list(rowUS,colUS))

  # Build Model with Confidence Intervals
  modelACEra <- mxModel( "ACE_5group_ra", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, estUS)

#
#   # RUN MODEL
#
  # Run ACEra Model - Qualitative (Ra) & Quantative Sex Differences ACE model
  fitACEra  <- mxTryHard(modelACEra)
  out <- fit_ACE_5group_sub_models(fitACEra, svPe = svPe, svPa = svPa)
  out$response_type <- x$response_type
  out$trait <- x$trait
  class(out) <- "ACE.5group"
  out

}



# Bivariate ---------------------------------------------------------------



#' @export
fit_ACE.prep.biv <- function(x, covs = c("Female", "Birth_year_first", "Birth_year_second")) {

  m <- umxACE(
    name = "ACE_biv",
    selDVs = c("X", "Y"),
    opt = "NPSOL",
    selCovs = c("Female", "Birth_year_first"),
    sep = "",
    dzData = as.data.frame(x$DZ),
    mzData = as.data.frame(x$MZ),
    tryHard = "yes",
    addCI = FALSE,
    autoRun = FALSE
  )

  ACE <- mxTryHard(m)

  out <- list(ACE = ACE, response_type = x$response_type, traitX = x$traitX, traitY = x$traitY)
  class(out) <- c("ACE.biv")
  out


}
