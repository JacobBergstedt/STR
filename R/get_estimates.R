

#' @export
get_estimates <- function(x, ...) {

  UseMethod("get_estimates")

}


#' @export
get_estimates.ACE.uni <- function(x) {


  ACE <- x$ACE


  # Variance components
  A <- as.numeric(mxEval(top.A, ACE))
  C <- as.numeric(mxEval(top.C, ACE))
  E <- as.numeric(mxEval(top.E, ACE))

  E_sd <- as.numeric(mxSE(top.E, ACE))
  A_sd <- as.numeric(mxSE(top.A, ACE))
  C_sd <- as.numeric(mxSE(top.C, ACE))


  # Heritability
  h2 <- as.numeric(mxEval(top.A / top.Vtot, ACE))
  h2_SD <- as.numeric(mxSE(top.A / top.Vtot, ACE))
  c2 <- as.numeric(mxEval(top.C / top.Vtot, ACE))
  c2_SD <- as.numeric(mxSE(top.C / top.Vtot, ACE))
  e2 <- as.numeric(mxEval(top.E / top.Vtot, ACE))
  e2_SD <- as.numeric(mxSE(top.E / top.Vtot, ACE))

  tibble(A = A,
         C = C,
         E = E,
         A_SD = A_sd,
         C_SD = C_sd,
         E_SD = E_sd,
         h2 = h2,
         h2_SD = h2_SD,
         c2 = c2,
         c2_SD = c2_SD,
         e2 = e2,
         e2_SD = e2_SD,
         Trait = x$trait,
         Response_type = x$response_type,
         Constrained = x$constrained)

}


#' @export
get_estimates.ACE.biv <- function(x) {

  ACE <- x$ACE

  # Covariance decomposition
  bivA <- mxEval(top.A / top.Vtot, ACE)[2, 1]
  bivA_sd <- mxSE(top.A / top.Vtot, ACE)[2, 1]

  bivC <- mxEval(top.C / top.Vtot, ACE)[2, 1]
  bivC_sd <- mxSE(top.C / top.Vtot, ACE)[2, 1]

  bivE <- mxEval(top.E / top.Vtot, ACE)[2, 1]
  bivE_sd <- mxSE(top.E / top.Vtot, ACE)[2, 1]

  # Correlation decomposition
  gen_cor <- mxEval(cov2cor(top.A), ACE)[2, 1]
  gen_cor_sd <- mxSE(cov2cor(top.A), ACE)[2, 1]

  shared_cor <- mxEval(cov2cor(top.C), ACE)[2, 1]
  shared_cor_sd <- mxSE(cov2cor(top.C), ACE)[2, 1]

  non_shared_cor <- mxEval(cov2cor(top.E), ACE)[2, 1]
  non_shared_cor_sd <- mxSE(cov2cor(top.E), ACE)[2, 1]

  # Variance components
  A1 <- mxEval(top.A, ACE)[1, 1]
  C1 <- mxEval(top.C, ACE)[1, 1]
  E1 <- mxEval(top.E, ACE)[1, 1]

  A1_sd <- mxSE(top.A, ACE)[1, 1]
  C1_sd <- mxSE(top.C, ACE)[1, 1]
  E1_sd <- mxSE(top.E, ACE)[1, 1]

  A2 <- mxEval(top.A, ACE)[2, 2]
  C2 <- mxEval(top.C, ACE)[2, 2]
  E2 <- mxEval(top.E, ACE)[2, 2]

  A2_sd <- mxSE(top.A, ACE)[2, 2]
  C2_sd <- mxSE(top.C, ACE)[2, 2]
  E2_sd <- mxSE(top.E, ACE)[2, 2]


  # Heritability
  h2_1 <- mxEval(top.A / top.Vtot, ACE)[1, 1]
  h2_1_sd <- mxSE(top.A / top.Vtot, ACE)[1, 1]
  h2_2 <- mxEval(top.A / top.Vtot, ACE)[2, 2]
  h2_2_sd <- mxSE(top.A / top.Vtot, ACE)[2, 2]


  tibble(traitX = x$traitX,
         traitY = x$traitY,
         bivA = bivA,
         bivC = bivC,
         bivE = bivE,
         bivA_sd = bivA_sd,
         bivC_sd = bivC_sd,
         bivE_sd = bivE_sd,
         cor_gen = gen_cor,
         cor_shared = shared_cor,
         cor_non_shared = non_shared_cor,
         cor_gen_SD = gen_cor_sd,
         cor_shared_SD = shared_cor_sd,
         cor_non_shared_SD = non_shared_cor_sd,
         var_comp_A1 = A1,
         var_comp_C1 = C1,
         var_comp_E1 = E1,
         var_comp_A1_sd = A1_sd,
         var_comp_C1_sd = C1_sd,
         var_comp_E1_sd = E1_sd,
         var_comp_A2 = A2,
         var_comp_C2 = C2,
         var_comp_E2 = E2,
         var_comp_A2_sd = A2_sd,
         var_comp_C2_sd = C2_sd,
         var_comp_E2_sd = E2_sd,
         h2_1 = h2_1,
         h2_2 = h2_2,
         h2_1_sd = h2_1_sd,
         h2_2_sd = h2_2_sd)

}
