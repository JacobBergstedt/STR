

#' @export
get_estimates <- function(x, ...) {

  UseMethod("get_estimates")

}


#' @export
get_estimates.ACE.uni <- function(x) {


  extract_from_model <- function(m) {


    # Variance components
    A <- as.numeric(mxEval(top.A, m))
    C <- as.numeric(mxEval(top.C, m))
    E <- as.numeric(mxEval(top.E, m))

    E_sd <- as.numeric(mxSE(top.E, m))
    A_sd <- as.numeric(mxSE(top.A, m))
    C_sd <- as.numeric(mxSE(top.C, m))


    # Heritability
    h2 <- as.numeric(mxEval(top.A / top.Vtot, m))
    h2_SD <- as.numeric(mxSE(top.A / top.Vtot, m))
    c2 <- as.numeric(mxEval(top.C / top.Vtot, m))
    c2_SD <- as.numeric(mxSE(top.C / top.Vtot, m))
    e2 <- as.numeric(mxEval(top.E / top.Vtot, m))
    e2_SD <- as.numeric(mxSE(top.E / top.Vtot, m))

    tibble(Model = m$name,
           A = A,
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
           e2_SD = e2_SD)


  }

  out_info <- tibble(Trait = x$trait,
                     Response_type = x$response_type,
                     Constrained = x$constrained)


  out <- x[!names(x) %in% c("response_type", "trait", "constrained")] %>%
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}


#' @export
get_estimates.ADE.uni <- function(x) {


  extract_from_model <- function(m) {


    # Variance components
    A <- as.numeric(mxEval(top.A, m))
    D <- as.numeric(mxEval(top.C, m))
    E <- as.numeric(mxEval(top.E, m))

    A_sd <- as.numeric(mxSE(top.A, m))
    D_sd <- as.numeric(mxSE(top.C, m))
    E_sd <- as.numeric(mxSE(top.E, m))




    # Heritability
    h2 <- as.numeric(mxEval(top.A / top.Vtot, m))
    h2_SD <- as.numeric(mxSE(top.A / top.Vtot, m))
    d2 <- as.numeric(mxEval(top.C / top.Vtot, m))
    d2_SD <- as.numeric(mxSE(top.C / top.Vtot, m))
    e2 <- as.numeric(mxEval(top.E / top.Vtot, m))
    e2_SD <- as.numeric(mxSE(top.E / top.Vtot, m))

    tibble(Model = m$name,
           A = A,
           D = D,
           E = E,
           A_SD = A_sd,
           D_SD = D_sd,
           E_SD = E_sd,
           h2 = h2,
           h2_SD = h2_SD,
           d2 = d2,
           d2_SD = d2_SD,
           e2 = e2,
           e2_SD = e2_SD)


  }

  out_info <- tibble(Trait = x$trait,
                     Response_type = x$response_type,
                     Constrained = x$constrained)


  out <- x[!names(x) %in% c("response_type", "trait", "constrained")] %>%
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


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
