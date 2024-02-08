

#' @export
compare <- function(x, ...) {

  UseMethod("compare")

}


#' @export
compare.ACE.uni <- function(x) {

  mxCompare(x$ACE, list(x$AE, x$CE)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type, Constrained = x$constrained)


}



#' @export
compare.ADE.uni <- function(x) {

  mxCompare(x$ADE, list(x$AE, x$DE)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type, Constrained = x$constrained)


}



#' @export
compare.ACE.5group <- function(x) {

  # Test Significance of Sources of Variance of ACEra/rc model with Qualitative and Quantitative Sex differences
  # Run AEra model
  out_ra <- mxCompare(x$ACEra, list(x$ACErc, x$ACEq, x$ACE)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))


  # Run CErc model
  # Print Comparative Fit Statistics
  out_ra_subs <- mxCompare(x$ACEra, list(x$ACErc, x$AEra, x$CErc)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))


  # Test Significance of Sources of Variance of ACEq model with Quantitative Sex differences
  # Run AEq model
  # Run CEq model
  # Print Comparative Fit Statistics
  out_q_subs <- mxCompare(x$ACEq, list(x$AEq, x$CEq)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))

  # Test Significance of Sources of Variance of ACE model without Sex differences
  # Run AE model
  # Print Comparative Fit Statistics
  out <- mxCompare(x$ACE, list(x$AE, x$CE)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))


  bind_rows(out_ra, out_ra_subs, out_q_subs, out) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)

}


#' @export
compare.ADE.5group <- function(x) {


  out1 <- mxCompare(x$ADEra, list(x$ADErd, x$ADEq, x$ADE)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)

  out2 <- mxCompare(x$ADEra, list(x$ADErd, x$AEra, x$AErd)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))

  # Print Comparative Fit Statistics
  out3 <- mxCompare(x$ADEq, list( x$AEq, x$Eq)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))


  # Print Comparative Fit Statistics
  out4 <- mxCompare(x$ADE, list(x$AE, x$E)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq))

  bind_rows(out1, out2, out3, out4) %>% mutate(Trait = x$trait, Response_type = x$response_type)

}



#' @export
compare.saturated.binary <- function(x) {

  mxCompare(x$sat, list(x$no_cov, x$ETO, x$ETOZ)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)


}

#' @export
compare.saturated.num <- function(x) {

  mxCompare(x$sat, list(x$no_cov, x$EMO, x$EMVO, x$EMVOZ)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)

}

#' @export
compare.saturated.5group.binary <- function(x) {


  mxCompare(x$sat, list(x$no_cov, x$no_sex_diff_cov, x$ETO, x$ETZ, x$ETP, x$ETS)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)


}

#' @export
compare.saturated.5group.num <- function(x) {

  mxCompare(x$sat, list(x$no_cov, x$no_sex_diff_cov, x$EMO, x$EMVO, x$EMVZ, x$EMVP, x$EMVS)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)

}

