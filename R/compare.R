

#' @export
compare <- function(x, ...) {

  UseMethod("compare")

}


#' @export
compare.ACE.uni <- function(x) {

  mxCompare(x$ACE, list(x$AE, x$CE, x$E)) %>%
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
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)


  # Run CErc model
  # Print Comparative Fit Statistics
  out_ra_subs <- mxCompare(x$ACEra, list(x$ACErc, x$AEra, x$CErc)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)


  # Test Significance of Sources of Variance of ACEq model with Quantitative Sex differences
  # Run AEq model
  # Run CEq model
  # Print Comparative Fit Statistics
  out_q_subs <- mxCompare(x$ACEq, list(x$AEq, x$CEq)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)

  # Test Significance of Sources of Variance of ACE model without Sex differences
  # Run AE model
  # Print Comparative Fit Statistics
  out <- mxCompare(x$ACE, list(x$AE, x$CE)) %>%
    as_tibble() %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)


  bind_rows(out_ra, out_ra_subs, out_q_subs, out)

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

