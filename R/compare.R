

OK_status <- c("infeasible linear constraint",
               "infeasible non-linear constraint",
               "iteration limit/blue",
               "not convex/red",
               "nonzero gradient/red",
               "bad deriv",
               "internal error",
               "infeasible start")

models_ran <- function(m1, m2) {

  grepl("OK", summary(m1)$statusCode) & grepl("OK", summary(m2)$statusCode)

}

compare_internal <- function(m1, m2) {

  if (models_ran(m1, m2)) {

    mxCompare(m1, m2) |>
      as_tibble() |>
      select(-(fit:SBchisq))


  }

}

#' @export
compare <- function(x, ...) {

  UseMethod("compare")

}


#' @export
compare.ACE.uni <- function(x, ...) {

  bind_rows(compare_helper(x$ACE, x$AE),
            compare_helper(x$ACE, x$CE))


}



#' @export
compare.ADE.uni <- function(x, ...) {

  bind_rows(compare_helper(x$ADE, x$AE),
            compare_helper(x$ADE, x$DE))


}



#' @export
compare.ACE.5group <- function(x, ...) {

  # Test Significance of Sources of Variance of ACEra/rc model with Qualitative and Quantitative Sex differences
  # Run AEra model

  bind_rows(compare_internal(x$ACEra, x$ACEq),
            compare_internal(x$ACEq, x$ACE),
            compare_internal(x$ACErc, x$ACEq),
            compare_internal(x$ACEq, x$ACE),
            compare_internal(x$ACEra, x$AEra),
            compare_internal(x$ACErc, x$CErc),
            compare_internal(x$ACEq, x$AEq),
            compare_internal(x$ACEq, x$CEq),
            compare_internal(x$ACE, x$AE),
            compare_internal(x$ACE, x$CE)) |>
    mutate(Trait = x$trait, Response_type = x$response_type, Constrained = x$constrained)

}


#' @export
compare.ADE.5group <- function(x) {

  out_ra <- mxCompare(x$ADEra, list(x$ADEq, x$ADE)) %>%
    as_tibble()

  out_rd <- mxCompare(x$ADErd, list(x$ADEq, x$ADE)) %>%
    as_tibble()

  out2 <- mxCompare(x$ADEra, x$AEra) %>%
    as_tibble()

  # Print Comparative Fit Statistics
  out3 <- mxCompare(x$ADEq, x$AEq) %>%
    as_tibble()

  # Print Comparative Fit Statistics
  out4 <- mxCompare(x$ADE, x$AE) %>%
    as_tibble()

  bind_rows(out_ra, out_rd,  out2, out3, out4) %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type)


}



#' @export
compare.saturated.binary <- function(x) {

  out_init <- x[!names(x) %in% c("response_type", "trait")] |>
    map(mxCompare) |>
    map_dfr(as_tibble)

  if ("no_cov" %in% names(x)) {

    out_no_cov <- mxCompare(x$sat, x$no_cov) %>%
      as_tibble()

    out_init <- bind_rows(out_init, out_no_cov)

  }

  out1 <- mxCompare(x$sat, x$ETO) %>%
    as_tibble()

  out2 <- mxCompare(x$ETO, x$ETOZ) %>%
    as_tibble()

  out <- bind_rows(out1, out2) %>%
    as_tibble() %>%
    mutate(Trait = x$trait, Response_type = x$response_type) |>
    filter(!is.na(comparison))

  bind_rows(out_init, out) |>
    select(-(fit:SBchisq))



}

#' @export
compare.saturated.num <- function(x) {

  out_init <- x[!names(x) %in% c("response_type", "trait")] |>
    map(mxCompare) |>
    map_dfr(as_tibble)


  out1 <- mxCompare(x$sat, x$EMO) %>%
    as_tibble()

  out2 <- mxCompare(x$EMO, x$EMVO) %>%
    as_tibble()

  out3 <- mxCompare(x$EMVO, x$EMVOZ) %>%
    as_tibble()

  out <- bind_rows(out1, out2, out3) %>%
    as_tibble() %>%
    mutate(Trait = x$trait, Response_type = x$response_type) |>
    filter(!is.na(comparison))

  bind_rows(out_init, out) |>
    select(-(fit:SBchisq))





}

#' @export
compare.saturated.5group.binary <- function(x) {


  out_init <- x[!names(x) %in% c("response_type", "trait")] |>
    map(mxCompare) |>
    map_dfr(as_tibble)


  out1 <- mxCompare(x$sat, list(x$no_cov, x$no_sex_diff_cov, x$equal_order)) %>%
    as_tibble()

  out2 <- mxCompare(x$equal_order, x$equal_order_zyg) %>%
    as_tibble()

  out3 <- mxCompare(x$equal_order_zyg, x$equal_order_zyg_ss_os) %>%
    as_tibble()

  out4 <- mxCompare(x$equal_order_zyg_ss_os, x$equal_order_zyg_ss_os_sex) %>%
    as_tibble()

  out <- bind_rows(out1, out2, out3, out4) %>%
    as_tibble() %>%
    mutate(Trait = x$trait, Response_type = x$response_type) |>
    filter(!is.na(comparison))

  bind_rows(out_init, out) |>
    select(-(fit:SBchisq))


}

#' @export
compare.saturated.5group.num <- function(x) {

  out_init <- x[!names(x) %in% c("response_type", "trait")] |>
    map(mxCompare) |>
    map_dfr(as_tibble)


  out1 <- mxCompare(x$sat, list(x$no_cov, x$no_sex_diff_cov, x$equal_mean_order)) %>%
    as_tibble()

  out2 <- mxCompare(x$equal_mean_order, x$equal_mean_var_order) %>%
    as_tibble()

  out3 <- mxCompare(x$equal_mean_var_order, x$equal_mean_var_order_zyg) %>%
    as_tibble()

  out4 <- mxCompare(x$equal_mean_var_order_zyg, x$equal_mean_var_order_zyg_ss_os) %>%
    as_tibble()

  out5 <- mxCompare(x$equal_mean_var_order_zyg_ss_os, x$equal_mean_var_order_zyg_ss_os_sex) %>%
    as_tibble()


  out <- bind_rows(out1, out2, out3, out4, out5) %>%
    mutate(Trait = x$trait, Response_type = x$response_type) |>
    filter(!is.na(comparison))

  bind_rows(out_init, out) |>
    select(-(fit:SBchisq))


}




#' @export
compare.ACE.biv.chol <- function(x, ...) {

  out_init <- x[!names(x) %in% c("traitX", "traitY")] |>
    map(mxCompare) |>
    map_dfr(as_tibble)

  out1 <- mxCompare(x$ACE, x$X_no_C) |> as_tibble()
  out2 <- mxCompare(x$ACE, x$Y_no_C) |> as_tibble()
  out3 <- mxCompare(x$X_no_C, x$AE) |> as_tibble()
  out4 <- mxCompare(x$Y_no_C, x$AE) |> as_tibble()

  out5 <- mxCompare(x$ACE, x$X_no_A) |> as_tibble()
  out6 <- mxCompare(x$ACE, x$Y_no_A) |> as_tibble()
  out7 <- mxCompare(x$X_no_A, x$CE) |> as_tibble()
  out8 <- mxCompare(x$Y_no_A, x$CE) |> as_tibble()

  out <- bind_rows(out1, out2, out3, out4, out5, out6) %>%
    select(-(fit:SBchisq)) %>%
    mutate(Trait = x$trait, Response_type = x$response_type) |>
    filter(!is.na(comparison))

  bind_rows(out_init, out) |>
    select(-(fit:SBchisq))


}
