

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
