
try_or_NA <- function(x) {

  tryCatch(x, error=function(err) NA)

}


#' @export
get_estimates <- function(x, ...) {

  UseMethod("get_estimates")

}


#' @export
get_estimates.ACE.uni <- function(x, ...) {


  extract_from_model <- function(m) {


    model = m$name

    if (model == "ACE") {

      # Variance components
      A <- as.numeric(mxEval(top.A, m))
      C <- as.numeric(mxEval(top.C, m))
      E <- as.numeric(mxEval(top.E, m))

      A_SD <- try_or_NA(as.numeric(mxSE(top.A, m)))
      C_SD <- try_or_NA(as.numeric(mxSE(top.C, m)))
      E_SD <- try_or_NA(as.numeric(mxSE(top.E, m)))


      # Heritability
      h2 <- as.numeric(mxEval(top.A / top.Vtot, m))
      h2_SD <- try_or_NA(as.numeric(mxSE(top.A / top.Vtot, m)))
      c2 <- as.numeric(mxEval(top.C / top.Vtot, m))
      c2_SD <- try_or_NA(as.numeric(mxSE(top.C / top.Vtot, m)))
      e2 <- as.numeric(mxEval(top.E / top.Vtot, m))
      e2_SD <- try_or_NA(as.numeric(mxSE(top.E / top.Vtot, m)))

      tibble(Param = c("A", "C", "E", "h2", "c2", "e2"),
             Value = c(A, C, E, h2, c2, e2),
             SD = c(A_SD, C_SD, E_SD, h2_SD, c2_SD, e2_SD),
             model = m$name)


    } else if (model == "AE") {


      # Variance components
      A <- as.numeric(mxEval(top.A, m))
      E <- as.numeric(mxEval(top.E, m))

      A_SD <- as.numeric(mxSE(top.A, m))
      E_SD <- as.numeric(mxSE(top.E, m))


      # Heritability
      h2 <- as.numeric(mxEval(top.A / top.Vtot, m))
      h2_SD <-try_or_NA(as.numeric(mxSE(top.A / top.Vtot, m)))
      e2 <- as.numeric(mxEval(top.E / top.Vtot, m))
      e2_SD <- try_or_NA(as.numeric(mxSE(top.E / top.Vtot, m)))

      tibble(Param = c("A", "E", "h2", "e2"),
             Value = c(A, E, h2, e2),
             SD = c(A_SD, E_SD, h2_SD, e2_SD),
             model = m$name)




    } else if (model == "CE") {


      # Variance components
      C <- as.numeric(mxEval(top.C, m))
      E <- as.numeric(mxEval(top.E, m))

      C_SD <- try_or_NA(as.numeric(mxSE(top.C, m)))
      E_SD <- try_or_NA(as.numeric(mxSE(top.E, m)))


      # Heritability
      c2 <- as.numeric(mxEval(top.C / top.Vtot, m))
      c2_SD <- try_or_NA(as.numeric(mxSE(top.C / top.Vtot, m)))
      e2 <- as.numeric(mxEval(top.E / top.Vtot, m))
      e2_SD <- try_or_NA(as.numeric(mxSE(top.E / top.Vtot, m)))

      tibble(Param = c("C", "E", "c2", "e2"),
             Value = c(C, E, c2, e2),
             SD = c(C_SD, E_SD, c2_SD, e2_SD),
             model = m$name)

    }





  }

  out_info <- tibble(Trait = x$trait,
                     Response_type = x$response_type,
                     Constrained = x$constrained)

  out <- x[!names(x) %in% c("response_type", "trait", "constrained")] %>%
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}


#' @export
get_estimates.ADE.uni <- function(x, ...) {


  extract_from_model <- function(m) {

    model <- m$name

    if (model == "ACE")

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
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}



#' @export
get_estimates.ACE.biv.chol <- function(x, ...) {


  extract_from_model <- function(m) {

    model <- m$name

    cov <- mxEval(top.Vtot, m)[2, 1]
    cov_SD <- try_or_NA(mxSE(top.Vtot, m)[2, 1])

    covA <- mxEval(top.A, m)[2, 1]
    covA_SD <- try_or_NA(mxSE(top.A, m)[2, 1])

    covC <- mxEval(top.C, m)[2, 1]
    covC_SD <- try_or_NA(mxSE(top.C, m)[2, 1])

    covE <- mxEval(top.E, m)[2, 1]
    covE_SD <- try_or_NA(mxSE(top.E, m)[2, 1])

    bivA <- mxEval(top.A / top.Vtot, m)[2, 1]
    bivA_SD <- try_or_NA(mxSE(top.A / top.Vtot, m)[2, 1])

    bivC <- mxEval(top.C / top.Vtot, m)[2, 1]
    bivC_SD <- try_or_NA(mxSE(top.C / top.Vtot, m)[2, 1])

    bivE <- mxEval(top.E / top.Vtot, m)[2, 1]
    bivE_SD <- try_or_NA(mxSE(top.E / top.Vtot, m)[2, 1])

    # Correlation decomposition
    rg <- mxEval(cov2cor(top.A), m)[2, 1]
    rg_SD <- try_or_NA(mxSE(cov2cor(top.A), m)[2, 1])

    rc <- mxEval(cov2cor(top.C), m)[2, 1]
    rc_SD <- try_or_NA(mxSE(cov2cor(top.C), m)[2, 1])

    re <- mxEval(cov2cor(top.E), m)[2, 1]
    re_SD <- try_or_NA(mxSE(cov2cor(top.E), m)[2, 1])

    # proportion unique additive genetic effect
    prop_unique_A_Y <- mxEval(a_r2c1 ^ 2 / (a_r2c1 ^ 2 + a_r2c2 ^ 2), m)
    prop_unique_A_Y_SD <- try_or_NA(as.numeric(mxSE(a_r2c1 ^ 2 / (a_r2c1 ^ 2 + a_r2c2 ^ 2), m)))

    # proportion unique shared environmental effect
    prop_unique_C_Y <- mxEval(c_r2c1 ^ 2 / (c_r2c1 ^ 2 + c_r2c2 ^ 2), m)
    prop_unique_C_Y_SD <- try_or_NA(as.numeric(mxSE(c_r2c1 ^ 2 / (c_r2c1 ^ 2 + c_r2c2 ^ 2), m)))

    # proportion unique random environmental effect
    prop_unique_E_Y <- mxEval(e_r2c1 ^ 2 / (e_r2c1 ^ 2 + e_r2c2 ^ 2), m)
    prop_unique_E_Y_SD <- try_or_NA(as.numeric(mxSE(e_r2c1 ^ 2 / (e_r2c1 ^ 2 + e_r2c2 ^ 2), m)))

    # proportion shared additive genetic effect
    prop_shared_A_Y <- mxEval(a_r2c2 ^ 2 / (a_r2c1 ^ 2 + a_r2c2 ^ 2), m)
    prop_shared_A_Y_SD <- try_or_NA(as.numeric(mxSE(a_r2c2 ^ 2 / (a_r2c1 ^ 2 + a_r2c2 ^ 2), m)))

    # proportion shared shared environmental effect
    prop_shared_C_Y <- mxEval(c_r2c2 ^ 2 / (c_r2c1 ^ 2 + c_r2c2 ^ 2), m)
    prop_shared_C_Y_SD <- try_or_NA(as.numeric(mxSE(c_r2c2 ^ 2 / (c_r2c1 ^ 2 + c_r2c2 ^ 2), m)))

    # proportion shared random environmental effect
    prop_shared_E_Y <- mxEval(e_r2c2 ^ 2 / (e_r2c1 ^ 2 + e_r2c2 ^ 2), m)
    prop_shared_E_Y_SD <- try_or_NA(as.numeric(mxSE(e_r2c2 ^ 2 / (e_r2c1 ^ 2 + e_r2c2 ^ 2), m)))


    # Variance components
    A1 <- mxEval(top.A, m)[1, 1]
    C1 <- mxEval(top.C, m)[1, 1]
    E1 <- mxEval(top.E, m)[1, 1]

    A1_SD <- try_or_NA(mxSE(top.A, m)[1, 1])
    C1_SD <- try_or_NA(mxSE(top.C, m)[1, 1])
    E1_SD <- try_or_NA(mxSE(top.E, m)[1, 1])

    A2 <- mxEval(top.A, m)[2, 2]
    C2 <- mxEval(top.C, m)[2, 2]
    E2 <- mxEval(top.E, m)[2, 2]

    A2_SD <- try_or_NA(mxSE(top.A, m)[2, 2])
    C2_SD <- try_or_NA(mxSE(top.C, m)[2, 2])
    E2_SD <- try_or_NA(mxSE(top.E, m)[2, 2])


    # Heritability
    h2_X <- mxEval(top.A / top.Vtot, m)[1, 1]
    h2_X_SD <- try_or_NA(mxSE(top.A / top.Vtot, m)[1, 1])
    h2_Y <- mxEval(top.A / top.Vtot, m)[2, 2]
    h2_Y_SD <- try_or_NA(mxSE(top.A / top.Vtot, m)[2, 2])

    tibble(Param = c("cov", "covA", "covC", "covE", "bivA", "bivC", "bivE", "rg", "rc", "re", "prop_unique_A_Y", "prop_unique_C_Y", "prop_unique_E_Y", "prop_shared_A_Y", "prop_shared_C_Y", "prop_shared_E_Y", "A1", "C1", "E1", "A2", "C2", "E2", "h2_X", "h2_Y"),
           Value = c(cov, covA, covC, covE, bivA, bivC, bivE, rg, rc, re, prop_unique_A_Y, prop_unique_C_Y, prop_unique_E_Y, prop_shared_A_Y, prop_shared_C_Y, prop_shared_E_Y, A1, C1, E1, A2, C2, E2, h2_X, h2_Y),
           SD = c(cov_SD, covA_SD, covC_SD, covE_SD, bivA_SD, bivC_SD, bivE_SD, rg_SD, rc_SD, re_SD, prop_unique_A_Y_SD, prop_unique_C_Y_SD, prop_unique_E_Y_SD, prop_shared_A_Y_SD, prop_shared_C_Y_SD, prop_shared_E_Y_SD, A1_SD, C1_SD, E1_SD, A2_SD, C2_SD, E2_SD, h2_X_SD, h2_Y_SD),
           model = m$name)

  }

  x[!names(x) %in% c("traitX", "traitY")] |>
    compact() |>
    map_dfr(extract_from_model) |>
    mutate(TraitX = x$TraitX, TraitY = x$TraitY)

}


#' @export
get_estimates.ACE.5group <- function(x, ...) {

  # Create Algebra for Variance Components
  # rowUS     <- rep('US', nv)
  # colUS     <- rep(c('VAf','VCf','VEf','SAf','SCf','SEf','VAm','VCm','VEm','SAm','SCm','SEm','rg','rc'), each=nv)
  # estUS     <- mxAlgebra(expression= cbind(VAf, VCf, VEf, VAf/Vf, VCf/Vf, VEf/Vf, VAm+VAms, VCm+VCms, VEm, (VAm+VAms)/Vm, (VCm+VCms)/Vm, VEm/Vm, rg, rc), name="US", dimnames=list(rowUS, colUS))

  extract_from_model <- function(m) {

    model <- m$name


    if (!model %in% c("ACE_5group", "CE_5group", "AE_5group")) {

      var_comp_A_fem <- as.numeric(mxEval(VAf11, m))
      var_comp_C_fem <- as.numeric(mxEval(VCf11, m))
      var_comp_E_fem <- as.numeric(mxEval(VEf11, m))

      var_comp_A_fem_SD <- try_or_NA(as.numeric(mxSE(VAf11, m)))
      var_comp_C_fem_SD <- try_or_NA(as.numeric(mxSE(VCf11, m)))
      var_comp_E_fem_SD <- try_or_NA(as.numeric(mxSE(VEf11, m)))

      var_comp_A_male <- as.numeric(mxEval(VAm11 + VAms11, m))
      var_comp_C_male <- as.numeric(mxEval(VCm11 + VCms11, m))
      var_comp_E_male <- as.numeric(mxEval(VEm11, m))

      var_comp_A_male_SD <- try_or_NA(as.numeric(mxSE(VAm11 + VAms11, m)))
      var_comp_C_male_SD <- try_or_NA(as.numeric(mxSE(VCm11 + VCms11, m)))
      var_comp_E_male_SD <- try_or_NA(as.numeric(mxSE(VEm11, m)))

      h2_fem <- as.numeric(mxEval(VAf11 / Vf, m))
      c2_fem <- as.numeric(mxEval(VCf11 / Vf, m))
      e2_fem <- as.numeric(mxEval(VEf11 / Vf, m))

      h2_fem_SD <- try_or_NA(as.numeric(mxSE(VAf11 / Vf, m)))
      c2_fem_SD <- try_or_NA(as.numeric(mxSE(VCf11 / Vf, m)))
      e2_fem_SD <- try_or_NA(as.numeric(mxSE(VEf11 / Vf, m)))

      h2_male <- as.numeric(mxEval((VAm11 + VAms11) / Vm, m))
      c2_male <- as.numeric(mxEval((VCm11 + VCms11) / Vm, m))
      e2_male <- as.numeric(mxEval(VEm11 / Vm, m))

      h2_male_SD <- try_or_NA(as.numeric(mxSE((VAm11 + VAms11) / Vm, m)))
      c2_male_SD <- try_or_NA(as.numeric(mxSE((VCm11 + VCms11) / Vm, m)))
      e2_male_SD <- try_or_NA(as.numeric(mxSE(VEm11 / Vm, m)))


      res <- tibble(Param = c("var_comp_A_fem", "var_comp_C_fem", "var_comp_E_fem", "var_comp_A_male", "var_comp_C_male", "var_comp_E_male", "h2_fem", "c2_fem", "e2_fem", "h2_male", "c2_male", "e2_male"),
                    Value = c(var_comp_A_fem, var_comp_C_fem, var_comp_E_fem, var_comp_A_male, var_comp_C_male, var_comp_E_male, h2_fem, c2_fem, e2_fem, h2_male, c2_male, e2_male),
                    SD = c(var_comp_A_fem_SD, var_comp_C_fem_SD, var_comp_E_fem_SD, var_comp_A_male_SD, var_comp_C_male_SD, var_comp_E_male_SD, h2_fem_SD, c2_fem_SD, e2_fem_SD, h2_male_SD, c2_male_SD, e2_male_SD),
                    model = m$name)

      if (model %in% c("ACE_5group_ra", "AE_5group_ra")) {

        rg_sex <- as.numeric(mxEval(rg, m))

        rg_sex_SD <- try_or_NA(as.numeric(mxSE(rg, m)))

        res1 <- tibble(Param = "rg_sex",
                       Value = rg_sex,
                       SD = rg_sex_SD,
                       model = model)

        res <- bind_rows(res, res1)


      }

      if (model %in% c("ACE_5group_rc", "CE_5group_rc")) {


        rc_sex <- as.numeric(mxEval(rc, m))
        rc_sex_SD <- try_or_NA(as.numeric(mxSE(rc, m)))

        res2 <- tibble(Param = "rc_sex",
                       Value = rc_sex,
                       SD = rc_sex_SD,
                       model = model)

        res <- bind_rows(res, res2)


      }

    } else  {

      # Variance components
      A <- as.numeric(mxEval(VA11, m))
      C <- as.numeric(mxEval(VC11, m))
      E <- as.numeric(mxEval(VE11, m))

      A_SD <- try_or_NA(as.numeric(mxSE(VA11, m)))
      C_SD <- try_or_NA(as.numeric(mxSE(VC11, m)))
      E_SD <- try_or_NA(as.numeric(mxSE(VE11, m)))

      # Heritability
      h2 <- as.numeric(mxEval(VA11 / (VA11 + VC11 + VE11), m))
      h2_SD <- try_or_NA(as.numeric(mxSE(VA11 / (VA11 + VC11 + VE11), m)))
      c2 <- as.numeric(mxEval(VC11 / (VA11 + VC11 + VE11), m))
      c2_SD <- try_or_NA(as.numeric(mxSE(VC11 / (VA11 + VC11 + VE11), m)))
      e2 <- as.numeric(mxEval(VE11 / (VA11 + VC11 + VE11), m))
      e2_SD <- try_or_NA(as.numeric(mxSE(VE11 / (VA11 + VC11 + VE11), m)))

      res <- tibble(Param = c("A", "C", "E", "h2", "c2", "e2"),
                    Value = c(A, C, E, h2, c2, e2),
                    SD = c(A_SD, C_SD, E_SD, h2_SD, c2_SD, e2_SD),
                    model = model)





    }


    res


  }


  out_info <- tibble(Trait = x$trait,
                     Response_type = x$response_type,
                     Constrained = x$constrained)

  out <- x[!names(x) %in% c("response_type", "trait", "constrained")] %>%
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}



# Get estimates from saturated model --------------------------------------
#' @export
get_estimates.saturated.binary <- function(x, ...) {




  extract_from_model <- function(m) {

    model <- m$name

    # Twin correlations
    rMZ <- as.numeric(mxEval(rMZ, m))
    rDZ <- as.numeric(mxEval(rDZ, m))


    rMZ_SD <- try_or_NA(as.numeric(mxSE(rMZ, m)))
    rDZ_SD <- try_or_NA(as.numeric(mxSE(rDZ, m)))

    res <- tibble(Param = c("rMZ", "rDZ"),
                  Value = c(rMZ, rDZ),
                  SD = c(rMZ_SD, rDZ_SD),
                  model = model)


  }





  out_info <- tibble(Trait = x$trait, Response_type = x$response_type)

  out <- x[!names(x) %in% c("response_type", "trait")] %>%
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out




}


#' @export
get_estimates.saturated.num <- function(x, ...) {




  extract_from_model <- function(m) {

    model <- m$name

    # Twin correlations
    rMZ <- as.numeric(mxEval(cMZ21, m))
    rDZ <- as.numeric(mxEval(cDZ21, m))


    rMZ_SD <- try_or_NA(as.numeric(mxSE(cMZ21, m)))
    rDZ_SD <- try_or_NA(as.numeric(mxSE(cDZ21, m)))

    res <- tibble(Param = c("rMZ", "rDZ"),
                  Value = c(rMZ, rDZ),
                  SD = c(rMZ_SD, rDZ_SD),
                  model = model)


  }





  out_info <- tibble(Trait = x$trait, Response_type = x$response_type)

  out <- x[!names(x) %in% c("response_type", "trait")] %>%
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}




#' @export
get_estimates.saturated.5group.num <- function(x, ...) {




  extract_from_model <- function(m) {

    model <- m$name

    # Twin correlations
    rMZf <- as.numeric(mxEval(cMZf21, m))
    rDZf <- as.numeric(mxEval(cDZf21, m))
    rMZm <- as.numeric(mxEval(cMZm21, m))
    rDZm <- as.numeric(mxEval(cDZm21, m))
    rDZo <- as.numeric(mxEval(cDZo21, m))

    rMZf_SD <- try_or_NA(as.numeric(mxSE(cMZf21, m)))
    rDZf_SD <- try_or_NA(as.numeric(mxSE(cDZf21, m)))
    rMZm_SD <- try_or_NA(as.numeric(mxSE(cMZm21, m)))
    rDZm_SD <- try_or_NA(as.numeric(mxSE(cDZm21, m)))
    rDZo_SD <- try_or_NA(as.numeric(mxSE(cDZo21, m)))

    tibble(Param = c("rMZf", "rDZf", "rMZm", "rDZm", "rDZo"),
           Value = c(rMZf, rDZf, rMZm, rDZm, rDZo),
           SD = c(rMZf_SD, rDZf_SD, rMZm_SD, rDZm_SD, rDZo_SD),
           model = model)


  }

  out_info <- tibble(Trait = x$trait, Response_type = x$response_type)

  out <- x[!names(x) %in% c("response_type", "trait")] %>%
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}


#' @export
get_estimates.saturated.5group.binary <- function(x, ...) {


  extract_from_model <- function(m) {

    model <- m$name

    # Twin correlations
    rMZf <- as.numeric(mxEval(r_mzf_21, m))
    rDZf <- as.numeric(mxEval(r_dzf_21, m))
    rMZm <- as.numeric(mxEval(r_mzm_21, m))
    rDZm <- as.numeric(mxEval(r_dzm_21, m))
    rDZo <- as.numeric(mxEval(r_dzo_21, m))

    rMZf_SD <- try_or_NA(as.numeric(mxSE(r_mzf_21, m)))
    rDZf_SD <- try_or_NA(as.numeric(mxSE(r_dzf_21, m)))
    rMZm_SD <- try_or_NA(as.numeric(mxSE(r_mzm_21, m)))
    rDZm_SD <- try_or_NA(as.numeric(mxSE(r_dzm_21, m)))
    rDZo_SD <- try_or_NA(as.numeric(mxSE(r_dzo_21, m)))

    tibble(Param = c("rMZf", "rDZf", "rMZm", "rDZm", "rDZo"),
           Value = c(rMZf, rDZf, rMZm, rDZm, rDZo),
           SD = c(rMZf_SD, rDZf_SD, rMZm_SD, rDZm_SD, rDZo_SD),
           model = model)


  }

  out_info <- tibble(Trait = x$trait, Response_type = x$response_type)

  out <- x[!names(x) %in% c("response_type", "trait")] |>
    compact() |>
    map_dfr(extract_from_model)

  out <- bind_cols(out, out_info)
  out


}

