

#' @import dplyr
#' @export
prep_uni_data <- function(db, db_cohort, trait, response_type, same_sex = TRUE) {

  db <- db %>%
    select(pairnnr, twinnr, all_of(trait), Female, Zyg, b_year) %>%
    rename(X = all_of(trait))

  db_cohort <- db_cohort %>%
    select(pairnnr, twinnr,  Zyg, Female, b_year)

  singletons <- db %>%
    group_by(pairnnr) %>%
    mutate(N = n()) %>%
    ungroup() %>%
    filter(N < 2)

  db_look_for_pairs <- db_cohort %>%
    filter(pairnnr %in% singletons$pairnnr) %>%
    group_by(pairnnr) %>%
    mutate(N = n()) %>%
    ungroup(pairnnr)

  db_found_other_twin <- db_look_for_pairs %>%
    filter(N == 2) %>%
    filter(!twinnr %in% singletons$twinnr)

  db_no_other_twin <- db_look_for_pairs %>%
    filter(N == 1)

  db_construct_other_twin <- tibble(pairnnr = db_no_other_twin$pairnnr,
                                    twinnr = paste0(pairnnr, ifelse(str_sub(db_no_other_twin$twinnr, start = -1) == 1, 2, 1)),
                                    Zyg = db_no_other_twin$Zyg,
                                    X = NA,
                                    b_year = db_no_other_twin$b_year,
                                    Female = case_when(
                                      Zyg == "MZ" ~ db_no_other_twin$Female,
                                      Zyg == "DZ_diff_sex" ~ 1 - db_no_other_twin$Female,
                                      Zyg == "DZ_same_sex" ~ db_no_other_twin$Female))


  db_complete <- bind_rows(db, db_found_other_twin, db_construct_other_twin) %>%
    arrange(pairnnr, twinnr) %>%
    mutate(Twin = str_sub(twinnr, start = -1)) %>%
    as_tibble() %>%
    select(-N)

  birth_year_poly <- db_complete %>%
    pull(b_year) %>%
    scale() %>%
    poly(degree = 2, raw = TRUE) %>%
    as_tibble() %>%
    rename(Birth_year_first = `1`, Birth_year_second = `2`)

  db_complete <- bind_cols(db_complete, birth_year_poly)

  datMZ <- db_complete %>%
    filter(Zyg == "MZ") %>%
    select(-Zyg, -twinnr) %>%
    pivot_wider(names_from = Twin, values_from = c(X, Female, Birth_year_first, Birth_year_second), names_sep = "") %>%
    select(pairnnr, X1, X2, Female1, Female2, contains("Birth_year")) %>%
    as_tibble()

  datDZ <- db_complete %>%
    {
      if (same_sex) filter(., Zyg == "DZ_same_sex") else filter(., Zyg %in% c("DZ_same_sex", "DZ_diff_sex"))

    } %>%
    group_by(pairnnr) %>%
    select(-Zyg, -twinnr) %>%
    pivot_wider(names_from = Twin, values_from = c(X, Female, Birth_year_first, Birth_year_second), names_sep = "") %>%
    select(pairnnr, X1, X2, Female1, Female2, contains("Birth_year")) %>%
    as_tibble()

  # Get rid of pairs where both are NAs
  exclude_MZ <- datMZ %>%
    rowwise() %>%
    mutate(Both_NA = is.na(X1) & is.na(X2)) %>%
    filter(Both_NA) %>%
    pull(pairnnr)


  exclude_DZ <- datDZ %>%
    rowwise() %>%
    mutate(Both_NA = is.na(X1) & is.na(X2)) %>%
    filter(Both_NA) %>%
    pull(pairnnr)

  datMZ <- datMZ %>% filter(!pairnnr %in% exclude_MZ)
  datDZ <- datDZ %>% filter(!pairnnr %in% exclude_DZ)

  if (response_type == "binary") {

    datMZ <- datMZ %>%
      mutate(X1 = mxFactor(X1, levels = sort(unique(X1))),
             X2 = mxFactor(X2, levels = sort(unique(X2))))

    datDZ <- datDZ %>%
      mutate(X1 = mxFactor(X1, levels = sort(unique(X1))),
             X2 = mxFactor(X2, levels = sort(unique(X2))))

    class <- c("prep.uni.bin", "prep.uni")

  } else class <- c("prep.uni.num", "prep.uni")


  out <- list(DZ = datDZ,
              MZ = datMZ,
              same_sex = same_sex,
              trait = trait,
              response_type = response_type)

  class(out) <- class
  out

}


#' @import dplyr
#' @export
prep_uni_data_non_expand <- function(db, trait, response_type, covs = NULL, same_sex = TRUE) {

  if (!is_null(covs)) covs_in_twin_frame <- paste0(covs, c(1, 2)) else covs_in_twin_frame <- NULL

  db <- db %>%
    select(pairnnr, twinnr, all_of(trait), Female, Zyg, b_year, all_of(covs)) %>%
    rename(X = all_of(trait)) %>%
    group_by(pairnnr) %>%
    mutate(n_twins = n()) %>%
    ungroup() %>%
    filter(n_twins > 1) %>%
    select(-n_twins) %>%
    mutate(Twin = str_sub(twinnr, start = -1))

  birth_year_poly <- db %>%
    pull(b_year) %>%
    scale(scale = FALSE) %>%
    poly(degree = 2, raw = TRUE) %>%
    as_tibble() %>%
    rename(Birth_year_first = `1`, Birth_year_second = `2`)

  db <- bind_cols(db, birth_year_poly)


  datMZ <- db %>%
    filter(Zyg == "MZ") %>%
    mutate(Twin = str_sub(twinnr, start = -1)) %>%
    select(-Zyg, -twinnr) %>%
    pivot_wider(names_from = Twin, values_from = c(X, Female, Birth_year_first, Birth_year_second, all_of(covs)), names_sep = "") %>%
    select(pairnnr, X1, X2, Female1, Female2, contains("Birth_year"), all_of(covs_in_twin_frame)) %>%
    filter(if_all(all_of(covs_in_twin_frame), ~ !is.na(.))) |>
    as_tibble()

  datDZ <- db %>%
    {
      if (same_sex) filter(., Zyg == "DZ_same_sex") else filter(., Zyg %in% c("DZ_same_sex", "DZ_diff_sex"))

    } %>%
    mutate(Twin = str_sub(twinnr, start = -1)) %>%
    select(-Zyg, -twinnr) %>%
    pivot_wider(names_from = Twin, values_from = c(X, Female, Birth_year_first, Birth_year_second, all_of(covs)), names_sep = "") %>%
    select(pairnnr, X1, X2, Female1, Female2, contains("Birth_year"), all_of(covs_in_twin_frame)) %>%
    filter(if_all(all_of(covs_in_twin_frame), ~ !is.na(.))) |>
    as_tibble()


  if (response_type == "binary") {

    datMZ <- datMZ %>%
      mutate(X1 = mxFactor(X1, levels = sort(unique(X1))),
             X2 = mxFactor(X2, levels = sort(unique(X2))))

    datDZ <- datDZ %>%
      mutate(X1 = mxFactor(X1, levels = sort(unique(X1))),
             X2 = mxFactor(X2, levels = sort(unique(X2))))

    class <- c("prep.uni.bin", "prep.uni")

  } else class <- c("prep.uni.num", "prep.uni")


  out <- list(DZ = datDZ, MZ = datMZ, same_sex = same_sex, trait = trait, response_type = response_type, covs = covs)
  class(out) <- class
  out

}


#' @import dplyr
#' @export
prep_5groups <- function(prep) {


  if (prep$same_sex) {stop("Doesn't make sense for same_sex = TRUE")}
  if (!is.null(prep$covs)) covs <- paste0(prep$covs, c(1, 2)) else covs <- NULL

  mzmData <- prep$MZ %>%
    filter(Female1 == 0) %>%
    select(pairnnr, X1, X2, contains("Birth_year"), all_of(covs))

  dzmData <- prep$DZ %>%
    filter(Female1 == 0, Female2 == 0) %>%
    select(pairnnr, X1, X2, contains("Birth_year"), all_of(covs))

  mzfData <- prep$MZ %>%
    filter(Female1 == 1) %>%
    select(pairnnr, X1, X2, contains("Birth_year"), all_of(covs))

  dzfData <- prep$DZ %>%
    filter(Female1 == 1, Female2 == 1) %>%
    select(pairnnr, X1, X2, contains("Birth_year"), all_of(covs))

  dzoData <- prep$DZ %>%
    filter(Female1 != Female2) %>%
    mutate(Xm = if_else(Female1 == 0, X1, X2), Xf = if_else(Female1 == 1, X1, X2)) %>%
    select(-X1, -X2) %>%
    rename(X1 = Xf, X2 = Xm) %>%
    select(pairnnr, X1, X2, contains("Birth_year"), all_of(covs))

  if (prep$response_type == "binary") {

    dzoData <- dzoData %>%
      mutate(X1 = mxFactor(X1, levels(X1)), X2 = mxFactor(X2, levels(X2)))

    class <- "prep.uni.5group.binary"

  } else {

    class <- "prep.uni.5group.num"

  }


  out <- list(mzm = mzmData,
              dzm = dzmData,
              mzf = mzfData,
              dzf = dzfData,
              dzo = dzoData,
              same_sex = FALSE,
              trait = prep$trait,
              response_type = prep$response_type,
              covs = prep$covs)


  class(out) <- class
  out

}

#' @export
prep_bivariate_data_non_expand <- function(db, traitX, traitY, covs = NULL, response_typeX, response_typeY, same_sex = TRUE) {


  if (!is_null(covs)) covs_in_twin_frame <- paste0(covs, c(1, 2)) else covs_in_twin_frame <- NULL

  X_is_factor <- is.factor(db[[traitX]])
  Y_is_factor <- is.factor(db[[traitY]])

  db <- db %>%
    select(pairnnr, twinnr, X = all_of(traitX), Y = all_of(traitY), Zyg, Female, b_year, all_of(covs))

  birth_year_poly <- db %>%
    pull(b_year) %>%
    scale() %>%
    poly(degree = 2, raw = TRUE) %>%
    as_tibble() %>%
    rename(Birth_year_first = `1`, Birth_year_second = `2`)

  db <- bind_cols(db, birth_year_poly)

  datMZ <- db %>%
    filter(Zyg == "MZ") %>%
    mutate(Twin = str_sub(twinnr, start = -1)) %>%
    select(-Zyg, -twinnr) %>%
    pivot_wider(id_cols = pairnnr, names_from = Twin, values_from = c(X, Y, Female, Birth_year_first, Birth_year_second, all_of(covs)), names_sep = "") %>%
    select(pairnnr, X1, X2, Y1, Y2, Female1, Female2, contains("Birth_year"), all_of(covs_in_twin_frame)) %>%
    as_tibble()

  datDZ <- db %>%
    {

      if (same_sex) filter(., Zyg == "DZ_same_sex") else filter(., Zyg %in% c("DZ_same_sex", "DZ_diff_sex"))

    } %>%
    mutate(Twin = str_sub(twinnr, start = -1)) %>%
    select(-Zyg, -twinnr) %>%
    pivot_wider(id_cols = pairnnr, names_from = Twin, values_from = c(X, Y, Female, Birth_year_first, Birth_year_second, all_of(covs)), names_sep = "") %>%
    select(pairnnr, X1, X2, Y1, Y2, Female1, Female2, contains("Birth_year"), all_of(covs_in_twin_frame)) %>%
    as_tibble()


  if (response_typeX == "binary") {

    datMZ <- datMZ %>%
      mutate(X1 = mxFactor(X1, levels = sort(unique(X1))),
             X2 = mxFactor(X2, levels = sort(unique(X2))))

    datDZ <- datDZ %>%
      mutate(X1 = mxFactor(X1, levels = sort(unique(X1))),
             X2 = mxFactor(X2, levels = sort(unique(X2))))


  }

  if (response_typeY == "binary") {

    datMZ <- datMZ %>%
      mutate(Y1 = mxFactor(Y1, levels = sort(unique(Y1))),
             Y2 = mxFactor(Y2, levels = sort(unique(Y2))))

    datDZ <- datDZ %>%
      mutate(Y1 = mxFactor(Y1, levels = sort(unique(Y1))),
             Y2 = mxFactor(Y2, levels = sort(unique(Y2))))


  }

  # Get rid of pairs where both are NAs
  these <- select(datMZ, X1, X2, Y1, Y2) %>%
    is.na() %>%
    rowSums()

  exclude_MZ <- datMZ$pairnnr[these == 4]

  these <- select(datDZ, X1, X2, Y1, Y2) %>%
    is.na() %>%
    rowSums()

  exclude_DZ <- datDZ$pairnnr[these == 4]

  datMZ <- datMZ %>% filter(!pairnnr %in% exclude_MZ)
  datDZ <- datDZ %>% filter(!pairnnr %in% exclude_DZ)


  out <- list(DZ = datDZ, MZ = datMZ, traitX = traitX, traitY = traitY, response_typeX = response_typeX, response_typeY = response_typeY)
  class(out) <- "prep.biv"
  out

}

