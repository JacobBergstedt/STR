


#' @export
transform_saturated_biv_output <- function(x) {

  spec <- tibble(X = c("X1", "X1", "X1", "Y1", "Y1", "X2"),
                 Y = c("Y1", "X2", "Y2", "X2", "Y2", "Y2"),
                 matrix_positions = c("c1r2", "c1r3", "c1r4", "c2r3", "c2r4", "c3r4"),
                 Labels = c("Correlation", "Within-trait X", "Cross-trait", "Cross-trait", "Within-trait Y", "correlation"))

  res <- summary(x$Saturated)$parameters |> select(name, rho = Estimate, se = Std.Error)


  out_MZ <- left_join(tibble(name = c("corMZ21", "corMZ31", "corMZ41", "corMZ32", "corMZ42", "corMZ43")), res) |>
    select(-name) |>
    bind_cols(spec) |>
    add_column(Zyg = "MZ")


  out_DZ <- left_join(tibble(name = c("corDZ21", "corDZ31", "corDZ41", "corDZ32", "corDZ42", "corDZ43")), res) |>
    select(-name) |>
    bind_cols(spec) |>
    add_column(Zyg = "DZ")

  bind_rows(out_DZ, out_MZ)

}

