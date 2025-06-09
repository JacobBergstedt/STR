
#' @export
fit_cor <- function(x, ...) {

  UseMethod("fit_cor")

}

#' @export
fit_cor.prep.biv <- function(x) {


  fit_polychor <- function(X, Y) {

    obj <- polycor::polychor(X, Y, ML = TRUE, std.err = TRUE)
    tibble(rho = obj$rho, se = sqrt(diag(obj$var)[1]))

  }

  fit_cor_for_zyg <- function(xz) {

    spec <- tibble(X = c("X1", "X1", "X1", "Y1", "Y1", "X2"),
                   Y = c("Y1", "X2", "Y2", "X2", "Y2", "Y2"),
                   matrix_positions = c("c1r2", "c1r3", "c1r4", "c2r3", "c2r4", "c3r4"),
                   Labels = c("Correlation", "Within-trait X", "Cross-trait", "Cross-trait", "Within-trait Y", "correlation"))


    map2_dfr(spec$X, spec$Y, ~ fit_polychor(xz[[.x]], xz[[.y]])) |>
      bind_cols(spec)


  }

  out_MZ <- fit_cor_for_zyg(x$MZ) |> add_column(Zyg = "MZ")
  out_DZ <- fit_cor_for_zyg(x$DZ) |> add_column(Zyg = "DZ")

  bind_rows(out_DZ, out_MZ)



}
