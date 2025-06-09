

.onLoad <- function(libname, pkgname) {
  # something to run

  OpenMx::mxOption(model = NULL, "Default optimizer", "NPSOL")

}
