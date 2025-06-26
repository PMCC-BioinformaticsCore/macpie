.onLoad <- function(libname, pkgname) {
  if (is.null(getOption("java.parameters"))) {
    options(java.parameters = "-Xmx3g")
  }
}