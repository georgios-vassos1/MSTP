# Internal environment for Julia session state
.mstp <- new.env(parent = emptyenv())
.mstp$loaded <- FALSE

.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}
