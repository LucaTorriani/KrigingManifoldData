.onAttach <- function(libname, pkgname)
  packageStartupMessage("Manifoldgstat 1.0.0 loaded\n")

.onUnload <- function(libpath){
  library.dynam.unload("Manifoldgstat",  libpath)
}
