.onAttach <- function(libname, pkgname)
  packageStartupMessage("Manifoldgstat 0.1.23 loaded\n")

.onUnload <- function(libpath){
  library.dynam.unload("Manifoldgstat",  libpath)
}
