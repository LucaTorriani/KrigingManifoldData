.onAttach <- function(libname, pkgname)
  packageStartupMessage("Manfidolgstat 0.1.23 loaded\n")

.onUnload <- function(libpath){
  library.dynam.unload("Manfidolgstat",  libpath)
}
