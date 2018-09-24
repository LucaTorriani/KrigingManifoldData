.onAttach <- function(libname, pkgname)
  packageStartupMessage("KrigingManifoldData 0.1.23 loaded\nCopyright Noi 2018")

.onUnload <- function(libpath)
  library.dynam.unload("KrigingManifoldData",  libpath)
