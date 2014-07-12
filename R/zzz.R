#If you have compiled code, add a .First.lib() function in the 'R' subdirectory to load it.
#.First.lib <- function(lib, pkg) library.dynam("modeest", pkg, lib)

.noGenerics <- TRUE

.onLoad <-
function(libname,
         pkgname)
{
  #! J'ai l'impression que dans les good practices,  packageStartupMessage doit plutôt être appelée dans .onAttach...
  packageStartupMessage("\nThis is package 'modeest' written by P. PONCET.\nFor a complete list of functions, use 'library(help = \"modeest\")' or 'help.start()'.\n")
}

.onUnload <-
function(libpath)
{
  library.dynam.unload("modeest", libpath)
}

