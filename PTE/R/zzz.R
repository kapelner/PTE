.onLoad = function(libname, pkgname){
	if (is.null(getOption("mc.cores"))){
		options(mc.cores = max(parallel::detectCores() - 1, 1))
	}
}

.onAttach = function(libname, pkgname){
  packageStartupMessage("Welcome to PTE v1.7 by Adam Kapelner, Alina Levine & Justin Bleich\n")
}