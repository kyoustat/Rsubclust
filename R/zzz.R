.pkgenv <- new.env(parent = emptyenv())

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version = utils::packageVersion("Rsubclust")
  
  ## Print on Screen
  packageStartupMessage("** -------------------------------------------------------- **")
  packageStartupMessage("**  Subspace Clustering Algorithms with 'Rsubclust' package  ")
  packageStartupMessage("**")
  packageStartupMessage("** Version    : ",this.version,"      (",this.year,")",sep="")
  packageStartupMessage("** Maintainer : Kisung You (kyoustat@gmail.com)")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
  packageStartupMessage("** -------------------------------------------------------- **")
}

.onUnload <- function(libpath) {
  library.dynam.unload("Rsubclust", libpath)
}
