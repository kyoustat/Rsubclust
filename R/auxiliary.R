## Auxiliary Functions for 'Rsubclust' Package
#  (1) rsc.d2subspace



# (1) rsc.d2subspace ------------------------------------------------------
#' @export
rsc.d2subspace <- function(X, basis, center){
  center = as.vector(center)
  if (is.vector(basis)){
    basis = matrix(basis, ncol=1)
  } else {
    basis = as.matrix(basis)
  }
  if (nrow(basis)!=ncol(X)){
    stop("* basis is invalid.")
  }
  
  projection = basis%*%t(basis)
  return(as.vector(rsc_d2subspace(X, projection, center)))
}
# library(Rdimtools)
# dat = as.matrix(scale(aux.gensamples(), center = TRUE, scale = FALSE))
# bb = rnorm(3)
# bb = bb/sqrt(sum(bb^2))
# 
# dd = rsc.d2subspace(dat, bb, rep(0,3))
