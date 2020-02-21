#' Prediction for Class Labels 
#' 
#' 
#' @export
scpredict <- function(X, foutput){
  #################################################
  # Needs to check 
  if ((!is.list(foutput))||("name"%in%names(foutput))){
    stop("* scpredict : 'foutput' should be a list output from functions in this package.")
  }
  allnames = c("kplane", "kpcf","kpcv")
  allnames = paste0("Rsubclust:",allnames)
  if (!(foutput$name %in% allnames)){
    stop("* scpredict : 'foutput' should be a list output from functions in this package.")
  }
  
  #################################################
  # Switching argument
  output = switch(foutput$name,
                  "Rsubclust:kplane" = predict_kplane(X, foutput),
                  "Rsubclust:kpcf"   = predict_kpc(X, foutput),
                  "Rsubclust:kpcv"   = predict_kpc(X, foutput))
  
  #################################################
  # Return 
  return(output)
}


# predict functions -------------------------------------------------------
# (1) predict_kplane
# (2) predict_kpc



# (2) predict_kpc ---------------------------------------------------------
#' @keywords internal
predict_kpc <- function(X, koutput){
  N = nrow(X)
  K = length(koutput$cluster)
  
  distmat = array(0,c(N,K))
  for (k in 1:K){
    k.basis  = koutput$basis[[k]]
    k.center = koutput$center[[k]]
    distmat[,k] = rsc.d2subspace(X, k.basis, k.center)
  }
  return(as.vector(apply(distmat, 1, which.min)))
}
# (1) predict_kplane ------------------------------------------------------
#' @keywords internal
predict_kplane <- function(X, kpoutput){
  ## simply apply cluster assignment
  K = ncol(kpoutput$w)
  pred.assign = kplane.assign.naive(X, K, kpoutput$w, kpoutput$gamma)
  return(pred.assign$cluster)
}