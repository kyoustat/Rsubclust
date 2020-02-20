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
  allnames = c("kplane", "kpc")
  allnames = paste0("Rsubclust:",allnames)
  if (!(foutput$name %in% allnames)){
    stop("* scpredict : 'foutput' should be a list output from functions in this package.")
  }
  
  #################################################
  # Switching argument
  output = switch(foutput$name,
                  "Rsubclust:kplane" = predict_kplane(X, foutput))
}


# predict functions -------------------------------------------------------
# (1) predict_kplane



# (1) predict_kplane ------------------------------------------------------
#' @keywords internal
predict_kplane <- function(X, kpoutput){
  ## simply apply cluster assignment
  K = ncol(kpoutput$w)
  pred.assign = kplane.assign.naive(X, K, kpoutput$w, kpoutput$gamma)
  return(pred.assign$cluster)
}