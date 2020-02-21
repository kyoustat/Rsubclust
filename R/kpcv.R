#' Projective k-Means Clustering with Varying Dimension
#' 
#' 
#' 
#' @examples 
#' ## generate a toy example of two line components
#' set.seed(19)
#' tester = simLP(n=100, nl=2, np=0)
#' data   = tester$data
#' label  = tester$class
#' 
#' ## do PCA for data reduction
#' proj = base::eigen(stats::cov(data))$vectors[,1:2]
#' dat2 = data%*%proj
#' 
#' ## run k-plane algorithm with K=2, 3, and 4
#' kpcv2 = kpcv(data, K=2)
#' kpcv3 = kpcv(data, K=3)
#' kpcv4 = kpcv(data, K=4)
#' 
#' ## extract clustering
#' finc2 = kpcv2$cluster
#' finc3 = kpcv3$cluster
#' finc4 = kpcv4$cluster
#' 
#' ## visualize
#' opar <- par(mfrow=c(3,4), no.readonly=TRUE)
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc2+1,main="K=2:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc2+1,main="K=2:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc3+1,main="K=3:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc3+1,main="K=3:Axis(2,3)")
#' 
#' plot(dat2[,1],dat2[,2],pch=19,cex=0.3,col=finc4+1,main="K=4:PCA")
#' plot(data[,1],data[,2],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(1,2)")
#' plot(data[,1],data[,3],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(1,3)")
#' plot(data[,2],data[,3],pch=19,cex=0.3,col=finc4+1,main="K=4:Axis(2,3)")
#' par(opar)
#' 
#' @references 
#' \insertRef{agarwal_k-means_2004}{Rsubclust}
#' 
#' @export
kpcv <- function(X, K=2, iter=496, init=c("kmeans","random"), print.progress=TRUE){
  #########################################################
  # Initialization : use the notation from the paper
  n = nrow(X)
  d = ncol(X)
  K = round(K)
  
  init = match.arg(init)
  if (all(init=="kmeans")){
    old.clust = round(stats::kmeans(X, center=K)$cluster)
  } else if (all(init=="random")){
    old.clust = base::sample(1:K, base::sample(1:K, n-K, replace=TRUE))
  }
  old.label = list()
  for (k in 1:K){
    old.label[[k]] = which(old.clust==k)
  }
  
  dims = sample(1:(d-1), K, replace = TRUE)
  run.old = kpc.compute.qflat(X, old.label, dims)
  
  #########################################################
  # Iterate
  for (it in 1:iter){
    # (a) computing optimal q-flats
    run.new = kpc.compute.qflat(X, old.label, dims)
    if ((length(run.new)==1)&&(run.new==FALSE)){
      if (print.progress){
        print(paste('* kpcv : iteration ', it,'/',iter,' terminated due to assignment collapse at ', Sys.time(),sep=""))  
      }
      break
    }
    dims  = run.new$dfixed
    
    # (b) re-assign points to nearest flat
    new.clust = kpc.assign.qflat(X, run.new)
    if ((length(new.clust)==1)&&(new.clust==FALSE)){
      if (print.progress){
        print(paste('* kpcv : iteration ', it,'/',iter,' terminated due to assignment collapse at ', Sys.time(),sep=""))  
      }
      break
    }
    new.label = list()
    for (k in 1:K){
      new.label[[k]] = which(new.clust==k)
    }
    
    # (c) Recompute dimension
    for (k in 1:K){
      xpart = X[new.label[[k]],]
      if (is.vector(xpart)){
        dims[k] = 1
      } else {
        dims[k] = max(min(round(Matrix::rankMatrix(xpart)[1]), d-1), 1)
      }
    }
    
    # (c) updating
    inc.clust = sum((old.clust-new.clust)^2)
    old.clust = new.clust
    old.label = list()
    run.old   = run.new
    for (k in 1:K){
      old.label[[k]] = which(old.clust==k)
    }
    
    if (it >= 5){
      if (inc.clust < 5){
        if (print.progress){
          print(paste('* kpcv : iteration ', it,'/',iter,' terminated at ', Sys.time(),sep=""))
        }
        break
      }
    }
    if (print.progress){
      if(it%%10==0){
        print(paste('* kpcv : iteration ', it,'/',iter,' complete at ', Sys.time(),sep=""))
      }
    }
  }
  
  #########################################################
  # Return output
  output = list()
  output$cluster = old.clust
  output$center  = run.old$centers
  output$basis   = run.old$basis
  output$dfixed  = dims
  output$name    = "Rsubclust:kpcv"
  return(output)
}
