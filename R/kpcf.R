#' Projective k-Means Clustering with Fixed Dimension
#' 
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
#' kpcf2 = kpcf(data, K=2, dims=c(1,1))
#' kpcf3 = kpcf(data, K=3, dims=c(1,1,1))
#' kpcf4 = kpcf(data, K=4, dims=c(1,1,1,1))
#' 
#' ## extract clustering
#' finc2 = kpcf2$cluster
#' finc3 = kpcf3$cluster
#' finc4 = kpcf4$cluster
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
kpcf<- function(X, K=2, dims=rep(1,K), iter=496, init=c("kmeans","random"), print.progress=TRUE){
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
  
  dims = round(dims)
  if ((length(dims)!=K)||(any(dims<1))||(any(dims>=d))){
    stop("* kpcf : 'dims' should be a vector of integer values in [1,ncol(X)).")
  }
  run.old = kpc.compute.qflat(X, old.label, dims)
  dims = run.old$dfixed
  
  #########################################################
  # Iterate
  for (it in 1:iter){
    # (a) computing optimal q-flats
    run.new = kpc.compute.qflat(X, old.label, dims)
    if ((length(run.new)==1)&&(run.new==FALSE)){
      if (print.progress){
        print(paste('* kpcf : iteration ', it,'/',iter,' terminated due to assignment collapse at ', Sys.time(),sep=""))  
      }
      break
    }
    dims  = run.new$dfixed

    # (b) re-assign points to nearest flat
    new.clust = kpc.assign.qflat(X, run.new)
    if ((length(new.clust)==1)&&(new.clust==FALSE)){
      if (print.progress){
        print(paste('* kpcf : iteration ', it,'/',iter,' terminated due to assignment collapse at ', Sys.time(),sep=""))  
      }
      break
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
          print(paste('* kpcf : iteration ', it,'/',iter,' terminated at ', Sys.time(),sep=""))
        }
        break
      }
    }
    if (print.progress){
      if(it%%10==0){
        print(paste('* kpcf : iteration ', it,'/',iter,' complete at ', Sys.time(),sep=""))
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
  output$name    = "Rsubclust:kpcf"
  return(output)
}

# auxiliary functions -----------------------------------------------------
#' @keywords internal
kpc.compute.qflat <- function(X, old.label, dims){
  K = length(dims)
  
  centers = list()
  baseonb = list()
  updated = rep(0,K)
  for (k in 1:K){
    xpart = X[old.label[[k]],]
    if (is.vector(xpart)){
      return(FALSE)
    }
    
    xcov = stats::cov(xpart)
    if (dims[k] >= nrow(xcov)){
      updated[k] = nrow(xcov)-1
      print(paste0("* kpcf : the number of observations in cluster ",k," is 
            smaller than the fixed dimension. Reduce the provided dimension."))
    } else {
      updated[k] = dims[k]
    }
    xeig = base::eigen(xcov)$vectors[,1:updated[k]]
    if (is.vector(xeig)){
      xeig = matrix(xeig, ncol = 1)
    }
    
    centers[[k]] = base::colMeans(xpart)
    baseonb[[k]] = xeig
  }
  
  output = list()
  output$center = centers
  output$basis  = baseonb
  output$dfixed = updated
  return(output)
}
#' @keywords internal
kpc.assign.qflat <- function(X, out.qflat){
  N = nrow(X)
  d = ncol(X)
  K = length(out.qflat$center)
  
  output = array(0,c(N,K))
  for (k in 1:K){
    k.basis  = out.qflat$basis[[k]]
    k.center = out.qflat$center[[k]]
    output[,k] = rsc.d2subspace(X, k.basis, k.center)
  }
  # for (k in 1:K){
  #   projU = base::diag(d)-(out.qflat$basis[[k]]%*%t(out.qflat$basis[[k]]))
  #   for (n in 1:N){
  #     xdiff = as.vector(projU%*%(as.vector(X[n,])-as.vector(out.qflat$center[[k]])))
  #     output[n,k] = sum(xdiff^2)
  #   }
  # }
  
  label = rep(0,N)
  for (n in 1:N){
    ntgt = as.vector(output[n,])
    label[n] = base::sample(which(ntgt <= min(ntgt)), 1)
  }
  if (length(unique(label))<K){
    return(FALSE)
  } else {
    return(label)
  }
}