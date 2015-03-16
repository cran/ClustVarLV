#' Print the CLV_kmeans results
#' 
#' Print the CLV_kmeans results
#' 
#' @param x an object of class clvkmeans
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso CLV_kmeans
#' @export print clvkmeans
print.clvkmeans =  function (x, ...) 
{
  if (!inherits(x, "clvkmeans")) 
    stop("non convenient object")
  p = x$param$p
  n = x$param$n
  EXTu=x$param$EXTu
  EXTr=x$param$EXTr
  method =  x$param$method
  strategy = x$param$strategy
  
  cat("\n")
  cat(paste("number of variables: ", p), sep = " ")
  cat("\n")
  cat(paste("number of observations: ", n), sep = " ")
  cat("\n")
  if (method==1) cat("measure of proximity: squared covariance")
  if (method==2) cat("measure of proximity: covariance")
  cat("\n")
  cat(paste("number of clusters: ", x$param$K), sep = " ")
  cat("\n")
  if (strategy=="sparselv") cat(paste("number of noise: ", length(which(unlist(x$sloading)==0))), sep = " ")
  if (strategy=="kplusone") cat(paste("number of noise: ", length(which(x$clusters[2,]==0))), sep = " ")
  cat("\n")
  cat("\n")
  cat("$tabres: clustering criterion")
  cat("\n")
  cat("$clusters: groups membership")
  cat("\n")
  cat("$comp: latent components of the clusters")
  if ((EXTu==1)|(EXTr==1)){
    cat("\n")
    cat("$loading: loadings of the external variables")
  }
  if (strategy=="sparselv"){
    cat("\n")
    cat("$sloading: sparse loading for latent variable")
  } 
  cat("\n")
  
}