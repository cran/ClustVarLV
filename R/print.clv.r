#' Print the CLV results
#' 
#' Print the CLV results
#' 
#' @param x an object of class clv
#' @param ... further arguments passed to or from other methods
#' 
#' @seealso CLV
#' @S3method print clv
print.clv =  function (x, ...) 
{
  if (!inherits(x, "clv")) 
    stop("non convenient object")
  p = x$param$p
  n = x$param$n
  EXTu=x$param$EXTu
  EXTr=x$param$EXTr
  method =  x$param$method
 
  cat("\n")
  cat(paste("number of variables: ", p), sep = " ")
  cat("\n")
  cat(paste("number of observations: ", n), sep = " ")
  cat("\n")
  if (method==1) cat("measure of proximity: squared covariance")
  if (method==2) cat("measure of proximity: covariance")
  cat("\n")
  cat(paste("consolidation for K in c(",x$param$nmax,":2)",sep = ""))
  cat("\n")
  cat("\n")
  cat("$tarbre: results of the clustering")
  cat("\n")
  cat("$partitionK or [[K]]: partition for K clusters")
  cat("\n")
  cat("     $clusters: groups membership (before and after consolidation)")
  cat("\n")
  cat("     $comp: latent components of the clusters (after consolidation)")
  if ((EXTu==1)|(EXTr==1)){
    cat("\n")
    cat("     $loading: loadings of the external variables (after consolidation)")
  }
  cat("\n")

}