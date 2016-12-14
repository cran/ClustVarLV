#' @title sparse loadings in each cluster when using the "sparselv" strategy
#' 
#' @description Applies only on CLV_kmeans output with strategy="sparselv".
#' 
#' @param resclv : result of CLV_kmeans() 
#'  
#' @return \item{sparse_loadings}{the loadings of the variables for each latent variables when the "sparselv strategy is used.}
#'        
#' @export
#' 
get_sparseload <-
  function(resclv)
  {
    
   if (!inherits(resclv, c("clv")))   stop("non convenient objects")
   if (!(is.null(resclv$param$nmax)&(resclv$param$strategy=="sparselv")))   
                   stop("only for output of strategy \"sparselv\" from CLV-kmeans")
  
      
   if(resclv$param$strategy=="sparselv") 
                   return(sparse_loadings=resclv$sload)
      
  }