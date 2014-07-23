#' K-means algorithm for the clustering of variables
#'
#' K-means algorithm for the clustering of variables. Directional or local groups may be defined. 
#' Each group of variables is associated with a latent component. 
#' Moreover external information collected on the observations or on the variables may be introduced.
#' 
#' The initalization can be made at random, repetitively, or can be defined by the user.
#'
#' @param X The matrix of the variables to be clustered
#' @param Xu The external variables associated with the columns of X
#' @param Xr The external variables associated with the rows of X
#' @param method The criterion to use in the cluster analysis.\cr 
#'        1 : the squared covariance is used as a measure of proximity (directional groups).\cr
#'        2 : the covariance is used as a measure of proximity (local groups)
#' @param sX TRUE/FALSE : standardization or not of the columns X (TRUE by default)\cr
#'        (predefined -> cX = TRUE : column-centering of X)
#' @param sXr TRUE/FALSE : standardization or not of the columns Xr (FALSE by default)\cr
#'        (predefined -> cXr    = TRUE : column-centering of Xr)
#' @param sXu TRUE/FALSE : standardization or not of the columns Xu (FALSE by default)\cr
#'        (predefined -> cXu= FALSE : no centering, Xu considered as a weight matrix)
#' @param init a number i.e.  the size of the partition, K,
#'        or  a vector of INTEGERS i.e. the group membership of each var in the initial parition (integer between 1 and K)
#' @param iter.max maximal number of iteration for the consolidation (20 by default)
#' @param nstart nb of random initialisations in the case of init = a number  (1 by default)
#' @param strategy "none", or "kplusone" (an additional cluster for the noise variables, default),
#'        or "sparselv" (zero loadings for the noise variables)
#' @param rho a threshold of correlation between 0 and 1 (0 by default)
#' 
#' @return \item{tabres}{ 
#'         The value of the clustering criterion at convergence.\cr
#'         The percentage of the explained initial criterion value.\cr
#'         The number of iterations in the partitioning algorithm.}
#'         \item{clusters}{ the group's membership}
#'         \item{comp}{ The latent components of the clusters}
#'         \item{loading}{ if there are external variables Xr or Xu :  The loadings of the external variables}
#' @seealso CLV, LCLV
#' @examples data(apples_sh)
#' #local groups with external variables Xr 
#' resclvkmYX <- CLV_kmeans(X = apples_sh$pref, Xr = apples_sh$senso, 
#'                          method = 2, sX = FALSE, sXr = TRUE, init = 2, nstart = 20)
CLV_kmeans <- function(X,Xu=NULL,Xr=NULL,method,sX=TRUE,sXr=FALSE,sXu=FALSE,
                       init, iter.max=20, nstart=1,strategy="none",rho=0)
{
  if(method!=1 & method!=2) stop("method should be 1 or 2")
  cX=TRUE
  cXr=TRUE
  cXu=FALSE
  X<- scale(X, center=cX, scale=sX)
  p <- ncol(X)
  n <- nrow(X)  

  if (is.null(Xr)) {
    EXTr<-0
  }
  else {
    EXTr<-1                                    
    Xr<- scale(Xr, center=cXr, scale=sXr)
    ntilde <- dim(Xr)[1]
    q<-dim(Xr)[2] 
    if (n != ntilde) stop("X and Xr must have the same number of observations")                              
  }

  if (is.null(Xu)) {
    EXTu<-0
  } 
  else {
    EXTu<-1                   
    Xu<- scale(Xu, center=cXu, scale=sXu) 
    ptilde <- dim(Xu)[1]  
    m<-dim(Xu)[2] 
    if (p != ptilde) {stop("X and Xu must be defined for the same number of
                           variables") }
    if (EXTr==1) {stop("this procedure doesn't allow Xr and Xu to be defined
                       simultenaously. Use l-clv instead")}
  }  
  
  
  crit<-crit_init(method,X,EXTr,Xr,EXTu,Xu)
  sbegin <- sum(crit)  
  
 if (missing(init))
     stop("'init' must be a number or a vector")
 if (length(init) == 1) {
     K <- init
     out<-mat_init(X,EXTr,Xr,EXTu,Xu,K)
     comp<-out$comp
     comp <- as.matrix(X[,sort(sample.int(p, K))]) # K columns of X chosen at random
     if (EXTr==1)  {  a<-out$a }
     if (EXTu==1)  {  u<-out$u }
     groupes <- as.factor(consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp,a,u))
  } else {
     nstart = 1
     if (!is.numeric(init))
         stop("init must be a vector of integers")
     groupes <- as.factor(init)
     K <- length(levels(groupes))
     if (p < K)
         stop("more cluster's centers than variables")
     if (length(which(init > K)) > 0)
         stop("clusters must be numbered from 1 to K")
     if (p != length(groupes))
         stop("the length of init must be equal to the number of variables")
     out<-mat_init(X,EXTr,Xr,EXTu,Xu,K)
     comp<-out$comp
     if (EXTr==1)  {a<-out$a}
     if (EXTu==1)  {u<-out$u}
 }
           
  if (length((intersect(strategy,c("kplusone","sparselv","none"))))==0) 
      stop("strategy must be either 'kplusone', 'sparselv', 'none'")
 
  if (strategy=="kplusone" & (EXTu!=0 | EXTr!=0) )
       stop(" 'k+1' strategy is not available with external variables, yet")
  if (strategy=="sparselv" & (EXTu!=0 | EXTr!=0) )
       stop(" 'Sparse LV' strategy is not available with external variables, yet")
 
  #####################################################################
  if (strategy=="none" | rho==0 ) 
    listcc = clvk_none(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart)
  if (strategy=="sparselv")       
    listcc = clvk_sparse(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart,rho)
  if (strategy=="kplusone")       
    listcc = clvk_kp1(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart,rho)
  #####################################################################


param<-list(method=method,n = n, p = p,K = K,nstart = nstart,EXTu=EXTu,EXTr=EXTr,
            sX=sX,sXr=sXr,cXu=cXu,sXu=sXu,strategy=strategy,rho=rho)
listcc= c(listcc, list(param=param)) 


class(listcc) = "clvkmeans"
# if (strategy=="sparselv") class(listcc) = "sparselv"
# if (strategy=="kplusone") class(listcc) = "kplusone"

return(listcc)
}
