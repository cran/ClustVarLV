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
CLV_kmeans <- function(X,Xu=NULL,Xr=NULL,method,sX=TRUE,sXr=FALSE,sXu=FALSE,init, iter.max=20, nstart=1)
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
    if (p != ptilde) {stop("X and Xu must be defined for the same number of variables") }
    if (EXTr==1) {stop("this procedure doesn't allow Xr and Xu to be defined simultenaously. Use l-clv instead")}
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
                   
if (K==1){
  initgroupes<-rep(1,p)
  lastgroupes<-initgroupes
  ind<-which(groupes == 1)    
  res<-consol_calcul(method,X,EXTr,Xr,EXTu,Xu,ind) 
  critere<-res$critere
  crit_trials<-sum(critere)
  comp<-as.matrix(res$comp)
  i<-1
  if(EXTr==1) a<-as.matrix(res$a)
  if(EXTu==1) u<-as.matrix(res$u)
}

if (K>1) {
cc_consol <- as.matrix(as.numeric(groupes))       
for (i in 1:iter.max) {
  critere <-rep(0,K)
  groupes_tmp<-cc_consol[,i]
  for (k in 1:K) { 
    ind<-which(groupes_tmp==k)
    if (length(ind) > 0) {          
      res <- consol_calcul(method,X,EXTr,Xr,EXTu,Xu,ind)  
      critere[k]<-res$critere
      comp[,k]<-res$comp
      if (EXTr==1)  a[,k]<-res$a
      if (EXTu==1)  u[,k]<-res$u
   }
  }
  groupes_tmp<-consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
  if (length(which((cc_consol[,i] == groupes_tmp) == FALSE, arr.ind = T)) == 0)    break
  cc_consol = cbind(cc_consol, groupes_tmp)
}
rownames(cc_consol) <- colnames(X)      
names(cc_consol) = NULL
initgroupes<-cc_consol[,1]
lastgroupes<-cc_consol[,ncol(cc_consol)]
iter<-i
crit_trials<-sum(critere)

 if (nstart >= 2) {
     best <- sum(critere)
      for (i in 2:nstart) {
          out<-mat_init(X,EXTr,Xr,EXTu,Xu,K)
          comp2<-out$comp
          comp2 <- X[,sort(sample.int(p, K))]
          if (EXTr==1)  {  a2<-out$a }
          if (EXTu==1)  {  u2<-out$u }
          groupes2 <- as.factor(consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp2,a2,u2))
          cc_consol2 <- as.matrix(as.numeric(groupes2))
          for (i in 1:iter.max) {
            critere2 <-rep(0,K)
            groupes_tmp2<-cc_consol2[,i]
            for (k in 1:K) {
              ind2<-which(groupes_tmp2==k)
              if (length(ind2) > 0) {
               res2 <- consol_calcul(method,X,EXTr,Xr,EXTu,Xu,ind2) 
               critere2[k]<-res2$critere
               comp2[,k]<-res2$comp
               if (EXTr==1)  a2[,k]<-res2$a
               if (EXTu==1)  u2[,k]<-res2$u
              }
            }
            groupes_tmp2<-consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp2,a2,u2)
            if (length(which((cc_consol2[, i] == groupes_tmp2) == FALSE, arr.ind = T)) == 0)    break
            cc_consol2 = cbind(cc_consol2, groupes_tmp2)
          }
          crit_trials<-c(crit_trials,sum(critere2))
          if ((zz <- sum(critere2)) > best) {
              comp<-comp2
              critere<-critere2
              if (EXTr==1) a<-a2
              if (EXTu==1) u<-u2
              initgroupes<-cc_consol2[,1]
              lastgroupes<-cc_consol2[,ncol(cc_consol2)]
              iter<-i
              best <- zz
          }
     }
 }
}

tabres<-as.matrix(t(c(sum(critere),(sum(critere)/sbegin)*100, i)))
colnames(tabres)<-c("clust.crit.cc","%S0expl.cc","nbiter")
clusters=rbind(initgroupes,lastgroupes)
colnames(clusters) <- colnames(X)
param<-list(method=method,n = n, p = p,K = K,nstart = nstart,EXTu=EXTu,EXTr=EXTr,sX=sX,sXr=sXr,cXu=cXu,sXu=sXu)
listcc<-list(tabres=tabres,clusters=clusters, comp=comp, trials=crit_trials,param=param)
if(EXTr==1) listcc= c(listcc, list(loading=a)) 
if(EXTu==1) listcc= c(listcc, list(loading=u)) 
  
class(listcc) = "clvkmeans"
return(listcc)
}
