#' Description of the clusters of variables
#' 
#' Complementary function that provides for each cluster their correlation with their own cluster latent component
#' and with the next neighbouring cluster latent component.
#' 
#' @param resclv resultat of CLV or CLV_kmeans
#' @param X the initial matrix used in CLV or CLV_kmeans
#' @param K the number of groups in the partition
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = 1, sX = TRUE, graph = TRUE)
#' descrip_gp(resclvX, X = apples_sh$senso, K = 4) 
descrip_gp <-
function(resclv, X, K=NULL) {
  if (!inherits(resclv, c("clv","clvkmeans"))) 
    stop("non convenient objects")
# resclv : resultat of clv.r or clv_keans.r
# X: the initial matrix 
# the number of groups in the partition

method<-resclv$param$method  
# group's membership of the variables  
if(is.null(resclv$param$K)) { 
    if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
       clusters<-resclv[[K]]$clusters[2,]
       latvar = resclv[[K]]$comp
       pk<- table(resclv[[K]]$clusters[2,])
} else {
       clusters<-resclv$clusters[2,] 
       K<-resclv$param$K
       latvar = resclv$comp
       pk<- table(resclv$clusters[2,])
} 
# pretreatment of X as in clv
X<- scale(X, center=T, scale=resclv$param$sX)
p <- dim(X)[2] 
n <- dim(X)[1]

# initialisation  
tab<-vector(length=K) 
correlation<-matrix(nrow=1,ncol=K)
colnames(correlation)<-paste("group",c(1:K))
groups<-list(NA)
#latvar<-resclv[[K]]$comp

    
for (k in 1:K) 
  {
  Xgroup<-as.matrix(X[,clusters==k])
  colnames(Xgroup)<-colnames(X)[which(clusters==k)]
  
  caract<- matrix(0,nrow=ncol(Xgroup),ncol=2)
  rownames(caract)<-colnames(Xgroup)
  colnames(caract)<-c("cor in group","cor next group")
    
  for (j in 1:ncol(Xgroup))    {
    veccov<-cov(Xgroup[,j],latvar) #covariance between var j and the latent variable of ist group (k)
    veccor<-cor(Xgroup[,j],latvar) #correlation between var j and the latent variable of ist group (k)
    if (method==1) ordrecov<-order(veccov^2,decreasing=T)
    if (method==2) ordrecov<-order(veccov,decreasing=T)
    verif<-(ordrecov[1]==k) 
    if (verif==F) {print (c(k,j)) }
    caract[j,1]<-round(veccor[ordrecov[1]],2) 
    caract[j,2]<-round(veccor[ordrecov[2]],2)
  }
  
  if (nrow(caract)==1) {
    groups[[k]]<-caract  
  }else{
    groups[[k]]<-caract[order(abs(caract[,1]),decreasing =T),] 
  }
 }

#   pk<- table(resclv[[K]]$clusters[2,]) 
  corlv<-round(cor(latvar),2)
  names(corlv) = NULL

return(list(number=pk,cormatrix=corlv,groups=groups))
  
}
