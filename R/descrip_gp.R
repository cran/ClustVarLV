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
       #pk<- table(resclv$clusters[2,])
} 



if (is.null(colnames(X)))  stop("Please set column names for the matrix X")
# pretreatment of X as in clv
X<- scale(X, center=T, scale=resclv$param$sX)
p <- dim(X)[2] 
n <- dim(X)[1]


clean_var<-NULL
if(resclv$param$strategy=="kplusone")  {
  clean_var<-which(clusters==0)
  pk<- table(clusters)
  if (length(pk)>K) names(pk)[1]<-"noise cluster"
  if (length(clean_var)==p)   stop ("All variables are considered as noise")
}
if(resclv$param$strategy=="sparselv"){
  names0=c()
  sloading = resclv$sloading
  for(k in 1:K) names0 = c(names0,dimnames(sloading[[k]])[[1]])
  names1 = names0[which(unlist(sloading)==0)]
  names2= dimnames(resclv$clusters)[[2]]
  clean_var = sort(match(names1,names2))
  names(clean_var) = colnames(X)[clean_var]  
  clusters[clean_var] = 0
  pk<- table(clusters)
  if (length(pk)>K)  names(pk)[1]<-"zero loading"
  if (length(clean_var)==p)   stop ("All variables have a zero loadings")
}
if(resclv$param$strategy=="none") pk<-table(clusters)

# initialisation  
tab<-vector(length=K) 
correlation<-matrix(nrow=1,ncol=K)
colnames(correlation)<-paste("group",c(1:K))
groups<-list(NA)
# latent var must be set as c'c=1
lvstd<-latvar/sqrt(matrix(apply(latvar,2,var)*(n-1),n,K,byrow=T))
    
    
for (k in 1:K) 
  {
  Xgroup<-as.matrix(X[,clusters==k])
  colnames(Xgroup)<-colnames(X)[which(clusters==k)]
  
  caract<- matrix(0,nrow=ncol(Xgroup),ncol=2)
  rownames(caract)<-colnames(Xgroup)
  colnames(caract)<-c("cor in group","cor next group")
      
  for (j in 1:ncol(Xgroup))    {
    veccov<-cov(Xgroup[,j],lvstd) #covariance between var j and the latent variable of its group (k)
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
    if (method==1) groups[[k]]<-caract[order(abs(caract[,1]),decreasing =T),] 
    if (method==2) groups[[k]]<-caract[order(caract[,1],decreasing =T),]
  }
  
  # sign modification if necessary
  if (method==1) {
    if (groups[[k]][1,1]<0) groups[[k]][,1]<-groups[[k]][,1]*(-1)
    if(K >1){
      if (groups[[k]][1,2]<0) groups[[k]][,2]<-groups[[k]][,2]*(-1)
    }
  }
  
  
}
  

  
  # percentage of the variation explained by the group's LV, within each group
  # percentage of the total variation explained by the K LV
  if (method==1) {
   prop_within<-matrix(0,ncol=K,nrow=1)
   colnames(prop_within)<-paste("Group.",1:K,sep="")
   CLVcp<-latvar
#   CLVcp<-matrix(0,nrow=n,ncol=K)
    for(k in 1:K) {
       Xgroup<-as.matrix(X[,clusters==k])
#      ressvd<-svd(Xgroup)
#      eigval<-ressvd$d^2
#      CLVcp[,k]<-ressvd$u[,1]*ressvd$d[1]
#      eigval<-svd(Xgroup)$d^2
#      prop_within[k]<-eigval[1]/sum(eigval)
       prop_within[k]<-var(CLVcp[,k])/sum(apply(Xgroup,2,var))
   }
   prop_tot<-sum(apply(CLVcp,2,var))/sum(apply(X,2,var))
  }

  
 
  corlv<-round(cor(latvar),2)
  names(corlv) = NULL
  
  if(resclv$param$strategy=="sparselv") {
     if (method==1)   return(list(number=pk,prop_within=round(prop_within,4), prop_tot=round(prop_tot,4), groups=groups,zero_loading=clean_var,cormatrix=corlv))
     if (method==2)   return(list(number=pk,groups=groups,zero_loading=clean_var,cormatrix=corlv))
  } else {
    if (method==1)   return(list(number=pk,prop_within=round(prop_within,4), prop_tot=round(prop_tot,4), groups=groups,set_aside=clean_var,cormatrix=corlv))
    if (method==2)   return(list(number=pk,groups=groups,set_aside=clean_var,cormatrix=corlv))
  }
}
