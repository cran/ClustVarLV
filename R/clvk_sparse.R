clvk_sparse <- function(X,n,p,sbegin,EXTr,EXTu,Xu,Xr,method,K,comp,groupes,a,u,iter.max,nstart,rho)
{                      
                                                             
# verification if there are NA values
  valmq=FALSE
  if (sum(is.na(X))>0)  {
    valmq=TRUE
    tauxNA=sum(is.na(X))/(n*p)
  }
epsil=0.00001 
sloading = list()
################################################
# when the variables form only one cluster  
if (K==1){
  cc_consol <- as.matrix(rep(1,p))
  for (i in 1:iter.max) {
    critere <-rep(0,K)
    groupes_tmp<-cc_consol[,i]
    ind<-which(groupes_tmp == 1)   
     res = consol_calcul_s(method,X,EXTr,Xr,EXTu,Xu,ind,rlevel = rho) 
  
    critere<-res$critere
    crit_trials<-sum(critere)
    comp<-as.matrix(res$comp)
    
    sloading[[1]] <-as.matrix(res$loading)
    
    if(EXTr==1) a<-as.matrix(res$a)
    if(EXTu==1) u<-as.matrix(res$u)
   
    # re-allocation of the variables

    groupes_tmp<-consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
    if (length(which((cc_consol[,i] == groupes_tmp) == FALSE, arr.ind = T)) == 0) break
    cc_consol = cbind(cc_consol, groupes_tmp)
  }
  crit_trials<-sum(critere)
  initgroupes<-cc_consol[,1]
  lastgroupes<-cc_consol[,ncol(cc_consol)]
  iter<-i
}

################################################
# for K (>1) clusters  
if (K>1) {
cc_consol <- as.matrix(as.numeric(groupes)) 
pcritav=0
                                                
for (i in 1:iter.max) {

  critere <-rep(0,K)
  groupes_tmp<-cc_consol[,i] 
 
  for (k in 1:K) { 
    ind<-which(groupes_tmp==k) 
                               
    if (length(ind) > 0) {    
        res <- consol_calcul_s(method,X,EXTr,Xr,EXTu,Xu,ind,rlevel=rho)  
        sloading[[k]] = res$loading
    
      critere[k]<-res$critere
      comp[,k]<-res$comp
       
      if (EXTr==1)  a[,k]<-res$a
      if (EXTu==1)  u[,k]<-res$u
   }
  }
   pcrit<-sum(critere)/sbegin
 
  # re-allocation of the variables
 
  groupes_tmp<-consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
  if (length(which((cc_consol[,i] == groupes_tmp) == FALSE, arr.ind = T)) == 0)    break
  if((pcrit-pcritav)<epsil) break
  cc_consol = cbind(cc_consol, groupes_tmp)
  pcritav<-pcrit
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
          
          sloading2 = list()
          
          if (EXTr==1)  {  a2<-out$a }
          if (EXTu==1)  {  u2<-out$u }
          groupes2 <- as.factor(consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp2,a2,u2))
          cc_consol2 <- as.matrix(as.numeric(groupes2))
          pcrit2av=0
          
          for (i in 1:iter.max) {
       
          critere2 <-rep(0,K)
           
            groupes_tmp2<-cc_consol2[,i]
            for (k in 1:K) {
              ind2<-which(groupes_tmp2==k)
              if (length(ind2) > 0) {
            
                 res2 <- consol_calcul_s(method,X,EXTr,Xr,EXTu,Xu,ind2,rlevel = rho) 
                 sloading2[[k]] = res2$loading     
         
               critere2[k]<-res2$critere
               comp2[,k]<-res2$comp
               if (EXTr==1)  a2[,k]<-res2$a
               if (EXTu==1)  u2[,k]<-res2$u
              }
            }
            pcrit2<-sum(critere2)/sbegin
            
            groupes_tmp2<-consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp2,a2,u2)
            if (length(which((cc_consol2[, i] == groupes_tmp2) == FALSE, arr.ind = T)) == 0)    break
            if((pcrit2-pcrit2av)<epsil) break
            cc_consol2 = cbind(cc_consol2, groupes_tmp2)
            pcrit2av<-pcrit2
          }
                                                      
          crit_trials<-c(crit_trials,sum(critere2))   
          if ((zz <- sum(critere2)) > best) {
              comp<-comp2
              critere<-critere2
              
              sloading = sloading2
              
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

tabres<-as.matrix(t(c(sum(critere),(sum(critere)/sbegin)*100, iter)))
colnames(tabres)<-c("clust.crit.cc","%S0expl.cc","nbiter")
colnames(comp)<-paste("Comp",c(1:K),sep="")
 
  
clusters=rbind(initgroupes,lastgroupes)
colnames(clusters) <- colnames(X)
listcc<-list(tabres=tabres,clusters=clusters, comp=comp, trials=crit_trials)

listcc= c(listcc, list(sloading=sloading))
if(EXTr==1) listcc= c(listcc, list(loading=a)) 
if(EXTu==1) listcc= c(listcc, list(loading=u)) 
  
return(listcc)
}
