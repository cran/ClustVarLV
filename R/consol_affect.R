consol_affect <-
function(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
{   
  p<-ncol(X)
  gtmp<-rep(0,p)
  n<-nrow(X)
 
  if (method == 1){
    if ((EXTu==0)&(EXTr==0)) {
      cova = t(comp) %*% X/(n-1)
      covaC<-  cova^2
      for (j in (1:p)) {
        maxj<-which(covaC[,j]==max(covaC[,j]))
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
        maxj<-which(covaC[,j]==max(covaC[,j]))
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      cova = t(comp) %*% X/(n-1)
      covaC <- cova^2
      for (j in (1:p)) {
        maxj<-which(covaC[,j]==max(covaC[,j]))
        gtmp[j] = maxj[1]   
      }
    }   
  }
  
  if (method==2){
    if ((EXTu==0)&(EXTr==0)) {
      cova = t(comp) %*% X/(n-1)
      covaC<-  cova
      for (j in (1:p)) {
        maxj<-which(covaC[,j]==max(covaC[,j]))
        gtmp[j] = maxj[1]   
      }
    } 
    if ((EXTu==0)&(EXTr==1)) {
      cova = t(comp) %*% X/(n-1)
      covaC <- cova
      for (j in (1:p)) {
        maxj<-which(covaC[,j]==max(covaC[,j]))
        gtmp[j] = maxj[1]   
      }
    }  
    if ((EXTu==1)&(EXTr==0)){          
      for (j in 1:p) {
        Pj<- as.matrix(X[,j])%*%Xu[j,]
        critind=diag(t(u)%*%t(Pj)%*%comp)          
        maxj=which(critind==max(critind))
        gtmp[j] = maxj[1]                    
      }
    }   
  } 
  
  return(gtmp)
}
