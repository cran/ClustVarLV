  #' linear model based on CLV : prediction of a response variable, y, based on clusters of explanatory variables, X.
  #'
  #' boosted-liked procedure for identifying groups of explanatory variables (and the associated latent variable) specifically correlated with the response variable, y.
  #' Directional groups are considered. 
  #' Discarding spurious variables is allowed using the strategy and rho parameters. 
  
  #' @param X : The matrix of the explanatory variables, to be clustered
  #' @param y : The response variable (numeric)
  #' @param method : The criterion to be use in the cluster analysis.  \cr
  #'        1 or "directional" : the squared covariance is used as a measure of proximity (directional groups). \cr    
  #'        2 or "local"       : the covariance is used as a measure of proximity (local groups)
  #' @param sX : TRUE/FALSE, i.e. standardization or not of the columns X (TRUE by default)\cr
  #'        (predefined -> cX = TRUE : column-centering of X)
  #' @param nbiter : maximum number of steps (by default nbiter=100)
  #' @param strategy : "none" (by default), or "kplusone" (an additional cluster for the noise variables),
  #'        or "sparselv" (zero loadings for the noise variables)
  #' @param rho : a threshold of correlation between 0 and 1 (used for kplusone or sparselv strategies, 0.3 by default)
  #' @param shrinkp : shrinkage  paramater used in the boosting (1 by default)
  #' @param validation TRUE/FALSE i.e. using a test set or not. By default no validation
  #' @param id.test : if validation==TRUE, the number of the observations used as test set
  #' 
  #' @return \item{Group}{a list of the groups of variables X iteratively extracted.}
  #' @return \item{Comp}{ a list of the latent variables associated with the groups of X variables extracted.}
  #' @return \item{Load}{ a list for the loadings of the X variables according to the associated latent variable.}
  #' @return \item{Beta}{ a list for the coefficients associated to the scaled predictors (X must to be scaled accordind to sX) regarding y .}
  #'
  #' @seealso CLV, CLV_kmeans
  #'
  #' @export
  #' 
lm_CLV = function(X,y,method="directional",sX=sX,nbiter=100, strategy="none",rho=0.3,shrinkp=1, validation=FALSE,id.test=NULL )
{  
 
  
  Group = list()
  Load=list()
  Comp =list()
  Beta=matrix(0,nbiter,ncol(X))
  n = dim(X)[1]
  p = dim(X)[2]
  gpmax=p

  
  if (!validation) {
    Xscale=scale(X,scale=sX)
    ycenter=scale(y,scale=FALSE)
  }
  if (validation) {
    test<-id.test
    train<-setdiff(1:n,test)
    meany=mean(y[train])
    scale.mean<-apply(X[train,],2,mean)
    Xtest=X[test,]
    ytest=y[test]
    if (sX==TRUE) {
      scale.scale<-apply(X[train,],2,sd)
    } else {
      scale.scale=rep(1,ncol(X))
    }
    X=X[train,]
    y=y[train]
    n<-length(train)
    Xscale=scale(X,scale=sX)
    ycenter=scale(y,scale=FALSE)
    ERR.val=matrix(NA,length(id.test),nbiter+1)
    ERR.val[,1]<-ytest-rep(meany,length(ytest))
  }
  ERR.cal=matrix(NA,n,nbiter+1)
  ERR.cal[,1]<-y-rep(mean(y),length(rep))
 
  # CLV  -----------------------------------------------
  # performed on X, one time only
  # ascendant hierarchy, with consolidation
  res.clv1_X=CLV(X,sX=sX,method=method,nmax=gpmax,maxiter=10)  # avec consolidation 
  
 
 # ------------------------------------------------------------
 # iterative search for groups of X variables (and associated LV) well correlated with y
  yt=ycenter
  for (iter in 1:nbiter) {
      tabgpy=matrix(0,p,gpmax)
      relate = list()
      for (k in 1:gpmax) {
        #gp=res.clv1_X[[k]]$clusters[1,]
        gp=get_partition(res.clv1_X,k)
        comp=get_comp(res.clv1_X,k)
        gpy=which.max(abs(cor(comp,yt)))
        tabgpy[which(gp==gpy),k]=1
        relate[[k]]=which(gp==gpy)
      }
   
      # identify the levels for which the X variables associated with y change
       nivparti=NULL
      for (k in 1:(length(relate)-1)) {
        respons=length(setdiff(relate[[k]],relate[[k+1]]))+length(setdiff(relate[[k+1]],relate[[k]]))
        if (respons>0) nivparti=c(nivparti,k) 
      }
      
      # consolidation, only for the levels in nivparti
      # CLV_kmeans according to "strategy". initialisation by cutting the hierarchy in res.clv1_X
      res.clv1Xy=NULL
      if ((strategy=="sparselv")|(strategy=="kplusone")) {
       for (kk in 1:length(nivparti)) {
         K=nivparti[kk]
         res.clv1Xy[[kk]]=CLV_kmeans(X,sX=sX,method="directional",clust=res.clv1_X[[K]]$clusters[1,],strategy=strategy,rho=rho)     
         names(res.clv1Xy)[[kk]]=paste("partition",nivparti[[kk]],sep="")
       }
      } 
      if (strategy=="none") {
        for (kk in 1:length(nivparti)) {
          K=nivparti[kk]
          res.clv1Xy[[kk]]=res.clv1_X[[K]]     
          names(res.clv1Xy)[[kk]]=paste("partition",nivparti[[kk]],sep="")
        }
      }
      res.clv1Xy$param=res.clv1_X$param
      res.clv1Xy$tabres=res.clv1_X$tabres[p-nivparti,]

       #----------------------------------------------------------------------
    
    # among the groups well associated with y, choice of one group
    # criterion here : modified Kaiser Guttman criterion on the standardized matrix of the X belonging to this group
      
      nbniv=length(nivparti)
      # let the X variables associated with y : is this group unidimensionnal ?
      rescrit=NULL
      tabgpyniv=matrix(0,p,nbniv)
      for(j in 1:nbniv){
        aclevel = j
        if ((strategy=="sparselv")|(strategy=="kplusone")) {
         gp = res.clv1Xy[[aclevel]]$clusters[2,]
         comp=res.clv1Xy[[aclevel]]$comp
         gpy=which.max(abs(cor(comp,yt)))
         tabgpyniv[which(gp==gpy),aclevel]=1
         indiy=which(gp==gpy)
        }
        if (strategy=="none") {
          tabgpyniv[,aclevel]=tabgpy[,nivparti[aclevel]]
          indiy=relate[[nivparti[aclevel]]]
        }
        #print(indiy)
        py=length(indiy)
        if (py==0) {break}
        
        Xj=as.matrix(X[,indiy])
        Xjsc=scale(Xj)
        svdXjsc=svd(Xjsc/sqrt(n-1))
        # KG
        l1=svdXjsc$d[1]^2
        l2=svdXjsc$d[2]^2
        seuill=1+ ( 2* sqrt((py-1)/(n-1)))
        KGX=(l1>seuill)&(l2<seuill)
   
        rescrit=rbind(rescrit,c(aclevel=nivparti[j],py=py,l1=l1,l2=l2,seuill=seuill,KGX=KGX))
     
     }  # fin boucle j
     
    rescrit=as.data.frame(rescrit)
    nf=rescrit[which(rescrit$KGX==1)[1],1];
                           #print(nf)
    
    if (is.na(nf)) {
      Group[[iter]]=NULL
      Coeff[[iter]]=NULL
      Comp[[iter]]=NULL
      iter=nbiter
    } else {
     ind= which(tabgpyniv[,which(nivparti==nf)]==1)
    if( strategy=="sparselv") {
      ind=intersect(ind,which(get_sparseload(res.clv1Xy[[which(nivparti==nf)]],"vector")!=0))
    }
     indname=colnames(X)[ind]
     Group[[iter]]=indname
                            #print(Group[[iter]])
     
     Xiter=as.matrix(Xscale[,ind])
     #peig=powerEigen(crossprod(Xiter))
     #load=peig$vectors
     vp <- eigen( t(t(Xiter)%*%Xiter)/(n-1))
     load=as.matrix(vp$vectors[,1])
     rownames(load)=colnames(Xiter)
     comp=Xiter%*%load
   
     Load[[iter]]=load
     Comp[[iter]]=comp
 
     # newt step unsing a boosted-like procedure
       beta=shrinkp*(cov(yt,comp)/var(comp))
       Beta[iter,which(apply(sapply(indname,FUN="==",colnames(X)),1,sum)==1)]=as.numeric(beta)*load
       yiter=comp%*%beta
       yt= as.vector( yt - yiter )
       # rq all the X variables are used as previoulsly, then the CLV results also
      
   }
        # part for prediction/ step "iter"
    if (iter==1) W=Beta[1,]
    if (iter>1)  W=apply(as.matrix( Beta[1:iter,]),2,sum)
    ypred=mean(y)+Xscale%*%W
    ERR.cal[,iter+1]=y-ypred
    if (validation){
      Xstest<-scale(Xtest,center=scale.mean,scale=scale.scale)
      ypredtest=meany+Xstest%*%W
      ERR.val[,iter+1]=ytest-ypredtest
    } 
    
  } # end of iter
                             #print(paste("nb iterations:" ,iter))
  
    if (!validation)  res=list(Group=Group,Load=Load,Comp=Comp,Beta=Beta,ERR.cal=ERR.cal)
    if (validation)   res=list(Group=Group,Load=Load,Comp=Comp,Beta=Beta,ERR.cal=ERR.cal,ERR.val=ERR.val)
    return(res)
}
