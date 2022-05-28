#'@title Boostrapping for assessing the stability of a CLV result
#'
#'@description
#'Bootstrapping on the individuals (in row) or the variables (in column) is performed. 
#'Choose the "row" option (the default), if the variables are measured on a random sample of individuals,
#'Choose the "column" option, if the variables are taken from a population of variables.
#'The first case is the more usual, 
#'but the second may occur, e.g. when variables are consumers assessing specific products. 
#'Each boostrapped data matrix is submitted to CLV in order to get partitions from 1 to nmax clusters.
#'For each number of clusters, K, the Rand Index, the adjusted Rand Index, 
#'as well as the cohesion and the isolation of the clusters
#'of the observed partition and the bootstrapped partitions are computed. 
#'These criteria are used for assessing the stability of the solution into K clusters.
#'Parallel computing is performed for time saving.
#' 
#' @param object : result of CLV()
#' @param case : "row" or "column" corresponding to the random effect of the design 
#' @param B : the number of bootstrap to be run (100 by default)
#' @param nmax : maximal size of the partitions to be considered (if NULL, the value of nmax used for the object is used)
#' 
#' @return \item{res}{a list of length 4 for the Rand Index, Adjusted Rand Index, Cohesion and Isolation of the partition \cr
#'                   results (matrix of size (B x nmax), respectively.}
#'
#' @seealso CLV
#' 
#' @export


 boot_clv =  function (object, case="row", B=100,nmax=NULL) 
 {
  
    resclv<-object   
    if (!inherits(resclv, "clv"))   stop("non convenient objects")
    
    if(!((case=="row")|(case=="column"))) stop('case parameter should be "row" or "column"')
    
    X<-resclv$param$X
    if(is.null(nmax)) nmax=resclv$param$nmax
    method<-resclv$param$method
    n<-resclv$param$n
    p<-resclv$param$p
    sX<-resclv$param$sX
    
  
    
    # boot_fx_row ------------------------------------------
    boot_fx_row <- function(b) {
      resari=c()
      resrand=c()
      rescohe=list()
      rescoheP=c(NA)
      resisol=list()
      resisolP=c(NA)
      bsamp=sample(1:n,replace=TRUE)
      if (p<=1500) {
        resboot=CLV(X[bsamp,],method=method,sX=sX,nmax=nmax)
        for (k in 1:nmax)  {
          parti_ref=get_partition(resclv,k)
          parti_boot=get_partition(resboot,k)
          resrand=c(resrand,ARI(parti_ref,parti_boot)[1])
          resari=c(resari,ARI(parti_ref,parti_boot)[2])
          rescic=cohe_isol_C(parti_ref,parti_boot)
          rescohe[[k]]=rescic$cohe
          resisol[[k]]=rescic$isol
          sizegp=as.vector(table(parti_ref))
          alphaC=sizegp*(rep(p,k)-sizegp)        # mC*(m-mC), nb pairs of objects in C and Cbar
          betaC=sizegp*(sizegp-1)/2              # mC * (mC-1)/2 , nb pairs of objects in C
          if (k>1) {
            rescoheP=c(rescoheP,sum(betaC*rescohe[[k]])/ sum(betaC))
            resisolP=c(resisolP,sum(alphaC*resisol[[k]]) /sum(alphaC))
          }
        } 
      } else {
        for (k in 1:nmax) {
          parti_ref=get_partition(resclv,k)
          resboot=CLV_kmeans(X[bsamp,],method=method,sX=sX,clust=k,nstart=50)
          parti_boot=get_partition(resboot)
          resrand=c(resrand,ARI(parti_ref,parti_boot)[1])
          resari=c(resari,ARI(parti_ref,parti_boot)[2])
          rescic=cohe_isol_C(parti_ref,parti_boot)
          rescohe[[k]]=rescic$cohe
          resisol[[k]]=rescic$isol
          sizegp=as.vector(table(parti_ref))
          alphaC=sizegp*(rep(p,k)-sizegp)        # mC*(m-mC), nb pairs of objects in C and Cbar
          betaC=sizegp*(sizegp-1)/2              # mC * (mC-1)/2 , nb pairs of objects in C
          if (k>1) {
            rescoheP=c(rescoheP,sum(betaC*rescohe[[k]])/ sum(betaC))
            resisolP=c(resisolP,sum(alphaC*resisol[[k]]) /sum(alphaC))
          }
        }
      }
      #    
      res=list(Rand=resrand, ARI=resari,CoheP=rescoheP,IsolP=resisolP)
      return(res)
 }  
 
    # boot_fx_col ------------------------------------------
    boot_fx_col<- function(b) {
      resari=c()
      resrand=c()
      rescohe=list()
      rescoheP=c(NA)
      resisol=list()
      resisolP=c(NA)
      bsamp=sample(1:p,replace=TRUE)
      bsamp=sort(bsamp)
      pbsamp=length(unique(bsamp))
      if (p<=1500) {
        resboot=CLV(X[,bsamp],method=method,sX=sX,nmax=min(nmax,pbsamp))
        for (k in 1:min(nmax,pbsamp))  {
          parti_ref=get_partition(resclv,k)
          parti_boot=get_partition(resboot,k)
          # parti_boot_reduced and parti_ref_reduced : partition without repeated variables in bootstrap sample
          parti_ref_reduced=parti_ref[unique(bsamp)]
          parti_boot_reduced=parti_boot[-which(bsamp[-1]==bsamp[-p])]
          # if a group disappear in the parti_ref_reduced
          tab=table(parti_ref_reduced)
          if (length(tab)<k) {
            mq=setdiff(1:k,unique(parti_ref_reduced))
            give=apply(cor(get_comp(resclv,k),get_comp(resboot,k))[,mq,drop=FALSE],2,which.max)
            for (i in 1:length(mq)) {
              whomove=which(parti_ref_reduced==give[i])
              if (length(whomove)>1)  {
                parti_ref_reduced[whomove[1]]=mq[i]
              }
            }
          }
          #
          out=ARI(parti_ref_reduced,parti_boot_reduced)
          resrand=c(resrand,out[1])
          resari=c(resari,out[2])
          rescic=cohe_isol_C(parti_ref_reduced,parti_boot_reduced)
          rescohe[[k]]=rescic$cohe
          resisol[[k]]=rescic$isol
          tab=table(parti_ref_reduced)
          sizegp=as.vector(tab)
          alphaC=sizegp*(rep(sum(tab),k)-sizegp) # mC*(m-mC), nb pairs of objects in C and Cbar
          betaC=sizegp*(sizegp-1)/2              # mC * (mC-1)/2 , nb pairs of objects in C
          if (k>1) {
            rescoheP=c(rescoheP,sum(betaC*rescohe[[k]],na.rm=TRUE)/ sum(betaC))
            resisolP=c(resisolP,sum(alphaC*resisol[[k]],na.rm=TRUE) /sum(alphaC))
          }
        }
      } else {
        for (k in 1:min(nmax,pbsamp)) {
          parti_ref=get_partition(resclv,k)
          resboot=CLV_kmeans(X[,bsamp],method=method,sX=sX,clust=k,nstart=50)
          parti_boot=get_partition(resboot)
          # parti_boot_reduced and parti_ref_reduced : partition without repeated variables in bootstrap sample
          parti_ref_reduced=parti_ref[unique(bsamp)]
          parti_boot_reduced=parti_boot[-which(bsamp[-1]==bsamp[-p])]
          # if a group disappear in the parti_ref_reduced
          tab=table(parti_ref_reduced)
          if (length(tab)<k) {
            mq=setdiff(1:k,unique(parti_ref_reduced))
            give=which.max(cor(get_comp(resclv,k),get_comp(resboot,k))[,mq,drop=FALSE])
            whomove=which(parti_ref_reduced==give)[1]
            parti_ref_reduced[whomove]=mq
          }
          #
          out=ARI(parti_ref_reduced,parti_boot_reduced)
          resrand=c(resrand,out[1])
          resari=c(resari,out[2])
          rescic=cohe_isol_C(parti_ref_reduced,parti_boot_reduced)
          rescohe[[k]]=rescic$cohe
          resisol[[k]]=rescic$isol
          tab=table(parti_ref_reduced)
          sizegp=as.vector(tab)
          alphaC=sizegp*(rep(sum(tab),k)-sizegp)  # mC*(m-mC), nb pairs of objects in C and Cbar
          betaC=sizegp*(sizegp-1)/2              # mC * (mC-1)/2 , nb pairs of objects in C
          if (k>1) {
            rescoheP=c(rescoheP,sum(betaC*rescohe[[k]],na.rm=TRUE)/ sum(betaC))
            resisolP=c(resisolP,sum(alphaC*resisol[[k]],na.rm=TRUE) /sum(alphaC))
            
          }
        }
      }
      #    
      res=list(Rand=resrand, ARI=resari,CoheP=rescoheP,IsolP=resisolP)
      return(res)
    }  
    
    # ARI -------------------------------------------------
    ARI<-function(P1,P2) {
      # compute the Adjusted Rand Index between two partitions
      tab <- table(as.vector(P1), as.vector(P2))
      totm1=apply(tab,1,sum)
      totm2=apply(tab,2,sum)
      tot=sum(tab)
      index=sum(choose(tab,2))
      expect=sum(choose(totm1,2))*sum(choose(totm2,2))/sum(choose(tot,2))
      maxindex=(sum(choose(totm1,2))+sum(choose(totm2,2)))/2
      ari=(index-expect)/(maxindex-expect)
      Rand=(choose(tot,2)+(2*index)-sum(choose(totm1,2))-sum(choose(totm2,2)))/choose(tot,2)
      return(c(Rand,ari))
    } 
    
    # Cohesion and isolation of a cluster------------------
    cohe_isol_C<-function(P1,P2) {
      # compute the values for cohesion and isolation assesment
      p1=length(P1)
      p2=length(P2)
      cohe<-rep(0,max(P1))
      isol<-rep(0,max(P1))
    
      cooc_P2<-matrix(0,nrow=p1,ncol=p1)
      for (i in 1:max(P2)) {
        #num=as.numeric(matrix(unlist(strsplit(names(which(P2==i)),"X")),nrow=2)[2,])
        num=match(names(which(P2==i)),names(P2))
        cooc_P2[num,num]<-1
      }
      #apply(cooc_P2,1,sum)
      diag(cooc_P2)<-1
    
      for (g in 1:max(P1)){
          mC=length(which(P1==g))
          data2=cooc_P2[which(P1==g),which(P1==g)]
          n11C=sum(data2[lower.tri(data2)]==1)
          data4=cooc_P2[which(P1==g),setdiff(1:p1,which(P1==g))]
          n00C=sum(data4==0)
          if (mC>1) cohe[g]<-n11C/(mC*(mC-1)/2)
          isol[g]=n00C / (mC*(p1-mC))
      }
      
      # rescohe[[k]]=cohe
      # resisol[[k]]=isol
      # tab=table(parti_ref_reduced)
      # sizegp=as.vector(tab)
      # alphaC=sizegp*(rep(p,k)-sizegp)        # mC*(m-mC), nb pairs of objects in C and Cbar
      # betaC=sizegp*(sizegp-1)/2              # mC * (mC-1)/2 , nb pairs of objects in C
      # if (k>1) {
      #    rescoheP=c(rescoheP,sum(betaC*rescohe[[k]],na.rm=TRUE)/ sum(betaC))
      #   resisolP=c(resisolP,sum(alphaC*resisol[[k]],na.rm=TRUE) /sum(alphaC))
      # }
      return(list(cohe=cohe,isol=isol))
    } 
    
    # Trimmed mean------------------
    trimmean<-function(x,f) {
      # compute the f fraction trimmed mean of the values in x
      xord=sort(x)
      n=length(x)
      k=n*f
      g=floor(k)
      r=k-g
      tmean=  (sum(xord[(g+2):(n-g-1)]) + ( (1-r)*(xord[g+1]+xord[n-g]) ))/(n-2*k)
      return(tmean)
    } 
  
  numCores <- parallel::detectCores()
  #  numCores
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  b<-NULL
  if (case=="row") {
     res<- foreach(b=1:B,.combine=list,.multicombine=TRUE,.packages="ClustVarLV") %dopar%
      {
           boot_fx_row(b)
      }
   }
  if (case=="column") {
    res<- foreach(b=1:B,.combine=list,.multicombine=TRUE,.packages="ClustVarLV") %dopar%
      {
        boot_fx_col(b)
      }
  }
  parallel::stopCluster(cl)   
  
   totres=array(data=NA,dim=c(nmax,4,B))
   dimnames(totres)[[1]]=1:nmax
   dimnames(totres)[[2]]=c("Rand","ARI","CoheP","IsolP")
   dimnames(totres)[[3]]=paste("boot", 1:B,sep=".")
   for (b in 1:B) totres[,,b]=matrix(unlist(res[[b]]),nmax,4)
   
   totres[1,,]=NA
   if (nmax==p) totres[nmax,,]=NA
   tmean=matrix(NA,nmax,4)
   for (crit in 1:4) {
     lmat<-as.list(data.frame(t(totres[,crit,])))
     tmean[,crit]=unlist(lapply(lmat,trimmean,f=0.2))
   }
   # meanres=apply(totres,c(1,2),mean,na.rm=TRUE)
   # medres=apply(totres,c(1,2),median,na.rm=TRUE)
   # madres=apply(totres,c(1,2),mad,na.rm=TRUE)
   # robsigres=1.4826*madres
   
  
#   boxplot(resari,ylim=c(mini,maxi),xlab="nb clusters",ylab="Adjusted Rand Index")
   dev.new()
   matplot(1:nmax,tmean,ylim=c(0,1),type="l",col=1:4,lty=1,xlab="nb clusters",main="2O% trimmed mean of the criteria")
   legend("bottomleft",col=1:4, lty=1,legend=c("Rand Index","Adjusted Rand Index", "Cohesion partition","Isolation partition"))
   dev.new()
   par(mfrow=c(2,2))
   boxplot(t(totres[,1,]), main="Rand Index")
   boxplot(t(totres[,2,]), main="Adjusted Rand Index")
   boxplot(t(totres[,3,]), main="Cohesion partition")
   boxplot(t(totres[,4,]), main="Isolation partition")
   par(mfrow=c(1,1))
  
   return(totres)  
 }
 