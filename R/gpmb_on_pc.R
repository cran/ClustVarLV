#' Representation of the variables and their group membership
#' 
#' Loading plot of the variables from a Principal Components Analysis. The group membership of the variables is superimposed.
#'
#' @param resclv results of CLV.r or CLV_kmeans.r
#' @param X the initial matrix
#' @param K the number of groups in the partition
#' @param axeh component number for the horizontal axis 
#' @param axev component number for the vertical axis 
#' @param v_colors default NULL
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = 1, sX = TRUE, graph = TRUE)
#' gpmb_on_pc(resclvX, X = apples_sh$senso, K = 4, axeh = 1, axev = 2)
gpmb_on_pc <-
function(resclv,X,K=NULL,axeh=1,axev=2,v_colors=NULL) {
  if (!inherits(resclv, c("clv","clvkmeans"))) 
    stop("non convenient objects")
  if(is.null(resclv$param$K)) { 
    if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
    clusters<-resclv[[K]]$clusters[2,]
} else {
    clusters<-resclv$clusters[2,] 
    K<-resclv$param$K
 }
  
  X<- scale(X, center=T, scale=resclv$param$sX)
  p <- dim(X)[2] 
  n <- dim(X)[1]
   
  if (is.null(v_colors)) {v_colors <- c("blue","red","green","black",
                                        "purple","orange","yellow","tomato","pink",
                                        "gold","cyan","turquoise","violet","green",
                                        "khaki","maroon","chocolate","burlywood")}
   # PCA of X
  A<-max(axeh,axev)
  if (n<p) {
     reseig<-eigen(X%*%t(X)/n)
     valp<-100*reseig$values[1:A]/sum(reseig$values)
     coordvar<-t(X)%*%reseig$vectors[,1:A]/sqrt(n)
  } else {
     reseig<-eigen(t(X)%*%X/(n))
     valp<-100*reseig$values[1:A]/sum(reseig$values)
     coordvar<-reseig$vectors[,1:A]%*%diag(sqrt(reseig$values[1:A]))
  }
 #re-orientation of the PC so that the maximal nb of var have positive coordinate alonf this axe
  for (a in 1:A) {
      if (sign(mean(coordvar[,a]))==(-1)) coordvar[,a]=coordvar[,a]*(-1)
  }
    
  
  dev.new() 
  par(pty="s")
  colpart<-NULL
  for (j in 1:p) {
    colpart[j]<-v_colors[clusters[j]]
  }
  plot(coordvar[,c(axeh,axev)],col=colpart,pch=20,cex.axis=0.7,cex.lab=0.7,
       xlab=paste("Dim ",axeh," (",round(valp[axeh],2),"%)"), 
       ylab=paste("Dim ",axev," (",round(valp[axev],2),"%)"))  #xlim=c(-1,1), ylim=c(-1,1),
  arrows(rep(0,p), rep(0,p), coordvar[,axeh], coordvar[,axev], col= colpart,angle=0)
  abline(h=0,v=0)
  #symbols(0, 0, circles = 1, inches = FALSE, add = TRUE) 
  legend("topleft",paste("G",1:K,sep=""),col=v_colors[1:K],title="Groupes",
         lty="solid",seg.len=0.6,cex=0.7,ncol=2)     
}
