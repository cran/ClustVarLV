#' @title Description of the clusters of variables
#' 
#' @description This function provides the list of the variables within each group and complementary informations.
#' Users will be asked to specify the number of clusters, 
#' 
#' @param resclv : result of CLV() or CLV_kmeans()
#' @param K : the number of clusters (unless if CLV_kmeans was used)
#' 
#'  
#' @details The ouputs include :
#' \itemize{
#' \item the size of the groups, \cr
#' \item the list of the variables within each group. FFor each cluster, the correlation of the each variable with its group latent component 
#'   and the correlation with the next neighbouring group latent component are given.  \cr
#' \item the proportion of the variance within each group explained by its latent variable, \cr
#'  \item the proportion of the whole dataset account by the group latent variables \cr
#'  \item the matrix of correlation between the latent variables. }
#' 
#' @examples 
#' data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = "directional", sX = TRUE)
#' summary_clv(resclvX, K = 4) 
#' @export
#' 

summary_clv <-
function(resclv,K=NULL)
  {
    if (!inherits(resclv, "clv"))   stop("non convenient objects")
  
    descrip_gp(resclv,K)
  }