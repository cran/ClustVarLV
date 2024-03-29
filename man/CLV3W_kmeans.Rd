% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CLV3W_kmeans.R
\name{CLV3W_kmeans}
\alias{CLV3W_kmeans}
\title{Partitioning algorithm of a set of variables (associated with mode 2) oh a  three-way array}
\usage{
CLV3W_kmeans(X,K,mode.scale=0,NN=FALSE,init=10,cp.rand=5)
}
\arguments{
\item{X}{: a three way array - variables of mode 2 will be clustered}

\item{K}{: number of clusters}

\item{mode.scale}{: scaling parameter applied to X, by default  centering of X (for mode 2 x mode 3) is done. By default no scaling (mode.scale=0) \cr
0 : no scaling only centering - the default \cr
1 : scaling with standard deviation of  (mode 2 x mode 3) elements \cr
2 : global scaling (each block i.e. each mode 2 slice will have the same inertia ) \cr
3 : global scaling (each block i.e. each mode 3 slice will have the same inertia )}

\item{NN}{: non Negativity constraint to be added on the loading coefficients. By default no constraint (NN=FALSE)   \cr
 TRUE : a non negativity constrained is applied on the loading coefficients to set them as positive values \cr
FALSE : loading coefficients may be either positive or negative}

\item{init}{: either the number of random starts i.e. partitions generated for the initialisation (By default init=10)}

\item{cp.rand}{: number of random starts associated with the one rank Candecomp/Parafac model (By default cp.rand=10)}
}
\value{
\item{results}{
         \itemize{
            \item {clusters}{:  in line 1, the groups membership in the initial partition; in line 2 the final groups membership }
            \item {comp}{: the latent components of the clusters associated with the first mode }
            \item {loading}{: the vector of loadings  associated with the second mode by cluster }
            \item {weigth}{: the vector of weights  associated with the third mode by cluster }
            \item {criterion}{: vector of loss giving for each cluster the residual amount between the sub-array and its reconstitution associated with the cluster one rank PARAFAC model}
            \item {niter}{: number of iterations of the partitioning alorithm}
         }}

\item{param}{contains the clustering parameters
         \itemize{
               \item {X}{: the scaled three-way array}
         }}

call : call of the method
}
\description{
Each group of variables is associated with a one-rank PARAFAC model (comp x loading x weight).
Moreover, a Non Negativity (NN) constraint may be added to the model, so that the loading coefficients have positive values.
Return an object of class clv3w.
}
\examples{
data(coffee)
## Cluster Analysis of coffee sensory descriptors with block scaling
## to set the assessors to the same footing
res.coffee <- CLV3W_kmeans(coffee,K=2,NN=TRUE,mode.scale=3,init=1,cp.rand=1)
summary(res.coffee)
get_partition(res.coffee)

}
\references{
Wilderjans, T. F., & Cariou, V. (2016). CLV3W: A clustering around latent variables approach to detect panel disagreement in three-way conventional sensory profiling data. Food quality and preference, 47, 45-53.

Cariou, V., & Wilderjans, T. F. (2018). Consumer segmentation in multi-attribute product evaluation by means of non-negatively constrained CLV3W. Food Quality and Preference, 67, 18-26.
}
\seealso{
summary.clv3W, print.clv3W
}
\author{
Veronique Cariou, \email{veronique.cariou@oniris-nantes.fr}
}
