% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CLV3W.R
\name{CLV3W}
\alias{CLV3W}
\title{Hierarchical clustering of variables (associated with mode 2 three-way array)  with consolidation}
\usage{
CLV3W(X,mode.scale=0,NN=FALSE,moddendoinertie=TRUE,gmax=20,graph=TRUE,cp.rand=10)
}
\arguments{
\item{X}{: a three way array - variables of mode 2 will be clustered}

\item{mode.scale}{: scaling parameter applied to X, by default  centering of X (for mode 2 x mode 3) is done. By default no scaling (mode.scale=0) \cr
0 : no scaling only centering - the default \cr
1 : scaling with standard deviation of  (mode 2 x mode 3) elements \cr
2 : global scaling (each block i.e. each mode 2 slice will have the same inertia ) \cr
3 : global scaling (each block i.e. each mode 3 slice will have the same inertia )}

\item{NN}{: non Negativity constraint to be added on the loading coefficients. By default no constraint (NN=FALSE)   \cr
 TRUE : a non negativity constrained is applied on the loading coefficients to set them as positive values \cr
FALSE : loading coefficients may be either positive or negative}

\item{moddendoinertie}{: dendrogram. By default it is based on the delta clustering criterion (moddendoinertie =TRUE) \cr
 TRUE : dendrogram associated with the clustering criterion delta \cr
FALSE : dendrogram associated with the the height (cumulative delta)}

\item{gmax}{: maximum number of partitions for which the consolidation will be done (default : gmax=11)}

\item{graph}{: boolean, if TRUE, the graphs associated with the dendrogram and the evolution of the aggregation criterion are displayed (default : graph=TRUE)}

\item{cp.rand}{: number of random starts associated with the one rank Candecomp/Parafac model (By default cp.rand=10)}
}
\value{
\item{tabres}{ Results of the hierarchical clustering algorithm.
        In each line you find the results of one specific step of the hierarchical clustering.
        \itemize{
               \item {Columns 1 and 2}{ : the numbers of the two groups which are merged}
               \item {Column 3}{ : name of the new cluster}
               \item {Column 4}{ : the value of the aggregation criterion for the Hierarchical Ascendant Clustering (delta) : delta loss}
               \item {Column 5}{ : the loss value of the clustering criterion for the HAC}
               \item {Column 6}{ : the percentage of explained inertia of the data array X}
               \item {Column 7}{ : the loss value of the clustering criterion  after consolidation}
               \item {Column 8}{ : the percentage of explained inertia of the data array X after consolidation }
               \item {Column 9}{ : number of iterations in the partitioning algorithm. \cr
               Remark : A zero in columns 7 to 9 indicates that no consolidation was done}
       }}

\item{hclust}{ contains the results of the HCA }

\item{partition K}{ contains a list for each number of clusters of the partition, K=1 to gmax with
         \itemize{
               \item {clusters}{ :  in line 1, the groups membership before consolidation; in line 2 the groups membership after consolidation}
               \item {comp}{ : the latent components of the clusters associated with the first mode (after consolidation)}
               \item {loading}{ : the vector of loadings  associated with the second mode by cluster (after consolidation)}
               \item {weigth}{ : the vector of weights  associated with the third mode by cluster (after consolidation)}
               \item {criterion}{ : vector of loss giving for each cluster the residual amount between the sub-array and its reconstitution associated with the cluster one rank PARAFAC model (after consolidation)}
         }}

\item{param}{ contains the clustering parameters
         \itemize{
               \item {gmax}{ :  maximum number of partitions for which the consolidation has been done}
               \item {X}{ : the scaled three-way array}
         }}

call : call of the method
}
\description{
Hierarchical Cluster Analysis of a set of variables (mode 2) given a three-way array with a further consolidation step.
Each group of variables is associated with a one-rank PARAFAC model (comp x loading x weight).
Moreover, a Non Negativity (NN) constraint may be added to the model, so that the loading coefficients have positive values.
Return an object of class clv3w.
}
\examples{
data(ciders)
## Cluster Analysis of cider sensory descriptors with block scaling
## to set the assessors to the same footing
res.cider<-CLV3W(ciders,mode.scale=3,NN=FALSE,moddendoinertie=FALSE,gmax=20,graph=FALSE,cp.rand=5)
plot(res.cider,type="delta")
plot(res.cider,type="dendrogram")
print(res.cider)
summary(res.cider,2)
get_comp(res.cider,2)
get_loading(res.cider,2)
get_weight(res.cider,2)


}
\references{
Wilderjans, T. F., & Cariou, V. (2016). CLV3W: A clustering around latent variables approach to detect panel disagreement in three-way conventional sensory profiling data. Food quality and preference, 47, 45-53.

Cariou, V., & Wilderjans, T. F. (2018). Consumer segmentation in multi-attribute product evaluation by means of non-negatively constrained CLV3W. Food Quality and Preference, 67, 18-26.
}
\seealso{
CLV3W_kmeans, get_comp, get_loading, get_partition, plot, plot_var.clv3w,
}
\author{
Veronique Cariou, \email{veronique.cariou@oniris-nantes.fr}
}
