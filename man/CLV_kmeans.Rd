\name{CLV_kmeans}
\alias{CLV_kmeans}
\title{K-means algorithm for the clustering of variables}
\usage{
  CLV_kmeans(X, Xu = NULL, Xr = NULL, method, sX = TRUE,
    sXr = FALSE, sXu = FALSE, init, iter.max = 20,
    nstart = 1)
}
\arguments{
  \item{X}{The matrix of the variables to be clustered}

  \item{Xu}{The external variables associated with the
  columns of X}

  \item{Xr}{The external variables associated with the rows
  of X}

  \item{method}{The criterion to use in the cluster
  analysis.\cr 1 : the squared covariance is used as a
  measure of proximity (directional groups).\cr 2 : the
  covariance is used as a measure of proximity (local
  groups)}

  \item{sX}{TRUE/FALSE : standardization or not of the
  columns X (TRUE by default)\cr (predefined -> cX = TRUE :
  column-centering of X)}

  \item{sXr}{TRUE/FALSE : standardization or not of the
  columns Xr (FALSE by default)\cr (predefined -> cXr =
  TRUE : column-centering of Xr)}

  \item{sXu}{TRUE/FALSE : standardization or not of the
  columns Xu (FALSE by default)\cr (predefined -> cXu=
  FALSE : no centering, Xu considered as a weight matrix)}

  \item{init}{a number i.e.  the size of the partition, K,
  or a vector of INTEGERS i.e. the group membership of each
  var in the initial parition (integer between 1 and K)}

  \item{iter.max}{maximal number of iteration for the
  consolidation (20 by default)}

  \item{nstart}{nb of random initialisations in the case of
  init = a number (1 by default)}
}
\value{
  \item{tabres}{ The value of the clustering criterion at
  convergence.\cr The percentage of the explained initial
  criterion value.\cr The number of iterations in the
  partitioning algorithm.} \item{clusters}{ the group's
  membership} \item{comp}{ The latent components of the
  clusters} \item{loading}{ if there are external variables
  Xr or Xu : The loadings of the external variables}
}
\description{
  K-means algorithm for the clustering of variables.
  Directional or local groups may be defined. Each group of
  variables is associated with a latent component. Moreover
  external information collected on the observations or on
  the variables may be introduced.
}
\details{
  The initalization can be made at random, repetitively, or
  can be defined by the user.
}
\examples{
data(apples_sh)
#local groups with external variables Xr
resclvkmYX <- CLV_kmeans(X = apples_sh$pref, Xr = apples_sh$senso,
                         method = 2, sX = FALSE, sXr = TRUE, init = 2, nstart = 20)
}
\seealso{
  CLV, LCLV
}

