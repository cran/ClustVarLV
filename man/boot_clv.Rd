% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot_clv.R
\name{boot_clv}
\alias{boot_clv}
\title{Boostrapping for assessing the stability of a CLV result}
\usage{
boot_clv(object, B = 100, nmax = NULL)
}
\arguments{
\item{object}{: result of CLV()}

\item{B}{: the number of bootstrap to be run (100 by default)}

\item{nmax}{: maximal size of the partitions to be considered (if NULL, the value of nmax used for the object is used)}
}
\value{
\item{matARI}{a matrix of the Adjusted Rand Index of size (B x nmax).}
}
\description{
Bootstrapping on the samples is performed. 
Each boostrapped data matrix is submitted to CLV in order to get partitions from 1 to nmax clusters.
For each number of clusters, K, the adjusted Rand Index between actual and the bootstrapped partitions are computed and used in order to assess the stability of the solution into K clusters.
Parallel computing is performed in order to save time.
}
\seealso{
CLV
}