\name{gpmb_on_pc}
\alias{gpmb_on_pc}
\title{Representation of the variables and their group membership}
\usage{
  gpmb_on_pc(resclv, X, K = NULL, axeh = 1, axev = 2,
    v_colors = NULL)
}
\arguments{
  \item{resclv}{results of CLV.r or CLV_kmeans.r}

  \item{X}{the initial matrix}

  \item{K}{the number of groups in the partition}

  \item{axeh}{component number for the horizontal axis}

  \item{axev}{component number for the vertical axis}

  \item{v_colors}{default NULL}
}
\description{
  Loading plot of the variables from a Principal Components
  Analysis. The group membership of the variables is
  superimposed.
}
\examples{
data(apples_sh)
resclvX <- CLV(X = apples_sh$senso, method = 1, sX = TRUE, graph = TRUE)
gpmb_on_pc(resclvX, X = apples_sh$senso, K = 4, axeh = 1, axev = 2)
}

