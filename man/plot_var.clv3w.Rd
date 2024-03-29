% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_var.clv3w.R
\name{plot_var.clv3w}
\alias{plot_var.clv3w}
\title{Scores plot from a Candecomp Parafac analysis. The group membership of the variables is superimposed.}
\usage{
plot_var.clv3w(
  resclv3w,
  K = NULL,
  axeh = 1,
  axev = 2,
  labels = FALSE,
  cex.lab = 1,
  v_colors = NULL,
  v_symbol = FALSE,
  beside = FALSE,
  mode3 = FALSE
)
}
\arguments{
\item{resclv3w}{the data matrix}

\item{K}{the number of groups in the partition (already defined if CLV3W_kmeans i used)}

\item{axeh}{component number for the horizontal axis}

\item{axev}{component number for the vertical axis}

\item{labels}{boolean to add variable' labels (label=TRUE) on the plot or not (label=FALSE). By default label=TRUE}

\item{cex.lab}{magnification to be used for labels (1 by default)}

\item{v_colors}{default NULL. If missing colors are given, by default}

\item{v_symbol}{symbols are given instead of colors for the identification of the groups/ =FALSE: no symbol (by default).}

\item{beside}{plot per cluster of variables, side-by-side\cr
=FALSE : an unique plot with all the variables with the identification of their group membership (by default).}

\item{mode3}{projection of the mode 3 elements onto the scores plot\cr
=FALSE : mode 3 elements are not represented (by default).}
}
\description{
Scores plot from a Candecomp Parafac analysis. The group membership of the variables is superimposed.
}
