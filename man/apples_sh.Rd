\name{apples_sh}
\alias{apples_sh}
\docType{data}
\title{
apples from southern hemisphere data set
}
\description{
Sensory characterization and consumers preference for 12 variables of apples
}
\usage{
data(apples_sh)
}
\format{
  A data frame with 12 observations and 2 blocks of variables.
  \describe{
    \item{senso}{43 sensory attributes}
    \item{pref}{hedonic scores given by a panel of 60 consumers}
  }
}


\references{
Daillant-Spinnler, B, MacFie, H.J.H, Beyts, P.K., Hedderley, D. (1996). Relationships between perceived sensory properties and major preference directions of 12 variaties of apples from the southern hemisphere. Food Quality and Preference, 7(2), 113-126.

}
\examples{
  data(apples_sh)
  names(apples_sh)
  apples_sh$senso
  apples_sh$pref
}
\keyword{datasets}
