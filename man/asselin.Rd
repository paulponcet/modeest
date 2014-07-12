\name{asselin}
%\alias{mlv.asselin}
\alias{asselin}
\alias{Asselin}

\title{ The Asselin de Beauville Mode Estimator }
\description{
  This mode estimator is based on the algorithm described in Asselin de Beauville (1978).
}
\usage{
  asselin(x, 
          bw = NULL, 
          ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. A number in \code{(0, 1]}. If \code{bw = 1}, the selected 'modal chain' may be too long. }
  \item{...}{ further arguments to be passed to the \code{\link{quantile}} function. }
}
\value{
A numeric value is returned, the mode estimate.
}
\note{
The user should preferentially call \code{asselin} through 
 \code{mlv(x, method = "asselin", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{
    \item Asselin de Beauville J.-P. (1978).
    Estimation non parametrique de la densite et du mode, exemple de la distribution Gamma.
    \emph{Revue de Statistique Appliquee}, \bold{26}(3):47-70.
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}} for general mode estimation }
\examples{
x <- rbeta(1000, shape1 = 2, shape2 = 5)

## True mode:
betaMode(shape1 = 2, shape2 = 5)

## Estimation:
asselin(x, bw = 1)
asselin(x, bw = 1/2)
M <- mlv(x, method = "asselin")
print(M)
plot(M)
}
\keyword{ univar }
\keyword{ nonparametric }
