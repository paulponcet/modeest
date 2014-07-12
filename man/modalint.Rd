\name{modalint}
\alias{modalint}

\title{ Modal Interval }
\description{
Modal interval
}
\usage{
  modalint(x,
           bw = seq(0.01, 0.99, 0.01),
           plot.it = FALSE,
           tie.action = "mean",
           tie.limit = 0.05)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The sequence of smoothing bandwidths to be used. }
  \item{plot.it}{ boolean.  }
  \item{tie.action}{ character. The action to take if a tie is encountered. }
  \item{tie.limit}{ numeric. A limit deciding whether or not a warning is given when a tie is encountered. }
}
\value{
\code{modalint} returns
}
\note{
The user should preferentially call \code{meanshift} through
 \code{mlv(x, method = "meanshift", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}}, \code{\link{tsybakov}} }
\examples{
x <- rweibull(1000, shape = 12, scale = 0.8)
modalint(x, plot.it = TRUE)
}
\keyword{ univar }
