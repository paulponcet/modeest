\name{MS_max_representation}
\alias{MS_max_representation}
\alias{MS_max_representation}

\title{ MS_max_representation }
\description{
MS_max_representation
}
\usage{
  modetree(x,
           bw = NULL,
           kernel = "gaussian",
           par = shorth(x),
           iter = 1000,
           eps = 1e-08)
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. The smoothing bandwidth to be used. }
  \item{kernel}{ character. The kernel to be used. Available kernels are
    \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"},
    \code{"gaussian"}, \code{"optcosine"}, \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}.
    See \code{\link{density.default}} for more details on some of these kernels. }
  \item{par}{ numeric. The initial value used in the meanshift algorithm. }
  \item{iter}{ numeric. Maximal number of iterations. }
  \item{eps}{ numeric. Stopping criteria. }
}
\value{
\code{meanshift} returns a numeric value, the mode estimate,
with an attribute \code{"iterations"}.
The number of iterations can be less than \code{iter}
if the stopping criteria specified by \code{eps} is reached.
}
\note{
The user should preferentially call \code{meanshift} through
 \code{mlv(x, method = "meanshift", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{
    \item Burger, T. and Dhorne T. (2009). 
    A Graphical Tool for the Detection of Modes in Continuous Data
    \emph{??}, \bold{??}.
  }
}
\author{ T. Burger and T. Dhorne }
\seealso{ \code{\link{mlv}}, \code{\link{tsybakov}} }
\examples{
boulons <- read.table("boulons.txt", header = TRUE, sep = " ");
data <- boulons[,2];
MS_max_representation(data, tailleNoyauMax = 0.1)

#Il est aussi possible de simuler un jeu de données bimodale:

data <- c(rnorm(80, mean = 990, sd = 5),rnorm(80, mean = 1000, sd = 2));
MS_max_representation(data, tailleNoyauMax = 0.1)
}
\keyword{ univar }
