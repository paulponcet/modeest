\name{tsybakov}
%\alias{mlv.tsybakov}
\alias{tsybakov}
\alias{Tsybakov}
%\alias{mlv.dmp}
%\alias{Djeddour}
%\name{mlv.wtp}
%\alias{mlv.wtp}
%\alias{mlv.mizoguchi}
%\alias{Mizoguchi}
%\alias{mlv.shimura}
%\alias{Shimura}

\title{ The Tsybakov Mode Estimator }
\description{
  This mode estimator is based on a gradient-like recursive algorithm.
  It includes the Mizoguchi-Shimura (1976) mode estimator, based on the 
  window training procedure. 
}
\usage{
  tsybakov(x, 
           bw = NULL, 
           a, 
           alpha = 0.9, 
           kernel = "triangular", 
           dmp = TRUE, 
           par = shorth(x))
}
\arguments{
  \item{x}{ numeric. Vector of observations. }
  \item{bw}{ numeric. Vector of length \code{length(x)} giving the sequence of smoothing bandwidths to be used. }
  \item{a}{ numeric. Vector of length \code{length(x)} used in the gradient algorithm. }
  \item{alpha}{ numeric. An alternative way of specifying \code{a}. See 'Details'. }
  \item{kernel}{ character. The kernel to be used. Available kernels are 
    \code{"biweight"}, \code{"cosine"}, \code{"eddy"}, \code{"epanechnikov"}, 
    \code{"gaussian"}, \code{"optcosine"}, \code{"rectangular"}, \code{"triangular"}, \code{"uniform"}.
    See \code{\link{density.default}} for more details on some of these kernels. }
  \item{dmp}{ logical. If \code{TRUE}, Djeddour et al. version of the estimate is used. }
  \item{par}{ numeric. Initial value in the gradient algorithm. 
              Default value is \code{\link{shorth}(x)}. }
}
\details{
If \code{bw} is missing, then \code{bw = (1:length(x))^(-1/7)}, which is the default value advised by Djeddour et al (2003).
If \code{a} is missing, then \code{a = (1:length(x))^(-alpha)} (with \code{alpha = 0.9} if \code{alpha} is missing), which is the default value advised by Djeddour et al (2003).
}
\value{
A numeric value is returned, the mode estimate.
}
\note{
The user should preferentially call \code{tsybakov} through 
 \code{mlv(x, method = "tsybakov", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\section{Warning}{
The Tsybakov mode estimate as it is presently computed does not work very well. 
The reasons of this inefficiency are under investigation.
}
\references{
  \itemize{
    \item Mizoguchi R. and Shimura M. (1976).
    Nonparametric Learning Without a Teacher Based on Mode Estimation.
    \emph{IEEE Transactions on Computers}, \bold{C25}(11):1109-1117.
    
    \item Tsybakov A. (1990).
    Recursive estimation of the mode of a multivariate distribution.
    \emph{Probl. Inf. Transm.}, \bold{26}:31-37.
    
    \item Djeddour K., Mokkadem A. et Pelletier M. (2003).
    Sur l'estimation recursive du mode et de la valeur modale d'une densite de probabilite.
    \emph{Technical report 105}.
    
    \item Djeddour K., Mokkadem A. et Pelletier M. (2003).
    Application du principe de moyennisation a l'estimation recursive du mode et de la valeur modale d'une densite de probabilite.
    \emph{Technical report 106}.
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}} for general mode estimation }
\examples{
x <- rbeta(1000, shape1 = 2, shape2 = 5)

## True mode:
betaMode(shape1 = 2, shape2 = 5)

## Estimation:
tsybakov(x, kernel = "triangular")
tsybakov(x, kernel = "gaussian", alpha = 0.99)
M <- mlv(x, method = "tsybakov", kernel = "gaussian", alpha = 0.99)
print(M)
plot(M)
}
\keyword{ univar }
\keyword{ nonparametric }
