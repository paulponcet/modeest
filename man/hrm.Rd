\name{hrm}
%\alias{mlv.hrm}
\alias{hrm}
\alias{HRM}
\alias{half.range.mode}
\alias{bickel}
\alias{Bickel}
%\alias{mlv.bickel}

\title{ Half-range Mode }
\description{
This function computes Bickel's half range mode estimator described in Bickel (2002).
}
\usage{
  hrm(x, 
      bw = NULL, 
      ...)
}
\arguments{
  \item{x}{ numeric. Vector of observations }
  \item{bw}{ numeric. The bandwidth to be used. Should belong to (0, 1]. 
  This gives the fraction of the observations to consider at each step of the iterative algorithm. }
  \item{...}{ further arguments. }
}
\details{
The mode estimator is computed by iteratively identifying densest half ranges. 
A densest half range is an interval whose width equals half the current range, and
which contains the maximal number of observations. 
The subset of observations falling in the selected densest half range is then used to compute
a new range, and the procedure is iterated.
}
\value{
A numeric value is returned, the mode estimate.
}
\note{
The user should preferentially call \code{hrm} through 
 \code{mlv(x, method = "hrm", bw)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{ 
    \item Bickel D.R. (2002).
    Robust estimators of the mode and skewness of continuous data.
    \emph{Computational Statistics and Data Analysis}, \bold{39}:153-163.
    
    \item Hedges S.B. and Shah P. (2003).
    Comparison of mode estimation methods and application in molecular clock analysis.
    \emph{BMC Bioinformatics}, \bold{4}:31-41.
    
    \item Bickel D.R. and Fruehwirth R. (2006).
    On a Fast, Robust Estimator of the Mode: Comparisons to Other Robust Estimators with Applications.
    \emph{Computational Statistics and Data Analysis}, \bold{50}(12):3500-3530.
  }
}
\author{ The C and R code are due to Richard Bourgon \email{bourgon@stat.berkeley.edu}, see package \pkg{genefilter}. 
The algorithm is described in Bickel (2002). }
\seealso{ \code{\link{mlv}} for general mode estimation;
\code{\link{hsm}} for the half sample mode; 
\code{\link{venter}} for the Venter mode estimate }
\examples{
# Unimodal distribution
x <- rgamma(1000, shape = 31.9)

## True mode
gammaMode(shape = 31.9)

## Estimate of the mode
hrm(x, bw = 0.4)
M <- mlv(x, method = "hrm", bw = 0.4)
print(M)
plot(M)
}
\keyword{ univar }
\keyword{ robust }
\keyword{ nonparametric }
