\name{mlv}
\alias{mlv}
\alias{mlv.default}
\alias{mlv.factor}
\alias{mlv.integer}
\alias{mlv.character}
\alias{mlv.density}
\alias{plot.mlv}
\alias{print.mlv}
\alias{as.numeric.mlv}

\title{ Estimation of the Mode }
\description{
\code{mlv} is a generic function which enables to compute an estimate of the mode of a univariate distribution.
Many different estimates (or methods) are provided: 
\itemize{
  \item \code{\link{mfv}}, 
  which returns the most frequent value(s) in a given numerical vector, 
  \item the \code{\link{Lientz}} mode estimator, which is the value minimizing the Lientz function estimate, 
  \item the Chernoff mode estimator, also called \code{\link{naive}} mode estimator, which is defined as the 
  center of the interval of given length containing the most observations, 
  \item the \code{\link{Venter}} mode estimator, including the \code{\link{shorth}}, i.e. the midpoint of the modal interval, 
  \item the \code{\link{Grenander}} mode estimator, 
  \item the half sample mode (\code{\link{HSM}}) and the half range mode (\code{\link{HRM}}), which are iterative versions of the Venter mode estimator, 
  \item \code{\link{Parzen}}'s kernel mode estimator, which is the value maximizing the kernel density estimate, 
  \item the \code{\link{Tsybakov}} mode estimator, based on a gradient-like recursive algorithm, 
  \item the \code{\link{Asselin}} de Beauville mode estimator, based on a algorithm detecting chains and holes in the sample, 
  \item the \code{\link{Vieu}} mode estimator, 
  \item the \code{\link{meanshift}} mode estimator. 
}

\code{mlv} can also be used to compute the mode of a given distribution, with \code{mlv.character}.

A 'plot' and a 'print' methods are provided. 
}
\usage{
mlv(x, 
    ...)
  
  \method{mlv}{default}(x, 
    bw = NULL, 
    method, 
    na.rm = FALSE, 
    boot = FALSE, 
    R = 100, 
    B = length(x), 
    ...)
      
  \method{mlv}{factor}(x, 
    ...)
      
  \method{mlv}{integer}(x, 
    na.rm = FALSE, 
    ...)
      
  \method{mlv}{character}(x, 
    ...)
      
  \method{mlv}{density}(x, 
    all = TRUE, 
    abc = FALSE, 
    ...)
  
  \method{plot}{mlv}(x, 
     ...)
       
  \method{print}{mlv}(x, 
      digits = NULL, 
      ...)
        
  \method{as.numeric}{mlv}(x, 
           ...)
}
\arguments{
  \item{x}{ numeric (vector of observations), or an object of class \code{"factor"}, \code{"integer"}, etc. 
  For the function \code{as.numeric}, an object of class \code{"mlv"}. }
  \item{bw}{ numeric. The bandwidth to be used. This may have different meanings regarding the \code{method} used. }
  \item{method}{ character. One of the methods available for computing the mode estimate. See 'Details'. }
  \item{na.rm}{ logical. Should missing values be removed? }
  %\item{dip.test}{ logical. Should Hartigan's \code{\link{DIP}} statistic for unimodality be computed? }
  \item{boot}{ logical. Should bootstrap resampling be done? }
  \item{R}{ numeric. If \code{boot = TRUE}, the number of bootstrap resampling rounds to use.}
  \item{B}{ numeric. If \code{boot = TRUE}, 
  the size of the bootstrap samples drawn from \code{x}. 
  Default is to use a sample which is the same size as data. 
  For large data sets, this may be slow and unnecessary. }
  \item{all}{ logical. }
  \item{abc}{ logical. If \code{FALSE} (the default), the estimate of the density function  
  is maximised using \code{\link{optim}}. }
  \item{digits}{ numeric. Number of digits to be printed. }
  \item{...}{ Further arguments to be passed to the function called for computation. 
  This function is related to the \code{method} argument. }
}
\details{
For the function \code{mlv.default}, available methods are \code{"mfv"}, \code{"lientz"}, \code{"naive"}, \code{"venter"},
\code{"grenander"}, \code{"hsm"}, \code{"hrm"}, \code{"parzen"}, \code{"tsybakov"}, and \code{"asselin"}. 
See the description above and the associated links. \cr
If \code{x} is of class \code{"factor"} or \code{"integer"}, the most frequent value found in \code{x} is returned.\cr
If \code{x} is of class \code{"character"}, \code{x} should be one of \code{"beta"}, \code{"cauchy"}, \code{"gev"}, etc. 
i.e. a character for which a function \code{'x'Mode} exists (for instance \code{betaMode}, \code{cauchyMode}, etc.). 
See \code{\link{distribMode}} for the available functions. The mode of the corresponding distribution is returned.\cr
If \code{x} is of class \code{"density"}, the value where the density is maximised is returned.\cr
For the S3 function \code{mlv.lientz}, see \code{\link{Lientz}} for more details.
}
\value{
\code{mlv} returns an object of \code{\link{class}} \code{"mlv"}. \cr
An object of class \code{"mlv"} is a list containing at least the following components:

  \item{M}{the value of the mode}
  \item{skewness}{Bickel's measure of \code{\link{skewness}}}
  \item{x}{the argument \code{x}}
  \item{method}{the argument \code{method}}
  \item{bw}{the bandwidth}
  %\item{dip.stat}{Hartigan's DIP statistic}
  \item{boot}{the argument \code{boot}}
  \item{boot.M}{if \code{boot = TRUE}, the resampled values of the mode}
  \item{call}{the call which produced the result}
}
\references{
See the references on mode estimation on the \code{\link{modeest-package}}'s page.\cr
%For the DIP test, see:
% 
%  \itemize{
%    \item Hartigan J.A. and Hartigan P.M. (1985).
%    The Dip Test of Unimodality.
%    \emph{Ann. Statist.}, \bold{13}:70-84.
%    
%    \item Hartigan P.M. (1985).
%    Computation of the Dip Statistic to Test for Unimodality.
%    \emph{Appl. Statist. (JRSS C)}, \bold{34}:320-325.
%  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mfv}},
          \code{\link{Lientz}},
          \code{\link{naive}},
          \code{\link{venter}},
          \code{\link{grenander}},
          \code{\link{hrm}},
          \code{\link{hsm}},
          \code{\link{parzen}},
          \code{\link{tsybakov}}, 
          \code{\link{skewness}} }
\examples{
# Unimodal distribution
x <- rbeta(1000,23,4)

## True mode
betaMode(23, 4)
# or
mlv("beta", 23, 4)

## Estimate of the mode
mlv(x, method = "lientz", bw = 0.2)
mlv(x, method = "naive", bw = 1/3)
mlv(x, method = "venter", type = "shorth")
mlv(x, method = "grenander", p = 4)
mlv(x, method = "hrm", bw = 0.3)
mlv(x, method = "hsm")
mlv(x, method = "parzen", kernel = "gaussian")
mlv(x, method = "tsybakov", kernel = "gaussian")
mlv(x, method = "asselin", bw = 2/3)
mlv(x, method = "vieu")
mlv(x, method = "meanshift")

## Bootstrap
M <- mlv(x, method = "kernel", boot = TRUE, R = 150)
print(M)
plot(M)
print(mean(M[["boot.M"]]))
}
\keyword{ univar }
\keyword{ robust }
\keyword{ nonparametric }
