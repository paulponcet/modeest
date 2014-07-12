\name{lientz}
\alias{lientz}
\alias{Lientz}
\alias{mlv.lientz}
\alias{plot.lientz}
\alias{print.lientz}

\title{ The Empirical Lientz Function and The Lientz Mode Estimator }
\description{
  The Lientz mode estimator is nothing but the value minimizing the empirical
  Lientz function. \cr
  A 'plot' and a 'print' methods are provided. 
}
\usage{
lientz(x, 
       bw = NULL)
  
  \method{mlv}{lientz}(x, 
    bw = NULL, 
    abc = FALSE, 
    par = shorth(x), 
    optim.method = "BFGS", 
    ...)
           
  \method{plot}{lientz}(x, 
     zoom = FALSE, 
     ...)
       
  \method{print}{lientz}(x, 
      digits = NULL, 
      ...)
}
\arguments{
  \item{x}{ numeric (vector of observations) or an object of class \code{"lientz"}. }
  \item{bw}{ numeric. The smoothing bandwidth to be used. 
  Should belong to (0, 1). Parameter 'beta' in Lientz (1970) function. }
  \item{abc}{ logical. If \code{FALSE} (the default), the Lientz empirical function 
  is minimised using \code{\link{optim}}. }
  \item{par}{ numeric. The initial value used in \code{\link{optim}}. }
  \item{optim.method}{ character. If \code{abc = FALSE}, the method used in \code{\link{optim}}. }  
  \item{zoom}{ logical. If \code{TRUE}, one can zoom on the graph created. }
  \item{digits}{ numeric. Number of digits to be printed. }
  \item{...}{ if \code{abc = FALSE}, further arguments to be passed to \code{\link{optim}}, or further arguments to be passed to \code{\link{plot.default}}.  }
}
\details{
The Lientz function is the smallest non-negative quantity \eqn{S(x,\beta)}{S(x,b)}, 
where \eqn{\beta}{b} = \code{bw}, such that 
\deqn{F(x+S(x,\beta)) - F(x-S(x,\beta)) \geq \beta.}{F(x+S(x,b)) - F(x-S(x,b)) >= b.}
Lientz (1970) provided a way to estimate \eqn{S(x,\beta)}{S(x,b)}; this estimate 
is what we call the empirical Lientz function. % a ajouter : le choix de beta : petit
}
\value{
\code{lientz} returns an object of class \code{c("lientz", "function")}; this is a 
function with additional attributes:
  \item{x}{ the \code{x} argument}
  \item{bw}{ the \code{bw} argument }
  \item{call}{ the call which produced the result }
  
\code{mlv.lientz} returns a numeric value, the mode estimate. If \code{abc = TRUE}, 
the \code{x} value minimizing the Lientz empirical function is returned. Otherwise, 
the \code{\link{optim}} method is used to perform minimization, and the  
attributes: 'value', 'counts', 'convergence' and 'message', coming from 
the \code{\link{optim}} method, are added to the result.
}
\note{
The user should preferentially call \code{mlv.lientz} through 
\code{mlv(x, method = "lientz", ...)}.
This returns an object of class \code{\link{mlv}}.
}
\references{
  \itemize{
    \item Lientz B.P. (1969).
    On estimating points of local maxima and minima of density functions.
    \emph{Nonparametric Techniques in Statistical Inference (ed. M.L. Puri, Cambridge University Press}, p.275-282.
     
    \item Lientz B.P. (1970).
    Results on nonparametric modal intervals.
    \emph{SIAM J. Appl. Math.}, \bold{19}:356-366.
    
    \item Lientz B.P. (1972).
    Properties of modal intervals.
    \emph{SIAM J. Appl. Math.}, \bold{23}:1-5.
  }
}
\author{ P. Poncet }
\seealso{ \code{\link{mlv}} for general mode estimation;
\code{\link{shorth}} for the shorth estimate of the mode }
\examples{
# Unimodal distribution
x <- rbeta(1000,23,4)

## True mode
betaMode(23, 4)

## Lientz object
f <- lientz(x, 0.2)
print(f)
plot(f)

## Estimate of the mode
mlv(f)              # optim(shorth(x), fn = f)
mlv(f, abc = TRUE) # x[which.min(f(x))]
M <- mlv(x, method = "lientz", bw = 0.2)
print(M)
plot(M)

# Bimodal distribution
x <- c(rnorm(1000,5,1), rnorm(1500, 22, 3))
f <- lientz(x, 0.1)
plot(f)
}
\keyword{ univar }