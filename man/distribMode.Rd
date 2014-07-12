\name{distribMode}
\alias{distribMode}
\alias{betaMode}
\alias{cauchyMode}
\alias{chisqMode}
\alias{expMode}
\alias{fMode}
\alias{fiskMode}
\alias{gammaMode}
\alias{normMode}
\alias{frechetMode}
\alias{gevMode}
\alias{ghMode}
\alias{gompertzMode}
\alias{gpdMode}
\alias{gumbelMode}
\alias{hypMode}
\alias{koenkerMode}
\alias{kumarMode}
\alias{laplaceMode}
\alias{logisMode}
\alias{lnormMode}
\alias{nigMode}
\alias{paralogisticMode}
\alias{paretoMode}
\alias{rayleighMode}
\alias{stableMode}
%\alias{symstbMode}
\alias{rweibullMode}
\alias{tMode}
\alias{unifMode}
\alias{weibullMode}
\alias{bernMode}
\alias{binomMode}
\alias{geomMode}
\alias{hyperMode}
\alias{nbinomMode}
\alias{poisMode}

\title{ Computing the Mode of Some Distributions }
\description{
These functions return the mode of the main probability distributions implemented in \R.
}

\usage{
  ## Continuous distributions
  
  # Beta
  betaMode(shape1, shape2, ncp = 0) 
  
  # Cauchy
  cauchyMode(location = 0, ...) 
  
  # Chisquare
  chisqMode(df, ncp = 0) 
  
  # Exponentiel
  expMode(...) 
  
  # F
  fMode(df1, df2) 
  
  # Fisk (package 'VGAM')
  fiskMode(shape1.a, scale = 1)
  
  # Frechet (package 'evd')
  frechetMode(loc = 0, scale = 1, shape = 1, ...) 
  
  # Gamma
  gammaMode(shape, rate = 1, scale = 1/rate) 
  
  # Normal (Gaussian)
  normMode(mean = 0, ...) 
  
  # Generalised Extreme Value (package 'evd')
  gevMode(loc = 0, scale = 1, shape = 0, ...) 
  
  # Generalised Hyperbolic (package 'fBasics')
  ghMode(alpha = 1, beta = 0, delta = 1, mu = 0, 
         lambda = -1/2, ...) 

  # Gompertz (package 'VGAM')
  gompertzMode(shape, scale = 1)
  
  # Generalised Pareto (package 'evd')
  gpdMode(loc = 0, scale = 1, shape = 0, ...) 
  
  # Gumbel (package 'evd')
  gumbelMode(loc = 0, ...) 
  
  # Hyperbolic (package 'fBasics')
  hypMode(alpha = 1, beta = 0, delta = 1, mu = 0, 
          pm = c(1, 2, 3, 4)) 
  
  # Koenker (package 'VGAM')
  koenkerMode(location = 0, ...)

  # Kumaraswamy (package 'VGAM')
  kumarMode(shape1, shape2)
  
  # Laplace (package 'VGAM')
  laplaceMode(location = 0, ...)
  
  # Logistic
  logisMode(location = 0, ...) 
  
  # Lognormal
  lnormMode(meanlog = 0, sdlog = 1) 
  
  # Normal Inverse Gaussian (package 'fBasics')
  nigMode(alpha = 1, beta = 0, delta = 1, mu = 0, ...) 
  
  # Stable (package 'stabledist')
  stableMode(alpha, beta, gamma = 1, delta = 0, pm = 0, ...) 
  
  # Negative Weibull (package 'evd')
  rweibullMode(loc = 0, scale = 1, shape = 1, ...) 

  # Paralogistic (package 'VGAM')
  paralogisticMode(shape1.a, scale = 1)
  
  # Pareto (package 'VGAM')
  paretoMode(location, ...)
  
  # Rayleigh (package 'VGAM')
  rayleighMode(scale = 1)

  # T (Student)
  tMode(df, ncp = 0) 
  
  # Uniform
  unifMode(min = 0, max = 1)

  # Weibull
  weibullMode(shape, scale = 1, ...)
  
  ## Discrete distributions
  
  # Bernoulli
  bernMode(prob) 
  
  # Binomial
  binomMode(size, prob) 
  
  # Geometric
  geomMode(...) 
  
  # Hypergeometric
  hyperMode(m, n, k, ...) 
  
  # Negative Binomial
  nbinomMode(size, prob, mu) 
  
  # Poisson
  poisMode(lambda) 
}
\arguments{
  \item{shape1}{ First positive parameter of the Beta and Kumaraswamy distributions. See the package \pkg{VGAM} for more details. }
  \item{shape2}{ Second positive parameter of the Beta and Kumaraswamy distributions. See the package \pkg{VGAM} for more details. }
  \item{shape1.a}{ Shape parameter of the Fisk and Paralogistic distributions. See the package \pkg{VGAM} for more details. }
  \item{ncp}{ Non-centrality parameter of the Beta, Chisquare, and Student distributions. }
  \item{location}{ Location parameter of the Cauchy, Koenker, Laplace, Logistic, and Pareto distributions. See the package \pkg{VGAM} for more details. }
  \item{df}{ Degree of freedom of the Chisquare and Student distributions. }
  \item{df1}{ First degree of freedom of the F distribution. }
  \item{df2}{ Second degree of freedom of the F distribution. }
  \item{loc}{ Location parameter of the Fr\'{e}chet, Generalized Extreme Value, Generalized Pareto, Gumbel, and Negative Weibull distributions. }
  \item{scale}{ Scale parameter of the Fisk, Fr\'{e}chet, Gamma, Generalized Extreme Value, Gompertz, Generalized Pareto, Negative Weibull, Paralogistic, Rayleigh, and Weibull distributions. See the packages \pkg{evd} and \pkg{VGAM} for more details. }
  \item{shape}{ Shape parameter of the Fr\'{e}chet, Gamma, Generalized Extreme Value, Gompertz, Generalized Pareto, Negative Weibull, and Weibull distributions. See the packages \pkg{evd} and \pkg{VGAM} for more details. }
  \item{rate}{ An alternative way to specify the scale of the Gamma distribution. }
  \item{mean}{ Mean of the Normal distribution. }
  \item{alpha}{ Parameter of the Hyperbolic, Generalised Hyperbolic, Stable, and Normal Inverse Gaussian distributions. See the packages \pkg{fBasics} and \pkg{stabledist} for more details. }
  \item{beta}{ Parameter of the Hyperbolic, Generalised Hyperbolic, Stable, and Normal Inverse Gaussian distributions. See the packages \pkg{fBasics} and \pkg{stabledist} for more details. }
  \item{delta}{ Parameter of the Hyperbolic, Generalised Hyperbolic, Stable, and Normal Inverse Gaussian distributions. See the packages \pkg{fBasics} and \pkg{stabledist} for more details. }
  \item{mu}{ Parameter of the Hyperbolic, Generalised Hyperbolic, Normal Inverse Gaussian, and Negative binomial distributions. See the package \pkg{fBasics} for more details. }
  \item{lambda}{ Vector of (non-negative) means of the Poisson distribution. }
  \item{pm}{ Integer value for the selection of the parameterization of the Hyperbolic and Stable distributions. See the packages \pkg{fBasics} and \pkg{stabledist} for more details. }
  \item{meanlog}{ Mean of the Lognormal distribution on the log scale. }
  \item{sdlog}{ Standard deviation of the Lognormal distribution on the log scale. }
  \item{gamma}{ Scale parameter of the Stable distribution. See the package \pkg{stabledist} for more details. }
  \item{min}{ Lower limit of the Uniform distribution. Must be finite. }
  \item{max}{ Upper limit of the Uniform distribution. Must be finite. }
  \item{prob}{ Probability of success on each trial (between 0 and 1), used in the Bernoulli, Binomial and Negative Binomial distributions. }
  \item{size}{ Number of trials (zero or more), used in the Binomial and Negative Binomial distributions. }
  \item{m}{ Number of white balls in the urn for the Hypergeometric distribution. }
  \item{n}{ Number of black balls in the urn for the Hypergeometric distribution. }
  \item{k}{ Number of balls drawn from the urn for the Hypergeometric distribution. }
  \item{...}{ Further arguments, which will be ignored. }
}
\value{
A numeric value is returned, the (true) mode of the distribution.
}
\author{ P. Poncet, 
except for \code{hypMode} written by David Scott, \code{nigMode} and \code{stableMode} written by Diethelm Wuertz; see packages \pkg{fBasics} and \pkg{stabledist}. }
\note{
Some functions like \code{normMode} or \code{cauchyMode}, which are related 
to symmetric distributions, are trivial, but are implemented for exhaustivity.
}
\seealso{ \code{\link{mlv}} for the estimation of the mode; 
the documentation of the related distributions \code{\link{Beta}}, \code{\link{GammaDist}}, etc. }
\examples{
layout(mat = matrix(1:2,1,2))

## Beta distribution
curve(dbeta(x, shape1 = 2, shape2 = 3.1), xlim = c(0,1), ylab = "Beta density")
M <- betaMode(shape1 = 2, shape2 = 3.1)
abline(v = M, col = 2)
mlv("beta", shape1 = 2, shape2 = 3.1)
 
## Lognormal distribution          
curve(dlnorm(x, meanlog = 3, sdlog = 1.1), xlim = c(0, 10), ylab = "Lognormal density")
M <- lnormMode(meanlog = 3, sdlog = 1.1)
abline(v = M, col = 2)
mlv("lnorm", meanlog = 3, sdlog = 1.1) 

## Poisson distribution
poisMode(lambda = 6)
poisMode(lambda = 6.1)
mlv("poisson", lambda = 6.1)

layout(mat = matrix(1,1,1)) 
}
\keyword{ univar }
\keyword{ distribution }
