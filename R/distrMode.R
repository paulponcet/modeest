#' @title 
#' Mode of some continuous and discrete distributions
#' 
#' @description 
#' These functions return the mode of the main probability 
#' distributions implemented in R.
#' 
#' @note 
#' Some functions like \code{normMode} or \code{cauchyMode}, which relate 
#' to symmetric distributions, are trivial, but are implemented for the sake of 
#' exhaustivity.
#' 
#' @param x
#' character. The name of the distribution to consider. 
#' 
#' @param ...
#' Additional parameters. 
#' 
#' @return 
#' A numeric value is returned, the (true) mode of the distribution.
#' 
#' @author 
#' \code{\link[fBasics]{ghMode}} and \code{\link[fBasics]{ghtMode}} are from 
#' package \pkg{fBasics}; 
#' \code{\link[fBasics]{hypMode}} was written by David Scott; 
#' \code{\link[fBasics]{gldMode}}, \code{\link[fBasics]{nigMode}} and 
#' \code{\link[stabledist]{stableMode}} were written by Diethelm Wuertz.
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}} for the estimation of the mode; 
#' the documentation of the related distributions 
#' \code{\link[stats]{Beta}}, \code{\link[stats]{GammaDist}}, etc.
#' 
#' @importFrom statip name2distr
#' @export
#' 
distrMode <- 
function(x, 
         ...)
{
  stopifnot(is.character(x))
  x <- match.arg(statip::name2distr(x), .distributionsList())
  do.call(paste0(x, "Mode"), list(...))
}


# Beta distribution
#' @inheritParams stats::dbeta
#' 
#' @export
#' @rdname distrMode
#' 
#' @examples 
#' ## Beta distribution
#' curve(dbeta(x, shape1 = 2, shape2 = 3.1), 
#'       xlim = c(0,1), ylab = "Beta density")
#' M <- betaMode(shape1 = 2, shape2 = 3.1)
#' abline(v = M, col = 2)
#' mlv("beta", shape1 = 2, shape2 = 3.1)
#' 
betaMode <-
function(shape1,
         shape2,
         ncp = 0)
{
  if (ncp == 0) {
    M <- (shape1-1)/(shape1+shape2-2)
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  return(M)
}


# Cauchy distribution
#' @inheritParams stats::dcauchy
#' 
#' @export
#' @rdname distrMode
#' 
cauchyMode <-
function(location = 0,
         ...)
{
  return(location)
}


# Chi-square distribution
#' @inheritParams stats::dchisq
#' 
#' @export
#' @rdname distrMode
#' 
chisqMode <-
function(df,
         ncp = 0)
{
  if (ncp == 0) {
    M <- max(df-2, 0)
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  M
}


# Dagum distribution
#' @inheritParams VGAM::ddagum
#' 
#' @export
#' @rdname distrMode
#' 
dagumMode <- 
function(scale = 1, shape1.a, shape2.p)
{
  scale*((shape1.a*shape2.p-1)/(shape1.a+1))^(1/shape1.a)
}


# Exponential distribution
#' @inheritParams stats::dexp
#' 
#' @export
#' @rdname distrMode
#' 
expMode <-
function(...)
{
  0
}


# F distribution
#' @inheritParams stats::df
#' 
#' @export
#' @rdname distrMode
#' 
fMode <-
function(df1,
         df2)
{ 
  if (df1 > 2) {
    M <- (1-2/df1)*(df2/(2+df2))
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  M
}


# Fisk distribution
#' @inheritParams VGAM::dfisk
#' 
#' @export
#' @rdname distrMode
#' 
fiskMode <- 
function(scale = 1, 
         shape1.a)
{
  scale*((shape1.a-1)/(shape1.a+1))^(1/shape1.a)
}


# Frechet distribution
#' @inheritParams VGAM::dfrechet
#' 
#' @export
#' @rdname distrMode
#' 
frechetMode <-
function(location = 0,
         scale = 1,
         shape = 1,
         ...)
{
  location + scale*(shape/(1+shape))^(1/shape)
}


# Gamma distribution
#' @inheritParams stats::dgamma
#' 
#' @export
#' @rdname distrMode
#' 
gammaMode <-
function(shape,
         rate = 1,
         scale = 1/rate)
{
  scale*(shape-1)
}


# Gaussian (normal) distribution
#' @inheritParams stats::dnorm
#' 
#' @export
#' @rdname distrMode
#' 
normMode <-
function(mean = 0,
         ...)
{
  return(mean)
}


# Generalized extreme value distribution
#' @inheritParams VGAM::dgev
#' 
#' @export
#' @rdname distrMode
#' 
gevMode <-
function(location = 0,
         scale = 1,
         shape = 0,
         ...)
{
  k <- pmax(0,(1+shape))^(-shape)-1
  shape[shape==0] <- Inf
  location + (scale/shape)*k
}


# Generalized hyperbolic distribution
#' @inheritParams fBasics::ghMode
#' 
#' @importFrom fBasics ghMode
#' @export
#' @rdname distrMode
#' 
ghMode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0,
         lambda = -1/2)
{
  fBasics::ghMode(alpha, beta, delta, mu, lambda)
}


# Generalized Hyperbolic Student-t
#' @inheritParams fBasics::ghtMode
#' 
#' @importFrom fBasics ghtMode
#' @export
#' @rdname distrMode
#' 
ghtMode <-
function(beta = 0.1, 
         delta = 1, 
         mu = 0, 
         nu = 10) 
{
  fBasics::ghtMode(beta, delta, mu, nu)
}


# Generalized Lambda Distribution
#' @inheritParams fBasics::gldMode
#' 
#' @importFrom fBasics gldMode
#' @export
#' @rdname distrMode
#' 
gldMode <-
function(lambda1 = 0, 
         lambda2 = -1, 
         lambda3 = -1/8, 
         lambda4 = -1/8) 
{
  fBasics::gldMode(lambda1, lambda2, lambda3, lambda4)
}


# Gompertz distribution
#' @inheritParams VGAM::dgompertz
#' 
#' @export
#' @rdname distrMode
#' 
gompertzMode <-
function(scale = 1, 
         shape)
{
  if (shape < scale) {
    M <- log(scale/shape)/scale
  } else {
    M <- 0
  }
  M
}


# Generalized Pareto distribution
#' @inheritParams VGAM::dgpd
#' 
#' @export
#' @rdname distrMode
#' 
gpdMode <-
function(location = 0,
         scale = 1,
         shape = 0)
{
  if (shape == -1) {
    warning("all values between 'loc' and 'loc+scale' are modes, only the mean value is returned")
    M <- location + scale/2  
  } else if (-2-1/shape > 0) {
    M <- location - scale/shape
  } else {
    warning("the density is not continuous at the mode.")
    M <- location
  }
  M
}


# Gumbel distribution
#' @inheritParams VGAM::dgumbel
#' 
#' @export
#' @rdname distrMode
#' 
gumbelMode <-
function(location = 0,
         ...)
{
  location
}


# Hyperbolic distribution
#' @inheritParams fBasics::hypMode
#' 
#' @importFrom fBasics hypMode
#' @export
#' @rdname distrMode
#' 
hypMode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0,
         pm = c(1, 2, 3, 4)) 
{
  fBasics::hypMode(alpha, beta, delta, mu, pm)
}


# Koenker distribution 
#' @inheritParams stats::dcauchy
#'
#' @export
#' @rdname distrMode
#' 
koenkerMode <- 
function(location = 0, 
         ...)
{
  location
}


# Kumaraswamy distribution
#' @inheritParams VGAM::dkumar
#' 
#' @export
#' @rdname distrMode
#' 
kumarMode <-
function(shape1,
         shape2)
{
  (shape1-1)/(shape1*shape2 - 1)^(1/shape1)
}


# Laplace distribution
#' @inheritParams VGAM::dlaplace
#' 
#' @export
#' @rdname distrMode
#' 
laplaceMode <-
function(location = 0,
         ...)
{
  return(location)
}


# Logistic distribution
#' @inheritParams stats::dlogis
#' 
#' @export
#' @rdname distrMode
#' 
logisMode <-
function(location = 0,
         ...)
{
  location
}


# Lognormal distribution
#' @inheritParams stats::dlnorm
#' 
#' @export
#' @rdname distrMode
#' 
#' @examples 
#' ## Lognormal distribution
#' curve(stats::dlnorm(x, meanlog = 3, sdlog = 1.1), 
#'       xlim = c(0, 10), ylab = "Lognormal density")
#' M <- lnormMode(meanlog = 3, sdlog = 1.1)
#' abline(v = M, col = 2)
#' mlv("lnorm", meanlog = 3, sdlog = 1.1)
#' 
lnormMode <-
function(meanlog = 0,
         sdlog = 1)
{
  exp(meanlog - sdlog^2)
}


# Lomax distribution
#' @inheritParams VGAM::dlomax
#' 
#' @export
#' @rdname distrMode
#' 
lomaxMode <- 
function(...)
{
  0
}


# Maxwell-Boltzmann distribution
#' @inheritParams VGAM::dmaxwell
#' 
#' @export
#' @rdname distrMode
#' 
maxwellMode <-
function(rate)
{
  sqrt(2/rate)
}


# Multivariate normal distribution
#' @inheritParams mvtnorm::dmvnorm
#' 
#' @export
#' @rdname distrMode
#' 
mvnormMode <- 
function(mean, 
         ...)
{
  mean
}


# Nakagami distribution
#' @inheritParams VGAM::dnaka
#' 
#' @export
#' @rdname distrMode
#'
nakaMode <- 
function(scale = 1, 
         shape)
{
  sqrt(2*(2*shape - 1)*scale/(4*shape))
}


# Normal Inverse Gaussian distribution
#' @inheritParams fBasics::nigMode
#' 
#' @importFrom fBasics nigMode
#' @export
#' @rdname distrMode
#' 
nigMode <-
function(alpha = 1,
         beta = 0,
         delta = 1,
         mu = 0)
{
  fBasics::nigMode(alpha, beta, delta, mu)
}


# Paralogistic distribution
#' @inheritParams VGAM::dparalogistic
#' 
#' @export
#' @rdname distrMode
#' 
paralogisticMode <- 
function(scale = 1, 
         shape1.a)
{
  if (shape1.a <= 1) {
    warning("the density is not continuous at the mode.")
    0
  } else {
    scale*((shape1.a-1)/(shape1.a^2+1))^(1/shape1.a)
  }
}


# Pareto distribution
#' @inheritParams VGAM::dpareto
#' 
#' @export
#' @rdname distrMode
#' 
#' @examples 
#' curve(VGAM::dpareto(x, scale = 1, shape = 1), xlim = c(0, 10))
#' abline(v = paretoMode(scale = 1), col = 2)
#' 
paretoMode <-
function(scale = 1,
         ...)
{
  warning("the density is not continuous at the mode.")  
  scale
}


# Rayleigh distribution       
#' @inheritParams VGAM::drayleigh
#' 
#' @export
#' @rdname distrMode
#' 
rayleighMode <-
function(scale = 1)
{
  return(scale)
}


# Stable distribution
#' @inheritParams stabledist::dstable
#' 
#' @importFrom stabledist dstable
#' @importFrom stats optimize
#' @export
#' @rdname distrMode
#' 
stableMode <-
function(alpha,
         beta,
         gamma = 1,
         delta = 0,
         pm = 0, 
         ...) 
{
  beta.max <- 1 - 1e-11
  tol <- .Machine$double.eps^0.25
  if (gamma == 1 & delta == 0 & pm == 0) {
    stopifnot(0 < alpha, alpha <= 2, length(alpha) == 1, -1 <= beta, beta <= 1, length(beta) == 1, length(beta.max) == 1)
    if (alpha * beta == 0) {
      M <- 0
    }
    if (beta > beta.max) {
      beta <- beta.max
    }
    M <- stats::optimize(stabledist::dstable, interval = c(-0.7, 0) * sign(beta), alpha = alpha, 
                         beta = beta, pm = 0, maximum = TRUE, tol = tol)$maximum    
    return(M)
  } else {
    warning("still to be done. 'NA' is returned")
    return(NA)
  }
}


# from package 'stable'
#' @inheritParams stable::stable.mode
#' 
#' @importFrom stable stable.mode
#' @export
#' @rdname distrMode
#' 
stableMode2 <- 
function(loc, 
         disp, 
         skew, 
         tail)
{
  stable::stable.mode(loc, disp, skew, tail)
}


# Weibull distribution (in the context of extreme value theory)
# 
# #' @export
# #' @rdname distrMode
# #' 
# nweibullMode <-
# function(loc = 0,
#          scale = 1,
#          shape = 1)
# {
#   #return(loc + (scale/-abs(shape))*((1-abs(shape))^(abs(shape))-1))
#   if (shape < 1) {
#     M <- loc
#   } else {
#     M <- loc - scale*((shape-1)/shape)^(1/shape)
#   }
#   M
# }


# Student distribution
#' @inheritParams stats::dt
#' 
#' @export
#' @rdname distrMode
#' 
tMode <-
function(df,
         ncp)
{
  if (ncp == 0) {
    M <- 0
  } else {
    warning("still to be done. 'NA' is returned")
    M <- NA
  }
  M
}


# Uniform distribution
#' @inheritParams stats::dunif
#' 
#' @export
#' @rdname distrMode
#' 
unifMode <-
function(min = 0,
         max = 1)
{
  warning("all values between 'min' and 'max' are modes, only the mean value is returned")
  (min+max)/2
}


# Weibull distribution
#' @inheritParams stats::dweibull
#' 
#' @export
#' @rdname distrMode
#' 
weibullMode <-
function(shape,
         scale = 1)
{
  scale*(1-1/shape)^(1/shape)
}


# Yules-Simon distribution
#' @inheritParams VGAM::dyules
#' 
#' @export
#' @rdname distrMode
#' 
yulesMode <- 
function(...)
{
  1
}


## Discrete distributions
#------------------------------------------------------

# Bernoulli distribution
#' @inheritParams statip::dbern
#' 
#' @export
#' @rdname distrMode
#' 
bernMode <-
function(prob)
{
  if (prob > 1 || prob < 0) return(NaN)
  q <- 1 - prob
  if (q > prob) return(0)
  if (q == prob) {
    c(0,1)
  } else {
    1
  }
}


# Binomial distribution
#' @inheritParams stats::dbinom
#' 
#' @export
#' @rdname distrMode
#' 
binomMode <-
function(size,
         prob)
{
  if (prob > 1 || prob < 0) return(NaN)
  if (prob == 0) {
    return(0)
  } else {
    if (prob == 1) {
      return(size)
    } else {
      x <- ceiling((size+1)*prob - 1)
      if (x == (size+1)*prob - 1) {
        return(c(x, x+1))
      } else {
        return(x)
      }      
    }
  }     
}


# Geometric distribution
#' @inheritParams stats::dgeom
#' 
#' @export
#' @rdname distrMode
#' 
geomMode <-
function(...)
{
  return(1)
}


# Hypergeometric distribution
#' @inheritParams stats::dhyper
#' 
#' @export
#' @rdname distrMode
#' 
hyperMode <-
function(m,
         n,
         k,
         ...)
{
  lambda <- (m+1)*(k+1)/(m+n+1)  
  if (lambda == 0) return(0)
  x <- floor(lambda)
  if (lambda == x) {
    return(c(lambda-1,lambda))
  } else {
    return(x)
  }
}


# Negative binomial distribution
#' @inheritParams stats::dnbinom
#' 
#' @export
#' @rdname distrMode
#' 
nbinomMode <-
function(size,
         prob,
         mu)
{
  if (prob > 1 || prob < 0) return(NaN)
  if (!missing(mu)) {
    prob <- size/(size+mu)    
  }
  if (size <= 1) {
    0
  } else {
    floor((size-1)*(1-prob)/prob)
  }
}


# Poisson distribution
#' @inheritParams stats::dpois
#' 
#' @export
#' @rdname distrMode
#' 
#' @examples 
#' ## Poisson distribution
#' poisMode(lambda = 6)
#' poisMode(lambda = 6.1)
#' mlv("poisson", lambda = 6.1)
#' 
poisMode <-
function(lambda)
{
  if (lambda < 0) return(NaN)
  if (lambda == 0) return(0)
  x <- floor(lambda)
  if (lambda == x) {
    c(lambda-1,lambda)
  } else {
    x
  }
}
