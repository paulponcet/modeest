\name{modeest}
\alias{modeest-package}
\alias{modeest}

\docType{package}
\title{
Mode Estimation 
}
\description{
This package intends to provide estimators of the mode of univariate 
unimodal (and sometimes multimodal) data and values of the modes of usual probability distributions. 
%It also includes computation of the density function, distribution function, and quantile 
%function of the Chernoff distribution, which is the limiting distribution of the
%Chernoff mode estimator.

For a complete list of functions, use \code{library(help = "modeest")} or \code{help.start()}.
}
\details{
\tabular{ll}{
Package: \tab modeest\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2012-10-15\cr
License: \tab GPL version 2 or newer\cr
}
}
\author{
P. Poncet 
}
\references{
  \itemize{
    \item Parzen E. (1962).
    On estimation of a probability density function and mode.
    \emph{Ann. Math. Stat.}, \bold{33}(3):1065-1076. 
    
    \item Chernoff H. (1964).
    Estimation of the mode.
    \emph{Ann. Inst. Statist. Math.}, \bold{16}:31-41.

    \item Huber P.J. (1964).
    Robust estimation of a location parameter.
    \emph{Ann. Math. Statist.}, \bold{35}:73-101.

    \item Dalenius T. (1965).
    The Mode - A Negleted Statistical Parameter.
    \emph{J. Royal Statist. Soc. A}, \bold{128}:110-117.

    \item Grenander U. (1965).
    Some direct estimates of the mode.
    \emph{Ann. Math. Statist.}, \bold{36}:131-138.

    \item Venter J.H. (1967).
    On estimation of the mode.
    \emph{Ann. Math. Statist.}, \bold{38}(5):1446-1455.
    
    \item Lientz B.P. (1969).
    On estimating points of local maxima and minima of density functions.
    \emph{Nonparametric Techniques in Statistical Inference (ed. M.L. Puri, Cambridge University Press)}, p.275-282.

    \item Lientz B.P. (1970).
    Results on nonparametric modal intervals.
    \emph{SIAM J. Appl. Math.}, \bold{19}:356-366.

    \item Wegman E.J. (1971).% a revoir !!
    A note on the estimation of the mode.
    \emph{Ann. Math. Statist.}, \bold{42}(6):1909-1915.

    \item Yamato H. (1971). % a revoir !!
    Sequential estimation of a continuous probability density function and mode.
    \emph{Bull. Math. Statist.}, \bold{14}:1-12.

    \item Ekblom H. (1972).
    A Monte Carlo investigation of mode estimators in small samples.
    \emph{Applied Statistics}, \bold{21}:177-184.

    \item Lientz B.P. (1972).
    Properties of modal intervals.
    \emph{SIAM J. Appl. Math.}, \bold{23}:1-5.

    \item Konakov V.D. (1973).
    On the asymptotic normality of the mode of multidimensional distributions.
    \emph{Theory Probab. Appl.}, \bold{18}:794-803.

    \item Robertson T. and Cryer J.D. (1974).
    An iterative procedure for estimating the mode.
    \emph{J. Amer. Statist. Assoc.}, \bold{69}(348):1012-1016.

    \item Kim B.K. and Van Ryzin J. (1975).
    Uniform consistency of a histogram density estimator and modal estimation.
    \emph{Commun. Statist.}, \bold{4}:303-315.

    \item Sager T.W. (1975).% a revoir !!
    Consistency in nonparametric estimation of the mode.
    \emph{Ann. Statist.}, \bold{3}(3):698-706.
   
    \item Stone C.J. (1975).
    Adaptive maximum likelihood estimators of a location parameter.
    \emph{Ann. Statist.}, \bold{3}:267-284.

    \item Mizoguchi R. and Shimura M. (1976).
    Nonparametric Learning Without a Teacher Based on Mode Estimation.
    \emph{IEEE Transactions on Computers}, \bold{C25}(11):1109-1117.

    \item Adriano K.N., Gentle J.E. and Sposito V.A. (1977).
    On the asymptotic bias of Grenander's mode estimator.
    \emph{Commun. Statist.-Theor. Meth. A}, \bold{6}:773-776.

    \item Asselin de Beauville J.-P. (1978).
    Estimation non parametrique de la densite et du mode, exemple de la distribution Gamma.
    \emph{Revue de Statistique Appliquee}, \bold{26}(3):47-70.

    \item Sager T.W. (1978).% a revoir !!
    Estimation of a multivariate mode.
    \emph{Ann. Statist.}, \bold{6}:802-812.

    \item Devroye L. (1979).% a revoir !!
    Recursive estimation of the mode of a multivariate density.
    \emph{Canadian J. Statist.}, \bold{7}(2):159-167.

    \item Sager T.W. (1979).% a revoir !!
    An iterative procedure for estimating a multivariate mode and isopleth.
    \emph{J. Amer. Statist. Assoc.}, \bold{74}(366):329-339.

    \item Eddy W.F. (1980).
    Optimum kernel estimators of the mode.
    \emph{Ann. Statist.}, \bold{8}(4):870-882.

    \item Eddy W.F. (1982).
    The Asymptotic Distributions of Kernel Estimators of the Mode.
    \emph{Z. Wahrsch. Verw. Gebiete}, \bold{59}:279-290.

    \item Hall P. (1982).
    Asymptotic Theory of Grenander's Mode Estimator.
    \emph{Z. Wahrsch. Verw. Gebiete}, \bold{60}:315-334.

    \item Sager T.W. (1983).
    Estimating modes and isopleths.
    \emph{Commun. Statist.-Theor. Meth.}, \bold{12}(5):529-557.

    \item Hartigan J.A. and Hartigan P.M. (1985).
    The Dip Test of Unimodality.
    \emph{Ann. Statist.}, \bold{13}:70-84.

    \item Hartigan P.M. (1985).
    Computation of the Dip Statistic to Test for Unimodality.
    \emph{Appl. Statist. (JRSS C)}, \bold{34}:320-325.

    \item Romano J.P. (1988).
    On weak convergence and optimality of kernel density estimates of the mode.
    \emph{Ann. Statist.}, \bold{16}(2):629-647.

    \item Tsybakov A. (1990).
    Recursive estimation of the mode of a multivariate distribution.
    \emph{Probl. Inf. Transm.}, \bold{26}:31-37.

    \item Hyndman R.J. (1996).
    Computing and graphing highest density regions.
    \emph{Amer. Statist.}, \bold{50}(2):120-126.

    \item Vieu P. (1996).
    A note on density mode estimation.
    \emph{Statistics \& Probability Letters}, \bold{26}:297--307.

    \item Leclerc J. (1997).
    Comportement limite fort de deux estimateurs du mode : le shorth et l'estimateur naif.
    \emph{C. R. Acad. Sci. Paris, Serie I}, \bold{325}(11):1207-1210.

    %\item Minnotte M. C. (1997).
    %Nonparametric testing of the existence of modes.
    %\emph{Ann. Statist.}, \bold{25}(4):1646-1660.
    %
    %\item Futschik A. (1999).
    %A new estimate of the mode based on the quantile density.
    %\emph{Statistics and Probability Letters}, \bold{43}:145-152.
    %
    \item Leclerc J. (2000).
    Strong limiting behavior of two estimates of the mode: the shorth and the naive estimator.
    \emph{Statistics and Decisions}, \bold{18}(4).

    %\item Groeneboom P.and Wellner J.A. (2001).
    %Computing Chernoff's distribution.
    %\emph{J. Comput. Graph. Statist.}, \bold{10}:388-400.
    %
    \item Shoung J.M. and Zhang C.H. (2001).
    Least squares estimators of the mode of a unimodal regression function.
    \emph{Ann. Statist.}, \bold{29}(3):648-665.

    \item Bickel D.R. (2002).
    Robust estimators of the mode and skewness of continuous data.
    \emph{Computational Statistics and Data Analysis}, \bold{39}:153-163.

    \item Abraham C., Biau G. and Cadre B. (2003).
    Simple Estimation of the Mode of a Multivariate Density.
    \emph{Canad. J. Statist.}, \bold{31}(1):23-34.

    \item Bickel D.R. (2003).
    Robust and efficient estimation of the mode of continuous data: The mode as a viable measure of central tendency.
    \emph{J. Statist. Comput. Simul.}, \bold{73}:899-912.

    \item Djeddour K., Mokkadem A. et Pelletier M. (2003).
    Sur l'estimation recursive du mode et de la valeur modale d'une densite de probabilite.
    \emph{Technical report 105}.

    \item Djeddour K., Mokkadem A. et Pelletier M. (2003).
    Application du principe de moyennisation a l'estimation recursive du mode et de la valeur modale d'une densite de probabilite.
    \emph{Technical report 106}.

    \item Hedges S.B. and Shah P. (2003).
    Comparison of mode estimation methods and application in molecular clock analysis.
    \emph{BMC Bioinformatics}, \bold{4}:31-41.

    %\item Ziegler K. (2003).
    %On the asymptotic normality of kernel regression estimators of the mode in the nonparametric random design model.
    %\emph{Journal of Statistical Planning and Inference}, \bold{115}:123-144.
    %
    \item Herrmann E. and Ziegler K. (2004).
    Rates of consistency for nonparametric estimation of the mode in absence of smoothness assumptions.
    \emph{Statistics and Probability Letters}, \bold{68}:359-368.

    \item Abraham C., Biau G. and Cadre B. (2004).
    On the Asymptotic Properties of a Simple Estimate of the Mode.
    \emph{ESAIM Probab. Stat.}, \bold{8}:1-11.

    \item Mokkadem A. and Pelletier M. (2005).
    Adaptive Estimation of the Mode of a Multivariate Density.
    \emph{J. Nonparametr. Statist.}, \bold{17}(1):83-105.

    \item Bickel D.R. and Fruehwirth R. (2006).
    On a Fast, Robust Estimator of the Mode: Comparisons to Other Robust Estimators with Applications.
    \emph{Computational Statistics and Data Analysis}, \bold{50}(12):3500-3530.
    
  }
} 

\keyword{ distribution  }
\keyword{ univar }
\keyword{ robust }
\keyword{ nonparametric }
\seealso{ \code{\link{mlv}} for general mode estimation%;
}
