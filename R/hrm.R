#' @title 
#' Bickel's half-range mode estimator
#' 
#' @description 
#' SINCE THIS FUNCTION USED TO DEPEND ON THE BIOCONDUCTOR PACKAGE 'GENEFILTER', 
#' IT IS CURRENTLY DEFUNCT.
#' 
#' This function computes Bickel's half range mode estimator 
#' described in Bickel (2002). It is a wrapper around the function 
#' \code{half.range.mode} from package \pkg{genefilter}. 
#'
#' @details 
#' The mode estimator is computed by iteratively identifying 
#' densest half ranges. A densest half range is an interval 
#' whose width equals half the current range, and which 
#' contains the maximal number of observations. 
#' The subset of observations falling in the selected 
#' densest half range is then used to compute a new range, 
#' and the procedure is iterated.
#'
#' @note 
#' The user may call \code{hrm} through 
#' \code{mlv(x, method = "hrm", bw, ...)}. 
#' 
#' @references 
#' \itemize{
#'   \item Bickel D.R. (2002). 
#'   Robust estimators of the mode and skewness of continuous data. 
#'   \emph{Computational Statistics and Data Analysis}, \bold{39}:153-163.
#'   
#'   \item Hedges S.B. and Shah P. (2003). 
#'   Comparison of mode estimation methods and application in molecular clock analysis. 
#'   \emph{BMC Bioinformatics}, \bold{4}:31-41.
#'   
#'   \item Bickel D.R. and Fruehwirth R. (2006). 
#'   On a Fast, Robust Estimator of the Mode: 
#'   Comparisons to Other Robust Estimators with Applications. 
#'   \emph{Computational Statistics and Data Analysis}, \bold{50}(12):3500-3530. 
#' }
#' 
#' @param x 
#' numeric. Vector of observations. 
#' 
#' @param bw
#' numeric. The bandwidth to be used. Should belong to (0, 1]. 
#' This gives the fraction of the observations to consider at 
#' each step of the iterative algorithm. 
#' 
#' @param ... 
#' Additional arguments. 
#' 
#' @return 
#' A numeric value is returned, the mode estimate.
#'
#' @author The C and R code are due to Richard Bourgon \email{bourgon@stat.berkeley.edu}, 
#' see package \pkg{genefilter}. The algorithm is described in Bickel (2002).
#' 
#' @seealso 
#' \code{\link[modeest]{mlv}()} for general mode estimation; 
#' \code{\link[modeest]{hsm}()} for the half sample mode;  
#' \code{\link[modeest]{venter}()} for the Venter mode estimate. 
#' 
# #' @importFrom genefilter half.range.mode
#' @export
#' @aliases HRM
#' 
#' @examples 
#' \dontrun{
#' # Unimodal distribution 
#' x <- rgamma(1000, shape = 31.9)
#' ## True mode
#' gammaMode(shape = 31.9)
#' 
#' ## Estimate of the mode
#' hrm(x, bw = 0.4)
#' mlv(x, method = "hrm", bw = 0.4)
#' }
#' 
hrm <-
function(x,
         bw = NULL,
         ...) # TODO: introduce a 'k' argument?
{
  .Defunct("genefilter::half.range.mode",
           msg = paste0("currently 'hrm()' is due to difficulties ",
                        "in installing automatically BioConductor's ",
                        "'genefilter' dependency"))
  if (is.null(bw)) stop("argument 'bw' is missing", call. = FALSE)
  if (bw <= 0 || bw > 1) {
    stop("argument 'bw' must belong to (0, 1]", 
         call. = FALSE)
  }
  
  y <- sort(x)
  #genefilter::half.range.mode(y, beta = bw)
} 
