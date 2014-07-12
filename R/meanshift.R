# Author: P. Poncet

meanshift <-
function(x,
         bw = NULL,
         kernel = "gaussian",
         par = shorth(x),
         iter = 1000,
         eps = 1e-08)
{
  if (is.null(bw)) bw <- bw.SJ(x)
  if (is.character(bw)) {
    if (nx < 2) 
      stop("need at least 2 points to select a bandwidth automatically")
    bw <- switch(tolower(bw), nrd0 = bw.nrd0(x), nrd = bw.nrd(x), 
                 ucv = bw.ucv(x), bcv = bw.bcv(x), sj = , 
                 "sj-ste" = bw.SJ(x, method = "ste"), 
                 "sj-dpi" = bw.SJ(x, method = "dpi"), 
                 stop("unknown bandwidth rule"))
  }
  s <- 0
  #th <- rep(0, iter)
  #Mrec <- c()
  x0 <- par
  for (j in 1:iter) {
    z <- (x - par)/bw
    k <- do.call(paste(".kernel.", kernel, sep = ""), list(z))$k
    M <- crossprod(x, k)/sum(k)
    #Mrec[j] <- M
    if (is.nan(M)) stop("sum(g) is zero in the meanshift function. Change the bandwidth 'bw' or the initial value 'par'.")
    #th[j] <- sum((M - par)^2)/sum(par^2)
    th <- sum((M - par)^2)/sum(par^2)
    #if (th[j] < eps) {
    if (th < eps) {
      s <- j
      break
    }
    par <- M
  }
  #attr(M, "meanshift.points") <- Mrec[1:s]
  #attr(M, "threshold.values") <- th[1:s]
  attr(M, "iterations") <- s
  #attr(M, "start") <- x0
  return(M)
}
