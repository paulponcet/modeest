
############## FUNCTIONS ############

MS_max_representation <-
function(x,                    #Vector of the mono-valued dataset
         percentMargin = 0.05, #Size of the margin, so that the extremal value are not stuck to the border of the image
         sizeKerMin = 0.001,   #Minimal value for the size of the kernel
         sizeKerMax = 1,       #Maximal value for the size of the kernel
         bwLen = 500,          #Number of convolutions with a different kernel. It corresponds to the number of lines in the display
         width = 500,          #Width of the display
         jitterOrHist = 0)     #Flag indicating the representation of the data in the lower part of the graphical representation. - 0 : automatic 1 : jittered density diagram 2 : histogram
{

  #définition des dimensions de la matrice image
  minx <- min(x);
  maxx <- max(x);
  rangex <- maxx - minx;
  newMinx <- minx-percentMargin*rangex; #on va jusqu'à un pourcentage au-dela des valeurs extremales
  newMaxx <- maxx+percentMargin*rangex;
  axeAbs <- seq(from = newMinx, to = newMaxx, length = width);
  precision <- axeAbs[2]-axeAbs;
  axeOrd <- seq(from = sizeKerMin , to = sizeKerMax , length = bwLen);
  MatConv <- matrix(rep(0, length(axeAbs)*(length(axeOrd))), ncol = length(axeAbs));

  # Remplissage de la matrice image
  for (ibw in 1:length(axeOrd)) {
    mode <- density(x, bw=axeOrd[ibw], kernel = "gaussian", n=length(axeAbs), from=newMinx, to=newMaxx);
    valueLine <- mode$y/max(mode$y);
    maxLine <- .localMode(valueLine);
    MatConv[ibw,] <- valueLine + maxLine ;
  }

  # Display
  modeInit <- density(x, bw=axeOrd[1], kernel = "gaussian", n=length(axeAbs), from=newMinx, to=newMaxx);
  xj <- jitter(x, factor=1, amount = NULL);
  axeInv <- seq(from=-5, to=5, length = bwLen)
  taillesKer <- round(length(axeOrd)/40)
  sizeAxis <- axeOrd[round(length(axeOrd))]
  mGaus1 <- dnorm(0, mean = 0, sd = 0.1)
  mGaus2 <- dnorm(0, mean = 0, sd = 1)
  mGaus3 <- dnorm(0, mean = 0, sd = 2.4)
  X11()
  nf <- layout(matrix(c(1,3,2,4),2,2,byrow=TRUE), c(3,1), c(3,1), TRUE)
  layout.show(nf)
  par(mar=c(3,3,1,1))
  image(t(MatConv), col = heat.colors(24), axes=F);
  par(mar=c(2,3,1,1))
  if(jitterOrHist == 0) { 
    if(length(x < 150) {
      jitterOrHist = 1
    } else{ 
      jitterOrHist = 2
    }
  }
  if(jitterOrHist == 1) { 
    plot(xj, rep(1,length(x)), type = "h", frame.plot=F, yaxt = "n", ylab="", ylim=c(0,1))
  }
  if(jitterOrHist == 2) { 
    hist(x, frame.plot=F, yaxt = "n", ylab="", main="")
  }
  par(mar=c(3,1,1,1))
  plot(axeInv, axeOrd, col="white", frame.plot=F, xaxt = "n")
  lines(axeInv, axeOrd[1]+ sizeAxis/(8*mGaus1)*dnorm(axeInv, mean = 0, sd = 0.1), type="l")
  lines(axeInv, sizeAxis*3/8+ sizeAxis/(8*mGaus2)*dnorm(axeInv, mean = 0, sd = 1), type="l")
  lines(axeInv, sizeAxis*7/8+ sizeAxis/(8*mGaus3)*dnorm(axeInv, mean = 0, sd = 2.4), type="l")
}

###
.localMode <-
function(y,
         threshold = 0.000001)
{
  ll <- length(y)
  ind <- (y[1:(ll-2)]< y[2:(ll-1)]) & (y[2:(ll-1)] > y[3:ll]) * (y[2:(ll-1)] > threshold);
  indnumResize <- c(0, as.numeric(ind), 0);
}
