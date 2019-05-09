#' Harmonic analysis for time series: fitting the harmonic regression model under the presence
#' of outliers 
#' 
#' iterative algorithm that computes amplitudes and phase angles in the 
#' typical harmonic regression framework. Based on these estimates a harmonic 
#' regression function is fitted. As part of the iterative algorithm, observations
#' are being excluded from the design matrix of the regression model if the distance
#' between them and the fitted curve exceeds the value of the parameter \code{fitErrorTol}.
#' \code{hants} is based on implementations of the same name written in Fortran 
#' and Matlab computer languages.
#' 
#' @param y               numeric vector containing time series on which harmonic 
#'                        regression will be fitted. Missing values are not allowed
#' @param ts              numeric vector of length \code{y} with the sampling
#'                        points for \code{y}. Default is \eqn{ts[i] = i, i=1,\ldots, 
#'                        \code{y}}
#' @param lenBasePeriod   numeric giving the length of the base period, reported in 
#'                        samples, e.g. days, dekads, months, years, etc.
#' @param numFreq         numeric indicating the total number of frequencies to be 
#'                        used in the harmonic regression
#' @param HiLo            character indicating whether high or low outliers must be rejected
#' @param low             numeric giving minimum valid value of fitted harmonic regression 
#'                        function
#' @param high            numeric giving maximum valid value of fitted harmonic regression
#'                        function
#' @param fitErrorTol     numeric giving maximum allowed distance between observations and fitted
#'                        curve; if difference between a given observation and its fitted value 
#'                        exceeds \code{fitErrorTol} then the observation is excluded from
#'                        next iteration of fitting algorithm
#' @param degreeOverDeter numeric; iteration stops when number of observations equals
#'                        number of observations for curve fitting plus \code{degreeOverDeter};
#'                        the latter in turns is by definition \code{y} minus
#'                        \eqn{min(2 * \code{numFreq} +1, \code{y})}.
#' @param delta           numeric giving a (small) regularization parameter to prevent 
#'                        non-invertible hat matrix (see \code{details}), probably caused by high
#'                        amplitudes.
#'                        
#' @export
#' 
#' @details This function does not allow missing data, missing values of the original
#' data \code{y} must be substituted by values considerably out of the range of the observations.
#' 
#' @examples
#' y <- c(5, 2, 5, 10, 12, 18, 20, 23, 27, 30, 40, 60, 66, 
#'       70, 90, 120, 160, 190, 105, 210, 104, 200, 90, 170, 
#'       50, 120, 80, 60, 50, 40, 30, 28, 24, 20, 15, 10)
#' fitLow <- hants(y = y, numFreq = 3, HiLo = "Lo", low = 0, 
#'           high = 255, fitErrorTol = 5, degreeOverDeter = 1, delta = 0.1)
#' fitHigh <- hants(y = y, numFreq = 3, HiLo = "Hi", low = 0, 
#'           high = 255, fitErrorTol = 5, degreeOverDeter = 1, delta = 0.1)
#' fit_harmR <- harmonicR(y = y, numFreq = 5, delta = 0.1)
#' plot(y, pch = 16, main = "hants fitting")
#' lines(fitLow$fitted, lty = 4, col = "red")
#' lines(fitHigh$fitted, lty = 2, col = "blue")
#' lines(fit_harmR$fitted, lty = 4, col = "green")
#' 
#' # A missing value should be substituted by a number out of the valid range
#' y1 <- y
#' y1[20] <- -10
#' fitLowMiss <- hants(y = y1, numFreq = 3, HiLo = "Lo", low = 0, high = 255, 
#'               fitErrorTol = 5, degreeOverDeter = 1, delta = 0.1)
#' fitHighMiss <- hants(y = y1, numFreq = 3, HiLo = "Hi", low = 0, 
#'                high = 255, fitErrorTol = 5, degreeOverDeter = 1, delta = 0.1)
#' fit_harmRMiss <- harmonicR(y = y1, numFreq = 3, delta = 0.1)
#' plot(y1, pch = 16, main = "hants fitting (missing values)")
#' lines(fitLowMiss$fitted, lty = 4, col = "red")
#' lines(fitHighMiss$fitted, lty = 2, col = "blue")
#' lines(fit_harmRMiss$fitted, lty = 4, col = "green")
#' 
#' @seealso \code{\link[geoTS]{harmonicR}}
#' 
#' @return A list containing:
#' \item{amplitude}{a numeric vector with amplitude estimates.}
#' \item{phase}{a numeric vector with phase estimates.}
#' \item{fitted}{a numeric vector with fitted values via harmonic regression.}
#' 
#' @references Roerink, G.J., Menenti, M., Verhoef, W. (2000).
#' \emph{Reconstructing clodfree NDVI composites using Fourier analysis of time series}, 
#' Int. J. Remote Sensing, \bold{21(9)}, 1911--1917.
#'   
#' @references The Matlab implementation of HANTS can be found 
#' \href{https://nl.mathworks.com/matlabcentral/fileexchange/38841-matlab-implementation-of-harmonic-analysis-of-time-series-hants}{here}.
#' 
hants <- function(y, ts = 1:length(y), lenBasePeriod = length(y), numFreq, 
                  HiLo = c("Hi", "Lo"), low, high, fitErrorTol, 
                  degreeOverDeter, delta){
  
  numImages <- length(y)
  mat <- matrix(0, nrow = min(2 * numFreq + 1, numImages), ncol = numImages)
  amp <- numeric(numFreq + 1)
  phi <- numeric(numFreq + 1)
  ySmooth <- numeric(numImages)

  HiLo <- match.arg(HiLo)
  sHiLo <- 0
  
  if(HiLo == "Hi"){
    sHiLo <- -1
  } else {
    sHiLo <- 1
  }
  
  numRows <- min(2 * numFreq + 1, numImages)
  numOutMax <- numImages - numRows - degreeOverDeter
  degree <- 180 / pi
  mat[1,] <- 1
  
  ang <- 2 * pi * (0:(lenBasePeriod-1)) / lenBasePeriod
  cs <- cos(ang)
  sn <- sin(ang)
  
  i <- 1:numFreq
  for( j in 1:numImages ){
    index <- 1 + (i * (ts[j]-1)) %% lenBasePeriod
    mat[2*i, j] <- cs[index]
    mat[2*i+1, j] <- sn[index]
  }
  
  p <- rep(1, numImages)
  
  p[y < low | y > high] <- 0
  numOut <- sum( p==0 )
  
  if(numOut > numOutMax){
    cat("numOut:", numOut, "numOutMax:", numOutMax)
    stop("Not enough data points")
  }
  
  ready <- F
  nLoop <- 0
  nLoopMax <- numImages
  
  while( !ready & nLoop < nLoopMax ){
    nLoop <- nLoop + 1
    za <- mat %*% ( p * y )
    
    A <- mat %*% diag(p) %*% t(mat)
    A <- A + diag(rep(1, numRows)) * delta
    A[1,1] <- A[1,1] - delta
    zr <- solve(A) %*% za
    
    ySmooth <- t(mat) %*% zr
    weightedResiduals <- sHiLo * (ySmooth - y)
    error <- p * weightedResiduals
    
    # rankVec <- sort(error)
    rankVec <- order(error)
    
    if( floor(rankVec[numImages]) == 0  ){
      maxError <- 0
    } else {
      maxError <- weightedResiduals[rankVec[numImages]]
    }
    
    ready <- (maxError < fitErrorTol) | (numOut == numOutMax)
    
    if(!ready){
      i <- numImages
      j <- rankVec[i]
      while( (p[j] * weightedResiduals[j] > maxError * 0.5) & (numOut < numOutMax) ){
        p[j] <- 0
        numOut <- numOut + 1
        i <- i - 1
        j <- rankVec[i]
      }
    }
  }
  
  amp[1] <- zr[1]
  phi[1] <- 0
  
  zr[numImages+1] <- 0
  
  i <- seq(2, numRows, by = 2)
  iFreq <- (i+2)/2
  ra <- zr[i]
  rb <- zr[i+1]
  amp[iFreq] <- sqrt( ra * ra + rb * rb )
  phase <- atan2(rb, ra) * degree
  phase[phase < 0] <- phase[phase <0] + 360
  phi[iFreq] <- phase
  
  list(amplitude = amp, phase = phi, fitted = ySmooth)  
}
