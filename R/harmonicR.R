#' Harmonic regression model
#' 
#' Computes amplitudes and phase angle in the typical harmonic regression framework.
#' Based on these estimates a harmonic regression function is fitted.
#' 
#' @param y             numeric vector containing time series on which harmonic 
#'                      regression will be fitted; missing values are not allowed
#' @param lenBasePeriod numeric
#' @param numFreq       numeric indicating the total number of frequencies to be 
#'                      used in the harmonic regression.
#' @param delta         numeric giving a (small) regularization parameter to prevent 
#'                      non-invertible hat matrix (see \code{details}), probably caused by high
#'                      amplitudes.
#'                      
#' @export
#' 
#' @details This function does not allow missing data, missing values of the original
#' data \code{y} must be substituted by values considerably out of the range of the observations.
#' 
#' @examples 
#' y <- c(5, 2, 5, 10, 12, 18, 20, 23, 27, 30, 40, 60, 66, 70, 90, 120, 
#'      160, 190, 105, 210, 104, 200, 90, 170, 50, 120, 80, 60, 50, 40, 
#'      30, 28, 24, 20, 15, 10)
#' y2 <- y[-c(19, 21, 23, 25)]
#' fit <- harmonicR(y = y2, numFreq = 5, delta = 0.01)
#' (fit)
#' plot(y2, pch = 16, col = "lightgray")
#' lines(fit$fitted)
#' 
#' y3 <- y[-c(16:18, 20, 22, 24, 26)]
#' fit1 <- harmonicR(y = y3, numFreq = 5, delta = 0.01)
#' fit2 <- harmonicR(y = y3, numFreq = 8, delta = 0.01)
#' fit3 <- harmonicR(y = y3, numFreq = 10, delta = 0.01)
#' plot(y3, pch = 16, col = "lightgray")
#' lines(fit1$fitted)
#' lines(fit2$fitted, lty = 4)
#' lines(fit3$fitted, lty = 2)
#' 
#' @seealso \code{\link[geoTS]{hants}}
#' 
#' @return A list containing:
#' \item{amplitude}{a numeric vector with amplitude estimates.}
#' \item{phase}{a numeric vector with phase estimates.}
#' \item{fitted}{a numeric vector with fitted values via harmonic regression.}
#' 
harmonicR <- function(y, lenBasePeriod = length(y), numFreq, delta){
  
  numImages <- length(y)
  # mat <- matrix(0, nrow = min(2 * numFreq + 1, numImages), ncol = numImages)
  amp <- numeric(numFreq + 1)
  phi <- numeric(numFreq + 1)
  ySmooth <- numeric(numImages)
  
  numRows <- min(2 * numFreq + 1, numImages)
  # numOutMax <- numImages - numRows - degreeOverDeter
  degree <- 180 / pi
  
  mat <- matrix(0, nrow = numImages, ncol = min(2 * numFreq + 1, numImages))
  vecTs <- 2 * pi * (0:(lenBasePeriod-1)) / lenBasePeriod
  mat[,1] <- 1
  mat[, seq(2, 2 * numFreq, 2)] <- sapply(1:numFreq, function(s) cos( s * vecTs ))
  mat[, seq(2, 2 * numFreq, 2)+1] <- sapply(1:numFreq, function(s) sin( s * vecTs ))
  
  za <- t(mat) %*% ( y )
  
  A <- t(mat) %*% mat
  A <- A + diag(rep(1, numRows)) * delta
  A[1,1] <- A[1,1] - delta
  zr <- solve(A) %*% za
  
  ySmooth <- (mat) %*% zr
  
  amp[1] <- zr[1]
  phi[1] <- 0
  
  zr[numImages+1] <- 0
  
  i <- seq(2, numRows, by = 2)
  iFreq <- (i+2)/2
  ra <- zr[i]
  rb <- zr[i+1]
  amp[iFreq] <- sqrt( ra * ra + rb * rb )
  phase <- atan2(rb, ra) * degree
  phase[phase < 0] <- phase[phase < 0] + 360
  phi[iFreq] <- phase
  
  list(amplitude = amp, phase = phi, fitted = as.numeric(ySmooth))  
}
