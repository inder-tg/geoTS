#' Harmonic analysis for time series
#' 
#' Fits harmonic regression (\code{harmR}) model, that is, computes amplitudes and phase angles 
#' in the typical harmonic regression framework. Based on these estimates a harmonic regression 
#' function is fitted. Also fits \code{hants}, a popular iterative algorithm that computes amplitudes and phase angles in the 
#' harmonic regression framework. As part of the iterative algorithm, observations
#' are being excluded from the design matrix of the regression model if the distance
#' between them and the fitted curve exceeds the value of the parameter \code{fitErrorTol}.
#' \code{hants} is based on implementations with the same name written in Fortran 
#' and Matlab computer languages.
#' 
#' @param              y  numeric vector containing time series on which harmonic 
#'                        regression will be fitted. Missing values are not allowed.
#' @param         method  character specifying algorithm to apply: \code{ols_harmR} (default),
#'                        \code{wls_harmR} (heteroscedastic model) or \code{hants}. 
#' @param          sigma  numeric vector of \code{length(y)} containing variance estimates.
#'                        Default set \code{NULL}.          
#' @param             ts  numeric vector of \code{length(y)} with the sampling
#'                        points for \code{y}. Default is \eqn{ts[i] = i, i=1,\ldots, 
#'                        \code{length(y)}}.
#' @param      lenPeriod  numeric giving the length of the base period, reported in 
#'                        samples, e.g. days, dekads, months, years, etc.
#' @param        numFreq  numeric indicating the total number of frequencies to be 
#'                        used in harmonic regression. For technical reasons, \code{2*numFreq+1}
#'                        must be lesser than \code{length(y)}. 
#' @param           HiLo  character indicating whether high or low outliers must be rejected
#'                        when \code{method=hants}.
#' @param            low  numeric giving minimum valid value of fitted harmonic regression 
#'                        function when \code{method=hants}.
#' @param           high  numeric giving maximum valid value of fitted harmonic regression
#'                        function when \code{method=hants}.
#' @param     fitErrorTol numeric giving maximum allowed distance between observations and fitted
#'                        curve; if difference between a given observation and its fitted value 
#'                        exceeds \code{fitErrorTol} then this observation will not be included
#'                        in the fitting procedure in the next iteration of the algorithm.
#' @param degreeOverDeter numeric; iteration stops when number of observations equals
#'                        number of observations for curve fitting plus \code{degreeOverDeter};
#'                        the latter in turns is by definition \code{length(y)} minus
#'                        \eqn{min(2 * \code{numFreq} +1, \code{length(y)})}.
#' @param           delta numeric (positive) giving a (small) regularization parameter to prevent 
#'                        non-invertible hat matrix (see \code{details}), probably caused by high
#'                        amplitudes.
#'                        
#' @export
#' 
#' @details Methods \code{ols_harmR} and \code{wls_harmR} do not allow missing values 
#' and utilizes parameters \code{y}, \code{lenBasePeriod}, \code{numFreq} and \code{delta} 
#' only. 
#' 
#' Method \code{hants} utilizes all the parameters presented above. This method 
#' does not allow missing values. Missing values in \code{y} must be substituted by values 
#' considerably out of observations range.
#' 
#' @note In previous versions: \code{method} allowed \code{harmR} which is equivalent 
#' to current \code{ols_harmR}); \code{lenBasePeriod} was replaced by \code{lenPeriod}.
#' 
#' @examples
#' y <- c(5, 2, 5, 10, 12, 18, 20, 23, 27, 30, 40, 60, 66,
#' 70, 90, 120, 160, 190, 105, 210, 104, 200, 90, 170,
#' 50, 120, 80, 60, 50, 40, 30, 28, 24, 20, 15, 10)
#' # -------------------------------------------------------------------------- 
#' fit_ols_harmR <- haRmonics(y = y, numFreq = 3, delta = 0.1)
#' fitLow_hants <- haRmonics(y = y, method = "hants", numFreq = 3, HiLo = "Lo", 
#'                          low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1, 
#'                          delta = 0.1)
#' fitHigh_hants <- haRmonics(y = y, method = "hants", numFreq = 3, HiLo = "Hi", 
#'                           low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1, 
#'                           delta = 0.1)
#' plot(y, pch = 16, main = "haRmonics fitting")
#' lines(fit_ols_harmR$fitted ,lty = 4, col = "green")
#' lines(fitLow_hants$fitted, lty = 4, col = "red")
#' lines(fitHigh_hants$fitted, lty = 2, col = "blue")
#' # -------------------------------------------------------------------------- 
#' # Substituting missing value by a number outside observations range
#' # -------------------------------------------------------------------------- 
#' y1 <- y
#' y1[20] <- -10
#' 
#' fitLow_hants_missing <- haRmonics(y = y1, method = "hants", numFreq = 3, HiLo = "Lo", 
#'                                  low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1, 
#'                                  delta = 0.1)
#' fitHigh_hants_missing <- haRmonics(y = y1, method = "hants", numFreq = 3, HiLo = "Hi", 
#'                                   low = 0, high = 255, fitErrorTol = 5, degreeOverDeter = 1, 
#'                                   delta = 0.1)
#' fit_ols_harmR_missing <- haRmonics(y = y1, numFreq = 3, delta = 0.1)
#' 
#' plot(y1, pch = 16, main = "haRmonics fitting (missing values)", ylim = c(-1,210))
#' lines(fitLow_hants_missing$fitted, lty = 4, col = "red")
#' lines(fitHigh_hants_missing$fitted, lty = 2, col = "blue")
#' lines(fit_ols_harmR_missing$fitted, lty = 4, col = "green")
#' 
#' @return A list containing:
#' \item{a.coef}{a numeric vector with estimates of cosine coefficients}
#' \item{b.coef}{a numeric vector with estimates of sine coefficients}
#' \item{amplitude}{a numeric vector with amplitude estimates.}
#' \item{phase}{a numeric vector with phase estimates.}
#' \item{fitted}{a numeric vector with fitted values via harmonic regression.}
#' 
#' @references Roerink, G.J., Menenti, M., Verhoef, W. (2000).
#' \emph{Reconstructing cloudfree NDVI composites using Fourier analysis of time series}, 
#' Int. J. Remote Sensing, \bold{21(9)}, 1911--1917.
#' 
#' @references Jakubauskas, M., Legates, D., Kastens, J. (2001).
#' \emph{Harmonic analysis of time-series AVHRR NDVI data},
#' Photogrammetric Engineering and Remote Sensing, \bold{67(4)}, 461--470.
#'   
#' @references The Matlab implementation of HANTS can be found 
#' \href{https://nl.mathworks.com/matlabcentral/fileexchange/38841-matlab-implementation-of-harmonic-analysis-of-time-series-hants}{here}.
#' 
haRmonics <- function(y, method = c("ols_harmR", "wls_harmR", "hants"), sigma=NULL,
                      ts = 1:length(y), lenPeriod = length(y), numFreq, HiLo = c("Hi", "Lo"), 
                      low, high, fitErrorTol, degreeOverDeter, delta){
  if(missing(y)){
    stop("y must be provided")
  }
  
  if(missing(numFreq) | numFreq <= 0){
    stop("A nonnegative numFreq must be provided")
  }
  
  if(missing(delta) | delta < 0){
    stop("A positive delta must be provided")
  }
  
  method <- match.arg(method)
  
  if(method == "ols_harmR"){
    output <- harmonicR(y = y, sigma = sigma, lenPeriod = lenPeriod, numFreq = numFreq, 
                        delta = delta)
  }
  
  if(method == "wls_harmR"){
    if(is.null(sigma)){
      stop("sigma must be provided. Check out geoTS::getSigma() for an alternative")
    }
    
    if( length(sigma) != lenPeriod ){
      stop("sigma must be a vector of length 'lenBasePeriod'. Check out geoTS::getSigma() for an alternative")
    }
    
    output <- harmonicR(y = y, sigma = sigma, lenPeriod = lenPeriod, 
                        numFreq = numFreq, delta = delta)
  }
  
  if(method == "hants"){
    if(missing(low) | missing(high)){
      stop("low and high must be provided")
    }
    
    if(missing(fitErrorTol)){
      stop("fitErrorTol must be provided")
    }
    
    if(missing(degreeOverDeter)){
      stop("degreeOverDeter must be provided")
    }
    
    output <- hants(y = y, ts = ts, lenPeriod = lenPeriod, numFreq = numFreq, 
                    HiLo = HiLo, low = low, high = high, fitErrorTol = fitErrorTol, 
                    degreeOverDeter = degreeOverDeter, delta = delta)
  }
  
output
}
# -----------------------------------------------------------------------------

