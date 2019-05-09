#' Time delay estimation for discrete signals
#' 
#' Computes the lag that maximizes the Pearson correlation between two \code{numeric}
#' time series (emitted and received signal). When the emmitted signal is sparse, this function
#' allows us to project this signal onto the space of columns of a sparse matrix whose entries
#' are defined by the values of the received signal. To get this projection this function allows 
#' for Lasso or SQRT Lassoregression.
#' 
#' @param x    numeric (tentatively sparse) vector.
#' @param y    numeric vector
#' @param method         character string specifying the type of time delay estimation 
#'                       to be performed. When \code{x} is sparse, a sparse-reprojection
#'                       is available via \code{lasso} or \code{sqrtLasso} regression.                        
#' @param corMethod      character string specifying the type of correlation to be 
#'                       computed.
#' @param searchGrid     numeric vector giving a grid on points on which the correlation 
#'                       will be maximized.
#' @param cvLambda       logical (FALSE by default). Used when \code{method = lasso} as it
#'                       indicates whether the regularization parameter is chosen via cross-validation or not.                  
#' @param prob           numeric. Used when \code{method = lasso} and \code{cvLambda = F}
#'                       as it corresponds to the quantile of the set of \code{lambda}s
#'                       (regularization parameter) calculated at the moment of 
#'                       representing \code{x} as a highly-sparse vector.
#' @param coefId         numeric. Used when \code{method = sqrtLasso} as it indicates which
#'                       sqrt-lasso regression parameter \eqn{\beta} should be chosen. See details.                    
#' @export 
#' 
#' @importFrom stats coef
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom stats cor.test
#' @importFrom Matrix Diagonal
#' 
#' @details For details on the computaion of cross-correlation see Section 2 of \cite{Tecuapetla-Gómez (2019)}.
#' 
#' When \code{method = sqrtLasso} a sqrt-lasso regression is performed via the 
#' \code{\link[flare]{slim}} function. By default this function returns, \eqn{\beta},
#' a matrix of regression estimates whose columns correspond to regularization parameters.
#' Thus \code{coefId} specifies the column of \eqn{\beta} used to project (predict)
#' \code{x} as a highly-sparse vector. By default \code{ncol(beta)} = 5.
#' 
#' @examples 
#' # Basic signals
#' createSignalPeak <- function(leftValue, rightValue, cp1, cp2, value) {
#' function(t) ifelse(t < cp1, rep(leftValue, length(t)),
#'                    ifelse(t < cp2, rep(value, length(t)), 
#'                          rep(rightValue, length(t))))
#' }
#' stepF1 <- createSignalPeak(0, 0, 3, 6, 1)
#' stepF2 <- createSignalPeak(0, 0, 6, 9, 1)
#' n <- 100
#' ind <- seq(0, 10, length = n)
#' signal1 <- stepF1(ind)
#' signal2 <- stepF2(ind)
#' timeDelay(signal1, signal2, method = "pearson")
#' timeDelay(signal1, signal2, method = "pearson", corMethod = "covariance")
#' timeDelay(signal1, signal2, method = "pearson", corMethod = "correlation")
#' timeDelay(signal1, signal2, method = "pearson", corMethod = "circularCov")
#' timeDelay(signal1, signal2, method = "pearson", corMethod = "circularCor")
#' 
#' @seealso \code{\link[stats]{cor}}, \code{\link[stats]{cor.test}}, \code{\link[glmnet]{glmnet}}, \code{\link[flare]{slim}}, 
#' 
#' @return A list containing:
#' \item{tau}{a numeric giving the estimated time delay between \code{x} and
#' \code{y}.}
#' \item{funCor}{a numeric vector containing the correlation between \code{x}
#' and \code{y}. See details.}
#' \item{pValue}{a numeric giving the p-value associated with the test for correlation
#' between \code{x} and \code{y} at \code{tau}. Default parameters
#' of function \code{\link[stats]{cor.test}} are used.}
#' \item{searchGrid}{the grid of points on which \code{funCor} is maximized.}
#' \item{xLASSO}{numeric vector giving the representation of \code{x} as a 
#' highly-sparse vector.}
#' \item{lambda}{a numeric, giving the value of the regularization parameter used to
#' calculate \code{xLASSO}.}
#' 
#' @references Tecuapetla-Gómez, I. (2019). \emph{Time delay estimation in 
#' satellite imagery time series of precipitation and NDVI: Pearson's cross correlation
#' revisited.}
#'                       
timeDelay <- function(x, y, 
                      method = c("lasso", "sqrtLasso", "pearson"),
                      corMethod = c("productMean", "covariance","correlation","circularCov", "circularCor"),
                      # alpha = 0.05,
                      searchGrid = unique(round((length(x) - length(x) *.4):(length(x) + length(x) * .4))),
                      cvLambda = F, prob = 0.1, coefId = 1){
  
  method <- match.arg(method)
  
  corMethod <- match.arg(corMethod)
  
  tau <- numeric(1)
  
  missingValuesX <- numeric(0)
  missingValuesY <- numeric(0)
  
  missingValues_x <- getMissingValues(x)
  missingValues_y <- getMissingValues(y)
  
  if( length(missingValues_y) == 0 & length(missingValues_x) == 0 ){
    output <- getTimeDelay(x, y, method = method, corMethod = corMethod,
                           searchGrid = searchGrid, prob = prob, cvLambda = cvLambda,
                           coefId = coefId)
    pValue <- getPValue(X = x, Y = y, h = output$tau, corMethod = corMethod)
  } else {
    stop("This version does not allow missing values")
  }

  # Below commented on Oct 23, 2018  
  # if( length(missingValues_x) != 0 & length(missingValues_y) == 0  ){
  #   missingValuesX <- missingValues_x
  #   missingValuesY <- missingValues_x
  # }
  # 
  # if( length(missingValues_y) != 0 & length(missingValues_x) == 0  ){
  #   missingValuesY <- missingValues_y
  # }
  # 
  # if( length(missingValues_x) !=0 & length(missingValues_y) != 0 ){
  #   missingValuesX <- missingValues_x
  #   missingValuesY <- sort(union( missingValues_x, missingValues_y ))
  # }
  # 
  # if( length(missingValuesX) == 0 & length(missingValuesY) == 0 ){
  #   output <- getTimeDelay(x, y, method = method, corMethod = corMethod,
  #                            searchGrid = searchGrid, prob = prob, cvLambda = cvLambda,
  #                            coefId = coefId)
  # } else {
  #   output <- getTimeDelay(x, y, method = method, corMethod = corMethod,
  #                            missingValuesX = missingValuesX, missingValuesY = missingValuesY, 
  #                            searchGrid = searchGrid, prob = prob, cvLambda = cvLambda,
  #                            coefId = coefId)
  # }
  # Above commented on Oct 23, 2018  
  
list(tau = output$tau, funCor = output$funCor, pValue = pValue,
     searchGrid = output$searchGrid, xLASSO = output$xLASSO, lambda = output$lambda)
}
