#' Estimation of heteroscedastic variance for remotely-sensed time series
#' 
#' @param         x numeric vector
#' @param lenPeriod numeric giving the number of observations per period.
#'                  Default, 23.
#' @param    method character specifying whether \code{standard} variance
#'                  or the median absolute deviation (\code{robust-mad})
#'                  method should be used
#' 
#' @export
#' 
#' @importFrom stats sd
#' @importFrom stats mad
#' @importFrom robustbase Qn
#' 
#' @details This function has been designed for time series of satellite imagery,
#' it is expected that \code{x} contains observations collected at pixels. 
#' That is, \code{length(x)} must be a multiple of \code{lenPeriod}. Default of 
#' \code{lenPeriod} corresponds to the temporal resolution of a MODIS product.
#' 
#' Method \code{standard} invokes \code{\link[stats]{sd}} whereas \code{robust-mad}
#' uses the median absolute deviation of \code{\link[stats]{mad}} and \code{robust-Qn}
#' utilizes the robust scale estimator implemented in \code{\link[robustbase]{Qn}}.
#' 
#' @return A numeric vector of length \code{lenPeriod}
#' 
hetervar <- function(x, lenPeriod=23, method=c("standard", "robust-mad", "robust-Qn")){
  
  if( length(x) %% lenPeriod != 0 ){
    stop("length(x) must be a multiple of lenPeriod")
  }
  
  method <- match.arg(method)
  
  m <- get_matrix_pixel(x, lenPeriod=lenPeriod)
  
  if( method=="standard" ){
    output <- apply(m, MARGIN = 2, FUN = sd)
  }
  
  if( method=="robust-mad" ){
    output <- apply(m, MARGIN = 2, FUN = mad)
  }
  
  if( method=="robust-Qn" ){
    output <- apply(m, MARGIN = 2, FUN = Qn)
  }
  
  # output <- apply(m, MARGIN = 2, FUN= ifelse(method=="standard", sd, 
  #                                            ifelse(method="robust-mad", mad, "robust-Qn")))

output
}