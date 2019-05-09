#' Get maximum lag of missing values
#' 
#' This function computes the maximum amount of consecutive missing values in a vector.
#' This quantity is also know as maximum lag, run, or record, and can be used as a
#' rough estimate of the quality of a dataset.
#' 
#' @param     x numeric vector
#' @param  type character specifying the type of missing value to consider. Default is 
#'              \code{type = "NA"}; when \code{type == "numeric"}, \code{value} must be provided
#' @param value numeric giving a figure to be considered as a missing value; often as part
#'              of a pre-processing missing values in a dataset (vector, time series, etc.)
#'              are fill in with pre-established values
#' @export
#' 
#' @examples
#' v <- c(NA, 0.12, 0.58, 0.75, NA, NA, NA, 0.46, 0.97, 0.39,
#'        NA, 0.13, 0.46, 0.95, 0.30, 0.98, 0.23, 0.98,
#'        0.68, NA, NA, NA, NA, NA, 0.11, 0.10, 0.79, 0.46, 0.27,
#'        0.44, 0.93, 0.20, 0.44, 0.66, 0.11, 0.88)
#' maxLagMissVal(x=v, type="NA")
#' 
#' w <- c(23,3,14,3,8,3,3,3,3,3,3,3,10,14,15,3,10,3,3,6)
#' maxLagMissVal(x=w, type = "numeric", value = 3)
#'        
#' @seealso \code{\link[base]{rle}}
#' 
#' @return A list containing:
#' \item{maxLag}{numeric giving the maximum lag of missing values in \code{x}}
#' \item{x}{numeric vector with the original data}
#' \item{value}{a numeric when \code{type == numeric}, \code{NA} otherwise}
#' 
maxLagMissVal <- function(x, type = c("NA", "numeric"), value){
  type <- match.arg(type)
  
  flag <- F
  
  if(type == "numeric"){
    if( missing(value) ){
      stop("When type is numeric, value must be passed")
    } else {
      flag <- T
      xCopy <- x
      x[x==value] <- NA
    }
  }
  
  is.na.rle <- rle(is.na(x))
  ind <- 1:length(is.na.rle$values)
  temp <- ind[is.na.rle$values == T]
  output <- ifelse(length(temp) != 0, max(is.na.rle$lengths[temp]), 0)
  
  if(flag == F){
    xOutput <- x
    valueOutput <- NA
  } else {
    xOutput <- xCopy
    valueOutput <- value
  }
  
  list(maxLag = output, x = xOutput, value = valueOutput)
}
