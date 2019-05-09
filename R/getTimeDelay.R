# @export
# @importFrom Matrix Diagonal

getTimeDelay <- function(timeSeriesX, timeSeriesY, 
                         method = c("lasso", "sqrtLasso", "pearson"),
                         corMethod = c("productMean", "covariance", "correlation", "circularCov", "circularCor"),
                         # alpha = 0.05,
                         searchGrid = unique(round((length(timeSeriesX) - length(timeSeriesX) *.4):(length(timeSeriesX) + length(timeSeriesX) * .4))),
                         # missingValuesX, missingValuesY, 
                         cvLambda = F, prob = 0.9, coefId = 1){
  

  if( length(timeSeriesX) != length(timeSeriesY) ){
    stop("length of time series must be equal")
  }
  
  r <- range(searchGrid)
  widthSearchGrid <- (r[2] - r[1]) / 2
  
  X <- timeSeriesX
  Y <- timeSeriesY
  
  funCor <- numeric(2 * length(timeSeriesX) - 1)
  
  method <- match.arg(method)
  
  corMethod <- match.arg(corMethod)
  
  if(method != "pearson"){
    # below commented on Oct 23
    # if( !missing(missingValuesX) | !missing(missingValuesY) ){
    #   Psi <- getPsiMissingValues(timeSeriesY, missingValuesY)
    #   x <- timeSeriesX
    #   
    #   if( length(missingValuesX) != 0 ){
    #     x[missingValuesX] <- 0
    #   }
    #   
    #   # s_lasso1 <- sXcorr( x = as.numeric(x), Psi = Psi, Phi = Diagonal( length(timeSeriesX) ) )
    #   sparseSol <- getSparse( x = as.numeric(x), Psi = Psi, Phi = Diagonal( length(timeSeriesX) ),
    #                           method = method, cvLambda = cvLambda )
    # } else {
    # above commented on Oct 23
      Psi <- getPsi(timeSeriesY)
      
      sparseSol <- getSparse( x = as.numeric(timeSeriesX), Psi = Psi, Phi = Diagonal( length(timeSeriesX) ),
                             method = method, cvLambda = cvLambda )
      # s_lasso1 <- sXcorr( x = as.numeric(timeSeriesX), Psi = Psi, Phi = Diagonal( length(timeSeriesX) ) ) # this lines was commented before Oct 23
    # }
    
    if(method == "lasso" & cvLambda == F){
      lambda_s <- quantile( sparseSol$lambda, probs = prob )
      temp <- as.numeric( coef(sparseSol, s = lambda_s) )
    }
    
    if(method == "lasso" & cvLambda == T){
      temp <- as.numeric( coef(sparseSol) )
    }
    
    if(method == "sqrtLasso"){
      temp <- getBeta(slim = sparseSol, id = coefId)
    }
    
    # lambda_s <- quantile( s_lasso1$lambda, probs = prob )
    # 
    # temp <- as.numeric( coef(s_lasso1, s = lambda_s) )
    
    X <- Psi %*% temp[2:length(temp)]
  }
  
  funCor <- getFunCor(X = X, Y = Y, corMethod = corMethod)
                      # missingValuesX = missingValuesX, missingValuesY = missingValuesY) # commmented on Oct 23, 2018
  
  # pValue <- getPValue(X = X, Y = Y, 
  #                     h = -widthSearchGrid + which.max(funCor[searchGrid]^2) - 1,
  #                     corMethod = corMethod)
  
  list(tau = -widthSearchGrid + which.max(funCor[searchGrid]^2) - 1, # pValue = pValue,
       funCor = funCor, searchGrid = searchGrid, xLASSO = X, 
       lambda = ifelse(method == "pearson", NA, sparseSol$lambda))
}
