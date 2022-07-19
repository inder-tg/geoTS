# -----------------------------------------------------------------------------
  hants <- function(y, ts = 1:length(y), lenPeriod = length(y), numFreq, 
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
    
    ang <- 2 * pi * (0:(lenPeriod-1)) / lenPeriod
    cs <- cos(ang)
    sn <- sin(ang)
    
    i <- 1:numFreq
    for( j in 1:numImages ){
      index <- 1 + (i * (ts[j]-1)) %% lenPeriod
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
    
    ready <- FALSE
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
    phase[phase < 0] <- phase[phase < 0] + 360
    phi[iFreq] <- phase
    
    list(a.coef = ra, b.coef = rb,
         amplitude = amp, phase = phi, fitted = ySmooth)  
  }
# -----------------------------------------------------------------------------
  harmonicR <- function(y, sigma, lenPeriod = length(y), numFreq, delta){
    
    numImages <- length(y)
    amp <- numeric(numFreq + 1)
    phi <- numeric(numFreq + 1)
    ySmooth <- numeric(numImages)
    
    numRows <- min(2 * numFreq + 1, numImages)
    degree <- 180 / pi
    
    mat <- matrix(0, nrow = numImages, ncol = min(2 * numFreq + 1, numImages))
    vecTs <- 2 * pi * (0:(lenPeriod-1)) / lenBasePeriod
    mat[,1] <- 1
    mat[, seq(2, 2 * numFreq, 2)] <- sapply(1:numFreq, function(s) cos( s * vecTs ))
    mat[, seq(2, 2 * numFreq, 2)+1] <- sapply(1:numFreq, function(s) sin( s * vecTs ))
    
    matCopy <- mat
    
    if( !is.null(sigma) ){
      y <- diag(1/sqrt(sigma)) %*% y
      mat <- diag(1/sqrt(sigma)) %*% mat
    }
    
    za <- t(mat) %*% ( y )
    
    A <- t(mat) %*% mat
    A <- A + diag(rep(1, numRows)) * delta
    A[1,1] <- A[1,1] - delta
    zr <- solve(A) %*% za
    
    ySmooth <- matCopy %*% zr
    
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
    
    list(a.coef = ra, b.coef = rb,
         amplitude = amp, phase = phi, fitted = as.numeric(ySmooth))
  }
# -----------------------------------------------------------------------------

globalVariables("j")
getParallelOutput <- function(raster, numCores, polygon_extent){
  closter <- parallel::makeCluster(spec = numCores, outfile = "")
  registerDoParallel(closter)
    
  output <- foreach(j = 1:nlayers(raster), .packages = c("raster")) %dopar% {
    return(crop(raster[[j]], polygon_extent))
  }
    
  stopCluster(closter)
    
output
}
  
# -----------------------------------------------------------------------------

# --- Added on July 18, 2022
  
# --- this function takes a time series -x- and returns a length(x)/lengthPeriod times lengthPeriod matrix
get_pixel_matrix <- function(x, lenPeriod=23){ 
  output <- matrix(nrow=length(x)/lenPeriod, ncol=lenPeriod)
  
  for(i in seq_len(nrow(output))){
    output[i,] <- x[((i-1) * lenPeriod + 1):(i * lenPeriod)]
  }
output
}

# -----------------------------------------------------------------------------
  

  