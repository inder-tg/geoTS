# -----------------------------------------------------------------------------
standardize <- function(v){
  (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
corFun <- function(x, y, method = c("productMean", "covariance", "correlation", 
                                    "circularCov", "circularCor") ){
  method <- match.arg(method)
  output <- numeric(1)
  
  if(method == "productMean"){
    output <- mean( x * y, na.rm = TRUE )
  }
  
  if(method == "covariance" | method == "circularCov"){
    output <- cov(x, y, method = "pearson")  # mean( x * y, na.rm = TRUE )
  }
  
  if(method == "correlation" | method == "circularCor"){
    output <- cor(x, y, method = "pearson") # mean( x * y, na.rm = TRUE ) / ( sd(x, na.rm = TRUE) * sd(y, na.rm = TRUE) ) # cor(x, y, method = "pearson")
  }
  
  output
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
getMissingValues <- function(x){
  l <- length(x)
  indices <- 1:l
  indices[is.na(x)]
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# getFunCor <- function(X, Y, corMethod, missingValuesX, missingValuesY){ # commented on Oct 23, 2018
getFunCor <- function(X, Y, corMethod){
  funCor <- numeric(2 * length(X) - 1)
  
  for( h in -(length(X) - 1):(length(X) - 1) ){
    t <- h + length(X) 
    
    if( h < 0 ){
      
      if( corMethod == "circularCov" | corMethod == "circularCor" ){
        x <- c(X[(1 + abs(h)):length(X)], X[1:abs(h)])
        y <- Y
      } else {
        x <- X[(1 + abs(h)):length(X)]
        y <- Y[1:(length(X) - abs(h))]
      }
    }
    
    if( h == 0 ){
      x <- X
      y <- Y
    }
    
    if( h > 0 ){
      
      if( corMethod == "circularCov" | corMethod == "circularCor" ){
        x <- X
        y <- c(Y[(h + 1):length(X)], Y[1:h])
      } else {
        x <- X[1:(length(X) - h)]
        y <- Y[(h + 1):length(X)]
      }
    }
    
    # Below commented on Oct 23, 2018
    # if( !missing(missingValuesX) | !missing(missingValuesY) ){
    #   x <- x[-missingValuesY]
    #   y <- y[-missingValuesY]
    # }
    # Above commented on Oct 23, 2018
    
    funCor[t] <- corFun(x, y, method = corMethod) #cor(x, y, method = "pearson")
  }
  
  funCor
}
# -----------------------------------------------------------------------------
getPsi <- function(p){
  na <- length(p)
  Psi <- matrix(0, nrow = na, ncol = 2 * na - 1)
  
  for(i in 1:na){
    Psi[,i] <- c(rep(0, na-i), p[1:i])
  }
  
  for(i in (na + 1):(2 * na - 1)){
    Psi[,i] <- c( p[(i + 1 - na):na], rep(0, (i - na) ) )
  }
  
  Psi  
}
# -----------------------------------------------------------------------------
# Commented on Oct 23, 2018
# getPsiMissingValues <- function(p, missingValues){
#   # p <- timeSeriesY
#   
#   p[missingValues] <- 0
#   
#   Psi <- getPsi(p)
#   
#   # Psi[missingValues, ] <- 0
#   
#   Psi
# }
# -----------------------------------------------------------------------------

#' @importFrom glmnet glmnet
#' @importFrom glmnet cv.glmnet
#' @importFrom flare slim
#' 
getSparse <- function(x, Psi, Phi, method = c("lasso", "sqrtLasso"), cvLambda = F){
  
  n <- length(x)
  M <- Phi %*% Psi
  y <- Phi %*% x
  
  method <- match.arg(method)
  
  if(method == "lasso"){
    sparseSol <- glmnet(x = M, y = y, family = "gaussian", intercept = F)
  }
  
  if(method == "sqrtLasso"){
    sparseSol <- slim(X = M, Y = y)
  }
  
  if(cvLambda == T){
    crossval <-  cv.glmnet(x = M, y = y)
    penalty <- crossval$lambda.min # optimal lambda
    sparseSol <- glmnet(x = M, y = y, family = "gaussian", lambda = penalty,
                        intercept = F)
  }
  
  # crossval <-  cv.glmnet(x = M, y = y)
  # plot(crossval)
  # penalty <- crossval$lambda.min #optimal lambda
  # penalty #minimal shrinkage
  # fit1 <-glmnet(x = X, y = y, alpha = 1, lambda = penalty ) #estimate the model with that
  # coef(fit1)
  
  # lasso_sol <- glmnet(x = M, y = y, family = "gaussian", intercept = F) # , lambda = penalty
  
  # lasso_sol <- glmnet(x = M, y = y, intercept = F) # , lambda = penalty
  
sparseSol
}
# -----------------------------------------------------------------------------
getBeta <- function(slim, id){
  beta <- numeric(nrow(slim$beta)+1)
  
  beta[1] <- slim$intercept[id]
  
  beta[2:length(beta)] <- slim$beta[,id]
  
beta
}
# -----------------------------------------------------------------------------
# @export
# @importFrom stats cor.test
# 
getPValue <- function(X, Y, h, corMethod){
  
  if( is.na(h) ){
    testCor <- NA
  } else {
    if( h < 0 ){
      
      if( corMethod == "circularCov" | corMethod == "circularCor" ){
        x <- c(X[(1 + abs(h)):length(X)], X[1:abs(h)])
        y <- Y
      } else {
        x <- X[(1 + abs(h)):length(X)]
        y <- Y[1:(length(X) - abs(h))]
      }
    }
    
    if( h == 0 ){
      x <- X
      y <- Y
    }
    
    if( h > 0 ){
      
      if( corMethod == "circularCov" | corMethod == "circularCor" ){
        x <- X
        y <- c(Y[(h + 1):length(X)], Y[1:h])
      } else {
        x <- X[1:(length(X) - h)]
        y <- Y[(h + 1):length(X)]
      }
    }
    
    testCor <- cor.test(x, y, alternative = "two.sided", method = "pearson")
  }
  
testCor$p.value  
}

# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
#' @importFrom glmnet glmnet
#' 
sXcorr <- function(x, Psi, Phi){
  
  n <- length(x)
  M <- Phi %*% Psi
  y <- Phi %*% x
  
  # crossval <-  cv.glmnet(x = M, y = y)
  # plot(crossval)
  # penalty <- crossval$lambda.min #optimal lambda
  # penalty #minimal shrinkage
  # fit1 <-glmnet(x = X, y = y, alpha = 1, lambda = penalty ) #estimate the model with that
  # coef(fit1)
  
  lasso_sol <- glmnet(x = M, y = y, family = "gaussian", intercept = F) # , lambda = penalty
  
  # lasso_sol <- glmnet(x = M, y = y, intercept = F) # , lambda = penalty
  
  lasso_sol
}
# -----------------------------------------------------------------------------