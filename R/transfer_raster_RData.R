#' Transfer the values of a RasterStack into an RData file
#'
#' Get the values of a Raster* and transfer them into an RData file in an 
#' \code{\link[base]{array}} format which allows for compatibility with multiple 
#' R functions as well as great portability.
#' 
#' @param inputFile  a character vector with full path name of input file
#' @param outputFile a character vector with full path name (where the \code{RData} file will be saved). 
#'                     Do not include the extension .RData
#' @param vmode      a character specifying the type of virtual storage mode \code{\link[ff]{vmode}} 
#'                   needed. Only \code{integer}, \code{single} and \code{double} are
#'                   allowed
#' 
#' @export
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom raster nrow
#' @importFrom raster ncol
#' @importFrom raster nlayers
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom ff delete
#' 
#' @details Prior to embark the user in a transfer that may not be successful due to the
#' lack of RAM, \code{transfer_raster_RData} provides an estimate of the amount of bytes to be used
#' in the transfer process. The estimate is obtained by multiplying the number of rows by the number of 
#' columns of the first layer of \code{inputFile} by its number of layers by the amount of 
#' bites used by \code{vmode} (32-bit float for \code{integer} or \code{single} and 64-bit float for \code{double}). 
#' Should the user decide not to continue with the importation \code{transfer_raster_RData} returns the message \code{"Did not transfer anything"}.
#' 
#' @seealso \code{\link[ff]{vmode}}
#' 
#' @return At the designated path (\code{outputFile}) the user will find an \code{.RData} file
#' 
transfer_raster_RData <- function(inputFile, outputFile, vmode = c("integer", "single", "double")){
  
  stack_TEMP <- stack(inputFile)
  
  # if(nlayers(stack_TEMP) == 1){
  #   stop("inputFile must be a RasterStack object")
  # }
  
  roughBites <- ncol(stack_TEMP) * nrow(stack_TEMP) * nlayers(stack_TEMP)
    
  vmode <- match.arg(vmode)
  
  if(vmode == "single" | vmode == "integer"){
    bytes <- roughBites * 4 / (2^20)
  }
  
  if(vmode == "double"){
    bytes <- roughBites * 8 / (2^20)
  }
  
  message(paste0("Need ", round(bytes, 3), " MB (approx.) to complete the transfer"))
  answer <- readline(prompt = "Would you like to continue (y/n): ")
  
  if( answer == "y" | answer == "Y" ){

    index <- ff::ff( vmode = vmode, dim = c(nrow(stack_TEMP), ncol(stack_TEMP), nlayers(stack_TEMP)),
                     filename = paste0(outputFile, ".ffdata") )
    
    # Sys.time()
    
    cat("Started at:", as.character(Sys.time()[1]), "\n")
    
    pBar <- txtProgressBar(min = 0, max = nlayers(stack_TEMP), style = 3)
    for(i in 1:nlayers(stack_TEMP)){
      Sys.sleep(0.1)
      index[, , i] <- stack_TEMP[[i]][]
      setTxtProgressBar(pBar, i)
    }
    close(pBar)
    
    temp <- 1:nlayers(stack_TEMP)
    
    index_array <- index[,,temp]  
    
    message("Please wait, your output file is being saved")
    
    save(index_array, file = paste0(outputFile, ".RData"))
    
    close(index)
    
    delete(index)
    
    # Sys.time()
    cat("Finished at:", as.character(Sys.time()[1]), "\n")
    
  } else {
    message("Did not transfer anything")
  }
}