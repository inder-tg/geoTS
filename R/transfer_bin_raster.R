#' Transfer values from a binary image file to a raster file
#'
#' Get the values of a binary file (in integer format) and transfer them to a raster file. All formats
#' considered in \code{\link[raster]{writeRaster}} are allowed.
#' 
#' @param inputPath  character with full path name of input file(s).
#' @param outputPath character with full path name (where the raster files will be saved). 
#' @param master     character with full path name of a raster file; extent and projection
#'                   of this file are applied to this function output.        
#' @param what       See \code{\link[base]{readBin}}. Default \code{integer()}.
#' @param signed     See \code{\link[base]{readBin}}. Default \code{TRUE}.
#' @param endian     See \code{\link[base]{readBin}}. Default \code{"little"}.
#' @param size       integer, number of bytes per element in the byte stream, default 2. See \code{\link[base]{readBin}}.
#' @param format     character, output file type. See \code{\link[raster]{writeFormats}}.
#' @param dataType   character, output data type. See \code{\link[raster]{dataType}}.
#' @param overwrite  logical, default \code{TRUE}, should the resulting raster be overwritten.
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' inputPath = system.file("extdata", package = "geoTS")
#' masterFile = system.file("extdata", "master.tif", package = "geoTS") 
#' transfer_bin_raster(inputPath = inputPath, outputPath = inputPath, 
#'                     master = masterFile, what = integer(),
#'                     signed = TRUE, endian = "little", size = 2,
#'                     format = "GTiff", dataType = "INT2S", overwrite = TRUE)
#' }
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom tools file_path_sans_ext
#' @importFrom raster nrow
#' @importFrom raster ncol
#' @importFrom raster nlayers
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' 
#' @return At the designated path (\code{outputPath}) the user will find \code{TIF} file(s).
#' 
transfer_bin_raster <- function(inputPath, outputPath, master, what = integer(),
                                signed = TRUE, endian = "little", size = 2,
                                format = "GTiff", dataType = "INT2S", overwrite = TRUE){
  
  if(missing(inputPath) | missing(outputPath)){
    stop("inputPath and outputPath must be provided")
  }
  
  listFilesBin <- list.files(path = inputPath, pattern = ".bin", full.names = TRUE)
  listFilesBinTemp <- list.files(path = inputPath, pattern = ".bin")
  
  if(missing(master)){
    stop("masterFile must be provided")
  } else {
    masterFile <- raster(master)
  }
  
  nRow <- nrow(masterFile)
  nCol <- ncol(masterFile)
  
  x <- 1:nCol
  y <- 1:nRow
  
  message(paste0("Started at: ", as.character(Sys.time()[1])))
  
  pBar <- txtProgressBar(min = 0, max = length(listFilesBin), style = 3)
  for(i in 1:length(listFilesBin)) {
    imagen <- matrix(NA, ncol = nCol, nrow = nRow, byrow = TRUE)
    
    to.read <- file(listFilesBin[i], "rb")
    for(r in 1:nRow){
      imagen[r,] <- readBin(to.read, what = what, n = nCol, signed = TRUE, size = size, endian = endian)
    }
    close(to.read)
    
    writeRaster(matrixToRaster(t(imagen), masterFile), 
                filename = paste0(outputPath, "/", file_path_sans_ext(listFilesBinTemp[i])), 
                format = format, datatype = dataType, overwrite = overwrite)
    setTxtProgressBar(pBar, i)
  }
  close(pBar)
  
  message(paste0("Finished at: ", as.character(Sys.time()[1])))
}
#
rotate <- function(x) t(apply(x, 2, rev))
