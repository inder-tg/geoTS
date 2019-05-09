#' Transfer values of a binary into a raster file
#'
#' Get the values of a binary file and transfer them into a raster file. All formats
#' considered in \code{\link[raster]{writeRaster}} are allowed.
#' 
#' @param inputPath  character with full path name of input file(s). A list of files is also allowed.
#' @param outputPath character with full path name (where the raster files will be saved). 
#' @param masterFile character with full path name of a raster file with basic specifications (such
#'                   as extent and projections) to be used in the transfer process.
#' @param what       See \code{\link[base]{readBin}}. Default \code{integer()}.
#' @param signed     See \code{\link[base]{readBin}}. Default \code{T}.
#' @param endian     See \code{\link[base]{readBin}}. Default \code{"little"}.
#' @param size       integer, number of bytes per element in the byte stream, default 2. See \code{\link[base]{readBin}}
#' @param format     character, output file type. See \code{\link[raster]{writeFormats}}.
#' @param dataType   character, output data type. See \code{\link[raster]{dataType}}.
#' @param overwrite  logical, default \code{T}, should the resulting raster be overwritten.
#' 
#' @export
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom tools file_path_sans_ext
#' @importFrom raster nrow
#' @importFrom raster ncol
#' @importFrom raster nlayers
#' @importFrom raster raster
#' @importFrom raster writeRaster
#' @importFrom raster rasterToPoints
#' @importFrom raster projection
#' 
#' @return At the designated path (\code{outputFile}) the user will find an \code{RData} file
#' 
transfer_bin_raster <- function(inputPath, outputPath, masterFile, what = integer(),
                                signed = T, endian = "little", size = 2,
                                format = "GTiff", dataType = "INT2S", overwrite = T){
  
  if(missing(inputPath)){
    stop("inputPath must be provided")
  }
  
  if(missing(outputPath)){
    stop("outputPath must be provided")
  }
  
  # if(length(inputPath) == 1){
  #   listFilesBin <- inputPath
  #   listFilesBinTemp <- inputPath
  # } 
  
  # inputPath = "/home/Desktop/posgradosUNAM/GEOGRAFIA/2018_TS/geoTS/data/site_BajaNorte_BIN"
  # outputPath = "/home/Desktop/posgradosUNAM/GEOGRAFIA/2018_TS/geoTS/data/site_BajaNorte_TIF"
  # dirBajaNorte <- "/home/itecuapetla/Desktop/posgradosUNAM/R/timeSeries/site_BajaNorte" #  paste(rootDir, "/data/site_QRoo", sep = "" )
  # masterFile <- raster(paste0(dirBajaNorte, "/master_BajaNorte.tif"))
  # masterFile = masterFile
  # size = 2
  
  # if(length(inputPath) > 1){
    listFilesBin <- list.files(path = inputPath, pattern = ".bin", full.names = T)
    listFilesBinTemp <- list.files(path = inputPath, pattern = ".bin")
  # }
  
  if(missing(masterFile)){
    stop("masterFile must be provided")
  } else {
    masterFile <- raster(masterFile)
  }
  
  nRow <- nrow(masterFile)
  nCol <- ncol(masterFile)
  
  x <- 1:nCol
  y <- 1:nRow
  
  cat("Started at:", as.character(Sys.time()[1]), "\n")
  
  pBar <- txtProgressBar(min = 0, max = length(listFilesBin), style = 3)
  for(i in 1:length(listFilesBin)) {
    imagen <- matrix(NA, ncol = nCol, nrow = nRow, byrow = TRUE)
    
    to.read <- file(listFilesBin[i], "rb")
    for(r in 1:nRow){
      imagen[r,] <- readBin(to.read, what = what, n = nCol, signed = T, size = size, endian = endian)
    }
    close(to.read)
    
    writeRaster(matrixToRaster(t(imagen), masterFile), 
                filename = paste0(outputPath, "/", file_path_sans_ext(listFilesBinTemp[i])), 
                format = format, datatype = dataType, overwrite = overwrite)
    setTxtProgressBar(pBar, i)
  }
  close(pBar)
  
  cat("Finished at:", as.character(Sys.time()[1]), "\n")
}
#
rotate <- function(x) t(apply(x, 2, rev))
#
matrixToRaster <- function(matrix, raster){
  rasterTable <- data.frame(rasterToPoints(raster)) 
  
  df <- data.frame(x = rasterTable$x, y = rasterTable$y, values = c(matrix))
  
  sp::coordinates(df) <- ~ x + y
  
  sp::gridded(df) <- TRUE
  
  raster_df <- raster(df)
  
  raster::projection(raster_df) <- projection(raster)
  
raster_df
}
