#' Creates a RasterLayer object from a matrix
#'
#' Transforms a \code{matrix} into a \code{RasterLayer} object. 
#' 
#' @param matrix     a matrix object. See Details.
#' @param raster     a \code{RasterLayer} object whose extent and projection will be used to
#'                   create a raster from \code{matrix}. 
#' @param projection a character vector providing a coordinate reference system. 
#'                   Required when \code{ncol(matrix)=3}.  
#'               
#' @export
#' 
#' @importFrom raster raster
#' @importFrom raster rasterToPoints
#' @importFrom raster projection
#' 
#' @details When \code{ncol(matrix)=3}, this function assumes that the first two
#' columns of argument \code{matrix} provide coordinates to create a \code{RasterLayer},
#' hence argument \code{projection} must be provided. When argument \code{matrix} has 
#' only 2 columns, then the argument \code{raster} must be provided because its 
#' \code{\link[sp]{coordinates}} and \code{\link[raster]{projection}} will be used 
#' to rasterize \code{matrix}.
#' 
#' @note In previous versions, \code{raster} argument was written in capital letters.
#' 
#' @seealso \code{\link[raster]{Raster-class}} 
#' 
#' @return A \code{RasterLayer} 
#' 
matrixToRaster <- function(matrix, raster=NULL, projection=NULL){
  
  if(ncol(matrix)==3){
    if(is.null(projection)){
      stop('projection must be provided')
    } 
    
    x <- matrix[,1]
    y <- matrix[,2]
    values <- matrix[,3]
    PROJ <- projection
    
  } else {
    
    if(is.null(raster)){
      stop('raster must be provided')
    }
    
    rasterTable <- data.frame(rasterToPoints(raster))
    x <- rasterTable$x
    y <- rasterTable$y
    values <- c(matrix)
    PROJ <- raster::projection(raster)
  }
  
  df <- data.frame(x=x, y=y, values=values)
  
  sp::coordinates(df) <- ~ x + y
  
  sp::gridded(df) <- TRUE
  
  raster_df <- raster(df)
  
  raster::projection(raster_df) <- PROJ
  
raster_df
}
