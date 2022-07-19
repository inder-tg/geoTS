#' Splits a Raster* object into smaller chunks and allows to replace cell values 
#'
#' This function will split a Raster* object into smaller chunks. The size of these chunks (number of cells) 
#' is controlled by \code{partPerSide}, \code{h} or \code{v}. Additionally, it allows to replace cell values (\code{valToReplace}) 
#' within Raster* object by another value of user's choice (\code{replacedBy}). When \code{save = TRUE}, 
#' the resulting \code{cellsToProcess} Raster* objects are saved in directory \code{outputPath}.
#' 
#' @param raster             Raster* object.
#' @param partPerSide        integer indicating number of cells in which \code{raster} will be split 
#'                           in each direction (horizontally and vertically). Use when \code{nrow(raster)} and 
#'                           \code{ncol(raster)} are multiples of \code{partPerSide}.
#' @param h                  integer indicating number of horizontal cells in which \code{raster} will be split.
#' @param v                  integer indicating number of vertical cells in which \code{raster} will be split.
#' @param outputPath         character with full path name where the resulting Raster* objects will be saved.
#' @param name               character with the name to assign to final products.
#' @param save               logical, should the output be saved, default is \code{TRUE}.
#' @param replace            logical, default \code{FALSE}, when \code{TRUE}, \code{valToReplace} and \code{replacedBy} must by specified.
#' @param valToReplace       indicates a value to be replaced across \code{raster} cells.
#' @param replacedBy         indicates the value by which \code{valToReplace} is replaced.
#' @param dataType           character, output data type. See \code{\link[raster]{dataType}}.
#' @param format             character, output file type, default \code{"GTiff"}. See \code{\link[raster]{writeFormats}}.
#' @param parallelProcessing logical, default \code{FALSE}, when \code{TRUE} raster splitting is done in parallel. See
#'                           \code{details}.
#' @param numCores           numeric indicating the number of cores used in parallel processing.
#' @param cellsToProcess     numeric vector indicating which smaller cells should be processed/saved. See \code{details}.
#' @param ...                additional arguments used by \code{\link[raster]{writeRaster}}. 
#' 
#' @export
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom raster nrow
#' @importFrom raster ncol
#' @importFrom raster ncell
#' @importFrom raster nlayers
#' @importFrom raster raster
#' @importFrom raster dataType
#' @importFrom raster stack
#' @importFrom raster reclassify
#' @importFrom raster crop
#' @importFrom raster rasterToPolygons
#' @importFrom raster writeRaster
#' @importFrom raster brick
#' @importFrom raster addLayer
#' @importFrom raster extent
#' @importFrom raster rasterOptions
#' @importFrom raster aggregate
#' @importFrom doParallel registerDoParallel 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel stopCluster
#' @importFrom utils globalVariables
#' 
#' @details Before processing any of the \code{cellsToProcess} the temporary raster 
#' directory is re-directed. Basically, prior to process the i-th cell, 
#' at \code{outputPath} a new subdirectory is created, which, in turn, is erased 
#' automatically once the i-th cell has been processed. As a result of several tests 
#' we found that this measure avoids memory overflow.
#' 
#' When \code{partPerSide} is used, \code{cellsToProcess = 1:(partPerSide^2)}. When \code{h}
#' and \code{v} are used, \code{cellsToProcess = 1:(ncells(raster)/(h*v))}. Since the code
#' assumes that \code{nrow(raster)} and \code{ncol(raster)} are multiples of \code{partPerSide}
#' or \code{h} and \code{v}, respectively, the user must be careful when selecting these
#' parameters.
#' 
#' For \code{parallelProcessing} the backend \code{\link[doParallel]{doParallel}} is employed.
#' 
#' @seealso \code{\link[raster]{writeRaster}}, \code{\link[raster]{aggregate}}, 
#' \code{\link[raster]{rasterOptions}}
#' 
#' @return At \code{outputPath} the user will find \code{length(cellsToProcess)} Raster* files
#' 
split_replace <- function(raster, partPerSide, h, v, outputPath, name, save = TRUE, 
                          replace = FALSE, valToReplace, replacedBy, 
                          dataType, format = "GTiff", 
                          parallelProcessing = FALSE, numCores = 20,
                          cellsToProcess, ...){
  
  if(missing(raster)){
    stop("raster must be provided")
  }
  
  if(missing(outputPath)){
    stop("outputPath must be provided")
  }
  
  if(missing(name)){
    stop("name must be provided")
  }
  
  if(replace){
    if(missing(valToReplace) | missing(replacedBy)){
      stop("When replace = TRUE, valToReplace and replacedBy must be specified")
    }
  }
  
  if(missing(partPerSide)){
    if(missing(h) | missing(v)){
      stop("h and v must be provided")
    } else {
      if(missing(cellsToProcess)){
        cellsToProcess <- 1:(ncell(raster)/(h*v))
      }
    }
    # stop("partPerSide must be provided")
  } else {
    h <- ceiling(ncol(raster)/partPerSide)
    
    v <- ceiling(nrow(raster)/partPerSide)
    
    if(missing(cellsToProcess)){
      cellsToProcess <- 1:(partPerSide^2)
    }
  }
  
  agg <- getAggregate(raster, h = h, v = v)
  
  agg[] <- 1:ncell(agg)
  
  agg_poly <- rasterToPolygons(agg)
  
  names(agg_poly) <- "polis"
  
  message(paste0("Started at: ", as.character(Sys.time()[1])))
  
  pb <- txtProgressBar(min = 0, max = length(cellsToProcess), style = 3)
  for(i in cellsToProcess){
    Sys.sleep(0.1)
    
    extent_polygon <- extent(agg_poly[agg_poly$polis == i,])
    
    # --- create temp dir
    dir.create(path = paste0(outputPath, "/temp_", name, "_", i),
               showWarnings = F)
    rasterOptions(tmpdir = paste0(outputPath, "/temp_", name, "_", i))
    # ---
    
    temp_r <- brick()  
    
    if( parallelProcessing ){

      output <- getParallelOutput(raster = raster, numCores = numCores, 
                                  polygon_extent = extent_polygon)
      
      for(k in 1:nlayers(raster)){
        temp_r <- addLayer(temp_r, output[[k]])
      }
      
    } else {
      for(k in 1:nlayers(raster)){
        aux_r <- crop(raster[[k]], extent_polygon)
        temp_r <- addLayer(temp_r, aux_r)
      }
    }
    
    if(replace == TRUE){
      temp_r <- reclassify( temp_r, cbind( valToReplace, replacedBy ) )
    }
    
    if(save == TRUE){
      writeRaster(temp_r, filename = paste0(outputPath, "/", name, "_", i),
                  format = format, datatype = dataType, overwrite = TRUE, ...)  
    }
    
    unlink(paste0(outputPath, "/temp_", name, "_", i), recursive = TRUE)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  message(paste0("Finished at: ", as.character(Sys.time()[1])))
  
}
#
getAggregate <- function(raster, h = 1, v = 1){
  
  if( nlayers(raster) > 1 ){
    raster <- raster::subset(raster,1)
  }
  
  agg <- aggregate(raster, fact = c(h,v))
  
  agg
}

