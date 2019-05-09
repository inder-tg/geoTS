#' Split a Raster* object and replace cell values (optional)
#'
#' This function will split a Raster* object into \code{partPerSide^2} parts. Additionally,
#' it allows to replace cell values (\code{valToReplace}) within Raster* object by another
#' value of user's choice (\code{replacedBy}). When \code{save = T}, the resulting \code{cellsToProcess} 
#' Raster* objects are saved in directory \code{outputPath}.
#' 
#' @param raster         Raster* object
#' @param partPerSide    numeric indicating the number of cells in which \code{raster} will be split
#' @param save           logical, should the output be saved, default is \code{TRUE}
#' @param replace        logical, default \code{FALSE}, when \code{TRUE}, \code{valToReplace} and \code{replacedBy} must by specified
#' @param valToReplace   indicates a value to be replaced across \code{raster} cells
#' @param replacedBy     indicates the value by which \code{valToReplace} is replaced
#' @param cellsToProcess numeric vector indicating which of the \code{partPerSide^2} should be processed/saved
#' @param format         character, output file type, default \code{"GTiff"}. See \code{\link[raster]{writeFormats}}
#' @param dataType       character, output data type. See \code{\link[raster]{dataType}}
#' @param outputPath     character with full path name where the resulting Raster* objects will be saved 
#' @param name           character with the name to assign to final products
#' @param ...            additional arguments used by \code{\link[raster]{writeRaster}} 
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
#' 
#' @details Before processing any of the \code{cellsToProcess} the temporary raster 
#' directory is re-directed. Basically, prior to process the i-th cell, 
#' at \code{outputPath} a new subdirectory is created, which, in turn, is erased 
#' automatically once the i-th cell has been processed. As a result of multiple testing 
#' we found that this measure avoids memory overflow.
#' 
#' @seealso \code{\link[raster]{writeRaster}}, \code{\link[raster]{aggregate}}, 
#' \code{\link[raster]{rasterOptions}}
#' 
#' @return At \code{outputPath} the user will find \code{length(cellsToProcess)} files
#' 
split_replace_raster <- function(raster, partPerSide, save = T, replace = F, 
                                 valToReplace, replacedBy, dataType,
                                 cellsToProcess = 1:(partPerSide^2), format = "GTiff",
                                 outputPath, name, ...){
  
  if(missing(raster)){
    stop("raster must be provided")
  }
  
  if(missing(partPerSide)){
    stop("partPerSide must be provided")
  }
  
  if(missing(outputPath)){
    stop("outputPath must be provided")
  }
  
  if(missing(name)){
    stop("outputPath must be provided")
  }
  
  # listFilesTEMP <- list.files(path = "/home/itecuapetla/Desktop/proyectoImputation/sim_dataCube/site_BajaNorte_TIF", full.names = T)
  # 
  # listFiles <- listFilesTEMP[-c(1:23)]
  # 
  # stackBajaNorte <- stack(listFiles)
  # 
  # raster <- stackBajaNorte
  # partPerSide <- 25
  # outputFile <- "/home/itecuapetla/Desktop/proyectoImputation/sim_dataCube/TIFs"
  # name <- "bajaNorte"
  
  if(replace){
    if(missing(valToReplace) | missing(replacedBy)){
      stop("When replace = T, valToReplace and replacedBy must be specified")
    }
  }
  
  h <- ceiling(ncol(raster)/partPerSide)
  
  v <- ceiling(nrow(raster)/partPerSide)
  
  agg <- getAggregate(raster, h = h, v = v)
  
  agg[] <- 1:ncell(agg)
  
  agg_poly <- rasterToPolygons(agg)
  
  names(agg_poly) <- "polis"
  
  cat("Started at:", as.character(Sys.time()[1]), "\n")
  
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
    
    for(j in 1:nlayers(raster)){
      aux_r <- crop(raster[[j]], extent_polygon)
      temp_r <- addLayer(temp_r, aux_r)
    }
    
    if(replace == T){
      temp_r <- reclassify( temp_r, cbind( valToReplace, replacedBy ) )
      # values(temp_r_copy)[values(temp_r_copy) == valToReplace] <- replacedBy
    }
    
    if(save == T){
      writeRaster(temp_r, filename = paste0(outputPath, "/", name, "_", i),
                  format = format, datatype = dataType, overwrite = TRUE, ...)  
    }
    
    unlink(paste0(outputPath, "/temp_", name, "_", i), recursive = TRUE)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  cat("Finished at:", as.character(Sys.time()[1]), "\n")
  
}
#
getAggregate <- function(raster, h = 1, v = 1){
  
  if( nlayers(raster) > 1 ){
    raster <- raster::subset(raster,1)
  }
  
  agg <- aggregate(raster, fact = c(h,v))
  
  agg
}

