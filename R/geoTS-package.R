#' Methods for Handling and Analyzing Time Series of Satellite Images
#' 
#' We provide tools for handling time series of satellite images as well
#' as some statistical methods for spatio-temporal analysis
#' 
#' @name geoTS-package
#' @author Tecuapetla-Gomez, I. \email{itecuapetla@@conabio.gob.mx}
#' 
#' @section Tools for handling time series of satellite images:
#' \code{\link{transfer_bin_raster}} transfers data from images originally
#' recorded in a binary format to images in any of the formats
#' allowed by the \code{\link{raster}} package. Similarly, 
#' \code{\link{transfer_raster_RData}} extracts the entries (numbers) of
#' images originally recorded as a \code{\link{tiff}} file, virtually storages them 
#' in an \code{\link{array}} object and, finally, this array is saved in an RData file.
#' \code{\link{split_replace}} allows us to split Raster* objects, which can
#' be arguably large, into smaller chunks. These chunks can be saved
#' in any of the formats allowed by \code{\link{writeRaster}}. Often, satellite
#' images come with missing values (or fill values assigned by other computer
#' programs), \code{\link{split_replace}} allows to replace these values by
#' values of users' convenience; see also \code{\link{reclassify}}.
#' 
#' @section Methods for analyzing time series of satellite images:
#' \code{\link{haRmonics}} allows us to fit harmonic regression models
#' to numeric vectors; the method \code{hants} is based on \cite{Roerink et al. (2000)}
#' whereas the method \code{ols_harmR} is based on \cite{Jakubauskas et al. (2001)}.
#' The \code{wls_harmR} is the weighted least squares method which requires pre-estimation
#' of heteroscedastic variance; \code{hetervar} allows for heteroscedastic variance
#' estimation for numeric vectors extracted from time series of satellite imagery. 
#' 
#' @keywords package
#' 
#' @references Roerink, G.J., Menenti, M., Verhoef, W. (2000).
#' \emph{Reconstructing clodfree NDVI composites using Fourier analysis of time series}, 
#' Int. J. Remote Sensing, \bold{21(9)}, 1911--1917.
#' 
#' @references Jakubauskas, M., Legates, D., Kastens, J. (2001).
#' \emph{Harmonic analysis of time-series AVHRR NDVI data},
#' Photogrammetric Engineering and Remote Sensing, \bold{67(4)}, 461--470.
#'   
#' @references The Matlab implementation of HANTS can be found 
#' \href{https://nl.mathworks.com/matlabcentral/fileexchange/38841-matlab-implementation-of-harmonic-analysis-of-time-series-hants}{here}. 
#' @docType package
NULL
# -----------------------------------------------------------------------------