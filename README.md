# geoTS

## Overview 
geoTS provides Time Series Methods for Analyzing Remote Sensing produtcs
such as

  - **split_replace**: Splits a ```Raster*``` object into smaller chunks and allows to replace cell values
  
  - **matrixToRaster**: Transforms a ```matrix``` into a ```RasterLayer``` object
  
  - **haRmonics**: Fits harmonic regression models (with homo- and heteroscedastic errors), 
  that is, computes amplitudes and phase angles. The popular method [HANTS](https://www.asprs.org/pers-archives-of-the-past?a_tag=2001&a_category=Apr&submit=GO)
  by Jakubauskas, Legates and Kastens is also implemented.
  
  - **hetervar**: Estimation of heteroscedastic variance for remotely-sensed time series

## Installation
``` r
# geoTS is available through CRAN, hence, you can type in
install.packages("geoTS")
```
