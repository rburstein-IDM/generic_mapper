## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Utility functions for generic_mapper
## ###########################################################

message(' <<< LOADING UTILITIES >>> ')

#' Insert Raster
#' Helper function stolen from seegMBG. Replace values in a raster with those in a column of a matrix.
#' @param raster Raster with target resolution and dimensions
#' @param new_vals New data matrix, each column representing a layer in the resulting raster brick
#'
#' @return a new raster or rasterBrick with the values supplies in new_vals
insertRaster <- function (raster, new_vals) {
  cellIdx <- function (x) 1:ncell(x)
  idx <- cellIdx(raster)
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }
  return (raster_new)
}



