#' Produce a stars object with zero value on a fixed resolution
#'
#' @param data spatial feature of the relevant locations
#' @param resolution cell size
#' @param values values to populate the stars values with, default to 0
#'
#' @export
stars_ref <- function(data, resolution, values = 0L) {

  stars::st_as_stars(sf::st_bbox(data), dx = resolution, dy = resolution, values = values)
}

#' Creates identifiers for (x,y) coordinates.
#' Used to work with crs 3035
#'
#' @param x matrix of coordinates (first column = x, second = y). Or a vector of x coordinates.
#' @param y default to NULL (if x is a matrix). Or a vector of y coordinates.
#' @param resolution default to 200. Size of raster cells which are identified.
#'
#' @export
coord2idINS <- function(x, y = NULL, resolution = 200) {

  if (!is.null(y)) x <- cbind(x, y)
  resolution <- round(resolution)

  res <- apply(x, MARGIN = 1:2, FUN = function(z) formatC(floor(z / resolution) * resolution, format = "d"))
  res <- paste0("r", resolution, "N", res[, 2], "E", res[, 1])

  nas <- apply(x, MARGIN = 1, FUN = anyNA)
  res[nas] <- NA

  res
}

#' Retrieves (x, y) coordinates from idINS.
#' Used to work with crs 3035
#'
#' @param ids character vector of idINS
#' @param stop_if_res_not_cst default to TRUE. Refuse to convert id with different resolutions
#'
#' @export
idINS2coord <- function(ids, stop_if_res_not_cst = TRUE) {

  ids <- strsplit(as.character(ids), "[rNE]")
  resolution <- lapply(ids, FUN = function(liste) liste[2]) |> unlist() |> as.numeric()
  if (stop_if_res_not_cst) stopifnot(length(unique(resolution)) == 1L)

  x <- lapply(ids, FUN = function(liste) liste[4]) |> unlist() |> as.numeric()
  y <- lapply(ids, FUN = function(liste) liste[3]) |> unlist() |> as.numeric()

  res <- matrix(c(x + resolution/2, y + resolution/2), ncol = 2)
  colnames(res) <- c("x", "y")

  res
}

#' Transform a data.table with idINS to a star object.
#' Assumes crs 3035
#'
#' @param data a data.table with an idINS column
#' @param idINS name of the idINS column. Default to idINS.
#'
#' @import data.table
#'
#' @export
dt2stars <- function(data, idINS = "idINS")
{
  if(!is.data.table(data)) setDT(data)
  xy <- idINS2coord(data[[idINS]]) |> as.data.table()

  data <- cbind(data[, .SD, .SDcol = !idINS], xy)

  #sf::st_as_sf(data, coords = c("x", "y"), crs = st_crs(3035))
  data <- stars::st_as_stars(data, dims = c("x", "y"))
  st_crs(data) <- sf::st_crs(3035)

  data
}

