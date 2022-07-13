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

#' Get the resolution of the data.frame or sf from the idINS column
#'
#' @param x a data.frame or sf object with a idINS column.
#' @param idINS name of the idINS column. default to idINS.
#'
#' @export
residINS <- function(x, idINS = "idINS") {

  if ("sf" %in% class(x)) x <- st_drop_geometry(x)

  na.omit(x[idINS])[1,1] |>
    gsub(x =_, pattern = "N.+", replacement = "") |>
    gsub(x =_, pattern = "r", replacement = "") |>
    as.numeric()
}

#' Transform a data.frame or a sf with an idINS to a stars object.
#' It doesn't take into account the geometry of the sf, but the idINS informations.
#' Assumes crs 3035
#'
#' @param x a data.frame or sf object with a idINS column.
#' @param crop a polygon of interest to crop the stars result.
#' @param idINS name of the idINS column. Default to idINS.
#' @param default_res if resolution isn't inscribed in idINS. Default to 200.
#'
#' @export
idINS2stars <- function(x, crop = NULL, idINS = "idINS", default_res = 200) {

  if (!is.null(crop)) x <- sf::st_crop(x, crop)
  if ("sf" %in% class(x)) x <- sf::st_drop_geometry(x)

  xy <- idINS2coord(x[idINS])

  xy_points <- xy |>
    sf::st_multipoint() |>
    sf::st_sfc(crs = st_crs(3035)) |>
    sf::st_cast(to = "POINT")

  res <- residINS(x)
  if (is.na(res)) res <- default_res

  bb <- sf::st_bbox(xy_points) + c(xmin = - res/2,
                                   ymin = - res/2,
                                   xmax = res/2,
                                   ymax = res/2)

  template <- stars::st_as_stars(bb, dx = res, dy = res)

  if (!"sf" %in% class(x)) x <- sf::st_as_sf(bind_cols(x, geometry = xy_points))

  stars::st_rasterize(x, template = template)

}
