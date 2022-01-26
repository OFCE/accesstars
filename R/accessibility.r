#' Calculation of accessibility on a stars object with an opportunity as attribute
#'
#' @param data Stars object of points with an information on opportunities
#' @param opportunities_colname Name of the opportunity attribute
#' @param departures_colname Name of an attribute to filter by, if the value is NA (Optional)
#' @param scope Maximal distance between origine and destination points, in km. If is NULL, scope is calculated with speed_max and duration_max
#' @param speed_max In km/h. By default, 60
#' @param duration_max in min. By default, 30 min which means scope is by default 30 km
#'
#' @import units
#' @import sf
#' @import stars
#' @import progress
#' @import data.table
#'
#' @export
#'
#' @section remarks : missing values on opportunity are excluded from any calculation, but not tiles with a zero value.
#' Transforming these values to NA before executing this function may be a good idea.
st_accessibility <- function(data,
                             opportunities_colname,
                             departures_colname = NULL,
                             scope = NULL, speed_max = 60, duration_max = 30,
                             r5r_core,
                             mode = "WALK", mode_egress = "WALK",
                             departure_datetime = Sys.time(),
                             time_window = 1L,
                             percentiles = 50L,
                             decay_function = "step",
                             cutoffs = 30L,
                             decay_value = 1,
                             max_walk_dist = Inf,
                             max_bike_dist = Inf,
                             max_trip_duration = 120L,
                             walk_speed = 3.6,
                             bike_speed = 12,
                             max_rides = 3,
                             max_lts = 2,
                             n_threads = Inf,
                             verbose = FALSE) {

  if (is.null(st_crs(data))) stop("Coordinates System Referential not defined.")

  resolution <- st_dimensions(data)$x$delta |> abs()
  if (is.na(resolution)) stop("Resolution is not defined in the stars dimensions")
  if (abs(st_dimensions(data)$y$delta) != resolution) stop("Tiles are not squares")

  # Creation of the idINS (crs 3035)
  if(st_crs(data) != st_crs(3035)) {
    idINS <- st_transform(data, crs = st_crs(3035)) |>
      st_coordinates() |>
      coord2idINS(resolution = resolution)
  } else {
    idINS <- st_coordinates(data) |>
      coord2idINS(resolution = resolution)
  }

  data[["idINS"]] <- matrix(idINS, nrow = dim(data)["x"], ncol = dim(data)["y"])
  message("Creation of idINS.\n")

  # Passage to WSG 84 mandatory to use r5r
  crs_ref <- st_crs(4326)
  if (st_crs(data) != crs_ref) data <- st_transform(data, crs = st_crs(4326))

  # Definition of scope
  if (is.null(scope)) {
    if (duration_max > 0 & speed_max > 0) {
      speed_max <- set_units(speed_max, km/h)
      duration_max <- set_units(duration_max, min)
      scope <- speed_max * duration_max
    } else stop("scope not defined")
  } else {
    scope <- set_units(scope, km)
  }

  # Summary on the area.
  # Assumes resolution in meters.
  side <- set_units(resolution, m)
  surf_area <- (st_dimensions(data)$x$to - st_dimensions(data)$x$from) * (st_dimensions(data)$y$to - st_dimensions(data)$y$from)
  surf_area <- surf_area * side^2
  message("Tile's dimension = ", side, " x ", side)
  message("Area's dimension = ", surf_area)
  message("Number of tiles on the area = ", round(surf_area/(side^2)))
  message("Scope = ", scope)
  message("Part of the area reachable (in percent) = ", round(100 * pi * scope^2 / surf_area), " %")

  # Selection of departures
  departs <- cbind(st_coordinates(data),
                   data[["idINS"]] |> as.vector()
  )
  names(departs) <- c("lon", "lat", "id")

  if (!is.null(departures_colname)) {
    exclure <- data[[departures_colname]] |> as.vector() |> is.na()
    departs <- departs[!exclure, ]
  }

  pb <- progress_bar$new(format = " [:bar] :percent in :elapsed eta: :eta",
                         total = nrow(departs),
                         clear = FALSE, width= 60)

  purrr::map_dfr(seq_along(idINS), function(i) {

    pb$tick()

    centre <- st_sfc(st_point(c(departs[i, 1], departs[i, 2])), crs = crs_ref)
    arrivees <- data[st_buffer(centre, scope)]

    arrivees <- cbind(st_coordinates(arrivees),
                      arrivees[["idINS"]] |> as.vector(),
                      arrivees[[opportunities_colname]] |> as.vector()
    )
    names(arrivees) <- c("lon", "lat", "id", opportunities_colname)
    arrivees <- arrivees[!is.na(opportunities_colname),]

    if (nrow(arrivees) == 0) return(NULL)

    r5r::accessibility(r5r_core = r5r_core,
                       origins = departs[i, ],
                       destinations = arrivees,
                       opportunities_colname = opportunities_colname,
                       mode = mode, mode_egress = mode_egress,
                       departure_datetime = departure_datetime,
                       time_window = time_window, percentiles = percentiles,
                       decay_function = decay_function, cutoffs = cutoffs, decay_value = decay_value,
                       max_walk_dist = max_walk_dist, max_bike_dist = max_bike_dist, max_trip_duration = max_trip_duration,
                       walk_speed = walk_speed, bike_speed = bike_speed, max_rides = max_rides,
                       max_lts = max_lts, n_threads = n_threads, verbose = verbose, progress = FALSE)

  })
}

#' Transform result of r5r::accessibiliy or st_accessibility from data.table to stars
#'
#' @param dt data.table of accessibility
#' @param cutoffs cutoffs of isochrones
#' @param percentiles percentiles of time travel
#' @param nom name of the opportunity which accessibility is evaluated
#'
#' @import data.table
acc_dt2stars <- function(dt, cutoffs, percentiles, nom = "accessibility") {

  grid_noms <- expand.grid(cutoffs, percentiles)
  names(grid_noms) <- c("cutoffs", "percentiles")
  if (length(percentiles) == 1) {
    if(length(cutoffs) == 1) {les_noms <- nom} else {les_noms <- paste0(nom, "_", cutoffs, "mn")}
  } else {
    if(length(cutoffs) == 1) {
      les_noms <- paste0(nom, "_", percentiles, "p")
    } else {
      les_noms <- paste0(nom, "_", grid_noms[, 1], "mn_", grid_noms[, 2], "p")
    }
  }

  res_stars <- purrr::imap(les_noms, ~ {
    res <- dt[percentile == grid_noms$percentiles[.y] & cutoff == grid_noms$cutoffs[.y], ][,
                                                                                           ':='(percentile = NULL, cutoff = NULL)] |>
      df2stars(idINS = "from_id")
    names(res) <- .x
    res
  })

  purrr::reduce(res_stars, c)

}


