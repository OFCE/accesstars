#' Calculation of accessibility on a stars object with an opportunity as attribute
#'
#' @param data Stars object of points with an information on opportunities
#' @param opportunities_colname Name of the opportunity attribute
#' @param departures_colname Name of an attribute to filter by, if the value is NA (Optional)
#' @param r5r_core a rJava object to connect with R5 routing engine
#' @param stars_result Defaults to TRUE for stars object. Otherwise a data.table.
#' @param mode string. Transport modes allowed for the trips. Defaults to
#'             "WALK". See details for other options.
#' @param mode_egress string. Transport mode used after egress from public
#'                    transport. It can be either 'WALK', 'BICYCLE', or 'CAR'.
#'                    Defaults to "WALK".
#' @param departure_datetime POSIXct object. If working with public transport
#'                           networks, please check \code{calendar.txt} within
#'                           the GTFS file for valid dates. See details for
#'                           further information on how datetimes are parsed.
#' @param time_window numeric. Time window in minutes for which r5r will
#'                    calculate travel times departing each minute. When using
#'                    frequency-based GTFS files, 5 Monte Carlo simulations will
#'                    be run for each minute in the time window. See details for
#'                    further information.
#' @param percentiles numeric vector. Defaults to '50', returning the accessibility
#'                    value for the median travel time computed for a given
#'                    time_window. If a numeric vector is passed, for example
#'                    c(25, 50, 75), the function will return accessibility
#'                    estimates for each percentile, by travel time cutoff. Only
#'                    the first 5 cut points of the percentiles are considered.
#'                    For more details, see R5 documentation at
#'                    'https://docs.conveyal.com/analysis/methodology#accounting-for-variability'
#' @param decay_function string. Choice of one of the following decay functions:
#'                       'step', 'exponential', 'fixed_exponential', 'linear',
#'                       and 'logistic'. Defaults to 'step', which yields
#'                       cumulative opportunities accessibility metrics.
#'                       More info in `details`.
#' @param cutoffs numeric. Cutoff times in minutes for calculating cumulative
#'                opportunities accessibility when using the 'step decay function'.
#'                This parameter has different effects for each of the other decay
#'                functions: it indicates the 'median' (or inflection point) of
#'                the decay curves in the 'logistic' and 'linear' functions, and
#'                the 'half-life' in the 'exponential' function. It has no effect
#'                when using the 'fixed exponential' function.
#' @param decay_value numeric. Extra parameter to be passed to the selected
#'                `decay_function`.
#' @param max_walk_dist numeric. Maximum walking distance (in meters) to access
#'                      and egress the transit network, or to make transfers
#'                      within the network. Defaults to no restrictions as long
#'                      as `max_trip_duration` is respected. The max distance is
#'                      considered separately for each leg (e.g. if you set
#'                      `max_walk_dist` to 1000, you could potentially walk up
#'                      to 1 km to reach transit, and up to _another_  1 km to
#'                      reach the destination after leaving transit). Obs: if you
#'                      want to set the maximum walking distance considering
#'                      walking-only trips you have to set the `max_trip_duration`
#'                      accordingly (e.g. to set a distance of 1 km assuming a
#'                      walking speed of 3.6 km/h you have to set `max_trip_duration = 1 / 3.6 * 60`).
#' @param max_bike_dist numeric. Maximum cycling distance (in meters) to access
#'                      and egress the transit network. Defaults to no
#'                      restrictions as long as `max_trip_duration` is respected.
#'                      The max distance is considered separately for each leg
#'                      (e.g. if you set `max_bike_dist` to 1000, you could
#'                      potentially cycle up to 1 km to reach transit, and up
#'                      to _another_ 1 km to reach the destination after leaving
#'                      transit). Obs: if you want to set the maximum cycling
#'                      distance considering cycling-only trips you have to set
#'                      the `max_trip_duration` accordingly (e.g. to set a
#'                      distance of 5 km assuming a cycling speed of 12 km/h you
#'                      have to set `max_trip_duration = 5 / 12 * 60`).
#' @param max_trip_duration numeric. Maximum trip duration in minutes. Defaults
#'                          to 120 minutes (2 hours).
#' @param walk_speed numeric. Average walk speed in km/h. Defaults to 3.6 km/h.
#' @param bike_speed numeric. Average cycling speed in km/h. Defaults to 12 km/h.
#' @param max_rides numeric. The max number of public transport rides allowed in
#'                  the same trip. Defaults to 3.
#' @param max_lts  numeric (between 1 and 4). The maximum level of traffic stress
#'                 that cyclists will tolerate. A value of 1 means cyclists will
#'                 only travel through the quietest streets, while a value of 4
#'                 indicates cyclists can travel through any road. Defaults to 2.
#'                 See details for more information.
#' @param n_threads numeric. The number of threads to use in parallel computing.
#'                  Defaults to use all available threads (Inf).
#' @param verbose logical. `TRUE` to show detailed output messages (the default).
#' @param progress logical. `TRUE` to show a progress counter. Only works when
#'                `verbose` is set to `FALSE`, so the progress counter does not
#'                interfere with R5's output messages. Setting `progress` to `TRUE`
#'                may impose a small penalty for computation efficiency, because
#'                the progress counter must be synchronized among all active
#'                threads.
#' @return A stars objects
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
                             stars_result = TRUE,
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
                             verbose = FALSE,
                             progress = TRUE) {

  if (is.null(st_crs(data))) stop("Coordinates System Referential not defined.")

  resolution <- st_dimensions(data)$x$delta |> abs()
  if (is.na(resolution)) stop("Resolution is not defined in the stars dimensions")
  if (abs(st_dimensions(data)$y$delta) != resolution) stop("Tiles are not squares")

  tictoc::tic()

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

  # Summary on the area.
  # Assumes resolution in meters.
  side <- set_units(resolution, m)
  surf_area <- (st_dimensions(data)$x$to - st_dimensions(data)$x$from) * (st_dimensions(data)$y$to - st_dimensions(data)$y$from)
  surf_area <- surf_area * side^2
  surf_area <- set_units(surf_area, km^2)
  message("Tile's dimension = ", side, "m x ", side, "m")
  message("Area's dimension = ", surf_area, " km^2")
  message("Number of tiles on the area ~ ", ceiling(surf_area/(side^2)))

  # Selection of origins
  departs <- cbind(st_coordinates(data),
                   data[["idINS"]] |> as.vector())

  names(departs) <- c("lon", "lat", "id")

  if (!is.null(departures_colname)) {
    exclure <- data[[departures_colname]] |> as.vector() |> is.na()
    departs <- departs[!exclure, ]
  }

  # Selection of destinations
  arrivees <- cbind(st_coordinates(data),
                    data[["idINS"]] |> as.vector(),
                    data[[opportunities_colname]] |> as.vector())

  names(arrivees) <- c("lon", "lat", "id", opportunities_colname)

  arrivees <- arrivees[!is.na(opportunities_colname),]

  resultat <- r5r::accessibility(r5r_core = r5r_core,
                                 origins = departs,
                                 destinations = arrivees,
                                 opportunities_colname = opportunities_colname,
                                 mode = mode, mode_egress = mode_egress,
                                 departure_datetime = departure_datetime,
                                 time_window = time_window, percentiles = percentiles,
                                 decay_function = decay_function, cutoffs = cutoffs, decay_value = decay_value,
                                 max_walk_dist = max_walk_dist, max_bike_dist = max_bike_dist, max_trip_duration = max_trip_duration,
                                 walk_speed = walk_speed, bike_speed = bike_speed, max_rides = max_rides,
                                 max_lts = max_lts, n_threads = n_threads, verbose = verbose, progress = progress)

  tictoc::toc()

  if (stars_result){
    acc_dt2stars(resultat, cutoffs = cutoffs, percentiles = percentiles, nom = opportunities_colname)
  } else {
    resultat
  }
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
              dt2stars(idINS = "from_id")
    names(res) <- .x
    res
  })

  purrr::reduce(res_stars, c)

}


