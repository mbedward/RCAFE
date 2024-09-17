#' Run an RCAFE simulation
#'
#' @param DB_path Path to a folder for the DuckDB database file to which
#'   simulation results will be recorded.
#'
#' @param tsf_init A matrix of of initial time since fire values for the
#'   landscape. The dimensions of this matrix will be used to define those of
#'   the simulation landscape. Values should be non-negative integers.
#'
#' @param regimes A list of one or definitions for wildfire and/or prescribed
#'   fire regime(s).
#'
#' @param n_rep (positive integer; default 1) Number of replicate simulations.
#'
#' @param n_iter (positive integer; default 100) Number of iterations (time
#'   steps) per simulation.
#'
#' @param n_burnin (non-negative integer; default 0) Number of initial time
#'   steps to run the simulation for prior to the \code{n_iter} times. No data
#'   are recorded for these initial time steps. This can be useful to evolve an
#'   initial, spatially-structured pattern of time since fire (see examples).
#'
#' @param tsf_breaks An integer vector specifying the breakpoints over which to
#'   record the proportions of the simulation landscape within time since fire
#'   intervals at each time step (year). For example, the default value of
#'   \code{tsf_breaks = c(5, 10, 20)} specifies to record landscape proportions
#'   in the intervals 0-5 years;  6-10 years; 10-20 years; and greater than 20
#'   years.
#'
#' @param tsf_map_times \strong{Not implemented yet}. An integer vector
#'   specifying the times at which to record a map of time since fire. Maps will
#'   only be recorded after any burn-in time steps (argument \code{n_burnin})
#'   have been completed. For example, the arguments
#'   \code{n_iter=500, n_burnin=100, tsf_map_times = c(1, seq(100, 500,
#'   50))} specify to run 100 burn-in simulations to evolve a time since fire
#'   pattern in the landscape, and then record a map in the first year and every
#'   hundred years of the simulation period. If this vector is missing,
#'   \code{NULL} or empty, no time since fire maps will be saved.
#'
#' @param map_fire_types (character; default \code{'none'}) Whether to save maps
#'   of fires in the output database. Can be one of the following: \code{'all'}
#'   (save both prescribed and wild fires); \code{'wild'} (save only wildfires);
#'   \code{'prescribed'} (save only prescribed fires); \code{'none'} (do not
#'   save fire maps). Can be abbreviated (e.g. \code{map_fire_types = 'w'}).
#'
#' @param display_progress (logical) Whether to display a progress bar.
#'
#' @return (invisibly) The path to the output database as provided via the
#'   \code{DB_path} argument.
#'
#' @examples
#' \dontrun{
#' # Function to relate cell burning probability to time since fire
#' fn_prob <- function(tsf) {
#'   1 / (1 + exp(1.65 - 0.16*tsf))
#' }
#'
#' # Visualize the function for a range of time since fire values
#' curve(fn_prob(x), from=0, to=40,
#'       ylim = c(0, 1), xlab = "TSF", ylab = "Burn probability")
#'
#' # Define a 100 x 100 landscape with uniform initial time since fire
#' # of 20 years
#' tsf0 <- matrix(20, nrow=100, ncol=100)
#'
#' # Path for the output DuckDB database file
#' dbpath <- "test_sim.duckdb"
#'
#' # Regime definition for a wildfire that (possibly) burns every 5th year of
#' # the simulation
#' regime_wild <- make_regime_wildfire("wildfire",
#'                                     fn_prob,
#'                                     fn_occur = function(x) {x %% 5 == 0})
#'
#' # Run the simulation for a single replicate of 500 years after an initial 100 year
#' # burn-in period
#' res <- rcafe_simulate(dbpath,
#'                       tsf_init = tsf0,
#'                       regimes = list(regime_wild),
#'                       n_rep = 1,
#'                       n_iter = 500,
#'                       n_burnin = 100)
#' }
#'
#' @export
#
rcafe_simulate <- function(DB_path,
                           tsf_init,
                           regimes,
                           n_rep = 1,
                           n_iter = 100,
                           n_burnin = 0,
                           tsf_breaks = c(5, 10, 20),
                           tsf_map_times = NULL,
                           map_fire_types = c('none', 'all', 'wild', 'prescribed'),
                           display_progress = TRUE) {

  checkmate::assert_path_for_output(DB_path)



  checkmate::assert_matrix(tsf_init, mode = "integerish", any.missing = FALSE)

  LandscapeDim <- dim(tsf_init)
  LandscapeNCells <- prod(LandscapeDim)

  checkmate::assert_list(regimes, any.missing = FALSE, min.len = 1)

  # TODO - check fire regime(s) for validity ?

  checkmate::assert_count(n_rep, positive = TRUE)
  checkmate::assert_count(n_iter, positive = TRUE)
  checkmate::assert_count(n_burnin)

  # breaks for monitoring TSF intervals
  checkmate::assert_integerish(tsf_breaks, lower = 0, any.missing = FALSE, min.len = 1, unique = TRUE)
  if (any(is.infinite(tsf_breaks))) stop("All values in tsf_breaks should be finite integers")

  tsf_breaks <- sort(tsf_breaks)

  map_fire_types <- match.arg(map_fire_types)

  map_prescribed <- map_fire_types %in% c('all', 'prescribed')
  map_wild <- map_fire_types %in% c('all', 'wild')

  # Whether to display a progress bar
  checkmate::assert_flag(display_progress)


  # Database connection for outputs (NOTE: PRESENTLY MUST BE A FILE)
  DB <- DBI::dbConnect(drv = duckdb::duckdb(dbdir = DB_path))

  ##### Create database table for fire regimes
  cmd <- glue::glue("CREATE TABLE regimes (
                       id INTEGER,
                       regime_type INTEGER,
                       name TEXT,
                       PRIMARY KEY(id)
                     );")

  DBI::dbExecute(DB, cmd)

  dat <- lapply(seq_along(regimes), function(i) {
    r <- regimes[[i]]
    data.frame(id = i,
               regime_type = r$regime_type,
               name = r$name)
  })
  dat <- do.call(rbind, dat)

  DBI::dbWriteTable(DB, name = "regimes", value = dat, append = TRUE)

  ##### Create database table for fire statistics
  cmd <- glue::glue("create table fires (
                       rep INTEGER,
                       time INTEGER,
                       regime_id INTEGER REFERENCES regimes(id),
                       size INTEGER,
                       PRIMARY KEY(rep, time, regime_id)
                     );")

  DBI::dbExecute(DB, cmd)

  ##### Create database table for the paths to the optional fire map raster files
  cmd <- glue::glue("CREATE TABLE firemaps (
                       rep INTEGER,
                       regime_id INTEGER REFERENCES regimes(id),
                       raster_path TEXT,
                       PRIMARY KEY(rep, regime_id)
                     );")

  DBI::dbExecute(DB, cmd)

  ##### Create database table for TSF interval monitoring
  s <- sprintf("tsf_leq_%g DOUBLE", tsf_breaks)
  s <- c(s, sprintf("tsf_gt_%g DOUBLE", tail(tsf_breaks, 1)))
  TSF_colnames <- paste(s, collapse = ", ")

  tsf_breaks_schema <- glue::glue("CREATE TABLE tsf_proportion (
                                     rep INTEGER,
                                     time INTEGER,
                                     {TSF_colnames},
                                     PRIMARY KEY(rep, time)
                                   );")

  DBI::dbExecute(DB, tsf_breaks_schema)

  # Create a string for a prepared statement to insert values into the table
  TSF_insert_stmt <- paste0("insert into tsf_proportion values (",
                           paste(rep('?', length(tsf_breaks) + 3), collapse = ", "),
                           ")")

  # Function to calculate TSF proportions
  fn_tsf_proportions <- function(tsf) {
    ii <- base::findInterval(c(tsf), c(tsf_breaks + 1, Inf))
    tabulate(ii + 1, nbins = length(tsf_breaks) + 1) / LandscapeNCells
  }

  # Path prefix for file name for fire rasters based on the output DB path
  FireRasterPrefix <- fs::path_ext_remove(DB_path)


  ##### Simulation starts here #####

  all_times <- c(-rev(seq_len(n_burnin)), seq_len(n_iter))

  if (display_progress) {
    if (n_rep == 1) bar_fmt <- "Time :current :percent [:bar]"
    else bar_fmt <- "Replicate :rep time :current :percent [:bar]"
  }

  for (irep in seq_len(n_rep)) {
    if (display_progress) {
      pbar <- progress::progress_bar$new(format = bar_fmt, total = n_iter, clear = TRUE)
      pbar$tick(0)
    }

    tsf = tsf_init

    res_replicate <- lapply(all_times, function(itime) {
      #
      # Simulate each regime
      #
      res_regimes <- lapply(seq_along(regimes), function(iregime) {
        regime <- regimes[[iregime]]
        res_cur_regime <- NULL

        if (regime$fn_occur(itime)) {
          fire_res <- doFire(tsf = tsf, regime)

          # Fire map (0 = unburnt, 1 = burnt)
          fire_footprint <- fire_res[["landscape"]]

          # Update the TSF matrix for burnt and unburnt cells
          tsf <<- tsf + 1

          # Set burnt cells to zero TSF
          tsf[ which(fire_footprint > 0) ] <<- 0

          # If there was a fire, and we are in the main simulation (after any burn-in phase)
          # record the fire summary and create a fire raster if requested
          #
          ncells_burnt <- fire_res[["ncells"]]

          if (ncells_burnt > 0 && itime > 0) {
            stmt <- DBI::dbSendStatement(DB, "insert into fires values (?, ?, ?, ?);")
            DBI::dbBind(stmt, list(irep, itime, iregime, ncells_burnt))
            DBI::dbClearResult(stmt)

            do_map <- (
              (regime$regime_type == .REGIME_TYPE_WILDFIRE && map_wild) ||
              (regime$regime_type == .REGIME_TYPE_PRESCRIBED_FIRE && map_prescribed)
            )

            if (do_map) {
              r <- terra::rast(fire_footprint)
              res_cur_regime <- list(rep = irep, time = itime, regime_index = iregime, layer = r)
            }
          }
        }

        # Return the regime map (or NULL if there isn't one)
        res_cur_regime
      })

      # Record TSF landscape proportions after all fires have been simulated for this year
      if (itime > 0) {
        if (display_progress) {
          tokens <- list()
          if (n_rep > 1) tokens <- list(rep = irep)
          pbar$tick(tokens = tokens)
        }

        p <- fn_tsf_proportions(tsf)
        vals <- c(list(irep, itime), p)

        stmt <- DBI::dbSendStatement(DB, TSF_insert_stmt)
        DBI::dbBind(stmt, vals)
        DBI::dbClearResult(stmt)
      }

      # Return list of regime result objects with fire maps (or NULLs if no maps were created)
      res_regimes
    })

    # Flatten the list of fire layers from regimes-within-time and write
    # rasters for each regime to a GeoTIFF file
    #
    res_all <- .flatten_list(res_replicate)

    for (iregime in seq_along(regimes)) {
      regime_name <- regimes[[iregime]]$name

      rlayers <- lapply(res_all, function(res) {
        if (!is.null(res) && res$regime_index == iregime) res$layer
      })
      rlayers <- base::Filter(.is_not_null, rlayers)

      times <- sapply(res_all, function(res) {
        if (!is.null(res) && res$regime_index == iregime) res$time
      })
      times <- unlist(times) # convert to vector and implicitly remove any NULLs

      if (length(rlayers) > 0) {
        regime_raster <- terra::rast(rlayers)
        names(regime_raster) <- sprintf("time%06d", times)

        out_path <- sprintf("%s_%s_%06d.tif",
                            FireRasterPrefix,
                            regime_name,
                            irep)

        terra::writeRaster(regime_raster,
                           filename = out_path,
                           overwrite = TRUE,
                           datatype = "INT1U",
                           gdal = c("COMPRESS=ZSTD", "ZSTD_LEVEL=1", "PREDICTOR=2"))

        # Store raster path in output database
        stmt <- DBI::dbSendStatement(DB, "insert into firemaps values (?, ?, ?);")
        DBI::dbBind(stmt, list(irep, iregime, out_path))
        DBI::dbClearResult(stmt)
      }
    }
  }

  # Close the database connection and return the path to the database file
  DBI::dbDisconnect(DB)

  invisible(DB_path)
}


# Private helper function to flatten a nested list
#
# See https://stackoverflow.com/a/77492151/40246
#
.flatten_list <- function(x, flatf = c) x |> Reduce (\(a,b) flatf(a,b) , x = _)


# Private helper function to test if a value is not NULL
.is_not_null <- function(x) !is.null(x)

