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
#' @param n_sim (positive integer) Maximum number of time steps for the
#'   simulation.
#'
#' @param n_burnin (non-negative integer; default 0) Number of initial time steps to
#'   run the simulation for prior to the \code{n_sim} times. No data are
#'   recorded for these initial time steps. This can be useful to evolve an
#'   initial, spatially-structured pattern of time since fire (see examples).
#'
#' @param tsf_breaks An integer vector specifying the breakpoints over which to
#'   record the proportions of the simulation landscape within time since fire
#'   intervals at each time step (year). For example, the default value of
#'   \code{tsf_breaks = c(5, 10, 20)} specifies to record landscape proportions
#'   in the intervals 0-5 years;  6-10 years; 10-20 years; and greater than 20
#'   years.
#'
#' @param tsf_map_times An integer vector specifying the times at which to record
#'   a map of time since fire. Maps will only be recorded after any burn-in time
#'   steps (argument \code{n_burnin}) have been completed. For example, the
#'   arguments \code{n_sim=500, n_burnin=100, tsf_map_times = c(1, seq(100, 500,
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
#' # Run the simulation for 500 years after an initial 100 year
#' # burn-in period
#' res <- rcafe_simulate(tsf0, n_sim=500, n_burnin=100,
#'                       fn_prob_tsf = fn_prob)
#' }
#'
#' @export
#
rcafe_simulate <- function(DB_path,
                           tsf_init,
                           regimes,
                           n_sim,
                           n_burnin = 0,
                           tsf_breaks = c(5, 10, 20),
                           tsf_map_times = NULL,
                           map_fire_types = 'none',
                           display_progress = TRUE) {

  checkmate::assert_path_for_output(DB_path)

  checkmate::assert_matrix(tsf_init, mode = "integerish", any.missing = FALSE)

  LandscapeDim <- dim(tsf_init)
  LandscapeNCells <- prod(LandscapeDim)

  checkmate::assert_list(regimes, any.missing = FALSE, min.len = 1)

  # TODO - check fire regime(s) for validity ?

  checkmate::assert_count(n_sim, positive = TRUE)
  checkmate::assert_count(n_burnin)

  # breaks for monitoring TSF intervals
  checkmate::assert_integerish(tsf_breaks, lower = 0, any.missing = FALSE, min.len = 1, unique = TRUE)
  if (any(is.infinite(tsf_breaks))) stop("All values in tsf_breaks should be finite integers")

  tsf_breaks <- sort(tsf_breaks)

  # Whether to display a progress bar
  checkmate::assert_flag(display_progress)


  # Database connection for outputs (NOTE: PRESENTLY MUST BE A FILE)
  DB <- DBI::dbConnect(drv = duckdb::duckdb(dbdir = DB_path))

  # Create table for fire regimes
  cmd <- glue::glue("CREATE TABLE regimes (
                       id INTEGER,
                       regime_type INTEGER,
                       name VARCHAR(32),
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


  # Create database table for fire statistics
  cmd <- glue::glue("create table fires (
                       time INTEGER,
                       regime_id INTEGER REFERENCES regimes(id),
                       size INTEGER,
                       PRIMARY KEY(time, regime_id)
                     );")

  DBI::dbExecute(DB, cmd)

  # Create database table for TSF interval monitoring
  s <- sprintf("tsf_leq_%g DOUBLE", tsf_breaks)
  s <- c(s, sprintf("tsf_gt_%g DOUBLE", tail(tsf_breaks, 1)))
  TSF_colnames <- paste(s, collapse = ", ")

  tsf_breaks_schema <- glue::glue("CREATE TABLE tsf_proportion (time INTEGER, {TSF_colnames});")

  DBI::dbExecute(DB, tsf_breaks_schema)

  # Create a string for a prepared statement to insert values into the table
  TSF_insert_stmt <- paste0("insert into tsf_proportion values (",
                           paste(rep('?', length(tsf_breaks) + 2), collapse = ", "),
                           ")")

  # Function to calculate TSF proportions
  fn_tsf_proportions <- function(tsf) {
    ii <- base::findInterval(c(tsf), c(tsf_breaks + 1, Inf))
    tabulate(ii + 1, nbins = length(tsf_breaks) + 1) / LandscapeNCells
  }

  ##### Simulation starts here #####
  tsf = tsf_init

  all_times <- c(-rev(seq_len(n_burnin)), seq_len(n_sim))

  progress_time_step <- (n_burnin + n_sim) / 10
  progress_times <- unique(round(seq(min(all_times), max(all_times), by = progress_time_step)))[-1]

  iSaveTSF <- 0

  for (itime in all_times) {
    if (display_progress && (itime %in% progress_times)) {
      cat("=", sep = "")
    }

    for (iregime in seq_along(regimes)) {
      regime <- regimes[[iregime]]

      if (regime$fn_occur(itime)) {
        fire_res <- doFire(tsf = tsf, regime)

        ncells_burnt <- fire_res[["ncells"]]
        if (ncells_burnt > 0 && itime > 0) {
          stmt <- DBI::dbSendStatement(DB, "insert into fires values (?, ?, ?);")
          DBI::dbBind(stmt, list(itime, iregime, ncells_burnt))
          DBI::dbClearResult(stmt)
        }

        landscape <- fire_res[["landscape"]]

        # Update the TSF matrix for burnt and unburnt cells
        tsf <- tsf + 1

        # Set burnt cells to zero TSF
        tsf[ which(landscape > 0) ] <- 0
      }
    }

    # Record TSF landscape proportions
    if (itime > 0) {
      p <- fn_tsf_proportions(tsf)
      vals <- c(list(itime), p)

      stmt <- DBI::dbSendStatement(DB, TSF_insert_stmt)
      DBI::dbBind(stmt, vals)
      DBI::dbClearResult(stmt)
    }
  }

  if (display_progress) cat("\n")

  # Close the database connection and return
  DBI::dbDisconnect(DB)
}

