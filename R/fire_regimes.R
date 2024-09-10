#' Define a wildfire regime
#'
#' @param name (character) The regime name to use in the simulation outputs. It
#'   is best to avoid having any spaces in the name.
#'
#' @param fn_prob_tsf A function that takes a single numeric argument for time
#'   since fire and returns a probability value.
#'
#' @param fn_occur A function that takes a single integer argument for
#'   simulation time step and returns a logical value indicating whether a fire
#'   based on this regime should be simulated. Note that when a burn-in period
#'   is specified for a simulation (via the \code{n_burnin} argument to
#'   \code{\link{rcafe_simulate}}) the burn-in time step values are negative
#'   integers, while the main simulation time step values are positive integers.
#'   You can take advantage of this to specify that a fire regime only operates
#'   during the burn-in period (e.g. \code{fn_occur = function(x) x < 0}); or
#'   that it only operates during the main simulation period (e.g.
#'   \code{fn_occur = function(x) x > 0}).
#'
#' @param diagonal (logical; default \code{FALSE}) Whether fire can spread from
#'   a burning cell to its diagonal neighbours. By default, only orthogonal
#'   neighbours are considered.
#'
#' @param single_test (logical) If \code{TRUE} (the default), a cell that
#'   neighbours currently burning cells will only be tested for ignition once.
#'   If the cell does not ignite it will be marked as unavailable for future
#'   spread iterations of the current fire. If \code{FALSE}, a cell that does
#'   not ignite when tested will still be removed from the list of candidate
#'   cells for fire spread, but may be added again later if one of its other
#'   neighbours ignites. The effect of such multiple testing is to elevate the
#'   probability of cell ignition beyond the base rate calculated by
#'   \code{fn_prob_tsf}.
#'
#' @param max_ignition_attempts (single positive integer) Specifies the maximum
#'   number of cells to test for ignition when attempting to start a wildfire.
#'   The default of 1 means that if the initially selected cell fails to ignite,
#'   no fire occurs.
#'
#' @return A list object with named elements and class \code{CAFEregime}.
#'
#' @export
#
make_regime_wildfire <- function(name, fn_prob_tsf, fn_occur, diagonal = FALSE, single_test = TRUE, max_ignition_attempts = 1) {
  regime <- .make_regime_base(regime_type = .REGIME_TYPE_WILDFIRE, name, fn_prob_tsf, fn_occur)

  checkmate::assert_flag(diagonal)
  checkmate::assert_count(max_ignition_attempts, positive = TRUE)

  regime <- c(regime, list(diagonal = diagonal, single_test = single_test, max_ignition_attempts = max_ignition_attempts))

  class(regime) <- "CAFEregime"
  regime
}


#' Define a prescribed fire regime.
#'
#' A prescribed fire is simulated by randomly selecting a specified proportion
#' of cells from the landscape and testing each for ignition.
#'
#' @param name (character) The regime name to use in the simulation outputs. It
#'   is best to avoid having any spaces in the name.
#'
#' @param fn_prob_tsf A function that takes a single numeric argument for time
#'   since fire and returns a probability value.
#'
#' @param fn_occur A function that takes a single integer argument for
#'   simulation time step and returns a logical value indicating whether a fire
#'   based on this regime should be simulated. Note that when a burn-in period
#'   is specified for a simulation (via the \code{n_burnin} argument to
#'   \code{\link{rcafe_simulate}}) the burn-in time step values are negative
#'   integers, while the main simulation time step values are positive integers.
#'   You can take advantage of this to specify that a fire regime only operates
#'   during the burn-in period (e.g. \code{fn_occur = function(x) x < 0}); or
#'   that it only operates during the main simulation period (e.g.
#'   \code{fn_occur = function(x) x > 0}).
#'
#' @param prop_landscape A single numeric value between 0 and 1 specifying what
#'   proportion of landscape cells to attempt to burn. Cells are randomly
#'   selected.
#'
#' @param replace (logical) Controls whether cells that fail to ignite are
#'   replaced by other candidates. The default behaviour (\code{replace=FALSE})
#'   is to simply test each cell in the randomly selected set without
#'   replacement. If \code{replace=TRUE}, any cell that fails to ignite is
#'   replaced by selecting a new, untested landscape cell at random making it
#'   more likely that the specified proportion of the landscape will be burnt.
#'
#' @export
#
make_regime_prescribed <- function(name, fn_prob_tsf, fn_occur, prop_landscape, replace = FALSE) {
  regime = .make_regime_base(regime_type = .REGIME_TYPE_PRESCRIBED_FIRE, name, fn_prob_tsf, fn_occur)

  checkmate::assert_number(prop_landscape, lower = 0, upper = 1)

  regime <- c(regime, list(prop_landscape = prop_landscape, replace = replace))

  class(regime) <- "CAFEregime"
  regime
}


#' Private helper function to check and list variables common to all regime types
#'
#' @param regime_type Integer constant: either .
#'
#' @noRd
#
.make_regime_base <- function(regime_type, name, fn_prob_tsf, fn_occur) {
  checkmate::assert_choice(regime_type, choices = c(.REGIME_TYPE_WILDFIRE, .REGIME_TYPE_PRESCRIBED_FIRE))
  checkmate::assert_string(name, min.chars = 1)
  checkmate::assert_function(fn_prob_tsf, nargs = 1)
  checkmate::assert_function(fn_occur, nargs = 1)

  # Check that the probability function returns sensible values for a range of input TSF values
  x <- sapply(0:100, fn_prob_tsf)
  ok <- is.numeric(x) && all(!is.na(x)) && all(x >= 0.0 & x <= 1.0)
  if (!ok) {
    msg <- glue::glue("Probability function supplied for regime {name} returns values that are
                       not valid probabilities when tested with input values 0:100")
    stop(msg)
  }

  # Check that the occurrence function returns single logical values
  x <- sapply(1:100, fn_occur)
  ok <- is.logical(x) && all(!is.na(x))
  if (!ok) {
    msg <- glue::glue("Fire occurrence function supplied for regime {name} returns either non-logical
                       or missing values when tested with input values 1:100")
    stop(msg)
  }

  regime <- list(regime_type = regime_type,
                 name = name,
                 fn_prob_tsf = fn_prob_tsf,
                 fn_occur = fn_occur)

  regime
}


#' @export
print.CAFEregime <- function(x, ...) {
  if (x$regime_type == .REGIME_TYPE_WILDFIRE) {
    rtype <- "wildfire"
  } else {
    rtype <- "prescribed fire"
  }

  print( glue::glue("CAFE {rtype} regime: {x$name}") )
}

