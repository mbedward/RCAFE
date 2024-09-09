#' Define a wildfire regime
#'
#' @export
#
make_regime_wildfire <- function(name, fn_prob_tsf, fn_occur, diagonal = FALSE, max_ignition_attempts = 1) {
  regime <- .make_regime_base(regime_type = .REGIME_TYPE_WILDFIRE, name, fn_prob_tsf, fn_occur)

  checkmate::assert_flag(diagonal)
  checkmate::assert_count(max_ignition_attempts, positive = TRUE)

  regime <- c(regime, list(diagonal = diagonal, max_ignition_attempts = max_ignition_attempts))

  class(regime) <- "CAFEregime"
  regime
}


#' Define a prescribed fire regime with random cell selection
#'
#' @export
#
make_regime_prescribed <- function(name, fn_prob_tsf, fn_occur, prop_landscape) {
  regime = .make_regime_base(regime_type = .REGIME_TYPE_PRESCRIBED_FIRE, name, fn_prob_tsf, fn_occur)

  checkmate::assert_number(prop_landscape, lower = 0, upper = 1)

  regime <- c(regime, list(prop_landscape = prop_landscape))

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

  glue::glue("CAFE {rtype} regime: {x$name}")
}

