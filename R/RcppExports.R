# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

get_wildfire_type <- function() {
    .Call(`_RCAFE_get_wildfire_type`)
}

TESTER <- function(n) {
    .Call(`_RCAFE_TESTER`, n)
}

#' @export
doFire <- function(tsf, regime) {
    .Call(`_RCAFE_doFire`, tsf, regime)
}

