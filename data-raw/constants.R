# Create constants that identify fire regime types.
# These are not exported by the package.

.REGIME_TYPE_WILDFIRE <- 1
.REGIME_TYPE_PRESCRIBED_FIRE <- 2

usethis::use_data(.REGIME_TYPE_WILDFIRE, .REGIME_TYPE_PRESCRIBED_FIRE, internal = TRUE, overwrite = TRUE)
