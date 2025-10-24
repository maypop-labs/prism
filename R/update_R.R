# =============================================================================
# update_R
#
# Updates R to the latest version. After updating R, it is recommended that you
# also update to the appropriate version of RTools. See the README.md file for
# more details. If you would like to also update RStudio, go to the following:
# https://posit.co/download/rstudio-desktop/
#
# =============================================================================

if (!requireNamespace("installr", dependencies = TRUE)) {
  install.packages("installr")
}
library(installr)
updateR()
update.packages(ask = FALSE, checkBuilt = TRUE)
message("Done!")
