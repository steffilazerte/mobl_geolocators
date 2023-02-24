
# Install packages and create folders if missing
install.packages("pak")

pak::pkg_install(c("dplyr", "tidyr", "lubridate", "ggplot2",
                   "sf", "rnaturalearth", "slisovski/SGAT", "slisovski/TwGeos"))

if(!dir.exists("MOBL_data/fixed")) dir.create("MOBL_data/fixed", recursive = TRUE)
if(!dir.exists("MOBL_data/final")) dir.create("MOBL_data/final", recursive = TRUE)


# Note that raw data files are expected to be found in /MOBL_data/RawData
