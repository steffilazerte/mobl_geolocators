library(SGAT)
library(TwGeos)
library(lubridate) # Better dates
library(ggplot2)   # Slower but easier plots


# Set details ---------------------------------------------
dir <- file.path("RMA4406 Deployed Data/")
out_dir <- "MOBL_data/fixed/"

# Kamloops, release location
lon_calib <- -120.33
lat_calib <- 50.67

# Offset
test_date <- as_datetime("2020-07-15", tz = "UTC")
test_date - force_tz(test_date, "America/Vancouver")

offset <- 7 + 12  # Difference between local and UTC


# ID000539 --------------------------------
ID  <- "ID000539"
file <- "ID000539_2022Oct05T142349.csv"

# Cut off problematic bits
date_range <- as_date(c("2020-07-15", "2021-05-29"))


## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, file))
names(raw) <- c("Date", "Light")
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]

# Filter to date range
raw <- filter(raw, Date >= date_range[1], Date <= date_range[2])

# Log-transform because these record the full light spectrum
# https://geolocationmanual.vogelwarte.ch/loadingData.html
raw$Light_orig <- raw$Light
raw$Light <- log10(raw$Light_orig)
head(raw)

#raw$Light[raw$Light > 100] <- 100

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Check Date range ----------------------

# Plot just a fraction of the points
ggplot(data = dplyr::slice_sample(raw, prop = 0.05), aes(x = Date, y = Light)) +
  geom_point()

# Look at the end of the range
ggplot(data = raw, aes(x = Date, y = Light)) +
  geom_point() +
  coord_cartesian(xlim = c(as_datetime("2021-05-28"), as_datetime("2021-06-01")))

# Use May 29th as cutoff (or best to use date before recovery if known)


## Checking times -----------------------------------------
# Plot
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# Adjustment
adjust <- 90

raw_summer <- filter(raw, Date < "2020-09-01")
raw_summer <- mutate(raw_summer, Date = Date + minutes(adjust))

lightImage(tagdata = raw_summer,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

tsimageDeploymentLines(raw_summer$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# Check sunrise/sunset
suntimes <- c("2020-07-15 05:07:00", "2020-07-16 05:08:00", "2020-07-17 05:09:00",
              "2020-07-15 21:07:00", "2020-07-16 21:06:00", "2020-07-17 21:05:00")
suntimes <- as_datetime(suntimes, tz = "America/Vancouver")

ggplot(data = raw[raw$Date >= as_date("2020-07-15") & raw$Date <= as_date("2020-07-19"), ],
       aes(x = Date, y = Light)) +
  geom_point() +
  geom_vline(xintercept = suntimes) + #Original
  geom_vline(xintercept = suntimes - minutes(adjust), colour = "orange") # Adjusted



# Using adjustment
raw_fix <- mutate(raw, Date = Date + minutes(adjust))
write.csv(raw_fix, file.path(out_dir, paste0(ID, "_fixed.csv")), row.names = FALSE)






# ID000547 --------------------------------
ID  <- "ID000547"
file <- "ID000547_2022Oct05T143422.csv"

# Cut off problematic bits
date_range <- as_date(c("2020-07-01", "2021-08-30"))

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, file))
names(raw) <- c("Date", "Light")
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]

# Filter to date range
raw <- filter(raw, Date >= date_range[1], Date <= date_range[2])

# Log-transform because these record the full light spectrum
# https://geolocationmanual.vogelwarte.ch/loadingData.html
raw$Light_orig <- raw$Light
raw$Light <- log10(raw$Light_orig)
head(raw)

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Check Date range ----------------------

# Plot just a fraction of the points
ggplot(data = dplyr::slice_sample(raw, prop = 0.05), aes(x = Date, y = Light)) +
  geom_point()

# Look at the end of the range
ggplot(data = raw, aes(x = Date, y = Light)) +
  geom_point() +
  coord_cartesian(xlim = c(as_datetime("2021-05-28"), as_datetime("2021-09-01")))


## Checking times -----------------------------------------
# Plot
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# Adjustment
adjust <- 140

raw_summer <- filter(raw, Date < "2020-09-01")
raw_summer <- mutate(raw_summer, Date = Date + minutes(adjust))

lightImage(tagdata = raw_summer,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

tsimageDeploymentLines(raw_summer$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# Check sunrise/sunset
suntimes <- c("2020-07-15 05:07:00", "2020-07-16 05:08:00", "2020-07-17 05:09:00",
              "2020-07-15 21:07:00", "2020-07-16 21:06:00", "2020-07-17 21:05:00")
suntimes <- as_datetime(suntimes, tz = "America/Vancouver")

ggplot(data = raw[raw$Date >= as_date("2020-07-15") & raw$Date <= as_date("2020-07-19"), ],
       aes(x = Date, y = Light)) +
  geom_point() +
  geom_vline(xintercept = suntimes) + #Original
  geom_vline(xintercept = suntimes - minutes(adjust), colour = "orange") # Adjusted



# Using adjustment
raw_fix <- mutate(raw, Date = Date + minutes(adjust))
write.csv(raw_fix, file.path(out_dir, paste0(ID, "_fixed.csv")), row.names = FALSE)

# ID000552 --------------------------------
ID  <- "ID000552"
file <- "ID000552_2022Oct05T152519.csv"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, file))
names(raw) <- c("Date", "Light")
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]

# Log-transform because these record the full light spectrum
# https://geolocationmanual.vogelwarte.ch/loadingData.html
raw$Light_orig <- raw$Light
raw$Light <- log10(raw$Light_orig)
head(raw)

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Check Date range ----------------------
date_range <- as_date(c("2020-06-01", "2021-05-25"))

# Plot just a fraction of the points
ggplot(data = dplyr::slice_sample(raw, prop = 0.05), aes(x = Date, y = Light)) +
  geom_point() +
  geom_vline(xintercept = as_datetime(date_range), colour = "red")

# Filter to date range
raw <- filter(raw, Date >= date_range[1], Date <= date_range[2])


## Checking times -----------------------------------------
# Plot
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# GOOD!

write.csv(raw, file.path(out_dir, paste0(ID, "_fixed.csv")), row.names = FALSE)


# ID000553 --------------------------------
ID  <- "ID000553"
file <- "ID000553_2022Oct05T144832.csv"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, file))
names(raw) <- c("Date", "Light")
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]

# Log-transform because these record the full light spectrum
# https://geolocationmanual.vogelwarte.ch/loadingData.html
raw$Light_orig <- raw$Light
raw$Light <- log10(raw$Light_orig)
head(raw)

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Check Date range ----------------------
date_range <- as_date(c("2020-06-01", "2021-07-01"))

# Plot just a fraction of the points
ggplot(data = dplyr::slice_sample(raw, prop = 0.05), aes(x = Date, y = Light)) +
  geom_point() +
  geom_vline(xintercept = as_datetime(date_range), colour = "red")

# Filter to date range
raw <- filter(raw, Date >= date_range[1], Date <= date_range[2])


## Checking times -----------------------------------------
# Plot
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# Adjustment
adjust <- 40

raw_summer <- filter(raw, Date > "2020-07-15", Date < "2020-09-01")
raw_summer <- mutate(raw_summer, Date = Date + minutes(adjust))

lightImage(tagdata = raw_summer,
           offset = offset,
           zlim = light_lim) # Use actual ranges in data

tsimageDeploymentLines(raw_summer$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

# Check sunrise/sunset
suntimes <- c("2020-07-15 05:07:00", "2020-07-16 05:08:00", "2020-07-17 05:09:00",
              "2020-07-15 21:07:00", "2020-07-16 21:06:00", "2020-07-17 21:05:00")
suntimes <- as_datetime(suntimes, tz = "America/Vancouver")

ggplot(data = raw[raw$Date >= as_date("2020-07-15") & raw$Date <= as_date("2020-07-19"), ],
       aes(x = Date, y = Light)) +
  geom_point() +
  geom_vline(xintercept = suntimes) + #Original
  geom_vline(xintercept = suntimes - minutes(adjust), colour = "orange") # Adjusted



# Using adjustment
raw_fix <- mutate(raw, Date = Date + minutes(adjust))
write.csv(raw_fix, file.path(out_dir, paste0(ID, "_fixed.csv")), row.names = FALSE)
