library(SGAT)
library(TwGeos)
library(MASS)
library(lubridate) # Better dates
library(ggplot2)   # Slower but easier plots
library(rnaturalearth) # Maps
library(sf) # better spatial
library(dplyr)
library(tidyr)

# Set details ---------------------------------------------
dir <- "MOBL_data/fixed/"
dir_out <- "MOBL_data/final"

# Offset
offset <- 7 + 12  # Difference between local and UTC

# Kamloops, release location
lon_calib <- -120.33
lat_calib <- 50.67

# Map for plots
north_america <- ne_states(country = c("United States of America", "Canada"),
                           returnclass = "sf") %>%
  filter(!name %in% c("Hawaii", "Alaska"))


## Movement models ---------------------
#https://geolocationmanual.vogelwarte.ch/SGAT.html#movement-model
# I think these can be set as a species-level effect
#
# CHOICE POINT - What is reasonable for this species?

beta  <- c(2.2, 0.08) # First is mean(ish), second is variation
matplot(0:100, dgamma(0:100, beta[1], beta[2]),
        type = "l", col = "orange",lty = 1,lwd = 2,ylab = "Density", xlab = "km/h")


# ID000539 ----------------------------------------------

ID  <- "ID000539"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
twl <- read.csv(file.path(dir, paste0(ID, "_twl.csv")))
head(twl)

twl$Twilight <- ymd_hms(twl$Twilight, tz = "UTC") # get the Twilight times back into the POSIX. class format
raw$Date <- ymd_hms(raw$Date, tz = "UTC")

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))


## Twilight tweaks ------------------------------
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 19, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"), )

twl_edits <- twilightEdit(
  twilights = twl,
  offset = offset,
  window = 4,           # two days before and two days after
  outlier.mins = 25,    # difference in mins
  stationary.mins = 60, # are the other surrounding twilights within 25 mins of one another
  plot = TRUE)

twl <- twl_edits
twl <- twl[!twl$Deleted,]

## Calibration --------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#calibration-1
# CHOICE POINTS
# - Which to use, during breeding or overwintering (alt)
# - Alt, which tol (tolerance to use)
# - For both, which dates to use

lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 16, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"))

tm.calib <- as.POSIXct(c("2020-07-15", "2020-09-01"), tz = "UTC")
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise,
                              lon_calib, lat_calib, method = "gamma")

zenith <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

# Compare with alternative
# https://geolocationmanual.vogelwarte.ch/SGAT.html#alternative-calibration
start_date <- "2020-12-21"
end_date <- "2021-04-21"
start <- min(which(as_date(twl$Twilight) == start_date))
end <- max(which(as_date(twl$Twilight) == end_date))

zenith_sd <- GeoLocTools::findHEZenith(twl, tol = 0.08, range=c(start, end))

# Let's try the alternate method, I think there are pretty substantial
# differences in light patterns, so perhaps a candidate for alt.

zenith <- zenith_sd
zenith0 <- zenith0 + abs(zenith-zenith_sd) # From end of Alt. Calib. section


## Initial path ------------------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#initial-path
#
# CHOICE POINT - tol (tolerance). Used to smooth over the equinox.
#                Don't want > 0.18, but need some smoothing.
#                Smaller is always better, but increase if too much noise

path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = 0.08)

# Steffi's workflow for visualizing
coords <- data.frame(time = with_tz(path$time, "UTC"),
                     lon = path$x[, 1],
                     lat = path$x[, 2]) %>%
  mutate(year = year(time))

path_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326) %>%
  mutate(equinox = time %within% interval("2020-09-11", "2020-09-30") |
           time %within% interval("2021-03-11", "2021-03-30"))
path_lines_sf <- path_sf %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot(data = path_sf) +
  geom_sf(data = north_america) +
  geom_sf(aes(colour = equinox)) +
  geom_sf(data = path_lines_sf)+
  facet_wrap(~year) +
  geom_sf(data = path_sf[1,], colour = "red") +
  geom_sf(data = st_sfc(st_point(c(lon_calib, lat_calib)), crs = 4326), colour = "blue") +
  labs(caption = "Blue dot is Kamloops; Red dot is first calculated location")


## Define known locations ------------------
# CHOICE POINT - What dates do you *know* they were on site for?
coords$fixed <- FALSE
coords$fixed[coords$time <= as_date("2020-08-15")] <- TRUE  # Before migration
coords$fixed[(nrow(coords) - 10):nrow(coords)] <- TRUE     # After, assuming 10 points in Kamloops

head(coords)

# Update the lat/lon for those points
coords$lat[coords$fixed] <- lat_calib
coords$lon[coords$fixed] <- lon_calib

## Model - Initial Run -------------------------

# NOTE: I skipped the land mask section

# Initial Run to get starting values
x0 <- as.matrix(coords[, c("lon", "lat")])
z0 <- trackMidpts(x0)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,  # From calibration
                        beta = beta,    # From movement model
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(x0))
proposal_z <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(z0))

fit <- estelleMetropolis(model, proposal_x, proposal_z, iters = 1000, thin = 20)

# https://geolocationmanual.vogelwarte.ch/SGAT.html#the-estelle-model
# "Once the chain meets the positivity constraint, the next step is to tune the proposal distributions. "
# WHAT DOES THIS MEAN? HOW DO WE KNOW?

## Model - Short runs to tune proposal ----------------------
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "Gamma",
                        alpha = alpha,
                        beta = beta,
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl))
proposal_z <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl) - 1)

for (k in 1:3) {
  fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                           z0 = chainLast(fit$z), iters = 300, thin = 20)

  proposal_x <- mvnorm(chainCov(fit$x), s = 0.2)
  proposal_z <- mvnorm(chainCov(fit$z), s = 0.2)
}

# Check samples

# CHOICE POINT - Is this okay?
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!coords$fixed, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!coords$fixed, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)


## Model - Final Run -----------------
# CHOICE POINTS
# - More iterations?

proposal_x <- mvnorm(chainCov(fit$x), s = 0.25)
proposal_z <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

## Results -------------------------------
# CHOICE POINTS
# - Get rid of wonky points?
# - Discard any chains
sm <- locationSummary(fit$z, time = fit$model$time)
head(sm)

summary(sm)

# Clean up - Try to omit some crazy variation in lat/lon
sm <- filter(sm, Lat.sd < 8, Lon.sd < 8)

xlim <- range(coords$lon + c(-35, 35))
ylim <- range(coords$lat + c(-15, 15))

r <- raster::raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
            xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = 4326)
s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- SGAT::slice(s, sliceIndices(s))

# Convert to starts to sf for plotting
sk_st <- stars::st_as_stars(sk) %>%
  st_as_sf() %>%
  mutate(layer = if_else(layer < 1, NA_real_, layer)) # If less than 1, don't plot

tracks <- dplyr::select(sm, lon = "Lon.50%", lat = "Lat.50%", Lat.sd, Lon.sd, "Time1") %>%
  mutate(year = year(Time1)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
paths <- tracks %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

g <- ggplot(sk_st) +
  theme_bw() +
  geom_sf(aes(fill = layer), colour = NA) +
  geom_sf(data = north_america, fill = NA, colour = "black", linewidth = 1) +
  scale_fill_viridis_c(na.value = NA)

g

g + geom_sf(data = tracks) +
  geom_sf(data = paths) +
  facet_wrap(~year)


# Lat/lon plots
sm_long <- sm %>%
  pivot_longer(cols = -c("Time1", "Time2"), names_to = c("Type", "Stat"),
               values_to = "Value", names_pattern = "(Lat|Lon)\\.(.+)") %>%
  pivot_wider(names_from = Stat, values_from = Value)

ggplot(data = sm_long, aes(x = Time1, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~Type, ncol = 1, scales = "free_y")


write.csv(sm, file.path(dir_out, paste0(ID, "_final.csv")), row.names = FALSE)
saveRDS(fit, file.path(dir_out, paste0(ID, "_model_fit.rds")))


# ID000547 FIX ME!!! ----------------------------------------------

ID  <- "ID000547"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
twl <- read.csv(file.path(dir, paste0(ID, "_twl.csv")))
head(twl)

twl$Twilight <- ymd_hms(twl$Twilight, tz = "UTC") # get the Twilight times back into the POSIX. class format
raw$Date <- ymd_hms(raw$Date, tz = "UTC")

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Twilight tweaks ------------------------------
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 19, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"), )

twl_edit <- twilightEdit(
  twilights = twl,
  offset = offset,
  window = 4,           # two days before and two days after
  outlier.mins = 25,    # difference in mins
  stationary.mins = 60, # are the other surrounding twilights within 25 mins of one another
  plot = TRUE)

twl <- twl_edit
twl <- twl[!twl$Deleted,]

## Calibration --------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#calibration-1
# CHOICE POINTS
# - Which to use, during breeding or overwintering (alt)
# - Alt, which tol (tolerance to use)
# - For both, which dates to use

lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 16, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"))

tm.calib <- as.POSIXct(c("2020-07-01", "2020-08-25"), tz = "UTC")
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise,
                              lon_calib, lat_calib, method = "gamma")

zenith <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

# Compare with alternative
# https://geolocationmanual.vogelwarte.ch/SGAT.html#alternative-calibration
start_date <- "2020-12-21"
end_date <- "2021-04-21"
start <- min(which(as_date(twl$Twilight) == start_date))
end <- max(which(as_date(twl$Twilight) == end_date))

zenith_sd <- GeoLocTools::findHEZenith(twl, tol = 0.08, range=c(start, end))

zenith_sd
zenith0

# Let's try using both methods, I think there are pretty substantial
# differences in light patterns, so perhaps a candidate for alt.

zenith <- zenith_sd
zenith0 <- zenith0 + abs(zenith-zenith_sd) # From end of Alt. Calib. section


## Initial path ------------------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#initial-path
#
# CHOICE POINT - tol (tolerance). Used to smooth over the equinox.
#                Don't want > 0.18, but need some smoothing.
#                Smaller is always better, but increase if too much noise

path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = 0.08)

# Steffi's workflow for visualizing
coords <- data.frame(time = with_tz(path$time, "UTC"),
                     lon = path$x[, 1],
                     lat = path$x[, 2]) %>%
  mutate(year = year(time))

path_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326) %>%
  mutate(equinox = time %within% interval("2020-09-11", "2020-09-30") |
           time %within% interval("2021-03-11", "2021-03-30"))
path_lines_sf <- path_sf %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot(data = path_sf) +
  geom_sf(data = north_america) +
  geom_sf(aes(colour = equinox)) +
  geom_sf(data = path_lines_sf)+
  facet_wrap(~year) +
  geom_sf(data = path_sf[1,], colour = "red") +
  geom_sf(data = st_sfc(st_point(c(lon_calib, lat_calib)), crs = 4326), colour = "blue") +
  labs(caption = "Blue dot is Kamloops; Red dot is first calculated location")


## Define known locations ------------------
# CHOICE POINT - What dates do you *know* they were on site for?
coords$fixed <- FALSE
coords$fixed[coords$time <= as_date("2020-08-15")] <- TRUE  # Before migration
coords$fixed[(nrow(coords) - 10):nrow(coords)] <- TRUE     # After, assuming 10 points in Kamloops

head(coords)

# Update the lat/lon for those points
coords$lat[coords$fixed] <- lat_calib
coords$lon[coords$fixed] <- lon_calib

## Model - Initial Run -------------------------

# NOTE: I skipped the land mask section

# Initial Run to get starting values
x0 <- as.matrix(coords[, c("lon", "lat")])
z0 <- trackMidpts(x0)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,  # From calibration
                        beta = beta,    # From movement model
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(x0))
proposal_z <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(z0))

fit <- estelleMetropolis(model, proposal_x, proposal_z, iters = 1000, thin = 20)

## Model - Short runs to tune proposal ----------------------
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "Gamma",
                        alpha = alpha,
                        beta = beta,
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl))
proposal_z <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl) - 1)

for (k in 1:3) {
  fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                           z0 = chainLast(fit$z), iters = 300, thin = 20)

  proposal_x <- mvnorm(chainCov(fit$x), s = 0.2)
  proposal_z <- mvnorm(chainCov(fit$z), s = 0.2)
}

# Check samples

# CHOICE POINT - Is this okay?
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!coords$fixed, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!coords$fixed, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)


## Model - Final Run -----------------
# CHOICE POINTS
# - More iterations?

proposal_x <- mvnorm(chainCov(fit$x), s = 0.25)
proposal_z <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

## Results -------------------------------
# CHOICE POINTS
# - Get rid of wonky points?
# - Discard any chains
sm <- locationSummary(fit$z, time = fit$model$time)
head(sm)

summary(sm)

# Clean up - Try to omit some crazy variation in lat/lon
sm <- filter(sm, Lat.sd < 8, Lon.sd < 8)

xlim <- range(coords$lon + c(-35, 35))
ylim <- range(coords$lat + c(-15, 15))

r <- raster::raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
                    xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = 4326)
s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- SGAT::slice(s, sliceIndices(s))

# Convert to starts to sf for plotting
sk_st <- stars::st_as_stars(sk) %>%
  st_as_sf() %>%
  mutate(layer = if_else(layer < 1, NA_real_, layer))  # If less than 1, don't plot

tracks <- dplyr::select(sm, lon = "Lon.50%", lat = "Lat.50%", Lat.sd, Lon.sd, "Time1") %>%
  mutate(year = year(Time1)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
paths <- tracks %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

g <- ggplot(sk_st) +
  theme_bw() +
  geom_sf(aes(fill = layer), colour = NA) +
  geom_sf(data = north_america, fill = NA, colour = "black", linewidth = 1) +
  scale_fill_viridis_c(na.value = NA)

g

g + geom_sf(data = tracks) +
  geom_sf(data = paths) +
  facet_wrap(~year)


# Lat/lon plots
sm_long <- sm %>%
  pivot_longer(cols = -c("Time1", "Time2"), names_to = c("Type", "Stat"),
               values_to = "Value", names_pattern = "(Lat|Lon)\\.(.+)") %>%
  pivot_wider(names_from = Stat, values_from = Value)

ggplot(data = sm_long, aes(x = Time1, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~Type, ncol = 1, scales = "free_y")


write.csv(sm, file.path(dir_out, paste0(ID, "_final.csv")), row.names = FALSE)
saveRDS(fit, file.path(dir_out, paste0(ID, "_model_fit.rds")))







# ID000552 ----------------------------------------------

ID  <- "ID000552"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
twl <- read.csv(file.path(dir, paste0(ID, "_twl.csv")))
head(twl)

twl$Twilight <- ymd_hms(twl$Twilight, tz = "UTC") # get the Twilight times back into the POSIX. class format
raw$Date <- ymd_hms(raw$Date, tz = "UTC")

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Twilight tweaks ------------------------------
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 19, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"), )

twl_edit <- twilightEdit(
  twilights = twl,
  offset = offset,
  window = 4,           # two days before and two days after
  outlier.mins = 25,    # difference in mins
  stationary.mins = 60, # are the other surrounding twilights within 25 mins of one another
  plot = TRUE)

twl <- twl_edit
twl <- twl[!twl$Deleted,]

## Calibration --------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#calibration-1
# CHOICE POINTS
# - Which to use, during breeding or overwintering (alt)
# - Alt, which tol (tolerance to use)
# - For both, which dates to use

lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 16, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"))

tm.calib <- as.POSIXct(c("2020-06-05", "2020-07-25"), tz = "UTC")
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise,
                              lon_calib, lat_calib, method = "gamma")

zenith <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

# Compare with alternative
# https://geolocationmanual.vogelwarte.ch/SGAT.html#alternative-calibration
start_date <- "2020-12-21"
end_date <- "2021-04-21"
start <- min(which(as_date(twl$Twilight) == start_date))
end <- max(which(as_date(twl$Twilight) == end_date))

zenith_sd <- GeoLocTools::findHEZenith(twl, tol = 0.08, range=c(start, end))

zenith_sd
zenith0

# Let's try using both methods, I think there are pretty substantial
# differences in light patterns, so perhaps a candidate for alt.

zenith <- zenith_sd
zenith0 <- zenith0 + abs(zenith-zenith_sd) # From end of Alt. Calib. section


## Initial path ------------------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#initial-path
#
# CHOICE POINT - tol (tolerance). Used to smooth over the equinox.
#                Don't want > 0.18, but need some smoothing.
#                Smaller is always better, but increase if too much noise

path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = 0.08)

# Steffi's workflow for visualizing
coords <- data.frame(time = with_tz(path$time, "UTC"),
                     lon = path$x[, 1],
                     lat = path$x[, 2]) %>%
  mutate(year = year(time))

path_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326) %>%
  mutate(equinox = time %within% interval("2020-09-11", "2020-09-30") |
           time %within% interval("2021-03-11", "2021-03-30"))
path_lines_sf <- path_sf %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot(data = path_sf) +
  geom_sf(data = north_america) +
  geom_sf(aes(colour = equinox)) +
  geom_sf(data = path_lines_sf)+
  facet_wrap(~year) +
  geom_sf(data = path_sf[1,], colour = "red") +
  geom_sf(data = st_sfc(st_point(c(lon_calib, lat_calib)), crs = 4326), colour = "blue") +
  labs(caption = "Blue dot is Kamloops; Red dot is first calculated location")


## Define known locations ------------------
# CHOICE POINT - What dates do you *know* they were on site for?
coords$fixed <- FALSE
coords$fixed[coords$time <= as_date("2020-08-15")] <- TRUE  # Before migration
coords$fixed[(nrow(coords) - 10):nrow(coords)] <- TRUE     # After, assuming 10 points in Kamloops

head(coords)

# Update the lat/lon for those points
coords$lat[coords$fixed] <- lat_calib
coords$lon[coords$fixed] <- lon_calib

## Model - Initial Run -------------------------

# NOTE: I skipped the land mask section

# Initial Run to get starting values
x0 <- as.matrix(coords[, c("lon", "lat")])
z0 <- trackMidpts(x0)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,  # From calibration
                        beta = beta,    # From movement model
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(x0))
proposal_z <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(z0))

fit <- estelleMetropolis(model, proposal_x, proposal_z, iters = 1000, thin = 20)

## Model - Short runs to tune proposal ----------------------
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "Gamma",
                        alpha = alpha,
                        beta = beta,
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl))
proposal_z <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl) - 1)

for (k in 1:3) {
  fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                           z0 = chainLast(fit$z), iters = 300, thin = 20)

  proposal_x <- mvnorm(chainCov(fit$x), s = 0.2)
  proposal_z <- mvnorm(chainCov(fit$z), s = 0.2)
}

# Check samples

# CHOICE POINT - Is this okay?
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!coords$fixed, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!coords$fixed, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)


## Model - Final Run -----------------
# CHOICE POINTS
# - More iterations?

proposal_x <- mvnorm(chainCov(fit$x), s = 0.25)
proposal_z <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

## Results -------------------------------
# CHOICE POINTS
# - Get rid of wonky points?
# - Discard any chains
sm <- locationSummary(fit$z, time = fit$model$time)
head(sm)

summary(sm)

# Clean up - Try to omit some crazy variation in lat/lon
sm <- filter(sm, Lat.sd < 8, Lon.sd < 8)

xlim <- range(coords$lon + c(-35, 35))
ylim <- range(coords$lat + c(-15, 15))

r <- raster::raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
                    xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = 4326)
s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- SGAT::slice(s, sliceIndices(s))

# Convert to starts to sf for plotting
sk_st <- stars::st_as_stars(sk) %>%
  st_as_sf() %>%
  mutate(layer = if_else(layer < 1, NA_real_, layer))  # If less than 1, don't plot

tracks <- dplyr::select(sm, lon = "Lon.50%", lat = "Lat.50%", Lat.sd, Lon.sd, "Time1") %>%
  mutate(year = year(Time1)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
paths <- tracks %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

g <- ggplot(sk_st) +
  theme_bw() +
  geom_sf(aes(fill = layer), colour = NA) +
  geom_sf(data = north_america, fill = NA, colour = "black", linewidth = 1) +
  scale_fill_viridis_c(na.value = NA)

g

g + geom_sf(data = tracks) +
  geom_sf(data = paths) +
  facet_wrap(~year)


# Lat/lon plots
sm_long <- sm %>%
  pivot_longer(cols = -c("Time1", "Time2"), names_to = c("Type", "Stat"),
               values_to = "Value", names_pattern = "(Lat|Lon)\\.(.+)") %>%
  pivot_wider(names_from = Stat, values_from = Value)

ggplot(data = sm_long, aes(x = Time1, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~Type, ncol = 1, scales = "free_y")


write.csv(sm, file.path(dir_out, paste0(ID, "_final.csv")), row.names = FALSE)
saveRDS(fit, file.path(dir_out, paste0(ID, "_model_fit.rds")))


# ID000553 ----------------------------------------------

ID  <- "ID000553"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
twl <- read.csv(file.path(dir, paste0(ID, "_twl.csv")))
head(twl)

twl$Twilight <- ymd_hms(twl$Twilight, tz = "UTC") # get the Twilight times back into the POSIX. class format
raw$Date <- ymd_hms(raw$Date, tz = "UTC")

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Twilight tweaks ------------------------------
lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

# Check the lat/lon where geos were deployed, can clearly see offset with migration
tsimageDeploymentLines(raw$Date, lon = lon_calib, lat = lat_calib,
                       offset = offset, lwd = 3, col = adjustcolor("orange", alpha.f = 0.5))

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 19, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"), )

twl_edit <- twilightEdit(
  twilights = twl,
  offset = offset,
  window = 4,           # two days before and two days after
  outlier.mins = 25,    # difference in mins
  stationary.mins = 60, # are the other surrounding twilights within 25 mins of one another
  plot = TRUE)

twl <- twl_edit
twl <- twl[!twl$Deleted,]

## Calibration --------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#calibration-1
# CHOICE POINTS
# - Which to use, during breeding or overwintering (alt)
# - Alt, which tol (tolerance to use)
# - For both, which dates to use

lightImage(tagdata = raw,
           offset = offset,
           zlim = light_lim)

tsimagePoints(twl$Twilight ,
              offset = offset,
              pch = 16, cex = 1.2,
              col = ifelse(twl$Rise, "firebrick", "cornflowerblue"))

tm.calib <- as.POSIXct(c("2020-07-01", "2020-08-25"), tz = "UTC")
abline(v = tm.calib, lwd = 2, lty = 2, col = "orange")

d_calib <- subset(twl, Twilight>=tm.calib[1] & Twilight<=tm.calib[2])
calib <- thresholdCalibration(d_calib$Twilight, d_calib$Rise,
                              lon_calib, lat_calib, method = "gamma")

zenith <- calib[1]
zenith0 <- calib[2]
alpha <- calib[3:4]

# Compare with alternative
# https://geolocationmanual.vogelwarte.ch/SGAT.html#alternative-calibration
start_date <- "2020-12-21"
end_date <- "2021-04-21"
start <- min(which(as_date(twl$Twilight) == start_date))
end <- max(which(as_date(twl$Twilight) == end_date))

zenith_sd <- GeoLocTools::findHEZenith(twl, tol = 0.08, range=c(start, end))

zenith_sd
zenith0

# Let's try using both methods, I think there are pretty substantial
# differences in light patterns, so perhaps a candidate for alt.

zenith <- zenith_sd
zenith0 <- zenith0 + abs(zenith-zenith_sd) # From end of Alt. Calib. section


## Initial path ------------------------------
# https://geolocationmanual.vogelwarte.ch/SGAT.html#initial-path
#
# CHOICE POINT - tol (tolerance). Used to smooth over the equinox.
#                Don't want > 0.18, but need some smoothing.
#                Smaller is always better, but increase if too much noise

path <- thresholdPath(twl$Twilight, twl$Rise, zenith = zenith, tol = 0.08)

# Steffi's workflow for visualizing
coords <- data.frame(time = with_tz(path$time, "UTC"),
                     lon = path$x[, 1],
                     lat = path$x[, 2]) %>%
  mutate(year = year(time))

path_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326) %>%
  mutate(equinox = time %within% interval("2020-09-11", "2020-09-30") |
           time %within% interval("2021-03-11", "2021-03-30"))
path_lines_sf <- path_sf %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot(data = path_sf) +
  geom_sf(data = north_america) +
  geom_sf(aes(colour = equinox)) +
  geom_sf(data = path_lines_sf)+
  facet_wrap(~year) +
  geom_sf(data = path_sf[1,], colour = "red") +
  geom_sf(data = st_sfc(st_point(c(lon_calib, lat_calib)), crs = 4326), colour = "blue") +
  labs(caption = "Blue dot is Kamloops; Red dot is first calculated location")


## Define known locations ------------------
# CHOICE POINT - What dates do you *know* they were on site for?
coords$fixed <- FALSE
coords$fixed[coords$time <= as_date("2020-08-15")] <- TRUE  # Before migration
coords$fixed[(nrow(coords) - 10):nrow(coords)] <- TRUE     # After, assuming 10 points in Kamloops

head(coords)

# Update the lat/lon for those points
coords$lat[coords$fixed] <- lat_calib
coords$lon[coords$fixed] <- lon_calib

## Model - Initial Run -------------------------

# NOTE: I skipped the land mask section

# Initial Run to get starting values
x0 <- as.matrix(coords[, c("lon", "lat")])
z0 <- trackMidpts(x0)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "ModifiedGamma",
                        alpha = alpha,  # From calibration
                        beta = beta,    # From movement model
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(x0))
proposal_z <- mvnorm(S = diag(c(0.0025, 0.0025)), n = nlocation(z0))

fit <- estelleMetropolis(model, proposal_x, proposal_z, iters = 1000, thin = 20)

## Model - Short runs to tune proposal ----------------------
x0 <- chainLast(fit$x)
z0 <- chainLast(fit$z)

model <- thresholdModel(twilight = twl$Twilight,
                        rise = twl$Rise,
                        twilight.model = "Gamma",
                        alpha = alpha,
                        beta = beta,
                        x0 = x0,
                        z0 = z0,
                        zenith = zenith0,
                        fixedx = coords$fixed)

proposal_x <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl))
proposal_z <- mvnorm(S = diag(c(0.005, 0.005)), n = nrow(twl) - 1)

for (k in 1:3) {
  fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                           z0 = chainLast(fit$z), iters = 300, thin = 20)

  proposal_x <- mvnorm(chainCov(fit$x), s = 0.2)
  proposal_z <- mvnorm(chainCov(fit$z), s = 0.2)
}

# Check samples

# CHOICE POINT - Is this okay?
opar <- par(mfrow = c(2, 1), mar = c(3, 5, 2, 1) + 0.1)
matplot(t(fit$x[[1]][!coords$fixed, 1, ]), type = "l", lty = 1, col = "dodgerblue", ylab = "Lon")
matplot(t(fit$x[[1]][!coords$fixed, 2, ]), type = "l", lty = 1, col = "firebrick", ylab = "Lat")
par(opar)


## Model - Final Run -----------------
# CHOICE POINTS
# - More iterations?

proposal_x <- mvnorm(chainCov(fit$x), s = 0.25)
proposal_z <- mvnorm(chainCov(fit$z), s = 0.25)

fit <- estelleMetropolis(model, proposal_x, proposal_z, x0 = chainLast(fit$x),
                         z0 = chainLast(fit$z), iters = 1000, thin = 20)

## Results -------------------------------
# CHOICE POINTS
# - Get rid of wonky points?
# - Discard any chains
sm <- locationSummary(fit$z, time = fit$model$time)
head(sm)

summary(sm)

# Clean up - Try to omit some crazy variation in lat/lon
sm <- filter(sm, Lat.sd < 8, Lon.sd < 8)

xlim <- range(coords$lon + c(-35, 35))
ylim <- range(coords$lat + c(-15, 15))

r <- raster::raster(nrows = 2 * diff(ylim), ncols = 2 * diff(xlim), xmn = xlim[1]-5,
                    xmx = xlim[2]+5, ymn = ylim[1]-5, ymx = ylim[2]+5, crs = 4326)
s <- slices(type = "intermediate", breaks = "week", mcmc = fit, grid = r)
sk <- SGAT::slice(s, sliceIndices(s))

# Convert to starts to sf for plotting
sk_st <- stars::st_as_stars(sk) %>%
  st_as_sf() %>%
  mutate(layer = if_else(layer < 1, NA_real_, layer))  # If less than 1, don't plot

tracks <- dplyr::select(sm, lon = "Lon.50%", lat = "Lat.50%", Lat.sd, Lon.sd, "Time1") %>%
  mutate(year = year(Time1)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
paths <- tracks %>%
  group_by(year) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

g <- ggplot(sk_st) +
  theme_bw() +
  geom_sf(aes(fill = layer), colour = NA) +
  geom_sf(data = north_america, fill = NA, colour = "black", linewidth = 1) +
  scale_fill_viridis_c(na.value = NA)

g

g + geom_sf(data = tracks) +
  geom_sf(data = paths) +
  facet_wrap(~year)


# Lat/lon plots
sm_long <- sm %>%
  pivot_longer(cols = -c("Time1", "Time2"), names_to = c("Type", "Stat"),
               values_to = "Value", names_pattern = "(Lat|Lon)\\.(.+)") %>%
  pivot_wider(names_from = Stat, values_from = Value)

ggplot(data = sm_long, aes(x = Time1, y = `50%`)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  geom_line() +
  geom_point(size = 1) +
  facet_wrap(~Type, ncol = 1, scales = "free_y")


write.csv(sm, file.path(dir_out, paste0(ID, "_final.csv")), row.names = FALSE)
saveRDS(fit, file.path(dir_out, paste0(ID, "_model_fit.rds")))


