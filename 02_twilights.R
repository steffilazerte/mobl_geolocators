library(SGAT)
library(TwGeos)
library(lubridate) # Better dates
library(ggplot2)   # Slower but easier plots


# Set details ---------------------------------------------
dir <- "MOBL_data/fixed/"

# Offset
offset <- 7 + 12  # Difference between local and UTC


# ID000539 ----------------------------------------------

ID  <- "ID000539"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]

# Get appropriate limits for light data
light_lim <- c(min(raw$Light), max(raw$Light))

## Twilights ------------------------
threshold <- 1.62

g <- ggplot(data = raw[2000:5000,], aes(x = Date, y = Light)) +
  theme_bw() +
  geom_line() +
  geom_point(size = 1, aes(colour = Light)) +
  geom_hline(yintercept = threshold) +
  scale_colour_viridis_c(option = "inferno", end = 0.8)

# Check that it looks good
g
g %+% raw[10000:12500,]
g %+% raw[30000:32500,]
g %+% raw[100000:102500,]
g %+% raw[110000:112500,]

# Process
twl <- preprocessLight(raw,
                       threshold = threshold,
                       offset = offset,
                       zlim = light_lim, stage = 2,
                       gr.Device = "x11")

# Save
write.csv(twl, file.path(dir, paste0(ID, "_twl.csv")))


# ID000547 ----------------------------------------------

ID  <- "ID000547"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]


## Twilights ------------------------
threshold <- 1.57

g <- ggplot(data = raw[2000:5000,], aes(x = Date, y = Light)) +
  theme_bw() +
  geom_line() +
  geom_point(size = 1, aes(colour = Light)) +
  geom_hline(yintercept = threshold) +
  scale_colour_viridis_c(option = "inferno", end = 0.8)

# Check that it looks good
g
g %+% raw[10000:12500,]
g %+% raw[30000:32500,]
g %+% raw[100000:102500,]
g %+% raw[110000:112500,]

# Process
twl <- preprocessLight(raw,
                       threshold = threshold,
                       offset = offset,
                       zlim = light_lim, stage = 2,
                       gr.Device = "x11")

# Save
write.csv(twl, file.path(dir, paste0(ID, "_twl.csv")))


# ID000552 ----------------------------------------------

ID  <- "ID000552"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]


## Twilights ------------------------
threshold <- 1.57

g <- ggplot(data = raw[2000:5000,], aes(x = Date, y = Light)) +
  theme_bw() +
  geom_line() +
  geom_point(size = 1, aes(colour = Light)) +
  geom_hline(yintercept = threshold) +
  scale_colour_viridis_c(option = "inferno", end = 0.8)

# Check that it looks good
g
g %+% raw[10000:12500,]
g %+% raw[30000:32500,]
g %+% raw[100000:102500,]
g %+% raw[110000:112500,]

# Process
twl <- preprocessLight(raw,
                       threshold = threshold,
                       offset = offset,
                       zlim = light_lim,
                       gr.Device = "x11")

# Save
write.csv(twl, file.path(dir, paste0(ID, "_twl.csv")))

# ID000553 ----------------------------------------------

ID  <- "ID000553"

## Load data --------------------------------------------

# Read data
raw <- read.csv(file.path(dir, paste0(ID, "_fixed.csv")))
head(raw)

# Dates
raw$Date <- ymd_hms(raw$Date, tz = "UTC")
raw$Date[1]


## Twilights ------------------------
threshold <- 1.57

g <- ggplot(data = raw[2000:5000,], aes(x = Date, y = Light)) +
  theme_bw() +
  geom_line() +
  geom_point(size = 1, aes(colour = Light)) +
  geom_hline(yintercept = threshold) +
  scale_colour_viridis_c(option = "inferno", end = 0.8)

# Check that it looks good
g
g %+% raw[10000:12500,]
g %+% raw[30000:32500,]
g %+% raw[100000:102500,]
g %+% raw[110000:112500,]

# Process
twl <- preprocessLight(raw,
                       threshold = threshold,
                       offset = offset,
                       zlim = light_lim,
                       gr.Device = "x11")

# Save
write.csv(twl, file.path(dir, paste0(ID, "_twl.csv")))

