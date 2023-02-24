# Mountain bluebird migration

This is a rough extraction of mountain bluebird migration tracks from geolocator
data.

Based on [Light level geolocation analyses](https://geolocationmanual.vogelwarte.ch/). 

See [Matt Reudink](https://www.mattreudink.com/) for more details on this study.

---

## Scripts
- 00_setup.R - Install packages and get folders setup
- 01_check_date_times.R - Fixes the times and date ranges
- 02_twilights.R - Get twilights
- 03_post_twilights.R - Run calibrations and models etc.
