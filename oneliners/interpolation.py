## author: Declan Finney
## date 14no2023
## purpose: To interpolate satellite (or could be applied to model) to FAAM lat,lon,time.
## TAGS: satellite ; UM ; coordinates ; interpolate ; FAAM ; aircraft

## all xarrays
## goes_rgd is whatever you want to interpolate, and ds_faam contains lats , lons and times of famm data.
goes_interp = goes_rgd.interp(lon=ds_faam.LON_GIN, lat=ds_faam.LAT_GIN, time=ds_faam.Time,method="linear")
