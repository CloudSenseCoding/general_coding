## author: Declan Finney
## date 14no2023
## purpose: nested unifed model grid is generally on rotated grid_latitude and grid_longitude. This needs to be unrotated to know the true lat and lons of the grid cells. This bit of code calculates those.
## TAGS: CASIM ; UM ; coordinates ; grid

sample_cube = iris.load_cube("20220730T0000Z_LMagda_km1p5set1_expt1_pz011.pp", ["toa_incoming_shortwave_flux"]) ## this can be anything containing grid_lats/loons and coord_system
## add lats and lons
# rotate coordinates
glons = sample_cube.coord("grid_longitude").points
glats = sample_cube.coord("grid_latitude").points
x, y = np.meshgrid(glons, glats)
cs = sample_cube.coord_system() # not sure where coord_system is in the xarray
lons_full, lats_full = iris.analysis.cartography.unrotate_pole(x, y, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

## You can add these to an xarray with 
ds = ds.assign_coords({"lon":(["grid_latitude","grid_longitude"], lons_full), "lat": (["grid_latitude","grid_longitude"], lats_full)})
