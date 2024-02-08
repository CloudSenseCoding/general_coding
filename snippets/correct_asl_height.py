## authors: Ezri Alkilani-Brown and Declan Finney
## date: 8 feb 2024
## purpose: pp model level height is not correct, this code creates the correct height (m) above sea level. It is 3d variable, which accounts for the sigma pressure coordinate and surface topography.
## TAGS: CASIM ; UM ; model level height

import iris
import xarray as xr
 
path = #<path to pp files>
pp_var = #<pp file containing 3d variable>
pp_alt = #<pp file containing surface_altitude >

cube_var = iris.load_cube(path + pp_var, '< 3d variable >')
cube_alt = iris.load_cube(path + pp_alt, 'surface_altitude’) # assuming surface topography is called 'surface_altitude’
 
var = xr.DataArray.from_iris(cube_var)
orog = xr.DataArray.from_iris(cube_alt)
 
# this is the 3d variable of correct height above sea level:
alt_ASL = var[‘level_height’] + var[‘sigma’] * orog # assuming pp file (should) contain level_height and sigma.
