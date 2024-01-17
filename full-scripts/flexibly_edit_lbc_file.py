# Erin Raif, Jan 2024
# code inspired by Paul Field set_ft_aerosol_lbc.py
# TAGS: LBCS, MULE, CASIM, NESTING SUITE

# Modifies the LBC file for the regional model that is used by the regional forecast process.

# Originally designed to create profiles (where a single value at each level is applied to each x,y point in the LBC)
# but contains functions to easily adapt to broader uses.

# NOTE: TO understand LBC file structure, use accompanying image (how_lbc_files_work.jpg) cause it's a nightmare!

# When using on Monsoon, must first perform:
# module load scitools
# module load um_tools

# Modifies the LBC file for the regional model that is used by the regional forecast process.

import numpy as np
import mule
from mule import ArrayDataProvider as ADP
from mule.operators import ScaleFactorOperator

#### TO EDIT ####

# I/O files - these must be different
input_file = 'small_lbc_test' # this is a minimal file containing two STASH variables at two timesteps with value 0 everywhere
                              # this example file is too big for github but available from Erin, email eeenr@leeds.ac.uk
output_file = 'small_lbc_test_results'

grid_x_pts = 800
grid_y_pts = 800
rimwidth   = 9 # UM default
LBC_x_pts  = 814
LBC_y_pts  = 814
model_levels = 70

# Optional input for creating 'profiles', where LBC has constant value at all x,y on a level (but different for each level):
# This should be an ASCII file in the format:
# STASH number 1
#   level1 value
#   level2 value
#   ...
#   top level value
# STASH number 2
#   as above
profile_file = 'example_lbc_profile_dump.txt' # this is a profile for one STASH only, 36003
profile_levels = 70
profile_stash_nos = 1

#### FUNCTIONS ####
def create_profiles(profile_file, profile_levels, profile_stash_nos, LBC_pts_per_lev):
  """Creates a dictionary containing profile data for a given stash code
  
  Parameters
  ----------
  profile_file : string
    file name of input file
  profile_levels : int
    number of levels in the profile
  profile_stash_nos: int
    number of stash entries to be edited
  LBC_pts_per_lev: int
    total number of points in the LBC file on each level

  Returns
  -------
  dict
    keys are stash numbers, values are profiles as ArrayDataProviders
  """
  f = open(profile_file,'r')
  stash_nums = []
  profile_data = {}
  for stash_index in range(profile_stash_nos):
    profile_values = np.zeros(profile_levels + 1) # extra level required for LBC
    for level_index in range(profile_levels + 1):
      if level_index == 0:
        stash = int(f.readline().strip())
        stash_nums.append(stash)
        # leave first LBC level as 0 (this does not correspond to a model level)
      else:
        val = float(f.readline().strip())
        profile_values[level_index] = val
    profile_data[stash] = profile_values

  LBC_profiles = {}
  for stash in profile_data:
    oned = profile_data[stash]
    twod = np.tile(oned[:,np.newaxis],(1,LBC_pts_per_lev)) # this applies the profile across all x,y points in the LBC, modify here if you want more customised boundaries.
    LBC_profiles[stash] = ADP(twod) # create ArrayDataProvider, object that allows exchange of data into the field
  return LBC_profiles

def check_valid_time(field, time_constraint):
  """Check that a field's time is within any supplied time constraint. Deliberately ugly to avoid np/pd reliance.
  
  Parameters
  ----------
  field:
    field from a mule LBC object
  time_constraint: None, string or 2-tuple (string, string)
    None if no time constraint, string YYYYMMDDTHHMMSS if must be at one time only, 2-tuple of strings in same format
    if the time should be >= start time and < end time.

  Returns
  -------
  bool
    is the field compliant with the time constraint?
  """
  if time_constraint == None:
    return True
  elif len(time_constraint) != 1:
    if (field.lbyr  == int(time_constraint[0:4]) and
      field.lbmon == int(time_constraint[4:6]) and 
      field.lbdat == int(time_constraint[6:8]) and
      field.lbhr  == int(time_constraint[9:11]) and
      field.lbmin == int(time_constraint[11:13]) and
      field.lbsec == int(time_constraint[13:15])):
      return True
    else:
      return False
  elif ((field.lbyr >= int(time_constraint[0][0:4]) and # after start time
        field.lbmon >= int(time_constraint[0][4:6]) and 
        field.lbdat >= int(time_constraint[0][6:8]) and
        field.lbhr  >= int(time_constraint[0][9:11]) and
        field.lbmin >= int(time_constraint[0][11:13]) and
        field.lbsec >= int(time_constraint[0][13:15])) and # before end time
        (field.lbyr >= int(time_constraint[1][0:4]) and
         field.lbmon < int(time_constraint[1][4:6]) and 
         field.lbdat < int(time_constraint[1][6:8]) and
         field.lbhr  < int(time_constraint[1][9:11]) and
         field.lbmin < int(time_constraint[1][11:13]) and
         field.lbsec < int(time_constraint[1][13:15])) ):
    return True
  else:
    return False

def create_custom_boundary_array(val, border, levels=None, arr=None):
  """Make a numpy array that will contain LBC data, or adjust one that already exists.
     Allows user to choose which of the boundaries to adjust, N, E, S or W or any combination of these.
     NOTE: To be able to apply this to a field, this array must be passed to an ArrayDataProvider object.

  Parameters
  ----------
  val: float
    a value to apply
  border: string
    string containing the cardinal directions corresponding to the boundaries to change, e.g. "NEW" to edit the
    Northern, Eastern and Western boundaries. This can be in any order. "NEW"="WNE"  Note that corners belong to
    the Northern and Southern boundaries. If corners need to belong to the Eastern/Western boundaries, this will
    need to be edited to access the co-ordinates described in the image that accompanies this code.
  levels: None or 2-tuple (int, int)
    None if values should apply to all levels. 2-tuple (min model level, max model level) if applies to a subset.
    To access LBC level 1 (not a model level), use 0.
  arr: None or ndarray
    None if no array already exists or array to be operated on. Must be the correct size for the LBC file, i.e.
    dimensions (model levels + 1 x total points in LBC)
  
  Returns
  -------
  ndarray
    Array containing wanted values.

  """
  if levels == None:
    levels = (0, model_levels+1)
  if arr is None:
    arr = np.zeros((model_levels+1, LBC_pts_per_lev)) # very occasionally there should not be an extra model level
  if 'N' in border:
    arr[levels[0]:levels[1], :N_boundary_end] = val
  if 'E' in border:
    arr[levels[0]:levels[1], N_boundary_end:E_boundary_end] = val
  if 'S' in border:
    arr[levels[0]:levels[1], E_boundary_end:S_boundary_end] = val
  if 'W' in border:
    arr[levels[0]:levels[1], S_boundary_end:W_boundary_end] = val

  return arr

def apply_profiles(df, LBC_data_prov, time_constraint = None, save=False):
  """Apply profiles to LBC and optionally save.
  
  Parameters
  ----------
  df: mule LBC object
    new df (copy of the original) that will be optionally saved in this function
  LBC_data_prov: dict
    dictionary where keys are stash numbers and values are profiles as ArrayDataProvider objects 
  time_constraint: None, string or 2-tuple (string, string)
    None if no time constraint, string YYYYMMDDTHHMMSS if must be at one time only, 2-tuple of strings in same format
    if the time should be >= start time and < end time.
  save: bool
    True if LBC is finalised after this and should be saved.

  Returns
  -------
  none
  """
  stash_nums = list(LBC_data_prov.keys())
  for field in df.fields:
    # loop through and update field if their stash number (lbuser4) and time meets the criteria
    if field.lbuser4 in stash_nums and check_valid_time(field, time_constraint):
      field.set_data_provider(LBC_data_prov[field.lbuser4])
      print('setting stash no ',field.lbuser4, 'at',field.lbyr,field.lbmon,field.lbdat,field.lbhr)
  if save:
    df.to_file(output_file)
  return df



#### CODE ####

halo_width = (LBC_x_pts - grid_x_pts) // 2
boundary_width  = rimwidth + halo_width
LBC_pts_per_lev = 2 * boundary_width * (LBC_x_pts + LBC_y_pts - 2 * boundary_width) # is this counting wrong?

N_boundary_end = LBC_x_pts * boundary_width # this is the index of the last point using the UM's index from 1. However, in python, it's the first of the next section.
E_boundary_end = N_boundary_end + (LBC_y_pts - 2 * boundary_width) * boundary_width
S_boundary_end = E_boundary_end + LBC_x_pts * boundary_width
W_boundary_end = LBC_pts_per_lev

boundary_ends = (N_boundary_end, E_boundary_end, S_boundary_end, W_boundary_end)


df = mule.LBCFile.from_file(input_file)
new_df = df.copy() # Create blank copy of the input file, containing no field data (one could probably do this more efficiently on the original but this is safer)
for field in df.fields:
  # Copy current fields to the new one
  new_df.fields.append(field)

#### ADD/REMOVE FUNCTIONS HERE ACCORDING TO NEED. THESE ARE USED AS EXAMPLES ####

# Use file to create profiles for 36003 and apply these at all times.
LBC_profiles = create_profiles(profile_file, profile_levels, profile_stash_nos, LBC_pts_per_lev) 
new_df = apply_profiles(new_df, LBC_profiles) 

# for stash no 36004 make the Northern boundary have value 2, Eastern and Southern 1 on levels 1 to 50.
# On all levels, have the Western as 3. 
arr = create_custom_boundary_array(3, 'W')
arr = create_custom_boundary_array(2, 'N', levels=(1,50), arr=arr)
arr = create_custom_boundary_array(1, 'ES', levels=(1,50), arr=arr)

# Apply these BCs only at 13:00
custom_LBCs = {}
custom_LBCs[36004] = ADP(arr)

new_df = apply_profiles(new_df, custom_LBCs, time_constraint = '20220328T130000', save=True)
