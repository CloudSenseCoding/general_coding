## a example ncview command to view an animation from a given expt
## ncview 20220730T0000Z_LMagda_km1p5set1_expt1_cloudprops01[5-9].nc 20220730T0000Z_LMagda_km1p5set1_expt1_cloudprops0[2-3]*.nc   &


import iris
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import numpy.ma as ma
import xarray as xr
import re
import gc ## garbage collection


## sticking with this. I think a cloud fraction condition will be more useful
cloudmass_cond = 1e-6 ## as suggested by paul for lowest value aircraft can measure
cf_cond = 0.1

suitedate = "da775_9nov2023/"
ini = "20220730T0000Z"
domain="LMagda"
grid_res="km1p5"
minfilenum = 12
maxfilenum = 35


directory = "/gws/nopw/j04/dcmex/users/dfinney/data/CASIM/"+suitedate
filelist = [ini+"_"+domain+"_"+grid_res+"set1_expt1",ini+"_"+domain+"_"+grid_res+"set1_expt2"
            ,ini+"_"+domain+"_"+grid_res+"set1_expt3",ini+"_"+domain+"_"+grid_res+"set1_expt4"
            ,ini+"_"+domain+"_"+grid_res+"set2_expt5",ini+"_"+domain+"_"+grid_res+"set2_expt6"
            ,ini+"_"+domain+"_"+grid_res+"set2_expt7",ini+"_"+domain+"_"+grid_res+"set2_expt8"]

filelist = [ini+"_"+domain+"_"+grid_res+"set1_expt2"
            ,ini+"_"+domain+"_"+grid_res+"set1_expt3",ini+"_"+domain+"_"+grid_res+"set1_expt4"
            ,ini+"_"+domain+"_"+grid_res+"set2_expt5",ini+"_"+domain+"_"+grid_res+"set2_expt6"
            ,ini+"_"+domain+"_"+grid_res+"set2_expt7",ini+"_"+domain+"_"+grid_res+"set2_expt8"]



##plotting
for ifn, fn in enumerate(filelist):
    print("==========")
    print(directory+fn)
    print("==========")
    #ds_list = [] ## to store for each time
    for it, t in enumerate(np.arange(minfilenum,maxfilenum+1)):
        print("==========")
        print(directory+fn,str(t))
        print("==========")
        
        
        expt_label = "expt"+fn[re.search("expt", fn).end():]
        # read other pp 
        cubelist_pp = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["mass_fraction_of_cloud_liquid_water_in_air","mass_fraction_of_rain_in_air","mass_fraction_of_cloud_ice_crystals_in_air","mass_fraction_of_cloud_ice_in_air","mass_fraction_of_graupel_in_air","m01s02i111","air_temperature","air_pressure"])
        cubelist_num = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["number_of_cloud_droplets_per_kg_of_air","number_of_rain_drops_per_kg_of_air","number_of_ice_particles_per_kg_of_air", "number_of_snow_aggregates_per_kg_of_air","number_of_graupel_particles_per_kg_of_air"])
        cubelist_cf = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["cloud_volume_fraction_in_atmosphere_layer","liquid_cloud_volume_fraction_in_atmosphere_layer","ice_cloud_volume_fraction_in_atmosphere_layer"])
        cubelist_rads = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["toa_incoming_shortwave_flux","toa_outgoing_longwave_flux","toa_outgoing_longwave_flux_assuming_clear_sky","toa_outgoing_shortwave_flux","toa_outgoing_shortwave_flux_assuming_clear_sky"])

        ## make rads xarrays. these will go straight into netcdf, no need for compositing.
        isr = xr.DataArray.from_iris(cubelist_rads[0])
        olr = xr.DataArray.from_iris(cubelist_rads[1])
        olr_clr = xr.DataArray.from_iris(cubelist_rads[2])
        osr = xr.DataArray.from_iris(cubelist_rads[3])
        osr_clr = xr.DataArray.from_iris(cubelist_rads[4])
        
        ## make xarrays
        xr_p = xr.DataArray.from_iris(cubelist_pp[-1])
        xr_t = xr.DataArray.from_iris(cubelist_pp[-2])
        xr_rho_full = xr.DataArray.from_iris(cubelist_pp[-3])
        xr_rho = xr_rho_full[:,xr_rho_full["level_height"] >=5.0,:,:]
        xr_rho = xr_rho.rename('Density used for rad layers (kg/m3)')
        ## rho times are from end of first timestep in 15min, whereas other variables times are from end of 15min!
        ## just going to say rho is the times of the variables. I looked and rho varied in level 30 by less than 0.1% between 15min outputs.
        xr_rho["time"] = xr_p["time"]
        
        xr_cloudmass = xr.DataArray.from_iris(cubelist_pp[0] + cubelist_pp[1] + cubelist_pp[2] + cubelist_pp[3] + cubelist_pp[4])
        ## also arrays for cloud top liq/ice num/mass
        xr_topdrop_num = xr.DataArray.from_iris(cubelist_num[0])
        xr_topdrop_mass = xr.DataArray.from_iris(cubelist_pp[0])
        xr_toprain_num = xr.DataArray.from_iris(cubelist_num[1])
        xr_toprain_mass = xr.DataArray.from_iris(cubelist_pp[1])
        xr_topcrys_num = xr.DataArray.from_iris(cubelist_num[2])
        xr_topcrys_mass = xr.DataArray.from_iris(cubelist_pp[2])
        xr_topsnow_num = xr.DataArray.from_iris(cubelist_num[3])
        xr_topsnow_mass = xr.DataArray.from_iris(cubelist_pp[3])
        xr_topgraupel_num = xr.DataArray.from_iris(cubelist_num[4])
        xr_topgraupel_mass = xr.DataArray.from_iris(cubelist_pp[4])
        ## also arrays for cloud fractions
        xr_cf = xr.DataArray.from_iris(cubelist_cf[0])
        xr_cfliq = xr.DataArray.from_iris(cubelist_cf[1])
        xr_cfice = xr.DataArray.from_iris(cubelist_cf[2])


        ## calculate model_height depths, and we'll get that for top level too.
        level_midpoints =  0.5*(xr_p.level_height[0:-1].values + xr_p.level_height[1:].values)
        level_midpoints = np.concatenate([[0.0], level_midpoints , [np.nan]])
        level_depths = level_midpoints[1:] - level_midpoints[0:-1] 
        # make xarray
        xr_depths = xr.DataArray(data=level_depths
                                 , coords=xr_p.level_height.coords
                                 , attrs=xr_p.level_height.attrs )

        
        ## having a few issues with jasmin killing the script. i guess memory.
        # trying to delete pp object, and tidy
        del cubelist_pp
        del cubelist_num
        del cubelist_cf
        del cubelist_rads
        gc.collect()

        ## make index or where xr_cloudmass > 1e-6
        ## value was suggested by Paul Field to Xinyi for the lowest value that can be measured by aircraft.
        index_cloud = xr.where((xr_cloudmass>cloudmass_cond) & (xr_cf>cf_cond),1,np.nan)

        ## add level_number value to cloudy cells
        modlev_cloud = index_cloud*xr_cloudmass.model_level_number

        ## find max and min
        ## need to subtract because model_level_number is >1. where python wants index from 0 obviously.
        min_modlev_cloud = modlev_cloud.min("model_level_number")-1
        max_modlev_cloud = modlev_cloud.max("model_level_number")-1

        ## need to load these to use them as indexes, and can't use nans
        min_modlev_cloud = min_modlev_cloud.fillna(-1).astype(int).load()
        max_modlev_cloud = max_modlev_cloud.fillna(-1).astype(int).load()

        ## find height, pressure and temperature
        ## because put in -1 instead of nan, need to remove this again. thanks chatgpt!
        cbh = xr_cloudmass.level_height.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        cth = xr_cloudmass.level_height.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        cbp = xr_p.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        ctp = xr_p.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        cbt = xr_t.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        ctt = xr_t.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        cbrho = xr_rho.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        ctrho = xr_rho.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)

        ## also get ice number in top level
        ctDROPnum = xr_topdrop_num.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctDROPmass = xr_topdrop_mass.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctRAINnum = xr_toprain_num.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctRAINmass = xr_toprain_mass.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctCRYSnum = xr_topcrys_num.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctCRYSmass = xr_topcrys_mass.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctSNOWnum = xr_topsnow_num.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctSNOWmass = xr_topsnow_mass.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctGRAUPELnum = xr_topgraupel_num.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctGRAUPELmass = xr_topgraupel_mass.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ## also cloud fractions
        ctCF = xr_cf.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctCFLIQ = xr_cfliq.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ctCFICE = xr_cfice.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)
        ## cloud top level depth
        ctLevDepth = xr_depths.isel(model_level_number=max_modlev_cloud).where(max_modlev_cloud != -1)

        ## cloud base droplets could be handy
        ## maybe w too, but don't have time to implement that now.
        cbDROPnum = xr_topdrop_num.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        cbDROPmass = xr_topdrop_mass.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        cbRAINnum = xr_toprain_num.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        cbRAINmass = xr_toprain_mass.isel(model_level_number=min_modlev_cloud).where(min_modlev_cloud != -1)
        

        ## trying to save memory
        del xr_p
        del xr_t
        del xr_rho_full
        del xr_rho
        del xr_cloudmass
        del xr_topdrop_num
        del xr_topdrop_mass
        del xr_toprain_num
        del xr_toprain_mass
        del xr_topcrys_num
        del xr_topcrys_mass
        del xr_topsnow_num
        del xr_topsnow_mass
        del xr_topgraupel_num
        del xr_topgraupel_mass
        del xr_cf
        del xr_cfliq
        del xr_cfice
        del xr_depths
        gc.collect()
        
        ## save data
        ## save names to metadata
        names_dict_rad = {"isr":"incoming shortwave TOA radiation W m-2"
                          , "olr":"outgoing longwave TOA radiation W m-2"
                          , "olr_clr":"outgoing longwave TOA radiation for clear sky W m-2"
                          , "osr":"outgoing shortwave TOA radiation W m-2"
                          , "osr_clr":"outgoing shortwave TOA radiation for clear sky W m-2"
                          }
        names_dict_time = {"cbh":"cloud base model level height", "cth":"cloud top model level height"
                           , "cbp":"cloud base pressure", "ctp":"cloud top pressure"
                           , "cbt":"cloud base temperature", "ctt":"cloud top temperature"
                           , "cbrho":"cloud base air density kg m-3", "ctrho":"cloud top air density kg m-3"
                           , "ctLevDepth":"cloud top model level depth"
                           , "ctDROPnum":"cloud top liq droplet number kg-1"
                           , "cbDROPmass":"cloud base liq droplet mass kg kg-1"
                           , "cbDROPnum":"cloud base liq droplet number kg-1"
                           , "ctDROPmass":"cloud top liq droplet mass kg kg-1"
                           , "ctRAINnum":"cloud top rain number kg-1"
                           , "ctRAINmass":"cloud top rain mass kg kg-1"
                           , "cbRAINnum":"cloud base rain number kg-1"
                           , "cbRAINmass":"cloud base rain mass kg kg-1"
                           , "ctCRYSnum":"cloud top ice crystal number kg-1"
                           , "ctCRYSmass":"cloud top ice crystal mass kg kg-1"
                           , "ctSNOWnum":"cloud top snow number kg-1"
                           , "ctSNOWmass":"cloud top snow mass kg kg-1"
                           , "ctGRAUPELnum":"cloud top graupel number kg-1"
                           , "ctGRAUPELmass":"cloud top graupel mass kg kg-1"
                           , "ctCF":"total cloud volume fraction"
                           , "ctCFLIQ":"liquid cloud volume fraction"
                           , "ctCFICE":"ice cloud volume fraction"
                           }

        # Merge DataArrays into a single Dataset for this time step
        ## aving radiation variables seperate.
        ## for some reason when doing the merge, they wipe the other variables
        ## it must be something to do with them being 2D to beging with and the others
        ## being composited from 3D. But I can't overcome it.
        ds_rad =  xr.merge([isr.rename("isr"), olr.rename("olr")
                            , olr_clr.rename("olr_clr"), osr.rename("osr")
                            , osr_clr.rename("osr_clr")
                            ], compat='override')
        ds_time = xr.merge([cbh.rename("cbh"), cth.rename("cth")
                            , cbp.rename("cbp"), ctp.rename("ctp")
                            , cbt.rename("cbt"), ctt.rename("ctt")
                            , cbrho.rename("cbrho"), ctrho.rename("ctrho")
                            , ctLevDepth.rename("ctLevDepth")
                            , ctDROPnum.rename("ctDROPnum")
                            , ctDROPmass.rename("ctDROPmass")
                            , ctRAINnum.rename("ctRAINnum")
                            , ctRAINmass.rename("ctRAINmass")
                            , cbDROPnum.rename("cbDROPnum")
                            , cbDROPmass.rename("cbDROPmass")
                            , cbRAINnum.rename("cbRAINnum")
                            , cbRAINmass.rename("cbRAINmass")
                            , ctCRYSnum.rename("ctCRYSnum")
                            , ctCRYSmass.rename("ctCRYSmass")
                            , ctSNOWnum.rename("ctSNOWnum")
                            , ctSNOWmass.rename("ctSNOWmass")
                            , ctGRAUPELnum.rename("ctGRAUPELnum")
                            , ctGRAUPELmass.rename("ctGRAUPELmass")
                            , ctCF.rename("ctCF")
                            , ctCFLIQ.rename("ctCFLIQ")
                            , ctCFICE.rename("ctCFICE")
                            ], compat='override')

        # Update attributes
        for var_name, long_name in names_dict_rad.items():
            ds_rad[var_name].attrs.update({"standard_name": var_name, "long_name": long_name})
        for var_name, long_name in names_dict_time.items():
            ds_time[var_name].attrs.update({"standard_name": var_name, "long_name": long_name})

    
        # Save to NetCDF
        ## whilst this is effectively just duplicating the radiation fields, it takes too long to open them from the 3D pz files, so this makes it easier, I should probably just save to a seperate STASH.
        output_filename_rad = f"{directory}/{fn}_rad_"+str(t).zfill(3)+".nc"
        ds_rad.to_netcdf(output_filename_rad)
        output_filename_time = f"{directory}/{fn}_cloudprops_CF"+str(int(cf_cond*100))+"pcMASS"+str(cloudmass_cond)+"_"+str(t).zfill(3)+".nc"
        ds_time.to_netcdf(output_filename_time)

        ## trying to save memory
        del ds_rad
        del ds_time
        del isr
        del olr
        del olr_clr
        del osr
        del osr_clr
        del ctt
        del cth
        del ctp
        del cbt
        del cbh
        del cbp
        del cbrho
        del ctrho
        del ctLevDepth
        del ctDROPnum
        del ctDROPmass
        del ctRAINnum
        del ctRAINmass
        del cbDROPnum
        del cbDROPmass
        del cbRAINnum
        del cbRAINmass
        del ctCRYSnum
        del ctCRYSmass
        del ctSNOWnum
        del ctSNOWmass
        del ctGRAUPELnum
        del ctGRAUPELmass
        del ctCF
        del ctCFLIQ
        del ctCFICE
        gc.collect()


