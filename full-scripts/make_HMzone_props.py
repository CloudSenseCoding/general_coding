## a example ncview command to view an animation from a given expt
## ncview 20220730T0000Z_LMagda_km1p5set1_expt1_cloudprops01[5-9].nc 20220730T0000Z_LMagda_km1p5set1_expt1_cloudprops0[2-3]*.nc   &

## dlf13nov2023
## this script is developed from make_cloudtopprops to instead calculate average properties with the HM zone -2.5 to -7.5C. In addition it will calculate properties in -12.5 to -17.5C. It additionally calculates mean and maxupdarught velocity in these zones


import iris
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import numpy.ma as ma
import xarray as xr
import re
import gc ## garbage collection


## conditions
cloudmass_cond = 1e-6 ## as suggested by paul for lowest value aircraft can measure
cf_cond = 0.1 ## arbitrary cloud fraction condition to ignore small cloud frac grid boxes

## temp ranges to average over
temp_ranges = [[265.65 , 270.65] , [255.65 , 260.65]]


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
filelist = [ini+"_"+domain+"_"+grid_res+"set1_expt3",ini+"_"+domain+"_"+grid_res+"set1_expt4"
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
        print(directory+fn, str(t))
        print("==========")
        
        
        expt_label = "expt"+fn[re.search("expt", fn).end():]
        # read other pp 
        cubelist_pp = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["mass_fraction_of_cloud_liquid_water_in_air","mass_fraction_of_rain_in_air","mass_fraction_of_cloud_ice_crystals_in_air","mass_fraction_of_cloud_ice_in_air","mass_fraction_of_graupel_in_air","upward_air_velocity", "m01s02i111","air_temperature","air_pressure"])
        cubelist_num = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["number_of_cloud_droplets_per_kg_of_air","number_of_rain_drops_per_kg_of_air","number_of_ice_particles_per_kg_of_air", "number_of_snow_aggregates_per_kg_of_air","number_of_graupel_particles_per_kg_of_air"])
        cubelist_cf = iris.load(directory+fn+"_pz"+str(t).zfill(3)+".pp", ["cloud_volume_fraction_in_atmosphere_layer","liquid_cloud_volume_fraction_in_atmosphere_layer","ice_cloud_volume_fraction_in_atmosphere_layer"])

        
        ## make xarrays
        xr_p = xr.DataArray.from_iris(cubelist_pp[-1])
        xr_t = xr.DataArray.from_iris(cubelist_pp[-2])
        xr_w = xr.DataArray.from_iris(cubelist_pp[-4])
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

        ## having a few issues with jasmin killing the script. i guess memory.
        # trying to delete pp object, and tidy
        del cubelist_pp
        del cubelist_num
        del cubelist_cf
        gc.collect()


        for Trange in temp_ranges:
            print("==========")
            print(directory+fn, str(t), str(Trange))
            print("==========")
            Tindex = xr.where((xr_t >= Trange[0]) & (xr_t <=Trange[1]), 1, np.nan)

            T_mn = xr_t.where(Tindex==1).mean("model_level_number")
            h_mn = xr_t.level_height.where(Tindex==1).mean("model_level_number")
            approxdepth = xr_t.level_height.where(Tindex==1).max("model_level_number") - xr_t.level_height.where(Tindex==1).min("model_level_number") ## doesn't account for bottom halflev or top halflev
            p_mn = xr_p.where(Tindex==1).mean("model_level_number")
            rho_mn = xr_rho.where(Tindex==1).mean("model_level_number")
            w_mn = xr_w.where(Tindex==1).mean("model_level_number")
            w_max = xr_w.where(Tindex==1).max("model_level_number")

            DROPnum_mn =  xr_topdrop_num.where(Tindex==1).mean("model_level_number")
            DROPmass_mn =  xr_topdrop_mass.where(Tindex==1).mean("model_level_number")
            RAINnum_mn =  xr_toprain_num.where(Tindex==1).mean("model_level_number")
            RAINmass_mn =  xr_toprain_mass.where(Tindex==1).mean("model_level_number")
            CRYSnum_mn = xr_topcrys_num.where(Tindex==1).mean("model_level_number")
            CRYSmass_mn = xr_topcrys_mass.where(Tindex==1).mean("model_level_number")
            SNOWnum_mn = xr_topsnow_num.where(Tindex==1).mean("model_level_number")
            SNOWmass_mn = xr_topsnow_mass.where(Tindex==1).mean("model_level_number")
            GRAUPELnum_mn = xr_topgraupel_num.where(Tindex==1).mean("model_level_number")
            GRAUPELmass_mn = xr_topgraupel_mass.where(Tindex==1).mean("model_level_number")

            CF_mn = xr_cf.where(Tindex==1).mean("model_level_number")
            CFLIQ_mn = xr_cfliq.where(Tindex==1).mean("model_level_number")
            CFICE_mn = xr_cfice.where(Tindex==1).mean("model_level_number")

            ## save data
            ## save names to metadata
            names_dict = {"T_mn":"mean temperature for levels used in mean", "h_mn":"mean height"
                          , "p_mn":"mean pressure", "rho_mn":"mean density"
                          , "w_mn":"mean updraught velocity", "w_max":"max updraught velocity"
                          , "approxdepth":"max level - min level height for levels used in mean"
                          , "DROPnum_mn":"cloud top liq droplet number kg-1"
                          , "DROPmass_mn":"cloud top liq droplet mass kg kg-1"
                          , "RAINnum_mn":"cloud top rain number kg-1"
                          , "RAINmass_mn":"cloud top rain mass kg kg-1"
                          , "CRYSnum_mn":"cloud top ice crystal number kg-1"
                          , "CRYSmass_mn":"cloud top ice crystal mass kg kg-1"
                          , "SNOWnum_mn":"cloud top snow number kg-1"
                          , "SNOWmass_mn":"cloud top snow mass kg kg-1"
                          , "GRAUPELnum_mn":"cloud top graupel number kg-1"
                          , "GRAUPELmass_mn":"cloud top graupel mass kg kg-1"
                          , "CF_mn":"total cloud volume fraction"
                          , "CFLIQ_mn":"liquid cloud volume fraction"
                          , "CFICE_mn":"ice cloud volume fraction"
                          }

            # Merge DataArrays into a single Dataset for this time step
            ds_time = xr.merge([T_mn.rename("T_mn"), h_mn.rename("h_mn")
                                , p_mn.rename("p_mn"), rho_mn.rename("rho_mn")
                                , w_mn.rename("w_mn"), w_max.rename("w_max")
                                , approxdepth.rename("approxdepth")
                                , DROPnum_mn.rename("DROPnum_mn")
                                , DROPmass_mn.rename("DROPmass_mn")
                                , RAINnum_mn.rename("RAINnum_mn")
                                , RAINmass_mn.rename("RAINmass_mn")
                                , CRYSnum_mn.rename("CRYSnum_mn")
                                , CRYSmass_mn.rename("CRYSmass_mn")
                                , SNOWnum_mn.rename("SNOWnum_mn")
                                , SNOWmass_mn.rename("SNOWmass_mn")
                                , GRAUPELnum_mn.rename("GRAUPELnum_mn")
                                , GRAUPELmass_mn.rename("GRAUPELmass_mn")
                                , CF_mn.rename("CF_mn")
                                , CFLIQ_mn.rename("CFLIQ_mn")
                                , CFICE_mn.rename("CFICE_mn")
                                ], compat='override')
            
            # Update attributes
            for var_name, long_name in names_dict.items():
                ds_time[var_name].attrs.update({"standard_name": var_name, "long_name": long_name})
    
                # Save to NetCDF
                output_filename = f"{directory}/{fn}_Trangeprops_"+str(Trange[0])+"_"+str(Trange[1])+"_"+str(t).zfill(3)+".nc"
                ds_time.to_netcdf(output_filename)


        ## trying to save memory
        del xr_p
        del xr_t
        del xr_rho_full
        del xr_rho
        del xr_w
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
        gc.collect()
        
