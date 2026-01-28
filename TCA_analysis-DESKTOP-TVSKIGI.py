# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:56:16 2024
Purpose: Compare SAR, CYGNSS, and JPL SMAP
@author: Ashley
"""

import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import glob
import h5py
import matplotlib
import matplotlib.pyplot as plt
from cartopy import crs as ccrs 
import cartopy.feature as cfeature
from shapely.geometry import Point
from shapely.prepared import prep
from scipy.stats import linregress
from density_plotter import density_plot


comp = 'ashle'
# one scene
# sar_files_list=[r'C:\Users\\'+comp+r'\\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\ALL_During\STAR_SAR_20231023183426_SH012024_01P_MERGED_FIX_3km.nc']
# cyg_files_list = [r'C:\Users\\'+comp+r'\\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\ALL_During\cyg.ddmi.s20231023-000000-e20231023-235959.l2.wind-mss.a32.d33.nc']
# SMAP_folder = r'C:\Users\\'+comp+r'\\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\SMAP_during'
SAR_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\sar full'
SMAP_folder = r'C:\Users\\'+comp+'\Documents\jpl smap full'
smap_files_list = np.asarray(glob.glob(SMAP_folder+ r"\*.h5"))
sar_files_list = glob.glob(SAR_folder + r'\*.nc')
cyg_version = 'l2_3.2' # 'noaa_1.2' OR 'l2_3.2' or 'storm_centric'
gap_smap_sar = 3 # smap and Cyg can only be x hours separate from each other
gap_mid_cyg=2 # cyg can be within x of the min/max of smap and sar. Ex, Smap is 3 hrs after SAR, so Cyg can be 3+1=4 hrs
res = 0.25 # the size of the grid
max_num_matchups = 200000
error_type = 'multi' # 'multi' or 'add'

cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\cyg_'+cyg_version
cyg_files_list = np.asarray(glob.glob(cyg_folder+ r"\*.nc"))

if cyg_version == 'storm_centric':
    cyg_files_list = [r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Data\cyg_storm_centric\cyg.ddmi.ADRIAN.ep.2023.01.l3.merge-grid-wind.archive.a32.d33.nc']
    
# make a list of cygnss dates and SMAP dates
cyg_file_dates = []
cyg_files_num = []
file_counter = 0
for cyg_file in cyg_files_list:
    if cyg_version == 'storm_centric':
        cyg_nc = nc.Dataset(cyg_file)
        date_format = "%y-%m-%d %H:%M:%S"
        cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
        cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
        cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(hours=int(x)))(cyg_nc.variables['time'][:])
        for i in cyg_datetime_array:
            cyg_file_dates.append(i.date())
            cyg_files_num.append(file_counter)
        file_counter += 1
    else: 
        start_ind = cyg_file.find('ddmi.s')
        cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
        cyg_file_dates.append(cyg_datetime_object.date())
    
    
# # have 2 options with SMAP - I could let it find the files, or I can make a netCDF of them all for it to slice into
smap_file_dates = []
for smap_file in smap_files_list:
    start_ind = smap_file.find('SSS')
    smap_datetime_object = datetime.strptime(smap_file[start_ind+10:start_ind+25], "%Y%m%dT%H%M%S")
    # start_ind = smap_file.find('ly_')
    # smap_datetime_object = datetime.strptime(smap_file[start_ind+3:start_ind+13], "%Y_%m_%d")
    smap_file_dates.append(smap_datetime_object.date())
smap_dates_unique = np.unique(smap_file_dates)    

#make an array to hold all data
matchups = np.empty((max_num_matchups,6))
counter = 0
num_events = 0

# go by SAR image
for SAR_file in sar_files_list:
    smap_AOI = False # to record is
    cyg_AOI = False
    
# SAR_file = sar_files_list[2]
#open the sar nc file and have the date
    sar_nc = nc.Dataset(SAR_file)
    sar_start_time = datetime(2000,1,1)
    seconds_since_start = float(sar_nc.variables['acquisition_time'][:])
    sar_time = datetime(2000,1,1) + timedelta(seconds=seconds_since_start)
    sar_date = sar_time.date()
    try:
        assert sar_date in smap_dates_unique and sar_date in cyg_file_dates
        
        sar_wind_w_land = sar_nc.variables['sar_wind'][:] ; sar_lat = sar_nc.variables['latitude'][:] ; sar_lon = sar_nc.variables['longitude'][:]
        
        # masking land surfaces
        if True:
            land_feature = cfeature.NaturalEarthFeature(
                'physical', 'land', '10m', edgecolor='black', facecolor='none')
            land_geoms = list(land_feature.geometries())
            prepared_land_geoms = [prep(geom) for geom in land_geoms]
            
            land_mask = np.zeros(sar_lon.shape, dtype=bool)
            
            # Iterate over each point in the grid
            for i in range(sar_lon.shape[0]):
                for j in range(sar_lon.shape[1]):
                    point = Point(sar_lon[i, j], sar_lat[i, j])
                    # Check if the point is within any of the land geometries
                    if any(geom.contains(point) for geom in prepared_land_geoms):
                        land_mask[i, j] = True
            combined_mask = np.logical_or(sar_wind_w_land.mask,land_mask)
            sar_wind = np.ma.masked_where(combined_mask, sar_wind_w_land)
        else: 
            sar_wind = sar_wind_w_land
        # define the AOI based on extent of SAR - assumed to be a simple box - maybe more complex options for if its star shaped
        AOI_lat_min,AOI_lat_max = np.min(sar_lat),np.max(sar_lat)
        AOI_lon_min,AOI_lon_max = np.min(sar_lon),np.max(sar_lon)
        
        
        SMAP_AOI_time_start,SMAP_AOI_time_end = sar_time-timedelta(hours=gap_smap_sar),sar_time+timedelta(hours=gap_smap_sar)
        
        #find and open the corresponding smap file
        if sar_time.hour > 24-gap_smap_sar: # if its possible there was similar measurement on day after
            smap_date_inds= np.where([sar_date == x or sar_date+timedelta(days =1) == x for x in smap_file_dates])[0]
        elif sar_time.hour < 0 + gap_smap_sar: # if its possible there was similar measurement on day before
            smap_date_inds= np.where([sar_date == x or sar_date-timedelta(days =1) == x for x in smap_file_dates])[0]
        else:
            smap_date_inds = np.where([sar_date == x for x in smap_file_dates])[0]
    
            
        # if there's multiple dates of measurements - maybe can update this to be like Cyg
        jpl_lats=[]; jpl_lons=[]; jpl_winds=[]; jpl_times=[]
        files_list = smap_files_list[smap_date_inds]      
        for files in files_list:
            file = h5py.File(files, 'r')
            ascat_wind = file["/smap_high_spd"][:] ; ascat_lats = file["/lat"][:] ; ascat_lons = file['lon'][:]; ascat_time = file['row_time'][:].astype(np.float64) ; ascat_uncertainty = file['smap_ambiguity_spd'][:]
            jpl_lats.append(ascat_lats.flatten()); jpl_lons.append(ascat_lons.flatten()); jpl_winds.append(ascat_wind.flatten()); jpl_times.append(ascat_time.flatten())
                # start_date = datetime(2015, 1, 1, 0, 0, 0, 0)
                # datetime_array_1D = np.vectorize(lambda x: start_date + timedelta(seconds=x))(ascat_time.data)
        
        jpl_lats_merged=np.concatenate(jpl_lats).flatten();jpl_lons_merged=np.concatenate(jpl_lons).flatten();jpl_winds_merged=np.concatenate(jpl_winds).flatten(); jpl_times_merged = np.repeat(np.concatenate(jpl_times), 76).flatten()
        jpl_lats_filtered = jpl_lats_merged[jpl_winds_merged>0]
        jpl_lons_filtered = jpl_lons_merged[jpl_winds_merged>0]
        jpl_winds_filtered = jpl_winds_merged[jpl_winds_merged>0]
        jpl_times_filtered = jpl_times_merged[jpl_winds_merged>0]
        
        
        #convert time to datetime format
        lats_boolean = (jpl_lats_filtered > AOI_lat_min) & (jpl_lats_filtered < AOI_lat_max)
        lons_boolean = (jpl_lons_filtered > AOI_lon_min) & (jpl_lons_filtered < AOI_lon_max)
        box_indices = np.where(lats_boolean & lons_boolean)[0]
        if len(box_indices)>0:
            smap_start_date = datetime(2015, 1, 1, 0, 0, 0, 0)
            smap_datetime_array = np.vectorize(lambda x: smap_start_date + timedelta(seconds=x))(jpl_times_filtered[box_indices])
            # print(np.where(~datetime_array.mask))
            time_boolean = (smap_datetime_array > SMAP_AOI_time_start) & (smap_datetime_array < SMAP_AOI_time_end) # error here
            true_indices = np.where(time_boolean)[0]
            if len(true_indices)>0:
                smap_wind_AOI = jpl_winds_filtered[box_indices[true_indices]] ;  smap_time_AOI = smap_datetime_array[true_indices] 
                smap_lats_AOI = jpl_lats_filtered[box_indices[true_indices]] ; smap_lons_AOI = jpl_lons_filtered[box_indices[true_indices]]
                smap_AOI = True
        
        
        if smap_AOI:
            mid_time = datetime.fromtimestamp((sar_time.timestamp() + smap_time_AOI[0].timestamp())/2) # convert to timestamps for average, then back to datetime
            cyg_AOI_time_start,cyg_AOI_time_end = mid_time-timedelta(hours=gap_mid_cyg),mid_time+timedelta(hours=gap_mid_cyg)
    
            #find and open the corresponding cygnss file - maybe need ot change this part to depend on SMAP 
            # if mid_time.hour > 24-gap_mid_cyg: # if its possible there was similar measurement on day after
            #     date_list = [mid_time.date,mid_time.date + timedelta(days=1)]
            
            #     # cyg_date_inds= np.where([sar_date == x or sar_date+timedelta(days =1) == x for x in cyg_file_dates])[0]
            # elif mid_time.hour < 0 + gap_mid_cyg: # if its possible there was similar measurement on day before
            #     date_list = [mid_time.date,mid_time.date - timedelta(days=1)]    
            #     # cyg_date_inds= np.where([sar_date == x or sar_date-timedelta(days =1) == x for x in cyg_file_dates])[0]
            # else:
            #     date_list = [mid_time.date]
            #     # cyg_date_inds = np.where([sar_date == x for x in cyg_file_dates])[0][0]
            # cyg_date_list = np.array(cyg_files_list[cyg_date_inds]) # may need to remove brackets here
            
            date_list = [cyg_AOI_time_start.date(),cyg_AOI_time_end.date()]
            if date_list[0]==date_list[1]:
                date_list.pop()
                
            cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_winds=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([])
            for date in date_list:
                assert date in cyg_file_dates
                if cyg_version == 'storm_centric':
                    cyg_file_path = cyg_files_list[cyg_files_num[np.where([date == x for x in cyg_file_dates])[0][0]]]
                    cyg_nc = nc.Dataset(cyg_file_path)
                    cygnss_lons = cyg_nc.variables['lon'][:] - 180 ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                    cygnss_time =cyg_nc.variables['time'][:] ; cygnss_time_offset = cyg_nc.variables['time_offset'][:]
                    lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                    lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                    box_indices = [0,0]
                    
                else:
                    cyg_file_path = cyg_files_list[np.where([date == x for x in cyg_file_dates])[0][0]]
                    cyg_nc = nc.Dataset(cyg_file_path)
                    cygnss_lons = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                    cygnss_qual = cyg_nc.variables['sample_flags']
                    cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                    lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                    lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                    box_indices = np.where(lats_boolean & lons_boolean)[0]
                    
                if cyg_version== 'l2_3.2':
                    cygnss_wind= cyg_nc.variables['preliminary_yslf_wind_speed'][:] ; 
                    inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['preliminary_yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
                    cygnss_qual = cyg_nc.variables['fds_sample_flags']
                #convert time to datetime format
                

                if len(box_indices)>0:
                    if cyg_version == 'storm_centric':
                        # find closest epoch to SAR time
                        date_format = "%y-%m-%d %H:%M:%S"
                        cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
                        cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                        cyg_measured_times = np.ma.MaskedArray(np.zeros(np.shape(cygnss_time_offset)))
                        for epoch in range(len(cygnss_time)):
                            cyg_measured_times[epoch] = cygnss_time_offset[epoch] + cygnss_time[epoch]
                        # Focus on the lat and lons
                        timedeltas = np.vectorize(lambda x: timedelta(hours=x))(cyg_measured_times[:,lats_boolean,:][:,:,lons_boolean])
                        # Add the base datetime to each timedelta to get the 3D array of datetime objects
                        cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + x)(timedeltas)                        
                        
                    else:
                        cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                        date_format = "%Y-%m-%d %H:%M:%S"
                        cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                        cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
                    
                    time_boolean = (cyg_datetime_array > cyg_AOI_time_start) & (cyg_datetime_array < cyg_AOI_time_end)
                    true_indices = np.where(time_boolean)[0]
                    if len(true_indices)>0:
                        if cyg_version == 'storm_centric':
                            cygnss_wind_AOI = cygnss_wind[:,lats_boolean,:][:,:,lons_boolean][time_boolean] ; cygnss_time_AOI = cyg_datetime_array[time_boolean]
                            cygnss_lats_AOI = cygnss_lats[lats_boolean][np.where(time_boolean)[1]] ; cygnss_lons_AOI = cygnss_lons[lons_boolean][np.where(time_boolean)[2]]
                        else:
                            cygnss_wind_AOI = cygnss_wind[box_indices[true_indices]] ; cygnss_time_AOI = cyg_datetime_array[true_indices] ; 
                            #   sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
                            cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]
                        cyg_AOI = True 
                        
                        cyg_lats=np.append(cyg_lats,cygnss_lats_AOI) ; cyg_lons = np.append(cyg_lons,cygnss_lons_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI); cyg_times=np.append(cyg_times,cygnss_time_AOI)
                        
                        
                    
        
                        
                    # true_lats = [] ; true_lons = []
                    # # Iterate through the 2D array and find indexes within the range
                    # for i, row in enumerate(datetime_array):
                    #     for j, dt in enumerate(row):
                    #         if type(dt) is datetime:
                    #             if AOI_time_start <= dt <= AOI_time_end:
                    #                 true_lats.append(i) ; true_lons.append(j)
                    # # true_indices = np.where(time_boolean)[0]
                    # if len(true_lons)>0:
                    #     smap_wind_AOI = smap_wind_clipped[true_lats,true_lons] ; smap_time_AOI = datetime_array[true_lats,true_lons]
                    #     smap_lats_AOI = smap_lats[lats_boolean][true_lats] ; smap_lons_AOI = smap_lons[lons_boolean][true_lons]
                    #     smap_AOI = True 
            
            
        
            
            if cyg_AOI and smap_AOI:
                # comparing the two datasets
                # analysis of .25 deg collocations
                # find how many boxes there are in AOI
                AOI_lat_grid = int((AOI_lat_max - AOI_lat_min)/res)
                AOI_lon_grid = int((AOI_lon_max - AOI_lon_min)/res)
                AOI_grid = np.empty((AOI_lat_grid,AOI_lon_grid))
                
                # Define spatial bins
                lat_bins = np.arange(AOI_lat_min, AOI_lat_max, res)
                lon_bins = np.arange(AOI_lon_min, AOI_lon_max, res)
                x_grid, y_grid = np.add(np.meshgrid(lon_bins, lat_bins),0.125)
                x_grid = np.ma.masked_array(x_grid)
                y_grid = np.ma.masked_array(y_grid)
                # Use digitize to assign each point to a bin
                sar_lat_bin_indices = np.digitize(sar_lat, lat_bins)
                sar_lon_bin_indices = np.digitize(sar_lon, lon_bins)
                cyg_lat_bin_indices = np.digitize(cyg_lats, lat_bins)
                cyg_lon_bin_indices = np.digitize(cyg_lons, lon_bins)
                smap_lat_bin_indices = np.digitize(smap_lats_AOI, lat_bins)
                smap_lon_bin_indices = np.digitize(smap_lons_AOI, lon_bins)
                
                
                # Create an empty array to store aggregated data, and a common mask
                aggregated_ref_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
                aggregated_cygnss_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
                aggregated_smap_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
        
                aggregated_ref_time = np.ma.MaskedArray(np.full((len(lat_bins), len(lon_bins)), datetime(2000, 1, 1, 12, 0), dtype=object))
                aggregated_cygnss_time = np.ma.MaskedArray(np.full((len(lat_bins), len(lon_bins)), datetime(2000, 1, 1, 12, 0), dtype=object))
                aggregated_smap_time = np.ma.MaskedArray(np.full((len(lat_bins), len(lon_bins)), datetime(2000, 1, 1, 12, 0), dtype=object))
        
                mask_array = np.zeros((len(lat_bins), len(lon_bins))).astype(datetime)
                
                
                # Aggregate ref data in each spatial bin
                for i in range(len(lat_bins)):
                    for j in range(len(lon_bins)):
                        indices_in_bin = np.where((sar_lat_bin_indices == i + 1) & (sar_lon_bin_indices == j + 1))[0]
                        if len(indices_in_bin) > 0:
                            aggregated_ref_data[i, j] = np.mean(sar_wind[indices_in_bin])
                            aggregated_ref_time[i, j] = sar_time
                        elif (len(indices_in_bin) == 0) :
                            mask_array[i, j] = 1
                        if aggregated_ref_data[i, j] is np.nan:
                            mask_array[i, j] = 1
                    
                
                # Aggregate cygnss data in each spatial bin
                for i in range(len(lat_bins)):
                    for j in range(len(lon_bins)):
                        indices_in_bin = np.where((cyg_lat_bin_indices == i + 1) & (cyg_lon_bin_indices == j + 1))[0]
                        if len(indices_in_bin) > 0:
                            aggregated_cygnss_data[i, j] = np.mean(cyg_winds[indices_in_bin])
                            cyg_timestamps = np.array([dt.timestamp() for dt in cyg_times[indices_in_bin]])
                            aggregated_cygnss_time[i,j] = np.nanmean(cyg_timestamps)
                            # aggregated_cygnss_time[i, j] = datetime.fromtimestamp(cyg_mean_timestamp)
                            if aggregated_cygnss_data.data[i, j] is np.nan or aggregated_cygnss_data.data[i, j]<0.1:
                                mask_array[i, j] = 1
                        else:
                            mask_array[i, j] = 1
                
                # Aggregate smap data in each spatial bin
                for i in range(len(lat_bins)):
                    for j in range(len(lon_bins)):
                        indices_in_bin = np.where((smap_lat_bin_indices == i + 1) & (smap_lon_bin_indices == j + 1))[0]
                        if len(indices_in_bin) > 0:
                            aggregated_smap_data[i, j] = np.mean(smap_wind_AOI[indices_in_bin])
                            # smap_timestamps = np.array([dt.timestamp() for dt in smap_time_AOI[indices_in_bin]])
                            # smap_mean_timestamp = np.nanmean(smap_timestamps)
                            # aggregated_smap_time[i, j] = datetime.fromtimestamp(smap_mean_timestampmp)
                            if aggregated_smap_data.data[i, j] is np.nan or aggregated_smap_data.data[i, j]<0.1:
                                    mask_array[i, j] = 1
                        else:
                            mask_array[i, j] = 1
                # convert masked array to boolean to mask the bad values
                bool_arr = mask_array.astype(bool)
                aggregated_ref_data.mask = bool_arr
                aggregated_cygnss_data.mask = bool_arr
                aggregated_cygnss_time.mask = bool_arr 
                aggregated_smap_data.mask = bool_arr
                y_grid.mask = bool_arr
                x_grid.mask = bool_arr
                
                # compressing the arrays to be only relevant values 
                ref_compressed = aggregated_ref_data.compressed()
                cyg_compressed = aggregated_cygnss_data.compressed()
                cyg_time_compressed = aggregated_cygnss_time.compressed()
                smap_compressed = aggregated_smap_data.compressed()
                lat_compressed = y_grid.compressed()
                lon_compressed = x_grid.compressed()
                
                # packaging it all up inot a list 
                values = [cyg_compressed,ref_compressed,smap_compressed,lat_compressed,lon_compressed,cyg_time_compressed, ]
                for i in range(len(values)):
                    matchups[counter:counter+len(ref_compressed),i] = values[i]
                counter += len(ref_compressed)
                num_events+=1
                # make an array of the differences 
                # diff = aggregated_ref_data - aggregated_cygnss_data
                # diff.mask = np.where(aggregated_cygnss_data < 0.1)
         
        
                plotting = True
                # if smap_AOI and cyg_AOI:
                if plotting:
                   
                    # Create a Cartopy PlateCarree projection (cylindrical projection)
                    projection = ccrs.PlateCarree()
                    
                    # Create a Matplotlib figure and axis
                    fig, ax = plt.subplots(subplot_kw={'projection': projection})
                    vmin = 0 ; vmax = 25
                    im = ax.pcolormesh(sar_lon, sar_lat, sar_wind, cmap='BuGn',transform=projection, vmin=vmin, vmax=vmax)
                    
                    im = ax.scatter(cygnss_lons_AOI, cygnss_lats_AOI, c=cygnss_wind_AOI , cmap='BuGn',transform=projection, marker='^',edgecolors = 'black',vmin=vmin, vmax=vmax, s=30)            
                    im = ax.scatter(smap_lons_AOI, smap_lats_AOI, c=smap_wind_AOI , cmap='BuGn',transform=projection, marker='s',edgecolors = 'black',vmin=vmin, vmax=vmax, s=30)            
                    # im = ax.scatter(cygnss_lons_AOI, cygnss_lats_AOI, c=cygnss_YSLF_AOI, cmap='viridis',transform=projection, vmin=vmin, vmax=vmax, marker='1')
                    cbar = plt.colorbar(im, ax=ax, orientation='vertical',location='left', shrink=0.8)
                    cbar.set_label('wind speed (m/s)')
                        
                    # if len(diff)>0:
                    #     im = ax.pcolormesh(x_grid, y_grid, diff, cmap='Reds_r', alpha=0.8)
                    #     cbar = plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.8)
                    #     cbar.set_label('wind speed diff (m/s)')
                
                    
                    # Add coastlines and gridlines for better context
                    ax.coastlines()
                    buffer = 0.07
                    ax.set_xlim(AOI_lon_min-buffer,AOI_lon_max+buffer)  # Set x-axis limits
                    ax.set_ylim(AOI_lat_min-buffer,AOI_lat_max+buffer)  # Set y-axis limits
                    plt.xlabel('Longitude')
                    plt.ylabel('Latitude')
                    gls = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
                    gls.xlocator = plt.FixedLocator(range(-180, 181, 5))  # Longitude lines every 10 degrees
                    gls.ylocator = plt.FixedLocator(range(-90, 91, 5))   # Latitude lines every 10 degrees
                    gls.top_labels=False   # suppress top labels
                    gls.right_labels=False # suppress right labels
                    ax.set_title(f'SMAP, SAR, & CYGNSS at {str(cyg_AOI_time_start)[11:16]}-{str(cyg_AOI_time_end)[11:16]}', loc='right')
                    
                    plt.subplots_adjust(left=0.55, right=0.65, bottom=0.1, top=0.9, wspace=0.5)
                    plt.tight_layout()
                    
                    # # saving file
                    directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\TCA{1}.png'.format(comp,num_events)
                    # fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
                    plt.show()
    except AssertionError:
        pass
            
    
final_matchups = matchups[:counter]

# Save the list to a file
import pickle
with open(SAR_folder+'\matchupsTCA_'+cyg_version+'.pkl', 'wb') as f:
    # filtered_array = np.array([[[x, y] for x, y in row if not np.isnan(x) and not np.isnan(y)] for row in matchups])
    pickle.dump(final_matchups, f, protocol=pickle.HIGHEST_PROTOCOL)




# Open the file and load the list
import pickle
with open(SAR_folder+'\matchupsTCA_'+cyg_version+'.pkl', 'rb') as f:
    filtered_array = pickle.load(f)

def tc_analysis(matchups,error_type): 
    # input: a array which had 3 or more 1d arrays, with a lat and lon 1d array    
    # output: RMSE, error_sd, and correlation to unknown
    Means = np.mean(matchups,axis=0)

    # all to be 2D arrays of values that have been aggregated to the grid
    # ensure that the pixels are already in Log form    
    if error_type == 'multi':
        matchups = np.log(matchups)
    
    cov_matrix=np.cov(matchups,rowvar=False)
    
    error_sd=np.empty(3)
    error_sd[0]=np.sqrt(cov_matrix[0,0]-(cov_matrix[0,1]*cov_matrix[0,2]/cov_matrix[1,2]))
    error_sd[1]=np.sqrt(cov_matrix[1,1]-(cov_matrix[0,1]*cov_matrix[1,2]/cov_matrix[0,2]))
    error_sd[2]=np.sqrt(cov_matrix[2,2]-(cov_matrix[0,2]*cov_matrix[1,2]/cov_matrix[0,1]))
    
    correlation = np.empty(3)
    correlation[0]=np.sqrt((cov_matrix[0,1]*cov_matrix[0,2])/(cov_matrix[0,0]*cov_matrix[1,2]))
    correlation[1]=np.sign(cov_matrix[0,2]*cov_matrix[1,2])* np.sqrt((cov_matrix[0,1]*cov_matrix[1,2])/(cov_matrix[1,1]*cov_matrix[0,2]))
    correlation[2]=np.sign(cov_matrix[0,1]*cov_matrix[1,2])* np.sqrt((cov_matrix[0,2]*cov_matrix[1,2])/(cov_matrix[2,2]*cov_matrix[0,1]))

    RMSE=np.empty(3)
    RMSE=error_sd*Means
                
    return error_sd, RMSE, correlation



# tc_analysis(filtered_array[:,0:3])
nan_mask = (filtered_array == 0) 

# Use np.any to find rows with any NaN values
rows_with_nan = np.any(nan_mask, axis=1)


# Use boolean indexing to select rows without NaN values
cleaned_array = filtered_array[~rows_with_nan]
matchups = cleaned_array[:,0:3]
ubRMSE, RMSE, cross_corr = tc_analysis(matchups,error_type)


# boot strapping

index_array = np.arange(len(matchups))

# bootstrapping
# Number of bootstrap samples
n_bootstraps = 1000

# Array to store bootstrap results
bootstrap_RMSE = np.zeros((n_bootstraps,3))
bootstrap_ubRMSE = np.zeros((n_bootstraps,3))
bootstrap_corr = np.zeros((n_bootstraps,3))

# Perform bootstrapping
for i in range(n_bootstraps):
    # Generate a bootstrap sample
    bootstrap_sample = np.random.choice(index_array, size=int(len(index_array)), replace=True)
    bootstrapped_matchups = matchups[bootstrap_sample]
    if np.isinf(bootstrapped_matchups).any():
        print("Data contains inf or -inf values.")
        
    # Calculate the function on the bootstrap sample
    bootstrap_ubRMSE[i], bootstrap_RMSE[i], bootstrap_corr[i] = tc_analysis(bootstrapped_matchups,error_type)
    
# Calculate the 99% confidence interval
alpha = 0.01
RMSE_CI = np.zeros((3,3))
ubRMSE_CI = np.zeros((3,3))
corr_CI = np.zeros((3,3))
for i in range(3):
    RMSE_CI[i,0] = np.percentile(bootstrap_RMSE[:,i], 100 * alpha / 2) # lower bound
    RMSE_CI[i,1] = np.percentile(bootstrap_RMSE[:,i], 100 * (1 - alpha / 2)) # upper bound
    RMSE_CI[i,2] = np.max((RMSE_CI[i,1] - RMSE[i],RMSE[i] - RMSE_CI[i,0])) # largest difference
    ubRMSE_CI[i,0] = np.percentile(bootstrap_ubRMSE[:,i], 100 * alpha / 2)
    ubRMSE_CI[i,1] = np.percentile(bootstrap_ubRMSE[:,i], 100 * (1 - alpha / 2))
    ubRMSE_CI[i,2] = np.max((ubRMSE_CI[i,1] - ubRMSE[i],ubRMSE[i] - ubRMSE_CI[i,0]))
    corr_CI[i,0] = np.percentile(bootstrap_corr[:,i], 100 * alpha / 2)
    corr_CI[i,1] = np.percentile(bootstrap_corr[:,i], 100 * (1 - alpha / 2))
    corr_CI[i,2] = np.max((corr_CI[i,1] - cross_corr[i],cross_corr[i] - corr_CI[i,0]))

cyg_vals = np.asarray([item[0] for item in matchups])
sar_vals = np.asarray([item[1] for item in matchups])
smap_vals = np.asarray([item[2] for item in matchups])

save_name = '\TCA_'+ cyg_version
save_dir = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots' 
save_location = save_dir + save_name
density_plot(sar_vals, 'SAR', cyg_vals, 'CYGNSS', save_location +'SAR_CYG.png')
density_plot(sar_vals, 'SAR', smap_vals, 'SMAP', save_location +'SAR_SMAP.png')
density_plot(smap_vals, 'SMAP',cyg_vals, 'CYGNSS', save_location +'SMAP_CYG.png')
