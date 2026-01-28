# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 09:58:31 2024
Matchups finder 
inputs: 2 datasets, 
output: plt of number of matchups 
@author: Ashley
"""

import numpy as np
from density_plotter import density_plot
import glob
# import h5py
import netCDF4 as nc
from datetime import datetime, timedelta 
import cartopy.feature as cfeature
from shapely.geometry import Point
from shapely.prepared import prep
import pickle
from scipy.stats import linregress


intervals = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
res = 0.25
max_num_matchups = 200000
cyg_version = 'noaa_1.2' # 'noaa_1.2' OR 'l2_3.2' or 'storm_centric'
comp = 'Ashley'

# one scene
SAR_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\sar full'
SMAP_folder = r'C:\Users\\'+comp+'\Documents\jpl smap full'
cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\cyg_'+cyg_version
smap_files_list = np.asarray(glob.glob(SMAP_folder+ r"\*.h5"))
sar_files_list = glob.glob(SAR_folder + r'\*.nc')
cyg_files_list = np.asarray(glob.glob(cyg_folder+ r"\*.nc"))



# make a list of cygnss dates and SMAP dates
cyg_file_dates = []
cyg_files_num = []
file_counter = 0
for cyg_file in cyg_files_list:
    start_ind = cyg_file.find('ddmi.')
    if cyg_version == 'storm_centric':
        end_ind = cyg_file.find('.202')
        storm_name = cyg_file[start_ind:end_ind]
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
        cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
        cyg_file_dates.append(cyg_datetime_object.date())
    
    
smap_file_dates = []
for smap_file in smap_files_list:
    start_ind = smap_file.find('SSS')
    smap_datetime_object = datetime.strptime(smap_file[start_ind+10:start_ind+25], "%Y%m%dT%H%M%S")
    # start_ind = smap_file.find('ly_')
    # smap_datetime_object = datetime.strptime(smap_file[start_ind+3:start_ind+13], "%Y_%m_%d")
    smap_file_dates.append(smap_datetime_object.date())
smap_dates_unique = np.unique(smap_file_dates)    

events_per_interval = []
meas_per_interval = []
R_per_interval = []
RMSE_per_interval = []

# for interval in intervals:
#     #make an array to hold all data
#     matchups = np.empty((max_num_matchups,2))
#     counter = 0
#     num_events = 0
    
#     # go by SAR image
#     for SAR_file in sar_files_list:
#         smap_AOI = False # to record is
#         cyg_AOI = False
        
#     # SAR_file = sar_files_list[2]
#     #open the sar nc file and have the date
#         sar_nc = nc.Dataset(SAR_file)
#         sar_start_time = datetime(2000,1,1)
#         seconds_since_start = float(sar_nc.variables['acquisition_time'][:])
#         sar_time = datetime(2000,1,1) + timedelta(seconds=seconds_since_start)
#         sar_date = sar_time.date()
#         try:
#             assert sar_date in smap_dates_unique
            
#             sar_wind_w_land = sar_nc.variables['sar_wind'][:] ; sar_lat = sar_nc.variables['latitude'][:] ; sar_lon = sar_nc.variables['longitude'][:]
            
#             # masking land surfaces
#             if True:
#                 land_feature = cfeature.NaturalEarthFeature(
#                     'physical', 'land', '10m', edgecolor='black', facecolor='none')
#                 land_geoms = list(land_feature.geometries())
#                 prepared_land_geoms = [prep(geom) for geom in land_geoms]
                
#                 land_mask = np.zeros(sar_lon.shape, dtype=bool)
                
#                 # Iterate over each point in the grid
#                 for i in range(sar_lon.shape[0]):
#                     for j in range(sar_lon.shape[1]):
#                         point = Point(sar_lon[i, j], sar_lat[i, j])
#                         # Check if the point is within any of the land geometries
#                         if any(geom.contains(point) for geom in prepared_land_geoms):
#                             land_mask[i, j] = True
#                 combined_mask = np.logical_or(sar_wind_w_land.mask,land_mask)
#                 sar_wind = np.ma.masked_where(combined_mask, sar_wind_w_land)
#             else: 
#                 sar_wind = sar_wind_w_land
#             # define the AOI based on extent of SAR - assumed to be a simple box - maybe more complex options for if its star shaped
#             AOI_lat_min,AOI_lat_max = np.min(sar_lat),np.max(sar_lat)
#             AOI_lon_min,AOI_lon_max = np.min(sar_lon),np.max(sar_lon)
            
            
#             SMAP_AOI_time_start,SMAP_AOI_time_end = sar_time-timedelta(hours=interval),sar_time+timedelta(hours=interval)
            
#             #find and open the corresponding smap file
#             if sar_time.hour > 24-interval: # if its possible there was similar measurement on day after
#                 smap_date_inds= np.where([sar_date == x or sar_date+timedelta(days =1) == x for x in smap_file_dates])[0]
#             elif sar_time.hour < 0 + interval: # if its possible there was similar measurement on day before
#                 smap_date_inds= np.where([sar_date == x or sar_date-timedelta(days =1) == x for x in smap_file_dates])[0]
#             else:
#                 smap_date_inds = np.where([sar_date == x for x in smap_file_dates])[0]
        
                
#             # if there's multiple dates of measurements - maybe can update this to be like Cyg
#             jpl_lats=[]; jpl_lons=[]; jpl_winds=[]; jpl_times=[]
#             files_list = smap_files_list[smap_date_inds]      
#             for files in files_list:
#                 with h5py.File(files, "r") as file:
#                     ascat_wind = file["/smap_high_spd"][:] ; ascat_lats = file["/lat"][:] ; ascat_lons = file['lon'][:]; ascat_time = file['row_time'][:].astype(np.float64) ; ascat_uncertainty = file['smap_ambiguity_spd'][:]
#                 jpl_lats.append(ascat_lats.flatten()); jpl_lons.append(ascat_lons.flatten()); jpl_winds.append(ascat_wind.flatten()); jpl_times.append(ascat_time.flatten())
#                         # start_date = datetime(2015, 1, 1, 0, 0, 0, 0)
#                         # datetime_array_1D = np.vectorize(lambda x: start_date + timedelta(seconds=x))(ascat_time.data)
                
#             jpl_lats_merged=np.concatenate(jpl_lats).flatten();jpl_lons_merged=np.concatenate(jpl_lons).flatten();jpl_winds_merged=np.concatenate(jpl_winds).flatten(); jpl_times_merged = np.repeat(np.concatenate(jpl_times), 76).flatten()
#             jpl_lats_filtered = jpl_lats_merged[jpl_winds_merged>0]
#             jpl_lons_filtered = jpl_lons_merged[jpl_winds_merged>0]
#             jpl_winds_filtered = jpl_winds_merged[jpl_winds_merged>0]
#             jpl_times_filtered = jpl_times_merged[jpl_winds_merged>0]
            
            
#             #convert time to datetime format
#             lats_boolean = (jpl_lats_filtered > AOI_lat_min) & (jpl_lats_filtered < AOI_lat_max)
#             lons_boolean = (jpl_lons_filtered > AOI_lon_min) & (jpl_lons_filtered < AOI_lon_max)
#             box_indices = np.where(lats_boolean & lons_boolean)[0]
#             if len(box_indices)>0:
#                 smap_start_date = datetime(2015, 1, 1, 0, 0, 0, 0)
#                 smap_datetime_array = np.vectorize(lambda x: smap_start_date + timedelta(seconds=x))(jpl_times_filtered[box_indices])
#                 # print(np.where(~datetime_array.mask))
#                 time_boolean = (smap_datetime_array > SMAP_AOI_time_start) & (smap_datetime_array < SMAP_AOI_time_end) # error here
#                 true_indices = np.where(time_boolean)[0]
#                 if len(true_indices)>0:
#                     smap_wind_AOI = jpl_winds_filtered[box_indices[true_indices]] ;  smap_time_AOI = smap_datetime_array[true_indices] 
#                     smap_lats_AOI = jpl_lats_filtered[box_indices[true_indices]] ; smap_lons_AOI = jpl_lons_filtered[box_indices[true_indices]]
#                     smap_AOI = True
                    
                    
#                     AOI_lat_grid = int((AOI_lat_max - AOI_lat_min)/res)
#                     AOI_lon_grid = int((AOI_lon_max - AOI_lon_min)/res)
#                     AOI_grid = np.empty((AOI_lat_grid,AOI_lon_grid))
                    
#                     # Define spatial bins
#                     lat_bins = np.arange(AOI_lat_min, AOI_lat_max, res)
#                     lon_bins = np.arange(AOI_lon_min, AOI_lon_max, res)
#                     x_grid, y_grid = np.add(np.meshgrid(lon_bins, lat_bins),0.125)
#                     x_grid = np.ma.masked_array(x_grid)
#                     y_grid = np.ma.masked_array(y_grid)
#                     # Use digitize to assign each point to a bin
#                     sar_lat_bin_indices = np.digitize(sar_lat, lat_bins)
#                     sar_lon_bin_indices = np.digitize(sar_lon, lon_bins)
#                     smap_lat_bin_indices = np.digitize(smap_lats_AOI, lat_bins)
#                     smap_lon_bin_indices = np.digitize(smap_lons_AOI, lon_bins)
                    
                    
#                     # Create an empty array to store aggregated data, and a common mask
#                     aggregated_ref_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
#                     aggregated_cygnss_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
#                     aggregated_smap_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
            
#                     mask_array = np.zeros((len(lat_bins), len(lon_bins))).astype(datetime)
                    
#                     # Aggregate ref data in each spatial bin
#                     for i in range(len(lat_bins)):
#                         for j in range(len(lon_bins)):
#                             indices_in_bin_sar = np.where((sar_lat_bin_indices == i + 1) & (sar_lon_bin_indices == j + 1))[0]
#                             indices_in_bin_smap = np.where((smap_lat_bin_indices == i + 1) & (smap_lon_bin_indices == j + 1))[0]
#                             if len(indices_in_bin_smap) > 0 and len(indices_in_bin_sar) > 0:
#                                 aggregated_smap_data[i, j] = np.mean(smap_wind_AOI[indices_in_bin_smap])
#                                 aggregated_ref_data[i, j] = np.mean(sar_wind[indices_in_bin_sar])
#                                 if aggregated_smap_data.data[i, j] is np.nan or aggregated_smap_data.data[i, j]<0.1:
#                                     mask_array[i, j] = 1
#                                 elif aggregated_ref_data[i, j] is np.nan or aggregated_ref_data[i, j]<15 :
#                                     mask_array[i, j] = 1
         
#                             else:
#                                 mask_array[i, j] = 1
                      
                            
                   
#                     # convert masked array to boolean to mask the bad values
#                     bool_arr = mask_array.astype(bool)
#                     aggregated_ref_data.mask = bool_arr
#                     aggregated_cygnss_data.mask = bool_arr
#                     aggregated_smap_data.mask = bool_arr
#                     y_grid.mask = bool_arr
#                     x_grid.mask = bool_arr
                    
#                     # compressing the arrays to be only relevant values 
#                     ref_compressed = aggregated_ref_data.compressed()
#                     cyg_compressed = aggregated_cygnss_data.compressed()
#                     smap_compressed = aggregated_smap_data.compressed()
#                     lat_compressed = y_grid.compressed()
#                     lon_compressed = x_grid.compressed()
                    
#                     # packaging it all up inot a list 
#                     values = [ref_compressed,smap_compressed]
#                     for i in range(len(values)):
#                         matchups[counter:counter+len(ref_compressed),i] = values[i]
#                     counter += len(ref_compressed)
#                     num_events+=1
                      

#         except AssertionError:
#             pass
#     final_matchups = matchups[:counter]
#     ref = final_matchups[:,0]
#     comparison = final_matchups[:,1]
#     R = np.around(np.corrcoef(ref,comparison)[0,1], decimals=2)
#     RMSE = np.around(np.sqrt(np.mean((comparison - ref) ** 2)), decimals=3)
#     R_per_interval.append(R)
#     RMSE_per_interval.append(RMSE)
#     events_per_interval.append(num_events)
#     meas_per_interval.append(len(final_matchups))
    
#     # Save the list to a file
#     with open(SAR_folder+'\matchups_SMAP_SAR_time_'+str(interval)+'.pkl', 'wb') as f:
#         # filtered_array = np.array([[[x, y] for x, y in row if not np.isnan(x) and not np.isnan(y)] for row in matchups])
#         pickle.dump(final_matchups, f, protocol=pickle.HIGHEST_PROTOCOL)
    

for interval in intervals:
    #make an array to hold all data
    matchups = np.empty((max_num_matchups,2))
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
            assert sar_date in smap_dates_unique
            
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
            
            cyg_AOI_time_start,cyg_AOI_time_end = sar_time-timedelta(hours=interval),sar_time+timedelta(hours=interval)

            date_list = [cyg_AOI_time_start.date(),cyg_AOI_time_end.date()]
            if date_list[0]==date_list[1]:
                date_list.pop()
                
            cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_winds=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([])
            for date in date_list:
                assert date in cyg_file_dates
                
                if cyg_version == 'storm_centric':
                    cyg_file_path = cyg_files_list[cyg_files_num[np.where([date == x for x in cyg_file_dates])[0][0]]]
                    cyg_nc = nc.Dataset(cyg_file_path)
                    cyg_storm_lat = cyg_nc.variables['best_track_storm_center_lat'][:] ; cyg_storm_lon = cyg_nc.variables['best_track_storm_center_lon'][:]
                    assert True in (cyg_storm_lat > AOI_lat_min) & True in (cyg_storm_lat < AOI_lat_max)
                    assert True in (cyg_storm_lon > AOI_lon_min) & True in (cyg_storm_lon < AOI_lon_max)

                    longitudes = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                    cygnss_lons = ((longitudes + 180) % 360) - 180

                    cygnss_time =cyg_nc.variables['time'][:] ; cygnss_time_offset = cyg_nc.variables['time_offset'][:]
                    lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                    lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                    box_indices = [0,0]
                    
                else:
                    cyg_file_path = cyg_files_list[np.where([date == x for x in cyg_file_dates])[0][0]]
                    cyg_nc = nc.Dataset(cyg_file_path)
                    longitudes = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                    cygnss_lons = ((longitudes + 180) % 360) - 180

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
                            # sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
                            cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]
                        cyg_AOI = True 
                        
                        
                        
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
                        smap_lat_bin_indices = np.digitize(cygnss_lats_AOI, lat_bins)
                        smap_lon_bin_indices = np.digitize(cygnss_lons_AOI, lon_bins)
                        
                        
                        # Create an empty array to store aggregated data, and a common mask
                        aggregated_ref_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
                        aggregated_smap_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
                
                        mask_array = np.zeros((len(lat_bins), len(lon_bins))).astype(datetime)
                        
                        # Aggregate ref data in each spatial bin
                        for i in range(len(lat_bins)):
                            for j in range(len(lon_bins)):
                                indices_in_bin_sar = np.where((sar_lat_bin_indices == i + 1) & (sar_lon_bin_indices == j + 1))[0]
                                indices_in_bin_smap = np.where((smap_lat_bin_indices == i + 1) & (smap_lon_bin_indices == j + 1))[0]
                                if len(indices_in_bin_smap) > 0 and len(indices_in_bin_sar) > 0:
                                    aggregated_smap_data[i, j] = np.mean(cygnss_wind_AOI[indices_in_bin_smap])
                                    aggregated_ref_data[i, j] = np.mean(sar_wind[indices_in_bin_sar])
                                    if aggregated_smap_data.data[i, j] is np.nan or aggregated_smap_data.data[i, j]<0.1:
                                        mask_array[i, j] = 1
                                    elif aggregated_ref_data[i, j] is np.nan or aggregated_ref_data[i, j]<15:
                                        mask_array[i, j] = 1
             
                                else:
                                    mask_array[i, j] = 1
                          
                        # convert masked array to boolean to mask the bad values
                        bool_arr = mask_array.astype(bool)
                        aggregated_ref_data.mask = bool_arr
                        aggregated_smap_data.mask = bool_arr
                        y_grid.mask = bool_arr
                        x_grid.mask = bool_arr
                        
                        # compressing the arrays to be only relevant values 
                        ref_compressed = aggregated_ref_data.compressed()
                        smap_compressed = aggregated_smap_data.compressed()
                        lat_compressed = y_grid.compressed()
                        lon_compressed = x_grid.compressed()
                        
                        # packaging it all up inot a list 
                        values = [ref_compressed,smap_compressed]
                        for i in range(len(values)):
                            matchups[counter:counter+len(ref_compressed),i] = values[i]
                        counter += len(ref_compressed)
                        num_events+=1
                      

        except AssertionError:
            pass
    final_matchups = matchups[:counter]
    events_per_interval.append(num_events)
    meas_per_interval.append(len(final_matchups))
        
    # Save the list to a file
    with open(SAR_folder+'\matchups_CYG'+cyg_version+'_SAR_time_'+str(interval)+'.pkl', 'wb') as f:
        # filtered_array = np.array([[[x, y] for x, y in row if not np.isnan(x) and not np.isnan(y)] for row in matchups])
        pickle.dump(final_matchups, f, protocol=pickle.HIGHEST_PROTOCOL)
        




# Open the file and load the list
import pickle
interval = 4
cyg_version = 'noaa_1.2' # 'noaa_1.2' OR 'l2_3.2' or 'storm_centric'
comparison_name = 'CYG' #'SMAP' or 'CYG'
if comparison_name =='CYG':
    comparison_name += cyg_version
with open(SAR_folder+r'\matchups_'+comparison_name+'_SAR_time_'+str(interval)+'.pkl', 'rb') as f:
    filtered_array = pickle.load(f)
    
ref = filtered_array[:,0] ; comparison = filtered_array[:,1]
slope, intercept, r_value, p_value, std_err = linregress(ref,comparison)
# R = np.around(np.corrcoef(ref,comparison)[0,1], decimals=2)
RMSE = np.around(np.sqrt(np.mean((comparison - ref) ** 2)), decimals=3)
num_meas = len(filtered_array)
print(r_value,p_value ,RMSE, num_meas)


# test linearity
Anom_ref = ref - np.mean(ref) ; Anom_comparison = comparison - np.mean(comparison) 
import matplotlib.pyplot as plt
if True:
    plt.scatter(Anom_ref,Anom_comparison)
    R = np.around(np.corrcoef(Anom_ref,Anom_comparison)[0,1], decimals=2)
    plt.text(0.05, 0.95, f'R: {R}', transform=plt.gca().transAxes, va='top', ha='left')
    plt.title(comparison_name + ' SAR')
    plt.show()
# density_plot(ref, 'SAR', comparison, comparison_name, interval)

# test stationarity
from statsmodels.tsa.stattools import adfuller
# sorted_indices = np.argsort(filtered_array[:, -1]) # sort based on time - unneisary as already it increases in time
# sorted_array = filtered_array[sorted_indices,1]

adf_result = adfuller(Anom_comparison) # sometimes gaps in timeseries though - not appropriate as not consistent
print('ADF Statistic:', adf_result[0])
print('p-value:', adf_result[1])
print('Critical Values:')
for key, value in adf_result[4].items():
    print(f'   {key}: {value}')
    
# can check visually of anomalies - see if there's trend overtime. 
plt.figure(figsize=(10, 5))
plt.plot(Anom_ref, label='SAR', color='blue',alpha=0.4)
plt.plot(Anom_comparison , label=comparison_name, color='red',alpha=0.4)
plt.plot([0, 0], [len(Anom_ref), 0], linestyle='--', color='gray')
plt.legend()
plt.show()