# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:34:20 2024

@author: ashle
"""
import numpy as np
from density_plotter import density_plot
import glob
import netCDF4 as nc
from datetime import datetime, timedelta 
import pickle


intervals = 2
res = 0.25
max_num_matchups = 10000   
num_variables = 5
cyg_version = 'l2_3.2' # 'noaa_1.2' OR 'l2_3.2' or 'storm_centric'
comp = 'ashle'

# one scene
SAR_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\sar full'
cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\cyg_'+cyg_version
sar_files_list = glob.glob(SAR_folder + r'\*.nc')
cyg_files_list = np.asarray(glob.glob(cyg_folder+ r"\*.nc"))

# make a list of cygnss dates and SMAP dates
cyg_file_dates = []
cyg_files_num = []
file_counter = 0
for cyg_file in cyg_files_list:
    start_ind = cyg_file.find('ddmi.s')
    cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
    cyg_file_dates.append(cyg_datetime_object.date())
    

for interval in intervals:
    #make an array to hold all data
    matchups = np.empty((max_num_matchups,2+num_variables))
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
                            sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
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
                        aggregated = np.ma.zeros((len(lat_bins), len(lon_bins),2+num_variables))

                        mask_array = np.zeros((len(lat_bins), len(lon_bins))).astype(datetime)
                        
                        # Aggregate ref data in each spatial bin
                        for i in range(len(lat_bins)):
                            for j in range(len(lon_bins)):
                                indices_in_bin_sar = np.where((sar_lat_bin_indices == i + 1) & (sar_lon_bin_indices == j + 1))[0]
                                indices_in_bin_smap = np.where((smap_lat_bin_indices == i + 1) & (smap_lon_bin_indices == j + 1))[0]
                                if len(indices_in_bin_smap) > 0 and len(indices_in_bin_sar) > 0:
                                    aggregated[i, j,1] = np.mean(cygnss_wind_AOI[indices_in_bin_smap])
                                    aggregated[i, j,0] = np.mean(sar_wind[indices_in_bin_sar])
                                    if aggregated.data[i,j,1] is np.nan or aggregated.data[i,j,1]<0.1:
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
                        
                        
                        