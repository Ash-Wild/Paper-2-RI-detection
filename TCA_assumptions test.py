# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 12:52:51 2024

@author: Ashley
"""

# Open the file and load the list
import pickle
import matplotlib.pyplot as plt
import cfgrib
# import adfuller
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import xarray as xr
import glob
import h5py
import cartopy.feature as cfeature
from shapely.geometry import Point
from shapely.prepared import prep
from density_plotter import density_plot


comp = 'Ashley'
file_dir = r'C:\Users\\'+comp+'\Documents\Jan_2023'
SAR_folder = file_dir +'\SAR'
SMAP_folder = file_dir +'\SMAP'
ERA5_folder = file_dir + r'\ERA5'
era5_files_list = glob.glob(ERA5_folder + r'\*.grib')
smap_files_list = np.asarray(glob.glob(SMAP_folder+ r"\*.h5"))
sar_files_list = glob.glob(SAR_folder + r'\*.nc')
cyg_version = r'\noaa_1.2' # 'noaa_1.2' OR 'storm_centric'
cyg_folder = file_dir+cyg_version
cyg_files_list = np.asarray(glob.glob(cyg_folder+ r"\*.nc"))

# make a list of cygnss dates and SMAP dates
cyg_file_dates = []
cyg_storm_names = []
cyg_files_num = []
file_counter = 0
for cyg_file in cyg_files_list:
    start_ind = cyg_file.find('ddmi.')
    if cyg_version == r'\storm_centric':
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
            cyg_storm_names.append(storm_name)
        file_counter += 1
    else:
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


# era5
counter = 0
for month in ['08','09','10']:
    u_wind_grib = cfgrib.open_dataset(ERA5_folder+r'\u_{0}.grib'.format(month))
    v_wind_grib = cfgrib.open_dataset(ERA5_folder+r'\v_{0}.grib'.format(month))
    
    if month == '08':
        era5_lat = v_wind_grib.latitude.values ; era5_lon  = v_wind_grib.longitude.values
        era5_winds=np.ma.MaskedArray(np.zeros([3000,len(era5_lat),len(era5_lon)])); era5_times=np.ma.MaskedArray(np.zeros(3000).astype('datetime64[h]'))
    
    u_wind= abs(u_wind_grib['u10'].values) ; v_wind = abs(v_wind_grib['v10'].values)
    era5_wind = np.hypot(u_wind, v_wind)
    era5_time = v_wind_grib.time.values.astype('datetime64[h]')
    era5_winds[counter:counter+len(era5_wind)] = era5_wind; era5_times[counter:counter+len(era5_wind)] = era5_time
    counter += len(era5_wind)
era5_winds = era5_winds[:counter] ; era5_times = era5_times[:counter]
with open(ERA5_folder+r'\ERA5_wind.pkl','wb') as f:
    pickle.dump(era5_winds, f, protocol=pickle.HIGHEST_PROTOCOL)
with open(ERA5_folder+r'\ERA5_times.pkl','wb') as f:
    pickle.dump(era5_times, f, protocol=pickle.HIGHEST_PROTOCOL)
with open(ERA5_folder+r'\ERA5_lats.pkl','wb') as f:
    pickle.dump(era5_lat, f, protocol=pickle.HIGHEST_PROTOCOL)
with open(ERA5_folder+r'\ERA5_lons.pkl','wb') as f:
    pickle.dump(era5_lon, f, protocol=pickle.HIGHEST_PROTOCOL)

    
with open(ERA5_folder+r'\ERA5_wind.pkl', 'rb') as f:
    era5_winds = pickle.load(f)
with open(ERA5_folder+r'\ERA5_times.pkl', 'rb') as f:
    era5_times = pickle.load(f)
with open(ERA5_folder+r'\ERA5_lons.pkl', 'rb') as f:
    longitudes = pickle.load(f)
    era5_lon = ((longitudes + 180) % 360) - 180
with open(ERA5_folder+r'\ERA5_lats.pkl', 'rb') as f:
    era5_lat = pickle.load(f)

#make an array to hold all data
matchups = np.empty((10000,7))
gap_smap_sar = 3 # smap and Cyg can only be x hours separate from each other
gap_mid_cyg=2 # cyg can be within x of the min/max of smap and sar. Ex, Smap is 3 hrs after SAR, so Cyg can be 3+1=4 hrs


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
            date_list = [cyg_AOI_time_start.date(),cyg_AOI_time_end.date()]
            if date_list[0]==date_list[1]:
                date_list.pop()
                
            cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_winds=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([])
            for date in date_list:
                assert date in cyg_file_dates
                if cyg_version == r'\storm_centric':
                    cyg_file_path = cyg_files_list[cyg_files_num[np.where([date == x for x in cyg_file_dates])[0][0]]]
                    cyg_nc = nc.Dataset(cyg_file_path)
                    cyg_storm_lat = cyg_nc.variables['best_track_storm_center_lat'][:] ; cyg_storm_lon = ((cyg_nc.variables['best_track_storm_center_lon'][:] + 180) % 360) - 180
                    
                    assert True in (cyg_storm_lon > AOI_lon_min) & True in (cyg_storm_lon < AOI_lon_max)
                    assert True in (cyg_storm_lat > AOI_lat_min) & True in (cyg_storm_lat < AOI_lat_max)

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
                    
                # if cyg_version== 'l2_3.2':
                #     cygnss_yslf_wind= cyg_nc.variables['preliminary_yslf_wind_speed'][:];
                #     merge_threshold = 17
                #     cygnss_wind = np.where(cygnss_wind < merge_threshold, cygnss_wind, cygnss_yslf_wind)
                #     inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['preliminary_yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
                #     # cygnss_qual = cyg_nc.variables['fds_sample_flags']
                #convert time to datetime format
                
                # assert date in cyg_file_dates
                # cyg_file_path = cyg_files_list[np.where([date == x for x in cyg_file_dates])[0][0]]
                # cyg_nc = nc.Dataset(cyg_file_path)
                # cygnss_lons = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                # cygnss_qual = cyg_nc.variables['sample_flags']
                # cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                # lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                # lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                # box_indices = np.where(lats_boolean & lons_boolean)[0]
                # if cyg_version== 'l2_3.2':
                #     cygnss_wind= cyg_nc.variables['preliminary_yslf_wind_speed'][:] ; 
                #     inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['preliminary_yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
                #     cygnss_qual = cyg_nc.variables['fds_sample_flags']

                if len(box_indices)>0:
                    if cyg_version == r'\storm_centric':
                        # find closest epoch to SAR time
                        # print('found storm')
                        date_format = "%y-%m-%d %H:%M:%S"
                        cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
                        cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                        cyg_measured_times = np.ma.MaskedArray(np.zeros(np.shape(cygnss_time_offset)))
                        for epoch in range(len(cygnss_time)):
                            cyg_measured_times[epoch] = cygnss_time_offset[epoch] + cygnss_time[epoch]
                        # Focus on the lat and lons
                        try:
                            timedeltas = np.vectorize(lambda x: timedelta(hours=x))(cyg_measured_times[:,lats_boolean,:][:,:,lons_boolean])
                        except ValueError:
                            print('error in time')
                            raise AssertionError
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
                        if cyg_version == r'\storm_centric':
                            cygnss_wind_AOI = cygnss_wind[:,lats_boolean,:][:,:,lons_boolean][time_boolean] ; cygnss_time_AOI = cyg_datetime_array[time_boolean]
                            cygnss_lats_AOI = cygnss_lats[lats_boolean][np.where(time_boolean)[1]] ; cygnss_lons_AOI = cygnss_lons[lons_boolean][np.where(time_boolean)[2]]
                        else:
                            cygnss_wind_AOI = cygnss_wind[box_indices[true_indices]] ; cygnss_time_AOI = cyg_datetime_array[true_indices] ; 
                            #   sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
                            cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]
                            cyg_qual_AOI = cygnss_qual[box_indices[true_indices]]
                        cyg_AOI = True 
                        
                        cyg_lats=np.append(cyg_lats,cygnss_lats_AOI) ; cyg_lons = np.append(cyg_lons,cygnss_lons_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI); cyg_times=np.append(cyg_times,cygnss_time_AOI)
                        

                
                # if len(box_indices)>0:
                #     cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                #     date_format = "%Y-%m-%d %H:%M:%S"
                #     cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                #     cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
                #     time_boolean = (cyg_datetime_array > cyg_AOI_time_start) & (cyg_datetime_array < cyg_AOI_time_end)
                #     true_indices = np.where(time_boolean)[0]
                #     if len(true_indices)>0:
                #         cygnss_wind_AOI = cygnss_wind[box_indices[true_indices]] ; cygnss_time_AOI = cyg_datetime_array[true_indices] ; 
                #         #   sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
                #         cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]
                #         cyg_AOI = True 
                #         cyg_lats=np.append(cyg_lats,cygnss_lats_AOI) ; cyg_lons = np.append(cyg_lons,cygnss_lons_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI); cyg_times=np.append(cyg_times,cygnss_time_AOI)
            
            if cyg_AOI and smap_AOI:
                # ERA5 analysis
                dt = sar_time
                if dt.minute >= 30:
                    dt += timedelta(hours=1)
                era5_AOI_time = np.datetime64(dt.replace(minute=0, second=0, microsecond=0)).astype('datetime64[h]')
                
                # find index when era5 is made to SAR time
                true_indices = np.where(era5_times==era5_AOI_time)
                assert len(true_indices) >0
                                
                # constrict era5 to lat lons
                lats_boolean = (era5_lat > AOI_lat_min) & (era5_lat < AOI_lat_max)
                lons_boolean = (era5_lon > AOI_lon_min) & (era5_lon < AOI_lon_max)
                
                era5_wind_AOI = era5_winds[true_indices][:,lats_boolean,:][:,:,lons_boolean] ;  era5_time_AOI = era5_times[true_indices] 
                era5_lats_AOI = era5_lat[lats_boolean] ; era5_lons_AOI = era5_lon[lons_boolean]

                # comparing the two datasets
                # analysis of .25 deg collocations
                
                # Define spatial bins
                lat_bins = era5_lats_AOI
                lon_bins = era5_lons_AOI
                x_grid, y_grid = np.meshgrid(lon_bins, lat_bins)
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
                era5_wind_AOI.mask = bool_arr
                y_grid.mask = bool_arr
                x_grid.mask = bool_arr
                
                # compressing the arrays to be only relevant values 
                ref_compressed = aggregated_ref_data.compressed()
                cyg_compressed = aggregated_cygnss_data.compressed()
                cyg_time_compressed = aggregated_cygnss_time.compressed()
                smap_compressed = aggregated_smap_data.compressed()
                era5_compressed = era5_wind_AOI.compressed()
                lat_compressed = y_grid.compressed()
                lon_compressed = x_grid.compressed()
                
                # packaging it all up inot a list 
                values = [cyg_compressed,ref_compressed,smap_compressed,era5_compressed, lat_compressed,lon_compressed,cyg_time_compressed ]
                for i in range(len(values)):
                    matchups[counter:counter+len(ref_compressed),i] = values[i]
                counter += len(ref_compressed)
                num_events+=1
    except AssertionError:
        pass
    
final_matchups = matchups[:counter]


# tc_analysis(filtered_array[:,0:3])
nan_mask = (final_matchups == 0) 

# Use np.any to find rows with any NaN values
rows_with_nan = np.any(nan_mask, axis=1)

# Use boolean indexing to select rows without NaN values
cleaned_array = final_matchups[~rows_with_nan]
end_matchups = cleaned_array[:,0:4]

with open(r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TCA_assumptions'+cyg_version+'.pkl', 'wb') as f:
    pickle.dump(end_matchups, f, protocol=pickle.HIGHEST_PROTOCOL)





# Open the file and load the list
import pickle
with open(r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TCA_assumptions'+cyg_version+'.pkl', 'rb') as f:
    matchups = pickle.load(f)
# # # Checking assumptions

Means= np.mean(matchups,axis=0)
anom_array = matchups -  Means

# #1. Linearity
cyg_vals = np.asarray([item[0] for item in matchups])
sar_vals = np.asarray([item[1] for item in matchups])
smap_vals = np.asarray([item[2] for item in matchups])
era5_vals = np.asarray([item[3] for item in matchups])
cyg_error_vals = cyg_vals- era5_vals
sar_error_vals = sar_vals - era5_vals
smap_error_vals = smap_vals - era5_vals


save_name = 'assumptions_matchups.'+ cyg_version[1:]
save_dir = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TCA_assumptions\\'
save_location = save_dir + save_name

density_plot(era5_vals, 'ERA5', cyg_vals, 'CYGNSS' + cyg_version, save_location + 'CYG.png' ,0,51,15,False)
density_plot(era5_vals, 'ERA5', sar_vals, 'SAR', save_location + 'SAR.png',0,51,15,False)
density_plot(era5_vals, 'ERA5',smap_vals, 'SMAP', save_location + 'SMAP.png',0,51,15,False)

density_plot(sar_error_vals, 'SAR - ERA5', cyg_error_vals, 'CYGNSS - ERA5' , save_location + 'Error SAR CYG.png' ,-15,15,15,False)
density_plot(sar_error_vals, 'SAR - ERA5', smap_error_vals, 'SMAP - ERA5', save_location + 'Error SAR SMAP.png',-15,15,15,False)
density_plot(smap_error_vals, 'SMAP - ERA5',cyg_error_vals, 'CYGNSS - ERA5' , save_location + 'Error SMAP CYG.png',-15,15,15,False)

# # if True:
# #     plt.scatter(anom_array[:,0],anom_array[:,1])
# #     R = np.around(np.corrcoef(anom_array[:,0],anom_array[:,1])[0,1], decimals=2)
# #     plt.text(0.05, 0.95, f'R: {R}', transform=plt.gca().transAxes, va='top', ha='left')
# #     plt.title('CYG SAR')
# #     plt.show()
# #     plt.scatter(anom_array[:,0],anom_array[:,2])
# #     R = np.around(np.corrcoef(anom_array[:,0],anom_array[:,2])[0,1], decimals=2)
# #     plt.text(0.05, 0.95, f'R: {R}', transform=plt.gca().transAxes, va='top', ha='left')
# #     plt.title('CYG SMAP')
# #     plt.show()
# #     plt.scatter(anom_array[:,2],anom_array[:,1])
# #     R = np.around(np.corrcoef(anom_array[:,2],anom_array[:,1])[0,1], decimals=2)
# #     plt.text(0.05, 0.95, f'R: {R}', transform=plt.gca().transAxes, va='top', ha='left')
# #     plt.title('SMAP SAR')
# #     plt.show()

# # 2. Stationarity
# # Perform the Augmented Dickey-Fuller test - inappropriate as not time checked
# # Open the file and load the list
# adf_result_cyg = adfuller(cleaned_array) # sometimes gaps in timeseries though - not appropriate as not consistent
# # can check visually of anomalies - see if there's trend overtime. 

# # Extract and print the results
# adf_result = adf_result_cyg
# print('ADF Statistic:', adf_result[0])
# print('p-value:', adf_result[1])
# print('Critical Values:')
# for key, value in adf_result[4].items():
#     print(f'   {key}: {value}')

# # can check visually of anomalies - see if there's trend overtime. 
# plt.figure(figsize=(10, 5))
# plt.plot(anom_array[:,1], label='SAR', color='blue',alpha=0.4)
# plt.plot(anom_array[:,0] , label='CYG', color='red',alpha=0.4)
# plt.plot(anom_array[:,2] , label='SMAP', color='black',alpha=0.4)
# plt.plot([0, len(anom_array[:,1])], [0, 0], linestyle='--', color='gray')
# plt.legend()
# plt.show()

# # 3. independence of errors

