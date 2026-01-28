# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 10:15:16 2024
Purpose: finding the DDMs of RI area
Input: 7 flighy machines datasets
Output: DDMs over RI area
@author: Ashley
"""
import glob
import netCDF4
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import h5py
import matplotlib.pyplot as plt
from cartopy import crs as ccrs 

# datasets
comp = 'ashle'
dir = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\\'
L2_dir_during = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\L2_during'
L2_dir_before = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\L2_before'
SMAP_before = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\SMAP_before'
SMAP_during = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\SMAP_during'
L1_dir_during = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\L1_during\others'
L1_dir_before = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\L1_before'
ASCAT_during = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\ASCAT during'
SMAP_RSS_during = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\SMAP_during'
# area of RI
# AOI_lat_min,AOI_lat_max = -13.1,-12.9
# AOI_lon_min,AOI_lon_max = 168.65,168.7
AOI_lat_min,AOI_lat_max = -14,-10
AOI_lon_min,AOI_lon_max = 168,170
timing_threshold = 4
approx_ref_time = datetime(2023, 10, 22, 18, 33, 33, 71000)
res = 0.25 # the size of the grid
# AOI_time_start,AOI_time_end = datetime(2023,10,22, 0,0,0),datetime(2023,10,22, 23,59,59)
AOI_time_start,AOI_time_end = approx_ref_time-timedelta(hours=timing_threshold),approx_ref_time+timedelta(hours=timing_threshold)
# AOI_time_start,AOI_time_end = datetime(2023,10,23, 22,38,10),datetime(2023,10,23, 22,38,59)

#say that initially there are no measurements
ascat_AOI = False
cygnss_AOI = False
SMAP_AOI = False
valid = []

datasets = [SMAP_before,L2_dir_before]
for dataset in datasets:
    if dataset in (SMAP_during,SMAP_before):
        files_list = glob.glob(dataset+ "\*.h5")
        for files in files_list:
            file = h5py.File(files, 'r')
            ascat_wind = file["/smap_high_spd"][:] ; ascat_lats = file["/lat"][:] ; ascat_lons = file['lon'][:]; ascat_time = file['row_time'][:].astype(np.float64) ; ascat_uncertainty = file['smap_ambiguity_spd'][:]
            start_date = datetime(2015, 1, 1, 0, 0, 0, 0)
            datetime_array_1D = np.vectorize(lambda x: start_date + timedelta(seconds=x))(ascat_time.data)
            datetime_array = np.tile(datetime_array_1D[:, np.newaxis], 76).T
            lats_boolean = (ascat_lats > AOI_lat_min) & (ascat_lats < AOI_lat_max)
            lons_boolean = (ascat_lons > AOI_lon_min) & (ascat_lons < AOI_lon_max)
            time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end)
            wind_boolean = (ascat_wind > 0)
            AOI_boolean = lats_boolean & lons_boolean & time_boolean & wind_boolean
            true_indices = np.where(AOI_boolean)
            if len(true_indices[0])>0:
                # mask any areas in the wind and time fields that are outside AOI
                ascat_wind_AOI = ascat_wind[true_indices]
                ascat_time_AOI = datetime_array[true_indices]
                ascat_lats_AOI = ascat_lats[true_indices]
                ascat_lons_AOI = ascat_lons[true_indices]
                ascat_uncertainty_AOI = np.ma.masked_equal(ascat_uncertainty[true_indices],-9999)
                SMAP_AOI = True
                    
    files_list = glob.glob(dataset+ "\*.nc")
    for file in files_list:
        nc = netCDF4.Dataset(file)
        if 'ASCAT' in file:
            ascat_lons = nc.variables['lon'][:] ; ascat_lats = nc.variables['lat'][:]; ascat_wind= nc.variables['wind_speed'][:] ; ascat_time =nc.variables['time'][:].astype(np.float64)
            #convert time to datetime format
            start_date = datetime(1990, 1, 1, 0, 0, 0, 0)
            datetime_array = np.vectorize(lambda x: start_date + timedelta(seconds=x))(ascat_time.data)
            # masking to AOI
            lats_boolean = (ascat_lats.data > AOI_lat_min) & (ascat_lats.data < AOI_lat_max)
            lons_boolean = (ascat_lons.data > AOI_lon_min) & (ascat_lons.data < AOI_lon_max)
            time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end)
            AOI_boolean = lats_boolean & lons_boolean & time_boolean
            
            true_indices = np.where(AOI_boolean)[0]
            if len(true_indices)>0:
                # mask any areas in the wind and time fields that are outside AOI
                wind_new_mask = np.logical_or(ascat_wind.mask, ~AOI_boolean)
                ascat_wind_AOI = np.ma.masked_array(ascat_wind.data, wind_new_mask).compressed()
                ascat_time_AOI = np.ma.masked_array(datetime_array.data, wind_new_mask).compressed()
                ascat_lats_AOI = np.ma.masked_array(ascat_lats.data, wind_new_mask).compressed()
                ascat_lons_AOI = np.ma.masked_array(ascat_lons.data, wind_new_mask).compressed()
                ascat_AOI = True
                # need to add in something to add in if multiple observations
        if 'RSS' in file:
            SMAP_wind, SMAP_lons, SMAP_lats, SMAP_time = nc.variables['wind'][:], nc.variables['lon'][:], nc.variables['lat'][:], nc.variables['minute'][:].astype(float)
            #convert SMAP time to datetime format
            start_date = datetime(2023, 10, 23, 0, 0, 0)
            SMAP_wind_asc, SMAP_wind_des = SMAP_wind[:,:,0], SMAP_wind[:,:,1]
            SMAP_time_asc, SMAP_time_des = SMAP_time[:,:,0], SMAP_time[:,:,1]
            lats_boolean = (SMAP_lats.data > AOI_lat_min) & (SMAP_lats.data < AOI_lat_max)
            lons_boolean = (SMAP_lons.data > AOI_lon_min) & (SMAP_lons.data < AOI_lon_max)
            lats_lons_boolean = lats_boolean[:, np.newaxis] & lons_boolean
            true_indices = np.where(lats_lons_boolean)
            if len(true_indices)>0:
                SMAP_wind_AOI = SMAP_wind[true_indices] 
                SMAP_lats_AOI = SMAP_lats[true_indices[0]] ; SMAP_lons_AOI = SMAP_lons[true_indices[1]] 
                SMAP_AOI = True
            
        if 'cyg' in file:
            if 'l1' in file:
                cygnss_lons = nc.variables['sp_lon'][:] ; cygnss_lats = nc.variables['sp_lat'][:]; cygnss_time_1D =nc.variables['ddm_timestamp_utc'][:]; nbrcs = nc.variables['brcs'][:]; 
                cygnss_time = np.tile(cygnss_time_1D[:, np.newaxis], 4)
                brcs = nc.variables['brcs']
            if 'l2' in file:
                cygnss_lons = nc.variables['lon'][:] ; cygnss_lats = nc.variables['lat'][:]; cygnss_wind= nc.variables['wind_speed'][:] ; cygnss_YSLF_wind= nc.variables['yslf_wind_speed'][:] ; cygnss_time =nc.variables['sample_time'][:]
                inc = nc.variables['incidence_angle']; rcg = nc.variables['range_corr_gain'][:]; uncertainty= nc.variables['yslf_wind_speed_uncertainty']; sc_num=nc.variables['spacecraft_num']
            #convert time to datetime format
            start_date = datetime(2023, 10, 22, 0, 0, 0, 499261) # need to change if changing start dates
            datetime_array = np.vectorize(lambda x: start_date + timedelta(seconds=x))(cygnss_time)
            
            lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
            lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
            time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end)
            AOI_boolean = lats_boolean & lons_boolean & time_boolean
            true_indices = np.where(AOI_boolean)
            if 'l2' in file:
                true_indices = np.where(AOI_boolean)[0]
            if len(true_indices)>0:
                if 'l2' in file:
                    cygnss_wind_AOI = cygnss_wind[true_indices] ; cygnss_YSLF_AOI = cygnss_YSLF_wind[true_indices]; cygnss_time_AOI = datetime_array[true_indices] ; sc_num_AOI = sc_num[true_indices]; inc_AOI = inc[true_indices] ; rcg_AOI = rcg[true_indices] ; uncertainty_AOI = uncertainty[true_indices]
                cygnss_lats_AOI = cygnss_lats[true_indices] ; cygnss_lons_AOI = cygnss_lons[true_indices] ; cygnss_time_AOI = datetime_array[true_indices]
                if 'l1' in file:
                    brcs_AOI = brcs[true_indices]
                    valid.append(brcs_AOI) # incase there's multiple DDMs to plot
                cygnss_AOI = True
            

# analysis of .25 deg collocations
# find how many boxes there are in AOI
AOI_lat_grid = int((AOI_lat_max - AOI_lat_min)/res)
AOI_lon_grid = int((AOI_lon_max - AOI_lon_min)/res)
AOI_grid = np.empty((AOI_lat_grid,AOI_lon_grid))

# Define spatial bins
lat_bins = np.arange(AOI_lat_min, AOI_lat_max, res)
lon_bins = np.arange(AOI_lon_min, AOI_lon_max, res)

# Use digitize to assign each point to a bin
lat_bin_indices = np.digitize(ascat_lats_AOI, lat_bins)
lon_bin_indices = np.digitize(ascat_lons_AOI, lon_bins)
cyg_lat_bin_indices = np.digitize(cygnss_lats_AOI, lat_bins)
cyg_lon_bin_indices = np.digitize(cygnss_lons_AOI, lon_bins)

# Create an empty array to store aggregated data, and a common mask
aggregated_ref_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
aggregated_cygnss_data = np.ma.zeros((len(lat_bins), len(lon_bins)))
aggregated_ref_time = np.full((len(lat_bins), len(lon_bins)), datetime(2000, 1, 1, 12, 0), dtype=object)
aggregated_cygnss_time = np.full((len(lat_bins), len(lon_bins)), datetime(2000, 1, 1, 12, 0), dtype=object)
mask_array = np.zeros((len(lat_bins), len(lon_bins))).astype(datetime)


# Aggregate ref data in each spatial bin
for i in range(len(lat_bins)):
    for j in range(len(lon_bins)):
        indices_in_bin = np.where((lat_bin_indices == i + 1) & (lon_bin_indices == j + 1))[0]
        if len(indices_in_bin) > 0:
            aggregated_ref_data[i, j] = np.mean(ascat_wind_AOI[indices_in_bin])
            timestamps = np.array([dt.timestamp() for dt in ascat_time_AOI[indices_in_bin]])
            mean_timestamp = np.mean(timestamps)
            aggregated_ref_time[i, j] = datetime.fromtimestamp(mean_timestamp)
        elif (len(indices_in_bin) == 0) :
            mask_array[i, j] = 1
        if aggregated_ref_data[i, j] is np.nan:
            mask_array[i, j] = 1
    
# 5. get CYGNSS around that measurement of threshold time interval
# Convert datetime values to timestamps
timestamps = np.array([dt.timestamp() for dt in ascat_time_AOI])
mean_timestamp = np.mean(timestamps)
average_datetime = datetime.fromtimestamp(mean_timestamp)

# Aggregate cygnss data in each spatial bin
for i in range(len(lat_bins)):
    for j in range(len(lon_bins)):
        indices_in_bin = np.where((cyg_lat_bin_indices == i + 1) & (cyg_lon_bin_indices == j + 1))[0]
        
        # Check if the time of each measurement is within x hours of the target time
        target_time = aggregated_ref_time[i, j]
        indices_in_time_window = np.where(np.abs(cygnss_time_AOI[indices_in_bin] - target_time) <= pd.Timedelta(hours=timing_threshold))[0]

        if len(indices_in_time_window) > 0:
            # condition to test if time collocates
            aggregated_cygnss_data[i, j] = np.mean(cygnss_YSLF_AOI[indices_in_time_window])
                # aggregated_time[i, j] = np.mean(ascat_time_AOI[indices_in_bin])
        elif (len(indices_in_time_window) == 0) :
            mask_array[i, j] = 1
        if aggregated_cygnss_data[i, j] is np.nan:
            mask_array[i, j] = 1

# convert masked array to boolean to mask the bad values
bool_arr = mask_array.astype(bool)
aggregated_ref_data.mask = bool_arr
aggregated_cygnss_data.mask = bool_arr

# make an array of the differences 
diff = aggregated_cygnss_data - aggregated_ref_data

# 6. stats of correlation, bias, KGE, RMSE
def stat_comparison(observed, simulated):
    r = np.corrcoef(observed.flatten(), simulated.flatten())[0, 1]
    alpha = np.std(simulated) / np.std(observed)
    beta = np.mean(simulated) / np.mean(observed)

    kge_value = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    
    bias = np.mean(observed -  simulated)
    RMSE = np.sqrt(np.mean((observed -  simulated)**2, axis=None))

    return kge_value,r,bias,RMSE

kge_value,r,bias,RMSE = stat_comparison(aggregated_cygnss_data,aggregated_ref_data)

plotting = True
if plotting:
    # # plotting
    import matplotlib.pyplot as plt
    from cartopy import crs as ccrs 
    # Create a Cartopy PlateCarree projection (cylindrical projection)
    projection = ccrs.PlateCarree()
    
    # Create a Matplotlib figure and axis
    fig, ax = plt.subplots(subplot_kw={'projection': projection})
     
    # im = ax.contourf(np.average(cygnss_wind,axis=0), origin='lower', extent=[cygnss_coverage[1], cygnss_coverage[3], cygnss_coverage[0], cygnss_coverage[2]], cmap='viridis',
    #                 transform=projection, aspect='auto')
    # im = ax.imshow(np.average(cygnss_wind_buoy,axis=0), origin='lower', extent=[buoy_area[1], buoy_area[3], buoy_area[0], buoy_area[2]], cmap='viridis',
    #                 transform=projection, aspect='auto')
    # im = ax.imshow(SMAP_wind_asc, origin='lower', extent=[SMAP_coverage[1], SMAP_coverage[3], SMAP_coverage[0], SMAP_coverage[2]], cmap='viridis',
    #                 transform=projection, aspect='auto')
    vmin = 10 ; vmax = 30
    if ascat_AOI:
        im = ax.scatter(ascat_lons_AOI, ascat_lats_AOI, c=ascat_wind_AOI, cmap='viridis',transform=projection, marker='s',vmin=vmin, vmax=vmax)
    
    if SMAP_AOI:
        im = ax.scatter(ascat_lons_AOI, ascat_lats_AOI, c=ascat_wind_AOI, cmap='viridis',transform=projection, marker='s',vmin=vmin, vmax=vmax)
    
    if cygnss_AOI:
        im = ax.scatter(cygnss_lons_AOI, cygnss_lats_AOI, c=cygnss_YSLF_AOI, cmap='viridis',transform=projection, vmin=vmin, vmax=vmax, marker='1')
        cbar = plt.colorbar(im, ax=ax, orientation='vertical',location='left', shrink=0.8)
        cbar.set_label('wind speed avg (m/s)')
        
    if len(diff)>0:
        x_grid, y_grid = np.add(np.meshgrid(lon_bins, lat_bins),0.12)
        im = ax.pcolormesh(x_grid, y_grid, diff, cmap='Reds_r', alpha=0.4)
        cbar = plt.colorbar(im, ax=ax, orientation='vertical', shrink=0.8)
        cbar.set_label('wind speed diff (m/s)')

    
    # Add coastlines and gridlines for better context
    ax.coastlines()
    buffer = 0.07
    ax.set_xlim(AOI_lon_min-buffer,AOI_lon_max+buffer)  # Set x-axis limits
    ax.set_ylim(AOI_lat_min-buffer,AOI_lat_max+buffer)  # Set y-axis limits
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    gls = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gls.top_labels=False   # suppress top labels
    gls.right_labels=False # suppress right labels
    ax.set_title(f'SMAP and CYGNSS at {str(AOI_time_start)[11:16]}-{str(AOI_time_end)[11:16]}', loc='right')
    
    plt.subplots_adjust(left=0.55, right=0.65, bottom=0.1, top=0.9, wspace=0.5)
    plt.tight_layout()
    
    # # saving file
    directory = r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Plots\DDM4.png'
    # fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    
    plt.show()
    
    # Plot the 2D array
    for DDM in valid:
        plt.imshow(valid[1][2][0], cmap='viridis', origin='lower', interpolation='none')
        plt.xlabel('doppler bins')
        plt.ylabel('delay bins')
        plt.show()

        # plt.savefig(directory,  dpi=1200, bbox_inches="tight")

