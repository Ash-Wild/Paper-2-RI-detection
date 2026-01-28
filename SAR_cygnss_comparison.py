# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:56:16 2024
Purpose: Compare SAR and CYGNSS for TC Lola
@author: Ashley
"""

import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import glob

comp = 'Ashley'
# SAR_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\SAR'
# cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\All_L2'
# SMAP_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\All_SMAP\RSS'
SAR_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\SAR'
cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\CYG L2'
SMAP_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\RSS WS'

time_threshold = 4
res = 0.25 # the size of the grid
kge_list,r_list,bias_list,RMSE_list = [],[],[],[]

# assumes that all dates of SAR images have a cygnss file to them - no dates are missed by cyg
sar_files_list = glob.glob(SAR_folder+ "\*.nc")
cyg_files_list = glob.glob(cyg_folder+ "\*.nc")
smap_files_list = glob.glob(SMAP_folder + "\*.nc")

def stat_comparison(observed, simulated):
    r = np.corrcoef(observed.flatten(), simulated.flatten())[0, 1]
    alpha = np.std(simulated) / np.std(observed)
    beta = np.mean(simulated) / np.mean(observed)

    kge_value = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    
    bias = np.mean(observed -  simulated)
    RMSE = np.sqrt(np.mean((observed -  simulated)**2, axis=None))

    return kge_value,r,bias,RMSE

# make a list of cygnss dates and SMAP dates
cyg_file_dates = []
for cyg_file in cyg_files_list:
    start_ind = cyg_file.find('ddmi.s')
    cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
    cyg_file_dates.append(cyg_datetime_object.date())
    
    
# have 2 options with SMAP - I could let it find the files, or I can make a netCDF of them all for it to slice into
smap_file_dates = []
for smap_file in smap_files_list:
    # start_ind = smap_file.find('SSS')
    # smap_datetime_object = datetime.strptime(smap_file[start_ind+10:start_ind+25], "%Y%m%dT%H%M%S")
    start_ind = smap_file.find('ly_')
    smap_datetime_object = datetime.strptime(smap_file[start_ind+3:start_ind+13], "%Y_%m_%d")
    smap_file_dates.append(smap_datetime_object.date())


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
    assert sar_date in smap_file_dates and sar_date in cyg_file_dates
    
    sar_wind = sar_nc.variables['sar_wind'][:] ; sar_lat = sar_nc.variables['latitude'][:] ; sar_lon = sar_nc.variables['longitude'][:]
    
    # define the AOI based on extent of SAR - assumed to be a simple box - maybe more complex options for if its star shaped
    AOI_lat_min,AOI_lat_max = np.min(sar_lat),np.max(sar_lat)
    AOI_lon_min,AOI_lon_max = np.min(sar_lon),np.max(sar_lon)
    AOI_time_start,AOI_time_end = sar_time-timedelta(hours=time_threshold),sar_time+timedelta(hours=time_threshold)
    
    
    #find and open the corresponding cygnss file
    cyg_date_ind = np.where([sar_time.date() == x for x in cyg_file_dates])[0]
    if sar_time.hour > 24-time_threshold:
        if cyg_date_ind != len(cyg_files_list)-1:
            cyg_date_list = [cyg_date_ind,cyg_date_ind+1]
    elif sar_time.hour < 0 + time_threshold: # if its possible there was similar CYG measurement on surrounding day
        if cyg_date_ind != 0:
            cyg_date_list = [cyg_date_ind,cyg_date_ind-1]
    else:
        cyg_date_list = [cyg_date_ind]
    
    # if there's multiple dates of measurements
    for cyg_ind in cyg_date_list:
        cyg_nc = nc.Dataset(cyg_files_list[cyg_ind[0]])
        cygnss_lons = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ; cygnss_YSLF_wind= cyg_nc.variables['yslf_wind_speed'][:] ; cygnss_time =cyg_nc.variables['sample_time'][:]
        inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
        
        #convert time to datetime format
    
        lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
        lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
        box_indices = np.where(lats_boolean & lons_boolean)[0]
        if len(box_indices)>0:
            cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
            date_format = "%Y-%m-%d %H:%M:%S"
            cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
            datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
            
            time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end)
            true_indices = np.where(time_boolean)[0]
            if len(true_indices)>0:
                cygnss_wind_AOI = cygnss_wind[box_indices[true_indices]] ; cygnss_YSLF_AOI = cygnss_YSLF_wind[box_indices[true_indices]]; cygnss_time_AOI = datetime_array[true_indices] ; sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
                cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]
                cyg_AOI = True 
            
    # need to add in code for whether there's another CYGNSS
    
    #find and open the corresponding smap file
    smap_date_ind = np.where([sar_date == x for x in smap_file_dates])[0]
    if sar_time.hour > 24-time_threshold:
        if smap_date_ind != len(smap_files_list)-1:
            smap_date_list = [smap_date_ind,smap_date_ind+1]
    elif sar_time.hour < 0 + time_threshold: # if its possible there was similar CYG measurement on surrounding day
        if smap_date_ind != 0:
            smap_date_list = [smap_date_ind,smap_date_ind-1]
    else:
        smap_date_list = [smap_date_ind]
        
    # if there's multiple dates of measurements
    for smap_ind in smap_date_list:
        smap_nc = nc.Dataset(smap_files_list[smap_ind[0]])
        smap_lons = smap_nc.variables['lon'][:]-180 ; smap_lats = smap_nc.variables['lat'][:]; smap_wind= np.mean(smap_nc.variables['wind'][:],axis=2) ; smap_time =np.mean(np.ma.masked_array(smap_nc.variables['minute'][:].data,mask=smap_nc.variables['wind'][:].mask), axis=2)
        
        #convert time to datetime format
        lats_boolean = (smap_lats.data > AOI_lat_min) & (smap_lats.data < AOI_lat_max)
        lons_boolean = (smap_lons.data > AOI_lon_min) & (smap_lons.data < AOI_lon_max)
        smap_wind_clipped = smap_wind[lats_boolean][:,lons_boolean]
        # print(np.where(smap_wind_clipped >0))
        
        smap_start_date = datetime.combine(smap_file_dates[smap_ind[0]], datetime.min.time())
        datetime_array = np.vectorize(lambda x: smap_start_date + timedelta(seconds=x*60))(smap_time[lats_boolean][:,lons_boolean])
        print(np.where(~datetime_array.mask))
        time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end) # error here
        true_lats = [] ; true_lons = []
        # Iterate through the 2D array and find indexes within the range
        for i, row in enumerate(datetime_array):
            for j, dt in enumerate(row):
                if type(dt) is datetime:
                    if AOI_time_start <= dt <= AOI_time_end:
                        true_lats.append(i) ; true_lons.append(j)
        # true_indices = np.where(time_boolean)[0]
        if len(true_lons)>0:
            smap_wind_AOI = smap_wind_clipped[true_lats,true_lons] ; smap_time_AOI = datetime_array[true_lats,true_lons]
            smap_lats_AOI = smap_lats[lats_boolean][true_lats] ; smap_lons_AOI = smap_lons[lons_boolean][true_lons]
            smap_AOI = True 
    
    if cyg_AOI:
        # comparing the two datasets
        # analysis of .25 deg collocations
        # find how many boxes there are in AOI
        AOI_lat_grid = int((AOI_lat_max - AOI_lat_min)/res)
        AOI_lon_grid = int((AOI_lon_max - AOI_lon_min)/res)
        AOI_grid = np.empty((AOI_lat_grid,AOI_lon_grid))
        
        # Define spatial bins
        lat_bins = np.arange(AOI_lat_min, AOI_lat_max, res)
        lon_bins = np.arange(AOI_lon_min, AOI_lon_max, res)
        
        # Use digitize to assign each point to a bin
        sar_lat_bin_indices = np.digitize(sar_lat, lat_bins)
        sar_lon_bin_indices = np.digitize(sar_lon, lon_bins)
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
                    aggregated_cygnss_data[i, j] = np.mean(cygnss_YSLF_AOI[indices_in_bin])
                    timestamps = np.array([dt.timestamp() for dt in cygnss_time_AOI[indices_in_bin]])
                    mean_timestamp = np.nanmean(timestamps)
                    aggregated_cygnss_time[i, j] = datetime.fromtimestamp(mean_timestamp)
                if (len(indices_in_bin) == 0):
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
        kge_value,r,bias,RMSE = stat_comparison(aggregated_cygnss_data,aggregated_ref_data)
        kge_list.append(kge_value); r_list.append(r); bias_list.append(bias); RMSE_list.append(RMSE)

# def tc_analysis(cyg,smap,sar): #def Tri_Col(pixels,Means):
#     # all to be 2D arrays of values that have been aggregated to the grid
#     #ensure that the pixels are already in Log form    
#     cov_matrix=np.cov(pixels,rowvar=False)
    
#     Log_RMSE=np.empty(3)
#     Log_RMSE[0]=cov_matrix[1,1]-(cov_matrix[1,2]*cov_matrix[1,3]/cov_matrix[2,3])
#     Log_RMSE[1]=cov_matrix[2,2]-(cov_matrix[1,2]*cov_matrix[2,3]/cov_matrix[1,3])
#     Log_RMSE[2]=cov_matrix[3,3]-(cov_matrix[1,3]*cov_matrix[2,3]/cov_matrix[1,2])
#     RMSE=np.empty(3)
#     for i in range(3):
#         RMSE[i]=Log_RMSE[i]*Means[i]
        
#     correlation = np.empty(3)
#     correlation[0]=(cov_matrix[1,2]*cov_matrix[1,3])/(cov_matrix[1,1]*cov_matrix[2,3])
#     correlation[1]=(cov_matrix[1,2]*cov_matrix[2,3])/(cov_matrix[2,2]*cov_matrix[1,3])
#     correlation[2]=(cov_matrix[1,3]*cov_matrix[2,3])/(cov_matrix[3,3]*cov_matrix[1,2])
    
#     return Log_RMSE, RMSE, correlation

        plotting = True
        if smap_AOI and cyg_AOI:
            # # plotting
            import matplotlib.pyplot as plt
            from cartopy import crs as ccrs 
            # Create a Cartopy PlateCarree projection (cylindrical projection)
            projection = ccrs.PlateCarree()
            
            # Create a Matplotlib figure and axis
            fig, ax = plt.subplots(subplot_kw={'projection': projection})
            vmin = 10 ; vmax = 80
            im = ax.pcolormesh(sar_lon, sar_lat, sar_wind, cmap='viridis',transform=projection, vmin=vmin, vmax=vmax)
            
            im = ax.scatter(smap_lons_AOI, smap_lats_AOI, c=smap_wind_AOI , cmap='viridis',transform=projection, marker='s',vmin=vmin, vmax=vmax, s=30)            
            # im = ax.scatter(cygnss_lons_AOI, cygnss_lats_AOI, c=cygnss_YSLF_AOI, cmap='viridis',transform=projection, vmin=vmin, vmax=vmax, marker='1')
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
            ax.set_title(f'SMAP, SAR, & CYGNSS at {str(AOI_time_start)[11:16]}-{str(AOI_time_end)[11:16]}', loc='right')
            
            plt.subplots_adjust(left=0.55, right=0.65, bottom=0.1, top=0.9, wspace=0.5)
            plt.tight_layout()
            
            # # saving file
            directory = r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Plots\DDM4.png'
            # fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
            
            plt.show()