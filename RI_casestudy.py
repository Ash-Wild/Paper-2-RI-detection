# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 08:43:49 2025

@author: ashle
"""
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import glob
from scipy.spatial.distance import cdist
import pandas as pd
import cfgrib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

comp = 'ashle'
event_list = [['Hinnamnor_2022',428]] # list of RI events from RI ,   ['Mocha_2023',1569]
# hinnamnor v5 '2022239N22150' 428
# kyarr 2019296N15066 210
# Bualoi 2019290N08169 200
interp_interval= 20 # minutes to interpolate best-track
MAX_DISTANCE = 80  # Maximum distance in kilometers
buffer = MAX_DISTANCE/100 # degrees buffer to search in storm_centric

variables = ['eff_scatter', 'brcs','sp_inc_angle', 'ddm_snr','gps_eirp', 'rx_to_sp_range','tx_to_sp_range','quality_flags','sv_num', 'ddm_ant','nst_att_status', 'ddm_nbrcs',]
counter = 0
num_valid = 0

def haversine_vectorized(latlon1, latlon2):
    """
    Compute the great-circle distance between two sets of (lat, lon) coordinates using vectorized operations.
    """
    R = 6371.0  # Earth radius in km
    lat1, lon1 = np.radians(latlon1[:, 0]), np.radians(latlon1[:, 1])
    lat2, lon2 = np.radians(latlon2[:, 0]), np.radians(latlon2[:, 1])

    dlat = lat2[:, None] - lat1[None, :]
    dlon = lon2[:, None] - lon1[None, :]

    a = np.sin(dlat / 2) ** 2 + np.cos(lat1[None, :]) * np.cos(lat2[:, None]) * np.sin(dlon / 2) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    return (R * c).T  # Transpose to match expected shape
# Function for vectorized linear interpolation
def linear_interpolation_vectorized(timestamps, variable_input, target_times): # this needs to be improved
    indices = np.searchsorted(timestamps, target_times) - 1
    if indices == -1: # if it is before the first speed, take the first speed
        return variable_input[0]
    if indices == len(timestamps)-1: # if it is after the last reading, take the last speed 
        return variable_input[-1]
    else: # for all else inbetween 
        before_times = timestamps[indices]
        time_diff = (target_times - before_times) / np.timedelta64(1, 's')
        total_time_interval = (timestamps[indices+1] - timestamps[indices]) / np.timedelta64(1, 's')
        variable_diff = variable_input[indices + 1] - variable_input[indices]
        interpolated_speed = variable_input[indices] + (time_diff / total_time_interval) * variable_diff
        if np.isnan(interpolated_speed):
            print('error interp')
        return interpolated_speed
    
for case_study in event_list:
    RI_file = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v4.nc'
    RI_nc = nc.Dataset(RI_file)
    time_int = RI_nc.variables['times'][case_study[1]]
    nans = ~np.isnan(time_int).data # remove nans from loaded nc
    time = np.datetime64('1970-01-01T00:00:00') + time_int.astype('timedelta64[s]') #- np.timedelta64(7,'D')
    tc_time = time[nans]
    tc_storm_id = RI_nc.variables['storm'][case_study[1]]
    tc_vmax = RI_nc.variables['vmax'][case_study[1]][nans]*0.514 # convert from knots to m/s
    tc_latitude = RI_nc.variables['latitude'][case_study[1]][nans]
    tc_longitude = RI_nc.variables['longitude'][case_study[1]][nans]
    ri_rd = RI_nc.variables['ri_rd'][case_study[1]]
    durations=RI_nc.variables['durations'][case_study[1]]
    intensities=RI_nc.variables['intensities'][case_study[1]]*0.514
    index_before=RI_nc.variables['index_before'][case_study[1]]
    index_after=RI_nc.variables['index_after'][case_study[1]]
    RI_nc.close()
    
    #create the interpolated time of the best-track
    interval = np.timedelta64(interp_interval, 'm')
    start = tc_time[0];    end = tc_time[-1]
    tc_time_interpolated = np.arange(start, end + interval, interval)
    
    # create arrays to store interpolated values
    tc_vmax_interpolated = np.zeros(np.size(tc_time_interpolated))
    tc_lat_interpolated = np.zeros(np.size(tc_time_interpolated))
    tc_lon_interpolated = np.zeros(np.size(tc_time_interpolated))

    for i in range(len(tc_time_interpolated)): 
        tc_vmax_interpolated[i] = linear_interpolation_vectorized(tc_time, tc_vmax, tc_time_interpolated[i])
        tc_lat_interpolated[i] = linear_interpolation_vectorized(tc_time, tc_latitude, tc_time_interpolated[i])
        tc_lon_interpolated[i] = linear_interpolation_vectorized(tc_time, tc_longitude, tc_time_interpolated[i])
    tc_coords = np.column_stack((tc_lat_interpolated, tc_lon_interpolated))  # Shape (M, 2)

    
    cyg_AOI_time_start,cyg_AOI_time_end = np.nanmin(tc_time),np.nanmax(tc_time)
    cyg_AOI_time_start = cyg_AOI_time_start.astype('datetime64[s]').astype('O')
    cyg_AOI_time_end = cyg_AOI_time_end.astype('datetime64[s]').astype('O')
    
    AOI_lat_min,AOI_lat_max = np.nanmin(tc_latitude)-buffer,np.nanmax(tc_latitude)+buffer
    AOI_lon_min,AOI_lon_max = np.nanmin(tc_longitude)-buffer,np.nanmax(tc_longitude)+buffer
    
    
    # cyg 
    case_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\Casestudy\\'+ case_study[0]
    cyg_files_list = np.asarray(glob.glob(case_folder+ r"\cyg_l1\*.nc"))
    
    coincident_meas= {i: np.ma.MaskedArray([]) for i in variables}
    cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([],dtype='datetime64') ; nearest_tc_coords=np.ma.MaskedArray([]);
    cyg_lats_all=np.ma.MaskedArray([]); cyg_lons_all=np.ma.MaskedArray([]); cyg_times_all=np.ma.MaskedArray([],dtype='datetime64') 
    for cyg_file_path in cyg_files_list:
        cyg_nc = nc.Dataset(cyg_file_path)
        longitudes = cyg_nc.variables['sp_lon'][:].flatten() ; cygnss_lats = cyg_nc.variables['sp_lat'][:].flatten()
        cygnss_lons = ((longitudes + 180) % 360) - 180
        cygnss_time =np.repeat(cyg_nc.variables['ddm_timestamp_utc'][:],4) ; cyg_start_time_str = cyg_nc.variables['ddm_timestamp_utc'].units[14:33]
        
    
        # find where the TC is for the date
        lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
        lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
        box_indices = np.where(lats_boolean & lons_boolean)[0] 
        if len(box_indices)>0:
            meas_lats = cygnss_lats[box_indices]
            meas_lons = cygnss_lons[box_indices]
            
            # temporal matching - time within x number of minutes
            date_format = "%Y-%m-%d %H:%M:%S"
            cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
            cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
            meas_times = np.array(cyg_datetime_array, dtype='datetime64')
            #save all current measures
            cyg_lats_all=np.append(cyg_lats_all,meas_lats) ; cyg_lons_all = np.append(cyg_lons_all,meas_lons) 
            cyg_times_all=np.append(cyg_times_all,meas_times)
            
            # Convert to NumPy arrays if not already
            meas_coords = np.column_stack((meas_lats, meas_lons))  # Shape (N, 2)
            # Compute the pairwise distance matrix
            distances = haversine_vectorized(meas_coords, tc_coords)
            # Compute pairwise time difference matrix (absolute time difference in hours)
            time_diffs = np.abs(meas_times[:, None] - tc_time_interpolated[None, :])
            # Find measurements where at least one track point is within 100 km & 30 minutes
            # print(np.timedelta64(np.min(time_diffs),'m'))
            valid_mask = (np.abs(distances) <= MAX_DISTANCE) & (time_diffs <= np.timedelta64(interp_interval, 'm'))
            # Select valid measurement indices
            valid_ind = np.where(valid_mask)
            selected_indices = np.unique(valid_ind[0]) # to remove duplicates
            
            
            # num_valid += len(selected_indices)
            
            if len(selected_indices) >0:
                #save meas that meet criteria
                cyg_lats=np.concatenate((cyg_lats,meas_lats[selected_indices])) ; cyg_lons = np.concatenate((cyg_lons,meas_lons[selected_indices]) )
                cyg_times=np.concatenate((cyg_times,meas_times[selected_indices]))
                # cyg_nearest = np.concatenate((cyg_nearest,np.avg(np.where(valid_mask)[1]) ))
                # selected_time_diffs = time_diffs[valid_mask]
                # selected_distances =  distances[valid_mask]
                
                for i in variables:
                    try:
                        if i in ['eff_scatter', 'brcs']:
                            var = cyg_nc.variables[i][:]
                            var = var.reshape(var.shape[0] * var.shape[1], var.shape[2], var.shape[3])
                        else:
                            var = cyg_nc.variables[i][:].flatten()
                        assert len(var) == len(cygnss_lons)
                        if len(coincident_meas[i]) == 0:
                            coincident_meas[i] = var[box_indices[selected_indices]]
                        else:
                            coincident_meas[i]=np.concatenate((coincident_meas[i],var[box_indices[selected_indices]]), axis=0) 
                    except AssertionError:
                        var=np.repeat(var,4)
                        assert len(var) == len(cygnss_lons)
                        coincident_meas[i]=np.concatenate((coincident_meas[i],var[box_indices][selected_indices]), axis=0) 
                        
                # surely faster way of doing this
                y_list=[]
                x_list=[]
                y_avg=[]
                for i in range(len(valid_ind[0])):
                    x_ind=valid_ind[0][i] ; y_ind = valid_ind[1][i]
                    if i == 0:
                        x_previous = x_ind
                        y_list.append(y_ind)
                        if i==len(valid_ind[0])-1:
                            x_list.append(x_previous)
                            y_avg.append(int(np.nanmean(y_list)))    
                    elif i != 0:    
                        if x_ind==x_previous:
                            y_list.append(y_ind)
                            if i==len(valid_ind[0])-1:
                                x_list.append(x_previous)
                                y_avg.append(int(np.nanmean(y_list)))                       
                        else:
                            x_list.append(x_previous)
                            y_avg.append(int(np.nanmean(y_list)))
                            y_list=[y_ind]
                            x_previous = x_ind    
                            if i==len(valid_ind[0])-1:
                                x_list.append(x_previous)
                                y_avg.append(int(np.nanmean(y_list)))    
                    else:
                        raise ValueError
                    
                cyg_valid_meas_inds = np.asarray(x_list) # same as selected ind
                nearest_tc_coord_ind = np.asarray(y_avg) 
                if len(nearest_tc_coord_ind) != len(meas_times[selected_indices]):
                    print('error')
                nearest_tc_coords=np.concatenate((nearest_tc_coords,nearest_tc_coord_ind))
        cyg_nc.close()

    # print(num_valid, len(cyg_lats))
     
    # # SWH 
    # cms_files_list = np.asarray(glob.glob(case_folder+ r"\cms\*.csv"))
    # cms_lats=np.ma.MaskedArray([]); cms_lons=np.ma.MaskedArray([]); cms_times=np.ma.MaskedArray([],dtype='datetime64') ; cms_swh=np.ma.MaskedArray([])
    # cms_lats_all=np.ma.MaskedArray([]); cms_lons_all=np.ma.MaskedArray([]); cms_times_all=np.ma.MaskedArray([],dtype='datetime64') 

    # for cms in cms_files_list:
    #     df = pd.DataFrame(pd.read_csv(cms,skiprows=5))
    #     cms_lat =np.asarray(df['latitude']) ; cms_lon = np.asarray(df['longitude']) ; 
    #     cms_time = np.asarray(df['time'],dtype='datetime64')
        
    #     # limitto box_indices
    #     lats_boolean = (cms_lat > AOI_lat_min) & (cms_lat < AOI_lat_max)
    #     lons_boolean = (cms_lon > AOI_lon_min) & (cms_lon < AOI_lon_max)
    #     box_indices = np.where(lats_boolean & lons_boolean)[0] 
    #     cms_lat = cms_lat[box_indices]
    #     cms_lon = cms_lon[box_indices]
    #     cms_time = cms_time[box_indices]
        
    #     cms_lats_all=np.append(cms_lats_all,cms_lat) ; cms_lons_all = np.append(cms_lons_all,cms_lon) 
    #     cms_times_all=np.append(cms_times_all,cms_time) 
        
    #     cms_coords = np.column_stack((cms_lat, cms_lon))  # Shape (N, 2)
    #     # Compute the pairwise distance matrix
    #     distances = haversine_vectorized(cms_coords, tc_coords)
    #     # Compute pairwise time difference matrix (absolute time difference in hours)
    #     time_diffs = np.abs(cms_time[:, None] - tc_time_interpolated[None, :])
    #     # Find measurements where at least one track point is within 100 km & 30 minutes
    #     valid_mask = np.any((distances <= MAX_DISTANCE + 100) & (time_diffs <= np.timedelta64(interp_interval, 'm')*3), axis=1)
    #     # valid_mask = distances <= MAX_DISTANCE + 500
        
    #     # Select valid measurement indices
    #     selected_indices = np.where(valid_mask)[0]
    #     if len(selected_indices)>0:
    #         cms_lats=np.append(cms_lats,cms_lat[selected_indices]) ; cms_lons = np.append(cms_lons,cms_lon[selected_indices]) 
    #         cms_times=np.append(cms_times,cms_time[selected_indices]) ; cms_swh = np.append(cms_swh, np.asarray(df['value'])[box_indices][selected_indices])
    # # COde for 2x2 L4 daily data
    # if len(cms_lats)==0:
    #     cms_files_list = np.asarray(glob.glob(cyg_folder+ r"\cms\*.nc"))
    #     for cms in cms_files_list:
    #         cms_nc = nc.Dataset(cms)
    #         cms_lat = cms_nc.variables['latitude'][:] ; cms_lon = cms_nc.variables['longitude'][:] ; cms_time = cms_nc.variables['time']
            
    #         cms_swh = cms_nc.variables['VAVH_DAILY_MAX']
    #         cms_nc.close
            
    #         lats_boolean = (cms_lat.data > AOI_lat_min) & (cms_lat.data < AOI_lat_max)
    #         lons_boolean = (cms_lon.data > AOI_lon_min) & (cms_lon.data < AOI_lon_max)
    #         box_indices = np.where(lats_boolean & lons_boolean)[0] 
         
    # IMERG
    imerg_files_list = np.asarray(glob.glob(case_folder+ r"\IMERG\*.nc4"))
    cyg_nc = nc.Dataset(imerg_files_list[0])
    imerg_lons = cyg_nc.variables['lon'][:] ; imerg_lats = cyg_nc.variables['lat'][:] ; 
    imerg_start_time_str = cyg_nc.variables['time'].units[14:33]
    cyg_nc.close
    imerg_times = np.ma.zeros((len(imerg_files_list))); imerg_precip = np.ma.zeros((len(imerg_files_list),len(imerg_lons),len(imerg_lats))); 
    imerg_precip_quality = np.ma.zeros((len(imerg_files_list),len(imerg_lons),len(imerg_lats)));
    for i,cyg_file_path in enumerate(imerg_files_list):
        cyg_nc = nc.Dataset(cyg_file_path)
        imerg_precip[i]= cyg_nc.variables['precipitation'][:] ; imerg_precip_quality[i]= cyg_nc.variables['precipitationQualityIndex'][:]
        imerg_times[i] =cyg_nc.variables['time'][0] 
        cyg_nc.close
    imerg_times = np.datetime64(imerg_start_time_str)+ imerg_times.astype('timedelta64[s]')
        
    # # ERA5
    # era5_dataset = cfgrib.open_dataset(case_folder+r'\1e6070f1609ff2b350fb377c94099b51.grib')
    # era5_lat = era5_dataset.latitude.values ; era5_lon  = era5_dataset.longitude.values
    # era5_time = era5_dataset.time.values.astype('datetime64[h]')
    # meanSea=era5_dataset.meanSea.values ; mslp=era5_dataset.msl.values
    # u_wind = era5_dataset.u10.values ; v_wind=era5_dataset.v10.values
    # height = era5_dataset.heightAboveGround.values
    # temp_2m= era5_dataset.t2m.values-273 ; surface_pressure= era5_dataset.sp.values
    # era5_wind = np.hypot(np.abs(u_wind), np.abs(v_wind))

    # analysing
    # classify measurements by quadrants
    # Example TC track interpolated positions (from your previous interpolation)
    tc_lats = tc_lat_interpolated
    tc_lons = tc_lon_interpolated
    # Example measurement points (sp_lat, sp_lon)
    meas_lats = cyg_lats
    meas_lons = cyg_lons
    meas_times = cyg_times
    # Compute TC motion vector (dx/dt, dy/dt) at each interpolated time step
    tc_dx = np.diff(tc_lons)
    tc_dy = np.diff(tc_lats)
    # Approximate TC heading angle (last step)
    tc_heading = np.arctan2(tc_dy, tc_dx)  # Angle in radians
    # Compute relative position of each measurement w.r.t. last TC position
    nearest_tc_coords = nearest_tc_coords.astype(int)-1
    rel_x = meas_lons - tc_lons[nearest_tc_coords] 
    rel_y = meas_lats - tc_lats[nearest_tc_coords]
    # Rotate coordinates to TC motion frame
    rot_x = rel_x * np.cos(tc_heading[nearest_tc_coords]) + rel_y * np.sin(tc_heading[nearest_tc_coords])  # Along TC motion direction
    rot_y = -rel_x * np.sin(tc_heading[nearest_tc_coords]) + rel_y * np.cos(tc_heading[nearest_tc_coords])  # Perpendicular to motion
    # Assign quadrants based on rotated coordinates
    quadrants = np.full(rot_x.shape, 'Unknown', dtype=object)  # Initialize array
    quadrants[(rot_x >= 0) & (rot_y >= 0)] = 'Front-Left'
    quadrants[(rot_x >= 0) & (rot_y < 0)] = 'Front-Right'
    quadrants[(rot_x < 0) & (rot_y >= 0)] = 'Rear-Left'
    quadrants[(rot_x < 0) & (rot_y < 0)] = 'Rear-Right'
    quadrant_colors = {'Front-Right': 'red', 'Front-Left': 'blue', 'Rear-Right': 'orange', 'Rear-Left': 'purple'}
    
    # # Plot of nearest setup
    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(figsize=(8, 6))
    # # Plot measurement points
    # ax.scatter(meas_lons, meas_lats, color='gray', label='Measurements')
    # # Plot TC track
    # ax.plot(tc_lons, tc_lats, color='blue', marker='o', label='TC Track')
    # # Draw lines from each measurement to its nearest TC coordinate
    # for i in range(len(meas_lats)):
    #     tc_i = nearest_tc_coords[i]
    #     ax.plot([meas_lons[i], tc_lons[tc_i]], [meas_lats[i], tc_lats[tc_i]], 
    #             color=quadrant_colors[quadrants[i]], linewidth=0.8, alpha=0.6)
    # # Labels and legend
    # ax.set_xlabel("Longitude")
    # ax.set_ylabel("Latitude")
    # ax.legend()
    # ax.set_title("Measurement to Nearest TC Position")
    # plt.show()
    
    # # number in each distance
    # # Convert to NumPy arrays if not already
    # meas_coords = np.column_stack((meas_lats, meas_lons))  # Shape (N, 2)
    # # Compute the pairwise distance matrix
    # distances = haversine_vectorized(meas_coords, tc_coords)
    # # Compute pairwise time difference matrix (absolute time difference in hours)
    # time_diffs = np.abs(meas_times[:, None] - tc_time_interpolated[None, :])
    # # Find measurements where at least one track point is within 100 km & 30 minutes
    # # print(np.timedelta64(np.min(time_diffs),'m'))
    # dist_increment = 10
    # for dist in range(0,100,dist_increment):
    #     valid_mask = ((np.abs(distances) > dist-dist_increment) & (np.abs(distances) <= dist)) & (time_diffs <= np.timedelta64(interp_interval, 'm'))
    #     # Select valid measurement indices
    #     selected_indices = np.unique(np.where(valid_mask)[0]) # to remove duplicates
    #     print(len(selected_indices),dist)
    #     # now need to remove distance rows identified
        
    # comparing quadrants
    # sort the arrays and dict
    valid_flag = coincident_meas['quality_flags'] % 2 == 0
    sorted_indices = np.argsort(cyg_times[valid_flag])
    meas_times = cyg_times[valid_flag][sorted_indices]
    meas_lats = cyg_lats[valid_flag][sorted_indices]
    meas_lons = cyg_lons[valid_flag][sorted_indices]
    quadrants = quadrants[valid_flag][sorted_indices]
    meas_data=coincident_meas.copy()
    for key in coincident_meas:
        meas_data[key] = coincident_meas[key][valid_flag][sorted_indices]
    brcs = meas_data['brcs']; eff_scatt = meas_data['eff_scatter']
    # Compute min and max along axis 0
    data_min = np.min(brcs, axis=(1, 2))
    data_max = np.max(brcs, axis=(1, 2))
    # Avoid division by zero
    range_ = data_max - data_min
    range_[range_ == 0] = np.nan  # or np.nan if you want to keep those locations masked
    # Normalize
    data_min = data_min[:, None, None]  # Shape becomes (359, 1, 1)
    range_   = range_[:, None, None]
    normalized_data = np.where(range_ == 0, np.nan, (brcs - data_min) / range_)
    mid_time = tc_time_interpolated[len(tc_time_interpolated) // 2]  # middle of the interpolated TC times
    ri_start_time = tc_time[index_before]
    ri_end_time = tc_time[len(tc_time)-index_after-1]
    assert ri_end_time - ri_start_time == np.timedelta64(int(durations),'h')
    plot_times = [cyg_AOI_time_start,ri_start_time,mid_time,ri_end_time,cyg_AOI_time_end]
    plot_times_names = ['12hrs before','Initial RI','Later RI','12hrs after']
    nbrcs = brcs; vmin = 0;  vmax = np.nanmax(np.mean(nbrcs))*5
    extra_var = meas_data['gps_eirp']
    extra_var2 = meas_data['sp_inc_angle']
    fig, axes = plt.subplots(len(plot_times)-1, 4, figsize=(13, (len(plot_times)-1)*4))
    img_handles = []
    for i, q in enumerate(quadrant_colors):
        for j in range(len(plot_times)-1):
            q_mask = quadrants == q
            time_segment_start=plot_times[j] ; time_segment_end=plot_times[j+1]
            mask = (meas_times > time_segment_start) & (meas_times <= time_segment_end) & q_mask
            ax = axes[j, i]
            avg_nbrcs = np.nanmean(nbrcs[mask], axis=0)
            im = ax.imshow(avg_nbrcs, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
            count = np.nansum(mask)
            ax.set_title(f'{q} {plot_times_names[j]}\nN={count}', color=quadrant_colors[q])
            img_handles.append(im)
            # Add inset histogram to the right of the subplot
            inset_ax = inset_axes(ax, width="30%", height="40%", loc='lower right',
                                  bbox_to_anchor=(0.55, 0.6, 1, 1),
                                  bbox_transform=ax.transAxes, borderpad=0)
            inset_ax.hist(extra_var[mask], bins=10, color='gray', edgecolor='black')
            # inset_ax.set_xticks([0,45,90]);   # inset_ax.set_yticks([])
            inset_ax.set_title("Inc. angle", fontsize=7)
            inset_ax = inset_axes(ax, width="30%", height="40%", loc='lower right',
                                  bbox_to_anchor=(0.55, 0., 1, 1),
                                  bbox_transform=ax.transAxes, borderpad=0)
            inset_ax.hist(extra_var2[mask], bins=10, color='gray', edgecolor='black')
            # inset_ax.set_xticks([0,750,1500]);   # inset_ax.set_yticks([])
            inset_ax.set_title("GPS EIRP", fontsize=7)
            after_mask = (meas_times >= mid_time) & q_mask
            
            # ax = axes[2, i]
            # change_ddm_mean = np.abs(avg_nbrcs_after - avg_nbrcs_before)
            # im = ax.imshow(change_ddm_mean, cmap='viridis', origin='lower', vmin=vmin, vmax=vmax)
            # ax.set_title(f'{q} (Change)', color=quadrant_colors[q])
            # img_handles.append(im)
            
    plt.tight_layout(rect=[0.03, 0.10, 0.97, 0.93])    # Adjust spacing to leave space at the bottom
    cbar = fig.colorbar(img_handles[0], cax=fig.add_axes([0.35, 0.08, 0.3, 0.02]), orientation='horizontal') # Manually add a smaller horizontal colorbar
    cbar.set_label('Average BRCS (m^2)')
    plt.show()
    
    # # quality flags    
    # vhex = np.vectorize(hex)
    # hex_array = vhex(coincident_meas['quality_flags'])
    # quality_flag_hex = np.array([int(h, 16) for h in hex_array])
    # flags = [    (0x00000001, 'poor_overall_quality'),    (0x00000002, 's_band_powered_up'),    (0x00000004, 'small_sc_attitude_err'),
    # (0x00000008, 'large_sc_attitude_err'),    (0x00000010, 'black_body_ddm'),    (0x00000020, 'ddmi_reconfgured'),
    # (0x00000040, 'spacewire_crc_invalid'),    (0x00000080, 'ddm_is_test_pattern'),    (0x00000100, 'channel_idle'),
    # (0x00000200, 'low_confdence_ddm_noise_foor'),    (0x00000400, 'sp_over_land'),    (0x00000800, 'sp_very_near_land'),
    # (0x00001000, 'sp_near_land'),    (0x00004000, 'large_step_lna_temp'),    (0x00008000, 'direct_signal_in_ddm'),
    # (0x00010000, 'low_confdence_gps_eirp_estimate'),    (0x00020000, 'rf_detected'),    (0x00040000, 'brcs_ddm_sp_bin_delay_error'),
    # (0x00080000, 'brcs_ddm_sp_bin_dopp_error'),    (0x00100000, 'neg_brcs_value_used_for_nbrcs'),
    # (0x00200000, 'gps_pvt_sp3_error'),    (0x00400000, 'sp_non_existent_error'),    (0x00800000, 'brcs_lut_range_error'),
    # (0x01000000, 'ant_data_lut_range_error'),    (0x02000000, 'bb_framing_error'),    (0x04000000, 'fsw_comp_shift_error'),
    # (0x08000000, 'low_quality_gps_ant_knowledge'),    (0x10000000, 'sc_altitude_out_of_nominal_range')]
    # # Function to decode a single hex value into list of flags
    # def decode_flags(hex_val):
    #     return [name for bit, name in flags if hex_val & bit]
    # decoded_flags = np.array([decode_flags(val) for val in quality_flag_hex], dtype=object)
    
    # # ax.scatter(cyg_lons_all,cyg_lats_all)
    # ax.scatter(cyg_lons,cyg_lats)
    # ax.scatter(tc_lon_interpolated,tc_lat_interpolated)
    # ax.scatter(cms_lons,cms_lats)
    # plt.show()      
    
    plotting = False
    if plotting:
        dlat = np.diff(tc_lat_interpolated, prepend=tc_lat_interpolated[0])
        dlon = np.diff(tc_lon_interpolated, prepend=tc_lon_interpolated[0])
        
        # Create plot
        import matplotlib.pyplot as plt
        import numpy as np
        import matplotlib.animation as animation
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors
        from matplotlib.patches import Circle
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        import matplotlib.patches as mpatches

        buffer=1.1
        latitudes = tc_lat_interpolated ; longitudes = tc_lon_interpolated; times = tc_time_interpolated
        winds = tc_vmax_interpolated
        fig, ax = plt.subplots(figsize=(16, 6))

        # Set plot limits
        ax.set_xlim(min(longitudes) - buffer, max(longitudes) + buffer)
        ax.set_ylim(min(latitudes) - buffer, max(latitudes) + buffer) 
        # ax.add_feature(cfeature.LAND, zorder=0, facecolor='lightgray')
        # ax.add_feature(cfeature.COASTLINE, zorder=1)
        # ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())                                    
        # Create a scalar mappable for the colorbar (needed for animation)
        norm = plt.Normalize(vmin=15, vmax=45)
        cmap = plt.cm.coolwarm  # Choose a colormap
        # Create colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # Required for colorbar
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', label='Wind Speed (m/s)')

        # Time label
        # Function to update frame
        def update_frame(frame):
            ax.clear()
        
            # Reset plot limits
            ax.set_xlim(min(longitudes) - buffer, max(longitudes) + buffer)
            ax.set_ylim(min(latitudes) - buffer, max(latitudes) + buffer)
            
            mask = (cyg_times_all <= times[frame]) & (cyg_times_all >= times[frame] - np.timedelta64(4,'h'))
            ax.scatter(cyg_lons_all[mask],cyg_lats_all[mask],c='grey', marker='x',s=10)
            # mask = (cms_times_all <= times[frame]) & (cms_times_all >= times[frame] - np.timedelta64(4,'h'))
            # ax.scatter(cms_lons_all[mask],cms_lats_all[mask],c='black', marker='o',s=10)
        
            # Plot arrows for wind direction with color mapped to wind speed
            ax.quiver(longitudes[:frame + 1], latitudes[:frame + 1], 
                                dlon[:frame + 1], dlat[:frame + 1], 
                                winds[:frame + 1],  # Use wind speed for color
                                cmap='coolwarm', width=0.005, norm=norm, angles='xy', scale_units='xy')
            
            # add rings around tc
            # for lat, lon in zip(longitudes[frame], latitudes[frame]):
            ring = Circle((longitudes[frame], latitudes[frame]), 0.95, fill=False, color='black', linestyle='--', linewidth=1)
            ax.add_patch(ring)
            ring = Circle((longitudes[frame], latitudes[frame]), 0.05, fill=False, color='black', linestyle='-', linewidth=0.5)
            ax.add_patch(ring)

            # valid cyg
            mask = (meas_times <= times[frame]) & (meas_times>= times[frame] - np.timedelta64(4,'h'))  
            color_list = [quadrant_colors[q] for q in quadrants[mask]]
            ax.scatter(meas_lons[mask], meas_lats[mask], c=color_list , marker='o',s=20)
                
            # # SWH
            # mask = (cms_times <= times[frame]) & (cms_times>= times[frame] - np.timedelta64(4,'h'))   # Select measurements up to this time
            # ax.scatter(cms_lons[mask], cms_lats[mask], c='blue', marker='*',s=20)
            # # ERA5
            # mask = (era5_time <= times[frame]) & (era5_time>= times[frame] - np.timedelta64(50,'m'))
            # mesh = ax.pcolormesh(era5_lon, era5_lat, np.squeeze(era5_wind[mask]), cmap='coolwarm', shading='auto', alpha=0.6)
            
            # Update time label
            ax.text(0.98, 0.95, f"Time: {times[frame]}", 
             transform=plt.gca().transAxes,  # Places text relative to axes
             fontsize=12, ha='right', color='black', bbox=dict(facecolor='white', alpha=0.6))

            # time_text = ax.text(0.1, 0.90, '', fontsize=12)
            # time_text.set_text(f'Time: {times[frame]} hours')
        
            # Titles & labels
            ax.set_title('ID: {0},Duration: {1}, Intensity: {2}'.format(tc_storm_id, durations, intensities))
            legend_handles = [mpatches.Patch(color=color, label=label) for label, color in quadrant_colors.items()]
            ax.legend(handles=legend_handles, loc='lower right', title='Quadrants')
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

        # Create animation
        pause_frames = 10 
        ani = animation.FuncAnimation(fig, update_frame, frames=list(range(len(times)))+[len(times) - 1] * pause_frames, repeat=True)
        
        # Save GIF
        ani.save(r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\gifs\L1_{1}_{2}.gif'.format(comp, tc_storm_id,case_study[1]), fps=2)
        
        # Show plot
        plt.show()