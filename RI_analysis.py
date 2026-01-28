# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:46:38 2024

@author: ashle
"""
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import pandas as pd
import glob
from Storm_timeline_visualiser import plot_timeline_refined, plot_distribution # for visualising timeline

# from scipy.spatial import cKDTree
comp = 'Ashley'

cyg_version = 'noaa_1.2' # 'noaa_1.2' or 'storm_centric'
buffer = 0.5 # degrees buffer to search in storm_centric
interp_interval= 20 # minutes to interpolate best-track
MAX_DISTANCE = 80  # Maximum distance in kilometers
case_study_ids = ['2022239N22150'] #['2019290N08169','2022239N22150','2019296N15066','2024279N21265','2023290N12256','2019077S13123']
case_study_inds = [428, 200, 210, 140, 590, 680]
variables = ['incidence_angle', 'rx_gain', 'snr', 'range_corr_gain','sample_flags', 'num_ddms_utilized','nbrcs_mean_corrected','nbrcs_mean','wind_speed_uncertainty','wind_speed','cyg_vmax_error','cyg_dist','time_since_ri'] # 'ddm_sample_index', 'ddm_channel',
noaa_cyg_folder = r'D:\Phd_data\cyg_noaa_1.2' # r'C:\Users\\'+comp+'\\Documents\CYGNSS_NOAA'  # 
merged_cyg_folder = r'D:\Phd_data\cyg_storm_centric' #r'C:\Users\\'+comp+'\\Documents\cyg_storm_centric' 
noaa_cyg_files_list = np.asarray(glob.glob(noaa_cyg_folder+ r"\*.nc"))
merged_cyg_files_list = np.asarray(glob.glob(merged_cyg_folder+ r"\*.nc"))

counter=0
RI_file = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v5.nc'
RI_nc = nc.Dataset(RI_file)

time_int = RI_nc.variables['times'][:]
time = np.datetime64('1970-01-01T00:00:00') + time_int.astype('timedelta64[s]')
storm_id = RI_nc.variables['storm'][:]
# nunique,counts=np.unique(storm_id, return_counts=True)
# print(len(nunique))
vmax = RI_nc.variables['vmax'][:]*0.514
latitude = RI_nc.variables['latitude'][:]
longitude = RI_nc.variables['longitude'][:]
ri_rd = RI_nc.variables['ri_rd'][:]
durations=RI_nc.variables['durations'][:]
intensities=RI_nc.variables['intensities'][:]*0.514
index_after=RI_nc.variables['index_after'][:]
index_before=RI_nc.variables['index_before'][:]
RI_nc.close()

merged_file_dates = []
noaa_file_dates =[]
cyg_storm_names = []
cyg_files_num = []
file_counter = 0
for cyg_file in merged_cyg_files_list:
    start_ind = cyg_file.find('ddmi.')
    end_ind = cyg_file.find('.20')
    storm_name = cyg_file[start_ind+5:end_ind-3]
    cyg_nc = nc.Dataset(cyg_file)
    date_format = "%y-%m-%d %H:%M:%S"
    cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
    cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(hours=int(x)))(cyg_nc.variables['time'][:])
    dates_only = np.array([dt.date() for dt in cyg_datetime_array])
    unique_dates = np.unique(dates_only).astype('datetime64[D]')
    for i in unique_dates:
        merged_file_dates.append(i)
        cyg_files_num.append(file_counter)
        cyg_storm_names.append(storm_name)
    file_counter += 1
for cyg_file in noaa_cyg_files_list:
    start_ind = cyg_file.find('ddmi.')
    cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
    noaa_file_dates.append(cyg_datetime_object.date())

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
        return interpolated_speed

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

flag_definitions = {1: 'poor_overall_quality',2: 'asc_node',4: 'block_IIF',8: 'block_IIR',16: 'block_IIRM',
    32: 'non_zero_nst',64: 'poor_quality_wind_nst',128: 'poor_quality_wind_sample'}

results_dt = pd.DataFrame(columns=['under/over','intensityBT','durationBT','storm_idBT','tc_eventBT','wind_change','time_diff','num_tracks'])
measurement_deets_dt = pd.DataFrame(columns=['Wind_speed','WS_uncertainty','Inc','nbrcs_corrected','nbrcs_raw','Sample_flag','SNR','Num_DDM','Rx','Gain','RCG','dist','vmax_error'])
noaa_df_list = []; tc_list= []; merged_list =[]
noaa_ri_change = [] ;merged_ri_change = [] ; noaa_wind_change=[] ; merged_wind_change=[] ;noaa_ri_durations=[] ;merged_ri_durations=[] ; ri_change_failed = [] ; cyg_wind_change_failed=[]
for tc_event in range(len(storm_id)): #Hin 428 Bua 200 Kya 210 case_study_inds: #
    if ri_rd[tc_event]==1: #and storm_id[tc_event]=='2024279N21265':
    # if time[tc_event][0] >np.datetime64('2022-12-31T23:59:59'):
        # need to know the length of storm and remove nans
        nans = ~np.isnan(latitude[tc_event]).data
        tc_time = time[tc_event,nans]
        tc_storm_id = storm_id[tc_event]
        
        tc_vmax = vmax[tc_event,nans]
        tc_latitude = latitude[tc_event,nans] ; tc_longitude = longitude[tc_event,nans]
        ri_start_time = np.datetime64(tc_time[index_before[tc_event]]); ri_end_time = np.datetime64(tc_time[-index_after[tc_event]])
        ri_vmax_end = tc_vmax[-index_after[tc_event]]
        ri_vmax_start= tc_vmax[index_before[tc_event]]
        ri_vmax_change = ri_vmax_end - ri_vmax_start

        # linearly interpolate event to X minutes 
        interval = np.timedelta64(interp_interval, 'm')
        start = tc_time[0]
        end = tc_time[-1]
        tc_time_interpolated = np.arange(start, end + interval, interval)
        time_since_ri = tc_time_interpolated-ri_start_time

        # create arrays to store interpolated values
        tc_vmax_interpolated = np.zeros(np.size(tc_time_interpolated))
        tc_lat_interpolated = np.zeros(np.size(tc_time_interpolated))
        tc_lon_interpolated = np.zeros(np.size(tc_time_interpolated))
        for i in range(len(tc_time_interpolated)): # issue here with generating nans
            tc_vmax_interpolated[i] = linear_interpolation_vectorized(tc_time, tc_vmax, tc_time_interpolated[i])
            tc_lat_interpolated[i] = linear_interpolation_vectorized(tc_time, tc_latitude, tc_time_interpolated[i])
            tc_lon_interpolated[i] = linear_interpolation_vectorized(tc_time, tc_longitude, tc_time_interpolated[i])

        # interpolate tc_vmax between known values
        # Copy the array to modify
        filled = tc_vmax.copy()
        # Iterate over masked elements
        for i in range(1, len(tc_vmax) - 1):
            if tc_vmax.mask[i] and not tc_vmax.mask[i - 1] and not tc_vmax.mask[i + 1]:
                # Take average of neighbors
                filled[i] = (tc_vmax[i - 1] + tc_vmax[i + 1]) / 2
                filled.mask[i] = False  # Unmask the value
        tc_vmax=filled
        
        cyg_AOI_time_start,cyg_AOI_time_end = np.nanmin(tc_time),np.nanmax(tc_time)
        cyg_AOI_time_start = cyg_AOI_time_start.astype('datetime64[s]').astype('O')
        cyg_AOI_time_end = cyg_AOI_time_end.astype('datetime64[s]').astype('O')
        date_list = np.arange(cyg_AOI_time_start, cyg_AOI_time_end+timedelta(days=1), timedelta(days=1) ,dtype='datetime64[D]')
        
        AOI_lat_min,AOI_lat_max = np.nanmin(tc_latitude)-buffer,np.nanmax(tc_latitude)+buffer
        AOI_lon_min,AOI_lon_max = np.nanmin(tc_longitude)-buffer,np.nanmax(tc_longitude)+buffer
        
        try:
            # if cyg_version == 'noaa_1.2':
            assert tc_vmax_interpolated[time_since_ri==np.timedelta64(0,'s')] < 28 # to filter for only the low or high initial speeds

            coincident_meas= {i: np.ma.MaskedArray([]) for i in variables}
            cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_winds=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([],dtype='datetime64') ; cyg_tracks=np.ma.MaskedArray([]); cyg_dist=np.ma.MaskedArray([]); cyg_vmax_error=np.ma.MaskedArray([]) 
            cyg_lats_all=np.ma.MaskedArray([]); cyg_lons_all=np.ma.MaskedArray([]); cyg_winds_all=np.ma.MaskedArray([]); cyg_times_all=np.ma.MaskedArray([],dtype='datetime64') ; cyg_tracks_all = np.ma.MaskedArray([])
            for date in date_list:
                try:
                    cyg_file_path = noaa_cyg_files_list[np.where([date == x for x in noaa_file_dates])[0][0]]
                except IndexError:
                    raise AssertionError
                cyg_nc = nc.Dataset(cyg_file_path)
                longitudes = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ; cygnss_tracks= cyg_nc.variables['track_id'][:]
                cygnss_lons = ((longitudes + 180) % 360) - 180
                cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                        
                # find where the TC is for the date
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = np.where(lats_boolean & lons_boolean)[0] 
                if len(box_indices)>0:
                    meas_lats = cygnss_lats[box_indices]  # Measurement latitudes
                    meas_lons = cygnss_lons[box_indices]  # Measurement longitudes
                    meas_winds = cygnss_wind[box_indices] 
                    # temporal matching - time within x number of minutes
                    date_format = "%Y-%m-%d %H:%M:%S"
                    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                    cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
                    meas_times = np.array(cyg_datetime_array, dtype='datetime64')
                                            
                    # Convert to NumPy arrays if not already
                    meas_coords = np.column_stack((meas_lats, meas_lons))  # Shape (N, 2)
                    tc_coords = np.column_stack((tc_lat_interpolated, tc_lon_interpolated))  # Shape (M, 2)
                    # Compute the pairwise distance matrix
                    distances = haversine_vectorized(meas_coords, tc_coords)
                    # Compute pairwise time difference matrix (absolute time difference in hours)
                    time_diffs = np.abs(meas_times[:, None] - tc_time_interpolated[None, :])
                    # Find measurements where at least one track point is within 100 km & 30 minutes
                    valid_mask = np.any((distances <= MAX_DISTANCE) & (time_diffs <= np.timedelta64(interp_interval, 'm')), axis=1)
                    # Select valid measurement indices
                    selected_indices = np.where(valid_mask)[0]
                    # Select valid measurements
                    if len(selected_indices)>0:
                        # # Find the minimum distance and the corresponding index of the nearest track point
                        # Get the filtered lat/lon values and their corresponding nearest track indices
                        filtered_distances = np.min(distances, axis=1)[selected_indices]
                        nearest_indices = np.argmin(distances, axis=1)[selected_indices]
                        filtered_lats = meas_lats[selected_indices]
                        filtered_lons = meas_lons[selected_indices]
                        filtered_times = meas_times[selected_indices]
                        
                        # cyg_lat_AOI = filtered_lats[valid] ; cyg_lon_AOI = filtered_lons[valid] ; cyg_times_AOI = filtered_times[valid]
                        # cyg_winds_AOI = meas_winds[selected_indices][valid] ; cyg_tracks_AOI = cygnss_tracks[box_indices][selected_indices][valid] # cygnss_tracks[box_indices][indices][valid_time]
                        cyg_lat_AOI = filtered_lats ; cyg_lon_AOI = filtered_lons ; cyg_times_AOI = filtered_times
                        cyg_winds_AOI = meas_winds[selected_indices] ; cyg_tracks_AOI = cygnss_tracks[box_indices][selected_indices] # cygnss_tracks[box_indices][indices][valid_time]
    
                        cyg_lats=np.append(cyg_lats,cyg_lat_AOI) ; cyg_lons = np.append(cyg_lons,cyg_lon_AOI) ; cyg_tracks = np.append(cyg_tracks, cyg_tracks_AOI)
                        cyg_winds = np.append(cyg_winds, cyg_winds_AOI); cyg_times=np.append(cyg_times,cyg_times_AOI)
                        cyg_dist = np.append(cyg_dist,filtered_distances)
                        
                        cyg_vmax_error = tc_vmax_interpolated[nearest_indices] - cyg_winds_AOI
                        
                        cyg_lats_all=np.append(cyg_lats_all,meas_lats) ; cyg_lons_all = np.append(cyg_lons_all,meas_lons) ; cyg_tracks_all = np.append(cyg_tracks_all, cygnss_tracks[box_indices])
                        cyg_winds_all = np.append(cyg_winds_all, meas_winds); cyg_times_all=np.append(cyg_times_all,meas_times)
                        for i in variables:
                            if i == 'cyg_dist':
                                var = filtered_distances
                            elif i =='cyg_vmax_error':
                                var = cyg_vmax_error
                            elif i =='time_since_ri':
                                var = filtered_times - ri_start_time
                            else: 
                                var = cyg_nc.variables[i][box_indices[selected_indices]]
                            assert len(var) == len(cyg_lat_AOI)
                            
                            # flags = [1 = poor_overall_quality: logical OR of poor_quality_wind_nst, poor_quality_wind_sample, 2 = asc_node: set during CYGNSS ascending node, 4 = block_IIF, 8 = block_IIR, 16 = block_IIRM, 32 = non_zero_nst: the Level 1 "nst_att_status" flag (i.e. nano star tracker) was nonzero, 64 = poor_quality_wind_nst: low confidence in the reported wind speed WHILE the star tracker flag (i.e. nst_att_status) is set. 128 = poor_quality_wind_sample: set when an unrealistic wind speed sample is detected.
                            if len(coincident_meas[i]) == 0:
                                coincident_meas[i] = var
                            else:
                                coincident_meas[i]=np.concatenate((coincident_meas[i],var), axis=0)         
                cyg_nc.close()
            
            if len(cyg_times) > 0:
                # Decode all flags
                # decoded_flags_list = []
                # decoded = []
                # for flag in coincident_meas['sample_flags']:
                #     if flag & 1: decoded.append("poor_overall_quality")
                #     if flag & 2: decoded.append("asc_node")
                #     if flag & 4: decoded.append("block_IIF")
                #     if flag & 8: decoded.append("block_IIR")
                #     if flag & 16: decoded.append("block_IIRM")
                #     if flag & 32: decoded.append("non_zero_nst")
                #     if flag & 64: decoded.append("poor_quality_wind_nst")
                #     if flag & 128: decoded.append("poor_quality_wind_sample")
                #     decoded_flags_list.append(decoded)
                # coincident_meas['sample_flags'] = np.array(decoded_flags_list, dtype=object)

                tc_list.append(pd.DataFrame({'time_since_ri': time_since_ri,'wind_speed' :tc_vmax_interpolated}))
                # calculating vmax based on max velocity of each track - assuming not parralel to path
                ri_mask=((cyg_times>ri_start_time) & (cyg_times<ri_end_time))
                ri_lats=cyg_lats[ri_mask]; ri_lons=cyg_lons[ri_mask]; ri_winds=cyg_winds[ri_mask]; ri_times=cyg_times[ri_mask] ; ri_tracks=cyg_tracks[ri_mask]; ri_dist = cyg_dist[ri_mask]     
                mid_time = ri_start_time + (ri_end_time - ri_start_time) / 2  # middle of the interpolated TC times
                noaa_df_list.append(pd.DataFrame(coincident_meas))

                unique_tracks = np.unique(ri_tracks)
                num_tracks = len(unique_tracks)
                if num_tracks > 1: # so 2 or more pass over
                    # cyg_lat_AOI = selected_lats[valid_time] ; cyg_lon_AOI = selected_lons[valid_time] ; cyg_times_AOI = cyg_datetime_array[valid_time]
                    # cyg_winds_AOI =cygnss_wind[box_indices][indices][valid_time]                
                    # Arrays to store results
                    vmax_time = np.zeros(num_tracks, dtype=object)
                    max_wind = np.zeros(num_tracks)
                    
                    # Process each track - breaking here
                    for i, track in enumerate(unique_tracks):
                        track_mask = (ri_tracks == track)  # Proper boolean mask for track filtering
                        track_times = ri_times[track_mask]  # Filter timestamps based on track mask
                        track_winds = ri_winds[track_mask]  # Filter wind speeds based on track mask
                    
                        # Ensure that we are not running into errors if the filtered array is empty
                        if track_times.size > 0:                         
                            # Store the values
                            vmax_index = np.argmax(track_winds)
                            vmax_time[i] = track_times[vmax_index]
                            max_wind[i] = track_winds[vmax_index]  # Maximum wind speed for the track
                                        
                    earliest_ind = np.argmin(vmax_time) ; latest_ind = np.argmax(vmax_time)
                    wind_change = np.max(max_wind) - max_wind[earliest_ind]
                    time_diff  = (vmax_time[np.argmax(max_wind)] - vmax_time[earliest_ind]).astype('timedelta64[h]')/np.timedelta64(1,'h')
                    if time_diff>=durations[tc_event]/2:
                        noaa_ri_change.append(ri_vmax_change)
                        noaa_wind_change.append(wind_change)
                        noaa_ri_durations.append(durations[tc_event])
                    else:
                        ri_change_failed.append(ri_vmax_change)
                        cyg_wind_change_failed.append(wind_change)
                    
                    # print(intensities[tc_event],durations[tc_event],wind_change,time_diff, cyg_AOI_time_start, cyg_AOI_time_end, cyg_lats[0],cyg_lons[0],ri_rd[tc_event])    
                    # now to filter for meeting criteria
                    plotting = False
                    if wind_change < intensities[tc_event]/2 and time_diff>=durations[tc_event]/2-2 and num_tracks>=durations[tc_event]/3:
                        new_row = ['under',intensities[tc_event],durations[tc_event],storm_id[tc_event], tc_event,int(wind_change),int(time_diff),num_tracks]
                        new_df = pd.DataFrame([new_row], columns=results_dt.columns)
                        results_dt = pd.concat([results_dt,new_df])
                        plotting=True
                    elif wind_change >= intensities[tc_event] and time_diff>=durations[tc_event]/2:# and num_tracks>=durations[tc_event]/3:
                        new_row = ['over',intensities[tc_event],durations[tc_event],storm_id[tc_event], tc_event,int(wind_change),int(time_diff),num_tracks]
                        new_df = pd.DataFrame([new_row], columns=results_dt.columns)
                        results_dt = pd.concat([results_dt,new_df])
                        plotting = True
                    if False:
                        dlat = np.diff(tc_lat_interpolated, prepend=tc_lat_interpolated[0])
                        dlon = np.diff(tc_lon_interpolated, prepend=tc_lon_interpolated[0])
                        
                        # Create plot
                        import matplotlib.pyplot as plt
                        import numpy as np
                        import matplotlib.animation as animation
                        import matplotlib.cm as cm
                        import matplotlib.colors as mcolors
                        buffer=1.1
                        latitudes = tc_lat_interpolated ; longitudes = tc_lon_interpolated; times = tc_time_interpolated
                        winds = tc_vmax_interpolated
                        fig, ax = plt.subplots(figsize=(8, 6))
    
                        # Set plot limits
                        ax.set_xlim(min(longitudes) - buffer, max(longitudes) + buffer)
                        ax.set_ylim(min(latitudes) - buffer, max(latitudes) + buffer)                                        
                        
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
                        
                            # Plot cyclone track up to current frame
                            # ax.plot(longitudes[:frame + 1], latitudes[:frame + 1], color='gray', linestyle='dashed')
                        
                            # Plot arrows for wind direction with color mapped to wind speed
                            quiver = ax.quiver(longitudes[:frame + 1], latitudes[:frame + 1], 
                                                dlon[:frame + 1], dlat[:frame + 1], 
                                                winds[:frame + 1],  # Use wind speed for color
                                                cmap='coolwarm', width=0.005, norm=norm, angles='xy', scale_units='xy')
                            
                            mask = (cyg_times_all <= times[frame]) & (cyg_times_all >= times[frame] - np.timedelta64(8,'h'))
                            ax.scatter(cyg_lons_all[mask],cyg_lats_all[mask],c='grey', marker='x',s=20)
                            # Add random measurements up to the current frame
                            mask = (cyg_times <= times[frame]) & (cyg_times >= times[frame] - np.timedelta64(8,'h'))
                            ax.scatter(cyg_lons[mask], cyg_lats[mask], 
                                       c=cyg_winds[mask], cmap=cmap, norm=norm, edgecolor='black', s=60)
        
                            # Update time label
                            if times[frame] < ri_start_time:
                                phase='Before RI'
                            elif times[frame] < mid_time:
                                phase='Early RI'
                            elif times[frame] <= ri_end_time:
                                  phase='Later RI'
                            elif times[frame] > ri_end_time:
                                  phase='After RI'
                            else:
                                phase='what happened?'
                            ax.text(0.98, 0.95, f"Phase: {phase}, Time: {times[frame]}", 
                             transform=plt.gca().transAxes,  # Places text relative to axes
                             fontsize=12, ha='right', color='black', bbox=dict(facecolor='white', alpha=0.6))
         
                            # time_text = ax.text(0.1, 0.90, '', fontsize=12)
                            # time_text.set_text(f'Time: {times[frame]} hours')
                        
                            # Titles & labels
                            ax.set_title('ID: {0},Duration: {1}, Intensity: {2}'.format(storm_id[tc_event],durations[tc_event],intensities[tc_event]))
                            ax.set_xlabel('Longitude')
                            ax.set_ylabel('Latitude')
        
                        # Create animation
                        pause_frames = 20 
                        ani = animation.FuncAnimation(fig, update_frame, frames=list(range(len(times)))+[len(times) - 1] * pause_frames, repeat=True)
                        
                        # Save GIF
                        ani.save(r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\gifs\NOAA_{1}_{2}.gif'.format(comp, tc_event, storm_id[tc_event]), fps=5)
                        
                        # Show plot
                        plt.show()
                        
                        # Time series plot
                        plot_timeline_refined(cyg_winds,cyg_times,cyg_dist,tc_time_interpolated,tc_vmax_interpolated,storm_id[tc_event]+'_'+str(tc_event),ri_start_time,ri_end_time)
                
            assert np.any(np.isin(date_list,  np.array(merged_file_dates)))
            # if cyg_version == 'storm_centric':
            valid_files = np.where(date_list[0] == np.array(merged_file_dates))[0]
            unique_storm_files=[]
            if len(valid_files)==1:
                unique_storm_files.append(merged_cyg_files_list[cyg_files_num[valid_files[0]]])
            else:
                for i in range(len(valid_files)-1):
                    if merged_cyg_files_list[cyg_files_num[valid_files[i]]] != merged_cyg_files_list[cyg_files_num[valid_files[i+1]]]:
                        unique_storm_files.append(merged_cyg_files_list[cyg_files_num[valid_files[i]]])
                    
            for cyg_file_path in unique_storm_files:
                cyg_nc = nc.Dataset(cyg_file_path)
                cyg_storm_lat = cyg_nc.variables['best_track_storm_center_lat'][:] ; cyg_storm_lon = ((cyg_nc.variables['best_track_storm_center_lon'][:] + 180) % 360) - 180
                assert True in (cyg_storm_lon > AOI_lon_min) & True in (cyg_storm_lon < AOI_lon_max)
                assert True in (cyg_storm_lat > AOI_lat_min) & True in (cyg_storm_lat < AOI_lat_max)
               
                epochs=np.where((cyg_storm_lon > AOI_lon_min)&(cyg_storm_lon < AOI_lon_max) & (cyg_storm_lat > AOI_lat_min) & (cyg_storm_lat < AOI_lat_max))[0]
                if len(epochs)>1: # to ensure we have two or more measurements over RI period
                    vmax_lat = np.round(cyg_nc.variables['cygnss_vmax_lat'][epochs],2) ; vmax_lon = np.round(((cyg_nc.variables['cygnss_vmax_lon'][epochs]+ 180) % 360) - 180,2)
                    assert True in (vmax_lon > AOI_lon_min) & True in (vmax_lon < AOI_lon_max)
                    assert True in (vmax_lat > AOI_lat_min) & True in (vmax_lat < AOI_lat_max)
                   
                    longitudes  = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][epochs] ;
                    cygnss_lons = np.round(((longitudes + 180) % 360) - 180,2)
                    cygnss_time =cyg_nc.variables['time'][epochs] ; cygnss_time_offset = cyg_nc.variables['time_offset'][epochs]
                    cygnss_time_units =cyg_nc.variables['time'].units[12:]
                    cyg_datetimes = np.datetime64(cygnss_time_units) + cygnss_time.astype('timedelta64[h]')
                    start_epoch = np.where(cyg_datetimes==ri_start_time)[0]
                    end_epoch = np.where(cyg_datetimes==ri_end_time)[0]
                    assert len(end_epoch)>0 and len(start_epoch)>0
                    cyg_vmax=np.zeros(len(epochs))
                    for j in range(len(epochs)):
                        cyg_vmax[j]= cygnss_wind[j,np.where(cygnss_lats == vmax_lat[j]),np.where(cygnss_lons == vmax_lon[j])]
                    change_vmax = cyg_vmax[end_epoch] - cyg_vmax[start_epoch]
                    # cyg_vmax_start = cygnss_wind[start_epoch,np.where(cygnss_lats == vmax_lat[start_epoch]),np.where(cygnss_lons == vmax_lon[start_epoch])]
                    # cyg_vmax_end = cygnss_wind[end_epoch,np.where(cygnss_lats == vmax_lat[end_epoch]),np.where(cygnss_lons == vmax_lon[end_epoch])]
                    # change_vmax = cyg_vmax_end - cyg_vmax_start 
                    merged_ri_change.append(ri_vmax_change)
                    merged_wind_change.append(change_vmax)
                    merged_ri_durations.append(durations[tc_event])
                    merged_list.append(pd.DataFrame({'time_since_ri': cyg_datetimes-ri_start_time,'wind_speed' :cyg_vmax}))
                    tc_list.append(pd.DataFrame({'time_since_ri': time_since_ri,'wind_speed' :tc_vmax_interpolated}))
                        
                    # if change_vmax >= intensities[tc_event]:
                    #     print(int(change_vmax),intensities[tc_event],durations[tc_event],storm_id[tc_event], tc_event,cyg_file_path[39:46], vmax_lat[0],vmax_lon[0],cyg_AOI_time_start)
            
        except AssertionError:
            print('assertion error - continuing')    
            pass
scatterplot=True
if scatterplot:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    ri_change = np.array(merged_ri_change) ; cyg_wind_change = np.array(merged_wind_change).flatten()
    lower = 15.3 ; upper = 23.1
    in_box      = ((ri_change >= lower) & (ri_change <= upper)) & ((cyg_wind_change >= lower) & (cyg_wind_change <= upper))
    x_high      = (ri_change > upper) & ((cyg_wind_change >= lower) & (cyg_wind_change <= upper))
    y_high      = ((ri_change >= lower) & (ri_change <= upper)) & (cyg_wind_change > upper)
    xy_high     = (ri_change > upper) & (cyg_wind_change > upper)
    x_low       = (ri_change < lower) & ((cyg_wind_change >= lower) & (cyg_wind_change <= upper))
    y_low       = ((ri_change >= lower) & (ri_change <= upper)) & (cyg_wind_change < lower)
    y_low_x_high = ((ri_change >= upper) & (cyg_wind_change < lower))
    xy_low      = (ri_change < lower) & (cyg_wind_change < lower)
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(ri_change, cyg_wind_change)
    plt.axvline(x=lower, color='grey', linestyle='--')
    plt.axvline(x=upper, color='grey', linestyle='--')
    plt.axhline(y=lower, color='grey', linestyle='--')
    plt.axhline(y=upper, color='grey', linestyle='--')
    plt.xlim(14,42) ; plt.ylim(-11, 39)
    # Add middle green box (in_box)
    plt.gca().add_patch(patches.Rectangle((15.3, 15.3), 23.1-15.3, 23.1-15.3,linewidth=0, facecolor='green', alpha=0.2))
    plt.gca().add_patch(patches.Rectangle((23.1, 23.1), plt.xlim()[1]-23.1, plt.ylim()[1]-23.1,linewidth=0, facecolor='green', alpha=0.2))
    total_15 = np.sum(in_box) + np.sum(y_high) + np.sum(y_low) ; total_23 = np.sum(y_low_x_high) + np.sum(xy_high) + np.sum(x_high)
    plt.text(16, 16, f"{np.sum(in_box)} pts  {np.sum(in_box)/total_15*100:.1f}%", fontsize=12, color='black')      # Middle box
    plt.text(23.5, 16, f"{np.sum(x_high)} pts  {np.sum(x_high)/total_23*100:.1f}%", fontsize=12, color='black')       # Right
    plt.text(16, 25.5, f"{np.sum(y_high)} pts  {np.sum(y_high)/total_15*100:.1f}%", fontsize=12, color='black')       # Above
    plt.text(23.5, 25.5, f"{np.sum(xy_high)} pts  {np.sum(xy_high)/total_23*100:.1f}%", fontsize=12, color='black')    # Top-right
    plt.text(16, 8, f"{np.sum(y_low)} pts  {np.sum(y_low)/total_15*100:.1f}%", fontsize=12, color='black')           # Below
    plt.text(23.5, 8, f"{np.sum(y_low_x_high)} pts  {np.sum(y_low_x_high)/total_23*100:.1f}%", fontsize=12, color='black')  # Lower-Right
    plt.axvline(x=15.3, color='grey', linestyle='--')
    plt.axvline(x=23.1, color='grey', linestyle='--')
    plt.axhline(y=15.3, color='grey', linestyle='--')
    plt.axhline(y=23.1, color='grey', linestyle='--')
    plt.xlabel('Best Track Vmax change (m/s)')
    plt.ylabel('CYGNSS Vmax change (m/s)')
    plt.title('merged')
    plt.savefig(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Plots\skill_merged_low_wind_speed.png' , dpi=500)
    plt.show()

results_dt.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\ri_results_v5.csv')
final_df = pd.concat(noaa_df_list, ignore_index=True)
final_tc_df = pd.concat(tc_list, ignore_index=True)
final_merged_df =pd.concat(merged_list, ignore_index=True)


final_df.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\NOAA_Measurements.csv')
final_merged_df.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\Merged_Measurements.csv')


plot_distribution(final_df,final_tc_df,final_merged_df)


bootstrapping=False
if bootstrapping:
    from sklearn.ensemble import RandomForestRegressor
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    # Prepare data
    df = final_df.select_dtypes(include=['number'])  # Keep only numeric columns
    df = df.dropna()

    y = df['cyg_vmax_error']
    X = df.drop(columns=['cyg_vmax_error', 'wind_speed', 'wind_speed_uncertainty', 'sample_flags'])
    X.columns = ['Incidence angle','Receiver gain','SNR','Range-corrected gain','# DDMs used',
                 'Avg.NBRCS corrected','Avg.NBRCS original','Distance to best-track','RI times']
    feature_names = X.columns
    # Parameters
    n_bootstrap = 500
    n_features = X.shape[1]
    importances_bootstrap = np.zeros((n_bootstrap, n_features))
    
    # Bootstrap loop
    for i in range(n_bootstrap):
        sample_indices = np.random.choice(len(X), size=len(X), replace=True)
        X_sample = X.iloc[sample_indices]
        y_sample = y.iloc[sample_indices]
        
        model = RandomForestRegressor(n_estimators=100, random_state=i)
        model.fit(X_sample, y_sample)
        importances_bootstrap[i, :] = model.feature_importances_
    
    # Compute mean and 95% confidence interval
    importance_mean = np.mean(importances_bootstrap, axis=0)
    importance_lower = np.percentile(importances_bootstrap, 5, axis=0)
    importance_upper = np.percentile(importances_bootstrap, 95, axis=0)
    
    # Create a DataFrame for plotting
    importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Mean Importance': importance_mean,
        'Lower CI': importance_lower,
        'Upper CI': importance_upper
    }).sort_values(by='Mean Importance', ascending=False)
    
    # Plot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=importance_df, x='Mean Importance', y='Feature', palette='viridis')
    plt.errorbar(importance_df['Mean Importance'], importance_df['Feature'],
                 xerr=[importance_df['Mean Importance'] - importance_df['Lower CI'],
                       importance_df['Upper CI'] - importance_df['Mean Importance']],
                 fmt='none', ecolor='black', capsize=3, label='95% CI')
    plt.xlabel('Bootstrapped Feature Importance')
    plt.tight_layout()
    plt.show()
                        # # extract the vmax at each epoch
                        # cyg_vmax=np.zeros(len(epochs))
                        # for j in range(len(epochs)):
                        #     cyg_vmax[j]= cygnss_wind[j,np.where(cygnss_lats == vmax_lat[j]),np.where(cygnss_lons == vmax_lon[j])]
                        # vmax_besttrack = cyg_nc.variables['best_track_vmax'][epochs]
                        # cyg_nc.close()
                        # diff = cyg_vmax - vmax_besttrack
                        # print(diff)
                        # start_ind = 0 ; end_ind = 1
                        # while end_ind < len(epochs):
                        #     if cygnss_time[end_ind] - cygnss_time[start_ind] < durations[tc_event]:
                        #         end_ind += 1
                        #     elif  cygnss_time[end_ind] - cygnss_time[start_ind] > durations[tc_event]:
                        #         start_ind += 1
                        #     elif cygnss_time[end_ind] - cygnss_time[start_ind] == durations[tc_event]:
                        #         assert cygnss_time[start_ind]
                        #         change_vmax = cyg_vmax[end_ind] - cyg_vmax[start_ind]
                        #         # change_vmax_besttrack = vmax_besttrack[end_ind] - vmax_besttrack[start_ind]
                        #         if change_vmax > intensities[tc_event]: # and change_vmax_besttrack > intensities[tc_event]:
                        #             print(int(change_vmax),intensities[tc_event],durations[tc_event],storm_id[tc_event], tc_event,cyg_storm_names[valid_file_ind], vmax_lat[0],vmax_lon[0],cyg_AOI_time_start,valid_file_ind)
                        #         start_ind += 1
                        #         end_ind += 1
                        #     else:
                        #         end_ind += 1

                            # new_row = ['over',intensities[tc_event],durations[tc_event],storm_id[tc_event], tc_event,int(change_vmax),int(time_diff),'NA']
                            # new_df = pd.DataFrame([new_row], columns=results_dt.columns)
                            # results_dt = pd.concat([results_dt,new_df])
            
              
        # plotting = False
        # if plotting:
        #     import matplotlib.pyplot as plt
        #     from cartopy import crs as ccrs 
        #     import matplotlib.patches as mpatches
        
        #     # Create a Cartopy PlateCarree projection (cylindrical projection)
        #     projection = ccrs.PlateCarree()
            
        #     # Create a Matplotlib figure and axis
        #     fig, ax = plt.subplots(subplot_kw={'projection': projection})
        #     vmin = 0 ; vmax = 25
        #     im = ax.scatter(cygnss_lons_AOI, cygnss_lats_AOI, s=cygnss_wind_AOI , c=cygnss_time_AOI, cmap='BuGn',transform=projection, marker='^',edgecolors = 'black',vmin=vmin, vmax=vmax)            
        #     if ri_rd[i] == 1:
        #         bbox = mpatches.Rectangle((longitude[i,0], latitude[i,0]), longitude[i,0] - longitude[i,1], latitude[i,0] - latitude[i,1],
        #                                   linewidth=1, edgecolor='red', facecolor='none', transform=ccrs.PlateCarree())
        #     else:
        #         bbox = mpatches.Rectangle((longitude[i,0], latitude[i,0]), longitude[i,0] - longitude[i,1], latitude[i,0] - latitude[i,1],
        #                                   linewidth=1, edgecolor='blue', facecolor='none', transform=ccrs.PlateCarree())
        #     ax.add_patch(bbox)
        #     # im = ax.scatter(cygnss_lons_AOI, cygnss_lats_AOI, c=cygnss_YSLF_AOI, cmap='viridis',transform=projection, vmin=vmin, vmax=vmax, marker='1')
        #     cbar = plt.colorbar(im, ax=ax, orientation='vertical',location='left', shrink=0.8)
        #     cbar.set_label('Time')
        
        #     # Add coastlines and gridlines for better context
        #     ax.coastlines()
        #     buffer = 0.07
        #     ax.set_xlim(AOI_lon_min-buffer,AOI_lon_max+buffer)  # Set x-axis limits
        #     ax.set_ylim(AOI_lat_min-buffer,AOI_lat_max+buffer)  # Set y-axis limits
        #     plt.xlabel('Longitude')
        #     plt.ylabel('Latitude')
        #     gls = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        #     gls.xlocator = plt.FixedLocator(range(-180, 181, 5))  # Longitude lines every 10 degrees
        #     gls.ylocator = plt.FixedLocator(range(-90, 91, 5))   # Latitude lines every 10 degrees
        #     gls.top_labels=False   # suppress top labels
        #     gls.right_labels=False # suppress right labels
        #     ax.set_title(f'CYGNSS {str(cyg_AOI_time_start)[0:10]}-{str(cyg_AOI_time_end)[11:16]}', loc='right')
            
        #     plt.subplots_adjust(left=0.55, right=0.65, bottom=0.1, top=0.9, wspace=0.5)
        #     plt.tight_layout()
            
        #     # # saving file
        #     directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\RI{1}{2}.png'.format(comp,i,cyg_version)
        #     fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        #     plt.show()
            
        
        # old distance/measurement aproaches
                        # measurement_lats = cygnss_lats[box_indices] ; measurement_lons = cygnss_lons[box_indices]
                        # EARTH_RADIUS = 6371.0  # Earth's radius in kilometers               
                        # # Function to convert lat/lon to Cartesian coordinates
                        # def latlon_to_cartesian(lat, lon):
                        #     """Convert latitude and longitude to Cartesian coordinates (x, y, z)."""
                        #     lat, lon = np.radians(lat), np.radians(lon)  # Convert to radians
                        #     x = EARTH_RADIUS * np.cos(lat) * np.cos(lon)
                        #     y = EARTH_RADIUS * np.cos(lat) * np.sin(lon)
                        #     z = EARTH_RADIUS * np.sin(lat)
                        #     return np.array([x, y, z]).T
                        # # Convert measurement points and track points to Cartesian coordinates
                        # measurement_cartesian = latlon_to_cartesian(measurement_lats, measurement_lons)
                        # track_cartesian = latlon_to_cartesian(tc_lat_interpolated, tc_lon_interpolated)
                        # # Build the KD-tree for measurement points
                        # tree = cKDTree(measurement_cartesian)
                        # # Query the tree for points within MAX_DISTANCE of each track point
                        # indices = set()  # Use a set to avoid duplicates
                        # for track_point in track_cartesian:
                        #     nearby_indices = tree.query_ball_point(track_point, MAX_DISTANCE)
                        #     indices.update(nearby_indices)
                            
                        # if len(indices)>0:
                        #     # Filter the measurement points based on the indices
                        #     indices = np.array(list(indices))
                        #     selected_lats = measurement_lats[indices]
                        #     selected_lons = measurement_lons[indices]
                            
                        #     # temporal matching - time within x number of minutes
                        #     cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                        #     date_format = "%Y-%m-%d %H:%M:%S"
                        #     cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                        #     cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices][indices])
                        #     measurement_times = np.array(cyg_datetime_array, dtype='datetime64')
                        #     time_diffs = np.abs(measurement_times - tc_time_interpolated)  # Shape (num_measurements, num_track_times)
                                                
                        #     # Find the minimum time difference for each measurement
                        #     min_time_diff = np.min(time_diffs, axis=1)
                            
                        #     # Define threshold of 30 minutes
                        #     threshold = np.timedelta64(interp_interval, 'm')
                            
                        #     # Find measurements within 30 minutes of any track time
                        #     valid_time = min_time_diff <= threshold
        
        
                        # # Convert lat/lon to Cartesian coordinates (approximate)
                        # def latlon_to_xyz(lat, lon):
                        #     R = 6371.0  # Earth radius in km
                        #     lat, lon = np.radians(lat), np.radians(lon)
                        #     x = R * np.cos(lat) * np.cos(lon)
                        #     y = R * np.cos(lat) * np.sin(lon)
                        #     z = R * np.sin(lat)
                        #     return np.vstack((x, y, z)).T
                        # # Build KDTree using track points
                        # track_xyz = latlon_to_xyz(tc_lat_interpolated, tc_lon_interpolated)
                        # tree = cKDTree(track_xyz)
                        # # Query closest track point for each measurement - issues here
                        # meas_xyz = latlon_to_xyz(meas_lats, meas_lons)
                        # distances, indices = tree.query(meas_xyz, distance_upper_bound=1) # this is broken
                        # inf_ind = np.where(distances != np.inf)[0]
                        # indices = indices[inf_ind] ; meas_times = meas_times[inf_ind]; distances=distances[inf_ind] 
                        # meas_lats = meas_lats[inf_ind] ;meas_lons = meas_lons[inf_ind] ; meas_winds = meas_winds[inf_ind]
                        
                        # Filter based on distance and time difference
                        # valid = np.abs(filtered_times - tc_time_interpolated[nearest_indices]) <= np.timedelta64(interp_interval, 'm')
                        
        
            #     lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
            #     lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
            #     box_indices = [0,0]
                
            # else:
            #     cyg_file_path = cyg_files_list[np.where([date == x for x in cyg_file_dates])[0][0]]
            #     cyg_nc = nc.Dataset(cyg_file_path)
            #     longitudes = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
            #     cygnss_lons = ((longitudes + 180) % 360) - 180
            #     cygnss_qual = cyg_nc.variables['sample_flags']
            #     cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
            #     lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
            #     lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
            #     box_indices = np.where(lats_boolean & lons_boolean)[0]
                
            
            
                        # # Haversine formula to calculate great-circle distance (in km)
                        # def haversine(lat1, lon1, lat2, lon2):
                        #     if None in (lat1, lon1, lat2, lon2):  # Skip if any value is None
                        #         return np.nan  # Return NaN instead of failing
                        #     R = 6371.0  # Earth radius in km
                        #     lat1, lon1, lat2, lon2 = np.radians([lat1, lon1, lat2, lon2])
                        #     dlat = lat2 - lat1
                        #     dlon = lon2 - lon1
                        #     a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
                        #     c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
                        #     return R * c  # Distance in km
                        
                        # # Compute pairwise distance matrix
                        # distances = np.array([[haversine(lat1, lon1, lat2, lon2) 
                        #                        for lat2, lon2 in zip(tc_lat_interpolated, tc_lon_interpolated)] 
                        #                       for lat1, lon1 in zip(meas_lats, meas_lons)])
                        
        
            # if len(box_indices)>0:
            #     if cyg_version == 'storm_centric':
            #         # find closest epoch to SAR time
            #         date_format = "%y-%m-%d %H:%M:%S"
            #         cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
            #         cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
            #         cyg_measured_times = np.ma.MaskedArray(np.zeros(np.shape(cygnss_time_offset)))
            #         for epoch in range(len(cygnss_time)):
            #             cyg_measured_times[epoch] = cygnss_time_offset[epoch] + cygnss_time[epoch]
            #         # Focus on the lat and lons
            #         try:
            #             timedeltas = np.vectorize(lambda x: timedelta(hours=x))(cyg_measured_times[:,lats_boolean,:][:,:,lons_boolean])
            #         except ValueError:
            #             raise AssertionError
            #         # Add the base datetime to each timedelta to get the 3D array of datetime objects
            #         cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + x)(timedeltas)                        
                    
            #     else:
            #         cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
            #         date_format = "%Y-%m-%d %H:%M:%S"
            #         cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
            #         cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
                
            #     time_boolean = (cyg_datetime_array > cyg_AOI_time_start) & (cyg_datetime_array < cyg_AOI_time_end)
            #     true_indices = np.where(time_boolean)[0]
            #     if len(true_indices)>0:
            #         if cyg_version == 'storm_centric':
            #             cygnss_wind_AOI = cygnss_wind[:,lats_boolean,:][:,:,lons_boolean][time_boolean] ; cygnss_time_AOI = cyg_datetime_array[time_boolean]
            #             cygnss_lats_AOI = cygnss_lats[lats_boolean][np.where(time_boolean)[1]] ; cygnss_lons_AOI = cygnss_lons[lons_boolean][np.where(time_boolean)[2]]
            #         else:
            #             cygnss_wind_AOI = cygnss_wind[box_indices[true_indices]] ; cygnss_time_AOI = cyg_datetime_array[true_indices] ; 
            #             #   sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
            #             cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]
            #             cyg_qual_AOI = cygnss_qual[box_indices[true_indices]]
            #         cyg_AOI = True                     
            #         cyg_lats=np.append(cyg_lats,cygnss_lats_AOI) ; cyg_lons = np.append(cyg_lons,cygnss_lons_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI); cyg_times=np.append(cyg_times,cygnss_time_AOI)
          
            
          # # Compute the change in wind speed - issue here, not just first and last
          # # Sort by time
          # sorted_indices = np.argsort(vmax_time)
          # sorted_times = vmax_time[sorted_indices]
          # sorted_wind_speeds = max_wind[sorted_indices]
          
          # # Initialize mask for significant wind changes
          # significant_change_mask = np.zeros_like(sorted_times, dtype=bool)
          
          # # Compare each measurement with future points within 24 hours
          # for i in range(len(sorted_times)):
          #     for j in range(i + 1, len(sorted_times)):
          #         # Stop checking if time difference is greater than 24 hours
          #         time_diff = (sorted_times[j] - sorted_times[i]).astype('timedelta64[m]')/np.timedelta64(1,'m')
          #         # if  time_diff < durations[tc_event]*60/2 - 60: # custon threshold of minimal separation
          #         #     break
                  
          #         # Check wind speed difference
          #         wind_change = sorted_wind_speeds[j] - sorted_wind_speeds[i]
          #         if wind_change > 12:
          #         # if wind_change > intensities[tc_event]:
          #             significant_change_mask[i] = True
          #             significant_change_mask[j] = True
          
          # # Select measurements with significant wind changes
          # changing_times = sorted_times[significant_change_mask]
          # changing_wind_speeds = sorted_wind_speeds[significant_change_mask]
          
          # if True in significant_change_mask:
                        