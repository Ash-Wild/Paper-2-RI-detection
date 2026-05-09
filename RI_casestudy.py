# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 08:43:49 2025

@author: ashle
"""
#%%
import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import glob
import os
import pandas as pd
import cfgrib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from cyg_L1_downloader import cyg_l1_downloader
from RI_casestudy_insitu_finder import haversine_vectorized

comp = 'S3987712' # 'S3987712' 'ashle'
drive = {'S3987712':'D:', 'ashle':'E:'}
event_list = [['Helene_2024',665]] # list of RI events from RI ,   
# hinnamnor v5 '2022239N22150' 428
# ['Hinnamnor_2022',428], ['Helene_2024',665]
# kyarr 2019296N15066 210
# Bualoi 2019290N08169 200
# Ian 2022266N12294 
interp_interval= 20 # minutes to interpolate best-track
MAX_DISTANCE = 80  # Maximum distance in kilometers
buffer = MAX_DISTANCE/100 # degrees buffer to search in storm_centric

variables = ['eff_scatter', 'brcs','sp_inc_angle', 'ddm_snr','gps_eirp', 'rx_to_sp_range', 'sc_vel_x_pvt','sc_vel_y_pvt','sc_vel_z_pvt','sp_theta_orbit','tx_to_sp_range','quality_flags','sv_num', 'ddm_ant','nst_att_status', 'ddm_nbrcs',]
counter = 0
num_valid = 0

# Loop through each case study and process the data
for case_study in event_list:
    RI_file = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Data', 'IBTrACS', 'RI_events_v5.nc')
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
    # Vectorized interpolation using seconds-since-epoch (faster than looping)
    tc_times_s = tc_time.astype('datetime64[s]').astype('int64')
    interp_s = tc_time_interpolated.astype('datetime64[s]').astype('int64')

    # Ensure numeric arrays (fill masked vmax with nan)
    tc_vmax_vals = np.ma.filled(tc_vmax, np.nan).astype(float)
    tc_lat_vals = np.asarray(tc_latitude).astype(float)
    tc_lon_vals = np.asarray(tc_longitude).astype(float)

    # Filter out NaN values from vmax for interpolation
    valid_idx = ~np.isnan(tc_vmax_vals)
    if np.sum(valid_idx) > 1:  # Need at least 2 valid points to interpolate
        tc_vmax_interpolated = np.interp(interp_s, tc_times_s[valid_idx], tc_vmax_vals[valid_idx])
    else:
        # If not enough valid vmax values, use all values and let np.interp handle NaNs
        tc_vmax_interpolated = np.interp(interp_s, tc_times_s, tc_vmax_vals)
    
    tc_lat_interpolated = np.interp(interp_s, tc_times_s, tc_lat_vals)
    tc_lon_interpolated = np.interp(interp_s, tc_times_s, tc_lon_vals)
    tc_coords = np.column_stack((tc_lat_interpolated, tc_lon_interpolated))  # Shape (M, 2)

    cyg_AOI_time_start,cyg_AOI_time_end = np.nanmin(tc_time),np.nanmax(tc_time)
    cyg_AOI_time_start = cyg_AOI_time_start.astype('datetime64[s]').astype('O')
    cyg_AOI_time_end = cyg_AOI_time_end.astype('datetime64[s]').astype('O')
    AOI_lat_min,AOI_lat_max = np.nanmin(tc_latitude)-buffer,np.nanmax(tc_latitude)+buffer
    AOI_lon_min,AOI_lon_max = np.nanmin(tc_longitude)-buffer,np.nanmax(tc_longitude)+buffer

    # cyg 
    case_folder = rf'{drive[comp]}:\Phd_data\cyg_l1\{case_study[0]}'
    try:
        cyg_files_list = np.asarray(glob.glob(case_folder+ r"\cyg*.nc"))
        assert len(cyg_files_list)>0
    except AssertionError:
        cyg_l1_downloader(AOI_lon_min,AOI_lon_max,AOI_lat_min,AOI_lat_max,cyg_AOI_time_start,cyg_AOI_time_end,case_study[0])
        cyg_files_list = np.asarray(glob.glob(case_folder+ r"\cyg*.nc"))
        assert len(cyg_files_list)>0, "No CYGNSS files found after download attempt. Please check the downloader function and file paths."
    
    coincident_meas= {i: np.ma.MaskedArray([]) for i in variables}
    cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([],dtype='datetime64') ; nearest_tc_coords=np.ma.MaskedArray([]);
    cyg_lats_all=np.ma.MaskedArray([]); cyg_lons_all=np.ma.MaskedArray([]); cyg_times_all=np.ma.MaskedArray([],dtype='datetime64') 
    unreadable_cyg_files = []
    num_cyg_files_loaded = 0
    for cyg_file_path in cyg_files_list:
        try:
            cyg_nc = nc.Dataset(cyg_file_path)
            num_cyg_files_loaded += 1
        except OSError as e:
            unreadable_cyg_files.append(cyg_file_path)
            print(f"Skipping unreadable CYGNSS file: {cyg_file_path} ({e})")
            continue
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
                    var = None
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
                            coincident_meas[i] = np.ma.concatenate((coincident_meas[i], var[box_indices[selected_indices]]), axis=0)
                    except AssertionError:
                        if var is not None:
                            var=np.repeat(var,4)
                            assert len(var) == len(cygnss_lons)
                            coincident_meas[i]=np.ma.concatenate((coincident_meas[i],var[box_indices][selected_indices]), axis=0) 
                        
                                # Group measurement-TC pairs and compute average TC index per measurement  
                meas_indices = valid_ind[0]  # measurement indices (x)
                tc_indices = valid_ind[1]    # TC track indices (y)
                
                # Use pandas groupby for efficient aggregation
                groupby_df = pd.DataFrame({'meas_idx': meas_indices, 'tc_idx': tc_indices})
                grouped = groupby_df.groupby('meas_idx')['tc_idx'].mean().astype(int)
                
                cyg_valid_meas_inds = np.asarray(grouped.index.values)  # unique measurement indices
                nearest_tc_coord_ind = np.asarray(grouped.values)        # averaged TC indices
                
                if len(nearest_tc_coord_ind) != len(meas_times[selected_indices]):
                    print('error')
                nearest_tc_coords = np.concatenate((nearest_tc_coords, nearest_tc_coord_ind))
        cyg_nc.close()

    if num_cyg_files_loaded == 0:
        raise OSError(
            f"No readable CYGNSS files found in {case_folder}. "
            f"Unreadable files: {len(unreadable_cyg_files)}"
        )

    # environmental data
    cms_buoy_data = None  # Initialize
    if case_study[0] == 'Helene_2024':
        ndbc_buoy_data = pd.read_csv(rf'{drive[comp]}\Phd_data\casestudy\Buoy\all_vars\\{tc_storm_id}_ndbc_insitu.csv')
        cms_buoy_data = pd.read_csv(rf'{drive[comp]}\Phd_data\casestudy\Buoy\all_vars\\{tc_storm_id}_{case_study[1]}_cms_insitu.csv')
    
    # if cms_l3 already exists
    cms_l3_csv_file = rf'{drive[comp]}\Phd_data\copernicus\{tc_storm_id}_{case_study[1]}_l3_swh.nc.csv'
    if os.path.exists(cms_l3_csv_file):
        cms_l3_df = pd.read_csv(cms_l3_csv_file)
        
    else:
        # Load CMS L3 files (.nc.csv preferred, .nc fallback) into pandas dataframe
        cms_l3_files = glob.glob(rf'{drive[comp]}\Phd_data\copernicus\{tc_storm_id}_{case_study[1]}_*_l3_swh.nc.csv')
        if len(cms_l3_files) == 0:
            cms_l3_files = glob.glob(rf'{drive[comp]}\Phd_data\copernicus\{tc_storm_id}_{case_study[1]}_*_l3_swh.nc')
        if len(cms_l3_files) > 0:
            cms_l3_data_list = []
            for l3_file in cms_l3_files:
                if l3_file.lower().endswith('.csv'):
                    file_df = pd.read_csv(l3_file)
                else:
                    with nc.Dataset(l3_file, 'r') as cms_nc:
                        # Extract all variables into a dictionary
                        file_data = {}
                        for var_name in cms_nc.variables:
                            var_obj = cms_nc.variables[var_name]
                            var_data = var_obj[:]
                            if var_name.lower() == 'time':
                                time_units = getattr(var_obj, 'units', None)
                                if time_units is not None:
                                    time_calendar = getattr(var_obj, 'calendar', 'standard')
                                    var_data = np.array(nc.num2date(var_data, time_units, time_calendar), dtype='datetime64[s]')
                            # Flatten multi-dimensional arrays if needed
                            if var_data.ndim > 1:
                                # For 2D+ arrays, flatten to 1D
                                file_data[var_name] = var_data.flatten()
                            else:
                                file_data[var_name] = var_data
                        
                        # Create dataframe from this file
                        # Handle different array lengths by repeating shorter dimensions
                        max_len = max(len(v) for v in file_data.values())
                        for key in file_data:
                            if len(file_data[key]) == 1:
                                file_data[key] = np.repeat(file_data[key], max_len)
                            elif len(file_data[key]) != max_len:
                                # If dimensions don't match, create a mesh grid
                                pass  # Handle case-by-case if needed
                        file_df = pd.DataFrame(file_data)

                rename_map = {}
                if 'longitude' in file_df.columns and 'lon' not in file_df.columns:
                    rename_map['longitude'] = 'lon'
                if 'latitude' in file_df.columns and 'lat' not in file_df.columns:
                    rename_map['latitude'] = 'lat'
                if rename_map:
                    file_df = file_df.rename(columns=rename_map)

                if 'time' in file_df.columns:
                    file_df['time'] = pd.to_datetime(file_df['time'], errors='coerce').values.astype('datetime64[s]')

                cms_l3_data_list.append(file_df)
            
            # Concatenate all dataframes
            cms_l3_df = pd.concat(cms_l3_data_list, ignore_index=True)
            # Vectorized nearest-track matching within 30 minutes
            cms_l3_df['dist_to_track'] = np.nan
            cms_l3_df['time_diff'] = np.timedelta64(0, 's')
            required_cols = {'lat', 'lon', 'time'}
            if required_cols.issubset(cms_l3_df.columns) and len(cms_l3_df) > 0 and len(tc_time_interpolated) > 0:
                cms_l3_df['time'] = pd.to_datetime(cms_l3_df['time'], errors='coerce').values.astype('datetime64[s]')

                cms_coords = cms_l3_df[['lat', 'lon']].to_numpy(dtype=float)
                track_coords = tc_coords
                cms_times = cms_l3_df['time'].to_numpy(dtype='datetime64[s]')
                track_times = tc_time_interpolated.astype('datetime64[s]')

                distances = np.abs(haversine_vectorized(cms_coords, track_coords))
                time_diffs = np.abs(cms_times[:, None] - track_times[None, :])

                time_window_mask = time_diffs <= np.timedelta64(30, 'm')
                valid_distances = np.where(time_window_mask, distances, np.inf)

                nearest_track_idx = np.argmin(valid_distances, axis=1)
                row_idx = np.arange(valid_distances.shape[0])
                min_dists = valid_distances[row_idx, nearest_track_idx]
                has_match = np.isfinite(min_dists)

                dist_out = np.full(len(cms_l3_df), np.nan, dtype=float)
                dist_out[has_match] = min_dists[has_match]

                time_out = np.full(len(cms_l3_df), np.timedelta64(0, 's'), dtype='timedelta64[s]')
                selected_time_diffs = time_diffs[row_idx, nearest_track_idx].astype('timedelta64[s]')
                time_out[has_match] = selected_time_diffs[has_match]

                cms_l3_df['dist_to_track'] = dist_out
                cms_l3_df['time_diff'] = time_out
            # assign True for rows with a match within 30 minutes and 100 km, False otherwise
            cms_l3_df['match_to_track'] = (cms_l3_df['dist_to_track'] <= MAX_DISTANCE) & (cms_l3_df['time_diff'] <= np.timedelta64(4, 'h'))
            # save to csv for future use
            cms_l3_df.to_csv(cms_l3_csv_file, index=False)
        else:
            cms_l3_df = pd.DataFrame()

    if len(cms_l3_df) > 0 and 'time' in cms_l3_df.columns:
        cms_l3_df = cms_l3_df.copy()
        cms_l3_df['time'] = pd.to_datetime(cms_l3_df['time'], errors='coerce').to_numpy(dtype='datetime64[s]')
        cms_l3_df = cms_l3_df.loc[~np.isnat(cms_l3_df['time'].to_numpy(dtype='datetime64[s]'))].reset_index(drop=True)
    
    # ERA5 and CMS
    from CMS_CDS_api import retrieve_vals
    filepath_swh = f"E:\\Phd_data\\copernicus\\{tc_storm_id}_SWH.nc"
    filepath_sss_sst = f"E:\\Phd_data\\copernicus\\{tc_storm_id}_SSS&SST.nc"
    filepath_era5 = f"E:\\Phd_data\\copernicus\\{tc_storm_id}_era5.grib"
    cms_res_SWH = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id,filepath_swh,['VHM0','VMDR'])
    cms_res_SSS = retrieve_vals(cyg_lats, cyg_lons, cyg_times,tc_storm_id,filepath_sss_sst,['thetao','so'])
    era5_res = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id, filepath_era5, ['msl', 'u10','v10','t2m','sp'])
    import cfgrib
    try:
        era5_dataset = cfgrib.open_dataset(filepath_era5, errors='ignore')
    except TypeError:
        # Fallback for cfgrib versions that do not expose the errors argument
        era5_dataset = cfgrib.open_dataset(filepath_era5)
    era5_time = era5_dataset.time.data ; era5_lat = era5_dataset.latitude.data ; era5_lon = era5_dataset.longitude.data
    era5_wind = np.sqrt(era5_dataset.u10.data**2 + era5_dataset.v10.data**2)
    era5_mslp = era5_dataset.msl.data
    era5_t2m = era5_dataset.t2m.data
    era5_dataset.close()
    swh_dataset = nc.Dataset(filepath_swh)
    swh_time_var = swh_dataset.variables['time']
    swh_lat = swh_dataset.variables['latitude'][:]
    swh_lon = swh_dataset.variables['longitude'][:]
    swh_time = nc.num2date(swh_time_var[:], swh_time_var.units, getattr(swh_time_var, 'calendar', 'standard'))
    swh_time = np.array(swh_time, dtype='datetime64[h]')
    swh_data = swh_dataset.variables['VHM0'][:]
    vmdr = swh_dataset.variables['VMDR'][:]
    swh_dataset.close()
    SSS_SST_dataset = nc.Dataset(filepath_sss_sst)
    sss_sst_time_var = SSS_SST_dataset.variables['time']
    sss_sst_lat = SSS_SST_dataset.variables['latitude'][:]
    sss_sst_lon = SSS_SST_dataset.variables['longitude'][:]
    sss_sst_time = nc.num2date(sss_sst_time_var[:], sss_sst_time_var.units, getattr(sss_sst_time_var, 'calendar', 'standard'))
    sss_sst_time = np.array(sss_sst_time, dtype='datetime64[h]')
    sss_sst_thetao = SSS_SST_dataset.variables['thetao'][:]
    sss_sst_so = SSS_SST_dataset.variables['so'][:]
    SSS_SST_dataset.close()
    # make sure all lons are -180 to 180
    swh_lon = np.where(swh_lon > 180, swh_lon - 360, swh_lon)
    sss_sst_lon = np.where(sss_sst_lon > 180, sss_sst_lon - 360, sss_sst_lon)
    era5_lon = np.where(era5_lon > 180, era5_lon - 360, era5_lon)

    # MSWEP precipitation
    mswep_casestudy_file = rf'{drive[comp]}\Phd_data\mswep\{tc_storm_id}_combined.csv'
    mswep_precip_list = np.array([])
    mswep_lats_list = np.array([])
    mswep_lons_list = np.array([])
    mswep_times_list = np.array([], dtype='datetime64[s]')
    if not os.path.exists(mswep_casestudy_file):
        mswep_folder = rf'{drive[comp]}\Phd_data\mswep'
        mswep_files = glob.glob(mswep_folder + r"\*.nc")
        mswep_precip_list = []
        mswep_lats_list = []
        mswep_lons_list = []
        mswep_times_list = []
        # find the day of year for the case study
        case_day_start = (np.datetime64(cyg_AOI_time_start, 'D') - np.datetime64(cyg_AOI_time_start,'Y')).astype(int)
        case_day_end = (np.datetime64(cyg_AOI_time_end, 'D') - np.datetime64(cyg_AOI_time_end,'Y')).astype(int) + 1
        case_year = cyg_AOI_time_start.year
        # Check if the string "case_year + day of year" exists in any of the MSWEP files
        mswep_files = np.array(mswep_files)
        # assert np.any(case_year + str(case_day_start) == mswep_files[:,21:28]), "MSWEP First day files for the case study period not found."
        # assert case_year + str(case_day_end) in mswep_files, "MSWEP Last day files for the case study period not found."

        for i in range(case_day_start, case_day_end):
            for hour in ['00','03','06','09','12','15','18','21']:
                mswep_nc = nc.Dataset(rf'{mswep_folder}\{case_year}{i}.{hour}.nc')
                # Extract precipitation variable and coordinates in the AOI and append to numpy arrays
                mswep_lats = mswep_nc.variables['lat'][:]
                mswep_lons = mswep_nc.variables['lon'][:]
                mswep_precip = mswep_nc.variables['precipitation'][:]
                # Find indices within AOI
                lat_mask = (mswep_lats >= AOI_lat_min) & (mswep_lats <= AOI_lat_max)
                lon_mask = (mswep_lons >= AOI_lon_min) & (mswep_lons <= AOI_lon_max)
                # Extract values within AOI as 2D (lat, lon)
                precip_aoi = mswep_precip[0][lat_mask][:, lon_mask]
                lats_aoi = mswep_lats[lat_mask]
                lons_aoi = mswep_lons[lon_mask]
                # Append to arrays
                if len(mswep_lats_list) == 0:
                    mswep_lats_list = lats_aoi
                    mswep_lons_list = lons_aoi
                mswep_precip_list.append(precip_aoi)
                mswep_times_list.append(
                    np.datetime64(f"{case_year:04d}-01-01T00:00:00")
                    + np.timedelta64(i - 1, 'D')
                    + np.timedelta64(int(hour), 'h')
                )

        mswep_times_list = np.asarray(mswep_times_list, dtype='datetime64[s]')
        mswep_precip_list = np.asarray(mswep_precip_list)

        # # Save combined MSWEP data to CSV
        # with open(mswep_casestudy_file, 'w', newline='') as f:
        #     writer = csv.writer(f)
        #     writer.writerow(['time', 'lat', 'lon', 'precipitation'])
        #     for i in range(len(mswep_times_list)):
        #         writer.writerow([mswep_times_list[i], mswep_lats_list[i], mswep_lons_list[i], mswep_precip_list[i]])    


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

    #%%
    # interactive app
    interactive_app = False
    if interactive_app:
        from RI_casestudy_app import (
            build_measurements_dataframe,
            build_track_dataframe,
            launch_case_study_app,
        )

        ri_start_time = tc_time[index_before]
        ri_end_time = tc_time[len(tc_time)-index_after-1]

        app_meas_df = build_measurements_dataframe(
            cyg_times=meas_times,
            cyg_lats=meas_lats,
            cyg_lons=meas_lons,
            nearest_tc_coords=nearest_tc_coords[valid_flag][sorted_indices],
            quadrants=quadrants,
            meas_data=meas_data,
            ri_start_time=ri_start_time,
            ri_end_time=ri_end_time,
            cms_res_swh=cms_res_SWH,
            era5_res=era5_res,
        )
        app_track_df = build_track_dataframe(
            tc_time_interpolated,
            tc_lat_interpolated,
            tc_lon_interpolated,
            tc_vmax_interpolated,
        )

        environmental_grids = {
            'SST': {
                'time': sss_sst_time,
                'lat': sss_sst_lat,
                'lon': sss_sst_lon,
                'data': sss_sst_thetao,
                'window': np.timedelta64(50, 'm'),
            },
            'SSS': {
                'time': sss_sst_time,
                'lat': sss_sst_lat,
                'lon': sss_sst_lon,
                'data': sss_sst_so,
                'window': np.timedelta64(50, 'm'),
            },
            'SWH': {
                'time': swh_time,
                'lat': swh_lat,
                'lon': swh_lon,
                'data': swh_data,
                'window': np.timedelta64(177, 'm'),
            },
            'ERA5 Wind': {
                'time': era5_time,
                'lat': era5_lat,
                'lon': era5_lon,
                'data': era5_wind,
                'window': np.timedelta64(50, 'm'),
            },
            'ERA5 MSLP': {
                'time': era5_time,
                'lat': era5_lat,
                'lon': era5_lon,
                'data': era5_mslp,
                'window': np.timedelta64(50, 'm'),
            },
            'ERA5 T2M': {
                'time': era5_time,
                'lat': era5_lat,
                'lon': era5_lon,
                'data': era5_t2m,
                'window': np.timedelta64(50, 'm'),
            },
        }

        if mswep_times_list.size > 0 and np.size(mswep_precip_list) > 0:
            environmental_grids['MSWEP Precip'] = {
                'time': mswep_times_list,
                'lat': mswep_lats_list,
                'lon': mswep_lons_list,
                'data': mswep_precip_list,
                'window': np.timedelta64(3, 'h'),
            }

        launch_case_study_app(
            tc_df=app_track_df,
            meas_df=app_meas_df,
            brcs=meas_data.get('brcs', None),
            eff_scatter=meas_data.get('eff_scatter', None),
            cms_buoy_data=cms_buoy_data,
            environmental_grids=environmental_grids,
            port=8050,
        )

    quadrant_plotting = False
    if quadrant_plotting:
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
        ri_end_time = tc_time[len(tc_time)-index_after]
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
            
        plt.tight_layout(rect=(0.03, 0.10, 0.97, 0.93))    # Adjust spacing to leave space at the bottom
        cbar = fig.colorbar(img_handles[0], cax=fig.add_axes((0.35, 0.08, 0.3, 0.02)), orientation='horizontal') # Manually add a smaller horizontal colorbar
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
    
    #%%
    animation = True
    if animation:
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
        latitudes = tc_lat_interpolated ; longitudes = tc_lon_interpolated; times = np.asarray(tc_time_interpolated, dtype='datetime64[s]')
        cyg_times_all = np.asarray(cyg_times_all, dtype='datetime64[s]')
        cms_l3_times = np.array([], dtype='datetime64[s]')
        if len(cms_l3_df) > 0 and 'time' in cms_l3_df.columns:
            cms_l3_times = cms_l3_df['time'].to_numpy(dtype='datetime64[s]')
        winds = tc_vmax_interpolated
        norm = mcolors.Normalize(vmin=15, vmax=45)
        cmap = cm.get_cmap('coolwarm')
        gif_output_folder = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\gifs'.format(comp)
        os.makedirs(gif_output_folder, exist_ok=True)

        background_layers = [
            {'key': 'sst', 'label': 'SST', 'unit': '°C', 'source': 'sss_sst', 'data': sss_sst_thetao, 'cmap': 'Reds', 'window': np.timedelta64(50, 'm'), 'vmin': 24, 'vmax': 32},
            {'key': 'sss', 'label': 'SSS', 'unit': 'psu', 'source': 'sss_sst', 'data': sss_sst_so, 'cmap': 'viridis', 'window': np.timedelta64(50, 'm'), 'vmin': 0, 'vmax': 40},
            {'key': 'swh', 'label': 'SWH', 'unit': 'm', 'source': 'swh', 'data': swh_data, 'cmap': 'Blues', 'window': np.timedelta64(177, 'm'), 'vmin': 0, 'vmax': 9},
            {'key': 'era5_wind', 'label': 'ERA5 Wind', 'unit': 'm/s', 'source': 'era5', 'data': era5_wind, 'cmap': 'coolwarm', 'window': np.timedelta64(50, 'm'), 'vmin': 15, 'vmax': 45},
            {'key': 'era5_mslp', 'label': 'ERA5 MSLP', 'unit': 'Pa', 'source': 'era5', 'data': era5_mslp, 'cmap': 'PuOr', 'window': np.timedelta64(50, 'm'), 'vmin': 950000, 'vmax': 1000000},
            {'key': 'era5_t2m', 'label': 'ERA5 T2M', 'unit': 'K', 'source': 'era5', 'data': era5_t2m, 'cmap': 'coolwarm', 'window': np.timedelta64(50, 'm'), 'vmin': np.nanmin(era5_t2m), 'vmax': np.nanmax(era5_t2m)},
        ]

        if mswep_times_list.size > 0 and np.size(mswep_precip_list) > 0:
            background_layers.append(
                {'key': 'mswep', 'label': 'MSWEP Precip', 'unit': 'mm/h', 'source': 'mswep', 'data': mswep_precip_list, 'cmap': 'Blues', 'window': np.timedelta64(3, 'h'), 'vmin': 0, 'vmax': 200}
            )

        # Quick switch: set to None for all layers, or a list like ['sst', 'swh', 'era5_mslp']
        enabled_background_layers = None

        if enabled_background_layers is None:
            layers_to_plot = background_layers
        else:
            enabled_set = set(enabled_background_layers)
            layers_to_plot = [layer for layer in background_layers if layer['key'] in enabled_set]
            missing_keys = sorted(enabled_set - {layer['key'] for layer in background_layers})
            if missing_keys:
                print(f"Warning: unknown background layer keys ignored: {missing_keys}")
            if len(layers_to_plot) == 0:
                raise ValueError('No valid background layers selected in enabled_background_layers')

        def _plot_background(ax, layer, frame):
            target_time = times[frame]
            mask = None
            if layer['source'] == 'swh':
                mask = (swh_time <= target_time) & (swh_time >= target_time - layer['window'])
                if np.any(mask):
                    return ax.pcolormesh(swh_lon, swh_lat, np.squeeze(layer['data'][mask]), cmap=layer['cmap'], vmin=layer['vmin'], vmax=layer['vmax'], shading='auto', alpha=0.6)
            elif layer['source'] == 'sss_sst':
                mask = (sss_sst_time <= target_time) & (sss_sst_time >= target_time - layer['window'])
                if np.any(mask):
                    return ax.pcolormesh(sss_sst_lon, sss_sst_lat, np.squeeze(layer['data'][mask]), cmap=layer['cmap'], vmin=layer['vmin'], vmax=layer['vmax'], shading='auto', alpha=0.6)
            elif layer['source'] == 'era5':
                mask = (era5_time <= target_time) & (era5_time >= target_time - layer['window'])
                if np.any(mask):
                    return ax.pcolormesh(era5_lon, era5_lat, np.squeeze(layer['data'][mask]), cmap=layer['cmap'], vmin=layer['vmin'], vmax=layer['vmax'], shading='auto', alpha=0.6)
            elif layer['source'] == 'mswep':
                mask = (mswep_times_list <= target_time) & (mswep_times_list >= target_time - layer['window'])
                if np.any(mask):
                    mswep_precip_plot = np.nanmean(layer['data'][mask], axis=0) 
                    return ax.pcolormesh(mswep_lons_list, mswep_lats_list, mswep_precip_plot, cmap=layer['cmap'], vmin=layer['vmin'], vmax=layer['vmax'], shading='auto', alpha=0.6)
            return None

        pause_frames = 10
        frames_to_plot = list(range(len(times))) + [len(times) - 1] * pause_frames

        lat_span = max(float(np.nanmax(latitudes) - np.nanmin(latitudes)), 0.5)
        lon_span = max(float(np.nanmax(longitudes) - np.nanmin(longitudes)), 0.5)
        fig_height = 6.0
        fig_width = float(np.clip(fig_height * (lon_span / lat_span), 10.0, 22.0))

        for layer in layers_to_plot:
            fig, ax = plt.subplots(figsize=(fig_width, fig_height))
            wind_sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            wind_sm.set_array([])
            plt.colorbar(wind_sm, ax=ax, orientation='vertical', label='Wind Speed (m/s)')

            bg_norm = mcolors.Normalize(vmin=layer['vmin'], vmax=layer['vmax'])
            bg_sm = plt.cm.ScalarMappable(cmap=layer['cmap'], norm=bg_norm)
            bg_sm.set_array([])
            plt.colorbar(bg_sm, ax=ax, orientation='vertical', label=f"{layer['label']} ({layer['unit']})")

            def update_frame(frame):
                ax.clear()
                ax.set_xlim(min(longitudes) - buffer, max(longitudes) + buffer)
                ax.set_ylim(min(latitudes) - buffer, max(latitudes) + buffer)

                mask = (cyg_times_all <= times[frame]) & (cyg_times_all >= times[frame] - np.timedelta64(4, 'h'))
                ax.scatter(cyg_lons_all[mask], cyg_lats_all[mask], c='grey', marker='x', s=10)

                # L3 SWH
                if len(cms_l3_df) > 0:
                    norm_cms = mcolors.Normalize(vmin=0, vmax=20) # vmin=min(dist_cyg_stat), vmax=max(dist_cyg_stat)

                    mask = (cms_l3_times <= times[frame]) & (cms_l3_times >= times[frame] - np.timedelta64(4, 'h'))
                    ax.scatter(cms_l3_df.loc[mask, 'lon'], cms_l3_df.loc[mask, 'lat'], norm=norm_cms, marker='o', s=10,cmap=cm.Reds_r, label='CMS L3 SWH')                    

                ax.quiver(
                    longitudes[:frame + 1],
                    latitudes[:frame + 1],
                    dlon[:frame + 1],
                    dlat[:frame + 1],
                    winds[:frame + 1],
                    cmap='coolwarm',
                    width=0.005,
                    norm=norm,
                    angles='xy',
                    scale_units='xy',
                )

                ring = Circle((longitudes[frame], latitudes[frame]), 0.95, fill=False, color='black', linestyle='--', linewidth=1)
                ax.add_patch(ring)
                ring = Circle((longitudes[frame], latitudes[frame]), 0.05, fill=False, color='black', linestyle='-', linewidth=0.5)
                ax.add_patch(ring)

                mask = (meas_times <= times[frame]) & (meas_times >= times[frame] - np.timedelta64(4, 'h'))
                color_list = [quadrant_colors[q] for q in quadrants[mask]]
                ax.scatter(meas_lons[mask], meas_lats[mask], c=color_list, marker='o', s=20)

                mesh = _plot_background(ax, layer, frame)
                if mesh is None:
                    ax.text(0.02, 0.95, 'No background data in window', transform=ax.transAxes, fontsize=10, ha='left', color='black', bbox=dict(facecolor='white', alpha=0.6))

                ax.text(
                    0.98,
                    0.95,
                    f"Time: {times[frame]}",
                    transform=ax.transAxes,
                    fontsize=12,
                    ha='right',
                    color='black',
                    bbox=dict(facecolor='white', alpha=0.6),
                )

                ax.set_title('ID: {0}, Duration: {1}, Intensity: {2}, Background: {3}'.format(tc_storm_id, durations, intensities, layer['label']))
                legend_handles = [mpatches.Patch(color=color, label=label) for label, color in quadrant_colors.items()]
                ax.legend(handles=legend_handles, loc='lower right', title='Quadrants')
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')

            ani = animation.FuncAnimation(fig, update_frame, frames=frames_to_plot, repeat=True)

            output_path = os.path.join(gif_output_folder, 'L1_{0}_{1}_{2}.gif'.format(tc_storm_id, case_study[1], layer['key']))
            ani.save(output_path, fps=2)
            plt.close(fig)
            print(f"Saved GIF: {output_path}")

#%%

# IMERG
# imerg_files_list = np.asarray(glob.glob(case_folder+ r"\IMERG\*.nc4"))
# cyg_nc = nc.Dataset(imerg_files_list[0])
# imerg_lons = cyg_nc.variables['lon'][:] ; imerg_lats = cyg_nc.variables['lat'][:] ; 
# imerg_start_time_str = cyg_nc.variables['time'].units[14:33]
# cyg_nc.close
# imerg_times = np.ma.zeros((len(imerg_files_list))); imerg_precip = np.ma.zeros((len(imerg_files_list),len(imerg_lons),len(imerg_lats))); 
# imerg_precip_quality = np.ma.zeros((len(imerg_files_list),len(imerg_lons),len(imerg_lats)));
# for i,cyg_file_path in enumerate(imerg_files_list):
#     cyg_nc = nc.Dataset(cyg_file_path)
#     imerg_precip[i]= cyg_nc.variables['precipitation'][:] ; imerg_precip_quality[i]= cyg_nc.variables['precipitationQualityIndex'][:]
#     imerg_times[i] =cyg_nc.variables['time'][0] 
#     cyg_nc.close
# imerg_times = np.datetime64(imerg_start_time_str)+ imerg_times.astype('timedelta64[s]')