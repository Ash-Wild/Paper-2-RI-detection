# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:46:38 2024

@author: ashle
"""
#%%
# notes: 

import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta 
import pandas as pd
import glob
import os
import pickle
import re
from typing import Optional, Tuple
from Storm_timeline_visualiser import plot_timeline_both, plot_distribution, plot_tc_animation # for visualising timeline
from CMS_CDS_api import cms_downloader, era5_downloader, retrieve_vals, expand_and_redownload
from RI_casestudy_insitu_finder import haversine_vectorized

# from scipy.spatial import cKDTree
comp = 'ashle'
case_study_ind = None # None 210 428 665 4

buffer = 1 # degrees buffer to search in storm_centric
interp_interval= 20 # minutes to interpolate best-track
noaa_time_threshold = 3 # hours - maximum time difference allowed between a CYGNSS measurement and a TC track point to be considered a match (used in both BallTree and fallback)
MAX_DISTANCE = 80  # Maximum distance in kilometers
load_checkpoint = False
save_checkpoint = False
skip_event_processing = False

# L2 sample-flag option: True: keep all L2 measurements, including flagged values.- False: keep only measurements where yslf_sample_flags == 0.
yslf_max_wind_uncertainty = None # 8.0 m/s
yslf_incidence_angle_range = None # (10, 70)
yslf_min_snr = None # 1.3
yslf_range_corr_gain_range = (35, 9000)
snr_rx_filter = False  # Enable combined SNR/Rx gain filter to remove low-quality signal/gain combinations
yslf_min_num_ddms_utilized = None # 2
variables = ['range_corr_gain','nbrcs_mean','wind_speed','cyg_vmax_error','time_since_ri'] 
# variables = ['incidence_angle', 'rx_gain', 'snr', 'range_corr_gain','num_ddms_utilized','nbrcs_mean_corrected','nbrcs_mean','wind_speed_uncertainty','wind_speed','cyg_vmax_error','cyg_dist','time_since_ri'] 
# #'sample_flags', 'ddm_sample_index', 'ddm_channel','SST','SSS','SWH','wave_dir',"mean_sea_level_pressure","era5_wind_speed","2m_temperature"]
cyg_match_var = 'cyg_dist'  # choose variable used where cyg_dist was used downstream (e.g., 'incidence_angle')
checkpoint_file = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Data', 'IBTrACS', 'RI_analysis_optimised_checkpoint_all_sources_RCG35.pkl')
noaa_cyg_folder = r'E:\Phd_data\cyg_noaa_1.2' # r'C:\Users\\'+comp+'\\Documents\CYGNSS_NOAA'  # 
l2_3p2_cyg_folder = r'E:\Phd_data\cyg_l2_3.2'
merged_cyg_folder = r'E:\Phd_data\cyg_storm_centric' # r'E:\Phd_data\cyg_storm_centric' r'C:\Users\\'+comp+'\\Documents\cyg_storm_centric' 
noaa_like_folders = {'l2_3.2': l2_3p2_cyg_folder,'noaa_1.2': noaa_cyg_folder} # 
noaa_like_files_map = {source_name: np.asarray(glob.glob(source_folder + r"\*.nc")) for source_name, source_folder in noaa_like_folders.items()}
# Fallback mappings for known schema differences (especially CYGNSS L2 v3.2 files).
var_aliases = {
    'rx_gain': ['port_gain_setting', 'starboard_gain_setting'],
    'snr': ['bit_ratio_lo_hi_port', 'bit_ratio_lo_hi_starboard', 'bit_ratio_lo_hi_zenith'],
    'nbrcs_mean_corrected': ['ddm_nbrcs'],
}
# Source-specific aliases let us choose algorithm-specific products.
# For L2 v3.2 we prefer YSLF fields instead of the default FDS/MV fields.
source_var_aliases = {
    'l2_3.2': {
        'wind_speed': ['preliminary_yslf_wind_speed'],
        'wind_speed_uncertainty': ['preliminary_yslf_wind_speed_uncertainty'],
        'sample_flags': ['yslf_sample_flags'], #quality_flags
    }
}
# For some derived fields, combine multiple instrument channels into one value.
# We compute an element-wise mean across available channels to keep a single
# analysis variable that is comparable with older NOAA-style naming.
avg_alias_fields = {
    'rx_gain': ['port_gain_setting', 'starboard_gain_setting'],
    'snr': ['bit_ratio_lo_hi_port', 'bit_ratio_lo_hi_starboard', 'bit_ratio_lo_hi_zenith'],
}
merged_cyg_files_list = np.asarray(glob.glob(merged_cyg_folder+ r"\*.nc"))

def parse_cyg_file_date(cyg_file_path):
    """Parse daily file date from NOAA 1.2 and CYGNSS L2 v3.2 filename patterns."""
    fname = os.path.basename(cyg_file_path)
    match = re.search(r"\.s(\d{8})-", fname)
    if not match:
        start_ind = fname.find('ddmi.')
        if start_ind >= 0:
            candidate = fname[start_ind+6:start_ind+14]
            if len(candidate) == 8 and candidate.isdigit():
                match = [None, candidate]
    if not match:
        return None
    date_str = match.group(1) if hasattr(match, 'group') else match[1]
    return np.datetime64(datetime.strptime(date_str, "%Y%m%d").date(), 'D')

def _extract_box_float_var(cached_vars, varname, box_indices):
    """Return selected variable values as float array aligned to box_indices."""
    var_full = cached_vars.get(varname, np.array([]))
    if getattr(var_full, 'size', 0) == 0:
        return np.array([])
    try:
        selected = np.asarray(np.ma.filled(var_full[box_indices], np.nan), dtype=float)
    except Exception:
        return np.array([])
    return selected

def _as_text(value):
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='ignore')
    if isinstance(value, np.bytes_):
        return value.astype(str)
    return str(value)

def _track_key(track_value):
    """Convert track IDs into a comparable hashable key."""
    if isinstance(track_value, np.ndarray):
        arr = np.asarray(track_value)
        if arr.ndim == 0:
            return (arr.item(),)
        return tuple(arr.tolist())
    if isinstance(track_value, list):
        return tuple(track_value)
    if isinstance(track_value, tuple):
        return track_value
    return (track_value,)

RI_file = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v5.nc'
RI_nc = nc.Dataset(RI_file)
time_int = RI_nc.variables['times'][:]
time = np.datetime64('1970-01-01T00:00:00') + time_int.astype('timedelta64[s]')
storm_id = RI_nc.variables['storm'][:]
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
noaa_file_dates = {source_name: [] for source_name in noaa_like_files_map}
cyg_storm_names = []
cyg_files_num = []
file_counter = 0
for cyg_file in merged_cyg_files_list:
    start_ind = cyg_file.find('ddmi.')
    end_ind = cyg_file.find('.20')
    storm_name = cyg_file[start_ind+5:end_ind-3]
    with nc.Dataset(cyg_file) as cyg_nc:
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
date_to_noaa_idx = {source_name: {} for source_name in noaa_like_files_map}
for source_name, cyg_files in noaa_like_files_map.items():
    for idx, cyg_file in enumerate(cyg_files):
        noaa_date = parse_cyg_file_date(cyg_file)
        if noaa_date is None:
            continue
        cyg_datetime_object = pd.Timestamp(noaa_date).to_pydatetime()
        date_to_noaa_idx[source_name][noaa_date] = idx
        noaa_file_dates[source_name].append(cyg_datetime_object.date())
noaa_cache = {source_name: {} for source_name in noaa_like_files_map}  # cache per source to avoid reopening files repeatedly

if cyg_match_var not in variables and cyg_match_var not in ('cyg_dist', 'cyg_vmax_error', 'wind_speed', 'time_since_ri'):
    variables.append(cyg_match_var)

results_dt = pd.DataFrame(columns=['under/over','intensityBT','durationBT','storm_idBT','tc_eventBT','wind_change','time_diff','num_tracks'])
noaa_df_list = []; tc_list= []; merged_list =[]; l2_df_list=[]
noaa_ri_change = [] ;merged_ri_change = [] ; noaa_wind_change=[] ; merged_wind_change=[] ;noaa_ri_durations=[] ;merged_ri_durations=[] ; noaa_ri_initial = [] ; merged_ri_initial = []
l2_3p2_ri_change = [] ; l2_3p2_wind_change = [] ; l2_3p2_ri_durations = [] ; l2_3p2_ri_initial = []


if load_checkpoint and os.path.exists(checkpoint_file):
    with open(checkpoint_file, 'rb') as f:
        checkpoint = pickle.load(f)
    results_dt = checkpoint.get('results_dt', results_dt)
    noaa_df_list = checkpoint.get('noaa_df_list', noaa_df_list)
    l2_df_list = checkpoint.get('l2_df_list', l2_df_list)
    tc_list = checkpoint.get('tc_list', tc_list)
    merged_list = checkpoint.get('merged_list', merged_list)
    noaa_ri_change = checkpoint.get('noaa_ri_change', noaa_ri_change)
    merged_ri_change = checkpoint.get('merged_ri_change', merged_ri_change)
    noaa_wind_change = checkpoint.get('noaa_wind_change', noaa_wind_change)
    merged_wind_change = checkpoint.get('merged_wind_change', merged_wind_change)
    noaa_ri_durations = checkpoint.get('noaa_ri_durations', noaa_ri_durations)
    merged_ri_durations = checkpoint.get('merged_ri_durations', merged_ri_durations)
    noaa_ri_initial = checkpoint.get('noaa_ri_initial', noaa_ri_initial)
    merged_ri_initial = checkpoint.get('merged_ri_initial', merged_ri_initial)
    l2_3p2_ri_change = checkpoint.get('l2_3p2_ri_change', l2_3p2_ri_change)
    l2_3p2_wind_change = checkpoint.get('l2_3p2_wind_change', l2_3p2_wind_change)
    l2_3p2_ri_durations = checkpoint.get('l2_3p2_ri_durations', l2_3p2_ri_durations)
    l2_3p2_ri_initial = checkpoint.get('l2_3p2_ri_initial', l2_3p2_ri_initial)
    skip_event_processing = True
    print(f"Loaded checkpoint: {checkpoint_file}")

for tc_event in range(0 if skip_event_processing else len(storm_id)): #Hin 428 Bua 200 Kya 210 Hele 665
    if ri_rd[tc_event]==1 and durations[tc_event] < 25: # and tc_event == 428:# if time[tc_event][0] >np.datetime64('2022-12-31T23:59:59'): '2024279N21265'
        if case_study_ind is not None and tc_event != case_study_ind:
            continue
        
        # need to know the length of storm and remove nans
        nans = ~np.isnan(latitude[tc_event]).data
        tc_time = time[tc_event,nans]
        tc_storm_id = storm_id[tc_event]
        
        tc_vmax = vmax[tc_event,nans]
        tc_durations = durations[tc_event]
        tc_intensities = intensities[tc_event]
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
        # create arrays to store interpolated values
        # Vectorized interpolation using seconds-since-epoch (faster than looping)
        tc_times_s = tc_time.astype('datetime64[s]').astype('int64')
        interp_s = tc_time_interpolated.astype('datetime64[s]').astype('int64')

        # Ensure numeric arrays (fill masked vmax with nan)
        tc_vmax_vals = np.ma.filled(tc_vmax, np.nan).astype(float)
        tc_lat_vals = np.asarray(tc_latitude).astype(float)
        tc_lon_vals = np.asarray(tc_longitude).astype(float)

        # Interpolate (np.interp extrapolates using endpoint values similar to previous logic)
        tc_vmax_interpolated = np.interp(interp_s, tc_times_s, tc_vmax_vals)
        tc_lat_interpolated = np.interp(interp_s, tc_times_s, tc_lat_vals)
        tc_lon_interpolated = np.interp(interp_s, tc_times_s, tc_lon_vals)

        # define AOI box
        AOI_lat_min,AOI_lat_max = np.nanmin(tc_latitude)-buffer,np.nanmax(tc_latitude)+buffer
        AOI_lon_min,AOI_lon_max = np.nanmin(tc_longitude)-buffer,np.nanmax(tc_longitude)+buffer
        cyg_AOI_time_start,cyg_AOI_time_end = np.nanmin(tc_time),np.nanmax(tc_time)
        cyg_AOI_time_start = cyg_AOI_time_start.astype('datetime64[s]').astype('O')
        cyg_AOI_time_end = cyg_AOI_time_end.astype('datetime64[s]').astype('O')
        date_list = np.arange(cyg_AOI_time_start, cyg_AOI_time_end+np.timedelta64(1, 'D'), np.timedelta64(1, 'D'), dtype='datetime64[D]')
        
        noaa_res = (False, 0, 0)
        l2_res = (False, 0, 0)
        noaa_track_count = 0
        l2_track_count = 0
        
        # Use Python lists for accumulating coincident measurements per variable
        coincident_meas = {i: [] for i in variables}
        coincident_meas['cyg_source'] = []
        # Use Python lists for incremental appends (faster than repeated np.append)
        cyg_lats = []
        cyg_lons = []
        cyg_winds = []
        cyg_times = []
        cyg_tracks = []
        cyg_dist = []
        cyg_match = []
        cyg_vmax_error = []
        cyg_sources = np.array([], dtype=object)
        store_all = False
        if case_study_ind is not None:
            store_all = True
            cyg_lats_all = []
            cyg_lons_all = []
            cyg_winds_all = []
            cyg_times_all = []
            cyg_tracks_all = []
            cyg_sources_all = []
        for source_name, source_files in noaa_like_files_map.items():
            source_date_to_idx = date_to_noaa_idx.get(source_name, {})
            source_cache = noaa_cache.get(source_name, {})

            for date in date_list:
                # Fast lookup to file path
                try:
                    cyg_file_path = str(source_files[source_date_to_idx[date]])
                except KeyError:
                    continue  # No file for this date, skip to next date

                # Load from cache if available; else open file once and cache necessary arrays
                if cyg_file_path in source_cache:
                    cached = source_cache[cyg_file_path]
                    cygnss_lons = cached['lons']
                    cygnss_lats = cached['lats']
                    cygnss_wind = cached['wind']
                    cygnss_tracks = cached['tracks']
                    cygnss_time = cached['times_full']
                else:
                    cyg_nc = nc.Dataset(cyg_file_path)
                    longitudes = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]
                    # Pick source-specific wind product when available.
                    # L2 v3.2 defaults to FDS in `wind_speed`, but we prefer YSLF.
                    source_aliases = source_var_aliases.get(source_name, {})
                    wind_candidates = source_aliases.get('wind_speed', []) + ['wind_speed']
                    wind_var_name = None
                    for candidate in wind_candidates:
                        if candidate in cyg_nc.variables:
                            wind_var_name = candidate
                            break
                    if wind_var_name is None:
                        raise KeyError('No wind speed variable found in CYGNSS file')
                    cygnss_wind = cyg_nc.variables[wind_var_name][:]
                    if 'track_id' in cyg_nc.variables:
                        cygnss_tracks = cyg_nc.variables['track_id'][:]
                    elif ('sv_num' in cyg_nc.variables) and ('spacecraft_num' in cyg_nc.variables):
                        sv_vals = np.asarray(np.ma.filled(cyg_nc.variables['sv_num'][:], -1), dtype=np.int64)
                        sc_vals = np.asarray(np.ma.filled(cyg_nc.variables['spacecraft_num'][:], -1), dtype=np.int64)
                        flat_tracks = np.empty(sv_vals.size, dtype=object)
                        sv_flat = sv_vals.ravel()
                        sc_flat = sc_vals.ravel()
                        for idx in range(sv_flat.size):
                            flat_tracks[idx] = (int(sv_flat[idx]), int(sc_flat[idx]))
                        cygnss_tracks = flat_tracks.reshape(sv_vals.shape)
                        
                    cygnss_lons = ((longitudes + 180) % 360) - 180
                    cygnss_time = cyg_nc.variables['sample_time'][:]
                    # parse sample_time units to get base datetime
                    date_format = "%Y-%m-%d %H:%M:%S"
                    cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                    # save full times array as numpy datetime64
                    times_full = np.datetime64(cyg_start_date) + cygnss_time.astype('timedelta64[s]')

                     # Cache all other variables we may need so we don't need an open netCDF handle later
                    vars_cache = {}
                    for varname in variables:
                        if varname in ('cyg_vmax_error', 'cyg_dist', 'time_since_ri'):
                            continue
                        try:
                            if varname in cyg_nc.variables:
                                var_data = cyg_nc.variables[varname][:]
                                # Only cache 1D variables for DataFrame compatibility
                                if np.asarray(var_data).ndim == 1:
                                    vars_cache[varname] = var_data
                                else:
                                    vars_cache[varname] = np.array([])
                            else:
                                # First apply source-specific aliases (e.g., L2 YSLF fields).
                                source_aliases = source_var_aliases.get(source_name, {})
                                source_candidates = source_aliases.get(varname, [])
                                source_found = None
                                for alias_name in source_candidates:
                                    if alias_name in cyg_nc.variables:
                                        source_found = alias_name
                                        break
                                if source_found is not None:
                                    var_data = cyg_nc.variables[source_found][:]
                                    if np.asarray(var_data).ndim == 1:
                                        vars_cache[varname] = var_data
                                    else:
                                        vars_cache[varname] = np.array([])
                                    continue

                                # For channelized L2 fields, average available channels
                                # (port/starboard/(zenith)) into one per-sample variable.
                                if varname in avg_alias_fields:
                                    alias_arrays = []
                                    for alias_name in avg_alias_fields[varname]:
                                        if alias_name in cyg_nc.variables:
                                            alias_vals = np.asarray(
                                                np.ma.filled(cyg_nc.variables[alias_name][:], np.nan),
                                                dtype=float,
                                            )
                                            alias_arrays.append(alias_vals)
                                    if len(alias_arrays) > 0:
                                        vars_cache[varname] = np.nanmean(np.stack(alias_arrays, axis=0), axis=0)
                                        continue

                                aliases = var_aliases.get(varname, [])
                                alias_found = None
                                for alias_name in aliases:
                                    if alias_name in cyg_nc.variables:
                                        alias_found = alias_name
                                        break
                                if alias_found is not None:
                                    var_data = cyg_nc.variables[alias_found][:]
                                    if np.asarray(var_data).ndim == 1:
                                        vars_cache[varname] = var_data
                                    else:
                                        vars_cache[varname] = np.array([])
                                else:
                                    # missing variable in this file, store empty array
                                    vars_cache[varname] = np.array([])
                        except Exception:
                            # On unexpected read errors, store empty array for alignment.
                            vars_cache[varname] = np.array([])

                    cyg_nc.close()
                    cached = {
                        'lons': cygnss_lons,
                        'lats': cygnss_lats,
                        'wind': cygnss_wind,
                        'tracks': cygnss_tracks, # pyright: ignore[reportPossiblyUnboundVariable]
                        'times_full': times_full,
                        'vars': vars_cache
                    }
                    source_cache[cyg_file_path] = cached
                    noaa_cache[source_name] = source_cache
                    cygnss_lons = cached['lons']
                    cygnss_lats = cached['lats']
                    cygnss_wind = cached['wind']
                    cygnss_tracks = cached['tracks']
                    cygnss_time = cached['times_full']

                # find where the TC is for the date (box selection)
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = np.where(lats_boolean & lons_boolean)[0]
                if len(box_indices)>0:
                    meas_lats = cygnss_lats[box_indices]  # Measurement latitudes
                    meas_lons = cygnss_lons[box_indices]  # Measurement longitudes
                    meas_winds = cygnss_wind[box_indices]
                    # temporal matching - time within x number of minutes
                    meas_times = cygnss_time[box_indices]

                    # Build a single quality mask for L2 v3.2 YSLF values.
                    if source_name == 'l2_3.2':
                        quality_mask = np.ones(len(box_indices), dtype=bool)
                        cached_vars = cached.get('vars', {})

                        if yslf_max_wind_uncertainty is not None:
                            meas_unc = _extract_box_float_var(cached_vars, 'wind_speed_uncertainty', box_indices)
                            if meas_unc.size > 0:
                                quality_mask &= np.isfinite(meas_unc) & (meas_unc <= yslf_max_wind_uncertainty)

                        if yslf_incidence_angle_range is not None:
                            inc_min, inc_max = yslf_incidence_angle_range # type: ignore
                            meas_inc = _extract_box_float_var(cached_vars, 'incidence_angle', box_indices)
                            if meas_inc.size > 0:
                                quality_mask &= np.isfinite(meas_inc) & (meas_inc >= inc_min) & (meas_inc <= inc_max)
                        
                        meas_snr = _extract_box_float_var(cached_vars, 'snr', box_indices)
                        meas_rcg = _extract_box_float_var(cached_vars, 'range_corr_gain', box_indices)

                        if yslf_min_snr is not None:
                            if meas_snr.size > 0:
                                quality_mask &= np.isfinite(meas_snr) & (meas_snr >= yslf_min_snr)

                        if snr_rx_filter:
                            # Remove low-quality signal/gain combinations:
                            # (Rx <= 3 and SNR < 9) or (SNR < 1 and Rx < 7) ((meas_rcg <= 50) & (meas_snr < 1.4))|((meas_snr < 1.2) & (meas_rcg < 100)) 
                            reject_rx_snr = (meas_rcg < 50)
                            quality_mask &= ~reject_rx_snr

                        if yslf_range_corr_gain_range is not None:
                            rcg_min, rcg_max = yslf_range_corr_gain_range # type: ignore
                            if meas_rcg.size > 0:
                                quality_mask &= np.isfinite(meas_rcg) & (meas_rcg >= rcg_min) & (meas_rcg <= rcg_max)

                        if yslf_min_num_ddms_utilized is not None:
                            meas_nddm = _extract_box_float_var(cached_vars, 'num_ddms_utilized', box_indices)
                            if meas_nddm.size > 0:
                                quality_mask &= np.isfinite(meas_nddm) & (meas_nddm >= yslf_min_num_ddms_utilized)

                        if not np.any(quality_mask):
                            continue

                        box_indices = box_indices[quality_mask]
                        meas_lats = meas_lats[quality_mask]
                        meas_lons = meas_lons[quality_mask]
                        meas_winds = meas_winds[quality_mask]
                        meas_times = meas_times[quality_mask]
                                            
                    # Convert to NumPy arrays if not already
                    meas_coords = np.column_stack((meas_lats, meas_lons))  # Shape (N, 2)
                    tc_coords = np.column_stack((tc_lat_interpolated, tc_lon_interpolated))  # Shape (M, 2)

                    # Use BallTree (haversine) to find nearby track points per measurement (memory-efficient)
                    # filtered_distances = np.array([])
                    # nearest_indices = np.array([], dtype=int)
                    use_fallback_distance_matrix = False
                    try:
                        from sklearn.neighbors import BallTree
                        # raise ValueError("Force fallback to avoid sklearn dependency for now")  # Remove this line when BallTree is available
                        meas_rad = np.radians(meas_coords)
                        track_rad = np.radians(tc_coords)
                        tree = BallTree(track_rad, metric='haversine')
                        radius = MAX_DISTANCE / 6371.0  # convert km to radians on unit sphere

                        neighbors = tree.query_radius(meas_rad, r=radius)

                        selected_indices_list = []
                        filtered_distances_list = []
                        nearest_indices_list = []
                        for mi, neigh in enumerate(neighbors):
                            if neigh.size == 0:
                                continue
                            # Check temporal condition among neighbor track indices
                            neigh_times = tc_time_interpolated[neigh]
                            time_diffs = np.abs(meas_times[mi] - neigh_times)
                            valid_time_mask = time_diffs <= np.timedelta64(noaa_time_threshold, 'h')
                            if not np.any(valid_time_mask):
                                continue
                            valid_neigh = neigh[valid_time_mask]
                            # Compute distances to valid neighbours (small arrays) and pick the nearest
                            meas_pt = np.array([[meas_lats[mi], meas_lons[mi]]])
                            neigh_pts = np.column_stack((tc_lat_interpolated[valid_neigh], tc_lon_interpolated[valid_neigh]))
                            dists = haversine_vectorized(meas_pt, neigh_pts).flatten()
                            idx_min = np.argmin(dists)
                            selected_indices_list.append(mi)
                            filtered_distances_list.append(dists[idx_min])
                            nearest_indices_list.append(valid_neigh[idx_min])

                        if len(selected_indices_list) > 0:
                            selected_indices = np.array(selected_indices_list, dtype=int)
                            filtered_distances = np.array(filtered_distances_list)
                            nearest_indices = np.array(nearest_indices_list, dtype=int)
                        else:
                            selected_indices = np.array([], dtype=int)
                            filtered_distances = np.array([])
                            nearest_indices = np.array([], dtype=int)

                    except Exception:
                        # Fallback: original full pairwise calculation (slower)
                        use_fallback_distance_matrix = True
                        distances = haversine_vectorized(meas_coords, tc_coords)
                        time_diffs = np.abs(meas_times[:, None] - tc_time_interpolated[None, :])
                        valid_mask = np.any((distances <= MAX_DISTANCE) & (time_diffs <= np.timedelta64(noaa_time_threshold, 'h')), axis=1)
                        selected_indices = np.where(valid_mask)[0]                              

                    # Select valid measurements
                    if len(selected_indices) > 0:
                        cyg_lat_AOI = meas_lats[selected_indices] ; cyg_lon_AOI = meas_lons[selected_indices] ; cyg_times_AOI = meas_times[selected_indices]
                        cyg_winds_AOI = meas_winds[selected_indices] ; cyg_tracks_AOI = cygnss_tracks[box_indices][selected_indices] 
                        cyg_winds_vals = np.ma.filled(cyg_winds_AOI, np.nan).astype(float)
                        if use_fallback_distance_matrix:
                            filtered_distances = np.min(distances, axis=1)[selected_indices]
                            nearest_indices = np.argmin(distances, axis=1)[selected_indices]
                        cyg_vmax_error_vals = cyg_winds_vals - np.asarray(tc_vmax_interpolated[nearest_indices], dtype=float)
                        

                        cyg_lats.extend(np.asarray(cyg_lat_AOI).tolist())
                        cyg_lons.extend(np.asarray(cyg_lon_AOI).tolist())
                        cyg_tracks.extend(np.asarray(cyg_tracks_AOI).tolist())
                        cyg_winds.extend(cyg_winds_vals.tolist())
                        cyg_times.extend(np.asarray(cyg_times_AOI).tolist())
                        cyg_dist.extend(np.asarray(filtered_distances).tolist())
                        cyg_vmax_error.extend(cyg_vmax_error_vals.tolist())
                        coincident_meas['cyg_source'].extend([source_name] * len(cyg_lat_AOI))
                        
                        if store_all:
                            cyg_lats_all.extend(np.asarray(meas_lats).tolist()) # pyright: ignore[reportPossiblyUnboundVariable]
                            cyg_lons_all.extend(np.asarray(meas_lons).tolist()) # pyright: ignore[reportPossiblyUnboundVariable]
                            cyg_tracks_all.extend(np.asarray(cygnss_tracks[box_indices]).tolist()) # pyright: ignore[reportPossiblyUnboundVariable]
                            cyg_winds_all.extend(np.asarray(meas_winds).tolist())
                            cyg_times_all.extend(np.asarray(meas_times).tolist())
                            cyg_sources_all.extend([source_name] * len(meas_times))
                        for i in variables:
                            if i == 'cyg_dist':
                                var = np.atleast_1d(filtered_distances)
                            elif i == 'cyg_vmax_error':
                                var = np.atleast_1d(cyg_vmax_error_vals)
                            elif i == 'wind_speed':
                                var = np.atleast_1d(cyg_winds_vals)
                            elif i == 'time_since_ri':
                                var = np.atleast_1d(cyg_times_AOI - ri_start_time)
                            else:
                                # Prefer cached variable data (no open netCDF handle required)
                                var_full = cached.get('vars', {}).get(i, None)
                                if var_full is None or getattr(var_full, 'size', 0) == 0:
                                    # variable not present in file; fill with NaN to keep lengths aligned
                                    var = np.full(len(cyg_lat_AOI), np.nan)
                                else:
                                    var_indexed = np.asarray(var_full)[box_indices[selected_indices]]
                                    # Ensure 1D by taking first element if indexing created extra dimensions
                                    var = np.atleast_1d(var_indexed) if var_indexed.ndim == 1 else var_indexed.ravel()[:len(cyg_lat_AOI)]
                            if len(var) != len(cyg_lat_AOI):
                                continue
                            # extend the list if present
                            if var is not None:
                                coincident_meas[i].extend(var.tolist())
        
        # # update the wind_speed values for yslf, so if nbrcs_mean is 13 or less, then apply correction
        # if 'wind_speed' in coincident_meas and 'nbrcs_mean' in coincident_meas:
        #     wind_speeds = np.array(coincident_meas['wind_speed'], dtype=float)
        #     cyg_vmax_errors = np.array(coincident_meas.get('cyg_vmax_error', []), dtype=float)
        #     nbrcs_means = np.array(coincident_meas['nbrcs_mean'], dtype=float)
        #     correction_mask = (nbrcs_means <= 13) 
        #     yslf_mask = (coincident_meas['cyg_source'] == 'l2_3.2')
        #     # correction_mask &= yslf_mask
        #     correction = -3.5992 * nbrcs_means[correction_mask] + 49.388
        #     wind_speeds[correction_mask] = wind_speeds[correction_mask] - correction
        #     cyg_vmax_errors[correction_mask] = cyg_vmax_errors[correction_mask] + correction
        #     coincident_meas['wind_speed'] = wind_speeds.tolist()
        #     coincident_meas['cyg_vmax_error'] = cyg_vmax_errors.tolist()
        #     cyg_winds = wind_speeds
        #     cyg_vmax_error = cyg_vmax_errors

        # Convert accumulated lists into numpy arrays for vectorized operations
        if len(cyg_times) > 0:
            cyg_lats = np.array(cyg_lats)
            cyg_lons = np.array(cyg_lons)
            cyg_tracks = np.array(cyg_tracks, dtype=object)
            cyg_winds = np.array(cyg_winds)
            cyg_times = np.array(cyg_times, dtype='datetime64[s]')
            cyg_dist = np.array(cyg_dist)
            cyg_vmax_error = np.array(cyg_vmax_error)
            cyg_sources = np.asarray(coincident_meas.get('cyg_source', []), dtype=object)
            if len(cyg_sources) != len(cyg_times):
                cyg_sources = np.full(len(cyg_times), 'noaa_1.2', dtype=object)
            if cyg_match_var == 'cyg_dist':
                cyg_match = cyg_dist
            elif cyg_match_var == 'cyg_vmax_error':
                cyg_match = cyg_vmax_error
            elif cyg_match_var == 'wind_speed':
                cyg_match = cyg_winds
            elif cyg_match_var == 'time_since_ri':
                cyg_match = np.array(cyg_times - ri_start_time)
            else:
                cyg_match = np.asarray(coincident_meas.get(cyg_match_var, []))
                if len(cyg_match) != len(cyg_times):
                    cyg_match = np.full(len(cyg_times), np.nan)
            if store_all:
                cyg_lats_all = np.array(cyg_lats_all)
                cyg_lons_all = np.array(cyg_lons_all)
                cyg_tracks_all = np.array(cyg_tracks_all, dtype=object)
                cyg_winds_all = np.array(cyg_winds_all)
                cyg_times_all = np.array(cyg_times_all, dtype='datetime64[s]')
                cyg_sources_all = np.asarray(cyg_sources_all, dtype=object)
            # Decode all flags
            if 'sample_flags' in variables:
                decoded_flags_list = []
                for flag in coincident_meas['sample_flags']:
                    decoded = []
                    if flag & 1: decoded.append("poor_overall_quality")
                    if flag & 2: decoded.append("asc_node")
                    if flag & 4: decoded.append("block_IIF")
                    if flag & 8: decoded.append("block_IIR")
                    if flag & 16: decoded.append("block_IIRM")
                    if flag & 32: decoded.append("non_zero_nst")
                    if flag & 64: decoded.append("poor_quality_wind_nst")
                    if flag & 128: decoded.append("poor_quality_wind_sample")
                    decoded_flags_list.append(decoded)
                # Keep as a plain list so pandas treats it as a 1D object column.
                coincident_meas['sample_flags'] = decoded_flags_list

            tc_list.append(pd.DataFrame({'time_since_ri': time_since_ri,'wind_speed' :tc_vmax_interpolated}))
            # calculating vmax based on max velocity of each track - assuming not parralel to path
            ri_mask=((cyg_times>ri_start_time-np.timedelta64(noaa_time_threshold,'h')) & (cyg_times<ri_end_time+np.timedelta64(noaa_time_threshold,'h')))
            ri_lats=cyg_lats[ri_mask]; ri_lons=cyg_lons[ri_mask]; ri_winds=cyg_winds[ri_mask]; ri_times=cyg_times[ri_mask] ; ri_tracks=cyg_tracks[ri_mask]; ri_sources=cyg_sources[ri_mask]
            mid_time = ri_start_time + (ri_end_time - ri_start_time) // 2  # middle of the interpolated TC times
            # plot_timeline_noaa(ri_winds, ri_times, ri_dist,tc_time_interpolated, tc_vmax_interpolated, tc_storm_id, ri_start_time, ri_end_time)
            # Convert accumulated lists to numpy arrays for consistent dtypes before creating DataFrame
            # coincident_meas_arr = {}
            # for k, v in coincident_meas.items():
            #     if len(v) == 0:
            #         coincident_meas_arr[k] = np.array([], dtype=object)
            #     else:
            #         coincident_meas_arr[k] = np.array(v)
            
            # add ERA5 and CMS to coincident_meas before creating DataFrame
            # check if file already exists
            cms_check = False
            if cms_check:
                filepath_swh = f"E:\\Phd_data\\copernicus\\{tc_storm_id}_SWH.nc"
                filepath_sss_sst = f"E:\\Phd_data\\copernicus\\{tc_storm_id}_SSS&SST.nc"
                if not os.path.exists(filepath_swh) or not os.path.exists(filepath_sss_sst):
                    # print(f"Files for {tc_storm_id} already exist. Skipping download.")  
                    cms_downloader(AOI_lon_min,AOI_lon_max,AOI_lat_min,AOI_lat_max,cyg_AOI_time_start,cyg_AOI_time_end,tc_storm_id)

                try:
                    cms_res_SWH = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id,filepath_swh,['VHM0','VMDR'])
                except ValueError as e:
                    try:
                        expand_and_redownload(e, filepath_swh, tc_storm_id, cms_downloader, file_type='netcdf')
                        cms_res_SWH = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id,filepath_swh,['VHM0','VMDR'])
                    except ValueError as e2:
                        expand_and_redownload(e2, filepath_swh, tc_storm_id, cms_downloader, file_type='netcdf')
                        cms_res_SWH = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id,filepath_swh,['VHM0','VMDR'])
                    
                try:
                    cms_res_SSS = retrieve_vals(cyg_lats, cyg_lons, cyg_times,tc_storm_id,filepath_sss_sst,['thetao','so'])
                except ValueError as e:
                    try:
                        expand_and_redownload(e, filepath_sss_sst, tc_storm_id, cms_downloader, file_type='netcdf')
                        cms_res_SSS = retrieve_vals(cyg_lats, cyg_lons, cyg_times,tc_storm_id,filepath_sss_sst,['thetao','so'])
                    except ValueError as e2:
                        expand_and_redownload(e2, filepath_sss_sst, tc_storm_id, cms_downloader, file_type='netcdf')
                        cms_res_SSS = retrieve_vals(cyg_lats, cyg_lons, cyg_times,tc_storm_id,filepath_sss_sst,['thetao','so'])
                    
                # add results into coincident_meas dictionary
                for k, v in cms_res_SWH.items():
                    if k not in coincident_meas:
                        coincident_meas[k] = []
                    coincident_meas[k].extend(v.tolist())

                for k, v in cms_res_SSS.items():
                    if k not in coincident_meas:
                        coincident_meas[k] = []
                    coincident_meas[k].extend(v.tolist())                        
                
                # check if file already exists
                import os
                filepath_era5 = f"E:\\Phd_data\\copernicus\\{tc_storm_id}_era5.grib"
                if not os.path.exists(filepath_era5):
                    # print(f"ERA5 file for {tc_storm_id}  already exists. Skipping download.")
                    era5_downloader(AOI_lon_min,AOI_lon_max,AOI_lat_min,AOI_lat_max,cyg_AOI_time_start,cyg_AOI_time_end,tc_storm_id)
                era5_vars = ['u10','v10','t2m','sp'] # 'msl',
                try:
                    era5_res = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id, filepath_era5, era5_vars)
                except ValueError as e:
                    try:
                        expand_and_redownload(e, filepath_era5, tc_storm_id, era5_downloader, file_type='grib')
                        era5_res = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id, filepath_era5, era5_vars)
                    except ValueError as e2:
                        expand_and_redownload(e2, filepath_era5, tc_storm_id, era5_downloader, file_type='grib')
                        era5_res = retrieve_vals(cyg_lats, cyg_lons, cyg_times, tc_storm_id, filepath_era5, era5_vars)

                era5_res['era5_ws'] = np.hypot(era5_res['u10'], era5_res['v10'])
                # era5_res['era5_dir'] = (np.arctan2(-era5_res['u10'], -era5_res['v10']) * 180 / np.pi) % 360
                # era5_res shape: (n_measurements, n_vars) -> iterate over columns
                for k, v in era5_res.items():
                    if k not in coincident_meas:
                        coincident_meas[k] = []
                    coincident_meas[k].extend(v.tolist())

            # Normalize arrays to 1D to avoid pandas construction errors.
            for k, v in list(coincident_meas.items()):
                if isinstance(v, np.ndarray) and v.ndim > 1:
                    if v.shape[-1] == 1:
                        coincident_meas[k] = v.reshape(-1).tolist()
                    else:
                        coincident_meas[k] = [np.asarray(row).tolist() for row in v]
            
            # After processing all files and accumulating measurements, analyze tracks within the RI window
            event_df = pd.DataFrame(coincident_meas)
            ri_track_keys = [_track_key(v) for v in ri_tracks]
            unique_tracks = list(dict.fromkeys(ri_track_keys))
            num_tracks = len(unique_tracks)
            if num_tracks > 1: # so 2 or more pass over
                # cyg_lat_AOI = selected_lats[valid_time] ; cyg_lon_AOI = selected_lons[valid_time] ; cyg_times_AOI = cyg_datetime_array[valid_time]
                # cyg_winds_AOI =cygnss_wind[box_indices][indices][valid_time]                
                # Arrays to store results
                vmax_time = np.zeros(num_tracks, dtype=object)
                max_wind = np.full(num_tracks, np.nan, dtype=float)
                track_sources = np.empty(num_tracks, dtype=object)
                valid_track_peak = np.zeros(num_tracks, dtype=bool)
                peak_global_indices = np.full(num_tracks, -1, dtype=int)
                ri_row_indices = np.where(ri_mask)[0]
                
                # Process each track 
                for i, track in enumerate(unique_tracks):
                    track_mask = np.array([k == track for k in ri_track_keys], dtype=bool)  # Proper boolean mask for track filtering
                    track_times = ri_times[track_mask]  # Filter timestamps based on track mask
                    track_winds = ri_winds[track_mask]  # Filter wind speeds based on track mask
                
                    # Ensure that we are not running into errors if the filtered array is empty
                    if track_times.size > 0:
                        track_winds_vals = np.asarray(np.ma.filled(track_winds, np.nan), dtype=float)
                        if not np.any(np.isfinite(track_winds_vals)):
                            continue

                        # Store the values
                        vmax_index = int(np.nanargmax(track_winds_vals))
                        vmax_time[i] = track_times[vmax_index]
                        max_wind[i] = track_winds_vals[vmax_index]  # Maximum wind speed for the track
                        track_indices_local = np.where(track_mask)[0]
                        peak_local_index = int(track_indices_local[vmax_index])
                        peak_global_indices[i] = int(ri_row_indices[peak_local_index])
                        track_src = ri_sources[track_mask]
                        if len(track_src) > 0:
                            unique_src, src_counts = np.unique(track_src, return_counts=True)
                            track_sources[i] = unique_src[np.argmax(src_counts)]
                        valid_track_peak[i] = True
                
                # remove any nan values (if any track had no measurements after filtering)
                # separate based on source
                noaa_track_mask = (track_sources == 'noaa_1.2') & valid_track_peak
                l2_track_mask = (track_sources == 'l2_3.2') & valid_track_peak
                noaa_track_count = int(np.sum(noaa_track_mask))
                l2_track_count = int(np.sum(l2_track_mask))
                vmax_time_noaa = vmax_time[noaa_track_mask] ; vmax_time_l2 = vmax_time[l2_track_mask]
                max_wind_noaa = max_wind[noaa_track_mask] ; max_wind_l2 = max_wind[l2_track_mask]

                # save track-level results into dataframe
                noaa_peak_rows = peak_global_indices[noaa_track_mask]
                l2_peak_rows = peak_global_indices[l2_track_mask]
                noaa_df_list.append(event_df.iloc[noaa_peak_rows].reset_index(drop=True))
                l2_df_list.append(event_df.iloc[l2_peak_rows].reset_index(drop=True))


                # find different max wind if max is the first one
                # Find indices where time is later than earliest
                def valid_RI_meas(vmax_time,max_wind):
                    # No measurements for this source/tracks.
                    if len(vmax_time) == 0 or len(max_wind) == 0:
                        return (False, 0, 0)

                    # the initial vmax value is the maximum wind speed either within 3 hours of first measurement or within 3 hours of RI start time, whichever is earlier. This is to account for cases where the first measurement is during the RI period but there were no measurements before it, so we want to look for the earliest maximum within a reasonable window to get a better estimate of the initial intensity before RI.
                    earlist_track_time = np.nanmin(vmax_time)
                    valid_initial_mask = (vmax_time <= earlist_track_time + np.timedelta64(noaa_time_threshold,'h'))&(vmax_time <= ri_start_time + np.timedelta64(noaa_time_threshold,'h'))

                    # Fallback when RI timing filter is too strict for sparse overpasses.
                    if not np.any(valid_initial_mask):
                        valid_initial_mask = (vmax_time <= earlist_track_time + np.timedelta64(noaa_time_threshold,'h'))
                    if not np.any(valid_initial_mask):
                        return (False, 0, 0)

                    valid_initial_indices = np.where(valid_initial_mask)[0]
                    inital_max_ind = valid_initial_indices[np.nanargmax(max_wind[valid_initial_indices])] # or swap with mean of tracks
                    vmax_time_arr = np.array(vmax_time, dtype='datetime64[s]')
                    # Calculate time differences in hours from the RI Vmax estimate to all other track maxima
                    time_diffs = np.round((vmax_time_arr - vmax_time_arr[inital_max_ind]) / np.timedelta64(1, 'h'))
                    valid_mask = (time_diffs >= tc_durations / 2) & (time_diffs <= tc_durations + 2)
                    
                    if np.any(valid_mask):
                        # Of valid indices, pick the one with the max wind
                        valid_indices = np.where(valid_mask)[0]
                        max_ind = valid_indices[np.argmax(max_wind[valid_indices])]
                        
                        wind_change = max_wind[max_ind] - max_wind[inital_max_ind]
                        time_diff = time_diffs[max_ind]
                        return (True, wind_change, time_diff)
                    return (False, 0, 0)
                l2_res = valid_RI_meas(vmax_time_l2, max_wind_l2)
                if l2_res[0]:    
                    l2_3p2_ri_change.append(ri_vmax_change)
                    l2_3p2_ri_initial.append(ri_vmax_start)
                    l2_3p2_wind_change.append(l2_res[1])
                    l2_3p2_ri_durations.append(tc_durations)
                noaa_res = valid_RI_meas(vmax_time_noaa, max_wind_noaa)
                if noaa_res[0]:
                    noaa_ri_change.append(ri_vmax_change)
                    noaa_ri_initial.append(ri_vmax_start)
                    noaa_wind_change.append(noaa_res[1])
                    noaa_ri_durations.append(tc_durations)
                
                    # print(tc_intensities,tc_durations,wind_change,time_diff, cyg_AOI_time_start, cyg_AOI_time_end, cyg_lats[0],cyg_lons[0],ri_rd[tc_event])

                # Initialize merged data variables for animation
                merged_vmax_lat = None
                merged_vmax_lon = None
                merged_datetimes = None
                merged_vmax_winds = None
                    

        # Merged 
        if not np.any(np.isin(date_list,  np.array(merged_file_dates))):
            continue
        # if cyg_version == 'storm_centric':
        valid_files = np.where(date_list[0] == np.array(merged_file_dates))[0]
        unique_storm_files=[]
        if len(valid_files)==1:
            unique_storm_files.append(merged_cyg_files_list[cyg_files_num[valid_files[0]]])
        else:
            for i in range(len(valid_files)-1):
                if merged_cyg_files_list[cyg_files_num[valid_files[i]]] != merged_cyg_files_list[cyg_files_num[valid_files[i+1]]]:
                    unique_storm_files.append(merged_cyg_files_list[cyg_files_num[valid_files[i]]])
        unique_storm_files = list(dict.fromkeys(unique_storm_files)) # remove duplicates while preserving order   
        found_already = False
        merged_valid_RI = False
    
        for cyg_file_path in unique_storm_files:
            with nc.Dataset(cyg_file_path) as cyg_nc:
                cyg_storm_lat = cyg_nc.variables['best_track_storm_center_lat'][:] ; cyg_storm_lon = ((cyg_nc.variables['best_track_storm_center_lon'][:] + 180) % 360) - 180
                cygnss_time =cyg_nc.variables['time'][:]
                cygnss_time_units =cyg_nc.variables['time'].units[12:]
                cyg_datetimes = np.datetime64(cygnss_time_units) + cygnss_time.astype('timedelta64[h]')
                # as long as at least 1 measurement falls within the AOI box during the RI period, we will consider it for comparison - to ensure we are not missing out on any potential matches due to the strict box selection
                if not (np.any((cyg_storm_lon > AOI_lon_min)&(cyg_storm_lon < AOI_lon_max) & (cyg_storm_lat > AOI_lat_min) & (cyg_storm_lat < AOI_lat_max)&(not found_already))):
                    continue
                found_already = True
                # print(cyg_file_path)
                epochs=np.where((cyg_storm_lon > AOI_lon_min)&(cyg_storm_lon < AOI_lon_max) & (cyg_storm_lat > AOI_lat_min) & (cyg_storm_lat < AOI_lat_max) & (cyg_datetimes >= ri_start_time- np.timedelta64(4, 'h')) & (cyg_datetimes <= ri_end_time+ np.timedelta64(4, 'h')))[0]
                if len(epochs)>1: # to ensure we have two or more measurements over RI period
                    vmax_lat = np.asarray(cyg_nc.variables['cygnss_vmax_lat'][epochs], dtype=float)
                    vmax_lon = np.asarray(((cyg_nc.variables['cygnss_vmax_lon'][epochs] + 180) % 360) - 180, dtype=float)
                    cyg_datetimes=cyg_datetimes[epochs]
                    longitudes  = cyg_nc.variables['lon'][:] ; cygnss_lats = np.asarray(cyg_nc.variables['lat'][:], dtype=float) ;
                    cygnss_lons = np.asarray(((longitudes + 180) % 360) - 180, dtype=float)
                    
                    merged_vmax = np.full(len(epochs), np.nan, dtype=float)
                    merged_vmax_time = np.full(len(epochs), np.datetime64('NaT'), dtype='datetime64[s]')

                    # finding the vmax vals at the vmax lat/lon for each epoch - this is to ensure we are comparing the same location for each epoch rather than just taking the maximum value within the AOI which may not be at the same location as the RI vmax and may not be consistent across epochs. This is especially important for storm-centric files where the vmax lat/lon can change significantly across epochs.
                    lat_res = np.nanmedian(np.abs(np.diff(cygnss_lats))) if len(cygnss_lats) > 1 else 0.25
                    lon_res = np.nanmedian(np.abs(np.diff(cygnss_lons))) if len(cygnss_lons) > 1 else 0.25
                    lat_tol = max(float(lat_res) * 0.6, 1e-4)
                    lon_tol = max(float(lon_res) * 0.6, 1e-4)
                    for j in range(len(epochs)):
                        if not np.isfinite(vmax_lat[j]) or not np.isfinite(vmax_lon[j]):
                            continue
                        lat_idx = int(np.nanargmin(np.abs(cygnss_lats - vmax_lat[j])))
                        lon_idx = int(np.nanargmin(np.abs(cygnss_lons - vmax_lon[j])))
                        if np.abs(cygnss_lats[lat_idx] - vmax_lat[j]) > lat_tol or np.abs(cygnss_lons[lon_idx] - vmax_lon[j]) > lon_tol:
                            continue
                        wind_val_epochs = cyg_nc.variables['wind_speed'][epochs, lat_idx, lon_idx]
                        wind_val = wind_val_epochs[j] 
                        merged_vmax[j] = np.nan if np.ma.is_masked(wind_val) else float(wind_val)

                        merged_vmax_time_offset_epochs = cyg_datetimes + np.array(cyg_nc.variables['time_offset'][epochs, lat_idx, lon_idx], dtype='timedelta64[s]')
                        merged_vmax_time[j] = merged_vmax_time_offset_epochs[j]
                        

                    # Find the epoch that corresponds to the RI start time (or the closest one within a reasonable window if an exact match is not found)
                    # find the closest epoch to the RI start time within a reasonable window (e.g., 3 hours)
                    time_diffs = np.abs(merged_vmax_time - ri_start_time)
                    time_diffs_hours = time_diffs / np.timedelta64(1, 'h')
                    valid_time_mask = (time_diffs_hours <= noaa_time_threshold) & (merged_vmax > 0)  # Only consider epochs with valid wind speeds
                    if not np.any(valid_time_mask):
                        continue
                    start_epoch = int(np.where(valid_time_mask)[0][np.argmin(time_diffs_hours[valid_time_mask])])
                

                    # Find epochs that are valid at or after tc_duration/2.
                    # Use a masked-array-safe conversion to elapsed hours.
                    time_deltas = np.ma.asarray(merged_vmax_time - cyg_datetimes[start_epoch])
                    time_diffs = np.ma.filled(time_deltas / np.timedelta64(1, 'h'), np.nan).astype(float)
                    valid_mask = (time_diffs >= tc_durations / 2) & (time_diffs <= tc_durations + 2) & (merged_vmax > 0)  # Only consider epochs with valid wind speeds
                    valid_indices = np.where(valid_mask)[0]
                    if len(valid_indices) == 0:
                        continue
                    end_epoch = valid_indices[np.argmax(merged_vmax[valid_indices])]
                    change_vmax = float(merged_vmax[end_epoch] - merged_vmax[start_epoch])

                    merged_ri_change.append(float(ri_vmax_change))
                    merged_ri_initial.append(float(ri_vmax_start))
                    merged_wind_change.append(change_vmax)
                    merged_ri_durations.append(float(tc_durations))
                    merged_valid_RI = True

                    merged_list.append(pd.DataFrame({'time_since_ri': merged_vmax_time-ri_start_time,'wind_speed' :merged_vmax}))
                    # tc_list.append(pd.DataFrame({'time_since_ri': time_since_ri,'wind_speed' :tc_vmax_interpolated}))

                    # Store merged data for animation
                    merged_vmax_lat = vmax_lat # cygnss_lats
                    merged_vmax_lon = vmax_lon # cygnss_lons
                    merged_datetimes = merged_vmax_time
                    merged_vmax_winds = merged_vmax # cygnss_wind 

        if case_study_ind is not None:   
            # Time series plot
            plot_timeline_both(ri_winds, ri_times, cyg_match, merged_vmax_winds, merged_datetimes, tc_time_interpolated, tc_vmax_interpolated, tc_storm_id+'_'+str(tc_event), ri_start_time, ri_end_time, cyg_sources=ri_sources)
        # Call animation function after all data processing
            # Creates 4 GIFs: NOAA, YSLF, Merged, and Combined
            plot_tc_animation(
                tc_lat_interpolated, tc_lon_interpolated, tc_time_interpolated, tc_vmax_interpolated,
                cyg_lats, cyg_lons, cyg_times, cyg_winds, cyg_sources,
                cyg_lats_all, cyg_lons_all, cyg_times_all, cyg_sources_all,
                merged_vmax_lat=merged_vmax_lat, merged_vmax_lon=merged_vmax_lon,
                merged_datetimes=merged_datetimes, merged_vmax_winds=merged_vmax_winds,
                tc_storm_id=tc_storm_id, tc_durations=tc_durations, tc_intensities=tc_intensities,
                tc_event=tc_event, ri_start_time=ri_start_time, ri_end_time=ri_end_time,
                mid_time=mid_time, comp=comp, MAX_DISTANCE=MAX_DISTANCE
            )
            output_dir = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Plots', 'case_studies')
            os.makedirs(output_dir, exist_ok=True)
            noaa_df = pd.DataFrame(noaa_df_list[0])
            noaa_csv = os.path.join(output_dir, f'{tc_event}_NOAA_measurements.csv')
            noaa_df.to_csv(noaa_csv, index=False)
            yslf_df = pd.DataFrame(l2_df_list[0])
            yslf_csv = os.path.join(output_dir, f'{tc_event}_YSLF_measurements.csv')
            yslf_df.to_csv(yslf_csv, index=False)
            if merged_valid_RI == True:
                merged_data = {
                    
                    'latitude': merged_vmax_lat,
                    'longitude': merged_vmax_lon,
                    'datetime': merged_datetimes,
                    'wind_speed': merged_vmax_winds
                }
                merged_df = pd.DataFrame(merged_data)
                merged_csv = os.path.join(output_dir, f'{tc_event}_Merged_measurements.csv')
                merged_df.to_csv(merged_csv, index=False)

        
        
# save the results
if (not skip_event_processing) and save_checkpoint:
    checkpoint_payload = {
        'results_dt': results_dt,
        'noaa_df_list': noaa_df_list,
        'l2_df_list': l2_df_list,

        'tc_list': tc_list,
        'merged_list': merged_list,
        'noaa_ri_change': noaa_ri_change,
        'merged_ri_change': merged_ri_change,
        'noaa_wind_change': noaa_wind_change,
        'merged_wind_change': merged_wind_change,
        'noaa_ri_durations': noaa_ri_durations,
        'merged_ri_durations': merged_ri_durations,
        'noaa_ri_initial': noaa_ri_initial,
        'merged_ri_initial': merged_ri_initial,
        'l2_3p2_ri_change': l2_3p2_ri_change,
        'l2_3p2_wind_change': l2_3p2_wind_change,
        'l2_3p2_ri_durations': l2_3p2_ri_durations,
        'l2_3p2_ri_initial': l2_3p2_ri_initial,
    }
    with open(checkpoint_file, 'wb') as f:
        pickle.dump(checkpoint_payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved checkpoint: {checkpoint_file}")
        
#%%
def product_skill_scatterplot(ri_change,cyg_wind_change,ri_initial, ri_durations, filename, shortname):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    if len(ri_change)==0:
        print(f"No data to plot for {shortname}. Skipping scatterplot.")
        return
    cyg_wind_change = np.array(cyg_wind_change) ; ri_change = np.array(ri_change) ; ri_initial = np.array(ri_initial) ; ri_durations = np.array(ri_durations)
    cyg_change_error = cyg_wind_change - ri_change

    # number of points for each duration
    unique_durations, durations_count = np.unique(ri_durations, return_counts=True)
    if len(unique_durations) < 3:
        # add 6 to start of unique_durations and 0 to durations_count to ensure we have 3 categories for the legend
        unique_durations = np.insert(unique_durations, 0, 6)
        durations_count = np.insert(durations_count, 0, 0)
    idx = np.argsort(unique_durations)[::-1]
    unique_durations = unique_durations[idx]
    durations_count = durations_count[idx] 
    colours = ['blue', 'orange', 'green']

    # filter for low initial speeds
    lower = 15.3 ; upper = 23.1
    in_box      = ((ri_change >= lower) & (ri_change <= upper)) & ((cyg_wind_change >= lower) & (cyg_wind_change <= upper))
    x_high      = (ri_change > upper) & ((cyg_wind_change >= lower) & (cyg_wind_change <= upper))
    y_high      = ((ri_change >= lower) & (ri_change <= upper)) & (cyg_wind_change > upper)
    xy_high     = (ri_change > upper) & (cyg_wind_change > upper)
    x_low       = (ri_change < lower) & ((cyg_wind_change >= lower) & (cyg_wind_change <= upper))
    y_low       = ((ri_change >= lower) & (ri_change <= upper)) & (cyg_wind_change < lower) & (cyg_wind_change >= 0) # ensure no negative values in y_low
    y_low_x_high = ((ri_change >= upper) & (cyg_wind_change < lower) & (cyg_wind_change >= 0)) # ensure no negative values in y_low_x_high
    xy_low      = (ri_change < lower) & (cyg_wind_change < lower) & (cyg_wind_change >= 0) # ensure no negative values in xy_low
    y_superlow_x_low = (ri_change < upper) & (cyg_wind_change < 0)
    y_superlow_x_high = (ri_change > upper) & (cyg_wind_change < 0)

    plt.figure(figsize=(8, 6))
    for i, duration in enumerate(unique_durations):
        mask = ri_durations == duration
        plt.scatter(ri_change[mask], cyg_wind_change[mask], c=colours[i],  s=60, alpha=0.7, marker='o', label=f'{duration} hr ({durations_count[i]}/{num_events_duration[i]} events)')

    plt.legend(title='Durations', loc='upper right')   # add legend

    plt.axvline(x=lower, color='grey', linestyle='--')
    plt.axvline(x=upper, color='grey', linestyle='--')
    plt.axhline(y=lower, color='grey', linestyle='--')
    plt.axhline(y=upper, color='grey', linestyle='--')
    plt.axhline(y=0, color='grey', linestyle='--')

    plt.xlim(14,42) ; plt.ylim(-11, 39)
    # Add middle green box (in_box)
    plt.gca().add_patch(patches.Rectangle((15.3, 15.3), 23.1-15.3, 23.1-15.3,linewidth=0, facecolor='green', alpha=0.2))
    plt.gca().add_patch(patches.Rectangle((23.1, 23.1), plt.xlim()[1]-23.1, plt.ylim()[1]-23.1,linewidth=0, facecolor='green', alpha=0.2))
    total_15 = np.sum(in_box) + np.sum(y_high) + np.sum(y_low) + np.sum(y_superlow_x_low)  ; total_23 = np.sum(y_low_x_high) + np.sum(xy_high) + np.sum(x_high) + np.sum(y_superlow_x_high)

    # Adding text annotations
    text = True
    if text:
        plt.text(16, 16, f"{np.sum(in_box)} pts  {np.sum(in_box)/total_15*100:.1f}%", fontsize=12, color='black')      # Middle box
        plt.text(23.5, 16, f"{np.sum(x_high)} pts  {np.sum(x_high)/total_23*100:.1f}%", fontsize=12, color='black')       # Right
        plt.text(16, 25.5, f"{np.sum(y_high)} pts  {np.sum(y_high)/total_15*100:.1f}%", fontsize=12, color='black')       # Above
        plt.text(23.5, 25.5, f"{np.sum(xy_high)} pts  {np.sum(xy_high)/total_23*100:.1f}%", fontsize=12, color='black')    # Top-right
        plt.text(16, 8, f"{np.sum(y_low)} pts  {np.sum(y_low)/total_15*100:.1f}%", fontsize=12, color='black')           # Below
        plt.text(23.5, 8, f"{np.sum(y_low_x_high)} pts  {np.sum(y_low_x_high)/total_23*100:.1f}%", fontsize=12, color='black')  # Lower-Right
        plt.text(16, -8, f"{np.sum(y_superlow_x_low)} pts {np.sum(y_superlow_x_low)/total_15*100:.1f}%", fontsize=12, color='black')           # Below x-axis
        plt.text(23.5, -8, f"{np.sum(y_superlow_x_high)} pts {np.sum(y_superlow_x_high)/total_23*100:.1f}%", fontsize=12, color='black')           # Below x-axis
        shortname = shortname + "_numbers"

        # add bias, RMSE, and R2 values to the plot
        bias = np.mean(cyg_change_error)
        rmse = np.sqrt(np.mean(cyg_change_error**2))
        from sklearn.metrics import r2_score
        r2 = r2_score(ri_change, cyg_wind_change)
        plt.text(16, -15, f'Bias: {bias:.2f} m/s\nRMSE: {rmse:.2f} m/s\nR²: {r2:.2f}', fontsize=10, color='black')

    plt.xlabel('Best Track Vmax change (m/s)')
    plt.ylabel('CYGNSS Vmax change (m/s)')
    plt.title(shortname)
    plt.savefig(f'{filename}_{shortname}.png' , dpi=500)
    plt.show(block=False)
    plt.pause(1)
    plt.close()

    
if case_study_ind is None:
    # find the total number of events for each duration 
    unique_durations, num_events_duration  = np.unique(durations, return_counts=True)
    # swap the order of unique_durations and num_events_duration so that they are in descending order of duration
    idx = np.argsort(unique_durations)[::-1]
    unique_durations = unique_durations[idx]
    num_events_duration = num_events_duration[idx]

    # create scatterplots for each product
    
    l2_ri_change = l2_3p2_ri_change
    l2_wind_change = l2_3p2_wind_change
    l2_ri_initial = l2_3p2_ri_initial
    l2_ri_durations = l2_3p2_ri_durations

    gif_file = r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Plots\skill_scatter'
    product_skill_scatterplot(noaa_ri_change,noaa_wind_change,noaa_ri_initial,noaa_ri_durations, gif_file, 'NOAA')
    product_skill_scatterplot(l2_ri_change,l2_wind_change,l2_ri_initial,l2_ri_durations, gif_file, 'YSLF')
    product_skill_scatterplot(merged_ri_change,merged_wind_change,merged_ri_initial, merged_ri_durations, gif_file, 'Merged')
    
    # # figure for error in change vs initial speed
    # plt.figure(figsize=(8, 6))
    # plt.scatter(noaa_ri_initial,cyg_change_error,c=ri_change, cmap='viridis', edgecolor='black', s=60, alpha=0.7, marker='o', vmin=0, vmax=40)
    # plt.xlabel('Initial Vmax (m/s)')
    # plt.ylabel('NOAA - IBTrACS Vmax change error (m/s)')
    # # insert a trendline 
    # z = np.polyfit(noaa_ri_initial, cyg_change_error, 1)
    # p = np.poly1d(z)
    # plt.plot(noaa_ri_initial,p(noaa_ri_initial),"r--")
    # cbar = plt.colorbar()
    # cbar.set_label('Vmax change over 24hr (m/s)')
    # plt.savefig(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Plots\skill_iniital_2026.png' , dpi=500)
    # plt.show()

    results_dt.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\ri_results_v5.csv')
    final_noaa_df = pd.concat(noaa_df_list, ignore_index=True)
    final_l2_df = pd.concat(l2_df_list, ignore_index=True)
    final_tc_df = pd.concat(tc_list, ignore_index=True)
    final_merged_df =pd.concat(merged_list, ignore_index=True)
    # final_noaa_df.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\NOAA_Measurements.csv')
    # final_l2_df.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\YSLF_Measurements.csv')
    # final_merged_df.to_csv(r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\Merged_Measurements.csv')
    plot_distribution(final_tc_df,final_noaa_df,final_merged_df,final_l2_df)

# Make a series of scatterplots showing the relationship between the error in RI change (CYGNSS - Best Track) to each of the variables in the final_df. This will help to identify if there are any relationships between the error and the variables, and if the initial Vmax has any influence on this relationship.
charts = False
if charts:
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    comp = 'ashle'
    for product in ['NOAA','YSLF']: # 
        df = pd.read_csv(r'C:\Users\{0}\OneDrive - RMIT University\PHD\Data\IBTrACS\{1}_Measurements_RCG35.csv'.format(comp,product))
        df = df.dropna(subset=['cyg_vmax_error'])
        for col in df.columns:
            if col not in ['Index', 'cyg_vmax_error', 'sample_flags', 'time_since_ri','ddm_sample_index','ddm_channel','VMDR','u10','v10','era5_dir','cyg_source']:
                print(col)
                plt.figure(figsize=(8, 6))
                sns.kdeplot(
                    data=df,
                    x=col,
                    y='cyg_vmax_error',
                    fill=True,
                    levels=30,
                    thresh=0.02,
                    cmap='viridis')
                plt.xlabel(col)
                plt.ylabel('CYGNSS - Best Track Vmax change error (m/s)')
                plt.title(f'{product} Error vs {col}')
                # add a quadratic trendline to the plot
                z_poly = np.polyfit(df[col], df['cyg_vmax_error'], 2)
                z_lin = np.polyfit(df[col], df['cyg_vmax_error'], 1)
                p = np.poly1d(z_poly)
                p_lin = np.poly1d(z_lin)

                x_sorted = np.sort(df[col])
                plt.plot(x_sorted, p(x_sorted), "r--", linewidth=2)
                plt.plot(x_sorted, p_lin(x_sorted), "b--", linewidth=2)
                # add R2 value for the goodness of fit of the trendline
                from sklearn.metrics import r2_score
                r2 = r2_score(df['cyg_vmax_error'], p(df[col]))
                plt.text(0.05, 0.97, f'R2 = {r2:.2f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
                # show best-fit line equation on the plot
                x='x'
                plt.text(0.05, 0.05, f'Error = {z_poly[0]:.6e}*{x}^2 + {z_poly[1]:.6e}*{x} + {z_poly[2]:.6e}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
                plt.text(0.05, 0.1, f'Error = {z_lin[0]:.6e}*{x} + {z_lin[1]:.6e}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')

                plt.savefig(r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\{1}_error_vs_{2}.png'.format(comp, product, col), dpi=500)
                plt.show(block=False)
                plt.pause(1)
                plt.close()

bootstrapping=False
if bootstrapping:
    import xgboost as xgb
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    # Check if df exists, if not create it by concatenating noaa_df_list and save to csv
    comp = 'ashle'
    for product in ['NOAA','YSLF']:
        df_filename = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Data\IBTrACS\{1}_Measurements_RCG35.csv'.format(comp, product)
        df = pd.read_csv(df_filename)

        # remove any rows where Y is empty or NaN
        df = df.dropna(subset=['cyg_vmax_error'])

        y = df['cyg_vmax_error']
        X = df.drop(columns=['Index', 'cyg_vmax_error', 'wind_speed_uncertainty', 'sample_flags', 'time_since_ri','ddm_sample_index','ddm_channel','VMDR','u10','v10','Unnamed: 0'], errors='ignore')
        X = X.select_dtypes(include=['number']).copy()
        X = X.dropna(axis=1, how='all')

        feature_label_map = {
            'incidence_angle': 'Incidence angle',
            'wind_speed': 'CYGNSS wind speed',
            'rx_gain': 'Receiver gain',
            'snr': 'SNR',
            'range_corr_gain': 'Range-corrected gain',
            'num_ddms_utilized': '# DDMs used',
            'nbrcs_mean_corrected': 'Avg.NBRCS corrected',
            'nbrcs_mean': 'Avg.NBRCS',
            'cyg_dist': 'Distance to best-track',
            'thetao': 'SST',
            'so': 'Salinity',
            'msl': 'MSLP',
            'sp': 'Surface Pressure',
            'VHM0': 'Significant wave height',
            'era5_ws': 'ERA5 Wind Speed',
            't2m': '2m Air Temp.',
            'SST': 'SST',
            'SSS': 'Salinity',
            'mean_sea_level_pressure': 'MSLP',
            '2m_temperature': '2m Air Temp.',
        }
        X = X.rename(columns=feature_label_map)
        feature_names = X.columns

        if X.shape[1] == 0:
            raise ValueError('No numeric predictor columns available for bootstrapping after filtering.')
        # Parameters
        n_bootstrap = 1000
        n_features = X.shape[1]
        importances_bootstrap = np.zeros((n_bootstrap, n_features))
        use_gpu = True
        
        # Bootstrap loop
        for i in range(n_bootstrap):
            sample_indices = np.random.choice(len(X), size=len(X), replace=True)
            X_sample = X.iloc[sample_indices].to_numpy(dtype=np.float32)
            y_sample = y.iloc[sample_indices].to_numpy(dtype=np.float32)

            model = xgb.XGBRegressor(
                n_estimators=100,
                random_state=i,
                tree_method='hist',
                device='cuda' if use_gpu else 'cpu',
                n_jobs=1,
            )
            try:
                model.fit(X_sample, y_sample)
            except xgb.core.XGBoostError:
                if use_gpu:
                    use_gpu = False
                    model = xgb.XGBRegressor(
                        n_estimators=100,
                        random_state=i,
                        tree_method='hist',
                        device='cpu',
                        n_jobs=1,
                    )
                    model.fit(X_sample, y_sample)
                else:
                    raise
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

        # Save importance_df to csv
        importance_df.to_csv(r'C:\Users\{0}\OneDrive - RMIT University\PHD\Data\IBTrACS\{1}_feature_importance_bootstrap.csv'.format(comp, product), index=False)
        
        # Plot
        plt.figure(figsize=(10, 6))
        sns.barplot(data=importance_df, x='Mean Importance', y='Feature', palette='viridis')
        plt.errorbar(importance_df['Mean Importance'], importance_df['Feature'],
                    xerr=[importance_df['Mean Importance'] - importance_df['Lower CI'],
                        importance_df['Upper CI'] - importance_df['Mean Importance']],
                    fmt='none', ecolor='black', capsize=3, label='95% CI')
        plt.xlabel('Bootstrapped Feature Importance')
        plt.tight_layout()
        plt.savefig(r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\{1}_feature_importance_bootstrap_v9.png'.format(comp,product), dpi=500)
        plt.show(block=False)
        plt.pause(1)
        plt.close()

# %%

    # OLD CODE
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
                        #     if cygnss_time[end_ind] - cygnss_time[start_ind] < tc_durations:
                        #         end_ind += 1
                        #     elif  cygnss_time[end_ind] - cygnss_time[start_ind] > tc_durations:
                        #         start_ind += 1
                        #     elif cygnss_time[end_ind] - cygnss_time[start_ind] == tc_durations:
                        #         assert cygnss_time[start_ind]
                        #         change_vmax = cyg_vmax[end_ind] - cyg_vmax[start_ind]
                        #         # change_vmax_besttrack = vmax_besttrack[end_ind] - vmax_besttrack[start_ind]
                        #         if change_vmax > tc_intensities: # and change_vmax_besttrack > tc_intensities:
                        #             print(int(change_vmax),tc_intensities,tc_durations,tc_storm_id, tc_event,cyg_storm_names[valid_file_ind], vmax_lat[0],vmax_lon[0],cyg_AOI_time_start,valid_file_ind)
                        #         start_ind += 1
                        #         end_ind += 1
                        #     else:
                        #         end_ind += 1

                            # new_row = ['over',tc_intensities,tc_durations,tc_storm_id, tc_event,int(change_vmax),int(time_diff),'NA']
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
          #         # if  time_diff < tc_durations*60/2 - 60: # custon threshold of minimal separation
          #         #     break
                  
          #         # Check wind speed difference
          #         wind_change = sorted_wind_speeds[j] - sorted_wind_speeds[i]
          #         if wind_change > 12:
          #         # if wind_change > tc_intensities:
          #             significant_change_mask[i] = True
          #             significant_change_mask[j] = True
          
          # # Select measurements with significant wind changes
          # changing_times = sorted_times[significant_change_mask]
          # changing_wind_speeds = sorted_wind_speeds[significant_change_mask]
          
          # if True in significant_change_mask:

# %%
