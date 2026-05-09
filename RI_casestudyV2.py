# -*- coding: utf-8 -*-
"""
RI Case Study Analysis V2
Performs detailed case study analysis on selected TC events for YSLF, Merged, and NOAA compared to IBTrACS reference.

@author: ashle
"""

import netCDF4 as nc
import numpy as np
from datetime import datetime, timedelta
import pandas as pd
import glob
import os
import re
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
from RI_casestudy_insitu_finder import haversine_vectorized

# Configuration
comp = 'ashle'
drive = {'S3987712': 'D:', 'ashle': 'E:'}

# Case study events (list of [storm_name, RI_index])
case_study_events = [
    ['Helene_2024', 665],
    ['Hinnamnor_2022', 428],
    ['Kyarr_2019', 210],
    ['Bualoi_2019', 290],
]

# Processing parameters
interp_interval = 20  # minutes to interpolate best-track
MAX_DISTANCE = 80  # Maximum distance in kilometers
noaa_time_threshold = 3  # hours - maximum time difference allowed
buffer = MAX_DISTANCE / 100  # degrees buffer to search in storm_centric

# Data source folders
noaa_cyg_folder = r'E:\Phd_data\cyg_noaa_1.2'
l2_3p2_cyg_folder = r'E:\Phd_data\cyg_l2_3.2'
merged_cyg_folder = r'E:\Phd_data\cyg_storm_centric'

# Variables to extract
variables = ['incidence_angle', 'rx_gain', 'snr', 'range_corr_gain', 'sample_flags',
             'num_ddms_utilized', 'nbrcs_mean_corrected', 'nbrcs_mean',
             'wind_speed_uncertainty', 'wind_speed']

# Variable aliases for different data sources
source_var_aliases = {
    'l2_3.2': {
        'wind_speed': ['preliminary_yslf_wind_speed'],
        'wind_speed_uncertainty': ['preliminary_yslf_wind_speed_uncertainty'],
        'sample_flags': ['yslf_sample_flags'],
    }
}

var_aliases = {
    'rx_gain': ['port_gain_setting', 'starboard_gain_setting'],
    'snr': ['bit_ratio_lo_hi_port', 'bit_ratio_lo_hi_starboard', 'bit_ratio_lo_hi_zenith'],
    'nbrcs_mean_corrected': ['ddm_nbrcs'],
}

avg_alias_fields = {
    'rx_gain': ['port_gain_setting', 'starboard_gain_setting'],
    'snr': ['bit_ratio_lo_hi_port', 'bit_ratio_lo_hi_starboard', 'bit_ratio_lo_hi_zenith'],
}

# Load IBTrACS RI events
RI_file = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Data', 'IBTrACS', 'RI_events_v5.nc')
RI_nc = nc.Dataset(RI_file)
time_int = RI_nc.variables['times'][:]
time = np.datetime64('1970-01-01T00:00:00') + time_int.astype('timedelta64[s]')
storm_id = RI_nc.variables['storm'][:]
vmax = RI_nc.variables['vmax'][:] * 0.514  # convert to m/s
latitude = RI_nc.variables['latitude'][:]
longitude = RI_nc.variables['longitude'][:]
ri_rd = RI_nc.variables['ri_rd'][:]
durations = RI_nc.variables['durations'][:]
intensities = RI_nc.variables['intensities'][:] * 0.514
index_after = RI_nc.variables['index_after'][:]
index_before = RI_nc.variables['index_before'][:]
RI_nc.close()

# Helper functions
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

def valid_RI_meas(vmax_time, max_wind, tc_durations, ri_start_time, noaa_time_threshold):
    """
    Find valid RI measurements following RI_analysis.py methodology.

    Returns:
        tuple: (has_valid_data, wind_change, time_diff_hours)
    """
    # No measurements for this source
    if len(vmax_time) == 0 or len(max_wind) == 0:
        return (False, 0, 0)

    # Find initial maximum: within 3 hours of earliest measurement or RI start, whichever is earlier
    earlist_track_time = np.nanmin(vmax_time)
    valid_initial_mask = ((vmax_time <= earlist_track_time + np.timedelta64(noaa_time_threshold, 'h')) &
                          (vmax_time <= ri_start_time + np.timedelta64(noaa_time_threshold, 'h')))

    # Fallback when RI timing filter is too strict
    if not np.any(valid_initial_mask):
        valid_initial_mask = (vmax_time <= earlist_track_time + np.timedelta64(noaa_time_threshold, 'h'))
    if not np.any(valid_initial_mask):
        return (False, 0, 0)

    valid_initial_indices = np.where(valid_initial_mask)[0]
    initial_max_ind = valid_initial_indices[np.nanargmax(max_wind[valid_initial_indices])]

    vmax_time_arr = np.array(vmax_time, dtype='datetime64[s]')
    # Calculate time differences in hours from initial to all other maxima
    time_diffs = np.round((vmax_time_arr - vmax_time_arr[initial_max_ind]) / np.timedelta64(1, 'h'))
    valid_mask = (time_diffs >= tc_durations / 2) & (time_diffs <= tc_durations + 2)

    if np.any(valid_mask):
        # Of valid indices, pick the one with max wind
        valid_indices = np.where(valid_mask)[0]
        final_max_ind = valid_indices[np.argmax(max_wind[valid_indices])]

        wind_change = max_wind[final_max_ind] - max_wind[initial_max_ind]
        time_diff = time_diffs[final_max_ind]
        return (True, wind_change, time_diff)
    return (False, 0, 0)

def filter_files_by_date(file_list, date_list):
    """Pre-filter files to only those overlapping with date range (huge speedup)."""
    filtered_files = []
    filtered_dates = []

    for file_path in file_list:
        file_date = parse_cyg_file_date(file_path)
        if file_date is not None and np.any(np.isin(date_list, np.array([file_date]))):
            filtered_files.append(file_path)
            filtered_dates.append(file_date)

    return np.array(filtered_files), np.array(filtered_dates)

# Get file dates for merged data
merged_cyg_files_list = np.asarray(glob.glob(merged_cyg_folder + r"\*.nc"))
merged_file_dates = []
cyg_files_num = []
cyg_storm_names = []
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
    cyg_nc.close()

# Setup file lookups for NOAA and L2 sources
noaa_like_folders = {
    'l2_3.2': l2_3p2_cyg_folder,
    'noaa_1.2': noaa_cyg_folder
}
noaa_like_files_map = {source_name: np.asarray(glob.glob(source_folder + r"\*.nc"))
                       for source_name, source_folder in noaa_like_folders.items()}
date_to_noaa_idx = {source_name: {} for source_name in noaa_like_files_map}
for source_name, cyg_files in noaa_like_files_map.items():
    for idx, cyg_file in enumerate(cyg_files):
        noaa_date = parse_cyg_file_date(cyg_file)
        if noaa_date is None:
            continue
        date_to_noaa_idx[source_name][noaa_date] = idx

noaa_cache = {source_name: {} for source_name in noaa_like_files_map}

# Process each case study event
for case_study in case_study_events:
    print(f"\n{'='*80}")
    print(f"Processing: {case_study[0]} (RI index: {case_study[1]})")
    print(f"{'='*80}\n")

    tc_event = case_study[1]

    if ri_rd[tc_event] != 1 or durations[tc_event] >= 25:
        print(f"Skipping event (not RI or duration too long)")
        continue

    # Extract TC event data
    nans = ~np.isnan(latitude[tc_event]).data
    tc_time = time[tc_event, nans]
    tc_storm_id = storm_id[tc_event]
    tc_vmax = vmax[tc_event, nans]
    tc_durations = durations[tc_event]
    tc_intensities = intensities[tc_event]
    tc_latitude = latitude[tc_event, nans]
    tc_longitude = longitude[tc_event, nans]
    ri_start_time = np.datetime64(tc_time[index_before[tc_event]])
    ri_end_time = np.datetime64(tc_time[-index_after[tc_event]])
    ri_vmax_end = tc_vmax[-index_after[tc_event]]
    ri_vmax_start = tc_vmax[index_before[tc_event]]
    ri_vmax_change = ri_vmax_end - ri_vmax_start

    # Interpolate TC track
    interval = np.timedelta64(interp_interval, 'm')
    start = tc_time[0]
    end = tc_time[-1]
    tc_time_interpolated = np.arange(start, end + interval, interval)

    # Handle masked vmax values
    filled = tc_vmax.copy()
    for i in range(1, len(tc_vmax) - 1):
        if tc_vmax.mask[i] and not tc_vmax.mask[i - 1] and not tc_vmax.mask[i + 1]:
            filled[i] = (tc_vmax[i - 1] + tc_vmax[i + 1]) / 2
            filled.mask[i] = False
    tc_vmax = filled

    # Vectorized interpolation
    tc_times_s = tc_time.astype('datetime64[s]').astype('int64')
    interp_s = tc_time_interpolated.astype('datetime64[s]').astype('int64')

    tc_vmax_vals = np.ma.filled(tc_vmax, np.nan).astype(float)
    tc_lat_vals = np.asarray(tc_latitude).astype(float)
    tc_lon_vals = np.asarray(tc_longitude).astype(float)

    tc_vmax_interpolated = np.interp(interp_s, tc_times_s, tc_vmax_vals)
    tc_lat_interpolated = np.interp(interp_s, tc_times_s, tc_lat_vals)
    tc_lon_interpolated = np.interp(interp_s, tc_times_s, tc_lon_vals)

    # Define AOI
    AOI_lat_min, AOI_lat_max = np.nanmin(tc_latitude) - buffer, np.nanmax(tc_latitude) + buffer
    AOI_lon_min, AOI_lon_max = np.nanmin(tc_longitude) - buffer, np.nanmax(tc_longitude) + buffer
    cyg_AOI_time_start = np.nanmin(tc_time)
    cyg_AOI_time_end = np.nanmax(tc_time)
    cyg_AOI_time_start_dt = cyg_AOI_time_start.astype('datetime64[s]').astype('O')
    cyg_AOI_time_end_dt = cyg_AOI_time_end.astype('datetime64[s]').astype('O')
    date_list = np.arange(cyg_AOI_time_start.astype('datetime64[D]'),
                          cyg_AOI_time_end.astype('datetime64[D]') + np.timedelta64(1, 'D'),
                          np.timedelta64(1, 'D'))

    # Storage for results from each source
    results_by_source = {
        'noaa_1.2': {'lats': [], 'lons': [], 'winds': [], 'times': [], 'tracks': [], 'data': {v: [] for v in variables}},
        'l2_3.2': {'lats': [], 'lons': [], 'winds': [], 'times': [], 'tracks': [], 'data': {v: [] for v in variables}},
        'merged': {'epochs': [], 'vmax': [], 'times': []}
    }

    # Process NOAA and L2 sources
    print("Processing NOAA 1.2 and L2 3.2 data...")
    for source_name, source_files in noaa_like_files_map.items():
        source_date_to_idx = date_to_noaa_idx.get(source_name, {})
        source_cache = noaa_cache.get(source_name, {})

        for date in date_list:
            try:
                cyg_file_path = str(source_files[source_date_to_idx[date]])
            except KeyError:
                continue

            # Load or retrieve from cache
            if cyg_file_path in source_cache:
                cached = source_cache[cyg_file_path]
                cygnss_lons = cached['lons']
                cygnss_lats = cached['lats']
                cygnss_wind = cached['wind']
                cygnss_tracks = cached['tracks']
                cygnss_time = cached['times_full']
            else:
                cyg_nc = nc.Dataset(cyg_file_path)
                longitudes = cyg_nc.variables['lon'][:]
                cygnss_lats = cyg_nc.variables['lat'][:]

                # Select appropriate wind speed variable
                source_aliases = source_var_aliases.get(source_name, {})
                wind_candidates = source_aliases.get('wind_speed', []) + ['wind_speed']
                wind_var_name = None
                for candidate in wind_candidates:
                    if candidate in cyg_nc.variables:
                        wind_var_name = candidate
                        break
                if wind_var_name is None:
                    cyg_nc.close()
                    continue

                cygnss_wind = cyg_nc.variables[wind_var_name][:]

                # Extract track IDs
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
                date_format = "%Y-%m-%d %H:%M:%S"
                cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                times_full = np.datetime64(cyg_start_date) + cygnss_time.astype('timedelta64[s]')

                # Cache variables
                vars_cache = {}
                for varname in variables:
                    try:
                        if varname in cyg_nc.variables:
                            vars_cache[varname] = cyg_nc.variables[varname][:]
                        else:
                            source_aliases_var = source_var_aliases.get(source_name, {})
                            source_candidates = source_aliases_var.get(varname, [])
                            source_found = None
                            for alias_name in source_candidates:
                                if alias_name in cyg_nc.variables:
                                    source_found = alias_name
                                    break
                            if source_found is not None:
                                vars_cache[varname] = cyg_nc.variables[source_found][:]
                            elif varname in avg_alias_fields:
                                alias_arrays = []
                                for alias_name in avg_alias_fields[varname]:
                                    if alias_name in cyg_nc.variables:
                                        alias_vals = np.asarray(np.ma.filled(cyg_nc.variables[alias_name][:], np.nan), dtype=float)
                                        alias_arrays.append(alias_vals)
                                if len(alias_arrays) > 0:
                                    vars_cache[varname] = np.nanmean(np.stack(alias_arrays, axis=0), axis=0)
                            else:
                                aliases = var_aliases.get(varname, [])
                                alias_found = None
                                for alias_name in aliases:
                                    if alias_name in cyg_nc.variables:
                                        alias_found = alias_name
                                        break
                                if alias_found is not None:
                                    vars_cache[varname] = cyg_nc.variables[alias_found][:]
                                else:
                                    vars_cache[varname] = np.array([])
                    except Exception:
                        vars_cache[varname] = np.array([])

                cyg_nc.close()
                cached = {
                    'lons': cygnss_lons,
                    'lats': cygnss_lats,
                    'wind': cygnss_wind,
                    'tracks': cygnss_tracks,
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

            # Box selection
            lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
            lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
            box_indices = np.where(lats_boolean & lons_boolean)[0]

            if len(box_indices) > 0:
                meas_lats = cygnss_lats[box_indices]
                meas_lons = cygnss_lons[box_indices]
                meas_winds = cygnss_wind[box_indices]
                meas_times = cygnss_time[box_indices]

                # Quality filtering for L2 3.2
                quality_mask = np.ones(len(box_indices), dtype=bool)
                if source_name == 'l2_3.2':
                    cached_vars = cached.get('vars', {})
                    meas_flags = _extract_box_float_var(cached_vars, 'sample_flags', box_indices)
                    if meas_flags.size > 0:
                        quality_mask &= np.isfinite(meas_flags) & (meas_flags == 0)
                    if not np.any(quality_mask):
                        continue
                    box_indices = box_indices[quality_mask]
                    meas_lats = meas_lats[quality_mask]
                    meas_lons = meas_lons[quality_mask]
                    meas_winds = meas_winds[quality_mask]
                    meas_times = meas_times[quality_mask]

                # Spatiotemporal matching
                meas_coords = np.column_stack((meas_lats, meas_lons))
                tc_coords = np.column_stack((tc_lat_interpolated, tc_lon_interpolated))

                distances = haversine_vectorized(meas_coords, tc_coords)
                time_diffs = np.abs(meas_times[:, None] - tc_time_interpolated[None, :])
                valid_mask = np.any((distances <= MAX_DISTANCE) & (time_diffs <= np.timedelta64(noaa_time_threshold, 'h')), axis=1)
                selected_indices = np.where(valid_mask)[0]

                if len(selected_indices) > 0:
                    cyg_winds_vals = np.ma.filled(meas_winds[selected_indices], np.nan).astype(float)
                    cyg_tracks_vals = cygnss_tracks[box_indices][selected_indices]

                    results_by_source[source_name]['lats'].extend(np.asarray(meas_lats[selected_indices]).tolist())
                    results_by_source[source_name]['lons'].extend(np.asarray(meas_lons[selected_indices]).tolist())
                    results_by_source[source_name]['winds'].extend(cyg_winds_vals.tolist())
                    results_by_source[source_name]['times'].extend(np.asarray(meas_times[selected_indices]).tolist())
                    results_by_source[source_name]['tracks'].extend(np.asarray(cyg_tracks_vals).tolist())

                    for varname in variables:
                        var_full = cached.get('vars', {}).get(varname, None)
                        if var_full is None or getattr(var_full, 'size', 0) == 0:
                            var = np.full(len(selected_indices), np.nan)
                        else:
                            var = np.atleast_1d(np.asarray(var_full)[box_indices[selected_indices]])
                        if len(var) == len(selected_indices):
                            results_by_source[source_name]['data'][varname].extend(var.tolist())

    # Process Merged data
    print("Processing Merged storm-centric data...")
    if np.any(np.isin(date_list, np.array(merged_file_dates))):
        valid_files = np.where(date_list[0] == np.array(merged_file_dates))[0]
        unique_storm_files = []
        if len(valid_files) == 1:
            unique_storm_files.append(merged_cyg_files_list[cyg_files_num[valid_files[0]]])
        else:
            for i in range(len(valid_files) - 1):
                if merged_cyg_files_list[cyg_files_num[valid_files[i]]] != merged_cyg_files_list[cyg_files_num[valid_files[i+1]]]:
                    unique_storm_files.append(merged_cyg_files_list[cyg_files_num[valid_files[i]]])
        unique_storm_files = list(dict.fromkeys(unique_storm_files))

        for cyg_file_path in unique_storm_files:
            cyg_nc = nc.Dataset(cyg_file_path)
            cyg_storm_lat = cyg_nc.variables['best_track_storm_center_lat'][:]
            cyg_storm_lon = ((cyg_nc.variables['best_track_storm_center_lon'][:] + 180) % 360) - 180

            if not np.any((cyg_storm_lon > AOI_lon_min) & (cyg_storm_lon < AOI_lon_max) &
                         (cyg_storm_lat > AOI_lat_min) & (cyg_storm_lat < AOI_lat_max)):
                cyg_nc.close()
                continue

            epochs = np.where((cyg_storm_lon > AOI_lon_min) & (cyg_storm_lon < AOI_lon_max) &
                            (cyg_storm_lat > AOI_lat_min) & (cyg_storm_lat < AOI_lat_max))[0]

            if len(epochs) > 1:
                cygnss_wind = cyg_nc.variables['wind_speed'][epochs]
                cygnss_time = cyg_nc.variables['time'][epochs]
                cygnss_time_units = cyg_nc.variables['time'].units[12:]
                cyg_datetimes = np.datetime64(cygnss_time_units) + cygnss_time.astype('timedelta64[h]')
                cyg_vmax = np.full(len(epochs), np.nan, dtype=float)

                # Extract wind speed at vmax location
                cygnss_lats = np.asarray(cyg_nc.variables['lat'][:], dtype=float)
                cygnss_lons = np.asarray(((cyg_nc.variables['lon'][:] + 180) % 360) - 180, dtype=float)
                vmax_lat = np.asarray(cyg_nc.variables['cygnss_vmax_lat'][epochs], dtype=float)
                vmax_lon = np.asarray(((cyg_nc.variables['cygnss_vmax_lon'][epochs] + 180) % 360) - 180, dtype=float)

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
                    wind_val = cygnss_wind[j, lat_idx, lon_idx]
                    cyg_vmax[j] = np.nan if np.ma.is_masked(wind_val) else float(wind_val)

                # Find RI initial and final epochs
                time_diffs = np.abs(cyg_datetimes - ri_start_time)
                time_diffs_hours = time_diffs / np.timedelta64(1, 'h')
                valid_time_mask = (time_diffs_hours <= noaa_time_threshold) & (cyg_vmax > 0)

                if np.any(valid_time_mask):
                    start_epoch = int(np.where(valid_time_mask)[0][np.argmin(time_diffs_hours[valid_time_mask])])

                    time_deltas = np.ma.asarray(cyg_datetimes - cyg_datetimes[start_epoch])
                    time_diffs_hr = np.ma.filled(time_deltas / np.timedelta64(1, 'h'), np.nan).astype(float)
                    valid_mask = (time_diffs_hr >= tc_durations / 2) & (time_diffs_hr <= tc_durations + 2) & (cyg_vmax > 0)
                    valid_indices = np.where(valid_mask)[0]

                    if len(valid_indices) > 0:
                        results_by_source['merged']['epochs'].append([start_epoch, cyg_vmax[start_epoch], cyg_datetimes[start_epoch]])
                        results_by_source['merged']['epochs'].append([valid_indices[np.argmax(cyg_vmax[valid_indices])],
                                                                      cyg_vmax[valid_indices[np.argmax(cyg_vmax[valid_indices])]],
                                                                      cyg_datetimes[valid_indices[np.argmax(cyg_vmax[valid_indices])]]])

            cyg_nc.close()

    print("\nCreating comparison summary...")

    # Summary statistics
    summary_data = {
        'source': [],
        'num_measurements': [],
        'num_tracks': [],
        'mean_wind': [],
        'max_wind': [],
        'min_wind': [],
        'initial_wind': [],
        'final_wind': [],
        'wind_change_estimate': [],
        'ri_valid': []
    }

    # Save measurement metadata for each source to CSV
    print(f"Saving measurement metadata for {tc_storm_id}...")
    for source_name in ['noaa_1.2', 'l2_3.2', 'merged']:
        if source_name == 'merged':
            if len(results_by_source['merged']['epochs']) > 0:
                merged_data = {
                    'epoch_index': [e[0] for e in results_by_source['merged']['epochs']],
                    'wind_speed': [e[1] for e in results_by_source['merged']['epochs']],
                    'time': [e[2] for e in results_by_source['merged']['epochs']],
                }
                merged_df = pd.DataFrame(merged_data)
                output_dir = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Plots', 'case_studies')
                os.makedirs(output_dir, exist_ok=True)
                merged_csv = os.path.join(output_dir, f'{tc_event}_Merged_measurements.csv')
                merged_df.to_csv(merged_csv, index=False)
                print(f"  Saved: {os.path.basename(merged_csv)}")
        else:
            if len(results_by_source[source_name]['winds']) > 0:
                meas_data = {
                    'latitude': results_by_source[source_name]['lats'],
                    'longitude': results_by_source[source_name]['lons'],
                    'wind_speed': results_by_source[source_name]['winds'],
                    'time': results_by_source[source_name]['times'],
                }

                # Add all extracted variables
                for varname in variables:
                    if varname in results_by_source[source_name]['data']:
                        meas_data[varname] = results_by_source[source_name]['data'][varname]

                meas_df = pd.DataFrame(meas_data)

                # Determine product name for filename
                product_name = 'NOAA' if source_name == 'noaa_1.2' else 'YSLF'
                output_dir = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Plots', 'case_studies')
                os.makedirs(output_dir, exist_ok=True)
                meas_csv = os.path.join(output_dir, f'{tc_event}_{product_name}_measurements.csv')
                meas_df.to_csv(meas_csv, index=False)
                print(f"  Saved: {os.path.basename(meas_csv)}")

    # Process NOAA and L2 sources with track-based analysis
    for source_name in ['noaa_1.2', 'l2_3.2']:
        if len(results_by_source[source_name]['winds']) > 0:
            winds = np.array(results_by_source[source_name]['winds'])
            times = np.array(results_by_source[source_name]['times'], dtype='datetime64[s]')
            tracks = np.array(results_by_source[source_name]['tracks'], dtype=object)

            # Group measurements by track
            track_keys = [_track_key(t) for t in tracks]
            unique_tracks = list(dict.fromkeys(track_keys))

            # Find peak wind for each track
            vmax_time = np.zeros(len(unique_tracks), dtype=object)
            max_wind = np.full(len(unique_tracks), np.nan, dtype=float)
            valid_track_peak = np.zeros(len(unique_tracks), dtype=bool)

            for i, track in enumerate(unique_tracks):
                track_mask = np.array([k == track for k in track_keys], dtype=bool)
                track_times = times[track_mask]
                track_winds = winds[track_mask]

                if track_times.size > 0:
                    track_winds_vals = np.asarray(np.ma.filled(track_winds, np.nan), dtype=float)
                    if np.any(np.isfinite(track_winds_vals)):
                        vmax_index = int(np.nanargmax(track_winds_vals))
                        vmax_time[i] = track_times[vmax_index]
                        max_wind[i] = track_winds_vals[vmax_index]
                        valid_track_peak[i] = True

            # Apply valid RI measurement logic
            valid_peaks = vmax_time[valid_track_peak]
            valid_winds = max_wind[valid_track_peak]

            if len(valid_peaks) > 0:
                ri_result = valid_RI_meas(valid_peaks, valid_winds, tc_durations, ri_start_time, noaa_time_threshold)

                summary_data['source'].append(source_name)
                summary_data['num_measurements'].append(len(winds))
                summary_data['num_tracks'].append(len(unique_tracks))
                summary_data['mean_wind'].append(np.nanmean(winds))
                summary_data['max_wind'].append(np.nanmax(winds))
                summary_data['min_wind'].append(np.nanmin(winds))
                summary_data['ri_valid'].append(ri_result[0])

                if ri_result[0]:
                    # Find initial and final peaks
                    earlist_track_time = np.nanmin(valid_peaks)
                    valid_initial_mask = ((valid_peaks <= earlist_track_time + np.timedelta64(noaa_time_threshold, 'h')) &
                                         (valid_peaks <= ri_start_time + np.timedelta64(noaa_time_threshold, 'h')))
                    if not np.any(valid_initial_mask):
                        valid_initial_mask = (valid_peaks <= earlist_track_time + np.timedelta64(noaa_time_threshold, 'h'))

                    valid_initial_indices = np.where(valid_initial_mask)[0]
                    initial_max_ind = valid_initial_indices[np.nanargmax(valid_winds[valid_initial_indices])]

                    vmax_time_arr = np.array(valid_peaks, dtype='datetime64[s]')
                    time_diffs = np.round((vmax_time_arr - vmax_time_arr[initial_max_ind]) / np.timedelta64(1, 'h'))
                    valid_mask = (time_diffs >= tc_durations / 2) & (time_diffs <= tc_durations + 2)
                    valid_indices = np.where(valid_mask)[0]
                    final_max_ind = valid_indices[np.argmax(valid_winds[valid_indices])]

                    initial_wind = valid_winds[initial_max_ind]
                    final_wind = valid_winds[final_max_ind]

                    summary_data['initial_wind'].append(initial_wind)
                    summary_data['final_wind'].append(final_wind)
                    summary_data['wind_change_estimate'].append(final_wind - initial_wind)
                else:
                    summary_data['initial_wind'].append(np.nan)
                    summary_data['final_wind'].append(np.nan)
                    summary_data['wind_change_estimate'].append(np.nan)

    # Process Merged source
    if len(results_by_source['merged']['epochs']) > 0:
        vmax_values = [e[1] for e in results_by_source['merged']['epochs']]
        vmax_times = np.array([e[2] for e in results_by_source['merged']['epochs']], dtype='datetime64[s]')

        summary_data['source'].append('merged')
        summary_data['num_measurements'].append(len(results_by_source['merged']['epochs']))
        summary_data['num_tracks'].append(1)  # Storm-centric, single "track"
        summary_data['mean_wind'].append(np.nanmean(vmax_values))
        summary_data['max_wind'].append(np.nanmax(vmax_values))
        summary_data['min_wind'].append(np.nanmin(vmax_values))

        if len(vmax_values) >= 2:
            initial_wind = vmax_values[0]
            final_wind = vmax_values[-1]
            summary_data['initial_wind'].append(initial_wind)
            summary_data['final_wind'].append(final_wind)
            summary_data['wind_change_estimate'].append(final_wind - initial_wind)
            summary_data['ri_valid'].append(True)
        else:
            summary_data['initial_wind'].append(np.nan)
            summary_data['final_wind'].append(np.nan)
            summary_data['wind_change_estimate'].append(np.nan)
            summary_data['ri_valid'].append(False)

    summary_df = pd.DataFrame(summary_data)

    # Add IBTrACS reference
    summary_df.loc[len(summary_df)] = {
        'source': 'IBTrACS',
        'num_measurements': 1,
        'num_tracks': 1,
        'mean_wind': (ri_vmax_start + ri_vmax_end) / 2,
        'max_wind': ri_vmax_end,
        'min_wind': ri_vmax_start,
        'initial_wind': ri_vmax_start,
        'final_wind': ri_vmax_end,
        'wind_change_estimate': ri_vmax_change,
        'ri_valid': True
    }

    print("\nCase Study Summary:")
    print(f"Storm ID: {tc_storm_id}")
    print(f"Event Duration: {tc_durations} hours")
    print(f"RI Period: {ri_start_time} to {ri_end_time}")
    print(f"IBTrACS Vmax Change: {ri_vmax_change:.2f} m/s ({ri_vmax_start:.2f} → {ri_vmax_end:.2f} m/s)")
    print("\nData Source Comparison:")
    print(summary_df.to_string(index=False))

    # Create comparison plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    sources = summary_df['source'].values
    colors = {'noaa_1.2': 'blue', 'l2_3.2': 'orange', 'merged': 'green', 'IBTrACS': 'red'}

    # Wind change comparison
    ax = axes[0, 0]
    for i, source in enumerate(sources):
        wind_change = summary_df.loc[summary_df['source'] == source, 'wind_change_estimate'].values[0]
        is_valid = summary_df.loc[summary_df['source'] == source, 'ri_valid'].values[0]
        marker = 'o' if is_valid else 'x'
        ax.scatter(i, wind_change, s=200, color=colors.get(source, 'gray'), alpha=0.7, marker=marker, edgecolors='black', linewidth=2)
    ax.axhline(y=ri_vmax_change, color='red', linestyle='--', linewidth=2, alpha=0.5, label='IBTrACS')
    ax.set_ylabel('Vmax Change (m/s)', fontsize=11, fontweight='bold')
    ax.set_title('Intensity Change Comparison', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(sources)))
    ax.set_xticklabels(sources, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)
    ax.legend()

    # Initial vs Final wind
    ax = axes[0, 1]
    for i, source in enumerate(sources):
        initial = summary_df.loc[summary_df['source'] == source, 'initial_wind'].values[0]
        final = summary_df.loc[summary_df['source'] == source, 'final_wind'].values[0]
        is_valid = summary_df.loc[summary_df['source'] == source, 'ri_valid'].values[0]
        if np.isfinite(initial) and np.isfinite(final):
            ax.arrow(i, initial, 0, final - initial, head_width=0.15, head_length=0.5,
                    fc=colors.get(source, 'gray'), ec=colors.get(source, 'gray'), alpha=0.7, linewidth=2)
            ax.scatter([i, i], [initial, final], s=100, color=colors.get(source, 'gray'),
                      alpha=0.7, edgecolors='black', linewidth=1.5, zorder=5)
    ax.set_ylabel('Wind Speed (m/s)', fontsize=11, fontweight='bold')
    ax.set_title('Initial → Final Peak Winds', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(sources)))
    ax.set_xticklabels(sources, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)

    # Measurement count
    ax = axes[1, 0]
    for i, source in enumerate(sources):
        num_meas = summary_df.loc[summary_df['source'] == source, 'num_measurements'].values[0]
        ax.bar(i, num_meas, color=colors.get(source, 'gray'), alpha=0.7, edgecolor='black', linewidth=1.5)
    ax.set_ylabel('Number of Measurements', fontsize=11, fontweight='bold')
    ax.set_title('Data Point Count', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(sources)))
    ax.set_xticklabels(sources, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)

    # Track count
    ax = axes[1, 1]
    for i, source in enumerate(sources):
        num_tracks = summary_df.loc[summary_df['source'] == source, 'num_tracks'].values[0]
        ax.bar(i, num_tracks, color=colors.get(source, 'gray'), alpha=0.7, edgecolor='black', linewidth=1.5)
    ax.set_ylabel('Number of Satellite Tracks', fontsize=11, fontweight='bold')
    ax.set_title('Satellite Track Count', fontsize=12, fontweight='bold')
    ax.set_xticks(range(len(sources)))
    ax.set_xticklabels(sources, rotation=45, ha='right')
    ax.grid(axis='y', alpha=0.3)

    plt.suptitle(f'{tc_storm_id} - RI Case Study Comparison\n{ri_start_time} to {ri_end_time}',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Save figure
    output_dir = os.path.join('C:\\Users', comp, 'OneDrive - RMIT University', 'PHD', 'Plots', 'case_studies')
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f'{tc_storm_id}_casestudy_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved comparison plot to: {output_file}")
    plt.show(block=False)
    plt.pause(1)
    plt.close(fig)

print("\n" + "="*80)
print("Case study analysis complete!")
print("="*80)
