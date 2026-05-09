# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 10:32:03 2025

# CMS api downloader

@author: ashle
"""

from datetime import datetime

from scipy import datasets


def cms_downloader(minimum_longitude,  maximum_longitude,  minimum_latitude,  maximum_latitude,start_datetime, end_datetime,storm_id,overide=False,buffer = 5):
    import copernicusmarine
    import datetime
    import os
    #  # check if files already exist
    # filepath_swh = f"E:\\Phd_data\\copernicus\\{storm_id}_SWH.nc"
    # filepath_sss_sst = f"E:\\Phd_data\\copernicus\\{storm_id}_SSS&SST.nc"
    # if overide == True:
    #     if os.path.exists(filepath_swh):
    #         os.remove(filepath_swh)
    #     if os.path.exists(filepath_sss_sst):
    #         os.remove(filepath_sss_sst)
    # if os.path.exists(filepath_swh) or os.path.exists(filepath_sss_sst):   
    #     print(f"SWH file for {storm_id}  already exists. Skipping download.")
    #     return  # Skip download if file already exists

    # SWH and wave direction 
    copernicusmarine.subset(
    dataset_id="cmems_mod_glo_wav_my_0.2deg_PT3H-i",
    dataset_version="202411",
    variables=["VHM0", "VMDR"],
    minimum_longitude=minimum_longitude-buffer,
    maximum_longitude=maximum_longitude+buffer,
    minimum_latitude=minimum_latitude-buffer,
    maximum_latitude=maximum_latitude+buffer,
    start_datetime=start_datetime,
    end_datetime=end_datetime +datetime.timedelta(hours=48),
    coordinates_selection_method="strict-inside",
    netcdf_compression_level=1,
    disable_progress_bar=True,
    username='asiedlecki',
    password='XsySq5LJ8TmZKs',
    output_filename = f"{storm_id}_SWH.nc",
    output_directory = r'E:\Phd_data\copernicus'
    )

    # SSS and SST
    if end_datetime > datetime.datetime(2022, 6, 1, 23, 59, 59):
        copernicusmarine.subset(
        dataset_id="cmems_mod_glo_phy_anfc_0.083deg_PT1H-m",
        variables=["thetao", "so"],
        minimum_longitude=minimum_longitude-buffer,
        maximum_longitude=maximum_longitude+buffer,
        minimum_latitude=minimum_latitude-buffer,
        maximum_latitude=maximum_latitude+buffer,
        start_datetime=start_datetime,
        end_datetime=end_datetime +datetime.timedelta(hours=48),
        minimum_depth=0.49402499198913574,
        maximum_depth=0.49402499198913574,
        username='asiedlecki',
        password='XsySq5LJ8TmZKs',
        output_filename = f"{storm_id}_SSS&SST.nc",
        output_directory = r'E:\Phd_data\copernicus'
        )
    else:
        copernicusmarine.subset(
        dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
        dataset_version="202311",
        variables=["thetao", "so"],
        minimum_longitude=minimum_longitude-buffer,
        maximum_longitude=maximum_longitude+buffer,
        minimum_latitude=minimum_latitude-buffer,
        maximum_latitude=maximum_latitude+buffer,
        start_datetime=start_datetime,
        end_datetime=end_datetime +datetime.timedelta(hours=48),
        minimum_depth=0.49402499198913574,
        maximum_depth=0.49402499198913574,
        coordinates_selection_method="strict-inside",
        netcdf_compression_level=1,
        disable_progress_bar=True,
        username='asiedlecki',
        password='XsySq5LJ8TmZKs',
        output_filename = f"{storm_id}_SSS&SST.nc",
        output_directory = r'E:\Phd_data\copernicus'
        )


# def cms_vals(meas_lats, meas_lons, meas_times,storm_id):
#     # a function that returns CMS values for each triplet of CYGNSS measurement lat, lon, time
#     import numpy as np
#     import netCDF4 as nc
#     import os
#     import datetime
#     import pandas as pd

#     # Helper to request CMS files covering the measurements (with buffers)
#     def _request_if_needed_box(min_lat, max_lat, min_lon, max_lon, start_dt, end_dt, storm_id, lat_buffer=0, lon_buffer=0, overide_first=True):
#         """Request CMS data for a given box. Handles dateline-crossing by splitting into two requests if needed.
#         min_*/max_* are in degrees (can be -180..180 or 0..360); start_dt/end_dt are datetimes.
#         If overide_first is True the first call will remove existing files (overide=True) and subsequent split calls will not.
#         """
#         # apply buffers (buffers may have already been applied upstream but do here for safety)
#         min_lon = float(min_lon) - lon_buffer
#         max_lon = float(max_lon) + lon_buffer
#         min_lat = float(min_lat) - lat_buffer
#         max_lat = float(max_lat) + lat_buffer

#         # convert to 0..360 for wrap-aware checks
#         min_lon360 = (min_lon + 360) % 360
#         max_lon360 = (max_lon + 360) % 360

#         # If the interval in 0..360 appears to wrap (max < min), split into two requests
#         if max_lon360 < min_lon360:
#             # Convert segments back to -180..180 for cms_downloader
#             seg1_min = ((min_lon360) + 180) % 360 - 180
#             seg1_max = 180.0
#             seg2_min = -180.0
#             seg2_max = ((max_lon360) + 180) % 360 - 180
#             # First call: overide (remove previous partial files)
#             cms_downloader(seg1_min, seg1_max, min_lat, max_lat, start_dt, end_dt, storm_id, overide=overide_first)
#             # Second call: do not overide to avoid deleting first segment
#             cms_downloader(seg2_min, seg2_max, min_lat, max_lat, start_dt, end_dt, storm_id, overide=False)
#         else:
#             # Normal contiguous box
#             # convert back to -180..180 for cms_downloader
#             box_min_lon = ((min_lon360) + 180) % 360 - 180
#             box_max_lon = ((max_lon360) + 180) % 360 - 180
#             cms_downloader(box_min_lon, box_max_lon, min_lat, max_lat, start_dt, end_dt, storm_id, overide=overide_first)

#     def _request_if_needed_arrays(lat_arr, lon_arr, time_arr, storm_id, lat_buffer=10, lon_buffer=10, hours_before=12, hours_after=48):
#         """Compute bounding box from arrays and request CMS files covering that (with buffers and time padding)."""
#         start_dt = pd.to_datetime(np.min(time_arr)).to_pydatetime() - datetime.timedelta(hours=hours_before)
#         end_dt = pd.to_datetime(np.max(time_arr)).to_pydatetime() + datetime.timedelta(hours=hours_after)
#         min_lon = float(np.min(lon_arr))
#         max_lon = float(np.max(lon_arr))
#         min_lat = float(np.min(lat_arr))
#         max_lat = float(np.max(lat_arr))
#         # Request the box (dateline-safe)
#         _request_if_needed_box(min_lat, max_lat, min_lon, max_lon, start_dt, end_dt, storm_id, lat_buffer=lat_buffer, lon_buffer=lon_buffer, overide_first=True)

#     def _closest_lon_idx(lon_array, lon):
#         # Compute wrap-aware distance by converting to 0..360 and using circular distance
#         lon_arr_360 = (np.array(lon_array, dtype=float) + 360) % 360
#         lon_360 = (float(lon) + 360) % 360
#         diffs = np.minimum(np.abs(lon_arr_360 - lon_360), 360 - np.abs(lon_arr_360 - lon_360))
#         return np.argmin(diffs)

#     def _min_lon_dist(lon_array, lon):
#         lon_arr_360 = (np.array(lon_array, dtype=float) + 360) % 360
#         lon_360 = (float(lon) + 360) % 360
#         d = np.minimum(np.abs(lon_arr_360 - lon_360), 360 - np.abs(lon_arr_360 - lon_360))
#         return np.min(d)

#     # Open the CMS files
#     cms_dataset = nc.Dataset(r'E:\Phd_data\copernicus\{}_SSS&SST.nc'.format(storm_id))
#     cms_swh_dataset = nc.Dataset(r'E:\Phd_data\copernicus\{}_SWH.nc'.format(storm_id))
#     cms1_lats = cms_dataset.variables['latitude'][:]
#     cms1_lons = cms_dataset.variables['longitude'][:]
#     tvar1 = cms_dataset.variables['time']
#     cms1_time = nc.num2date(tvar1[:], tvar1.units, getattr(tvar1, 'calendar', 'standard'))
#     cms1_time = np.array(cms1_time, dtype='datetime64[s]')
#     cms_sst = cms_dataset.variables['thetao'][:]  # sea surface temperature
#     cms_sss = cms_dataset.variables['so'][:]      # sea surface salinity

#     cms2_lats = cms_swh_dataset.variables['latitude'][:]
#     cms2_lons = cms_swh_dataset.variables['longitude'][:]
#     tvar2 = cms_swh_dataset.variables['time']
#     cms2_time = nc.num2date(tvar2[:], tvar2.units, getattr(tvar2, 'calendar', 'standard'))
#     cms2_time = np.array(cms2_time, dtype='datetime64[s]')
#     cms_swh = cms_swh_dataset.variables['VHM0'][:]  # significant wave height
#     cms_dir = cms_swh_dataset.variables['VMDR'][:]  # mean wave direction

#     # Normalize longitudes to -180..180 to handle dateline conventions (CMS may use 0..360)
#     def _normalize_to_pm180(lon_arr):
#         lon_arr = np.array(lon_arr, dtype=float)
#         # If longitudes appear in 0..360 range (max > 180), convert to -180..180
#         if np.nanmax(lon_arr) > 180:
#             lon_arr = ((lon_arr + 180) % 360) - 180
#         return lon_arr

#     cms1_lons = _normalize_to_pm180(cms1_lons)
#     cms2_lons = _normalize_to_pm180(cms2_lons)

#     # Ensure measurement longitudes are in same convention as CMS lon arrays
#     meas_lons = np.array(meas_lons, dtype=float)
#     # If cms lons are in -180..180 but meas_lons might not be, normalize meas_lons too
#     if np.nanmax(cms1_lons) <= 180:
#         meas_lons = ((meas_lons + 180) % 360) - 180

#     # Quick diagnostic when coverage fails (commented out by default)
#     # print(f"CMS1 lon range: {np.nanmin(cms1_lons):.2f} to {np.nanmax(cms1_lons):.2f}; meas lon range: {np.nanmin(meas_lons):.2f} to {np.nanmax(meas_lons):.2f}")

#     # Check coverage and automatically attempt a larger re-download if necessary
#     try:
#         # run the existing per-point checks once to surface any out-of-coverage points
#         for i, (lat, lon, time) in enumerate(zip(meas_lats, meas_lons, meas_times)):
#             if min(abs(cms1_lats - lat)) > 0.5 or _min_lon_dist(cms1_lons, lon) > 0.5 or min(abs(cms1_time - time)) > np.timedelta64(24,'h'):
#                 raise ValueError("cms1_coverage")
#             if min(abs(cms2_lats - lat)) > 0.5 or _min_lon_dist(cms2_lons, lon) > 0.5 or min(abs(cms2_time - time)) > np.timedelta64(3,'h'):
#                 raise ValueError("cms2_coverage")
#     except ValueError:
#         # Attempt to expand coverage and re-download once
#         cms_dataset.close(); cms_swh_dataset.close()
#         # Expand buffers and request files using the union of measurement and existing CMS extents (dateline-safe)
#         _request_if_needed_arrays(
#             np.concatenate([np.asarray(meas_lats), np.asarray(cms1_lats), np.asarray(cms2_lats)]),
#             np.concatenate([np.asarray(meas_lons), np.asarray(cms1_lons), np.asarray(cms2_lons)]),
#             np.concatenate([np.asarray(meas_times), np.asarray(cms1_time), np.asarray(cms2_time)]),
#             storm_id, lat_buffer=5, lon_buffer=5, hours_before=12, hours_after=48)
#         # Re-open files and refresh variables
#         cms_dataset = nc.Dataset(r'E:\Phd_data\copernicus\{}_SSS&SST.nc'.format(storm_id))
#         cms_swh_dataset = nc.Dataset(r'E:\Phd_data\copernicus\{}_SWH.nc'.format(storm_id))
#         cms1_lats = cms_dataset.variables['latitude'][:]
#         cms1_lons = cms_dataset.variables['longitude'][:]
#         tvar1 = cms_dataset.variables['time']
#         cms1_time = nc.num2date(tvar1[:], tvar1.units, getattr(tvar1, 'calendar', 'standard'))
#         cms1_time = np.array(cms1_time, dtype='datetime64[s]')
#         cms_sst = cms_dataset.variables['thetao'][:]  # sea surface temperature
#         cms_sss = cms_dataset.variables['so'][:]      # sea surface salinity

#         cms2_lats = cms_swh_dataset.variables['latitude'][:]
#         cms2_lons = cms_swh_dataset.variables['longitude'][:]
#         tvar2 = cms_swh_dataset.variables['time']
#         cms2_time = nc.num2date(tvar2[:], tvar2.units, getattr(tvar2, 'calendar', 'standard'))
#         cms2_time = np.array(cms2_time, dtype='datetime64[s]')
#         cms_swh = cms_swh_dataset.variables['VHM0'][:]  # significant wave height
#         cms_dir = cms_swh_dataset.variables['VMDR'][:]  # mean wave direction

#         # Re-normalize longitudes after re-opening to keep lon conventions consistent
#         cms1_lons = _normalize_to_pm180(cms1_lons)
#         cms2_lons = _normalize_to_pm180(cms2_lons)
#         # Ensure measurement longitudes are in same convention as CMS lon arrays
#         if np.nanmax(cms1_lons) <= 180:
#             meas_lons = ((np.asarray(meas_lons) + 180) % 360) - 180

#         # If still out of coverage then raise informative error
#         for i, (lat, lon, time) in enumerate(zip(meas_lats, meas_lons, meas_times)):
#             if min(abs(cms1_lats - lat)) > 0.5 or _min_lon_dist(cms1_lons, lon) > 0.5 or min(abs(cms1_time - time)) > np.timedelta64(24,'h'):
#                 raise ValueError(
#                     f"Measurement point ({lat}, {lon}) is outside the CMS1 SST/SSS data coverage even after expanding the domain. "
#                     f"CMS1 lon range: {np.nanmin(cms1_lons):.2f}–{np.nanmax(cms1_lons):.2f}, lat range: {np.nanmin(cms1_lats):.2f}–{np.nanmax(cms1_lats):.2f}, "
#                     f"time range: {np.nanmin(cms1_time)}–{np.nanmax(cms1_time)}")
#             if min(abs(cms2_lats - lat)) > 0.5 or _min_lon_dist(cms2_lons, lon) > 0.5 or min(abs(cms2_time - time)) > np.timedelta64(3,'h'):
#                 raise ValueError(
#                     f"Measurement point ({lat}, {lon}) is outside the CMS2 SWH data coverage even after expanding the domain. "
#                     f"CMS2 lon range: {np.nanmin(cms2_lons):.2f}–{np.nanmax(cms2_lons):.2f}, lat range: {np.nanmin(cms2_lats):.2f}–{np.nanmax(cms2_lats):.2f}, "
#                     f"time range: {np.nanmin(cms2_time)}–{np.nanmax(cms2_time)}")

#     # load up the values from the netcdf that match lat/lon/time
#     cms_vals = []
#     for i, (lat, lon, time) in enumerate(zip(meas_lats, meas_lons, meas_times)):
#         # Find the closest CMS grid point - Can improve here as CMS uses regular grid. 
#         lat_idx = np.argmin(np.abs(cms1_lats - lat))
#         lon_idx = _closest_lon_idx(cms1_lons, lon)
#         time_idx = np.argmin(np.abs(cms1_time - time))

#         # Extract the CMS values for this point and time
#         sst_val = cms_sst[time_idx, 0, lat_idx, lon_idx]
#         sss_val = cms_sss[time_idx, 0, lat_idx, lon_idx]

#         # Find the closest CMS SWH grid point
#         lat_idx2 = np.argmin(np.abs(cms2_lats - lat))
#         lon_idx2 = _closest_lon_idx(cms2_lons, lon)
#         time_idx2 = np.argmin(np.abs(cms2_time - time))

#         # Extract the CMS SWH value for this point and time
#         swh_val = cms_swh[time_idx2, lat_idx2, lon_idx2]
#         dir_val = cms_dir[time_idx2, lat_idx2, lon_idx2]

#         cms_vals.append((sst_val, sss_val, swh_val, dir_val))
#     return np.array(cms_vals)


def era5_downloader(minimum_longitude,  maximum_longitude,  minimum_latitude,  maximum_latitude,start_datetime, end_datetime,storm_id,overide=False):
    # ERA5 data downloader
    import cdsapi
    import datetime
    import os
     # check if file already exists
    filepath_era5 = f"E:\\Phd_data\\copernicus\\{storm_id}_era5.grib"
    if overide == True:
        # delete existing file
        if os.path.exists(filepath_era5):
            os.remove(filepath_era5)
    if os.path.exists(filepath_era5):
        print(f"ERA5 file for {storm_id}  already exists. Skipping download.")
        return  # Skip download if file already exists
    # Generate list of dates
    date_list = []
    current_date = start_datetime
    while current_date <= end_datetime+datetime.timedelta(days=2):
        date_list.append(current_date.strftime("%Y-%m-%d"))
        current_date += datetime.timedelta(days=1)

    # round coordinate to integer values as required by CDS API
    minimum_longitude = int(round(minimum_longitude))-5
    maximum_longitude = int(round(maximum_longitude))+5
    minimum_latitude = int(round(minimum_latitude))-5
    maximum_latitude = int(round(maximum_latitude))+5 

    # ERA5 data downloader
    dataset = "reanalysis-era5-single-levels"
    request = {
        "product_type": ["reanalysis"],
        "variable": [
            "10m_u_component_of_wind",
            "10m_v_component_of_wind",
            "2m_temperature",
            "mean_sea_level_pressure",
            "surface_pressure",
            "total_precipitation",
        ],
        "date": date_list,
        "time": [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ],
        "data_format": "grib",
        "download_format": "unarchived",
        "area": [maximum_latitude, minimum_longitude, minimum_latitude, maximum_longitude]
    }

    client = cdsapi.Client()
    client.retrieve(dataset, request, filepath_era5)

def retrieve_vals(meas_lats, meas_lons, meas_times,storm_id,file_name,variables):
    # a function that returns ERA5 values for each triplet of CYGNSS measurement lat, lon, time
    import numpy as np
    import cfgrib
    import netCDF4 as nc
    from cfgrib.dataset import DatasetBuildError
    
    if '.grib' in file_name:
        file_type = 'ERA5'
    elif '.nc' in file_name:
        file_type = 'CMS' 
    # Open the file
    if file_type == 'ERA5':
        try:
            dataset = cfgrib.open_dataset(file_name)
        except DatasetBuildError as e:
            # raise a concise error to avoid verbose cfgrib internal tracebacks
            raise RuntimeError(f"Failed to open ERA5 GRIB for {storm_id}: {e}") from None
        lat = dataset.latitude.values ; lon  = dataset.longitude.values
        time = dataset.time.values.astype('datetime64[h]')
        data = {}
        for var in variables:
            data[var] = dataset[var].values # I dont think this will work. 
            # mslp=era5_dataset.msl.values
            # u_wind = era5_dataset.u10.values ; v_wind=era5_dataset.v10.values
            # temp_2m= era5_dataset.t2m.values-273 ; surface_pressure= era5_dataset.sp.values
            # era5_wind = np.hypot(np.abs(u_wind), np.abs(v_wind))
    elif file_type == 'CMS':
        try:
            dataset = nc.Dataset(file_name)
        except DatasetBuildError as e:
            # raise a concise error to avoid verbose cfgrib internal tracebacks
            raise RuntimeError(f"Failed to open CMS nc for {storm_id}: {e}") from None
        lat = dataset.variables['latitude'][:]
        lon = dataset.variables['longitude'][:]
        tvar = dataset.variables['time']
        time = nc.num2date(tvar[:], tvar.units, getattr(tvar, 'calendar', 'standard'))
        time = np.array(time, dtype='datetime64[h]')
        data = {}
        for var in variables:
            # Extract variable data (assuming dimensions are time, lat, lon)
            values = dataset.variables[var][:]
            # remove the single depth dimension if it exists (CMS SSS/SST files have a depth dimension of size 1)
            if values.ndim == 4 and values.shape[1] == 1:
                values = values[:, 0, :, :]
            data[var] = values
    dataset.close()

    # ensure lon is in -180..180 to match CMS convention and handle dateline-crossing correctly
    lon = ((lon + 180) % 360) - 180

    # Quick coverage check and one retry with expanded buffer if needed
    # Use dataset time step to set a reasonable tolerance (CMS can be daily)
    if ['VHM0', 'VMDR'] == variables:  # CMS SWH files have 3h time steps
        time_tol = np.timedelta64(3, 'h')
    elif ['thetao', 'so'] == variables:  # CMS SSS/SST files have daily time steps
        time_tol = np.timedelta64(24, 'h')
    else:  # ERA5 has hourly time steps
        time_tol = np.timedelta64(2, 'h')
    meas_lats_arr = np.asarray(meas_lats)
    meas_lons_arr = np.asarray(meas_lons)
    meas_times_arr = np.asarray(meas_times)

    # Per-measurement minima, then max across measurements
    lat_min_list = np.min(np.abs(lat[:, np.newaxis] - meas_lats_arr[np.newaxis, :]), axis=0)
    lon_min_list = np.min(np.abs(lon[:, np.newaxis] - meas_lons_arr[np.newaxis, :]), axis=0)
    time_min_list = np.min(np.abs(time[:, np.newaxis] - meas_times_arr[np.newaxis, :]), axis=0)
    lat_min = np.max(lat_min_list)
    lon_min = np.max(lon_min_list)
    time_min = np.max(time_min_list).astype('timedelta64[h]')
    if lat_min > 0.5 or lon_min > 0.5 or time_min > time_tol:
        raise ValueError(f'Dataset not fully covering measurements. '
                       f'lat_min={lat_min}, lon_min={lon_min}, time_min={time_min}, '
                       f'file={file_name}, storm_id={storm_id}. '
                       f'Out of bounds: lat={lat_min > 0.5}, lon={lon_min > 0.5}, time={time_min > time_tol}')

    # Vectorized: compute all indices at once instead of looping through each measurement
    meas_lats = meas_lats_arr
    meas_lons = meas_lons_arr
    meas_times = meas_times_arr
    
    # Find closest indices for all measurements at once using broadcasting
    lat_idx = np.argmin(np.abs(lat[:, np.newaxis] - meas_lats[np.newaxis, :]), axis=0)
    lon_idx = np.argmin(np.abs(lon[:, np.newaxis] - meas_lons[np.newaxis, :]), axis=0)
    time_idx = np.argmin(np.abs(time[:, np.newaxis] - meas_times[np.newaxis, :]), axis=0)
    
    # Extract values for all measurements at once
    result = {}
    for var in variables:
        result[var] = data[var][time_idx, lat_idx, lon_idx]

    return result

def expand_and_redownload(error, filepath, storm_id, downloader_func, file_type='netcdf'):
    """
    Parse a ValueError from retrieve_vals, read existing file bounds, expand only the
    limiting dimensions, and re-download the file.
    
    Parameters:
    -----------
    error : ValueError
        The error raised by retrieve_vals
    filepath : str
        Path to the existing file
    storm_id : str
        Storm identifier
    downloader_func : callable
        Function to call for re-downloading (cms_downloader or era5_downloader)
    file_type : str
        'netcdf' for CMS files or 'grib' for ERA5 files
        
    Returns:
    --------
    tuple : (new_lon_min, new_lon_max, new_lat_min, new_lat_max, new_time_start, new_time_end)
    """
    import os
    from datetime import datetime, timedelta
    import netCDF4 as nc
    import numpy as np
    
    error_msg = str(error)
    lat_oob = 'lat=True' in error_msg
    lon_oob = 'lon=True' in error_msg
    time_oob = 'time=True' in error_msg
    
    print(f"{file_type.upper()} coverage issue for {storm_id}: lat_out_of_bounds={lat_oob}, lon={lon_oob}, time={time_oob}")
    
    # Read existing file bounds
    if file_type == 'netcdf':
        dataset = nc.Dataset(filepath)
        existing_lat_min = float(dataset.variables['latitude'][:].min())
        existing_lat_max = float(dataset.variables['latitude'][:].max())
        existing_lon_min = float(dataset.variables['longitude'][:].min())
        existing_lon_max = float(dataset.variables['longitude'][:].max())
        tvar = dataset.variables['time']
        existing_times = nc.num2date(tvar[:], tvar.units, getattr(tvar, 'calendar', 'standard'))
        existing_times = np.array(existing_times, dtype='datetime64[s]')
        existing_time_start = np.min(existing_times)
        existing_time_end = np.max(existing_times)
        # Convert numpy datetime64 to Python datetime for timedelta arithmetic
        existing_time_start = datetime.utcfromtimestamp(
            existing_time_start.astype('datetime64[s]').astype('int64')
        )
        existing_time_end = datetime.utcfromtimestamp(
            existing_time_end.astype('datetime64[s]').astype('int64')
        )
        dataset.close()
    else:  # grib (ERA5)
        import cfgrib
        import pandas as pd
        dataset = cfgrib.open_dataset(filepath)
        existing_lat_min = float(dataset.latitude.values.min())
        existing_lat_max = float(dataset.latitude.values.max())
        existing_lon_min = float(dataset.longitude.values.min())
        existing_lon_max = float(dataset.longitude.values.max())
        existing_times = dataset.time.values.astype('datetime64[s]')
        existing_time_start = pd.Timestamp(existing_times.min()).to_pydatetime()
        existing_time_end = pd.Timestamp(existing_times.max()).to_pydatetime()
        dataset.close()
    
    # Expand only the limiting dimensions
    new_lat_min = existing_lat_min - 4 if lat_oob else existing_lat_min
    new_lat_max = existing_lat_max + 4 if lat_oob else existing_lat_max
    new_lon_min = existing_lon_min - 4 if lon_oob else existing_lon_min
    new_lon_max = existing_lon_max + 4 if lon_oob else existing_lon_max
    new_time_start = existing_time_start - timedelta(hours=24) if time_oob else existing_time_start
    new_time_end = existing_time_end + timedelta(hours=48) if time_oob else existing_time_end
    
    print(f"Expanding {file_type.upper()} bounds: "
          f"lat [{existing_lat_min:.2f}, {existing_lat_max:.2f}] -> [{new_lat_min:.2f}, {new_lat_max:.2f}], "
          f"lon [{existing_lon_min:.2f}, {existing_lon_max:.2f}] -> [{new_lon_min:.2f}, {new_lon_max:.2f}], "
          f"time [{existing_time_start}] -> [{new_time_start}]")
    
    # Delete and re-download
    if os.path.exists(filepath):
        os.remove(filepath)
    
    downloader_func(new_lon_min, new_lon_max, new_lat_min, new_lat_max, new_time_start, new_time_end, storm_id, overide=True)
    
    return new_lon_min, new_lon_max, new_lat_min, new_lat_max, new_time_start, new_time_end

# Wave data downloader
def SWH_L3_searcher(minimum_longitude,  maximum_longitude,  minimum_latitude,  maximum_latitude,start_datetime, end_datetime,storm_id, tc_event, overide=False):
    # download the SWH L3 data for the given storm and tc event. The tc_event is used to determine which dataset to download as the datasets have changed over time. The function checks if the file already exists and skips download if it does unless overide is True. The function handles both old and new datasets based on the start_datetime of the storm event.
    import copernicusmarine
    import datetime
    import os

    datasets_new=['cfo','c2','h2b','h2c','j3','al','s3a','s3b','s6a','swon']
    datasets_old = ['cfo','c2','en','j1','j2','j3','al']
    # datasets=['cfo']
    
    if start_datetime < datetime.datetime(2021,1,1):
        datasets = datasets_old
    else:
        datasets = datasets_new
    for data in datasets:
        # check if file already exists
        output_filename = f"{storm_id}_{tc_event}_{data}_l3_swh.nc"
        if not overide and os.path.exists(os.path.join(r'E:\Phd_data\copernicus\\', output_filename)):
            continue
        if data == 'cfo' and start_datetime < datetime.datetime(2021,1,1):
            product_code = 'cmems_obs-wave_glo_phy-swh_my_'+data+'-l3_PT1S'
        elif data in ['c2','en','j1','j2','j3','al'] and start_datetime < datetime.datetime(2021,1,1):
            product_code = 'cci_obs-wave_glo_phy-swh_my_'+data+'-l3_PT1S'
        else:
            product_code = 'cmems_obs-wave_glo_phy-swh_nrt_'+data+'-l3_PT1S'
        try:
            copernicusmarine.subset(
                dataset_id=product_code,
                variables=["VAVH","VAVH_UNFILTERED","WIND_SPEED"],
                username='asiedlecki',
                password='XsySq5LJ8TmZKs',
                minimum_longitude=minimum_longitude,
                maximum_longitude=maximum_longitude,
                minimum_latitude=minimum_latitude,
                maximum_latitude=maximum_latitude,
                start_datetime=start_datetime,
                end_datetime=end_datetime,
                netcdf_compression_level=1,
                output_filename = output_filename,
                output_directory = r'E:\Phd_data\copernicus'
            )
        except copernicusmarine.core_functions.exceptions.CoordinatesOutOfDatasetBounds:
            continue

# # ERA5 data downloader
# import cdsapi

# dataset = "reanalysis-era5-single-levels"
# request = {
#     "product_type": ["reanalysis"],
#     "variable": [
#         "10m_u_component_of_wind",
#         "10m_v_component_of_wind",
#         "2m_temperature",
#         "mean_sea_level_pressure",
#         "mean_wave_period",
#         "significant_height_of_combined_wind_waves_and_swell",
#         "surface_pressure",
#         "total_precipitation",
#         "ocean_surface_stress_equivalent_10m_neutral_wind_speed",
#         "significant_height_of_wind_waves",
#         "maximum_total_precipitation_rate_since_previous_post_processing"
#     ],
#     "year": ["2022"],
#     "month": ["08"],
#     "day": ["28", "29", "30"],
#     "time": [
#         "00:00", "01:00", "02:00",
#         "03:00", "04:00", "05:00",
#         "06:00", "07:00", "08:00",
#         "09:00", "10:00", "11:00",
#         "12:00", "13:00", "14:00",
#         "15:00", "16:00", "17:00",
#         "18:00", "19:00", "20:00",
#         "21:00", "22:00", "23:00"
#     ],
#     "data_format": "grib",
#     "download_format": "unarchived",
#     "area": [28, 142, 26, 147]
# }

# client = cdsapi.Client()
# client.retrieve(dataset, request).download()