import netCDF4 as nc
import numpy as np
import glob
import os
import pandas as pd
from ndbc_api import NdbcApi


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

def cms_insitu_download(vars_list,cms_path,lon_min, lon_max, lat_min, lat_max, time_start, time_end):
            copernicusmarine.subset(
                dataset_id="cmems_obs-ins_glo_phybgcwav_mynrt_na_irr",
                dataset_part="monthly",
                dataset_version="202311",
                variables=vars_list,
                minimum_longitude=lon_min,
                maximum_longitude=lon_max,
                minimum_latitude=lat_min,
                maximum_latitude=lat_max,
                start_datetime=time_start,
                end_datetime=time_end,
                minimum_depth=-20,
                maximum_depth=20,
                coordinates_selection_method="strict-inside",
                disable_progress_bar=False,
                username='asiedlecki',
                password='XsySq5LJ8TmZKs',
                output_filename = cms_path,
                output_directory = r'E:\Phd_data\casestudy'
            )

if __name__ == "__main__":
    # find case events with in-situ data and list the events with the most data to prioritise which ones to analyse first
    comp = 'ashle'
    buffer = 1 
    interp_interval = 20
    phase = 1 # phase 1 is creating the list of events with some data. Phase 2 is downloading and comparing data for those events to identify which ones are the best for case studies. 

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

    # create txt files to save when a case study is completed and which events have in-situ data available for future reference
    insitu_file = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_with_insitu.txt'
    # create a list of file names in this folder
    gifs_folder = r"C:\Users\\" + comp + r"\OneDrive - RMIT University\PHD\Plots\gifs"
    existing_gifs = glob.glob(gifs_folder + r"\NOAA_*.gif")
    # select the index in the middle of each file name between the '_'s and save index in a list
    existing_gif_indices = np.unique([int(os.path.basename(gif).split('_')[1].split('.')[0]) for gif in existing_gifs])
    # existing_gif_indices = [200,210]

    if phase == 1:
        for tc_event in existing_gif_indices:
            if ri_rd[tc_event]==1 and durations[tc_event]>13 :
                tc_latitude = latitude[tc_event]
                tc_longitude = longitude[tc_event]
                tc_time = time[tc_event]
                valid_tc_time = tc_time[~np.isnat(tc_time)]
                if valid_tc_time.size == 0:
                    continue

                # linearly interpolate event to X minutes 
                interval = np.timedelta64(interp_interval, 'm')
                start = valid_tc_time[0] ; end = valid_tc_time[-1]
                tc_time_interpolated = np.arange(start, end + interval, interval)
                        # Vectorized interpolation using seconds-since-epoch (faster than looping)
                tc_times_s = tc_time.astype('datetime64[s]').astype('int64')
                interp_s = tc_time_interpolated.astype('datetime64[s]').astype('int64')

                # Ensure numeric arrays (fill masked vmax with nan)
                tc_lat_vals = np.asarray(tc_latitude).astype(float)
                tc_lon_vals = np.asarray(tc_longitude).astype(float)

                # Interpolate (np.interp extrapolates using endpoint values similar to previous logic)
                tc_lat_interpolated = np.interp(interp_s, tc_times_s, tc_lat_vals)
                tc_lon_interpolated = np.interp(interp_s, tc_times_s, tc_lon_vals)

                cyg_AOI_time_start, cyg_AOI_time_end = np.min(valid_tc_time), np.max(valid_tc_time)
                cyg_AOI_time_start = cyg_AOI_time_start.astype('datetime64[s]').astype('O')
                cyg_AOI_time_end = cyg_AOI_time_end.astype('datetime64[s]').astype('O')
                AOI_lat_min,AOI_lat_max = np.nanmin(tc_latitude)-buffer,np.nanmax(tc_latitude)+buffer
                AOI_lon_min,AOI_lon_max = np.nanmin(tc_longitude)-buffer,np.nanmax(tc_longitude)+buffer
                tc_coords = np.column_stack((tc_lat_interpolated, tc_lon_interpolated))  # Shape (M, 2)
                
                # in-situ data
                # in-situ, load from the all vars folder, and place all data from both cms and ndbc into single pandas dataframe
                insitu_dfs = []  
                api = NdbcApi()
                # specify desired latitude, longitude, radius, and units
                # find the station IDs of all NDBC stations within the radius
                # calculate radius in km between maximum values of the TC track 
                radius = np.max(haversine_vectorized(tc_coords, tc_coords))
                nearby_stations_df = api.radial_search(lat=np.nanmean(tc_latitude), lon=np.nanmean(tc_longitude), radius=radius, units='km')
                if len(nearby_stations_df)!=0:
                    # check if any of the buoys are within 100 km of the TC track
                    stations_lats = nearby_stations_df['Lat'].values ; stations_lons = nearby_stations_df['Lon'].values
                    stations_coords = np.column_stack((np.asarray(stations_lats), np.asarray(stations_lons)))  # Shape (N, 2)
                    distances_to_track = haversine_vectorized(stations_coords, tc_coords)
                    min_distances = np.min(distances_to_track, axis=1)
                    nearby_stations_df['Min_Distance_to_Track_km'] = min_distances  

                    # select which station to download data from (e.g., the ones within 100km distance)
                    nearby_stations_df = nearby_stations_df[nearby_stations_df['Min_Distance_to_Track_km'] < 100]
                    stations_list = nearby_stations_df['Station'].tolist()
                    if len(stations_list) > 0:
                        stdmet_df = api.get_data(
                            station_ids=stations_list,
                            modes=['stdmet', 'cwind','ocean','spec'],
                            start_time=cyg_AOI_time_start,
                            end_time=cyg_AOI_time_end,
                            as_df=True
                        )
                        # save the data to a csv file for future reference and to avoid having to re-download it, and add the event to the insitu file
                        if isinstance(stdmet_df, dict):
                            stdmet_df = pd.DataFrame(stdmet_df)
                        insitu_dfs.append(stdmet_df)
                        
                        # save the insitu data to a file for future reference
                        noaa_file = r'E:\Phd_data\casestudy\all_vars\\' + f"{storm_id[tc_event]}_ndbc_insitu.csv"
                        stdmet_df.to_csv(noaa_file, index=False)


                # copernicus in-situ data with retry: try small var list then full list
                if tc_time[0] < np.datetime64('2020-01-01'):
                    continue  # Skip events before 2020 due to data availability issues
                import copernicusmarine
                cms_small_path =  f"{storm_id[tc_event]}_{tc_event}_cms_insitu.csv"
                cms_all_path = f"all_vars\\{storm_id[tc_event]}_{tc_event}_cms_insitu.csv"
                cms_dir = r'E:\Phd_data\casestudy\\'
                small_vars = ["ATMS", "DRYT"]
                all_vars = ["ATMP", "ATMS", "ATPT", "DEWT", "DRYT", "GDIR", "GSPD", "NRAD", "RELH", "SLEV", "TEMP", "VGTA", "VHZA", "VTMX", "VTZM", "WDIR", "WETT", "WSPD", "WSPN", "VCMX", "PRRT", "PSAL", "DENS", "EWCT", "NSCT", "VHM0", "VTPK", "VMDR"]

                # look for small var list first as a quick check for data, then if not try the full list (some events have issues with some variables but may still have some valid data)
                try:
                    if not glob.glob(cms_dir+cms_small_path):
                        cms_insitu_download(small_vars,cms_small_path,AOI_lon_min, AOI_lon_max, AOI_lat_min, AOI_lat_max, cyg_AOI_time_start, cyg_AOI_time_end)
                    cms_df = pd.read_csv(cms_dir+cms_small_path)
                    if cms_df.empty or 'latitude' not in cms_df.columns or 'longitude' not in cms_df.columns:
                        raise pd.errors.EmptyDataError
                except (pd.errors.EmptyDataError, FileNotFoundError, OSError):
                    try:
                        if os.path.exists(cms_dir+cms_small_path):
                            os.remove(cms_dir+cms_small_path)
                    except Exception:
                        pass
                    try:
                        if not os.path.exists(cms_dir+cms_all_path):
                            cms_insitu_download(all_vars, cms_all_path, AOI_lon_min, AOI_lon_max, AOI_lat_min, AOI_lat_max, cyg_AOI_time_start, cyg_AOI_time_end)
                        cms_df = pd.read_csv(cms_dir+cms_all_path)
                        if cms_df.empty or 'latitude' not in cms_df.columns or 'longitude' not in cms_df.columns:
                            raise Exception('No valid CMS data')
                    except Exception:
                        print(f"No CMS data found for {storm_id[tc_event]}; skipping this event.")
                        continue

                cms_lats = cms_df['latitude'].values ; cms_lons = cms_df['longitude'].values
                cms_coords = np.column_stack((np.asarray(cms_lats), np.asarray(cms_lons)))  # Shape (N, 2)
                distances_to_track = haversine_vectorized(cms_coords, tc_coords)
                min_distances = np.min(distances_to_track, axis=1)
                cms_df['Min_Distance_to_Track_km'] = min_distances
                cms_df = cms_df[cms_df['Min_Distance_to_Track_km'] < 100]
                if len(cms_df)>0:
                    # check if cms_all_path already exists. If not, download the full var list and save the data, if it does exist, it means the small var list had some valid data but we want to get the full var list for the events with some data available
                    if not os.path.exists(cms_dir+cms_all_path):
                        cms_insitu_download(all_vars, cms_all_path, min(cms_df['longitude']), max(cms_df['longitude']), min(cms_df['latitude']), max(cms_df['latitude']), min(cms_df['time']), max(cms_df['time']))
                    cms_df = pd.read_csv(cms_dir+cms_all_path)
                    cms_lats = cms_df['latitude'].values ; cms_lons = cms_df['longitude'].values
                    cms_coords = np.column_stack((np.asarray(cms_lats), np.asarray(cms_lons)))  # Shape (N, 2) 
                    distances_to_track = haversine_vectorized(cms_coords, tc_coords)
                    min_distances = np.min(distances_to_track, axis=1)
                    cms_df['Min_Distance_to_Track_km'] = min_distances
                    cms_df = cms_df[cms_df['Min_Distance_to_Track_km'] < 100]
                    insitu_dfs.append(cms_df)
                    # with open(insitu_file, 'a') as f:
                    #     f.write(f"{storm_id[tc_event]}, {tc_event}\n")
                if len(insitu_dfs)>0:
                    combined_insitu_df = pd.concat(insitu_dfs, ignore_index=True)
                    combined_insitu_df.to_csv(r'E:\Phd_data\casestudy\\' + f"{storm_id[tc_event]}_buoy_data.csv", index=False)


                # copernicus L3 SWH data via API
                download_swh = False
                if download_swh:
                    from CMS_CDS_api import SWH_L3_searcher
                    SWH_L3_searcher(AOI_lon_min,  AOI_lon_max,  AOI_lat_min,  AOI_lat_max,cyg_AOI_time_start, cyg_AOI_time_end,storm_id[tc_event],tc_event, overide=False)
                    L3_files = glob.glob(r'E:\Phd_data\casestudy\L3_SWH\\' + f"{storm_id[tc_event]}_*_swh.nc")
                    #open each of the files and append the data to a single dataframe, then save the dataframe as a csv file
                    L3_dfs = []
                    for L3_file in L3_files:
                        L3_nc = nc.Dataset(L3_file)
                        L3_time_int = L3_nc.variables['time'][:]
                        L3_time = np.datetime64('1970-01-01T00:00:00') + L3_time_int.astype('timedelta64[s]')
                        L3_latitude = L3_nc.variables['latitude'][:]
                        L3_longitude = L3_nc.variables['longitude'][:]
                        L3_swh = L3_nc.variables['swh'][:]
                        L3_df = pd.DataFrame({'time': L3_time, 'latitude': L3_latitude, 'longitude': L3_longitude, 'swh': L3_swh})
                        # add a column for the data source
                        L3_df['source'] = L3_file.split('\\')[-1].split('_')[1]  # extract the source from the file name
                        L3_dfs.append(L3_df)
                        # Combine all dataframes into a single dataframe
                    if len(L3_dfs)>0:
                        combined_L3_df = pd.concat(L3_dfs, ignore_index=True)  
                        # find measurements that are within 100 km of TC track and measured within 3 hours of the TC time, and add a column for whether the measurement is within these thresholds or not
                        meas_coords = np.column_stack((combined_L3_df['latitude'].values, combined_L3_df['longitude'].values))  # Shape (N, 2)
                        # Compute the pairwise distance matrix
                        distances = haversine_vectorized(meas_coords, tc_coords)
                        # Compute pairwise time difference matrix (absolute time difference in hours)
                        time_diffs = np.abs(combined_L3_df['time'].values[:, None] - tc_time[None, :])
                        # Find measurements where at least one track point is within 100 km & 3 hours
                        # print(np.timedelta64(np.min(time_diffs),'m'))
                        valid_mask = (np.abs(distances) <= 100) & (time_diffs <= np.timedelta64(180, 'm'))
                        # Select valid measurement indices
                        valid_ind = np.where(valid_mask)
                        selected_indices = np.unique(valid_ind[0]) # to remove duplicates
                        
                        # Add a column to indicate whether the measurement is within the thresholds
                        combined_L3_df['within_thresholds'] = False
                        combined_L3_df.loc[selected_indices, 'within_thresholds'] = True
                        
                        # Save the combined dataframe with the new column
                        combined_L3_df.to_csv(r'E:\Phd_data\casestudy\\' + f"{storm_id[tc_event]}_swh_with_thresholds.csv", index=False)
                        
                        # print number of valid measurements within the thresholds
                        num_valid_measurements = combined_L3_df['within_thresholds'].sum()
                        print(f"{storm_id[tc_event]} has {num_valid_measurements} valid SWH measurements within 100 km and 3 hours of the TC track.")


    # evaluating each pf the case studies based on their data availability and proximity to the track, to prioritise which ones to analyse first (e.g., which ones have the most data within 100km of the track, which ones have data closest to the track, which ones have data closest to the time of RI, etc.)
    if phase == 2:
        # now that data has been downloaded, want to analyse to see which has the most meaningful data 
        cms_all_folder = r'E:\Phd_data\casestudy\\all_vars\\' 
        cms_files = glob.glob(cms_all_folder + r"\*_cms_insitu.csv")
        noaa_files = glob.glob(r'E:\Phd_data\casestudy\all_vars\*_ndbc_insitu.csv')
        for cms_file in cms_files:
            cms_df = pd.read_csv(cms_file)
            if cms_df.empty or 'latitude' not in cms_df.columns or 'longitude' not in cms_df.columns:
                continue
            # check if there are any valid data points within 100km of the track, and if so, how many variables have valid data at those points (some events have a lot of data but most of it is missing values, so want to prioritise events with more valid data)
            var_count = cms_df['variable'].unique().shape[0]  # Count unique variables with valid data
            min_distance = cms_df['Min_Distance_to_Track_km'].min()
            print(f"{os.path.basename(cms_file)} has {var_count} unique variables and a minimum distance of {min_distance} km.")

        for noaa_file in noaa_files:
            noaa_df = pd.read_csv(noaa_file)
            if noaa_df.empty or 'Lat' not in noaa_df.columns or 'Lon' not in noaa_df.columns:
                continue
            # check how many data points there are within 100km of the track, and how many variables have valid data at those points
            var_count = noaa_df['var'].unique().shape[0]  # Count unique variables with valid data
            min_distance = noaa_df['Min_Distance_to_Track_km'].min()
            print(f"{os.path.basename(noaa_file)} has {var_count} unique variables and a minimum distance of {min_distance} km.")