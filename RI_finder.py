# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 07:21:02 2024
Purpose: determine best RI events in best-track data
Input: IBTrACS dataset
Output: List of RI events, time and lat and lon. 
@author: ashle
"""
import netCDF4 as nc
import numpy as np 

# User identifier for file path construction
comp = 'S3987712'

# mode - either create or analyse. Create will run the code to find RI events and save to a new netCDF file. Analyse will load the saved file and run some basic analysis on it.
mode = 'create'# 'create' # 'analyse'

version = '6'

# Time intervals (in hours) to check for RI events
time_thresholds = [6, 12, 24]
# Wind speed thresholds (in knots) for identifying RI intensity changes
wind_thresholds = [45, 30]

# Lists to store identified RI event data
storm_ids = []
storm_names = []
inds_before = []  # Index before the RI event for context
inds_after = []   # Index after the RI event for context
times = []        # Timestamps for each RI event
lats = []         # Latitude coordinates
lons = []         # Longitude coordinates
winds = []        # Wind speeds throughout the event
ri_or_rd = []     # Classification: 1 for RI (intensification), -1 for RD (decay)
duration = []     # Duration of the event in hours
intensity = []    # Wind speed threshold met for the event

path = r'C:\Users\\' +comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS'
file_name = r'\IBTrACS.since1980.v04r01.nc' #'\IBTrACS.last3years.v04r00.nc'
RI_file = r'C:\Users\\' + comp + rf'\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v{version}.nc'

if mode == 'create':
    IB_nc = nc.Dataset(path + file_name)
    WMO_wind = IB_nc.variables['wmo_wind'][:]
    wmo_agency = IB_nc.variables['wmo_agency']
    storm_id = IB_nc.variables['sid'][:].astype(str)
    storm_name = IB_nc.variables['name'][:].astype(str)
    dt_str = IB_nc.variables['iso_time'][:].astype('U2')
    lat = IB_nc.variables['lat'][:]
    lon = IB_nc.variables['lon'][:]

    rows, cols, depth = dt_str.shape
    dt_concatenated = np.empty((rows, cols))
    storm_id_concatenated = np.zeros((rows)).astype(str) 

    for i in range(rows):
        storm_id_concatenated[i] = np.array([''.join(storm_id[i,:])])[0]
        
        for j in range(cols):
            if not type(dt_str[i,j,0]) == np.ma.core.MaskedConstant:
                dt_concatenated[i,j]= np.array([''.join(dt_str[i,j,:])],dtype = 'datetime64[s]')[0]

    # find locations of RI
    from datetime import datetime, timedelta
    dt= dt_concatenated.astype(datetime)

    for time_threshold in time_thresholds:
        for tc in range(rows):
            start_ind = 0 ; end_ind = 1
            counter = 0
            year = np.datetime64('1970-01-01T00:00:00') + np.timedelta64(int(dt[tc,end_ind]), 's')
            if np.datetime64(year,'Y') > np.datetime64('2018'): # only run for Valid CYGNSS years. and np.datetime64(year,'Y') < np.datetime64('2027')
                while dt[tc,end_ind] > 1 and counter <cols*2: 
                    counter+=1
                    time_difference = (dt[tc,end_ind] - dt[tc,start_ind])/60/60 # convert from seconds to hours
                    if time_difference < time_threshold:
                        end_ind += 1
                        if end_ind >= cols: # if it makes it to the end of array
                            break 
                    elif time_difference > time_threshold:
                        start_ind =+ 1
            
                    elif time_difference == time_threshold and np.nan not in [WMO_wind[tc,end_ind], WMO_wind[tc,start_ind]]:
                        # if (storm_id_concatenated[tc] not in storm_ids): # or (storm_id_concatenated[tc] in storm_ids and dt[tc,start_ind] > times[-1][-1]): 
                        # now to compare speeds to see if threshold is met 
                        WMO_wind_dif = WMO_wind[tc,end_ind] - WMO_wind[tc,start_ind]
                        for wind_threshold in wind_thresholds:
                            if WMO_wind_dif >= wind_threshold: # or 
                                storm_ids.append(storm_id_concatenated[tc])  
                                
                                # update storm_name to be a simple string of the name instead of an array of characters
                                storm_name_clean = ''.join(
                                    str(char) for char in storm_name[tc]
                                    if not isinstance(char, np.ma.core.MaskedConstant)
                                ).strip()
                                storm_names.append(storm_name_clean)
                                # end_ind += 1 # to make sure we include the end measurement
                                start_ind_before = np.abs(dt[tc]-(dt[tc,start_ind]-43200)).argmin() # np.where(dt[tc]==dt[tc,start_ind]-43200)[0][0] # maybe adjust to for loop if not work and add in as many 
                                end_ind_after = np.abs(dt[tc]-(dt[tc,end_ind]+43200)).argmin() #np.where(dt[tc]==dt[tc,end_ind]+43200)[0][0]
                                ind_before = start_ind - start_ind_before ; ind_after = end_ind_after - end_ind
                                inds_before.append(ind_before); inds_after.append(ind_after)
                                times.append(dt[tc,start_ind_before:end_ind_after])
                                lats.append(lat[tc,start_ind_before:end_ind_after])
                                lons.append(lon[tc,start_ind_before:end_ind_after])
                                duration.append(time_difference)
                                intensity.append(wind_threshold)
                                winds.append(WMO_wind[tc,start_ind_before:end_ind_after])
                                
                                if WMO_wind_dif > 0: 
                                    ri_or_rd.append(1)
                                else:
                                    ri_or_rd.append(-1)
                                break
                                # if end_ind-start_ind >9:
                                #     print((tc,start_ind,end_ind))
                        start_ind += 1
                            
                    else:
                        start_ind += 1
        
    # saving the details

    def padding(data,max_len):
        # input data list of numpy arrays and adds numpy nan to the end so all the arrays are same length
        import numpy as np
        padded_data=np.ma.MaskedArray(np.full((len(data), max_len), np.nan))
        # padded_data.data =   # Create 2D array filled with NaN
        for i, arr in enumerate(data):
            padded_data[i, :len(arr)] = arr  # Fill in values from each array
            
        return padded_data
            
            
    import netCDF4 as nc
    with nc.Dataset(RI_file, 'w') as ncfile:
        max_len = max(len(arr) for arr in times)
        # Pad arrays to the same length with NaN
        times_padded = padding(times,max_len)  
        lats_padded = padding(lats,max_len)  
        lons_padded = padding(lons,max_len)  
        winds_padded = padding(winds,max_len)
        new_mask = np.isnan(winds_padded.data) | (winds_padded.data < 0)     # Create new condition mask: values < 0 or NaN
        winds_padded.mask = winds_padded.mask | new_mask     # Combine old and new masks
        # Create dimensions
        length = len(storm_ids)
        ncfile.createDimension('events', length)
        ncfile.createDimension('values', max_len)
        
        # Create variables
        time = ncfile.createVariable('times', 'f8', ('events','values'))
        storm = ncfile.createVariable('storm', str, ('events',))
        latitude = ncfile.createVariable('latitude', 'f4', ('events','values'))
        longitude = ncfile.createVariable('longitude', 'f4', ('events','values'))
        vmax = ncfile.createVariable('vmax', 'f4', ('events','values'))
        ri_rd = ncfile.createVariable('ri_rd', 'i4', ('events',))
        durations = ncfile.createVariable('durations', 'i4', ('events',))
        intensities = ncfile.createVariable('intensities', 'i4', ('events',))
        index_before = ncfile.createVariable('index_before', 'i4', ('events',))
        index_after = ncfile.createVariable('index_after', 'i4', ('events',))


        # Add data to variables
        time[:] = times_padded
        storm[:] = np.array(storm_ids,dtype="object")
        latitude[:] = lats_padded
        longitude[:] = lons_padded
        vmax[:] = winds_padded
        ri_rd[:] = ri_or_rd
        durations[:] = duration
        intensities[:] = intensity
        index_after[:] = inds_after
        index_before[:] = inds_before
    
    # save also as a csv file for easier access
    import pandas as pd
    df = pd.DataFrame({
        'storm_id': storm_ids,
        'storm_name': storm_names,
        'times': times,
        'latitude': lats,
        'longitude': lons,
        'vmax': winds,
        'ri_rd': ri_or_rd,
        'durations': duration,
        'intensities': intensity,
        'index_before': inds_before,
        'index_after': inds_after
    })
    df.to_csv(path + f'RI_events_v{version}.csv', index=False)

if mode == 'analyse':
    RI_nc = nc.Dataset(RI_file)
    times = RI_nc.variables['times'][:]
    storm = RI_nc.variables['storm'][:]
    latitude = RI_nc.variables['latitude'][:]
    longitude = RI_nc.variables['longitude'][:]
    vmax = RI_nc.variables['vmax'][:]
    ri_rd = RI_nc.variables['ri_rd'][:]
    durations = RI_nc.variables['durations'][:]
    intensities = RI_nc.variables['intensities'][:]
    index_before = RI_nc.variables['index_before'][:]
    index_after = RI_nc.variables['index_after'][:]

    # print the number of total events and the number of RI vs RD events
    total_events = len(storm)
    ri_events = np.sum(ri_rd == 1)
    rd_events = np.sum(ri_rd == -1)
    print(f'Total events: {total_events}, RI events: {ri_events}, RD events: {rd_events}')

    # print the average duration and intensity of RI vs RD events
    avg_duration_ri = np.mean(durations[ri_rd == 1])
    avg_duration_rd = np.mean(durations[ri_rd == -1])
    avg_intensity_ri = np.mean(intensities[ri_rd == 1])
    avg_intensity_rd = np.mean(intensities[ri_rd == -1])
    print(f'Average duration of RI events: {avg_duration_ri} hours, Average duration of RD events: {avg_duration_rd} hours')
    print(f'Average intensity of RI events: {avg_intensity_ri} knots, Average intensity of RD events: {avg_intensity_rd} knots')

    # print the number of unique storms that experienced RI vs RD events at each time and intensity threshold
    for time_threshold in time_thresholds:
        for wind_threshold in wind_thresholds:
            ri_storms = set(storm[ri_rd == 1][(durations[ri_rd == 1] == time_threshold) & (intensities[ri_rd == 1] == wind_threshold)])
            rd_storms = set(storm[ri_rd == -1][(durations[ri_rd == -1] == time_threshold) & (intensities[ri_rd == -1] == wind_threshold)])
            print(f'Time threshold: {time_threshold} hours, Wind threshold: {wind_threshold} knots, Unique storms with RI events: {len(ri_storms)}')

    # determine how many storms cross the longitude of 180 and 360 during an RI event
    storms_crossing_180 = 0
    storms_crossing_360 = 0

    def check_crossing(longitudes):
        # Check if the storm crosses the 180 or 360 longitude
        crosses_180 = np.any(longitudes > 180) and np.any(longitudes < -180)
        crosses_360 = np.any(longitudes > 360) and np.any(longitudes < 0)
        return crosses_180, crosses_360

    for i in range(total_events):
        if ri_rd[i] == 1:
            crosses_180, crosses_360 = check_crossing(longitude[i,:])
            if crosses_180:
                storms_crossing_180 += 1
            if crosses_360:
                storms_crossing_360 += 1
    print(f'Storms crossing longitude 180 during RI events: {storms_crossing_180}')
    print(f'Storms crossing longitude 360 during RI events: {storms_crossing_360}')