# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 08:13:49 2024
BOM AWS and PSLGM analysis
Inputs: Not sure, CYGNSS L2 files
Outputs: R, RMSE etc 
@author: ashle
"""
import glob
from datetime import datetime, timedelta 
import numpy as np
import pandas as pd
import pickle
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
#from BOM_high_speed_finder import wind_thresh, time_thresh

# input files from PSLGM_high_speed_finder
comp = 'Ashley'
cyg_version = 'l2_3.2' # 'noaa_1.2' OR 'l2_3.2' 
data_folder = '\Aus 1hr' # '\BOM Stations' OR ' PSLGM' OR 'Aus 1hr'
csv_folder  = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\BOM stations' + data_folder
search_distance = 0.3 # number of gedrees for searching for CYGNSS
cyg_time_thresh = 1 # number of hours to look for CYGNSS outside of event 

max_num_matchups = 300000


cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\cyg_'+cyg_version
cyg_files_list = glob.glob(cyg_folder+ "\*.nc")

# find the dates of cyg 
cyg_file_dates = []
for cyg_file in cyg_files_list:
    start_ind = cyg_file.find('ddmi.s')
    cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
    cyg_file_dates.append(cyg_datetime_object.date())

if data_folder == '\PSLGM':
    csv_files_list = glob.glob(csv_folder+ "\*.csv")
else:
    csv_files_list = glob.glob(csv_folder+ "\HC06D_Data_*.txt")

location_dict = {
    'Cook Islands': {'lat': -21.199, 'lon': -159.786},
    'Federated States of Micronesia': {'lat': 6.978, 'lon': 158.197},
    'FSM': {'lat': 6.978, 'lon': 158.197},
    ' Fiji': {'lat': -17.605, 'lon': 177.438},
    'Kiribati': {'lat': 1.362, 'lon': 172.93},
    'Marshall Islands': {'lat': 7.108, 'lon': 171.371},
    'Nauru': {'lat': -0.532, 'lon': 166.909},
    'Papua New Guinea': {'lat': -2.036, 'lon': 147.375},
    'Samoa': {'lat': -13.819, 'lon': -171.756},
    'Solomon Islands': {'lat': -9.422, 'lon': 159.955},
    'Tonga': {'lat': -21.14, 'lon': -175.179},
    'Tuvalu': {'lat': -8.503, 'lon': 179.209},
    'Vanuatu': {'lat': -17.761, 'lon': 168.293},
    'Niue': {'lat': -19.053, 'lon': -169.921}
}

# Function for vectorized linear interpolation
def linear_interpolation_vectorized(timestamps, wind_speeds, target_times): # this needs to be improved
    indices = np.searchsorted(timestamps, target_times) - 1
    if indices == -1: # if it is before the first speed, take the first speed
        return wind_speeds[0]
    if indices == len(timestamps)-1: # if it is after the last reading, take the last speed 
        return wind_speeds[-1]
    else: # for all else inbetween 
        before_times = timestamps[indices]
        after_times = before_times + np.timedelta64(1,'h')
        time_diff = (target_times - before_times) / np.timedelta64(1, 's')
        if abs(time_diff) > 3600: # to see if there's some error of more than 1 hr between measurements
            return np.nan
        total_time_interval = (after_times - before_times) / np.timedelta64(1, 's')
        speed_diff = wind_speeds[indices + 1] - wind_speeds[indices]
        interpolated_speed = wind_speeds[indices] + (time_diff / total_time_interval) * speed_diff
        return interpolated_speed
    
    
matchups = np.empty((max_num_matchups,5))
counter = 0

for csv in csv_files_list: # go through each file and link up CYGNSS
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv)
    
    # Convert the date-time column to a pandas DateTimeIndex
    if data_folder == '\PSLGM':
        df['Datetime2'] = pd.to_datetime(df[' Date & UTC Time'])
    else:
        df['Datetime2'] = pd.to_datetime(df[['Year', 'Month', 'Day', 'Hour']])
    df['Datetime'] = df['Datetime2']
    df['Date_column'] = df['Datetime'].dt.date

    df.set_index('Datetime2', inplace=True)

    if data_folder == '\PSLGM':
        country = df.columns[-3]
        gauge_lat,gauge_lon = location_dict[country]['lat'],location_dict[country]['lon']
        wind_var = 'Wind Gust'  # 'Wind Speed' or 
        winds_array = np.array(df[wind_var])    

    elif data_folder == '\BOM Stations' or data_folder == '\Aus 1hr':
        country = [df['Latitude'][0],df['Longitude'][0]]
        gauge_lat,gauge_lon = float(country[0]),float(country[1])
        wind_var = 'Wind_speed_float'
        df[wind_var] = pd.to_numeric(df['Wind speed measured in m/s'], errors='coerce')
        winds_array = np.array(df[wind_var])    
        
    AOI_lat_min,AOI_lat_max = gauge_lat-search_distance,gauge_lat+search_distance
    AOI_lon_min,AOI_lon_max = gauge_lon-search_distance,gauge_lon+search_distance # bounding box around buoy location

    # Convert data to numpy arrays for vectorized operations
    times_array = np.array(df['Datetime'])  
    dates_array = np.array(df['Date_column'])
    AOI_time_start,AOI_time_end = np.min(times_array), np.max(times_array)

        
    # Generate a range of dates from start to end (exclusive)
    date_range = np.arange(AOI_time_start, AOI_time_end, dtype='datetime64[D]')

    # # Calculate the number of days between the dates
    # num_days = (date2 - date1).days
    # # Generate a list of dates between the start and end timestamps
    # date_list = [AOI_time_start + np.timedelta64(i,'D') for i in range(num_days + 1)]
    
    for date in date_range:
        # isolate the day's station values
        daily_bool = dates_array == date
            #date_measurement for date_measurement in times_array if np.datetime64(date_measurement, 'D') == date)
        daily_index = np.where(daily_bool)[0]
        if len(daily_index)>0:
            daily_times = times_array[daily_index]
            daily_winds = winds_array[daily_index]
            
            #find and open the corresponding cygnss file
            cyg_date_ind = np.where([date == x for x in cyg_file_dates])[0]
            try:
                # check if there is an available date 
                if len(cyg_date_ind) > 0:  
                    cyg_nc = nc.Dataset(cyg_files_list[cyg_date_ind[0]])
                else:
                    #print(date)
                    raise AssertionError
                
                # load in cyg variables
                cygnss_lons = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                cygnss_qual = cyg_nc.variables['sample_flags']
                cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
    
                if cyg_version== 'l2_3.2':
                    cygnss_YSLF_wind= cyg_nc.variables['preliminary_yslf_wind_speed'][:] ; 
                    inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['preliminary_yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
                    cygnss_qual = cyg_nc.variables['fds_sample_flags']

                #convert time to datetime format        
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = np.where(lats_boolean & lons_boolean)[0]
                if len(box_indices)>0:
                    date_format = "%Y-%m-%d %H:%M:%S"
                    cyg_start_date = np.datetime64(cyg_start_time_str)
                    datetime_array = np.vectorize(lambda x: cyg_start_date + np.timedelta64(int(x),'s'))(cygnss_time[box_indices])
                    
                    #time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end)
                    #true_indices = np.where(time_boolean)[0]
                    # if len(true_indices)>0:
                    cygnss_time_AOI = datetime_array ; 
                    
                    for i in range(len(cygnss_time_AOI)):                  
                        #interpolation of winds
                        
                        if type(cygnss_wind[box_indices][i]) is not np.ma.core.MaskedConstant and cygnss_qual[box_indices][i]%2 == 0:
                            station_wind_interp = linear_interpolation_vectorized(daily_times, daily_winds, cygnss_time_AOI[i])
                        
                            # now to check the values are valid then append to list
                            if station_wind_interp > 0: 
                                if cyg_version== 'l2_3.2' and type(cygnss_YSLF_wind[box_indices][i]) is not np.ma.core.MaskedConstant:
                                    matchups[counter] = [station_wind_interp,cygnss_wind[box_indices][i],cygnss_YSLF_wind[box_indices][i],cygnss_qual[box_indices][i],cygnss_time_AOI[i]] # variables to be saved
                                elif cyg_version== 'noaa_1.2':
                                    matchups[counter] = [station_wind_interp,cygnss_wind[box_indices][i],0,cygnss_qual[box_indices][i],cygnss_time_AOI[i]]
                                
                                counter+=1
                                if station_wind_interp > 40:
                                    print('weirdness')
                        else:
                            pass
            except AssertionError:
                pass

final_matchups = matchups[:counter]
if data_folder == '\PSLGM':
    area = 'Pac'
else:
    area = 'Aus'
    
# Save the list to a file
with open(cyg_folder+'\matchups_'+area+cyg_version+'.pkl', 'wb') as f:
    # filtered_array = np.array([[[x, y] for x, y in row if not np.isnan(x) and not np.isnan(y)] for row in matchups])
    pickle.dump(final_matchups, f, protocol=pickle.HIGHEST_PROTOCOL)

# Open the file and load the list
import pickle
with open(cyg_folder+'\matchups_'+area+cyg_version+'.pkl', 'rb') as f:
    filtered_array = pickle.load(f)

# Extract x and y values from the data
station_vals = np.asarray([item[0] for item in filtered_array])
cyg_FDS_vals = np.asarray([item[1] for item in filtered_array])
cyg_YSLF_vals = np.asarray([item[2] for item in filtered_array])
qual_flags = np.asarray([item[3] for item in filtered_array])

## flitering
# x_bool = (station_vals > 0)
# y_bool = (~np.isnan(cyg_vals))
# filtered_indices = np.where(x_bool & y_bool)[0]
# x_filter = station_vals[filtered_indices]
# y_filter = cyg_vals[filtered_indices]

# Plot the scatterplot/desnity plot for each value
# plt.scatter(station_vals, cyg_vals) # can be a box plot - want to see for 
plt.hist2d(station_vals, cyg_FDS_vals, bins=50, cmap='Blues', norm=matplotlib.colors.LogNorm())

plt.xlabel('Station (m/s)')
plt.ylabel('CYGNSS (m/s)')

cbar = plt.colorbar()
cbar.set_label('log(# meas.)')

# Calculate the minimum and maximum values of x and y
min_val = min(min(station_vals), min(cyg_FDS_vals))
max_val = max(max(station_vals), max(cyg_FDS_vals))

# Plot the dashed line along the 1:1 diagonal
plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='gray')

# # saving file
if cyg_version == 'l2_3.2':
    plt.title('Area:{0}, Alg:{1},FDS'.format(area,cyg_version))
    directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots\Station_CYG_{0}.{1},FDS.png'.format(cyg_version,area)
if cyg_version == 'noaa_1.2':
    plt.title('Area:{0}, Alg:{1}'.format(area,cyg_version))
    directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots\Station_CYG_{0}.{1}.png'.format(cyg_version,area)

plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)

if cyg_version == 'l2_3.2':
    # Plot the scatterplot for each value
    # plt.scatter(station_vals, cyg_vals) # can be a box plot - want to see for 
    plt.hist2d(station_vals, cyg_YSLF_vals, bins=50, cmap='Blues', norm=matplotlib.colors.LogNorm())
    
    plt.xlabel('Station (m/s)')
    plt.ylabel('CYGNSS (m/s)')
    plt.title('Area:{0}, Alg:{1},YSLF'.format(area,cyg_version))
    cbar = plt.colorbar()
    cbar.set_label('log(# meas.)')
    
    # Calculate the minimum and maximum values of x and y
    min_val = min(min(station_vals), min(cyg_YSLF_vals))
    max_val = max(max(station_vals), max(cyg_YSLF_vals))
    
    # Plot the dashed line along the 1:1 diagonal
    plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='gray')
    
    # # saving file
    directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots\Station_CYG_{0}.{1},YSLF.png'.format(cyg_version,area)
    plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
