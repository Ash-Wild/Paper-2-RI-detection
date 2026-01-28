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
import pickle
import netCDF4 as nc
import matplotlib.pyplot as plt
from bisect import bisect_left
import pandas as pd
# from BOM_high_speed_finder import wind_thresh, time_thresh

# input files from PSLGM_high_speed_finder
comp = 'Ashley'
cyg_version = 'l2_3.2' # 'noaa_1.2' OR 'l2_3.2' or 'storm_centric'
data_folder = '\Aus 1hr' # data_folder = '\BOM Stations' OR ' PSLGM' OR 'Aus 1hr'
csv_folder  = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\BOM stations' + data_folder
search_distance = 0.5 # number of gedrees for searching for CYGNSS
cyg_time_thresh = 1 # number of hours to look for CYGNSS outside of event 
cyg_AOI = False # will switch to true if there's some overlap

# Open the file and load the list
with open(csv_folder+'\events_list.pkl', 'rb') as f:
    events_list = pickle.load(f)

# interpolation
# make a numpy array based on number of events, 

# For each event, open relevant CYGNSS file and clip to .5 deg and time
# location_dict = {'Cook Islands':[21° 11' 58" S    159° 47' 10" W],
# 'Federated States of Micronesia':[06° 58' 42" N    158° 11' 50" E],
# 'Fiji':[17° 36' 19" S    177° 26' 17" E],
# 'Kiribati' : [1° 21' 45" N    172° 55' 48" E],  
# 'Marshall Islands': [7° 06' 27" N    171° 22' 15" E],
# 'Nauru':[0° 31' 55" S    166° 54' 33" E],
# 'Papua New Guinea':[2° 02' 10" S    147° 22' 31" E],
# 'Samoa':[13° 49'09" S    171° 45' 21" W],
# 'Solomon Islands':[9° 25' 18" S    159° 57' 19" E],
# 'Tonga':[21° 08' 25" S    175° 10' 45" W],
# 'Tuvalu':[08° 30' 10" S    179° 12' 33" E],
# 'Vanuatu':[17° 45' 41" S    168° 17' 35" E],
# 'Niue':[19° 03' 10" S    169° 55' 15" W]
#     }
location_dict = {
    'Cook Islands': {'lat': -21.199, 'lon': -159.786},
    'Federated States of Micronesia': {'lat': 6.978, 'lon': 158.197},
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
def linear_interpolation_vectorized(timestamps, wind_speeds, target_times):
    
    indices = np.searchsorted(timestamps, target_times) - 1
    # index = bisect_left(np.flip(timestamps), target_times)

    if indices == -1: # if it is before the first speed, take the first speed
        return wind_speeds[0]
    if indices == len(timestamps)-1: # if it is after the last reading, take the last speed 
        return wind_speeds[-1]
    else: # for all else inbetween 
        before_times = timestamps[indices]
        if indices == len(timestamps)-1:
            after_times =  before_times + timedelta(hours=1)
        else:
            after_times =  timestamps[indices+1]
            
        time_diff = min((target_times - before_times) / np.timedelta64(1, 's'),(target_times - after_times) / np.timedelta64(1, 's'))
        if abs(time_diff) > 3600: # to see if there's some error of more than 1 hr between measurements
            print('something wrong')
        total_time_interval = (after_times - before_times) / np.timedelta64(1, 's')
        speed_diff = wind_speeds[indices + 1] - wind_speeds[indices]
        interpolated_speed = wind_speeds[indices] + (time_diff / total_time_interval) * speed_diff
        return interpolated_speed

def get_matchups(cyg_files_list,start_time,end_time,matchups_list,times_array, winds_array,AOI):
    # Extract date components from the timestamps
    date_start_storm = start_time.date()
    date_end_storm = end_time.date()

    # Calculate the number of days between the dates
    num_days_storm = (date_end_storm - date_start_storm).days


    # Generate a list of dates between the start and end timestamps
    date_list_storm = [date_start_storm + timedelta(days=i) for i in range(num_days_storm + 1)]

    AOI_lat_min,AOI_lat_max,AOI_lon_min,AOI_lon_max = AOI

    # cyg_lats=np.ma.MaskedArray([]); cyg_lons=np.ma.MaskedArray([]); cyg_winds=np.ma.MaskedArray([]); cyg_times=np.ma.MaskedArray([])

    for date in date_list_storm:
        #find and open the corresponding cygnss file
        try:
            assert date in cyg_file_dates
            if cyg_version == 'storm_centric':
                cyg_file_path = cyg_files_list[cyg_files_num[np.where([date == x for x in cyg_file_dates])[0][0]]]
                cyg_nc = nc.Dataset(cyg_file_path)
                cyg_storm_lat = cyg_nc.variables['best_track_storm_center_lat'][:] ; cyg_storm_lon = ((cyg_nc.variables['best_track_storm_center_lon'][:] + 180) % 360) - 180
                # assert True in (cyg_storm_lon > AOI_lon_min) & True in (cyg_storm_lon < AOI_lon_max)
                # assert True in (cyg_storm_lat > AOI_lat_min) & True in (cyg_storm_lat < AOI_lat_max)
    
                longitudes  = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                cygnss_lons = ((longitudes + 180) % 360) - 180
                cygnss_time =cyg_nc.variables['time'][:] ; cygnss_time_offset = cyg_nc.variables['time_offset'][:]
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = [0,0]
                
            else:
                cyg_file_path = cyg_files_list[np.where([date == x for x in cyg_file_dates])[0][0]]
                cyg_nc = nc.Dataset(cyg_file_path)
                longitudes = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                cygnss_lons = ((longitudes + 180) % 360) - 180
                cygnss_qual = cyg_nc.variables['sample_flags']
                cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = np.where(lats_boolean & lons_boolean)[0]
                
            # if cyg_version== 'l2_3.2':
            #     cygnss_yslf_wind= cyg_nc.variables['preliminary_yslf_wind_speed'][:];
            #     merge_threshold = 17
            #     cygnss_wind = np.where(cygnss_wind < merge_threshold, cygnss_wind, cygnss_yslf_wind)
            #     inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['preliminary_yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
            #     # cygnss_qual = cyg_nc.variables['fds_sample_flags']
            #convert time to datetime format
            
    
            if len(box_indices)>0:
                if cyg_version == 'storm_centric':
                    # find closest epoch to SAR time
                    date_format = "%y-%m-%d %H:%M:%S"
                    cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
                    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                    cyg_measured_times = np.ma.MaskedArray(np.zeros(np.shape(cygnss_time_offset)))
                    for epoch in range(len(cygnss_time)):
                        cyg_measured_times[epoch] = cygnss_time_offset[epoch] + cygnss_time[epoch]
                    # Focus on the lat and lons
                    try:
                        timedeltas = np.vectorize(lambda x: timedelta(hours=x))(cyg_measured_times[:,lats_boolean,:][:,:,lons_boolean])
                    except ValueError:
                        raise AssertionError
                    # Add the base datetime to each timedelta to get the 3D array of datetime objects
                    cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + x)(timedeltas)                        
                    
                else:
                    cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
                    date_format = "%Y-%m-%d %H:%M:%S"
                    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                    cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
                
                # time_boolean = (cyg_datetime_array > cyg_AOI_time_start) & (cyg_datetime_array < cyg_AOI_time_end)
                # true_indices = np.where(time_boolean)[0]
                # if len(true_indices)>0:
                if cyg_version == 'storm_centric':
                    cygnss_wind_AOI = cygnss_wind[:,lats_boolean,:][:,:,lons_boolean] ; cygnss_time_AOI = cyg_datetime_array
                    cygnss_lats_AOI = cygnss_lats[lats_boolean] ; cygnss_lons_AOI = cygnss_lons[lons_boolean]
                else:
                    cygnss_wind_AOI = cygnss_wind[box_indices] ; cygnss_time_AOI = cyg_datetime_array ; 
                    #   sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
                    cygnss_lats_AOI = cygnss_lats[box_indices] ; cygnss_lons_AOI = cygnss_lons[box_indices]
                    cyg_qual_AOI = cygnss_qual[box_indices]
                cyg_AOI = True 
                
                # cyg_lats=np.append(cyg_lats,cygnss_lats_AOI) ; cyg_lons = np.append(cyg_lons,cygnss_lons_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI); cyg_times=np.append(cyg_times,cygnss_time_AOI)
                
                for i in range(len(cygnss_time_AOI)):
                    # # Calculate the differences between the timestamp and each time in the list
                    # differences = [abs((time - cygnss_time_AOI[i]).total_seconds()) for time in times_array]
                    # # Find the index of the time with the minimum difference
                    # closest_time_index = differences.index(min(differences))
                    
                    #interpolation
                    if cyg_version == 'storm_centric':
                        for j in range(len(cygnss_time_AOI[0,:,0])):
                            for k in range(len(cygnss_time_AOI[0,0,:])):
                                if isinstance(cygnss_time_AOI[i,j,k], (datetime, np.datetime64)):
                                    if cygnss_time_AOI[i,j,k] > start_time and cygnss_time_AOI[i,j,k] < end_time:
                                        station_wind_interp = linear_interpolation_vectorized(times_array, winds_array, pd.Timestamp(cygnss_time_AOI[i,j,k]))
                                        matchups_list.append([station_wind_interp, cygnss_wind_AOI[i,j,k]])
                    else:
                        if cygnss_time_AOI[i] > start_time and cygnss_time_AOI[i] < end_time:
                            station_wind_interp = linear_interpolation_vectorized(times_array, winds_array, pd.Timestamp(cygnss_time_AOI[i]))
                            matchups_list.append([station_wind_interp, cygnss_wind_AOI[i]])
                    # conditions_all.append([inc_AOI[i], rcg_AOI[i], uncertainty_AOI[i],sc_num_AOI[i]])
                    # if station_wind_interp > 30:
                    #     continue
                    # if station_wind_interp - cygnss_wind_AOI[i] < -20:
                    #     conditions_all.pop()
                        # conditions_biased.append([inc_AOI[i], rcg_AOI[i], uncertainty_AOI[i],sc_num_AOI[i]])
        except AssertionError:
            pass
                        
                        
# maybe I need ot interpolate wind data to minutes to match up?
        
cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\cyg_'+cyg_version
cyg_files_list = glob.glob(cyg_folder+ "\*.nc")

# find the dates of cyg 
cyg_file_dates = []
cyg_storm_names = []
cyg_files_num = []
file_counter = 0 
for cyg_file in cyg_files_list:
    start_ind = cyg_file.find('ddmi.')
    if cyg_version == 'storm_centric':
        end_ind = cyg_file.find('.202')
        storm_name = cyg_file[start_ind:end_ind]
        cyg_nc = nc.Dataset(cyg_file)
        date_format = "%y-%m-%d %H:%M:%S"
        cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
        cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
        cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(hours=int(x)))(cyg_nc.variables['time'][:])
        for i in cyg_datetime_array:
            cyg_file_dates.append(i.date())
            cyg_files_num.append(file_counter)
            cyg_storm_names.append(storm_name)
        file_counter += 1
    else:
        cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
        cyg_file_dates.append(cyg_datetime_object.date())

num_measured_events = 0
matchups_storm = []
matchups_before = []
matchups_after = []

conditions_all = []
conditions_biased = []

for country,dates,wind_gust_before,wind_gust_after in events_list:
    if data_folder == '\PSLGM':
        gauge_lat,gauge_lon = location_dict[country]['lat'],location_dict[country]['lon']
    elif data_folder == '\BOM Stations' or data_folder == '\Aus 1hr':
        gauge_lat,gauge_lon = float(country[0]),float(country[1])
        
    AOI_lat_min,AOI_lat_max = gauge_lat-search_distance,gauge_lat+search_distance
    AOI_lon_min,AOI_lon_max = gauge_lon-search_distance,gauge_lon+search_distance # bounding box around buoy location
    AOI = (AOI_lat_min,AOI_lat_max,AOI_lon_min,AOI_lon_max )
    AOI_time_start_storm,AOI_time_end_storm = dates[0][0]-timedelta(hours=1), dates[-1][0]+timedelta(hours=1)
    AOI_time_start_before,AOI_time_end_before = wind_gust_before[-1][0]-timedelta(hours=cyg_time_thresh), wind_gust_before[0][0]
    AOI_time_start_after,AOI_time_end_after = wind_gust_after[0][0], wind_gust_after[-1][0]+timedelta(hours=cyg_time_thresh)

    # Convert data to numpy arrays for vectorized operations
    timestamps_storm, wind_speeds_storm = zip(*dates)
    times_array_storm = np.array(timestamps_storm)
    winds_array_storm = np.array(wind_speeds_storm)    
    
    timestamps_before, wind_speeds_before = zip(*wind_gust_before)
    times_array_before = np.flip(np.array(timestamps_before))
    winds_array_before = np.flip(np.array(wind_speeds_before))
    
    timestamps_after, wind_speeds_after = zip(*wind_gust_after)
    times_array_after = np.array(timestamps_after)
    winds_array_after = np.array(wind_speeds_after)
    
    get_matchups(cyg_files_list, AOI_time_start_storm, AOI_time_end_storm,matchups_storm,times_array_storm,winds_array_storm,AOI)
    get_matchups(cyg_files_list, AOI_time_start_before, AOI_time_end_before,matchups_before,times_array_before,winds_array_before,AOI)
    get_matchups(cyg_files_list, AOI_time_start_after, AOI_time_end_after,matchups_after,times_array_after,winds_array_after,AOI)




# Extract x and y values from the data
station_vals_storm = [item[0] for item in matchups_storm]
Wind_Vals_Storm = [item[1] for item in matchups_storm]

station_vals_before = [item[0] for item in matchups_before]
Wind_Vals_before = [item[1] for item in matchups_before]

station_vals_after = [item[0] for item in matchups_after]
Wind_Vals_after = [item[1] for item in matchups_after]

station_all = np.array(station_vals_before+station_vals_storm+station_vals_after)
wind_all = np.array(Wind_Vals_before+Wind_Vals_Storm+Wind_Vals_after)

mask = ~np.isnan(station_all) & ~np.isnan(wind_all)
station_all_masked = station_all[mask]
wind_all_masked = wind_all[mask]

R=np.around(np.corrcoef(station_all_masked,wind_all_masked)[0,1], decimals=2)
RMSE = np.around(np.sqrt(np.mean((wind_all_masked - station_all_masked) ** 2)), decimals=2)
bias = np.around(np.mean(wind_all_masked - station_all_masked),decimals=2)

# Plot the scatterplot for each value
# import matplotlib.pyplot as plt
# import matplotlib.patches as ptc
# import numpy as np
# from shapely.geometry import Point
# from shapely.ops import cascaded_union

# size = 0.02
# alpha = 0.5

# polygons1 = [Point(station_vals_storm[i], YSLF_vals_storm[i]).buffer(size) for i in range(len(matchups_storm))]
# polygons2 = [Point(station_vals_before[i], YSLF_vals_before[i]).buffer(size) for i in range(len(matchups_before))]
# polygons3 = [Point(station_vals_after[i], YSLF_vals_after[i]).buffer(size) for i in range(len(matchups_after))]

# polygons1 = cascaded_union(polygons1)
# polygons2 = cascaded_union(polygons2)
# polygons3 = cascaded_union(polygons3)


# fig = plt.figure(figsize=(4,4))
# ax = fig.add_subplot(111, title="Test scatter")
# for polygon1 in polygons1:
#     polygon1 = ptc.Polygon(np.array(polygon1.exterior), facecolor="blue", lw=0, alpha=alpha)
#     ax.add_patch(polygon1)
# for polygon2 in polygons2:
#     polygon2 = ptc.Polygon(np.array(polygon2.exterior), facecolor="grey", lw=0, alpha=alpha)
#     ax.add_patch(polygon2)
# for polygon3 in polygons3:
#     polygon3 = ptc.Polygon(np.array(polygon3.exterior), facecolor="red", lw=0, alpha=alpha)
#     ax.add_patch(polygon2)
# ax.axis([-0.2, 1.2, -0.2, 1.2])
marker = 'o'
size = 16
min_speed = 0
max_speed = 31
plt.figure(figsize=(5,5))
plt.scatter(station_vals_before, Wind_Vals_before, edgecolor='grey', facecolors='none', label='Before', marker = marker, s = size) 
plt.scatter(station_vals_after, Wind_Vals_after, edgecolor='red', facecolors='none', label = 'After', marker = marker, s = size) 
plt.scatter(station_vals_storm, Wind_Vals_Storm, edgecolor='blue', facecolors='none', label = 'Storm', marker = marker, s = size) # can be a box plot - want to see for 
plt.xlabel('Station (m/s)')
plt.ylabel('CYGNSS (m/s)')
plt.xlim(min_speed, max_speed)  # Change these values to your desired x bounds
plt.ylim(min_speed, max_speed)  # Change these values to your desired y bounds

plt.legend(loc='upper right')

if data_folder == '\PSLGM':
    area = 'Pac'
else:
    area = 'Aus'
    
# plt.title('Area:{0}, Alg: {1}, Speed:{2} Hrs:{3}'.format(area, cyg_version, wind_thresh,time_thresh))
plt.text(0.05, 0.95, f'R: {R}', transform=plt.gca().transAxes, va='top', ha='left')
plt.text(0.05, 0.90, f'RMSE: {RMSE}', transform=plt.gca().transAxes, va='top', ha='left')
plt.text(0.05, 0.85, f'Bias: {bias}', transform=plt.gca().transAxes, va='top', ha='left')
plt.text(0.05, 0.80, f'# Cyg Meas.: {len(station_vals_storm)}', transform=plt.gca().transAxes, va='top', ha='left')

# plt.text(0.05, 0.75, f'# Measured Events: {num_measured_events}', transform=plt.gca().transAxes, va='top', ha='left')



# Calculate the minimum and maximum values of x and y
min_val = min(min(station_vals_before), min(station_vals_after),min(Wind_Vals_before), min(Wind_Vals_after))
max_val = max(max(station_vals_storm), max(Wind_Vals_Storm),max(station_vals_after), max(Wind_Vals_after),max(station_vals_before), max(Wind_Vals_before))

# Plot the dashed line along the 1:1 diagonal
plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='gray')

# # saving file
directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots\Station_CYG_{0}.{1}W17T3.png'.format(cyg_version,area)
plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)



# # analysis results
# # Convert the list of lists to NumPy arrays
# data_array1 = np.array(conditions_biased)
# data_array2 = np.array(conditions_all)

# # Extract columns from both datasets
# data_columns1 = [data_array1[:, i] for i in range(data_array1.shape[1])]
# data_columns2 = [data_array2[:, i] for i in range(data_array2.shape[1])]
# focus_variable = 3
# if focus_variable>-1:
#     col4_data1  = data_array1[:, focus_variable]
#     col4_data2  = data_array2[:, focus_variable]
    
#     # Create a list to hold the data for boxplot
#     data_columns = [col4_data1, col4_data2]
    
#     # Define boxplot positions
#     positions = [1, 2]
    
#     # Create box and whisker plots
#     plt.figure(figsize=(8, 6))
#     box1 = plt.boxplot([col4_data1], positions=[1], patch_artist=True, widths=0.6)
#     box2 = plt.boxplot([col4_data2], positions=[2], patch_artist=True, widths=0.6)

# else:
#     # Create a list to hold all boxplot data for combined plotting
#     combined_data = []
    
#     # Combine the data for each column side by side
#     for col1, col2 in zip(data_columns1, data_columns2):
#         combined_data.append(col1)
#         combined_data.append(col2)
    
#     # Define boxplot positions
#     positions = []
#     for i in range(1, len(data_columns1) + 1):
#         positions.append(i * 2 - 1)  # Position for dataset 1
#         positions.append(i * 2)      # Position for dataset 2
    
#     # Create box and whisker plots
#     plt.figure(figsize=(12, 6))
#     box1 = plt.boxplot(data_columns1, positions=positions[::2], patch_artist=True, widths=0.6)
#     box2 = plt.boxplot(data_columns2, positions=positions[1::2], patch_artist=True, widths=0.6)

# # Color the boxplots
# colors1 = 'blue'
# colors2 = 'red'

# for box in box1['boxes']:
#     box.set(facecolor=colors1)

# for box in box2['boxes']:
#     box.set(facecolor=colors2)

# # Customize the plot
# plt.title("Comparison of Two Datasets with Box and Whisker Plots")
# plt.xlabel("Column Number")
# plt.ylabel("Value")
# plt.xticks(range(1, len(data_columns1) * 2 + 1, 2), ['Column 1', 'Column 2', 'Column 3', 'Column 4'])

# # Add a legend
# plt.legend([box1["boxes"][0], box2["boxes"][0]], ['Dataset 1', 'Dataset 2'], loc='upper right')

# # Show the plot
# plt.show()