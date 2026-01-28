# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:56:06 2024
Analyse buoy winds 
Input: buoy nc files
Output: 
@author: ashle
"""

import netCDF4 as nc
import numpy as np
from datetime import datetime
from datetime import timedelta 

comp = 'Ashley'
file =  r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\Buoy\w4n23w_hr.cdf'
folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\Buoy'

wind_threshold = 10
analysis_start_date = datetime(2021,1,1)
time_threshold = 2
search_area = 1 # how far in degrees to look around buoy

buoy_nc = nc.Dataset(file)
#time minutes since 2018-06-21 18:00:00
buoy_time = buoy_nc.variables['time'][:].astype('float'); buoy_speed = buoy_nc.variables['WS_401'][:,0,0,0]; buoy_lat = buoy_nc.variables['lat'][:]; buoy_lon = buoy_nc.variables['lon'][:]; buoy_depth = buoy_nc.variables['depth'][:]
if 'hr.cdf' in file:
    start_time_str = buoy_nc.variables['time'].units[12:]
elif '10m.cdf' in file:
    start_time_str = buoy_nc.variables['time'].units[14:]
# Define the format of the date and time in the string
date_format = "%Y-%m-%d %H:%M:%S"
# Convert the string to a datetime object
buoy_start_date = datetime.strptime(start_time_str, date_format)
# Calculate the timedelta between the two datetime objects
time_difference =  analysis_start_date - buoy_start_date
# Convert the time difference to minutes
if 'hr.cdf' in file:
    hrs_difference = time_difference.total_seconds()/60/60
    if hrs_difference >0: # need to shorten array to relevant times
        if buoy_time[-1] > hrs_difference: # need to make sure it has values in analysis
            ind= np.where(buoy_time>hrs_difference)[0][0]
            buoy_time_clipped=buoy_time[ind:]
            buoy_speed_clipped=buoy_speed[ind:]
elif '10m.cdf' in file:
    minutes_difference = time_difference.total_seconds() / 60 / 10 
    if minutes_difference >0: # need to shorten array to relevant times
        if buoy_time[-1] > minutes_difference: # need to make sure it has values in analysis
            ind= np.where(buoy_time>minutes_difference)[0][0]
            buoy_time_clipped=buoy_time[ind:]
            buoy_speed_clipped=buoy_speed[ind:]
 
# focus on high wind events after 2021
high_speed_inds = np.where(buoy_speed_clipped >wind_threshold)[0]
buoy_high_speed=buoy_speed_clipped[high_speed_inds]
buoy_high_time = buoy_time_clipped[high_speed_inds]
if 'hr.cdf' in file:
    buoy_datetime_array = np.vectorize(lambda x: buoy_start_date + timedelta(hours=x))(buoy_high_time)
elif '10m.cdf' in file:
    buoy_datetime_array = np.vectorize(lambda x: buoy_start_date + timedelta(minutes=x*10))(buoy_high_time)

# extract dates of high wind events
# Use a set to store unique days
unique_days = set()
# Iterate through the datetime objects and extract unique days
for dt in buoy_datetime_array:
    unique_days.add(dt.date())
# Convert the set to a list if needed
unique_days_list = list(unique_days)

# need to go through whole dataset to get values and save at netCDF

# # comparison to CYGNSS and SMAP
# for day in unique_days_list:
#     d = str(day.day); m=str(day.month); y = str(day.year)
#     if len(d) == 1 :
#         d= '0'+d
#     if len(m) == 1:
#         m = '0'+m
#     cyg_file_name = 'cyg.ddmi.s{}{}{}-000000-e{}{}{}-235959.l2.wind-mss.a31.d32.nc'.format(y,m,d,y,m,d)
#     SMAP_file_name = 'SMAP_L2B_SSS_46590_20231022T040843_R18290_V5.0'
    
#     # get the values of the particular days
#     buoy_day_ind = np.where(buoy_datetime_array == day)[0]
#     buoy_high_speed_day = buoy_high_speed[buoy_day_ind]
#     buoy_high_time_day = buoy_high_time[buoy_day_ind]
    
    
# for each day load in CYGNSS and SMAP
temp_cyg = r'C:\Users\\' +comp+ '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\L2_before\cyg.ddmi.s20231022-000000-e20231022-235959.l2.wind-mss.a31.d32.nc'
nc = nc.Dataset(temp_cyg)
cygnss_lons = nc.variables['lon'][:] ; cygnss_lats = nc.variables['lat'][:]; cygnss_wind= nc.variables['wind_speed'][:] ; cygnss_YSLF_wind= nc.variables['yslf_wind_speed'][:] ; cygnss_time =nc.variables['sample_time'][:]
inc = nc.variables['incidence_angle']; rcg = nc.variables['range_corr_gain'][:]; uncertainty= nc.variables['yslf_wind_speed_uncertainty']; sc_num=nc.variables['spacecraft_num']
#convert time to datetime format
cyg_start_time_str = nc.variables['sample_time'].units[14:33]
# Define the format of the date and time in the string
date_format = "%Y-%m-%d %H:%M:%S"
cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
# # Calculate distances
# points = np.zeros((2,len(cygnss_lats)))
# points[0] = cygnss_lats; points[1] = cygnss_lons
# points = points.T
# see how many observation in radius within x hours
# def euclidean_distance(array, point):
#     # Calculate Euclidean distance between two points
#     return np.sqrt(np.sum((point1 - point2)**2))
# distances = np.array([euclidean_distance(target_point, p) for p in points]) # takes a long time
# Select points within a 1-degree radius or box
buoy_location = [buoy_lat,buoy_lon]  # Latitude and longitude of San Francisco, CA
buoy_area = [buoy_location[0]-search_area,buoy_location[1]-search_area,
             buoy_location[0]+search_area,buoy_location[1]+search_area] # bounding box around buoy location

ind_in_radius = np.where((cygnss_lats > buoy_area[0]) & (cygnss_lons > buoy_area[1]) & (cygnss_lats < buoy_area[2]) & (cygnss_lons < buoy_area[3]) )
# ind_in_radius = np.where((cygnss_lats > buoy_area[0]))

cyg_times_in_radius = cygnss_time[ind_in_radius]
#now with the set of times, subset to be within range
cyg_datetime_array_day = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cyg_times_in_radius)
cyg_YSLF_in_radius = cygnss_YSLF_wind[ind_in_radius]
# need to filter by time
buoy_day_ind = np.where([x.date() == datetime(2023,1,20).date() for x in buoy_datetime_array])
buoy_datetime_array_day = buoy_datetime_array[buoy_day_ind]
buoy_high_speed_day = buoy_high_speed[buoy_day_ind]

# ind_in_time_radius = ind_in_radius[np.where(abs(cyg_datetime_array-buoy_datetime_array_day)<time_threshold)]
# points_in = points[ind_in_time_radius]
# wind_in = cygnss_YSLF_wind[ind_in_time_radius]

# save matching values 
# combine all different matchups 
# stats overall and by station

# plotting
import matplotlib.pyplot as plt

# plot timeseries
fig, ax = plt.subplots()

ax.scatter(buoy_datetime_array_day, buoy_high_speed_day, 'o--')
ax.scatter(cyg_datetime_array_day, cyg_YSLF_in_radius, marker = 'o')

ax.tick_params(rotation=45)

plt.show()

