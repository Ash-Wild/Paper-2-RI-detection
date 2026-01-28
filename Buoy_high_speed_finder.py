# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:21 2024
Purpose: find the buoys with wind speed above a threshold
input: folder of buoy.cdf files, and wind threshold
return: list of dates with wind above threshold after particular date
@author: Ashley
"""

import netCDF4 as nc
import numpy as np
from datetime import datetime
from datetime import timedelta 
import glob

comp = 'Ashley'
# file =  r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\Buoy\w4n23w_hr.cdf'
folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\Buoy'

wind_threshold = 10
analysis_start_date = datetime(2021,1,1)
valid_buoy_files = []
valid_days = []

files_list = glob.glob(folder+ "\*.cdf")
for file in files_list: 
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
                buoy_time=buoy_time[ind:]
                buoy_speed=buoy_speed[ind:]
    elif '10m.cdf' in file:
        minutes_difference = time_difference.total_seconds() / 60 / 10
        if minutes_difference >0: # need to shorten array to relevant times
            if buoy_time[-1] > minutes_difference: # need to make sure it has values in analysis
                ind= np.where(buoy_time>minutes_difference)[0][0]
                buoy_time=buoy_time[ind:]
                buoy_speed=buoy_speed[ind:]
    
    # create list 
    if np.max(buoy_speed) > wind_threshold:
        valid_buoy_files.append(file)
        