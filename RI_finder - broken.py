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
comp = 'ashle'
time_thresholds = [6,12,24]
wind_thresholds = [45,30] # units of knots
max_len = 21
storm_ids=np.array([],dtype=str)
storm_names=np.array([],dtype=str)
times=np.empty((0, max_len))  # 0 rows, 2 columns
lats=np.empty((0, max_len))
lons=np.empty((0, max_len))
winds=np.empty((0, max_len))
ri_or_rd=[]
duration=np.array([])
intensity=np.array([])

path = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS'
file_name = '\IBTrACS.since1980.v04r01.nc' #'\IBTrACS.last3years.v04r00.nc'
IB_nc = nc.Dataset(path + file_name)
WMO_wind = IB_nc.variables['wmo_wind'][:]
USA_wind = IB_nc.variables['usa_wind']
storm_id = IB_nc.variables['sid'][:].astype(str)
storm_name = IB_nc.variables['name'][:].astype(str)
dt_str = IB_nc.variables['iso_time'][:].astype('U2')
lat = IB_nc.variables['lat'][:]
lon = IB_nc.variables['lon'][:]

rows, cols, depth = dt_str.shape
dt_concatenated = np.empty((rows, cols))
storm_id_concatenated = np.zeros((rows)).astype(str) 

# need ot fix here for name


for i in range(rows):
    storm_id_concatenated[i] = np.array([''.join(storm_id[i,:])])[0]
    
    for j in range(cols):
        if not type(dt_str[i,j,0]) == np.ma.core.MaskedConstant:
            dt_concatenated[i,j]= np.array([''.join(dt_str[i,j,:])],dtype = 'datetime64[s]')[0]

def padding(data,max_len):
    # input data list of numpy arrays and adds numpy nan to the end so all the arrays are same length
    padded_data = np.full(max_len, np.nan)  # Create 2D array filled with NaN
    padded_data[0:len(data)] = data  # Fill in values from each array        
    return padded_data

# find locations of RI
from datetime import datetime
dt= dt_concatenated.astype(datetime)
for time_threshold in time_thresholds:
    for wind_threshold in wind_thresholds:
        for tc in range(rows):
            start_ind = 0 ; end_ind = 1
            counter = 0
            year = np.datetime64('1970-01-01T00:00:00') + np.timedelta64(int(dt[tc,end_ind]), 's')
            if np.datetime64(year,'Y') > np.datetime64('2018'): # only run for Valid CYGNSS years
                while dt[tc,end_ind] > 1 and counter <cols*2: 
                    counter+=1
                    time_difference = (dt[tc,end_ind] - dt[tc,start_ind])/60/60 # convert from seconds to hours
                    if time_difference < time_threshold:
                        end_ind += 1
                        if end_ind >= cols: # if it makes it to the end of array
                            break 
                    elif time_difference > time_threshold:
                        start_ind =+ 1
            
                    elif time_difference == time_threshold:
                        valid = True
                        # if (storm_id_concatenated[tc] not in storm_ids):
                        #     valid = True
                        # elif wind_threshold not in np.atleast_1d(intensity[storm_ids==storm_id_concatenated[tc]]):
                        #     valid = True    
                        if valid:
                            # now to compare speeds to see if threshold is met 
                            WMO_wind_dif = WMO_wind[tc,end_ind] - WMO_wind[tc,start_ind]
                            USA_wind_dif = USA_wind[tc,end_ind] - USA_wind[tc,start_ind]
                            if WMO_wind_dif >= wind_threshold or USA_wind_dif >= wind_threshold: # may not be comparable
                                storm_ids=np.concatenate((storm_ids,[storm_id_concatenated[tc]]))  
                                # storm_names= np.concatenate((storm_names,storm_name[tc]))
                                end_ind += 1 # to make sure we include the end measurement
                                # Pad arrays to the same length with NaN
                                times_padded = padding(dt[tc,start_ind:end_ind],max_len)  
                                lats_padded = padding(lat[tc,start_ind:end_ind],max_len)  
                                lons_padded = padding(lon[tc,start_ind:end_ind],max_len)  
                                times=np.concatenate((times,[times_padded]))
                                lats=np.concatenate((lats,[lats_padded]))
                                lons=np.concatenate((lons,[lons_padded]))
                                duration=np.concatenate((duration,[time_difference]))
                                intensity=np.concatenate((intensity,[wind_threshold]))
                                if abs(USA_wind_dif) >= abs(WMO_wind_dif) or type(WMO_wind_dif) != np.int16:
                                    winds_padded = padding(USA_wind[tc,start_ind:end_ind],max_len)
                                    winds=np.concatenate((winds,[winds_padded]))
                                else:
                                    winds_padded = padding(WMO_wind[tc,start_ind:end_ind],max_len)
                                    winds=np.concatenate((winds,[winds_padded]))
                                
                                if WMO_wind_dif > 0 or USA_wind_dif > 0:
                                    ri_or_rd.append(1)
                                else:
                                    ri_or_rd.append(-1)
                                # break # to ensure that there is only 1 event saved for each interval
                                # if end_ind-start_ind >9:
                                #     print((tc,start_ind,end_ind))
                            start_ind += 1
                            # end_ind += 1
                            
                    else:
                        start_ind += 1
    
# saving the details      
        
import netCDF4 as nc
RI_file = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v3.nc'
with nc.Dataset(RI_file, 'w') as ncfile:
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
    vmax.missing_value = -9999.0
    ri_rd = ncfile.createVariable('ri_rd', 'i4', ('events',))
    durations = ncfile.createVariable('durations', 'i4', ('events',))
    intensities = ncfile.createVariable('intensities', 'i4', ('events',))


    # Add data to variables
    time[:] = times
    storm[:] = storm_ids
    latitude[:] = lats
    longitude[:] = lons
    vmax[:] = winds
    ri_rd[:] = ri_or_rd
    durations[:]=duration
    intensities[:]=intensity
    
