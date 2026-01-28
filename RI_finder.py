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
storm_ids=[]
storm_names=[]
inds_before=[]
inds_after=[]
times=[]
lats=[]
lons=[]
winds=[]
ri_or_rd=[]
duration=[]
intensity=[]

path = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS'
file_name = '\IBTrACS.since1980.v04r01.nc' #'\IBTrACS.last3years.v04r00.nc'
IB_nc = nc.Dataset(path + file_name)
WMO_wind = IB_nc.variables['wmo_wind'][:]
USA_wind = IB_nc.variables['usa_wind']
wmo_agency = IB_nc.variables['wmo_agency']
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

# find locations of RI
from datetime import datetime, timedelta
dt= dt_concatenated.astype(datetime)

for time_threshold in time_thresholds:
    for tc in range(rows):
        start_ind = 0 ; end_ind = 1
        counter = 0
        year = np.datetime64('1970-01-01T00:00:00') + np.timedelta64(int(dt[tc,end_ind]), 's')
        if np.datetime64(year,'Y') > np.datetime64('2018') and np.datetime64(year,'Y') < np.datetime64('2025'): # only run for Valid CYGNSS years
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
                    # USA_wind_dif = USA_wind[tc,end_ind] - USA_wind[tc,start_ind]
                    for wind_threshold in wind_thresholds:
                        if WMO_wind_dif >= wind_threshold: # or USA_wind_dif >= wind_threshold: # may not be comparable
                            if storm_id_concatenated[tc] == '2022239N22150':
                                pass
                            storm_ids.append(storm_id_concatenated[tc])  
                            storm_names.append(storm_name[tc])
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
                            # if abs(USA_wind_dif) >= abs(WMO_wind_dif) or type(WMO_wind_dif) != np.int16:
                            #     winds.append(USA_wind[tc,start_ind_before:end_ind_after])
                            # else:
                            winds.append(WMO_wind[tc,start_ind_before:end_ind_after])
                            
                            if WMO_wind_dif > 0: #or USA_wind_dif > 0:
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
RI_file = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v5.nc'
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

