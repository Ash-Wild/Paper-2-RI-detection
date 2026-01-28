# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 15:25:40 2023

Purpose: to plot multiple values as a timeseries at particular location

@author: Ashley
"""
import netCDF4
import datetime
import numpy as np


dir = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Data\Oct 17 2023 comparison'
# SMAP
SMAP = netCDF4.Dataset(dir + '\RSS_smap_wind_daily_2023_10_17_NRT_v01.0.nc')
SMAP_wind, SMAP_lons, SMAP_lats, SMAP_time = SMAP.variables['wind'][:], SMAP.variables['lon'][:], SMAP.variables['lat'][:], SMAP.variables['minute'][:]
# #try to update SMAP mask - slow and SMAP has strange values
# new_mask = np.ones(np.shape(SMAP_wind))
# for orbit in range(np.shape(SMAP_wind)[2]):
#     for i in range(np.shape(SMAP_wind)[0]):
#         for j in range(np.shape(SMAP_wind)[1]):
#             new_mask[i,j,orbit]=new_mask[i,j,orbit]-(0 < SMAP_wind.data[i,j,orbit] <60)
# new_mask= new_mask.astype(dtype=bool)
# SMAP_wind = np.ma.masked_array(SMAP_wind.data,mask=new_mask)            
SMAP_wind_asc, SMAP_wind_des = SMAP_wind[:,:,0], SMAP_wind[:,:,1]
SMAP_time_asc, SMAP_time_des = SMAP_time[:,:,0], SMAP_time[:,:,1]

# buoy 
buoy_location = [2,345] # deg north and east
search_area = 1 # how far in degrees to look aroudn buoy
buoy_area = [buoy_location[0]-search_area,buoy_location[1]-search_area,
             buoy_location[0]+search_area,buoy_location[1]+search_area] # bounding box aroudn buoy location


import pandas as pd
buoy_data = pd.read_excel(dir + '\Buoy NOAA.xlsx')
times = buoy_data.columns
timeseries = []

for hour in range(len(buoy_data)):
    for time in times[1:7]: 
        timeseries.append(float(buoy_data[time][hour][0:4]))

timeseries_correct = timeseries.reverse() # error here

start_time = datetime.datetime(2023,10,17,0,10)
end_time = datetime.datetime(2023,10,17,23,0)
buoy_times = pd.date_range(start_time,end_time,freq='10t')



#cygnss
cygnss = netCDF4.Dataset(dir + '\cyg.ddmi.s20231017-003000-e20231017-233000.l3.grid-wind.a31.d32.nc')
cygnss_lons = cygnss.variables['lon'][:] ; cygnss_lats = cygnss.variables['lat'][:]; cygnss_wind= cygnss.variables['wind_speed'][:] ; cygnss_time =cygnss.variables['time'][:]


#ascat
import glob
path = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Data\Oct 17 2023 comparison\ASCAT\\'
file_name = 'ascat_20231017_055400_metopc_25647_eps_o_250_3301_ovw.l2.nc'
files_list = glob.glob(path+ "\*.nc")
ascat = netCDF4.Dataset(path + file_name)
ascat_lons = ascat.variables['lon'][:] ; ascat_lats = ascat.variables['lat'][:]; ascat_wind= ascat.variables['wind_speed'][:] ; ascat_time =ascat.variables['time'][:]


# clipping to area around buoy
SMAP_coverage = [-90,0,90,360] #change this depending on dataset, global would be -90,-180 (or 0), 90, 180 (or 360).
cygnss_coverage = [-40,0,40,360]
SMAP_resolution = 0.25
cygnss_resolution = 0.2

SMAP_difference = np.subtract(buoy_area,SMAP_coverage)
cygnss_difference = np.subtract(buoy_area,cygnss_coverage)

SMAP_difference_indexes = np.divide(SMAP_difference,SMAP_resolution)
cygnss_difference_indexes = np.divide(cygnss_difference,cygnss_resolution)

SMAP_subset_float_indexes = [SMAP_difference_indexes[0],SMAP_difference_indexes[1],
                             SMAP_difference_indexes[2]+SMAP_lats.shape[0],SMAP_difference_indexes[3]+SMAP_lons.shape[0]]
SMAP_subset_int_indexes =  np.array(SMAP_subset_float_indexes).astype(int)       # Convert all items in the list to integers

cygnss_subset_float_indexes = [cygnss_difference_indexes[0],cygnss_difference_indexes[1],
                             cygnss_difference_indexes[2]+cygnss_lats.shape[0],cygnss_difference_indexes[3]+cygnss_lons.shape[0]]
cygnss_subset_int_indexes =  np.array(cygnss_subset_float_indexes).astype(int)       # Convert all items in the list to integers


# make the arrays for lats, lons and data
SMAP_lat_array, SMAP_lon_array = SMAP_lats[SMAP_subset_int_indexes[0]:SMAP_subset_int_indexes[2]], SMAP_lons[SMAP_subset_int_indexes[1]:SMAP_subset_int_indexes[3]]
SMAP_wind_asc_buoy = SMAP_wind_asc[SMAP_subset_int_indexes[0]:SMAP_subset_int_indexes[2],SMAP_subset_int_indexes[1]:SMAP_subset_int_indexes[3]]
SMAP_wind_des_buoy = SMAP_wind_des[SMAP_subset_int_indexes[0]:SMAP_subset_int_indexes[2],SMAP_subset_int_indexes[1]:SMAP_subset_int_indexes[3]]
SMAP_time_asc_buoy = SMAP_time_asc[SMAP_subset_int_indexes[0]:SMAP_subset_int_indexes[2],SMAP_subset_int_indexes[1]:SMAP_subset_int_indexes[3]]
SMAP_time_des_buoy = SMAP_time_des[SMAP_subset_int_indexes[0]:SMAP_subset_int_indexes[2],SMAP_subset_int_indexes[1]:SMAP_subset_int_indexes[3]]

#try to update SMAP mask


cygnss_lat_array, cygnss_lon_array = cygnss_lats[cygnss_subset_int_indexes[0]:cygnss_subset_int_indexes[2]], cygnss_lons[cygnss_subset_int_indexes[1]:cygnss_subset_int_indexes[3]]
cygnss_wind_buoy = cygnss_wind[:,cygnss_subset_int_indexes[0]:cygnss_subset_int_indexes[2],cygnss_subset_int_indexes[1]:cygnss_subset_int_indexes[3]]
# cygnss_time_buoy = cygnss_time[cygnss_subset_int_indexes[0]:cygnss_subset_int_indexes[2],cygnss_subset_int_indexes[1]:cygnss_subset_int_indexes[3]]

# get timeseries
cygnss_hourly_avg_buoy = np.average(np.average(cygnss_wind_buoy,axis=1),axis=1)[0:23] # remove the last value due to buoy missing last hour
cygnss_times_buoy = times = pd.date_range(start_time,end_time,freq='h')


# plotting
import matplotlib.pyplot as plt

# # plot timeseries
# fig, ax = plt.subplots()

# ax.plot_date(buoy_times, timeseries, 'o--')
# ax.plot_date(cygnss_times_buoy, cygnss_hourly_avg_buoy, 'o:')

# ax.tick_params(rotation=45)

# plt.show()

# plot cygnss
from cartopy import crs as ccrs 
# Create a Cartopy PlateCarree projection (cylindrical projection)
projection = ccrs.PlateCarree()

# Create a Matplotlib figure and axis
fig, ax = plt.subplots(subplot_kw={'projection': projection})
 
# im = ax.contourf(np.average(cygnss_wind,axis=0), origin='lower', extent=[cygnss_coverage[1], cygnss_coverage[3], cygnss_coverage[0], cygnss_coverage[2]], cmap='viridis',
#                 transform=projection, aspect='auto')
# im = ax.imshow(np.average(cygnss_wind_buoy,axis=0), origin='lower', extent=[buoy_area[1], buoy_area[3], buoy_area[0], buoy_area[2]], cmap='viridis',
#                 transform=projection, aspect='auto')
# im = ax.imshow(SMAP_wind_asc, origin='lower', extent=[SMAP_coverage[1], SMAP_coverage[3], SMAP_coverage[0], SMAP_coverage[2]], cmap='viridis',
#                 transform=projection, aspect='auto')
im = ax.contourf(ascat_lons, ascat_lats, ascat_wind, cmap='viridis',transform=projection)


# Add coastlines and gridlines for better context
ax.coastlines()
ax.gridlines()

# Add a colorbar for the plot
cbar = plt.colorbar(im, ax=ax, orientation='vertical')
cbar.set_label('wind speed avg (m/s)')

# saving file
directory = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Plots\MM2.png'
# fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.savefig(directory,  dpi=1200, bbox_inches="tight")
plt.show()

