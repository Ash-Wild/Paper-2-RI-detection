# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 11:28:02 2024
input: Storm files and files
output: list of matchups CYG storm files and Stations
purpose: see if there are matchups
@author: ashle
"""
import glob 
import pandas as pd
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

comp = 'ashle'
# data_folder = '\BOM Stations' OR ' PSLGM' OR 'Aus 1hr'
data_folder = r'\Aus 1hr'
csv_folder  = r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\BOM stations' + data_folder
cyg_folder = r'C:\Users\\'+comp+r'\OneDrive - RMIT University\PHD\Data\cyg_storm_centric'
sar_folder = r'C:\Users\\'+comp+r'\OneDrive - RMIT University\PHD\Data\sar full'
cyg_files_list = glob.glob(cyg_folder+ r"\*.nc")
csv_files_list = glob.glob(csv_folder+ r"\HC06D_Data_*.txt")
sar_files_list = glob.glob(sar_folder + r'\*.nc')
search_distance = 0.3 # number of degrees for searching for CYGNSS
cyg_time_thresh = 1 # number of hours to look for CYGNSS outside of event 

location_dict = {
    'Federated States of Micronesia': {'lat': 6.978, 'lon': 158.197},
    'Kiribati': {'lat': 1.362, 'lon': 172.93},
    'Marshall Islands': {'lat': 7.108, 'lon': 171.371},
    'Nauru': {'lat': -0.532, 'lon': 166.909},
    'Tuvalu': {'lat': -8.503, 'lon': 179.209},
}

station_coords = np.zeros((len(location_dict)+len(csv_files_list),2))
i = 0
for csv in csv_files_list:
    df = pd.read_csv(csv)
    station_coords[i,0],station_coords[i,1] = df['Latitude'][0],df['Longitude'][0]
    i+=1
    
for location in location_dict:
    station_coords[i,0],station_coords[i,1] = location_dict[location]['lat'],location_dict[location]['lon']
    i+=1

latitudes =  station_coords[0:len(csv_files_list),0]
longitudes = station_coords[0:len(csv_files_list),1]
latitudes2 =  station_coords[len(csv_files_list)::,0]
longitudes2 = station_coords[len(csv_files_list)::,1]

# Create a map
plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
ax.set_extent([110, 185, -40, 10], crs=ccrs.PlateCarree())  # [min_lon, max_lon, min_lat, max_lat]
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS)
ax.coastlines('50m')

# ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgreen')
# ax.add_feature(cfeature.OCEAN, edgecolor='black', facecolor='lightblue')
ax.stock_img()

# Plot the points on the map
plt.scatter(longitudes, latitudes, color='red', marker='o', transform=ccrs.PlateCarree())
plt.scatter(longitudes2, latitudes2, color='yellow', marker='o', transform=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.xlocator = plt.FixedLocator(range(-180, 181, 10))  # Longitude lines every 10 degrees
gl.ylocator = plt.FixedLocator(range(-90, 91, 10))   # Latitude lines every 10 degrees

gl.top_labels = False
gl.right_labels = False
gl.xlabel_style = {'size': 10, 'color': 'black'}
gl.ylabel_style = {'size': 10, 'color': 'black'}
    
directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots\station_domain.png'

plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()


# cyg_station_match = []
# for cyg_file in cyg_files_list[0:1]:    
#     cyg_nc = nc.Dataset(cyg_file)
#     cyg_lons = cyg_nc.variables['lon'][:] ; cyg_lats = cyg_nc.variables['lat'][:];
    
    
#     cyg_time = cyg_nc.variables['time'][:] ; cyg_time_offset = cyg_nc.variables['time_offset'][:]
#     cyg_wind= cyg_nc.variables['wind_speed'][:] ; cyg_uncertainty = cyg_nc.variables['wind_speed_uncertainty'][:]