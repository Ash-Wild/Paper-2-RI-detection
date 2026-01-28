# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 07:32:00 2023

@author: ashle
"""
# change this to be where the files are kept
wdir='C:/Users/ashle/OneDrive/Documents/Uni/PHD'

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import geopandas
from cartopy import crs as ccrs # maybe not needed?
import glob
from Griddify_points import gridify

path = r'C:\Users\ashle\OneDrive\Documents\Uni\PHD\Data\CyGNSS NOAA L2\\'
# file_name = 'cyg.ddmi.s20170501-000002-e20170501-235959.l2.wind_trackgridsize25km_NOAAv1.2_L1a21.d21.nc'
files_list = glob.glob(path+ "\*.nc")
   
# inputs - dataset as netCDF, lat, lon, search_distance, time in hours, search_duration
# output - list of measurements within that criteria
def extract_points(dataset, lat, lon, search_distance, time = 12, search_duration = 12):
    assert (-40 < lat < 40) & (0 < lon < 360) 
    
    lons = dataset.variables['lon'][:] ; lats = dataset.variables['lat'][:] ; times= dataset.variables['sample_time'][:]

#Code from here: https://stackoverflow.com/questions/41818927/how-to-subset-data-using-multidimensional-coordinates-using-python-xarray

    lat_bnds = [lat-search_distance, lat+search_distance]
    lon_bnds = [lon-search_distance, lon+search_distance]
    time_bnds = [time*3600 - search_duration*3600, time*3600 + search_duration*3600] # convert from hours to seconds

    lat_inds = np.where((lats > lat_bnds[0]) & (lats < lat_bnds[1]))[0]
    lon_inds = np.where((lons > lon_bnds[0]) & (lons < lon_bnds[1]))[0]
    time_inds = np.where((times > time_bnds[0]) & (times < time_bnds[1]))[0]

    # need to be merged all together for 1 list of values
    array_inds = np.concatenate((lat_inds, lon_inds, time_inds))
# https://stackoverflow.com/questions/11528078/determining-duplicate-values-in-an-array
    s = np.sort(array_inds, axis=None)
    rep_el = s[:-1][s[1:] == s[:-1]] # removes duplicated values
    # rep_el = s[1:][np.diff(s) == 0]
    rep_el2 =  rep_el[:-1][rep_el[1:] == rep_el[:-1]] # removes duplicates again - so now only list of all 3
    
    if len(rep_el2) == 0 :        
        raise ValueError('No values found. Please enter different search parameters')
        
    return rep_el2
    #return dataset.variables['wind_speed'][:,lat_inds,lon_inds, time_inds]
       

for file in files_list:
    noaa = netCDF4.Dataset(file)

lat, lon, search_distance = 35, 10, 10
inds = extract_points(noaa, lat, lon, search_distance,10,5)

   # plotting
   
import pandas as pd
import geopandas
import shapely

lonvals = noaa.variables['lon'][inds] ; latvals = noaa.variables['lat'][inds] ; wind= noaa.variables['wind_speed'][inds]
df = pd.DataFrame({'Wind': wind, 'Longitude': lonvals, 'Latitude': latvals})
crs = "EPSG:4326"    

gdf = geopandas.GeoDataFrame(
df, geometry=geopandas.points_from_xy(x=df.Longitude, y=df.Latitude, crs=crs))
# Info on plotting:
    # https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.plot.html
    # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
ax1= gdf.plot(column='Wind', figsize=(15,6), vmin=0,vmax=20,markersize=0.25, cmap = 'viridis', legend_kwds={"fmt": ","})

# add in a box https://stackoverflow.com/questions/66464565/how-to-update-a-shape-on-matplotlib
x_box = np.array([lon-search_distance, lon-search_distance, lon+search_distance, lon+search_distance, lon-search_distance])
y_box = np.array([lat-search_distance, lat+search_distance, lat+search_distance, lat-search_distance, lat-search_distance])
plt.plot(x_box, y_box, c='r') 
# world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
# world.to_crs(gdf.crs).plot(ax=ax1,color='none',edgecolor = 'black')

   
# lon = noaa.variables['lon'][inds] ; lat = noaa.variables['lat'][inds] ; wind= noaa.variables['wind_speed'][inds]
# gridded_wind = gridify(lon, lat, wind) # one issue is that 

# fig, ax = plt.subplots(1,1)
# ax.set_aspect('equal')
# plt.figure(figsize = (20,6))
# plt.contourf(gridded_wind)
# # world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
# # world.to_crs(crs).plot(ax=ax,color='none',edgecolor = 'black')

    
    