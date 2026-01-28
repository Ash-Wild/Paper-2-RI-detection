# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 15:34:40 2023

@author: ashle
"""

# change this to be where the files are kept
wdir='C:/Users/ashle/OneDrive/Documents/Uni/PHD'

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cartopy import crs as ccrs # maybe not needed?
import glob

path = r'C:\Users\ashle\OneDrive\Documents\Uni\PHD\Data\L1 V3.0\\'
path = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Data\L1 V3.0\\'
# file_name = 'cyg.ddmi.s20170501-000002-e20170501-235959.l2.wind_trackgridsize25km_NOAAv1.2_L1a21.d21.nc'
files_list = glob.glob(path+ "\*.nc")

crs = "EPSG:4326"
dopp_res=[]

for file in files_list:
    noaa = netCDF4.Dataset(file)
    lon = noaa.variables['sp_lon'][:] ; lat = noaa.variables['sp_lat'][:] ;  snr=noaa.variables['ddm_snr'][:]; brcs = noaa.variables['brcs'][:]; angle = noaa.variables['sp_inc_angle'][:]
    lon[lon>180]= lon[lon>180]-360
    dopp_res.append(noaa.variables['dopp_resolution'][:])

    # df = pd.DataFrame({'SNR': snr[0:1000].flatten(), 'Longitude': lon[0:1000].flatten(), 'Latitude': lat[0:1000].flatten()})
    # gdf = geopandas.GeoDataFrame(
    #     df, geometry=geopandas.points_from_xy(x=df.Longitude, y=df.Latitude, crs=crs))
    
    # Info on plotting:
        # https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.plot.html
        # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
    # ax1= gdf.plot(column='SNR', figsize=(11,6),markersize=5, legend = True, cmap = 'viridis', legend_kwds={"fmt": ","})
    # world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
    # world.to_crs(gdf.crs).plot(ax=ax1,color='none',edgecolor = 'black')

