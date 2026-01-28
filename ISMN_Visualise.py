# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:44:18 2023
To visualise ISMN data 
@author: ashle
"""
# https://ismn.readthedocs.io/en/latest/examples/interface.html
from ismn.interface import ISMN_Interface
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt  
import cartopy.crs as ccrs
import os

# Specify the directory path
directory_path = r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Data\ISMN\Aus_only'

# Get a list of folders (subdirectories) in the specified directory
stations_list = [folder for folder in os.listdir(directory_path+'\OZNET') if os.path.isdir(os.path.join(directory_path+'\OZNET', folder))]


ismn_data = ISMN_Interface(directory_path)
location = 'Aus'
network = 'OZNET'
station = stations_list[0]
nu_days_in_ISMN = 610

# # examples
ids = ismn_data.get_dataset_ids(variable = 'soil_moisture', filter_meta_dict={'station':station})
ts_hourly, meta = ismn_data.read(ids, return_meta=True)
ts_masked = ts_hourly.mask(ts_hourly == -9999)
ts_masked.columns = ts_masked.columns.droplevel()
ts_xarray = xr.DataArray(
       data=ts_masked['soil_moisture'],
       dims=["time"],
       coords=dict(
           time=pd.date_range("2020-01-01", periods=14617, freq=pd.DateOffset(hours=1)),
       ),
       name='ISMN',
   )
ts_masked_daily = ts_xarray.resample(time='D').mean()
# ts_masked_hourly = ts_masked.drop(['soil_moisture_flag','soil_moisture_orig_flag'], axis = 1)

# ax = ts_masked.plot(figsize=(12,4), title=f'Time series for ID {ids}', ylim = (0,0.5), xlabel="Time [year]", ylabel="Soil Moisture [$m^3 m^{-3}$]")

# conditions = (ismn_data.metadata['elevation'].val > 3200) 
# ismn_data.metadata[conditions].index.to_list()

# #plot available station on a map
# fig, axs = plt.subplots(1, 1, figsize=(16,10), subplot_kw={'projection': ccrs.Robinson()})
# ismn_data.plot_station_locations('soil_moisture',  ax=axs, markersize=5, text_scalefactor=2)

# location_dictionary = {'Aus':(-38,110,-10,154), 'Vic':(-38,140.5,-34,150) }
# bounds=location_dictionary[location]
# axs.set_extent([bounds[1], bounds[3], bounds[0], bounds[2]])
# plt.show()

# comparing to CYGNSS. First extract the lat lon of station
metadata_dic = ismn_data[network][station ].metadata.to_dict()  #  Get metadata for the station
ISMN_lat,ISMN_lon = metadata_dic['latitude'][0][0],metadata_dic['longitude'][0][0]

# load in cygnss
cygnuss_path = r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Data\UCAR SM\Results'
Cyg_lats,Cyg_lons = np.load(cygnuss_path+'\Aus_UCAR_SM_lat',allow_pickle=True),np.load(cygnuss_path+'\Aus_UCAR_SM_lon',allow_pickle=True)
Cyg_datacube = np.load(cygnuss_path+'\Aus_daily_UCAR_SM_seasonal',allow_pickle=True)

def find_closest_index(lst, target):
    return min(range(len(lst)), key=lambda i: abs(lst[i] - target))

closest_index_lat = find_closest_index(Cyg_lats, ISMN_lat)
closest_index_lon = find_closest_index(Cyg_lons, ISMN_lon)

Cyg_datacube_nearby = Cyg_datacube[0:nu_days_in_ISMN,closest_index_lat-1:closest_index_lat+1, closest_index_lon-1:closest_index_lon+1]

# remove nan values to calculate correlation of known values
nan_indices=np.where(np.isnan(ts_masked_daily.data))
ts_without_nan = ts_masked_daily.data[~np.isnan(ts_masked_daily.data)]
mask = np.ones(ts_masked_daily.data.shape, dtype=bool)
mask[nan_indices] = False    
Cyg_datacube_nearby_without_nan = Cyg_datacube_nearby[mask]


import scipy.stats as stats
correlation_coefficient1, p_value1 = stats.pearsonr(ts_without_nan, Cyg_datacube_nearby_without_nan[:,0,0])
correlation_coefficient2, p_value2 = stats.pearsonr(ts_without_nan, Cyg_datacube_nearby_without_nan[:,1,0])
correlation_coefficient3, p_value3 =stats.pearsonr(ts_without_nan, Cyg_datacube_nearby_without_nan[:,0,1])
correlation_coefficient4, p_value4= stats.pearsonr(ts_without_nan, Cyg_datacube_nearby_without_nan[:,1,1])


