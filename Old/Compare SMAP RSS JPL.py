# """"""
# Created on Tue Mar 26 15:56:37 2024
# Purpose: compare SMAP from JPL and RSS

# @author: Ashley
# """"""
import netCDF4 as nc
import numpy as np
import glob 
import h5py
from datetime import timedelta,datetime
import matplotlib.pyplot as plt
from cartopy import crs as ccrs 

comp = 'ashle'

rss_nc = nc.Dataset(r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\All_SMAP\RSS\RSS_smap_wind_daily_2023_10_23_v01.0.nc')
rss_lons = rss_nc.variables['lon'][:]-180 ; rss_lats = rss_nc.variables['lat'][:]; rss_wind= np.mean(rss_nc.variables['wind'][:],axis=2) ; rss_time =np.mean(np.ma.masked_array(rss_nc.variables['minute'][:].data,mask=rss_nc.variables['wind'][:].mask), axis=2)


L2_dir_before = r'C:\Users\\' +comp + '\OneDrive - RMIT University\PHD\Data\TC Lola Oct 23 2023\SMAP_during'
files_list = glob.glob(L2_dir_before+ "\*.h5")
jpl_lats=[]; jpl_lons=[]; jpl_winds=[]
for files in files_list:
    file = h5py.File(files, 'r')
    ascat_wind = file["/smap_high_spd"][:] ; ascat_lats = file["/lat"][:] ; ascat_lons = file['lon'][:]; ascat_time = file['row_time'][:].astype(np.float64) ; ascat_uncertainty = file['smap_ambiguity_spd'][:]
    jpl_lats.append(ascat_lats.flatten()); jpl_lons.append(ascat_lons.flatten()); jpl_winds.append(ascat_wind.flatten())
    # start_date = datetime(2015, 1, 1, 0, 0, 0, 0)
    # datetime_array_1D = np.vectorize(lambda x: start_date + timedelta(seconds=x))(ascat_time.data)
jpl_lats_merged=np.concatenate(jpl_lats);jpl_lons_merged=np.concatenate(jpl_lons);jpl_winds_merged=np.concatenate(jpl_winds)
jpl_lats_filtered = jpl_lats_merged[jpl_winds_merged>0]
jpl_lons_filtered = jpl_lons_merged[jpl_winds_merged>0]
jpl_winds_filtered = jpl_winds_merged[jpl_winds_merged>0]


# plotting RSS
projection = ccrs.PlateCarree()
fig, ax = plt.subplots(subplot_kw={'projection': projection})
vmin = 10 ; vmax = 80

im = ax.scatter(jpl_lons_filtered, jpl_lats_filtered, c=jpl_winds_filtered , cmap='viridis',transform=projection, marker='s',vmin=vmin, vmax=vmax, s=1)            

# im = ax.pcolormesh(rss_lons, rss_lats, rss_wind, cmap='viridis',transform=projection, vmin=vmin, vmax=vmax)
cbar = plt.colorbar(im, ax=ax, orientation='vertical',location='left', shrink=0.8)
cbar.set_label('wind speed avg (m/s)')

# #convert time to datetime format
# lats_boolean = (smap_lats.data > AOI_lat_min) & (smap_lats.data < AOI_lat_max)
# lons_boolean = (smap_lons.data > AOI_lon_min) & (smap_lons.data < AOI_lon_max)
# smap_wind_clipped = smap_wind[lats_boolean][:,lons_boolean]
# # print(np.where(smap_wind_clipped >0))

# smap_start_date = datetime.combine(smap_file_dates[smap_ind[0]], datetime.min.time())
# datetime_array = np.vectorize(lambda x: smap_start_date + timedelta(seconds=x*60))(smap_time[lats_boolean][:,lons_boolean])
# print(np.where(~datetime_array.mask))
# time_boolean = (datetime_array > AOI_time_start) & (datetime_array < AOI_time_end) # error here
