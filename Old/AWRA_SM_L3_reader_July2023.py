# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 09:44:12 2023

Purpose: Read AWRA file and dispay
Input: file location of AWRA
Output: Datacube of AWRA
@author: ashle
"""


# change this to be where the files are kept - helps see then
wdir='C:/Users/Ashley/OneDrive - RMIT University/PHD//'

import netCDF4
import numpy as np
import glob

path = r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Data\AWRA\Monthly AWRA relative'
# file = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Data\AWRA\WRALv7_dd_2023.nc'
files_list = glob.glob(path+ "\*.nc")

# store all the locations of interest for sub-setting data. Lat lower, lon lower, lat upper, lon upper.
location_dictionary = {'Aus':(-44,112,-10,154), 'Vic':(-40,140.5,-34,150) }
data_coverage = (-44,112,-10,154) #change this depending on dataset, global would be -90,-180 (or 0), 90, 180 (or 360).
data_resolution = (0.05,0.05,0.05,0.05) # maybe not need this?


location = 'Aus' # make selection of location
location_bounds = location_dictionary[location]

difference = tuple(x - y for x, y in zip(location_dictionary[location], data_coverage))
difference_indexes = tuple(y/ x for x, y in zip(data_resolution, difference)) 

        
file_count = 0
for file in files_list:
    dataset = netCDF4.Dataset(file)
    key_var= dataset.variables['s0_pct']; time =dataset.variables['time']
    # Possibly subset it for region
    
    if file_count == 0:
        lon = dataset.variables['longitude'] ; lat = dataset.variables['latitude'] 
        subset_float_indexes = (difference_indexes[0],difference_indexes[1],difference_indexes[2]+lat.shape[0],difference_indexes[3]+lon.shape[0])
        subset_int_indexes = tuple(int(item) for item in subset_float_indexes)        # Convert all items in the tuple to integers
        lats = lat[subset_int_indexes[0]:subset_int_indexes[2]]
        lons = lon[subset_int_indexes[1]:subset_int_indexes[3]]
        datacube=np.zeros((36,lats.shape[0],lons.shape[0]))

    datacube =key_var[-42:-6,subset_int_indexes[0]:subset_int_indexes[2],subset_int_indexes[1]:subset_int_indexes[3]] # allocates the key variable to that layer of the datacube
    # datacube[file_count]=key_var[0,subset_int_indexes[2]*(-1):subset_int_indexes[0]*(-1),subset_int_indexes[1]:subset_int_indexes[3]] # allocates the key variable to that layer of the datacube
    file_count += 1
    
# Set the values that are -9999 to NaN (mask them out)
datacube_masked = np.ma.masked_equal(datacube, -999)

# separate into averages
import xarray as xr
import pandas as pd
# date_list = np.arange(start_date, end_date + np.timedelta64(1, 'D'), dtype='datetime64[D]').tolist() #this had an issue of wrong formatting
monthly_datacube = xr.DataArray(
    data=datacube_masked,
    dims=["time","y", "x"],
    coords=dict(
        lat=(["y"], lats),
        lon=(["x"], lons),
        time=pd.date_range("2020-01-01", periods=36, freq='M'), 
    ),
    attrs=dict(
        description="Top layer SM.",
        units="units",
    ),
    name='SM',
)

# Monthly and Seasonal averages
total_average = monthly_datacube.mean(axis=0)
monthly_average = monthly_datacube.groupby('time.month').mean(dim='time')
seasonal_average = monthly_datacube.groupby('time.season').mean(dim='time')

#convert to ease-2 grid
import pyproj
# Define CRS for source (lat/lon) and target (EASE-2 grid)
source_crs = pyproj.CRS.from_string('+proj=longlat +datum=WGS84')
target_crs = pyproj.CRS.from_string('+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs')  

# Create a coordinate transformer
transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)

# Reproject coordinates
ease2_x, ease2_y = transformer.transform(lons, lats)



# new lists of lats and lons to the index of ease-2 grid
from ease_lonlat import EASE2GRID, SUPPORTED_GRIDS
grid = EASE2GRID(name='EASE2_G36km', **SUPPORTED_GRIDS['EASE2_G36km'])
col_first, row_first = grid.lonlat2rc(lon=lons[0], lat=lats[0])
col_last, row_last = grid.lonlat2rc(lon=lons[-1], lat=lats[-1])
EASE2_array=np.zeros((col_last - col_first+1,row_last-row_first+1))

#calculate the number of AWRA cells in each cell of the EASE2 array
for la in lats:
    for lo in lons:
        col, row = grid.lonlat2rc(lon=lo, lat=la)
        EASE2_array[col-col_first,row-row_first]+=1



# # analysis
# ms = datacube.mean(axis=0)
# vs = datacube.var(axis=0)

# # try plot
# import matplotlib.pyplot as plt
# from cartopy import crs as ccrs 
# # Create a Cartopy PlateCarree projection (cylindrical projection)
# projection = ccrs.PlateCarree()

# # Create a Matplotlib figure and axis
# fig, ax = plt.subplots(subplot_kw={'projection': projection})

# # Plot the 2D array using imshow
# # You can customize the colormap, extent, etc. based on your data
# im = ax.imshow(ms, origin='upper', extent=[location_bounds[1], location_bounds[3], location_bounds[0], location_bounds[2]], cmap='viridis',
#                 transform=projection, aspect='auto')

# # Add coastlines and gridlines for better context
# ax.coastlines()
# ax.gridlines()

# # Add a colorbar for the plot
# cbar = plt.colorbar(im, ax=ax, orientation='vertical')
# cbar.set_label('relative SM (%)')

# # Set the title for the plot
# ax.set_title('Title')
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')

# # # saving file
# # directory = r'C:/Users/Ashley/OneDrive - RMIT University/PHD/Plots/MM.png'
# # # fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
# # plt.savefig(directory,  dpi=1200, bbox_inches="tight")

# # Show the plot
# plt.show()

