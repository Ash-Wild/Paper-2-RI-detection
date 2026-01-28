# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 11:12:04 2024
Purpose: plotting locations of RI
Input: NetCDF of RI locations
Output: Map of areas
@author: Ashley
"""

import netCDF4 as nc
import datetime
import numpy as np
comp = 'ashle'

RI_file = r'C:\Users\\' + comp + r'\OneDrive - RMIT University\PHD\Data\IBTrACS\RI_events_v2.nc'
RI_nc = nc.Dataset(RI_file)

time = RI_nc.variables['times'][:]
storm = RI_nc.variables['storm'][:]
latitude = RI_nc.variables['latitude'][:]
longitude = RI_nc.variables['longitude'][:]
ri_rd = RI_nc.variables['ri_rd'][:]
RI_nc.close()

nans = ~np.isnan(latitude).data # remove nans from loaded nc

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.patches as mpatches

# Create a Cartopy plot
fig, ax = plt.subplots() # subplot_kw={'projection': ccrs.PlateCarree()}

# Add map features
# ax.coastlines()
# ax.gridlines(draw_labels=True)

# Set the extent of the plot
# ax.set_extent([-179, 180, -80, 80]) # , crs=ccrs.PlateCarree()

# Add a bounding box
for i in range(len(latitude)):
    bbox = mpatches.Rectangle((longitude[i][nans[i]][0], latitude[i][nans[i]][0]), longitude[i][nans[i]][0] - longitude[i][nans[i]][-1], latitude[i][nans[i]][0] - latitude[i][nans[i]][-1],
                                  linewidth=1, edgecolor='red')
    ax.add_patch(bbox)

# Show the plot
plt.show()

# saving file
directory =r'C:\Users\ashle\OneDrive - RMIT University\PHD\Plots\RI_locations_v2.png'
fig.savefig(directory, format='png', dpi=600, bbox_inches='tight', pad_inches=0)