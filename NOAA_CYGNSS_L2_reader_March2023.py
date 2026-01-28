# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:29:56 2023

@author: ashle
The filenaming convention is as follows:
cyg.ddmi.sYYYYMMDD-HHmmSS-eYYYYMMDD-HHmmSS.l2.wind_trackgridsize25km_NOAAv1.x_L1a21.d21.nc

cyg: Spacecraft/Mission identifier, where "cyg" refers to the CYGNSS constellation and NASA EV mission title.
ddmi: Primary remote sensing subsystem known as the "Delay Doppler Mapping Instrument", or "ddmi".
sYYYYMMDD-HHmmSS: Date/time stamp representing the UTC date and time of the first data point in the file in ISO-8601 format. YYYYMMDD represents the 4-digit year, 2-digit month, and 2-digit day of month (respectively). The "s" prefix indicates this is the "starting" point representing the first retrieved measurement in the data file. HHmmSS represents the 2-digit hour, 2-digit minute, and 2-digit second (respectively).
eYYYYMMDD-HHmmSS: Date/time stamp representing the UTC date and time of the last data point in the file in ISO-8601 format. YYYYMMDD represents the 4-digit year, 2-digit month, and 2-digit day of month (respectively). The "e" prefix indicates this is the "ending" point representing the last retrieved measurement in the data file. HHmmSS represents the 2-digit hour, 2-digit minute, and 2-digit second (respectively).
l2: Data Processing Level 2.
wind: Primary Measurement, in which the "wind" represents the availability of wind speed data in the data files.
trackgridsize25km: Indicates the horizontal grid resolution of each grid cell, which is 25x25 kilometers.
NOAAv1.x: Indicates the provider (NOAA) as well as the version identifier for the processing of this primary data product (v1.x = version 1.x <x = 1 or 2>).
L1aAA: Algorithm version of the Level 1 (L1) input data, where AA is the numerical algorithm version identifier. E.g., "a21" = Algorithm Version 2.1. Note: the Algorithm version is numerically de-coupled from the Dataset Version (see below).
dVV: Dataset version of the Level 1 (L1) input data, where VV is the numerical dataset version identifier. E.g., "d21" = Dataset Version 2.1. Note: the Dataset version is numerically de-coupled from the Algorithm Version (see above).
nc: File type: "nc" = netCDF.

"""
# change this to be where the files are kept - helps see then
wdir='C:/Users/ashle/OneDrive/Documents/Uni/PHD'

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
# from cartopy import crs as ccrs # maybe not needed?
import glob
from Griddify_points import gridify

path = r'C:\Users\Ashley\OneDrive\Documents\Uni\PHD\Data\CyGNSS NOAA L2\\'
# file_name = 'cyg.ddmi.s20170501-000002-e20170501-235959.l2.wind_trackgridsize25km_NOAAv1.2_L1a21.d21.nc'
files_list = glob.glob(path+ "\*.nc")

crs = "EPSG:4326"

sizes = []
        
for file in files_list:
    noaa = netCDF4.Dataset(file)
    lon = noaa.variables['lon'] ; lat = noaa.variables['lat'] ; wind= noaa.variables['wind_speed']
    gridded_wind = gridify(lon, lat, wind) # one issue is that the sizes change for some reason
    # maybe need to load more variables as data flags etc
#     sizes.append(np.shape(gridded_wind))
    
    
# fig, ax = plt.subplots(1,1)
# ax.set_aspect('equal')
# plt.figure(figsize = (20,6))
# plt.contourf(gridded_wind)
# world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
# world.to_crs(crs).plot(ax=ax,color='none',edgecolor = 'black')


# # saving file
# directory = 'Downloads/MM.png'
# # fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
# plt.savefig(directory,  dpi=1200, bbox_inches="tight")
# plt.show()
# plt.clf()

# # Old approach for plotting - makes a grid of vectors which isn't good for analysis
# import pandas as pd
# import geopandas
# import shapely

# df = pd.DataFrame({'Wind': wind[:], 'Longitude': lonvals, 'Latitude': latvals})

# gdf = geopandas.GeoDataFrame(
#     df, geometry=geopandas.points_from_xy(x=df.Longitude, y=df.Latitude, crs=crs))

# # Info on plotting:
#     # https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.plot.html
#     # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.plot.html
# ax1= gdf.plot(column='Wind', figsize=(15,6), vmin=0,vmax=20,markersize=0.0025, cmap = 'viridis', legend_kwds={"fmt": ","})
# world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))
# world.to_crs(gdf.crs).plot(ax=ax1,color='none',edgecolor = 'black')

# # making geopandas grid from: https://james-brennan.github.io/posts/fast_gridding_geopandas/
# xmin, ymin, xmax,ymax = gdf.total_bounds
# cell_size = 0.25
# grid_cells = []
# lon_array = np.arange(xmin, xmax+cell_size, cell_size)
# lat_array = np.arange(ymin, ymax+cell_size, cell_size)

# for x0 in lon_array:
#     for y0 in  lat_array:
#         # bounds
#         x1 = x0-cell_size
#         y1 = y0+cell_size
#         grid_cells.append( shapely.geometry.box(x0, y0, x1, y1)  )
# cell = geopandas.GeoDataFrame(grid_cells, columns=['geometry'], crs=crs)

# # merge 2 GeoPandasDataframes
# merged = geopandas.sjoin(gdf, cell, how='left')

# # calculate average wind in each cell
# merged['mean_wind']=-9999 # make new field
# dissolve = merged.dissolve(by="index_right", aggfunc={'Wind':'mean'})
# cell.loc[dissolve.index, 'Wind'] = dissolve.Wind.values
# ax2 = cell.plot(column = 'Wind',figsize=(15,6), vmin=0,vmax=20, cmap = 'viridis')





