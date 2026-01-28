# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 11:55:29 2023

@author: ashle

"""

import numpy as np
import gdal
import geopandas
import pandas as pd
import scipy

# approach 2 https://stackoverflow.com/questions/54842690/how-to-efficiently-convert-large-numpy-array-of-point-cloud-data-to-downsampled
def f_pp(pcloud_np, resolution):
    xy = pcloud_np.T[:2] # array of just lat and lon
    xy = ((xy + resolution / 2) // resolution).astype(int) # gives an index for each place
    mn, mx = xy.min(axis=1), xy.max(axis=1) # get the min and max indexes
    sz = mx + 1 - mn # total size of the array, size of z (Wind for us)
    flatidx = np.ravel_multi_index(xy-mn[:, None], sz)
    histo = np.bincount(flatidx, pcloud_np[:, 2], sz.prod()) / np.maximum(1, np.bincount(flatidx, None, sz.prod()))
    return (histo.reshape(sz), *(xy * resolution))

def do_kdtree(combined_x_y_arrays,points):
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes

def gridify(lon, lat, wind, resolution = 0.25):
    # extract lat/lon values (in degrees) to numpy arrays
    lonvals = np.append(lon[:],[-180,180]); latvals = np.append(lat[:],[-38,38]);  windvals = np.append(wind[:],[np.NaN,np.NaN])
    # shift the end axis 
    lonvals[lonvals>180]= lonvals[lonvals>180]-360
    
        
    # EASE2_gdf = gdf.to_crs("EPSG:6933")
    
# approach 4 - EASE 3km grid
    # INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 6933, 'EPSG', 6933, '+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs', 'PROJCS["WGS 84 / NSIDC EASE-Grid 2.0 Global",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Cylindrical_Equal_Area"],PARAMETER["standard_parallel_1",30],PARAMETER["central_meridian",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","6933"]]');
    from ease_grid import EASE2_grid
    egrid = EASE2_grid(3000)
    # these two attributes contain the longitude and latitude coordinate dimension
    xs = egrid.londim
    ys = egrid.latdim
    grid = np.meshgrid(xs,ys)
    combined_x_y_arrays = np.dstack([grid[0].ravel(),grid[1].ravel()])[0]
    points = np.dstack([lonvals, latvals])
    results = do_kdtree(combined_x_y_arrays,points)

    
    # # Step 3: Determine grid cell for each point
    # points['grid_cell'] = points.geometry.apply(lambda p: egrid.loc[grid.contains(p)].iloc[0]['grid_id'])
    
    # # Step 4: Perform the join
    # joined_data = grid.merge(points, left_on='grid_id', right_on='grid_cell')

    df = pd.DataFrame({'Wind': windvals, 'Longitude': lonvals, 'Latitude': latvals})
    gdf = geopandas.GeoDataFrame(
        df, geometry=geopandas.points_from_xy(x=df.Longitude, y=df.Latitude, crs="EPSG:4326"))
        
    
    # # Approach 3    # https://gis.stackexchange.com/questions/396995/using-geopandas-geodataframe-in-gdal-grid-for-spatial-interpolation-viz-idw-nea
    # GDAL_dataset = gdal.OpenEx(gdf.to_json(), gdal.OF_VECTOR)
    # GDAL_gridded = gdal.Grid('out_grid.tif',GDAL_dataset) # ERROR 1: Size and resolutions are missing
    # # gdal.Rasterize('out_grid2.tif',GDAL_dataset, )
    # grid_wind = np.array(GDAL_gridded.GetRasterBand(1).ReadAsArray())
    # return grid_wind

    # pcloud_np=np.transpose(np.array([lonvals, latvals, Windvals])) # rotated to be 3 columns
    # grid_wind, x, y = f_pp(pcloud_np, resolution)
    # grid_wind[grid_wind==0]= np.NaN # a crude way of applying the mask from previous data
    # grid_wind_t = np.transpose(grid_wind)
    # return grid_wind_t



# # approach 1 for griddifying: https://earthscience.stackexchange.com/questions/12057/how-to-interpolate-scattered-data-to-a-regular-grid-in-python

# # target grid to interpolate to
# xi = np.arange(-180,180,resolution)
# yi = np.arange(-36,36,resolution)
# xi,yi = np.meshgrid(xi,yi)

# # interpolate
# zi = griddata((lonvals,latvals),Windvals,(xi,yi),method='linear') # issue here that it cannot cope with masked

# # plot
# fig = plt.figure()
# ax = fig.add_subplot(111)
# plt.contourf(xi,yi,zi)
# plt.plot(lonvals,latvals,data=Windvals)
# plt.xlabel('xi',fontsize=16)
# plt.ylabel('yi',fontsize=16)
# plt.savefig('interpolated.png',dpi=100)
# plt.close(fig)