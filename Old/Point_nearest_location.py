# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:11:28 2023

@author: ashle
"""
# a function to find the index of the point closest pt
# (in squared distance) to give lat/lon value.
desired_coord =(50,200)
# nearest_index, nearest_lat, nearest_lon = getclosest_ij(latvals, lonvals, desired_coord[0], desired_coord[1])

def getclosest_ij(lats,lons,latpt,lonpt):
    # find squared distance of every point on grid
    dist_sq = (lats-latpt)**2 + (lons-lonpt)**2  
    # 1D index of minimum dist_sq element
    minindex_flattened = dist_sq.argmin()    
    # Get 2D index for latvals and lonvals arrays from 1D index
    return (minindex_flattened, lats[minindex_flattened],lons[minindex_flattened])
