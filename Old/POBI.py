# -*- coding: utf-8 -*-
"""
Created on Fri May  5 11:25:06 2023

@author: ashle
"""
import scipy.stats as sp
import numpy as np


def POBI(today_grid,cygnss_grid,ref_grid,radius):
    assert np.size(cygnss_grid) == np.size(ref_grid)
    assert np.size(today_grid) == np.size(cygnss_grid[0,:,:])
    assert type(radius) is int
    for i in range(radius, np.size(cygnss_grid)[1]-radius):
        for j in range(radius, np.size(cygnss_grid)[2]-radius):
            # needed to avoid border cases of no comparisons
            if today_grid[i,j] is np.nan: # keep original values, only do for times not covered
                sum_1 = 0
                sum_2 = 0
                for x in range(i-radius,i+radius):
                    for y in range(j-radius,j+radius):
                        if today_grid[x,y] is not np.nan:
                            #pearson = sp.pearsonr(cygnss_grid[:,i,j],ref_grid[:,x,y])
                            pearson = np.corrcoef(cygnss_grid[:,i,j],ref_grid[:,x,y])**2
                            # https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
                            a,b= np.linalg.lstsq(cygnss_grid[:,i,j],ref_grid[:,x,y])
                            sum_1 += pearson*(a*today_grid[x,y]+b)
                            sum_2 += pearson
                if sum_2 != 0:
                    today_grid[i,j] = sum_1/sum_2
    return today_grid
            