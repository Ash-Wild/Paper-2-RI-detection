# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 09:57:21 2024
purpose: create 2d histogram of layers
inputs: 2 arrays, and names of datasets, save directory 
output: saving maps 
@author: Ashley
"""

def density_plot(ref, ref_name, comparison, comparison_name, save_location,min_speed,max_speed,colourbar_upper,cbar_show):
    import matplotlib
    import matplotlib.pyplot as plt
    from cartopy import crs as ccrs 
    import cartopy.feature as cfeature
    import numpy as np
    from scipy.stats import linregress

    comp = 'Ashley'
    num_bins = 50
    colourbar_lower = 1
    colourbar_upper = colourbar_upper # TCA: 100, gauge: 1000
        
    # Perform linear regression to get the trendline
    slope, intercept, r_value, p_value, std_err = linregress(ref, comparison)
    R = np.around(np.corrcoef(ref,comparison)[0,1], decimals=2)
    RMSE = np.around(np.sqrt(np.mean((comparison - ref) ** 2)), decimals=2)
    bias = np.around(np.mean(comparison - ref), decimals=2)

    hist, xedges, yedges = np.histogram2d(ref, comparison, bins=num_bins)
    plt.figure(figsize=(5,5))
    plt.hist2d(ref, comparison, bins=num_bins, cmap='Blues', norm=matplotlib.colors.LogNorm(vmax=colourbar_upper,vmin=colourbar_lower), range = [[min_speed, max_speed], [min_speed, max_speed]])
    # plt.hist2d(ref, comparison, bins=num_bins, cmap='Blues', vmax=colourbar_upper,vmin=colourbar_lower, range = [[min_speed, max_speed], [min_speed, max_speed]])
    
    # # Add horizontal lines at the average value of each column
    # avg_y_per_column = [np.mean(comparison[(ref >= xedges[i]) & (ref < xedges[i+1])]) for i in range(len(xedges)-1)]
    # for i, avg in enumerate(avg_y_per_column):
    #     if not np.isnan(avg):
    #         plt.axhline(y=avg, color='green', linestyle='-', xmin=(xedges[i]-xedges.min())/(xedges.ptp()), xmax=(xedges[i+1]-xedges.min())/(xedges.ptp()))
    
    # Generate trendline values
    x_trend = np.linspace(ref.min(), ref.max(), 100)
    y_trend = slope * x_trend + intercept
    plt.plot(x_trend, y_trend, color='gray', linestyle='-', linewidth=2, label='Trendline')
    # Add black lines at x=0 and y=0
    plt.axhline(0, color='black', linewidth=1)  # Horizontal line at y=0
    plt.axvline(0, color='black', linewidth=1)  # Vertical line at x=0

    plt.xlabel(ref_name + ' (m/s)')
    plt.ylabel(comparison_name + ' (m/s)')
    plt.xlim(min_speed, max_speed)  # Change these values to your desired x bounds
    plt.ylim(min_speed, max_speed)  # Change these values to your desired y bounds
    
    plt.text(0.05,0.95, f'# Matches: {len(ref)}',transform=plt.gca().transAxes, va='top', ha='left')
    plt.text(0.05, 0.90, f'R: {R}', transform=plt.gca().transAxes, va='top', ha='left')
    plt.text(0.05, 0.85, f'RMSD: {RMSE}', transform=plt.gca().transAxes, va='top', ha='left')
    plt.text(0.05,0.80, f'Bias: {bias}',transform=plt.gca().transAxes, va='top', ha='left')
    # Plot the dashed line along the 1:1 diagonal
    plt.plot([min_speed, max_speed], [min_speed, max_speed], linestyle='--', color='black')
    
    if cbar_show:
        cbar = plt.colorbar(orientation='horizontal')
        cbar.set_label('# meas.')
    
    # # saving file
    plt.savefig(save_location, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()