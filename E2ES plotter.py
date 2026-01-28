# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 11:03:42 2024

@author: Ashley
"""

import netCDF4 as nc
import glob
import matplotlib.pyplot as plt
import numpy as np

comp = 'Ashley'
E2ES_nc_results = r"\\wsl.localhost\Ubuntu-20.04\home\awild\e2es_160615b\run.20100825\output"
files_list = glob.glob(E2ES_nc_results+ "\*.nc")

nc = nc.Dataset(files_list[0])
DDMs = nc.variables['DDMs']

max_count = np.max(DDMs)
min_count = 6e-18
# min_count = np.min(DDMs[0,0,:,:])

nrows= np.shape(DDMs)[0]
ncols=np.shape(DDMs)[1]//2+1
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(30, 4))

for i in range(nrows):
    for j in range(ncols):
        ax=axes[i,j]
        im = ax.imshow(DDMs[i][j], cmap='viridis', origin='lower', interpolation='none', vmin=min_count, vmax=max_count)
        ax.set_title(('{0} {1}').format(i,j))

cbar_ax = fig.add_axes([0.96, 0.1, 0.02, 0.8])  # [left, bottom, width, height]
cbar = fig.colorbar(im, cax=cbar_ax)

fig.supylabel('doppler bins')
fig.supxlabel('delay bins')
plt.tight_layout(rect=[0.01, 0, 0.96, 0.96])  # Adjust layout to make space for the common labels

# # saving file
directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Plots\E2ES DDMs.png'
plt.savefig(directory, format='png', dpi=300,bbox_inches='tight', pad_inches=0)

plt.show()

