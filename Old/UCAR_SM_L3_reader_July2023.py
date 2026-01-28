"""
Code to visualise UCAR L3 SM data
Input: folder containing data
Output: Data cube of values
"""

# change this to be where the files are kept - helps see then
wdir='C:/Users/Ashley/OneDrive - RMIT University/PHD//'

import netCDF4
import numpy as np
import glob
from datetime import date, timedelta
import xarray as xr
import pandas as pd

path = r'C:\Users\Ashley\OneDrive - RMIT University\PHD\Data\UCAR SM'

#select changes
location = 'Aus' # make selection of location
start_date = date(2020,1,1)
end_date = date(2022,12,31)


def process_UCAR_SM(path,location,start_date,end_date):
    '''
    inputs: the path to the daily data files, and the location which we want to cut the datacube to. 
    outputs: saves the file 
    '''
    # file_name = 'cyg.ddmi.s20170501-000002-e20170501-235959.l2.wind_trackgridsize25km_NOAAv1.2_L1a21.d21.nc'
    files_list = glob.glob(path+ "\*.nc4")
    # store all the locations of interest for sub-setting data. Lat lower, lon lower, lat upper, lon upper.
    location_dictionary = {'Aus':(-38,110,-10,154), 'Vic':(-38,140.5,-34,150) }
    data_coverage = (-38,-135,38,164) #change this depending on dataset, global would be -90,-180 (or 0), 90, 180 (or 360).
    data_resolution = (0.3,0.37,0.3,0.37) # maybe not need this?
        
    difference = tuple(x - y for x, y in zip(location_dictionary[location], data_coverage))
    difference_indexes = tuple(y/ x for x, y in zip(data_resolution, difference)) # feel iffy as it's EASE grid not normal proj
    
    reference_date = date(1970,1,1)
    times = pd.date_range(start_date,end_date,freq='d')
    file_count = 0
    for file in files_list:
        dataset = netCDF4.Dataset(file)
        key_var= dataset.variables['SM_daily']
        # subset it for region
        
        if file_count == 0:
            lon = dataset.variables['longitude'] ; lat = dataset.variables['latitude'] 
            subset_float_indexes = (difference_indexes[0],difference_indexes[1],difference_indexes[2]+lat.shape[0],difference_indexes[3]+lat.shape[1])
            subset_int_indexes = tuple(int(item) for item in subset_float_indexes)        # Convert all items in the tuple to integers
            lat_array, lon_array = lat[subset_int_indexes[0]:subset_int_indexes[2],0], lon[0,subset_int_indexes[1]:subset_int_indexes[3]]
            previous_date = date(2019,12,31)
            datacube=np.zeros((len(times),subset_int_indexes[2]-subset_int_indexes[0],subset_int_indexes[3]-subset_int_indexes[1]))
        
        # need to fill in days with no data as NaN
        current_date = reference_date + timedelta(days=int(dataset.variables['time'][:][0]))
        
        if previous_date + timedelta(days=1) != current_date:
            print('days skipped ')
            # Calculate the timedelta between the two dates
            time_difference = current_date - previous_date -timedelta(days=1)
            days_skipped = time_difference.days
            for i in range(days_skipped):
                datacube[file_count]=np.nan
                file_count += 1
        datacube[file_count]=key_var[0,subset_int_indexes[0]:subset_int_indexes[2],subset_int_indexes[1]:subset_int_indexes[3]] # allocates the key variable to that layer of the datacube
        file_count += 1
        previous_date = current_date
        
    # Set the values that are -9999 to NaN (mask them out)
    datacube_masked = np.ma.masked_equal(datacube, -9999)
    datacube_masked.dump(path +'/Results/' + location+'_daily_UCAR_SM_seasonal')
    lat_array.dump(path +'/Results/' + location+'_UCAR_SM_lat')
    lon_array.dump(path +'/Results/' + location+'_UCAR_SM_lon')
    lat[subset_int_indexes[0]:subset_int_indexes[2]].dump(path +'/Results/' + location+'_UCAR_SM_lat_full')
    lon[subset_int_indexes[1]:subset_int_indexes[3]].dump(path +'/Results/' + location+'_UCAR_SM_lon_full')
 

    # separate into averages
    # date_list = np.arange(start_date, end_date + np.timedelta64(1, 'D'), dtype='datetime64[D]').tolist() #this had an issue of wrong formatting
    xr_datacube = xr.DataArray(
        data=datacube_masked,
        dims=["time","y", "x"],
        coords=dict(
            lat=(["y"], lat_array),
            lon=(["x"], lon_array),
            time=times,
        ),
        attrs=dict(
            description="Top layer SM.",
            units="units",
        ),
        name='SM',
    )
    
    monthly_datacube = xr_datacube.resample(time='1M').mean()
    monthly_datacube.to_netcdf(path +'/Results/' + location+'_monthly_UCAR_SM_seasonal.nc')
    monthly_datacube.values.dump(path +'/Results/' + location+'_monthly_UCAR_SM_seasonal')

    # # remove seasonality
    # import statsmodels.api as sm
    # # Convert the datacube to a pandas DataFrame with multi-index (lat, lon)
    # data_df = monthly_datacube.to_dataframe()
    # # Perform seasonal decomposition using statsmodels on each pixel
    # seasonal_decomposition = data_df.groupby(['lat', 'lon'])['SM'].apply(
    #     lambda x: sm.tsa.seasonal_decompose(x, model='additive', period=12)
    # )
    # # Remove the seasonal component to get the deseasonalized data
    # data_df['deseasonalized'] = data_df['variable'] - seasonal_decomposition['seasonal']
    # # Convert the deseasonalized data back to a datacube
    # deseasonalized_datacube = data_df['deseasonalized'].unstack().to_xarray()
      
    # deseasonalized_datacube.to_netcdf(path +'/' + location+'_monthly_UCAR_SM_deseasonalised.nc')
    # deseasonalized_datacube.values.dump(path +'/Results/' + location+'_monthly_UCAR_SM_deseasonalised')
    
    return monthly_datacube

monthly_datacube = process_UCAR_SM(path, location,start_date,end_date)


# Monthly and Seasonal averages
total_average = monthly_datacube.mean(axis=0)
monthly_average = monthly_datacube.groupby('time.month').mean(dim='time')
seasonal_average = monthly_datacube.groupby('time.season').mean(dim='time')


# # try plot
# import matplotlib.pyplot as plt
# from cartopy import crs as ccrs 
# # Create a Cartopy PlateCarree projection (cylindrical projection)
# projection = ccrs.PlateCarree()

# # Create a Matplotlib figure and axis
# fig, ax = plt.subplots(subplot_kw={'projection': projection})
 
# # Create a figure with subplots for each season
# fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)

# # Plot each season on a different subplot
# for i, season in enumerate(seasonal_average['season']):
#     row = i // 2
#     col = i % 2
#     seasonal_average.sel(season=season).plot(ax=axs[row, col], marker='o', linestyle='-')
#     axs[row, col].set_title(season)
#     axs[row, col].grid(True)

# # Set common labels for x and y axes
# fig.text(0.5, 0.04, 'Months', ha='center')
# fig.text(0.04, 0.5, 'Average Value', va='center', rotation='vertical')

# plt.suptitle('Seasonal Averages', fontsize=16)
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.show()

# # Plot the 2D array using imshow
# # You can customize the colormap, extent, etc. based on your data
# im = ax.imshow(datacube_masked[0], origin='lower', extent=[location_bounds[1], location_bounds[3], location_bounds[0], location_bounds[2]], cmap='viridis',
#                 transform=projection, aspect='auto')

# # Add coastlines and gridlines for better context
# ax.coastlines()
# ax.gridlines()

# # Add a colorbar for the plot
# cbar = plt.colorbar(im, ax=ax, orientation='vertical')
# cbar.set_label('SM (cm3/cm3)')

# # Set the title for the plot
# ax.set_title('Title')


# # saving file
# directory = r'C:\Users\ashle\OneDrive - RMIT University\PHD\Plots\MM.png'
# # fig.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
# plt.savefig(directory,  dpi=1200, bbox_inches="tight")

# # Show the plot
# plt.show()