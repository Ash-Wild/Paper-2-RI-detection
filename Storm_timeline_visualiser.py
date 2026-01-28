"""
Created on Tue Jul 23 17:48:38 2024
Timeline visualiser
Input: times and values of station, and CYGNSS files
output: map showing timeline shift
@author: Ashley
"""

import glob
from datetime import datetime, timedelta 
import numpy as np
import pickle
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# from BOM_high_speed_finder import wind_thresh, time_thresh

# input files from PSLGM_high_speed_finder
comp = 'ashle'  # computer username
cyg_version = 'noaa_1.2' # 'noaa_1.2' OR 'l2_3.2' or 'storm_centric'
data_folder = '\Aus 1hr' # data_folder = '\BOM Stations' OR ' PSLGM' OR 'Aus 1hr'
csv_folder  = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\BOM stations' + data_folder
search_distance = 0.5 # number of gedrees for searching for CYGNSS
cyg_time_thresh = 1 # number of hours to look for CYGNSS outside of event 
cyg_AOI = False # will switch to true if there's some overlap

# Open the file and load the list
with open(csv_folder+'\events_list.pkl', 'rb') as f:
    events_list = pickle.load(f)

location_dict = {
    'Cook Islands': {'lat': -21.199, 'lon': -159.786},
    'Federated States of Micronesia': {'lat': 6.978, 'lon': 158.197},
    ' Fiji': {'lat': -17.605, 'lon': 177.438},
    'Kiribati': {'lat': 1.362, 'lon': 172.93},
    'Marshall Islands': {'lat': 7.108, 'lon': 171.371},
    'Nauru': {'lat': -0.532, 'lon': 166.909},
    'Papua New Guinea': {'lat': -2.036, 'lon': 147.375},
    'Samoa': {'lat': -13.819, 'lon': -171.756},
    'Solomon Islands': {'lat': -9.422, 'lon': 159.955},
    'Tonga': {'lat': -21.14, 'lon': -175.179},
    'Tuvalu': {'lat': -8.503, 'lon': 179.209},
    'Vanuatu': {'lat': -17.761, 'lon': 168.293},
    'Niue': {'lat': -19.053, 'lon': -169.921}
}

cyg_folder = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\cyg_'+cyg_version
cyg_files_list = glob.glob(cyg_folder+ "\*.nc")

# find the dates of cyg 
# make a list of cygnss dates and SMAP dates
cyg_file_dates = []
cyg_files_num = []
file_counter = 0
for cyg_file in cyg_files_list:
    if cyg_version == 'storm_centric':
        cyg_nc = nc.Dataset(cyg_file)
        date_format = "%y-%m-%d %H:%M:%S"
        cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
        cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
        cyg_datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(hours=int(x)))(cyg_nc.variables['time'][:])
        for i in cyg_datetime_array:
            cyg_file_dates.append(i.date())
            cyg_files_num.append(file_counter)
        file_counter += 1
    else: 
        start_ind = cyg_file.find('ddmi.s')
        cyg_datetime_object = datetime.strptime(cyg_file[start_ind+6:start_ind+14], "%Y%m%d")
        cyg_file_dates.append(cyg_datetime_object.date())


def plot_distribution(cyg_df, ref_df,merged_df):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    plt.rcParams.update({
    'font.size': 14,            # Applies to most text elements
    'axes.titlesize': 16,       # Title size
    'axes.labelsize': 14,       # Axis label size
    'xtick.labelsize': 12,      # X tick label size
    'ytick.labelsize': 12,      # Y tick label size
    'legend.fontsize': 12       # Legend font size
    })
    plt.figure(figsize=(10, 6))

    # Setup: define each dataset with label and color
    datasets = [[ref_df, 'IBTrACS', 'blue'],[cyg_df, 'NOAA', 'red'],[merged_df,'Merged','brown']]

    # Loop through each dataset
    for df, label, color in datasets:
        df['time_since_RI'] = df['time_since_ri'].dt.total_seconds() / 3600

        # Create 30-minute bins from -12 to 36 hours
        bins = np.arange(-12, 36, 2)
        labels_bin = bins[:-1] + 0.25  # center of each bin

        # Assign bins
        df['time_bin'] = pd.cut(df['time_since_RI'], bins=bins, labels=labels_bin)

        # Group and compute stats
        grouped = df.groupby('time_bin')
        mean_ws = grouped['wind_speed'].mean()
        lower_ws = grouped['wind_speed'].quantile(0.05)
        upper_ws = grouped['wind_speed'].quantile(0.95)
        
        if label=='Merged':
            # Interpolate stats to full bin range
            mean_ws = mean_ws.reindex(labels_bin).interpolate(method='linear')
            lower_ws = lower_ws.reindex(labels_bin).interpolate(method='linear')
            upper_ws = upper_ws.reindex(labels_bin).interpolate(method='linear')

        # Plot
        plt.plot(labels_bin, mean_ws, label=f'{label} mean', color=color)
        plt.fill_between(labels_bin, lower_ws, upper_ws, color=color, alpha=0.2, label=f'{label} 95% range')

    # Set consistent x and y limits
    plt.xlim(-12, 36)
    plt.ylim(0, 66)

    # Add vertical lines and annotation
    plt.axvline(0, color='grey', linestyle='--')
    plt.axvline(24, color='grey', linestyle='--')
    plt.text(12, plt.ylim()[1]*0.95, 'RI period', color='grey', ha='center', fontsize=20, fontstyle='italic')

    # Final plot labels
    plt.xlabel('Hours since RI onset')
    plt.ylabel('Wind speed (m/s)')
    # plt.title('Windspeed Evolution Relative to RI Onset')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\Timeseries\Merged_NOAA_IBtracs_comparison.png'.format(comp)
    plt.savefig(directory, format='png', dpi=500, bbox_inches='tight', pad_inches=0)
    plt.show()

def plot_timeline_refined(cyg_wind,cyg_time,cyg_distances,ref_times, ref_winds,name,ri_start,ri_end):
    # arrays needed for final results
    cyg_times = cyg_time
    cyg_winds = cyg_wind
    dist_cyg_stat = cyg_distances
 
    # Normalize dist_cyg_stat to be between 0 and 1 for color mapping
    norm = mcolors.Normalize(vmin=0, vmax=80) # vmin=min(dist_cyg_stat), vmax=max(dist_cyg_stat)
    
    # Start plotting
    plt.figure(figsize=(10, 5))
    # Use `c` instead of `alpha` for color mapping
    sc = plt.scatter(cyg_times, cyg_winds, label='CYGNSS', c=dist_cyg_stat, norm=norm, cmap=cm.Reds_r, alpha=0.7)
    plt.plot(ref_times, ref_winds, label='IBTrACS', color='blue')
    plt.xticks(rotation=45)
    plt.xlabel('Time')
    plt.ylabel('Wind speed (m/s)')
    plt.legend()
    # Add grey dashed vertical lines
    plt.axvline(ri_start, color='grey', linestyle='--')
    plt.axvline(ri_end, color='grey', linestyle='--')
    # Add label in between
    midpoint = ri_start + (ri_end - ri_start) / 2
    plt.text(midpoint, plt.ylim()[1]*0.15, 'RI period', color='grey',
             ha='center', va='top', fontsize=10, fontstyle='italic')
    # Add colorbar
    cbar = plt.colorbar(sc, orientation='vertical')  # Associate the colorbar with the scatter plot
    # cbar.set_label('Degrees from gauge (°)')
    cbar.set_label('Distance from centre (km)')
    # from sklearn.metrics import mean_squared_error
    # rmse = np.sqrt(mean_squared_error(cyg_winds, ref_winds))
    # plt.title(f"{name} – CYGNSS vs IBTrACS Wind - RMSE: {rmse:.2f} m/s")

    # cbar.ax.invert_yaxis()  # Flip the colorbar

    # Save the figure
    directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\Timeseries\{1}_{2}.png'.format(
        comp, cyg_version, name)
    plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()


def plot_timeline(cyg_files_list,start_time,end_time,times_array, winds_array,country):
    # extract lat and lon
    if data_folder == '\PSLGM':
        gauge_lat,gauge_lon = location_dict[country]['lat'],location_dict[country]['lon']
    elif data_folder == '\BOM Stations' or data_folder == '\Aus 1hr':
        gauge_lat,gauge_lon = float(country[0]),float(country[1])
    
    # arrays needed for final results
    cyg_times = np.ma.MaskedArray([])
    cyg_winds = np.ma.MaskedArray([])
    dist_cyg_stat = np.ma.MaskedArray([])
    
    # Extract date components from the timestamps
    date_start_storm = start_time.date()
    date_end_storm = end_time.date()

    # Calculate the number of days between the dates
    num_days_storm = (date_end_storm - date_start_storm).days


    # Generate a list of dates between the start and end timestamps
    date_list_storm = [date_start_storm + timedelta(days=i) for i in range(num_days_storm + 1)]



    for date in date_list_storm:
        #find and open the corresponding cygnss file
        try: 
            if cyg_version == 'storm_centric':
                valid_place = np.where([date == x for x in cyg_file_dates])[0]
                if len(valid_place)==0: 
                    raise AssertionError
                cyg_file_path = cyg_files_list[cyg_files_num[valid_place[0]]]
                cyg_nc = nc.Dataset(cyg_file_path)
                cygnss_lons = cyg_nc.variables['lon'][:] - 180 ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ;
                cygnss_time =cyg_nc.variables['time'][:] ; cygnss_time_offset = cyg_nc.variables['time_offset'][:]
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = [0,0]
                
            else:
                cyg_date_ind = np.where([date == x for x in cyg_file_dates])[0]
                if len(cyg_date_ind)==0: 
                    raise AssertionError
                cyg_nc = nc.Dataset(cyg_files_list[cyg_date_ind[0]])
                cygnss_lons = cyg_nc.variables['lon'][:] ; cygnss_lats = cyg_nc.variables['lat'][:]; cygnss_wind= cyg_nc.variables['wind_speed'][:] ; 
                cygnss_qual = cyg_nc.variables['sample_flags']
                cygnss_time =cyg_nc.variables['sample_time'][:] ; cyg_start_time_str = cyg_nc.variables['sample_time'].units[14:33]
        
                if cyg_version== 'l2_3.2':
                    cygnss_wind= cyg_nc.variables['preliminary_yslf_wind_speed'][:] 
                    inc = cyg_nc.variables['incidence_angle']; rcg = cyg_nc.variables['range_corr_gain'][:]; uncertainty= cyg_nc.variables['preliminary_yslf_wind_speed_uncertainty']; sc_num=cyg_nc.variables['spacecraft_num']
                    cygnss_qual = cyg_nc.variables['fds_sample_flags']
                
                lats_boolean = (cygnss_lats.data > AOI_lat_min) & (cygnss_lats.data < AOI_lat_max)
                lons_boolean = (cygnss_lons.data > AOI_lon_min) & (cygnss_lons.data < AOI_lon_max)
                box_indices = np.where(lats_boolean & lons_boolean)[0]
            
            if len(box_indices)>0:
                if cyg_version == 'storm_centric':
                    # find closest epoch to SAR time
                    date_format = "%y-%m-%d %H:%M:%S"
                    cyg_start_time_str = cyg_nc.variables['time'].units[14:31]
                    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                    cyg_measured_times = np.ma.MaskedArray(np.zeros(np.shape(cygnss_time_offset)))
                    for epoch in range(len(cygnss_time)):
                        cyg_measured_times[epoch] = cygnss_time_offset[epoch] + cygnss_time[epoch]
                    # Focus on the lat and lons
                    timedeltas = np.vectorize(lambda x: timedelta(hours=x))(cyg_measured_times[:,lats_boolean,:][:,:,lons_boolean])
                    # Add the base datetime to each timedelta to get the 3D array of datetime objects
                    datetime_array = np.vectorize(lambda x: cyg_start_date + x)(timedeltas)                        
                    
                else:
                    date_format = "%Y-%m-%d %H:%M:%S"
                    cyg_start_date = datetime.strptime(cyg_start_time_str, date_format)
                    datetime_array = np.vectorize(lambda x: cyg_start_date + timedelta(seconds=x))(cygnss_time[box_indices])
                    
                time_boolean = (datetime_array > start_time) & (datetime_array < end_time)
                true_indices = np.where(time_boolean)[0]
                if len(true_indices)>0:
                    if cyg_version == 'storm_centric':
                        cygnss_wind_AOI = cygnss_wind[:,lats_boolean,:][:,:,lons_boolean][time_boolean] ; cygnss_time_AOI = datetime_array[time_boolean]
                        cyg_times=np.append(cyg_times,cygnss_time_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI) 
                        cygnss_lats_AOI = cygnss_lats[lats_boolean][np.where(time_boolean)[1]] ; cygnss_lons_AOI = cygnss_lons[lons_boolean][np.where(time_boolean)[2]]
                    else:
                        cygnss_wind_AOI = cygnss_wind[box_indices[true_indices]] ; cygnss_time_AOI = datetime_array[true_indices] 
                        cyg_times=np.append(cyg_times,cygnss_time_AOI) ; cyg_winds = np.append(cyg_winds, cygnss_wind_AOI) 
                        cygnss_lats_AOI = cygnss_lats[box_indices[true_indices]] ; cygnss_lons_AOI = cygnss_lons[box_indices[true_indices]]

                    # Initialize a Geod object (WGS84 ellipsoid is the default)
                    from cartopy.geodesic import Geodesic
                    
                    # Calculation using  a direct time methods of a single geo distances %}
                    # Create a Geodesic object
                    geod = Geodesic()
                    
                    # Create an array of (lon, lat) pairs for the input points
                    points = np.array(list(zip(cygnss_lons_AOI, cygnss_lats_AOI)))
                    
                    # Reference point as a single (lon, lat) pair
                    ref_point = np.array([(gauge_lon, gauge_lat)])
                    
                    # Use the geodesic inverse method to compute distances
                    distances = geod.inverse(ref_point, points)[:,0]/1000 # Distance is in the last column

                    dist_cyg_stat = np.append(dist_cyg_stat,distances)  # Convert meters to kilometers

                    
                    # # alternative distance calculation using lat lon point
                    # distance = np.sqrt((cygnss_lats_AOI - gauge_lat)**2 + (cygnss_lons_AOI - gauge_lon)**2)
                    # max_dist = (np.sqrt((search_distance)**2 + (search_distance)**2))
                    # dist_normalised = abs(distance-max_dist) / max_dist
                    # dist_cyg_stat = np.append(dist_cyg_stat, dist_normalised)
                    
                    # sc_num_AOI = sc_num[box_indices[true_indices]]; inc_AOI = inc[box_indices[true_indices]] ; rcg_AOI = rcg[box_indices[true_indices]] ; uncertainty_AOI = uncertainty[box_indices[true_indices]]
        except AssertionError:
            pass
    # remove invalid values 
    cyg_times = cyg_times[cyg_winds>0]
    cyg_winds = cyg_winds[cyg_winds>0]
    dist_cyg_stat = dist_cyg_stat[cyg_winds>0]
    
    counter = 0
    if len(cyg_times) >0:     
        # plot        
        # plt.figure(figsize=(10, 5))
        # plt.scatter(cyg_times, cyg_winds, label='CYGNSS', color='red', alpha = dist_cyg_stat)
        # plt.plot(times_array, winds_array, label='Station', color='blue')
        # plt.title(str(country))
        # plt.xlabel('Time')
        # plt.ylabel('Wind speed (m/s)')
        # cbar = plt.colorbar(sm, orientation='vertical')
        # cbar.set_label('Alpha Value (dist_cyg_stat)')
        # plt.legend()
        # directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\Timeseries\{1}{2}{3}.png'.format(comp,cyg_version,str(country[0:3]),str(counter))
        # plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        # plt.show()
        
                
        # Normalize dist_cyg_stat to be between 0 and 1 for color mapping
        norm = mcolors.Normalize(vmin=min(dist_cyg_stat), vmax=max(dist_cyg_stat))
        cmap = cm.Reds_r  # Choose a colormap
        
        # Start plotting
        plt.figure(figsize=(10, 5))
        
        # Use `c` instead of `alpha` for color mapping
        sc = plt.scatter(cyg_times, cyg_winds, label='CYGNSS', c=dist_cyg_stat, cmap=cmap, alpha=0.7)
        plt.plot(times_array, winds_array, label='Station', color='blue')
        plt.xticks(rotation=45)

        plt.title(str(country))
        plt.xlabel('Time')
        plt.ylabel('Wind speed (m/s)')
        plt.legend()
        
        # Add colorbar
        cbar = plt.colorbar(sc, orientation='vertical')  # Associate the colorbar with the scatter plot
        # cbar.set_label('Degrees from gauge (°)')
        cbar.set_label('Distance from gauge (km)')
        cbar.ax.invert_yaxis()  # Flip the colorbar

        # Save the figure
        directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\Timeseries\{1}{2}{3}.png'.format(
            comp, cyg_version, str(country[0:3]), str(counter)
        )
        plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
        plt.show()
        counter += 1
        
        




for country,dates,wind_gust_before,wind_gust_after in events_list:
    if data_folder == '\PSLGM':
        gauge_lat,gauge_lon = location_dict[country]['lat'],location_dict[country]['lon']
    elif data_folder == '\BOM Stations' or data_folder == '\Aus 1hr':
        gauge_lat,gauge_lon = float(country[0]),float(country[1])
        
    AOI_lat_min,AOI_lat_max = gauge_lat-search_distance,gauge_lat+search_distance
    AOI_lon_min,AOI_lon_max = gauge_lon-search_distance,gauge_lon+search_distance # bounding box around buoy location
    
    # AOI_time_start_storm,AOI_time_end_storm = dates[0][0]-timedelta(hours=1), dates[-1][0]+timedelta(hours=1)
    # AOI_time_start_before,AOI_time_end_before = wind_gust_before[-1][0]-timedelta(hours=cyg_time_thresh), wind_gust_before[0][0]
    # AOI_time_start_after,AOI_time_end_after = wind_gust_after[0][0], wind_gust_after[-1][0]+timedelta(hours=cyg_time_thresh)
    AOI_time_start_before = wind_gust_before[-1][0]-timedelta(hours=cyg_time_thresh)
    AOI_time_end_after = wind_gust_after[-1][0]+timedelta(hours=cyg_time_thresh)
    

    # Convert data to numpy arrays for vectorized operations
    timestamps_storm, wind_speeds_storm = zip(*dates)
    times_array_storm = np.array(timestamps_storm)
    winds_array_storm = np.array(wind_speeds_storm)    
    
    timestamps_before, wind_speeds_before = zip(*wind_gust_before)
    times_array_before = np.array(timestamps_before)
    winds_array_before = np.array(wind_speeds_before)
    
    timestamps_after, wind_speeds_after = zip(*wind_gust_after)
    times_array_after = np.array(timestamps_after)
    winds_array_after = np.array(wind_speeds_after)
    
    # combine the 3 segments, dates, before, and after. 
    times_array = np.concatenate((np.flip(times_array_before),times_array_storm,times_array_after))
    winds_array = np.concatenate((np.flip(winds_array_before),winds_array_storm,winds_array_after))
    
    plot_timeline(cyg_files_list, AOI_time_start_before, AOI_time_end_after,times_array,winds_array,country)