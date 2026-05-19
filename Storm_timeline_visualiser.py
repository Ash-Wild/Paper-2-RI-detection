"""
Created on Tue Jul 23 17:48:38 2024
Timeline visualiser
Input: times and values of station, and CYGNSS files
output: map showing timeline shift
@author: Ashley
"""
def plot_distribution(ref_df,noaa_df, merged_df, l2_df):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    comp = 'ashle'
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
    datasets = [[ref_df, 'IBTrACS', 'blue'],[noaa_df, 'NOAA', 'red'],[merged_df,'Merged','orange'],[l2_df,'L2 3.2', 'green']]

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
        lower_ws = grouped['wind_speed'].quantile(0.1)
        upper_ws = grouped['wind_speed'].quantile(0.9)
        
        if label=='Merged':
            # Interpolate stats to full bin range
            mean_ws = mean_ws.reindex(labels_bin).interpolate(method='linear')
            lower_ws = lower_ws.reindex(labels_bin).interpolate(method='linear')
            upper_ws = upper_ws.reindex(labels_bin).interpolate(method='linear')

        # Plot
        plt.plot(labels_bin, mean_ws, label=f'{label} mean', color=color)
        plt.fill_between(labels_bin, lower_ws, upper_ws, color=color, alpha=0.2, label=f'{label} 90% range')

    # Set consistent x and y limits
    plt.xlim(-12, 36)
    plt.ylim(0, 69)

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
    plt.show(block=False)
    plt.pause(1)
    plt.close()

def plot_timeline_both(cyg_wind,cyg_time,cyg_distances,merged_wind, merged_times, ref_times, ref_winds,name,ri_start,ri_end,cyg_sources=None):
    import matplotlib.pyplot as plt
    import numpy as np
    comp = 'ashle'

    # Start plotting
    plt.figure(figsize=(10, 5))

    # Plot NOAA-like points as separate colors per source so all three datasets
    # are visually distinct on one figure (noaa_1.2, l2_3.2, and storm_centric).
    cyg_time_arr = np.asarray(cyg_time)
    cyg_wind_arr = np.asarray(cyg_wind)
    if cyg_sources is None:
        cyg_sources_arr = np.full(len(cyg_time_arr), 'noaa_1.2', dtype=object)
    else:
        cyg_sources_arr = np.asarray(cyg_sources, dtype=object)
        if len(cyg_sources_arr) != len(cyg_time_arr):
            cyg_sources_arr = np.full(len(cyg_time_arr), 'noaa_1.2', dtype=object)

    source_style = [
        ('l2_3.2', 'YSLF', '#1a9850', 'o'),  # circle
        ('noaa_1.2', 'NOAA', '#d73027', 'o'),  # square
    ]
    for source_key, source_label, source_color, marker_style in source_style:
        mask = (cyg_sources_arr == source_key)
        if np.any(mask):
            plt.scatter(cyg_time_arr[mask], cyg_wind_arr[mask], label=source_label,
                       color=source_color, alpha=0.6, marker=marker_style)

    # Plot storm-centric (merged) as the third dataset color.
    plt.scatter(merged_times, merged_wind, label='Merged', color="#bbbe0a", alpha=1, marker='s')
    plt.plot(ref_times, ref_winds, label='IBTrACS', color='blue')
    plt.xticks(rotation=45)
    plt.xlabel('Time')
    plt.ylabel('Wind speed (m/s)')
    plt.legend()
    # Add grey dashed vertical lines
    plt.axvline(ri_start, color='grey', linestyle='--')
    plt.axvline(ri_end, color='grey', linestyle='--')
    
    # Set consistent x and y limits
    plt.ylim(0, 70)

    # Add label in between
    midpoint = ri_start + (ri_end - ri_start) / 2
    plt.text(midpoint, plt.ylim()[1]*0.95, 'RI period', color='grey',
             ha='center', va='top', fontsize=10, fontstyle='italic')
    # from sklearn.metrics import mean_squared_error
    # rmse = np.sqrt(mean_squared_error(cyg_winds, ref_winds))
    # plt.title(f"{name} – CYGNSS vs IBTrACS Wind - RMSE: {rmse:.2f} m/s")

    # cbar.ax.invert_yaxis()  # Flip the colorbar

    # Save the figure
    directory = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\Timeseries\both_{1}.png'.format(
        comp, name)
    plt.savefig(directory, format='png', dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show(block=False)
    plt.pause(1)  # Brief pause to display
    plt.close()  # Properly close the figure to prevent hanging
                    

def plot_tc_animation(tc_lat_interpolated, tc_lon_interpolated, tc_time_interpolated, tc_vmax_interpolated,
                      cyg_lats, cyg_lons, cyg_times, cyg_winds, cyg_sources,
                      cyg_lats_all, cyg_lons_all, cyg_times_all, cyg_sources_all,
                      merged_vmax_lat=None, merged_vmax_lon=None,
                      merged_datetimes=None, merged_vmax_winds=None,
                      tc_storm_id='', tc_durations=0, tc_intensities=0, tc_event=0,
                      ri_start_time=None, ri_end_time=None, mid_time=None,
                      comp='ashle', MAX_DISTANCE=80):
    """
    Create 4 animations showing TC track with CYGNSS measurements:
    - NOAA 1.2 only
    - L2 3.2 YSLF only
    - Merged measurements only
    - Combined (all three products together)

    Parameters
    ----------
    tc_lat_interpolated, tc_lon_interpolated, tc_time_interpolated, tc_vmax_interpolated : arrays
        Interpolated track data
    cyg_lats, cyg_lons, cyg_times, cyg_winds, cyg_sources : arrays
        CYGNSS measurements within RI window
    cyg_lats_all, cyg_lons_all, cyg_times_all, cyg_sources_all : arrays
        All CYGNSS measurements before filtering
    merged_vmax_lat, merged_vmax_lon : arrays, optional
        Latitude/longitude of maximum wind for merged measurements
    merged_datetimes, merged_vmax_winds : arrays, optional
        Datetime and wind speed for merged measurements
    """
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import warnings
    from matplotlib.patches import Circle
    import numpy as np
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    dlat = np.diff(tc_lat_interpolated, prepend=tc_lat_interpolated[0])
    dlon = np.diff(tc_lon_interpolated, prepend=tc_lon_interpolated[0])

    # Filter out zero vectors to avoid divide by zero in quiver
    magnitude = np.sqrt(dlon**2 + dlat**2)
    zero_mask = magnitude < 1e-6
    dlon_safe = np.where(zero_mask, 0.01, dlon)
    dlat_safe = np.where(zero_mask, 0.01, dlat)

    buffer = 1.1
    latitudes = tc_lat_interpolated
    longitudes = tc_lon_interpolated
    times = tc_time_interpolated
    winds = tc_vmax_interpolated

    max_distance_deg = MAX_DISTANCE / 111.0
    inner_ring_deg = max(max_distance_deg * (0.03 / 0.99), 1e-4)

    # Create colorbar setup
    norm = plt.Normalize(vmin=15, vmax=45)
    cmap = plt.cm.coolwarm

    # Animation source configurations: (name, show_noaa, show_yslf, show_merged)
    animation_configs = [
        ('NOAA', True, False, False),
        ('YSLF', False, True, False),
        ('Merged', False, False, True),
        ('Combined', True, True, True),
    ]

    # Create animations for each configuration
    for anim_name, show_noaa, show_yslf, show_merged in animation_configs:
        fig, ax = plt.subplots(figsize=(8, 6))

        # Set plot limits
        ax.set_xlim(min(longitudes) - buffer, max(longitudes) + buffer)
        ax.set_ylim(min(latitudes) - buffer, max(latitudes) + buffer)

        # Create colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', label='Wind Speed (m/s)')

        def update_frame(frame,
                        show_noaa_data=show_noaa,
                        show_yslf_data=show_yslf,
                        show_merged_data=show_merged):
            ax.clear()

            # Reset plot limits
            ax.set_xlim(min(longitudes) - buffer, max(longitudes) + buffer)
            ax.set_ylim(min(latitudes) - buffer, max(latitudes) + buffer)

            # Plot arrows for track movement
            quiver = ax.quiver(longitudes[:frame + 1], latitudes[:frame + 1],
                                dlon_safe[:frame + 1], dlat_safe[:frame + 1],
                                winds[:frame + 1],
                                cmap='coolwarm', width=0.005, norm=norm, angles='xy', scale_units='xy')

            # Track if we've added any measurements for legend
            has_measurements = False

            # Plot NOAA 1.2 measurements as squares
            if show_noaa_data:
                mask_noaa = (cyg_times <= times[frame]) & (cyg_times >= times[frame] - np.timedelta64(8,'h')) & (cyg_sources == 'noaa_1.2')
                if np.any(mask_noaa):
                    ax.scatter(cyg_lons[mask_noaa], cyg_lats[mask_noaa],
                                c=cyg_winds[mask_noaa], cmap=cmap, norm=norm,
                                edgecolor='black', s=60, marker='s', label='NOAA 1.2')
                    has_measurements = True

            # Plot L2 3.2 YSLF measurements as circles
            if show_yslf_data:
                mask_yslf = (cyg_times <= times[frame]) & (cyg_times >= times[frame] - np.timedelta64(8,'h')) & (cyg_sources == 'l2_3.2')
                if np.any(mask_yslf):
                    ax.scatter(cyg_lons[mask_yslf], cyg_lats[mask_yslf],
                                c=cyg_winds[mask_yslf], cmap=cmap, norm=norm,
                                edgecolor='black', s=40, marker='o', label='YSLF (L2 3.2)')
                    has_measurements = True

            # Plot merged measurements (as diamonds)
            if show_merged_data and merged_vmax_lat is not None and merged_vmax_lon is not None and merged_datetimes is not None and merged_vmax_winds is not None:
                mask_merged = (merged_datetimes <= times[frame]) & (merged_datetimes >= times[frame] - np.timedelta64(8,'h'))
                if np.any(mask_merged):
                    if len(merged_vmax_lon) == len(merged_datetimes):
                    # Plot merged measurements at their vmax lat/lon locations with wind speed coloring
                        ax.scatter(merged_vmax_lon[mask_merged], merged_vmax_lat[mask_merged],
                                c=merged_vmax_winds[mask_merged], cmap=cmap, norm=norm,
                                edgecolor='black', s=70, marker='d', label='Merged')
                    else:
                        # lat and lon are a grid and want to show the epoch of winds at the time in mask
                        ax.scatter(merged_vmax_lon, merged_vmax_lat,
                                c=merged_vmax_winds[mask_merged], cmap=cmap, norm=norm,
                                edgecolor='black', s=20, marker='s', label='Merged')
                    has_measurements = True

            ax.legend(loc='upper left', fontsize=8)

            # Update phase label
            if times[frame] < ri_start_time:
                phase = 'Before RI'
            elif times[frame] < mid_time:
                phase = 'Early RI'
            elif times[frame] <= ri_end_time:
                phase = 'Later RI'
            elif times[frame] > ri_end_time:
                phase = 'After RI'
            else:
                phase = 'Unknown'

            # Add rings
            ring = Circle((longitudes[frame], latitudes[frame]), max_distance_deg,
                         fill=False, color='black', linestyle='--', linewidth=1)
            ax.add_patch(ring)
            ring = Circle((longitudes[frame], latitudes[frame]), inner_ring_deg,
                         fill=False, color='black', linestyle='-', linewidth=0.5)
            ax.add_patch(ring)

            ax.text(0.98, 0.95, f"Phase: {phase}, Time: {times[frame]}",
                    transform=ax.transAxes, fontsize=12, ha='right', color='black',
                    bbox=dict(facecolor='white', alpha=0.6))

            # Titles & labels
            ax.set_title('ID: {0}, Duration: {1}, Intensity: {2}\n{3}'.format(
                tc_storm_id, tc_durations, tc_intensities, anim_name))
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')

        # Create animation
        pause_frames = 20
        ani = animation.FuncAnimation(fig, update_frame,
                                     frames=list(range(len(times))) + [len(times) - 1] * pause_frames,
                                     repeat=True)

        # Save GIF
        save_dir = r'C:\Users\{0}\OneDrive - RMIT University\PHD\Plots\gifs'.format(comp)
        import os
        os.makedirs(save_dir, exist_ok=True)
        ani.save(os.path.join(save_dir, f'{anim_name}_{tc_event}_{tc_storm_id}.gif'), fps=5)

        # Show plot (non-blocking)
        plt.show(block=False)
        plt.pause(1)
        plt.close(fig)


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
        
    