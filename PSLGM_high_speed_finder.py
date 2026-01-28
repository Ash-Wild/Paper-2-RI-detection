
"""
Created on Thu Mar 21 08:13:49 2024
BOM AWS and PSLGM analysis
Inputs: PSLGM Bom files
Outputs: file of extreme events in measurement
@author: ashle
"""
import pandas as pd
import glob
from datetime import timedelta 
import pickle

# input files 
comp = 'Ashley'
csv_folder  = r'C:\Users\\' + comp + '\OneDrive - RMIT University\PHD\Data\BOM stations\PSLGM'
wind_var = 'Wind Gust' # 'Wind Gust'  or 'Wind Speed'
wind_thresh = 17
time_thresh = 3
num_hrs_before_after = 24

# event finder
csv_files_list = glob.glob(csv_folder+ "\*.csv")
events_details = []
for csv in csv_files_list:
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv)
    
    # Convert the date-time column to a pandas DateTimeIndex
    df['Datetime'] = pd.to_datetime(df[' Date & UTC Time'])
    df['Datetime_copy'] = df['Datetime']
    df.set_index('Datetime', inplace=True)
    
    # Filter the DataFrame for wind gust higher than threshold m/s
    df_high_wind = df[df[wind_var] > wind_thresh]
    if len(df_high_wind) > time_thresh:
        hours_in_period = 0
        event_durations = [] # a list of event end times, and the duration before them.
        for i in range(len(df_high_wind)-1):
            if df_high_wind.iloc[i]['Datetime_copy'] == df_high_wind.iloc[i+1]['Datetime_copy'] - timedelta(hours=1):
                hours_in_period += 1
                if i == len(df_high_wind)-2 and hours_in_period >= time_thresh:
                    event_durations.append([df_high_wind.iloc[i+1]['Datetime_copy'],hours_in_period])
            else:
                if hours_in_period >= time_thresh:
                    event_durations.append([df_high_wind.iloc[i]['Datetime_copy'],hours_in_period])
                hours_in_period = 0 # reset timer
                    
    # # Calculate the number of consecutive hours where wind gust exceeds 20 m/s
    # consecutive_hours = df[wind_var].rolling(window=time_thresh, min_periods=time_thresh).sum()
    
    # # Identify time periods with at least time thresh consecutive hours of high wind gust
    # time_periods = consecutive_hours[consecutive_hours >= wind_thresh*time_thresh].index # need to also save the wind and gust
    
    # # now to join if there's sequential periods
    # if len(time_periods) >0:
    #     if len(time_periods) ==1:
    #         event_durations = [[time_periods[0],time_thresh]]
    #     else:
    #         event_durations = []
    #         hours_in_period = time_thresh
    #         for i in range(0,len(time_periods)-1):
    #             if time_periods[i] == time_periods[i+1] - timedelta(hours=1):
    #                 hours_in_period += 1
    #                 if i == len(time_periods)-2:
    #                     event_durations.append([time_periods[i+1],hours_in_period])
    #             else:
    #                 event_durations.append([time_periods[i],hours_in_period])
    #                 hours_in_period = time_thresh
    #                 if i == len(time_periods)-2:
    #                     event_durations.append([time_periods[i+1],hours_in_period])
                
           
        if len(event_durations) > 0:        
            # Filter wind gust data for the identified time periods
            # wind_gust_at_time_periods = df_high_wind[df_high_wind['Datetime'].isin(time_periods)]
            # Iterate over the identified time periods
            for end_time,duration in event_durations:
                wind_gust_before = []
                wind_gust_at_time_periods = []
                wind_gust_after = []
                
                try:
                    for counter in range(num_hrs_before_after): # add in the measurements 24 hours before event
                        time_selected_before = end_time - timedelta(hours=duration+counter+1)
                        wind_gust_value = df.loc[time_selected_before][wind_var]
                        wind_gust_before.append([time_selected_before, wind_gust_value])
                        
                        time_selected_after = end_time + timedelta(hours=counter)
                        wind_gust_value = df.loc[time_selected_after][wind_var]
                        wind_gust_after.append([time_selected_after, wind_gust_value])
                        
                    # iterate for the number of measurements required
                    for counter in range(duration):
                        # check that the measurements are sequential
                        time_selected = end_time - timedelta(hours=duration-counter)
                        # Filter wind gust data for the current time period
                        # wind_gust_data = df_high_wind[df_high_wind['Datetime_copy'] == time_selected]
                        # # Extract timestamp and wind gust value for the current time period
                        # timestamp = wind_gust_data['Datetime_copy'].iloc[0]  # Assuming there's only one row for each timestamp
                        # wind_gust_value = wind_gust_data[wind_var].iloc[0]
                        # get the wind gust, and can get other variables if desired to save into list
                        wind_gust_value = df_high_wind.loc[time_selected][wind_var]
        
                        # Append timestamp and wind gust value as a sublist
                        wind_gust_at_time_periods.append([time_selected, wind_gust_value])
                
                    # add them to the list
                    events_details.append([df.columns[-2],wind_gust_at_time_periods,wind_gust_before,wind_gust_after])
                except KeyError:
                    pass
                    

# Save the list to a file
with open(csv_folder+'\events_list.pkl', 'wb') as f:
    pickle.dump(events_details, f, protocol=pickle.HIGHEST_PROTOCOL)