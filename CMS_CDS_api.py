# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 10:32:03 2025

# CMS api downloader

@author: ashle
"""

# import copernicusmarine
# comp = 'ashle'

# datasets=['cfo','c3','h2b','h2c','j3','al','s3a','s3b','s6a','swon']
# datasets=['cfo']
# for data in datasets:
#     copernicusmarine.subset(
#       dataset_id='cmems_obs-wave_glo_phy-swh_nrt_'+data+'-l3_PT1S',
#       variables=["VAVH"],
#       username='asiedlecki',
#       password='XsySq5LJ8TmZKs',
#       minimum_longitude=140,
#       maximum_longitude=147,
#       minimum_latitude=26,
#       maximum_latitude=29,
#       start_datetime="2022-08-28T00:00:00",
#       end_datetime="2022-08-29T23:59:59",
#       output_filename = data+"_swh.nc",
#       output_directory = r'C:\Users\\'+comp+'\OneDrive - RMIT University\PHD\Data\Casestudy\copernicus-data'
#     )

import cdsapi

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_temperature",
        "mean_sea_level_pressure",
        "mean_wave_period",
        "significant_height_of_combined_wind_waves_and_swell",
        "surface_pressure",
        "total_precipitation",
        "ocean_surface_stress_equivalent_10m_neutral_wind_speed",
        "significant_height_of_wind_waves",
        "maximum_total_precipitation_rate_since_previous_post_processing"
    ],
    "year": ["2022"],
    "month": ["08"],
    "day": ["28", "29", "30"],
    "time": [
        "00:00", "01:00", "02:00",
        "03:00", "04:00", "05:00",
        "06:00", "07:00", "08:00",
        "09:00", "10:00", "11:00",
        "12:00", "13:00", "14:00",
        "15:00", "16:00", "17:00",
        "18:00", "19:00", "20:00",
        "21:00", "22:00", "23:00"
    ],
    "data_format": "grib",
    "download_format": "unarchived",
    "area": [28, 142, 26, 147]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()


client = cdsapi.Client()
client.retrieve(dataset, request).download()
