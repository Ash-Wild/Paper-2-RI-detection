def cyg_l1_downloader(minimum_longitude,  maximum_longitude,  minimum_latitude,  maximum_latitude,start_datetime, end_datetime,storm_id,overide=False):
    import earthaccess
    import xarray as xr
    EARTHDATA_USERNAME ='asiedlecki'
    EARTHDATA_PASSWORD =r'RAPK8AKnA%Ghz0F1'
    
    # Login first, before search
    auth = earthaccess.login(strategy="netrc")
    # auth = earthaccess.login(EARTHDATA_USERNAME, EARTHDATA_PASSWORD)
    if not auth.authenticated:
        raise RuntimeError("Failed to authenticate with Earthdata. Check credentials.")
    
# Raw IF
#     results = earthaccess.search_data(
#         daac='PODAAC',
#         short_name="CYGNSS_L1_RAW_IFC2A0", # CYGNSS_L1_CAL_RAW_IF_V1.0
#     bounding_box=(minimum_longitude, minimum_latitude, maximum_longitude, maximum_latitude),  
#         # Time period: One week in January 2024 (times are in UTC)
#     temporal=(start_datetime, end_datetime),
#     count=-1,
#     )

# # Full DDM
#     results = earthaccess.search_data(
#         daac='PODAAC',
#         short_name="CYGNSS_L1_FULL_DDM_V3.0C2A0",
#     bounding_box=(minimum_longitude, minimum_latitude, maximum_longitude, maximum_latitude),  
#         # Time period: One week in January 2024 (times are in UTC)
#     temporal=(start_datetime, end_datetime),
#     count=-1,
#     )
# L1
    results = earthaccess.search_data(
        short_name="CYGNSS_L1_V3.2",
        bounding_box=(minimum_longitude, minimum_latitude, maximum_longitude, maximum_latitude),  
        # Time period: One week in January 2024 (times are in UTC)
        temporal=(start_datetime, end_datetime),
        # count=-1,
        )
    
    # Download data
    earthaccess.download(granules=results, local_path=f"E:\\Phd_data\\cyg_l1\\{storm_id}", show_progress=True)

    

    # files = earthaccess.open(results) # 
    # file_p = files[0]
    # refl = xr.open_dataset(file_p)
    # wvl = xr.open_dataset(file_p, group="sensor_band_parameters")
    # loc = xr.open_dataset(file_p, group="location")
    # ds = xr.merge([refl, loc])
    # ds = ds.assign_coords(
    #     {
    #         "downtrack": (["downtrack"], refl.downtrack.data),
    #         "crosstrack": (["crosstrack"], refl.crosstrack.data),
    #         **wvl.variables,
    #     }
    # )
    # ds = xr.open_dataset(files[0], engine="h5netcdf", backend_kwargs={"use_cftime": True})
    # ds

    
    
    # for result in results:
    #     vds =  earthaccess.open_virtual_dataset(result, access="indirect")
    #     vds.load()  # Load the dataset into memory
    #     vds

    # vds = earthaccess.open_virtual_mfdataset(results, access="indirect", load=False, concat_dim="time", coords="minimal", compat="override", combine_attrs="drop_conflicts")
    

...