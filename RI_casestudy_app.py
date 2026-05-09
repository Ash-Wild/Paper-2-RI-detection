import numpy as np
import pandas as pd


def build_measurements_dataframe(
    cyg_times,
    cyg_lats,
    cyg_lons,
    nearest_tc_coords,
    quadrants,
    meas_data,
    ri_start_time,
    ri_end_time,
    cms_res_swh=None,
    era5_res=None,
    mswep_res=None,
):
    df = pd.DataFrame(
        {
            "time": pd.to_datetime(np.asarray(cyg_times)),
            "lat": np.asarray(cyg_lats, dtype=float),
            "lon": np.asarray(cyg_lons, dtype=float),
            "track_idx": np.asarray(nearest_tc_coords, dtype=int),
            "quadrant": np.asarray(quadrants, dtype=object),
        }
    )

    for key, value in meas_data.items():
        arr = np.asarray(value)
        if arr.ndim == 1 and len(arr) == len(df):
            df[key] = arr

    if cms_res_swh:
        if "VHM0" in cms_res_swh:
            swh_vals = np.asarray(cms_res_swh["VHM0"], dtype=float)
            if len(swh_vals) == len(df):
                df["cms_swh"] = swh_vals
    if era5_res:
        if "swh" in era5_res:
            swh_vals = np.asarray(era5_res["swh"], dtype=float)
            if len(swh_vals) == len(df):
                df["era5_swh"] = swh_vals
        elif "sp" in era5_res:
            sp_vals = np.asarray(era5_res["sp"], dtype=float)
            if len(sp_vals) == len(df):
                df["era5_sp"] = sp_vals
    if mswep_res is not None:
        mswep_arr = np.asarray(mswep_res)
        if mswep_arr.ndim == 1 and len(mswep_arr) == len(df):
            df["mswep_precip"] = mswep_arr

    ri_start = pd.to_datetime(ri_start_time)
    ri_end = pd.to_datetime(ri_end_time)
    df["phase"] = np.where(
        df["time"] < ri_start,
        "Before RI",
        np.where(df["time"] <= ri_end, "During RI", "After RI"),
    )

    df["meas_id"] = np.arange(len(df), dtype=int)
    return df


def build_track_dataframe(tc_time_interpolated, tc_lat_interpolated, tc_lon_interpolated, tc_vmax_interpolated):
    return pd.DataFrame(
        {
            "time": pd.to_datetime(np.asarray(tc_time_interpolated)),
            "lat": np.asarray(tc_lat_interpolated, dtype=float),
            "lon": np.asarray(tc_lon_interpolated, dtype=float),
            "vmax": np.asarray(tc_vmax_interpolated, dtype=float),
        }
    )


def launch_case_study_app(tc_df, meas_df, brcs=None, eff_scatter=None, cms_buoy_data=None, environmental_grids=None, port=8050):
    import dash
    from dash import Dash, Input, Output, dcc, html
    import plotly.graph_objects as go

    tc_df = tc_df.copy()
    meas_df = meas_df.copy()

    tc_df["time"] = pd.to_datetime(tc_df["time"])
    meas_df["time"] = pd.to_datetime(meas_df["time"])

    if "meas_id" not in meas_df.columns:
        meas_df["meas_id"] = np.arange(len(meas_df), dtype=int)

    if "phase" not in meas_df.columns:
        meas_df["phase"] = "Unknown"
    if "quadrant" not in meas_df.columns:
        meas_df["quadrant"] = "Unknown"

    # Build layer options dynamically
    layer_to_column = {"None": None}
    layer_to_grid = {}

    if environmental_grids is None:
        environmental_grids = {}

    for layer_name, layer_spec in environmental_grids.items():
        if not isinstance(layer_spec, dict):
            continue
        if not all(k in layer_spec for k in ["lat", "lon", "data"]):
            continue
        layer_to_grid[layer_name] = layer_spec
    if "ddm_snr" in meas_df.columns:
        layer_to_column["DDM SNR"] = "ddm_snr"
    if "sp_inc_angle" in meas_df.columns:
        layer_to_column["Incidence Angle"] = "sp_inc_angle"
    if "gps_eirp" in meas_df.columns:
        layer_to_column["GPS EIRP"] = "gps_eirp"
    if "cms_swh" in meas_df.columns:
        layer_to_column["CMS SWH"] = "cms_swh"
    if "mswep_precip" in meas_df.columns:
        layer_to_column["MSWEP Precip"] = "mswep_precip"
    if "msl" in meas_df.columns:
        layer_to_column["MSLP"] = "msl"
    elif "era5_sp" in meas_df.columns:
        layer_to_column["MSLP"] = "era5_sp"

    n_times = len(meas_df["time"].sort_values().unique())
    sorted_times = np.array(sorted(meas_df["time"].unique()))

    def _add_background_field(fig, source_df, value_col, layer_title):
        vals = pd.to_numeric(source_df[value_col], errors="coerce")
        valid = source_df[np.isfinite(vals.values)].copy()
        if len(valid) < 8:
            return

        valid_vals = pd.to_numeric(valid[value_col], errors="coerce").to_numpy(dtype=float)
        valid = valid[np.isfinite(valid_vals)]
        valid_vals = valid_vals[np.isfinite(valid_vals)]
        if len(valid_vals) < 8:
            return

        n_bins = int(np.clip(np.sqrt(len(valid_vals)), 12, 55))
        lon_edges = np.linspace(valid["lon"].min(), valid["lon"].max(), n_bins + 1)
        lat_edges = np.linspace(valid["lat"].min(), valid["lat"].max(), n_bins + 1)
        if (lon_edges[-1] - lon_edges[0]) <= 0 or (lat_edges[-1] - lat_edges[0]) <= 0:
            return

        count_grid, _, _ = np.histogram2d(valid["lat"], valid["lon"], bins=[lat_edges, lon_edges])
        sum_grid, _, _ = np.histogram2d(valid["lat"], valid["lon"], bins=[lat_edges, lon_edges], weights=valid_vals)
        mean_grid = np.divide(sum_grid, count_grid, out=np.full_like(sum_grid, np.nan), where=count_grid > 0)

        lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
        lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
        lon_mesh, lat_mesh = np.meshgrid(lon_centers, lat_centers)

        mask = np.isfinite(mean_grid)
        if not np.any(mask):
            return

        fig.add_trace(
            go.Scattergeo(
                lon=lon_mesh[mask],
                lat=lat_mesh[mask],
                mode="markers",
                marker={
                    "symbol": "square",
                    "size": 8,
                    "color": mean_grid[mask],
                    "colorscale": "Plasma",
                    "showscale": True,
                    "colorbar": {"title": layer_title, "x": -0.12},
                    "opacity": 0.55,
                    "line": {"width": 0},
                },
                hovertemplate=layer_title + ": %{marker.color:.3f}<extra></extra>",
                name=f"{layer_title} background",
                showlegend=False,
            )
        )

    def _add_gridded_background(fig, layer_name, layer_spec, selected_time):
        lats = np.asarray(layer_spec["lat"], dtype=float)
        lons = np.asarray(layer_spec["lon"], dtype=float)
        data = np.asarray(layer_spec["data"])
        times = layer_spec.get("time", None)
        window = layer_spec.get("window", np.timedelta64(2, "h"))

        if data.ndim == 3 and times is not None and selected_time is not None:
            time_vals = np.asarray(times).astype("datetime64[s]")
            target_t = np.datetime64(pd.Timestamp(selected_time).to_datetime64(), "s")
            mask = (time_vals <= target_t) & (time_vals >= target_t - window)
            if np.any(mask):
                grid2d = np.nanmean(data[mask], axis=0)
            else:
                nearest_idx = np.abs(time_vals - target_t).argmin()
                grid2d = data[nearest_idx]
        elif data.ndim == 3:
            grid2d = np.nanmean(data, axis=0)
        elif data.ndim == 2:
            grid2d = data
        else:
            return

        grid2d = np.asarray(np.ma.filled(grid2d, np.nan), dtype=float)
        if grid2d.shape != (len(lats), len(lons)):
            if grid2d.shape == (len(lons), len(lats)):
                grid2d = grid2d.T
            else:
                return

        lon_mesh, lat_mesh = np.meshgrid(lons, lats)
        valid = np.isfinite(grid2d)
        if not np.any(valid):
            return

        max_points = 15000
        total_points = valid.size
        stride = int(np.ceil(np.sqrt(total_points / max_points))) if total_points > max_points else 1
        if stride > 1:
            valid = valid[::stride, ::stride]
            lon_plot = lon_mesh[::stride, ::stride]
            lat_plot = lat_mesh[::stride, ::stride]
            val_plot = grid2d[::stride, ::stride]
        else:
            lon_plot = lon_mesh
            lat_plot = lat_mesh
            val_plot = grid2d

        valid_points = np.isfinite(val_plot) & valid
        if not np.any(valid_points):
            return

        fig.add_trace(
            go.Scattergeo(
                lon=lon_plot[valid_points],
                lat=lat_plot[valid_points],
                mode="markers",
                marker={
                    "symbol": "square",
                    "size": 7,
                    "color": val_plot[valid_points],
                    "colorscale": "Plasma",
                    "showscale": True,
                    "colorbar": {"title": layer_name, "x": -0.12},
                    "opacity": 0.5,
                    "line": {"width": 0},
                },
                hovertemplate=layer_name + ": %{marker.color:.3f}<extra></extra>",
                name=f"{layer_name} background",
                showlegend=False,
            )
        )

    app = Dash(__name__)

    app.layout = html.Div(
        [
            html.H3("TC Case Study Explorer"),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Background layer"),
                            dcc.Dropdown(
                                id="layer-dropdown",
                                options=[{"label": k, "value": k} for k in ["None", *layer_to_grid.keys(), *[k for k in layer_to_column.keys() if k != "None"]]],
                                value="None",
                                clearable=False,
                            ),
                        ],
                        style={"width": "24%", "display": "inline-block", "paddingRight": "1%"},
                    ),
                    html.Div(
                        [
                            html.Label("RI phase"),
                            dcc.Dropdown(
                                id="phase-dropdown",
                                options=[{"label": p, "value": p} for p in sorted(meas_df["phase"].dropna().unique())],
                                value=sorted(meas_df["phase"].dropna().unique()),
                                multi=True,
                            ),
                        ],
                        style={"width": "24%", "display": "inline-block", "paddingRight": "1%"},
                    ),
                    html.Div(
                        [
                            html.Label("Quadrant"),
                            dcc.Dropdown(
                                id="quadrant-dropdown",
                                options=[{"label": q, "value": q} for q in sorted(meas_df["quadrant"].dropna().unique())],
                                value=sorted(meas_df["quadrant"].dropna().unique()),
                                multi=True,
                            ),
                        ],
                        style={"width": "24%", "display": "inline-block", "paddingRight": "1%"},
                    ),
                    html.Div(
                        [
                            html.Label("Only valid quality"),
                            dcc.Checklist(
                                id="quality-check",
                                options=[{"label": "Use quality_flags % 2 == 0", "value": "valid"}],
                                value=["valid"] if "quality_flags" in meas_df.columns else [],
                            ),
                        ],
                        style={"width": "24%", "display": "inline-block"},
                    ),
                ]
            ),
            html.Div(
                [
                    html.Label("Time (measurements ±2 hours)"),
                    dcc.Slider(
                        id="time-slider",
                        min=0,
                        max=max(n_times - 1, 0),
                        step=1,
                        value=max(n_times - 1, 0) // 2,
                        tooltip={"always_visible": False},
                    ),
                ],
                style={"marginTop": "10px", "marginBottom": "10px"},
            ),
            dcc.Graph(id="track-map", style={"height": "55vh"}),
            html.Div(
                [
                    html.Div(dcc.Graph(id="ddm-brcs"), style={"width": "49%", "display": "inline-block"}),
                    html.Div(dcc.Graph(id="ddm-eff"), style={"width": "49%", "display": "inline-block"}),
                ]
            ),
            html.Div(
                [
                    html.Div(
                        [html.H4("Selected Measurement"), html.Pre(id="point-metadata", style={"whiteSpace": "pre-wrap", "fontSize": "11px"})],
                        style={"width": "49%", "display": "inline-block", "paddingRight": "1%", "verticalAlign": "top"}
                    ),
                    html.Div(
                        [html.H4("CMS Buoy @ Time"), html.Pre(id="buoy-data", style={"whiteSpace": "pre-wrap", "fontSize": "11px"})],
                        style={"width": "49%", "display": "inline-block", "verticalAlign": "top", "paddingLeft": "1%"}
                    ),
                ],
                style={"marginTop": "10px"}
            ),
        ],
        style={"padding": "12px"},
    )

    @app.callback(
        Output("track-map", "figure"),
        Input("layer-dropdown", "value"),
        Input("phase-dropdown", "value"),
        Input("quadrant-dropdown", "value"),
        Input("time-slider", "value"),
        Input("quality-check", "value"),
    )
    def update_map(layer_name, phases, quadrants, time_idx, quality_check):
        filtered = meas_df.copy()
        selected_time = None

        if phases:
            filtered = filtered[filtered["phase"].isin(phases)]
        if quadrants:
            filtered = filtered[filtered["quadrant"].isin(quadrants)]

        # Filter for ±2 hours around selected time
        if len(sorted_times) > 0:
            selected_time = sorted_times[int(time_idx)]
            time_window = pd.Timedelta(hours=2)
            filtered = filtered[(filtered["time"] >= selected_time - time_window) & (filtered["time"] <= selected_time + time_window)]

        if quality_check and "valid" in quality_check and "quality_flags" in filtered.columns:
            filtered = filtered[(filtered["quality_flags"] % 2) == 0]

        fig = go.Figure()
        
        # Add TC track colored by Vmax with colorbar (line shows connections)
        fig.add_trace(
            go.Scattergeo(
                lon=tc_df["lon"],
                lat=tc_df["lat"],
                mode="lines+markers",
                line={"color": "rgba(200, 200, 200, 0.5)", "width": 2},
                marker={
                    "size": 6,
                    "color": tc_df["vmax"].values,
                    "colorscale": "RdYlBu_r",
                    "showscale": True,
                    "colorbar": {"title": "Vmax (m/s)", "x": 1.05},
                    "line": {"width": 0.5, "color": "white"},
                },
                name="TC track",
                hovertemplate="TC<br>Time: %{text}<br>Lat: %{lat:.2f}<br>Lon: %{lon:.2f}<br>Vmax: %{marker.color:.1f}<extra></extra>",
                text=tc_df["time"].astype(str),
            )
        )
        
        # Add TC position at selected time
        if len(sorted_times) > 0:
            selected_time = sorted_times[int(time_idx)]
            closest_idx = (tc_df["time"] - selected_time).abs().argmin()
            tc_pos = tc_df.iloc[closest_idx]
            fig.add_trace(
                go.Scattergeo(
                    lon=[tc_pos["lon"]],
                    lat=[tc_pos["lat"]],
                    mode="markers",
                    marker={"size": 20, "color": "darkred", "symbol": "star", "line": {"width": 2, "color": "white"}},
                    name="TC @ selected time",
                    hovertemplate="TC Position<br>Time: " + str(tc_pos["time"]) + "<br>Lat: %{lat:.2f}<br>Lon: %{lon:.2f}<br>Vmax: " + f"{tc_pos['vmax']:.1f}" + "<extra></extra>",
                )
            )

        color_col = layer_to_column.get(layer_name)

        if layer_name in layer_to_grid:
            _add_gridded_background(fig, layer_name, layer_to_grid[layer_name], selected_time)
        elif color_col and color_col in filtered.columns:
            _add_background_field(fig, filtered, color_col, layer_name)

        if color_col and color_col in filtered.columns:
            marker = {
                "size": 8,
                "color": filtered[color_col],
                "colorscale": "Plasma",
                "showscale": True,
                "colorbar": {"title": layer_name, "x": -0.15},
                "opacity": 0.9,
                "line": {"width": 0.5, "color": "white"},
            }
        else:
            marker = {"size": 8, "color": "steelblue", "opacity": 0.8, "line": {"width": 0.5, "color": "white"}}

        fig.add_trace(
            go.Scattergeo(
                lon=filtered["lon"],
                lat=filtered["lat"],
                mode="markers",
                marker=marker,
                name="Measurements (click to inspect)",
                customdata=filtered["meas_id"].values,
                hovertemplate=(
                    "Meas %{customdata}<br>Time: %{text}<br>Lat: %{lat:.3f}°<br>Lon: %{lon:.3f}°<br><b>Click to see DDM</b><extra></extra>"
                ),
                text=filtered["time"].astype(str),
            )
        )

        # Add CMS Buoy data if available
        if cms_buoy_data is not None and not cms_buoy_data.empty:
            buoy_df = pd.DataFrame(cms_buoy_data).copy()
            if 'time' in buoy_df.columns:
                buoy_df['time'] = pd.to_datetime(buoy_df['time'])
            
            if len(buoy_df) > 0:
                fig.add_trace(
                    go.Scattergeo(
                        lon=buoy_df["lon"] if "lon" in buoy_df.columns else [],
                        lat=buoy_df["lat"] if "lat" in buoy_df.columns else [],
                        mode="markers",
                        marker={"size": 15, "color": "gold", "symbol": "diamond", "line": {"width": 2, "color": "darkgoldenrod"}},
                        name="CMS Buoy",
                        hovertemplate="CMS Buoy<br>Lat: %{lat:.3f}°<br>Lon: %{lon:.3f}°<extra></extra>"
                    )
                )

        fig.update_layout(
            title="TC Track + CYGNSS Measurements",
            template="plotly_white",
            uirevision="constant",
            legend={"orientation": "v", "yanchor": "top", "y": 0.99, "xanchor": "left", "x": 0.01},
            geo={
                "projection_type": "mercator",
                "fitbounds": "locations",
                "showland": True,
                "landcolor": "rgb(243, 243, 243)",
                "coastlinecolor": "rgb(100, 100, 100)",
                "showocean": True,
                "oceancolor": "rgb(204, 229, 255)",
                "showcountries": True,
                "countrycolor": "rgb(180, 180, 180)",
            },
        )
        return fig

    @app.callback(
        Output("ddm-brcs", "figure"),
        Output("ddm-eff", "figure"),
        Output("point-metadata", "children"),
        Input("track-map", "clickData"),
    )
    def show_selected_measurement(click_data):
        empty_fig = go.Figure()
        empty_fig.update_layout(template="plotly_white", title="Click a measurement point on the map")
        empty_msg = "Click a measurement point on the map to inspect DDM and metadata."

        if not click_data or not click_data.get("points") or len(click_data.get("points", [])) == 0:
            return empty_fig, empty_fig, empty_msg

        try:
            point = click_data["points"][0]
            # Extract customdata - it can be a scalar, list, or nested list
            customdata = point.get("customdata", None)
            
            if customdata is None:
                return empty_fig, empty_fig, "No customdata in click (try clicking on measurement points only)."
            
            # Handle various formats of customdata
            if isinstance(customdata, (list, tuple)):
                if len(customdata) > 0:
                    selected_id = customdata[0]
                else:
                    return empty_fig, empty_fig, "Empty customdata list."
            else:
                selected_id = customdata
            
            # Ensure it's numeric
            try:
                selected_id = int(selected_id)
            except (ValueError, TypeError):
                return empty_fig, empty_fig, f"Invalid measurement ID: {selected_id}"
                
        except Exception as e:
            return empty_fig, empty_fig, f"Error processing click data: {str(e)}"
        
        if selected_id is None:
            return empty_fig, empty_fig, "No measurement ID found in click payload."

        selected = meas_df[meas_df["meas_id"] == int(selected_id)]
        if len(selected) == 0:
            return empty_fig, empty_fig, f"Measurement {selected_id} not found."

        row = selected.iloc[0]
        brcs_fig = go.Figure()
        eff_fig = go.Figure()

        if brcs is not None and int(selected_id) < len(brcs):
            brcs_fig.add_trace(go.Heatmap(z=np.asarray(brcs[int(selected_id)]), colorscale="Viridis"))
            brcs_fig.update_layout(title=f"BRCS DDM (meas_id={selected_id})", template="plotly_white")
        else:
            brcs_fig = empty_fig

        if eff_scatter is not None and int(selected_id) < len(eff_scatter):
            eff_fig.add_trace(go.Heatmap(z=np.asarray(eff_scatter[int(selected_id)]), colorscale="Cividis"))
            eff_fig.update_layout(title=f"Eff. Scatter DDM (meas_id={selected_id})", template="plotly_white")
        else:
            eff_fig = empty_fig

        preferred_cols = [
            "time",
            "lat",
            "lon",
            "phase",
            "quadrant",
            "ddm_snr",
            "sp_inc_angle",
            "gps_eirp",
            "cms_swh",
            "era5_swh",
            "era5_sp",
            "mswep_precip",
            "quality_flags",
            "track_idx",
        ]
        cols = [c for c in preferred_cols if c in selected.columns]
        metadata_lines = [f"{c}: {row[c]}" for c in cols]
        metadata = "\n".join(metadata_lines)

        return brcs_fig, eff_fig, metadata

    @app.callback(
        Output("buoy-data", "children"),
        Input("time-slider", "value"),
    )
    def show_buoy_data(time_idx):
        if cms_buoy_data is None or cms_buoy_data.empty:
            return "No CMS buoy data available."
        
        try:
            buoy_df = pd.DataFrame(cms_buoy_data).copy()
            if 'time' not in buoy_df.columns:
                return "Buoy data missing 'time' column."
            
            # Convert buoy times to datetime, remove timezone if present
            buoy_df['time'] = pd.to_datetime(buoy_df['time']).dt.tz_localize(None)
            
            if len(sorted_times) == 0:
                return "No measurement times available."
            
            selected_time = pd.Timestamp(sorted_times[int(time_idx)]).tz_localize(None)
            
            # Find closest buoy measurement to selected time
            time_diffs = np.abs(buoy_df['time'] - selected_time)
            closest_idx = time_diffs.idxmin()
            closest_time_diff = time_diffs.iloc[closest_idx]
            
            if closest_time_diff > pd.Timedelta(hours=6):
                return f"No buoy data within 6 hours of selected time.\nClosest: {closest_time_diff.total_seconds() / 3600:.1f} hours away"
            
            buoy_row = buoy_df.iloc[closest_idx]
            
            # Format buoy data display
            buoy_lines = [f"Time: {buoy_row['time']} (Δ: {closest_time_diff.total_seconds() / 3600:.2f} hrs)"]
            
            # Display all numeric columns
            for col in buoy_df.columns:
                if col != 'time' and col in buoy_row.index:
                    value = buoy_row[col]
                    if pd.notna(value):
                        if isinstance(value, (int, np.integer)):
                            buoy_lines.append(f"{col}: {value}")
                        elif isinstance(value, (float, np.floating)):
                            buoy_lines.append(f"{col}: {value:.3f}")
                        else:
                            buoy_lines.append(f"{col}: {value}")
            
            return "\n".join(buoy_lines)
        
        except Exception as e:
            return f"Error processing buoy data: {str(e)}"

    app.run(debug=True, port=port)
