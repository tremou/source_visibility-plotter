import streamlit as st
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris, AltAz, get_sun, get_body_barycentric, get_body
from astropy.time import Time
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.io import curdoc
from bokeh.models import DatetimeTickFormatter




st.set_page_config(
    page_title="Astronomical Object Visibility Plotter",  # Optional: Set the title that appears in the browser tab
    page_icon="📡",  # Optional: Set the favicon
    layout="wide",  # Set the layout to wide
    initial_sidebar_state="expanded", # Optional: Set the initial state of the sidebar
    menu_items={  # Optional: Customize the menu
        'Get Help': 'https://gitlab.nrao.edu/etremou/visibility-plotter',
        'Report a bug': 'https://gitlab.nrao.edu/etremou/visibility-plotter',
        'About': '# etremou@nrao.edu'
    }
)

# Theming (Light Theme)
st.markdown(
    """
    <style>
    body {
        color: #000000;  /* Main text color */
        background-color: #FFFFFF; /* Main background color */
    }
    .main .block-container {
        max-width: 90%;
    }
    </style>
    """,
    unsafe_allow_html=True,
)





# Default location: Very Large Array (VLA)
default_latitude = 34.3031 * u.deg
default_longitude = -107.2206 * u.deg
default_elevation = 2124 * u.m
vla_location = EarthLocation(lat=default_latitude, lon=default_longitude, height=default_elevation)

# Set up Streamlit app
st.title("Astronomical Object Visibility Plotter")
st.write('etremou@nrao.edu')
# Location input (with default)
st.sidebar.subheader("Location (default: VLA location)")
# Convert Quantity to float for Streamlit input, but keep the Quantity for calculations
default_lat_float = default_latitude.to_value(u.deg) # or default_latitude.value
default_lon_float = default_longitude.to_value(u.deg) # or default_longitude.value
default_elev_float = default_elevation.to_value(u.m) # or default_elevation.value

latitude = st.sidebar.number_input("Latitude (degrees)", value=default_lat_float, step=0.01)
longitude = st.sidebar.number_input("Longitude (degrees)", value=default_lon_float, step=0.01)
elevation = st.sidebar.number_input("Elevation (meters)", value=default_elev_float, step=1.0)
  # Convert to value

# Object input (with name lookup and default)
st.sidebar.subheader("Object")
object_name = st.sidebar.selectbox("Object Name", ["Custom", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]) # Dropdown
if object_name == "Custom":
    custom_name = st.sidebar.text_input("Object Name (e.g., M31, Vega, 3C sources, A-team sources)", value="3C286")
    ra = st.sidebar.number_input("Right Ascension (degrees)", value=0.0, step=0.01)
    dec = st.sidebar.number_input("Declination (degrees)", value=0.0, step=0.01)
elif object_name in ["Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]:
    ra = None #set to None so that the next if statement is skipped.
    dec = None
else:
    ra = None
    dec = None



# Time range input
st.sidebar.subheader("Time Range")
start_date = st.sidebar.date_input("Start Date")
end_date = st.sidebar.date_input("End Date")
time_step = st.sidebar.number_input("Time Step (hours)", value=0.5, step=0.5)

# Calculate times (Improved, same as before)
try:
    start_time = Time(str(start_date))
    end_time = Time(str(end_date))+1
    if start_time >= end_time:
        st.error("Invalid time range: Start date must be before end date.")
        st.stop()
    time_step_days = time_step / 24.0
    times = np.arange(start_time.jd, end_time.jd + time_step_days, time_step_days)
    times = Time(times, format='jd')

except Exception as e:
    st.error(f"Error creating time array: {e}")
    st.stop()


# Observer location (use input or default)
observer_location = EarthLocation(lat=latitude * u.deg, lon=longitude * u.deg, height=elevation * u.m) # Units are essential here

# Object coordinates (with name lookup and error handling)
try:
    if object_name == "Sun":
        object_coordinates = get_sun(times)
    elif object_name == "Moon":
        from astropy.coordinates import get_moon
        object_coordinates = get_moon(times, observer_location)  # Observer location is important for Moon!
    elif object_name in ["Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]:
        object_coordinates = get_body(object_name, times)
    elif object_name == "Custom":
        if ra == 0.0 and dec == 0.0:  # Allow name lookup for custom objects
            object_coordinates = SkyCoord.from_name(custom_name)
        else:
            object_coordinates = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
    else:  # Default to custom if no object name matches
        if ra == 0.0 and dec == 0.0:  # Allow name lookup for custom objects
            object_coordinates = SkyCoord.from_name(custom_name)
        else:
            object_coordinates = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
except Exception as e:
    st.error(f"Error resolving object: {e}")
    st.stop()


# Calculate Alt/Az (Check for empty times)
if len(times) > 0:
    try:
        altaz = object_coordinates.transform_to(AltAz(obstime=times, location=observer_location))

        # Check for empty altaz and times
        if len(altaz) > 0 and len(times) > 0:
            # Create Pandas DataFrame for Bokeh
            df = pd.DataFrame({
                'time': times.utc.datetime,
                'altitude': altaz.alt.to_value(u.deg)  # Correct way to get degrees
                #'azimuth': altaz.az.to_value(u.deg)   # Correct way to get degrees
            })
            source = ColumnDataSource(df)
            
            # Create Bokeh plot
            if object_name== "Custom":
                p = figure(x_axis_type="datetime", title=f"Visibility of {custom_name} from {latitude}, {longitude}",
                       x_axis_label="UTC Time",
                       y_axis_label="Altitude (degrees)",
                       height=600, width=800, tools="hover,pan,box_zoom,reset,save", toolbar_location="above")
            else:
                p = figure(x_axis_type="datetime", title=f"Visibility of {object_name} from {latitude}, {longitude}",
                       x_axis_label="UTC Time",
                       y_axis_label="Altitude (degrees)",
                       height=600, width=800, tools="hover,pan,box_zoom,reset,save", toolbar_location="above")
            p.background_fill_color = "beige"
            p.background_fill_alpha = 0.3
            p.outline_line_color = "black"
            p.outline_line_width = 1.5

            # Night/Day Shading (Improved - Correct Bokeh Layering)
            sun_altaz = get_sun(times).transform_to(AltAz(obstime=times, location=observer_location))
            night_mask = sun_altaz.alt < 0 * u.deg

            # The fix is here:
            night_mask = night_mask.astype(np.bool_) # Convert to np.bool_ explicitly. Or just let it be.
            formatter1 = DatetimeTickFormatter(days="%H-%M-%S")
            if night_mask.size > 0:  # Check if it's empty *after* conversion
                for i in range(len(times) - 1):
                    if night_mask[i] != night_mask[i + 1]:
                        if object_name== "Custom": 
                            p.scatter(x='time', y='altitude', source=source, color="blue", size=10, legend_label=custom_name) 
                            p.line([df['time'].min(), df['time'].max()], [10, 10], color="orange", line_width=3, line_dash="dashed", legend_label="horizon limit: 10 deg")  # Span the entire x-axis
                            p.xaxis.major_label_text_font_size = "14pt"
                            p.yaxis.major_label_text_font_size = "14pt" 
                            p.xaxis.axis_label_text_font_size = "16pt"
                            p.yaxis.axis_label_text_font_size = "16pt"
                            p.xaxis.major_tick_in = 10
                            p.xaxis.minor_tick_in = 3  
                            p.yaxis.major_tick_in = 10
                            p.yaxis.minor_tick_in = 3
                            
                            #p.xaxis.formatter = formatter1                  
                        else: 
                            p.scatter(x='time', y='altitude', source=source, color="blue", size=10, legend_label=object_name) 
                            p.line([df['time'].min(), df['time'].max()], [10, 10], color="orange", line_width=3, line_dash="dashed", legend_label="horizon limit: 10 deg")  # Span the entire x-axis
                            p.xaxis.major_label_text_font_size = "14pt"
                            p.yaxis.major_label_text_font_size = "14pt"
                            p.xaxis.axis_label_text_font_size = "16pt"
                            p.yaxis.axis_label_text_font_size = "16pt"
                            p.xaxis.major_tick_in = 10
                            p.xaxis.minor_tick_in = 3  
                            p.yaxis.major_tick_in = 10
                            p.yaxis.minor_tick_in = 3
                            #p.xaxis.formatter = formatter1
                        #p.line(x='time', y='azimuth', source=source, legend_field="Azimuth", color="blue", line_width=2)

            # Hover Tool
            hover = p.select(tools=dict(type=HoverTool))
            hover.mode = "vline"
            hover.tooltips = [
                ("Time", "@time{%Y-%m-%d %H:%M:%S}"),
                ("Altitude", "@altitude{0.2f} degrees")
                #("Azimuth", "@azimuth{0.2f} degrees"),
            ]

            st.bokeh_chart(p, use_container_width=True)

        else:
            st.warning("Could not calculate or plot altitude. Check your inputs.")

    except Exception as e:
        st.error(f"Error calculating Alt: {e}")

else:
    st.warning("No times were generated. Check your date range.")

# Display some additional information:
if 'object_coordinates' in locals(): #check if object_coordinates is defined to avoid errors.
    st.write(f"**Object Coordinates (ICRS)**: RA={object_coordinates.ra}, Dec={object_coordinates.dec}")

st.write(f"**Observer Location**: Latitude={latitude}, Longitude={longitude}, Elevation={elevation}")