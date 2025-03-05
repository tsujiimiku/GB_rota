import argparse
from astropy.coordinates import EarthLocation, AltAz, get_body, Latitude, Longitude
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
from collections import defaultdict
from datetime import datetime, timedelta
import os
import calendar
from concurrent.futures import ProcessPoolExecutor

"""
This script calculates the times when celestial bodies (Moon, Jupiter, Venus, Mars, Sun)
reach specific altitudes as seen from GroundBIRD:

Usage:
  python script.py --date YYYY/MM [--moon ALT] [--jupiter ALT] [--venus ALT] 
                   [--mars ALT] [--sun ALT] [--output DIR]

Arguments:
  --date YYYY/MM    Target year and month (required)
  --moon ALT        Altitude threshold for the Moon (default: 70°)
  --jupiter ALT     Altitude threshold for Jupiter (default: 70°)
  --venus ALT       Altitude threshold for Venus (default: 70°)
  --mars ALT        Altitude threshold for Mars (default: 70°)
  --sun ALT         Altitude threshold for the Sun (default: 50°)
  --output DIR      Output directory (default: current directory)

Output:
  CSV file with rising and setting times for each celestial body crossing the specified altitude for each day of the month
"""

def calculate_altitude_crossing_times(body_name, altitude, latitude, longitude, height, start_time, end_time):
    """
    Calculate the times when a specified celestial body reaches a specific altitude
    
    Parameters:
        body_name (str): Name of the celestial body ('moon', 'jupiter', etc.)
        altitude (float): Target altitude (degrees)
        latitude (Latitude): Observer's latitude
        longitude (Longitude): Observer's longitude
        height (Quantity): Observer's elevation
        start_time (str): Start date for calculation (YYYY-MM-DD format)
        end_time (str): End date for calculation (YYYY-MM-DD format)
        
    Returns:
        pandas.DataFrame: DataFrame containing rising and setting times for each date
    """
    location = EarthLocation(lat=latitude, lon=longitude, height=height)
    t_start = Time(start_time)
    t_end = Time(end_time) + timedelta(days=1)
    sample_count = 2000 # More points increase accuracy but extend calculation time. 1000 is sufficient for this script.
    times = Time(np.linspace(t_start.mjd, t_end.mjd, num=sample_count), format='mjd')
    
    body_altaz = get_body(body_name, times).transform_to(AltAz(obstime=times, location=location))
    body_alt = body_altaz.alt.deg

    # Function to find altitude crossing time using linear interpolation
    def find_crossing_time(idx, alt, alts, times):
        """
        Estimate the exact crossing time using linear interpolation between two points
        
        Parameters:
            idx (int): Index of the later time point
            alt (float): Target altitude
            alts (array): Array of altitudes
            times (Time): Array of times
            
        Returns:
            Time: Estimated crossing time
        """
        # Initial estimation
        t_prev, t_curr = times[idx-1].mjd, times[idx].mjd
        alt_prev, alt_curr = alts[idx-1], alts[idx]
        
        # Calculate crossing time using linear interpolation
        t_crossing_mjd = t_prev + (t_curr - t_prev) * (alt - alt_prev) / (alt_curr - alt_prev)
        
        return Time(t_crossing_mjd, format='mjd')

    def crossing_time(alt, alts, times):
        is_upper_prev = None
        ret_list = []
        for i in range(1, len(alts)):
            is_upper = alts[i] > alt
            if is_upper_prev is not None and is_upper != is_upper_prev:
                crossing = find_crossing_time(i, alt, alts, times)
                ret_list.append((crossing, is_upper))
            is_upper_prev = is_upper
        return ret_list

    crossing_times = crossing_time(altitude, body_alt, times)
    daily_times = defaultdict(lambda: {'rising': None, 'setting': None})
    
    for time, is_rising in crossing_times:
        date = time.datetime.date()
        time_str = time.datetime.time().strftime('%H:%M:%S')
        if is_rising:
            daily_times[date]['rising'] = time_str
        else:
            daily_times[date]['setting'] = time_str

    data = []
    for date, times in sorted(daily_times.items()):
        data.append({
            'Date': date, 
            'Rising time': times['rising'], 
            'Setting time': times['setting']
        })
    
    df = pd.DataFrame(data)

    if not df.empty:
        df['Date'] = pd.to_datetime(df['Date'])
    else:
        df = pd.DataFrame(columns=['Date', 'Rising time', 'Setting time'])
        df['Date'] = pd.date_range(start=start_time, end=end_time)

    return df

def calculate_body_times(args):
    """
    Wrapper function for parallel processing
    
    Takes a tuple of arguments for ProcessPoolExecutor parallel execution
    and returns a pair of body name and calculation results.
    
    Parameters:
        args (tuple): (body_name, altitude, latitude, longitude, height, start_time, end_time)
        
    Returns:
        tuple: (body_name, result DataFrame)
    """
    body_name, altitude, latitude, longitude, height, start_time_str, end_time_str = args
    return body_name, calculate_altitude_crossing_times(
        body_name, altitude, latitude, longitude, height, start_time_str, end_time_str
    )

def main():
    """
    Main execution function
    
    Parses command line arguments, calculates altitude crossing times for
    specified celestial bodies, and outputs results to a CSV file.
    Uses parallel processing to reduce computation time for multiple bodies.
    """
    parser = argparse.ArgumentParser(description="Optimized Celestial Body Altitude Crossing Calculator")
    parser.add_argument("--moon", type=float, default=70.0, help="Altitude threshold for the Moon (default: 70°)")
    parser.add_argument("--jupiter", type=float, default=70.0, help="Altitude threshold for Jupiter (default: 70°)")
    parser.add_argument("--venus", type=float, default=70.0, help="Altitude threshold for Venus (default: 70°)")
    parser.add_argument("--mars", type=float, default=70.0, help="Altitude threshold for Mars (default: 70°)")
    parser.add_argument("--sun", type=float, default=50.0, help="Altitude threshold for the Sun (default: 70°)")
    parser.add_argument("--date", type=str, required=True, help="Target year and month (YYYY/MM format)")
    parser.add_argument("--output", type=str, default="", help="Output directory (default: current directory)")

    args = parser.parse_args()

    latitude = Latitude('28d18m00s', unit='deg')
    longitude = Longitude('-16d30m35s', unit='deg')
    height = 2390 * u.m

    # Calculate the first and last days of the specified month
    year, month = map(int, args.date.split('/'))
    start_time = datetime(year, month, 1)
    _, last_day = calendar.monthrange(year, month)
    end_time = datetime(year, month, last_day)

    start_time_str = start_time.strftime('%Y-%m-%d')
    end_time_str = end_time.strftime('%Y-%m-%d')

    # Prepare for parallel processing
    # Bundle calculation parameters for each celestial body
    calculation_args = [
        ('moon', args.moon, latitude, longitude, height, start_time_str, end_time_str),
        ('jupiter', args.jupiter, latitude, longitude, height, start_time_str, end_time_str),
        ('venus', args.venus, latitude, longitude, height, start_time_str, end_time_str),
        ('mars', args.mars, latitude, longitude, height, start_time_str, end_time_str),
        ('sun', args.sun, latitude, longitude, height, start_time_str, end_time_str)
    ]

    # Run calculations for each body in parallel
    print("Calculating altitude crossing times for each body in parallel...")
    results = {}
    with ProcessPoolExecutor(max_workers=5) as executor:
        for body_name, df in executor.map(calculate_body_times, calculation_args):
            results[body_name] = df
            print(f"  Completed calculations for {body_name.capitalize()}")

    # Prepare the DataFrame
    combined_df = pd.DataFrame({'Date': pd.date_range(start=start_time_str, end=end_time_str)})
    combined_df = combined_df.set_index('Date')

    # Merge all DataFrames
    name_map = {'moon': 'Moon', 'jupiter': 'Jupiter', 'venus': 'Venus', 'mars': 'Mars', 'sun': 'Sun'}
    alt_map = {'moon': args.moon, 'jupiter': args.jupiter, 'venus': args.venus, 'mars': args.mars, 'sun': args.sun}
    
    for body_name, df in results.items():
        name = name_map[body_name]
        alt = alt_map[body_name]
        if not df.empty:
            combined_df = combined_df.join(df.set_index('Date').rename(columns={
                'Rising time': f'{name} Rising ({alt}°)',
                'Setting time': f'{name} Setting ({alt}°)'
            }))
        else:
            combined_df[f'{name} Rising ({alt}°)'] = None
            combined_df[f'{name} Setting ({alt}°)'] = None

    combined_df = combined_df.reset_index()

    # Set up output directory
    output_dir = args.output if args.output else os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    file_name = f"{year}-{month:02d}.csv"
    file_path = os.path.join(output_dir, file_name)

    combined_df.to_csv(file_path, index=False)
    print(f"Results saved to {file_path}")

if __name__ == '__main__':
    main()