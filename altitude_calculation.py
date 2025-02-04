import argparse
from astropy.coordinates import EarthLocation, AltAz, get_moon, get_sun, get_body, Latitude, Longitude
from astropy.time import Time
import astropy.units as u
import numpy as np
import pandas as pd
from collections import defaultdict
from datetime import datetime, timedelta
import os
import calendar

def calculate_altitude_60deg_times(body_func, altitude, latitude, longitude, height, start_time, end_time):
    location = EarthLocation(lat=latitude, lon=longitude, height=height)
    t_start = Time(start_time)
    t_end = Time(end_time) + timedelta(days=1)
    times = Time(np.linspace(t_start.mjd, t_end.mjd, num=10000), format='mjd')

    body_altaz = body_func(times).transform_to(AltAz(obstime=times, location=location))
    body_alt = body_altaz.alt.deg

    def crossing_time(alt, alts, times):
        is_upper_prev = None
        ret_list = []
        for i in range(1, len(alts)):
            is_upper = alts[i] > alt
            if is_upper_prev is not None and is_upper != is_upper_prev:
                alt_prev, alt_curr = alts[i-1], alts[i]
                t_prev, t_curr = times[i-1].mjd, times[i].mjd
                t_crossing_mjd = (alt_curr - alt) / (alt_curr - alt_prev) * t_prev + \
                                 (alt - alt_prev) / (alt_curr - alt_prev) * t_curr
                ret_list.append((Time(t_crossing_mjd, format='mjd'), is_upper))
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

    df = pd.DataFrame([
        {'Date': date, 'Rising time': times['rising'], 'Setting time': times['setting']}
        for date, times in sorted(daily_times.items())
    ])

    if not df.empty:
        df['Date'] = pd.to_datetime(df['Date'])
    else:
        df = pd.DataFrame(columns=['Date', 'Rising time', 'Setting time'])
        df['Date'] = pd.date_range(start=start_time, end=end_time)

    return df

def add_months(date, months):
    month = date.month - 1 + months
    year = date.year + month // 12
    month = month % 12 + 1
    day = min(date.day, calendar.monthrange(year, month)[1])
    return date.replace(year=year, month=month, day=day)

def str_to_datetime(date_str):
    """
    'YYYY/MM' 形式の文字列を datetime オブジェクトに変換する
    """
    year, month = map(int, date_str.split('/'))
    return datetime(year, month, 1, 0, 0, 0)

def main():
    parser = argparse.ArgumentParser(description="Calculate times when celestial bodies reach a specified altitude")
    parser.add_argument("--moon", type=float, default=70.0, help="Altitude for the Moon (default 70°)")
    parser.add_argument("--jupiter", type=float, default=70.0, help="Altitude for Jupiter (default 70°)")
    parser.add_argument("--venus", type=float, default=70.0, help="Altitude for venus (default 70°)")
    parser.add_argument("--mars", type=float, default=70.0, help="Altitude for Mars (default 70°)")
    parser.add_argument("--sun", type=float, default=70.0, help="Altitude for the Sun (default 70°)")
    parser.add_argument("--date", type=str, required=True, help="Year and month in YYYY/MM format")
    parser.add_argument("--output", type=str, default="", help="Output directory for CSV file (default: current directory)")

    args = parser.parse_args()

    latitude = Latitude((+28, 18, 00), unit='deg')
    longitude = Longitude((-16, 30, 35), unit='deg')
    height = 2390 * u.m

    # 指定された年月の初日と最終日を計算
    year, month = map(int, args.date.split('/'))
    start_time = datetime(year, month, 1)
    _, last_day = calendar.monthrange(year, month)
    end_time = datetime(year, month, last_day)

    start_time_str = start_time.strftime('%Y-%m-%d')
    end_time_str = end_time.strftime('%Y-%m-%d')

    moon_df = calculate_altitude_60deg_times(get_moon, args.moon, latitude, longitude, height, start_time_str, end_time_str)
    jupiter_df = calculate_altitude_60deg_times(lambda time: get_body('jupiter', time), args.jupiter, latitude, longitude, height, start_time_str, end_time_str)
    venus_df = calculate_altitude_60deg_times(lambda time: get_body('venus', time), args.venus, latitude, longitude, height, start_time_str, end_time_str)
    mars_df = calculate_altitude_60deg_times(lambda time: get_body('mars', time), args.mars, latitude, longitude, height, start_time_str, end_time_str)
    sun_df = calculate_altitude_60deg_times(get_sun, args.sun, latitude, longitude, height, start_time_str, end_time_str)

    # 合体するデータフレームの準備
    combined_df = pd.DataFrame({'Date': pd.date_range(start=start_time_str, end=end_time_str)})
    combined_df = combined_df.set_index('Date')

    # 各データフレームをマージ
    for df, name, alt in [(moon_df, 'Moon', args.moon),
                          (jupiter_df, 'Jupiter', args.jupiter),
                          (venus_df, 'Venus', args.venus),
                          (mars_df, 'Mars', args.mars),
                          (sun_df, 'Sun', args.sun)]:
        if not df.empty:
            combined_df = combined_df.join(df.set_index('Date').rename(columns={
                'Rising time': f'{name} Rising ({alt}°)',
                'Setting time': f'{name} Setting ({alt}°)'
            }))
        else:
            combined_df[f'{name} Rising ({alt}°)'] = None
            combined_df[f'{name} Setting ({alt}°)'] = None

    combined_df = combined_df.reset_index()

    # 出力ディレクトリの設定
    output_dir = args.output if args.output else os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    file_name = f"{year}-{month:02d}.csv"
    file_path = os.path.join(output_dir, file_name)

    combined_df.to_csv(file_path, index=False)
    print(f"Results saved to {file_path}")

if __name__ == '__main__':
    main()
