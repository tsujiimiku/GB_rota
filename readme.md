# Celestial Body Altitude Calculator

This tool calculates the rising and setting times of celestial bodies (Moon, Jupiter, Venus, Mars, and Sun) when they reach specified altitudes at a given location.

## Features

- Calculates rising and setting times for:
  - Moon
  - Jupiter
  - Venus
  - Mars
  - Sun
- Outputs results in CSV format
- Configurable altitude thresholds for each celestial body
- Batch processing capabilities using the provided shell script

## Requirements

- Python 3.x
- Required Python packages:
  - astropy
  - numpy
  - pandas

## Installation

1. Clone this repository:
```bash
git clone https://github.com/YOUR_USERNAME/celestial-altitude-calculator.git
cd celestial-altitude-calculator
```

2. Install required packages:
```bash
pip install astropy numpy pandas
```

## Usage

### Single Month Calculation

```bash
python altitude_calculation.py --date YYYY/MM [options]
```

Options:
- `--moon`: Altitude threshold for Moon (default: 70°)
- `--jupiter`: Altitude threshold for Jupiter (default: 70°)
- `--venus`: Altitude threshold for Venus (default: 70°)
- `--mars`: Altitude threshold for Mars (default: 70°)
- `--sun`: Altitude threshold for Sun (default: 70°)
- `--output`: Output directory for CSV file (default: current directory)

Example:
```bash
python altitude_calculation.py --date 2024/01 --moon 70 --sun 50 --output ./results
```

### Batch Processing

The included `test.sh` script can process multiple months:

```bash
chmod +x test.sh
./test.sh
```

By default, it will process dates from 2024 to 2027 with the following altitude thresholds:
- Moon: 70°
- Jupiter: 70°
- Venus: 70°
- Mars: 70°
- Sun: 50°

## Output Format

The script generates CSV files named `YYYY-MM.csv` containing the following columns:
- Date
- Moon Rising/Setting times
- Jupiter Rising/Setting times
- Venus Rising/Setting times
- Mars Rising/Setting times
- Sun Rising/Setting times

## Location Settings

The calculations are currently set for the following location:
- Latitude: 28° 18' 00" N
- Longitude: 16° 30' 35" W
- Height: 2390m

To change the location, modify the corresponding variables in the `altitude_calculation.py` script.

## License

[Add your chosen license here]
