# TPPSS (TopopoSunsun)

TPPSS computes sunrise / sunset times taking into account local topography provided by the user as a DEM (it can also simply compute an horizon based on the DEM). It provides a library and a command line script to compute sunrise / sunset for a location for a single day or a whole year.

# Install

The tool requires Python 3.9+.

## Library

To install just the library, launch:

```console
pip install tppss
```

### RasterIO

Depending on your Python version and platform , the `rasterIO` library may not be available in binary form on PyPI. It is recommanded to use Conda to first create an environment then install TPPSS with:

```console
conda install -c conda-forge rasterio
pip install tppss
```

The second command will pick up the `rasterIO` installed by Conda.

## CLI

The `cli` extra includes a command-line script, also named `tppss`, whose functionality is described below. It has a few additional dependencies.

Launch:

```console
pip install tppss[cli]
```

# Library documentation

The library documentation is available [here](http://example.com). (*not yet*)

# CLI instructions

The `tppss` tool has 2 subcommands:

- `day`
- `year` 

## Day

The `day` subcommand prints to the console the sunrise and sunset times for the chosen day.

```console
~$ tppss day --help
Usage: tppss day [OPTIONS]

  Compute sunset / sunrise time for a single day

Options:
  -p, -position LATLON       Latitude and Longitude of the location to
                             consider (eg 45.2315,5.8389)  [required]

  -m, --dem FILE             DEM in TIFF format and Geographic CRS (eg WGS4)
                             [required]

  -j, --day DATE             Date to consider (in YYYY-MM-DD format)
                             [required]

  -t, --timezone TEXT        Timezone for the result [Default: Timezone of the
                             local machine]

  --distance INTEGER         Distance from the position to consider when
                             computing the horizon (in KM)  [default: 25]

  --angle-precision INTEGER  Precision of horizon angles (for each degree)
                             [default: 1]

  --time-precision INTEGER   Precision of times (for each hour)  [default: 60]
  --help                     Show this message and exit.
```


### Example

```
tppss day -m "C:\Users\gvellut\Documents\projects\sunlight\rge74_wgs84_cubic.tif" -j
2021-11-01 -p "46.179317,6.234106"
```

It outputs:

```
Timezone set to local: 'Romance Daylight Time'
Compute horizon...
Compute sunrise / sunset...
Sunrise: 2021-11-01 07:47:00+01:00 / Sunset: 2021-11-01 16:21:00+01:00
```

The first line is a warning indicating which timezone was selected, since none was indicated on the command-line. 'Romance Daylight Time' is the name of the local timezone on my Windows computer.

In case there is no sunset or sunrise (maybe in a location surrounded by high mountains or above the arctic circle), the last line is either: `Light all day!` or `Night all day!`.

## Year

The `year` subcommand outputs to a CSV file the sunrise and sunset times for each day of the chosen year.

```console
~$ tppss year --help
Usage: tppss year [OPTIONS]

  Compute sunset / sunrise time for a whole year

Options:
  -p, -position LATLON       Latitude and Longitude of the location to
                             consider (eg 45.2315,5.8389)  [required]

  -m, --dem FILE             DEM in TIFF format and Geographic CRS (eg WGS4)
                             [required]

  -y, --year YEAR            Year to take into account  [required]
  -t, --timezone TEXT        Timezone for the result [Default: Timezone of the
                             local machine]

  -o, --csv FILE             CSV to export  [required]
  --distance INTEGER         Distance from the position to consider when
                             computing the horizon (in KM)  [default: 25]

  --angle-precision INTEGER  Precision of horizon angles (for each degree)
                             [default: 1]

  --time-precision INTEGER   Precision of times (for each hour)  [default: 60]
  --help                     Show this message and exit.
```

### Example

```
tppss year -m "C:\Users\gvellut\Documents\projects\sunlight\rge74_wgs84_cubic.tif" -y
2021 -p "46.179317,6.234106" -t Europe/Paris -o ss2021.csv
```

It outputs to the console:

```
Compute horizon...
Compute sunrise / sunset for year 2021...
Writing results to C:\Users\gvellut\Documents\projects\github\calcsoleil\ss2021.csv...
```

The CSV is as follows:

```
DAY,SUNRISE,SUNSET
2021-01-01,08:48+0100,13:59+0100
2021-01-02,08:47+0100,13:59+0100
2021-01-03,08:47+0100,14:01+0100
2021-01-04,08:46+0100,14:04+0100
2021-01-05,08:46+0100,14:04+0100
2021-01-06,08:44+0100,14:05+0100
....
```

The times for sunset and sunrise also indicate the offset from UTC (which can change in the year with DST).

If there is no sunset or sunrise, the second and third columns have value `NA`.

## Some notes

The Latitude and Longitude must be in the same CRS as the DEM.

The DEM (Digital Elevation Model) raster must contain heights above the ellipsoid. If the DEM contains instead the altitude above the geoid (for example, EGM96 if SRTM is used), it should first be transformed into heights. However, since only differences of altitudes are considered, it shouldn't matter much depending on the specific purpose (I use this tool for planning photography outings and I am fine with the precision I get). If needed, the ellipsoidal heights can also be computed using GDAL (`gdalwarp`: See the example at the bottom of [this page](https://gdal.org/programs/gdalwarp.html)).

The value for the timezone option is something like `Europe/Paris` or `MST`. If not present, it is taken from the local machine. If the timezone has DST, the change is reflected in the times computed for surises and sunsets.

The `--distance` option indicates how far away from the position should heights be extracted from the DEM when computing the horizon.

# Notes

## COG

### Transform base GeoTIFF

#### gdal

```sh
# -of COG: Use the COG driver (GDAL 3.1+)
# -co COMPRESS=LZW: Use LZW or DEFLATE compression (DEFLATE is smaller, LZW is faster)
# -co OVERVIEWS=IGNORE_EXISTING: Recompute fresh overviews

gdal_translate input_dem.tif output_cog.tif \
    -of COG \
    -co COMPRESS=LZW \
    -co OVERVIEWS=IGNORE_EXISTING \
    -co RESAMPLING=BILINEAR
```

#### rio-cogeo

For DEM :

```sh
rio cogeo create input_dem.tif output_dem_cog.tif \
  --profile deflate \
  --overview-resampling bilinear
  ```

For sat iamge use, `--profile webp --web-optimized`.

##### Use with uvx

`uvx --from rasterio --with rio-cogeo rio cogeo ...`

### Validate

#### gdal

`gdalinfo ....tif`

- If it says Driver: GTiff/GeoTIFF, look for Metadata: layout=COG.
- If the file is hosted remotely (gs://), look for: Band 1 Block=512x512. If you see Block=25000x1 (strip layout), it is not a COG.

#### rio-cogeo 

`rio cogeo validate output_dem_cog.tif`

# Acknowledgements

The horizon computation is based on this Matlab code (adapted for Python + Numpy):

Benjamin Pillot (2021). DEM-based topography horizon model 
(https://www.mathworks.com/matlabcentral/fileexchange/59421-dem-based-topography-horizon-model),
MATLAB Central File Exchange. Retrieved October 19, 2021. 


The sun position computation is adapted from this code:

John Clark Craig. Python Sun Position for Solar Energy and Research
(https://levelup.gitconnected.com/python-sun-position-for-solar-energy-and-research-7a4ead801777)


# TODO 

- Example that draws horizon + sun course through the sky or just sunrise / sunset position
- Sample: Document how the input DEM is obtained from the RGE 
- numba or faster computations
- correct issue when pos lat/lon near border of DEM cell : if too steep, will seem higher than it should
- web interface for deployment on Cloud Run
- Generate doc
- Support projected CRS (use vertical compensation see https://www.usna.edu/Users/oceano/pguth/md_help/html/demb30q0.htm)