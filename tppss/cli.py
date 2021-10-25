from datetime import datetime
import logging
import sys
import traceback

import click
import colorama
from colorama import Fore, Style
from dateutil import tz as dutz
import rasterio

from . import horizon, KM, sunrise_sunset, sunrise_sunset_year

colorama.init()

logger = logging.getLogger(__package__)

# specify colors for different logging levels
LOG_COLORS = {logging.ERROR: Fore.RED, logging.WARNING: Fore.YELLOW}


class ColorFormatter(logging.Formatter):
    def format(self, record, *args, **kwargs):
        if record.levelno in LOG_COLORS:
            record.msg = "{color_begin}{message}{color_end}".format(
                message=record.msg,
                color_begin=LOG_COLORS[record.levelno],
                color_end=Style.RESET_ALL,
            )
        return super().format(record, *args, **kwargs)


def setup_logging(is_debug):
    global logger
    if is_debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    formatter = ColorFormatter("%(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


class LatLonParamType(click.ParamType):
    name = "latlon"

    def convert(self, value, param, ctx):
        try:
            parts = value.split(",")
            lat = float(parts[0].strip())
            lon = float(parts[1].strip())
            return lat, lon
        except Exception:
            self.fail(
                f"{value!r} is not a valid Lat Lon (eg '45.235555,5.83890')", param, ctx
            )


class CatchAllExceptionsCommand(click.Command):
    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except Exception as ex:
            raise UnrecoverableJNCEPError(str(ex), sys.exc_info())


class UnrecoverableJNCEPError(click.ClickException):
    def __init__(self, message, exc_info):
        super().__init__(message)
        self.exc_info = exc_info

    def show(self):
        logger.error("*** An unrecoverable error occured ***")
        logger.error(self.message)
        logger.debug("".join(traceback.format_exception(*self.exc_info)))


def colored(s, color):
    return f"{color}{s}{Fore.RESET}"


position_option = click.option(
    "-p",
    "-position",
    "latlon",
    help="Latitude and Longitude of the location to consider (eg 45.2315,5.8389)",
    required=True,
    type=LatLonParamType(),
)
dem_option = click.option(
    "-m",
    "--dem",
    "dem_filepath",
    help="DEM in TIFF format and Geographic CRS (eg WGS4)",
    required=True,
    type=click.Path(exists=True, resolve_path=True, dir_okay=False),
)
timezone_option = click.option(
    "-t",
    "--timezone",
    "timezone",
    help="Timezone for the result [Default: Timezone of the local machine]",
)
distance_option = click.option(
    "--distance",
    "distance",
    help="Distance from the position to consider when computing the horizon (in KM)",
    default=25,
    type=int,
    show_default=True,
)
angle_option = click.option(
    "--angle-precision",
    "angle_precision",
    help="Precision of horizon angles (for each degree)",
    default=1,
    type=int,
    show_default=True,
)
time_option = click.option(
    "--time-precision",
    "time_precision",
    help="Precision of times (for each hour)",
    default=60,
    type=int,
    show_default=True,
)


@click.group()
@click.option(
    "-d",
    "--debug",
    "is_debug",
    is_flag=True,
    help=("Flag to activate debug mode"),
    required=False,
)
@click.pass_context
def main(ctx, is_debug):
    """Computes sunset / sunrise time taking into account local topography"""
    setup_logging(is_debug)
    # special attribute of context
    ctx.obj = {"DEBUG": is_debug}


@main.command(
    "day",
    help="Compute sunset / sunrise time for a single day",
    cls=CatchAllExceptionsCommand,
)
@position_option
@dem_option
@click.option(
    "-j",
    "--day",
    "day",
    metavar="DATE",
    help="Date to consider (in YYYY-MM-DD format)",
    required=True,
    type=click.DateTime(formats=["%Y-%m-%d"]),
)
@timezone_option
@distance_option
@angle_option
@time_option
@click.pass_context
def tppss_day(
    ctx, latlon, dem_filepath, day, timezone, distance, angle_precision, time_precision
):
    tz = dutz.gettz(timezone)
    if timezone is None:
        zone_name = datetime.now(tz).tzname()
        logger.warning(f"Timezone set to local: '{zone_name}'")
    with rasterio.open(dem_filepath) as dataset:
        logger.info("Compute horizon...")
        horizon_ = horizon(
            latlon, dataset, distance=distance * KM, height=2, precision=angle_precision
        )

        logger.info("Compute sunrise / sunset...")
        res = sunrise_sunset(latlon, horizon_, day, tz, precision=time_precision)

        sunrise, sunset, is_light_all_day = res
        if is_light_all_day is None:
            text = f"Sunrise: {sunrise} / Sunset: {sunset}"
        else:
            if is_light_all_day:
                text = "Light all day!"
            else:
                text = "Night all day!"

        logger.info(colored(text, Fore.GREEN))


@main.command(
    "year",
    help="Compute sunset / sunrise time for a whole year",
    cls=CatchAllExceptionsCommand,
)
@position_option
@dem_option
@click.option(
    "-y",
    "--year",
    "year",
    metavar="YEAR",
    help="Year to take into account",
    required=True,
    type=int,
)
@timezone_option
@click.option(
    "-o",
    "--csv",
    "csv_filepath",
    help="CSV to export",
    type=click.Path(resolve_path=True, dir_okay=False, writable=True),
    required=True,
)
@distance_option
@angle_option
@time_option
@click.pass_context
def tppss_year(
    ctx,
    latlon,
    dem_filepath,
    year,
    timezone,
    csv_filepath,
    distance,
    angle_precision,
    time_precision,
):
    if not 1901 <= year <= 2099:
        logger.warning(
            "Sun position computation may not be accurate outside years 1901 to 2099!"
        )

    tz = dutz.gettz(timezone)
    if timezone is None:
        zone_name = datetime.now(tz).tzname()
        logger.warning(f"Timezone set to local: '{zone_name}'")

    with rasterio.open(dem_filepath) as dataset:
        logger.info("Compute horizon...")
        horizon_ = horizon(
            latlon, dataset, distance=distance * KM, height=2, precision=angle_precision
        )

        logger.info(f"Compute sunrise / sunset for year {year}...")
        sunsuns = sunrise_sunset_year(
            latlon, horizon_, year, tz, precision=time_precision
        )

        logger.info(f"Write results to {csv_filepath}...")
        with open(csv_filepath, "w", encoding="utf-8") as f:
            print_output(f, year, sunsuns)


def print_output(file_, year, sunsuns):
    format_day = "%Y-%m-%d"
    format_sunsun = "%H:%M:%S%z"
    file_.write("DAY,SUNRISE,SUNSET\n")
    for i in range(len(sunsuns)):
        day, sunrise, sunset, is_light_or_night_all_day = sunsuns[i]
        if is_light_or_night_all_day is None:
            file_.write(
                f"{day.strftime(format_day)},"
                f"{sunrise.strftime(format_sunsun)},"
                f"{sunset.strftime(format_sunsun)}\n"
            )
        else:
            file_.write(f"{day.strftime(format_day)},NA,NA\n")


if __name__ == "__main__":
    main()
