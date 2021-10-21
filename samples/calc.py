from datetime import date
import logging

import pytz
import rasterio

from tppss import horizon, sunrise_sunset, sunrise_sunset_year

FORMAT = "%(asctime)-15s %(levelname)s %(name)s %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

logger = logging.getLogger(__package__)


latlon = (45.837491, 6.218549)
# imperial (45.903562, 6.144632)  # (45.903627, 6.144714)
tz = pytz.timezone("Europe/Paris")
day = date(2021, 10, 20)
with rasterio.open(
    "C:\\Users\\gvellut\\Documents\\projects\\sunlight\\rge74_wgs84_cubic.tif"
) as dataset:
    logger.info("Computing horizon...")
    horizon = horizon(latlon, dataset, height=2, precision=2)
    logger.info("Computing sunrise / sunset...")
    res = sunrise_sunset(latlon, horizon, day, tz, precision=60)
    sunrise, sunset, is_light_all_day = res
    logger.info(f"Result: Sunrise: {sunrise} Sunset: {sunset}")

    logger.info("Computing sunrise / sunset for year 2021...")
    sunsuns = sunrise_sunset_year(latlon, horizon, 2021, tz, 60)
    logger.info("End")
    format_day = "%Y-%m-%d"
    format_sunsun = "%H:%M:%S"
    for i in range(280, 320):
        sunsun = sunsuns[i]
        day, sunrise, sunset, _ = sunsun
        print(
            f"{day.strftime(format_day)},"
            f"{sunrise.strftime(format_sunsun)},"
            f"{sunset.strftime(format_sunsun)}"
        )
