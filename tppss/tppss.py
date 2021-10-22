from collections import namedtuple
from datetime import date, datetime, time, timedelta
from fractions import Fraction
import logging

import numpy as np

from .sunpos import sunpos

logger = logging.getLogger(__name__)

SunriseSunset = namedtuple(
    "SunriseSunset", "sunrise_time sunset_time is_light_or_night_all_day"
)


def sunrise_sunset(latlon, horizon, day, tz, precision=1):

    times = times_in_day(day, tz, precision)
    above_horizon_indices, _ = above_horizon(latlon, times, horizon)

    if len(above_horizon_indices) == 0:
        # Night all day
        # TODO test in Norway for example
        return SunriseSunset(None, None, False)
    elif len(above_horizon_indices) == len(times):
        # Light all day
        return SunriseSunset(None, None, True)

    # take first for sunrise
    # last for sunset
    # possible that it is hidden again (between peaks):
    # TODO return all transitions instead of just first ?
    # TODO longest sequence ?

    sunrise_time_index = above_horizon_indices[0]
    sunrise_time = times[sunrise_time_index]
    sunset_time_index = above_horizon_indices[-1]
    sunset_time = times[sunset_time_index]

    return SunriseSunset(sunrise_time, sunset_time, None)


def times_in_day(day, tz, precision=1):
    # 12 am
    base_time = datetime.combine(day, time(0, 0, 0))
    base_time = base_time.replace(tzinfo=tz)

    # fraction of an hour
    freq = Fraction(1, precision)

    times = []
    for i in range(24 * precision + 1):
        step = i * freq
        t = base_time + timedelta(hours=float(step))
        times.append(t)

    return times


def above_horizon(latlon, times, horizon):
    helevations, _, hazimuths = horizon

    azimuths = []
    elevations = []
    for i in range(len(times)):
        azimuth, elevation = sunpos(times[i], latlon)
        azimuths.append(azimuth)
        elevations.append(elevation)

    # take nearest
    centers = (hazimuths[1:] + hazimuths[:-1]) / 2
    indices = np.digitize(azimuths, centers)
    # terrain elevation at sun positions for each time step
    helevations_by_time = helevations[indices]
    above_horizon_indices = np.nonzero(elevations - helevations_by_time > 0)
    above_horizon_indices = above_horizon_indices[0]

    return above_horizon_indices, helevations_by_time


def sunrise_sunset_year(latlon, horizon, year, tz, precision=1):
    day_first = date(year, 1, 1)
    day_end = date(year, 12, 31)
    delta = day_end - day_first

    sunsuns = []
    for i in range(delta.days + 1):
        day = day_first + timedelta(days=i)
        sunsun = sunrise_sunset(latlon, horizon, day, tz, precision)
        sunsuns.append((day, *sunsun))

    return sunsuns
