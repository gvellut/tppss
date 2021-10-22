from datetime import datetime, timezone
import math

# from
# https://levelup.gitconnected.com/python-sun-position-for-solar-energy-
# and-research-7a4ead801777


class NoTimeZoneException(Exception):
    pass


def sunpos(when: datetime, location, refraction=False):
    if not when.tzinfo:
        raise NoTimeZoneException(f"{when} doesn't have a timezone !")

    latitude, longitude = location
    # Math typing shortcuts
    rad, deg = math.radians, math.degrees
    sin, cos, tan = math.sin, math.cos, math.tan
    asin, atan2 = math.asin, math.atan2
    # Convert latitude and longitude to radians
    rlat = rad(latitude)
    rlon = rad(longitude)
    # Decimal hour of the day at Greenwich
    when_utc = when.astimezone(timezone.utc)
    dec_utc = when_utc.hour + when_utc.minute / 60 + when_utc.second / 3600
    # Days from J2000, accurate from 1901 to 2099
    daynum = (
        367 * when_utc.year
        - 7 * (when_utc.year + (when_utc.month + 9) // 12) // 4
        + 275 * when_utc.month // 9
        + when_utc.day
        - 730531.5
        + dec_utc / 24
    )
    # Mean longitude of the sun
    mean_long = daynum * 0.01720279239 + 4.894967873
    # Mean anomaly of the Sun
    mean_anom = daynum * 0.01720197034 + 6.240040768
    # Ecliptic longitude of the sun
    eclip_long = (
        mean_long
        + 0.03342305518 * sin(mean_anom)
        + 0.0003490658504 * sin(2 * mean_anom)
    )
    # Obliquity of the ecliptic
    obliquity = 0.4090877234 - 0.000000006981317008 * daynum
    # Right ascension of the sun
    rasc = atan2(cos(obliquity) * sin(eclip_long), cos(eclip_long))
    # Declination of the sun
    decl = asin(sin(obliquity) * sin(eclip_long))
    # Local sidereal time
    sidereal = 4.894961213 + 6.300388099 * daynum + rlon
    # Hour angle of the sun
    hour_ang = sidereal - rasc
    # Local elevation of the sun
    elevation = asin(sin(decl) * sin(rlat) + cos(decl) * cos(rlat) * cos(hour_ang))
    # Local azimuth of the sun
    azimuth = atan2(
        -cos(decl) * cos(rlat) * sin(hour_ang),
        sin(decl) - sin(rlat) * sin(elevation),
    )
    # Convert azimuth and elevation to degrees
    azimuth = _into_range(deg(azimuth), 0, 360)
    elevation = _into_range(deg(elevation), -180, 180)
    # Refraction correction (optional)
    if refraction:
        targ = rad((elevation + (10.3 / (elevation + 5.11))))
        elevation += (1.02 / tan(targ)) / 60
    # Return azimuth and elevation in degrees
    return (round(azimuth, 2), round(elevation, 2))


def _into_range(x, range_min, range_max):
    shiftedx = x - range_min
    delta = range_max - range_min
    return (((shiftedx % delta) + delta) % delta) + range_min
