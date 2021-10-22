import logging
import math

import numpy as np
from pyproj import CRS
import rasterio

logger = logging.getLogger(__name__)

KM = 1000.0
MILES = 1609.34

# TODO
# optional use of rasterio data ; extract proj / ellipsoid from it

# TODO handle nodata (0 dans le TIFF IGN) ;
# TIFF IGN EPSG:2154 - RGF93 / Lambert-93 - Projected

# TODO support negative elevation angle


def horizon(latlon, raster, distance=25 * KM, precision=1, height=0):
    if raster.crs.is_projected:
        raise Exception("Only geographic CRS are supported")

    crs = CRS.from_wkt(raster.crs.to_wkt())
    ellipsoid = crs.ellipsoid

    logger.debug("Extracting data...")
    study_area, study_transform = _extract_data(latlon, raster, distance, ellipsoid)

    lats, lons = _pixel_positions(study_area, study_transform)

    y_obs, x_obs = rasterio.transform.rowcol(study_transform, latlon[1], latlon[0])

    logger.debug("Computing azimuths...")
    # Azimuth of all pixel corners seen from observer point
    # bottom right corners
    lat_res, lon_res = _resolution(study_transform)
    azimuth = _azimuth(latlon, lats - lat_res / 2, lons + lon_res / 2, ellipsoid)

    # Corresponding elevation angle
    # add observer height
    z_obs = study_area[y_obs, x_obs] + height
    lon_grid, lat_grid = np.meshgrid(lons, lats, sparse=False)

    logger.debug("Computing elevations...")
    elevation = _elevation_angle(
        z_obs,
        study_area,
        latlon,
        lat_grid,
        lon_grid,
        ellipsoid,
    )

    logger.debug(f"Computing sun mask... precision={precision}")
    helevation = _compute_mask(x_obs, y_obs, study_area, precision, azimuth, elevation)

    # Finalization
    # Set horizon
    # TODO negative values should be allowed (when located at peak of mountain)
    helevation[helevation < 0] = 0
    hzenith = 90 - helevation
    # Store azimuth from 0 to 360 (for using with solar azimuth functions)
    hazimuth = np.linspace(0, 360, len(helevation))

    return helevation, hzenith, hazimuth


def _compute_mask(x_obs, y_obs, study_area, precision, azimuth, elevation):
    # Elevation vector length
    # degrees
    length_elevation = precision * 90 + 1
    height, width = study_area.shape

    # Specific azimuth values for NE (-180 to -90) and NW (90 to 180) areas
    # so bin edges are ordered (in the computation below)
    azimuthNE = azimuth.copy()
    azimuthNE[:y_obs, x_obs - 1] = azimuthNE[:y_obs, x_obs - 1] - 360
    azimuthNW = azimuth.copy()
    azimuthNW[:y_obs, x_obs] = azimuthNW[:y_obs, x_obs] + 360

    # Initialization
    north_distance = y_obs
    east_distance = width - x_obs
    south_distance = height - y_obs
    west_distance = x_obs

    # TODO initialize with -inf
    # TODO process NO DATA
    elevationNE = np.zeros((north_distance, length_elevation))
    elevationE = np.zeros((east_distance, 2 * length_elevation - 1))
    elevationS = np.zeros((south_distance, 2 * length_elevation - 1))
    elevationW = np.zeros((west_distance, 2 * length_elevation - 1))
    elevationNW = np.zeros((north_distance, length_elevation))

    azNE = np.linspace(-180, -90, length_elevation)
    azE = np.linspace(-180, 0, 2 * length_elevation - 1)
    azS = np.linspace(-90, 90, 2 * length_elevation - 1)
    azW = np.linspace(0, 180, 2 * length_elevation - 1)
    azNW = np.linspace(90, 180, length_elevation)

    # Main computation
    # North divided into 2 sections : -180 to -90 in W ; 90 to 180 in E
    # Retrieve all elevation angles for iso-azimuth lines (loxodromes / rhumb lines)
    for isoline in range(north_distance):
        k = np.digitize(azNE, azimuthNE[isoline, x_obs - 1 :])
        valid_k = (k != 0) & (k != east_distance + 1)
        elevationNE[isoline, valid_k] = elevation[isoline, x_obs - 1 + k[valid_k]]

        k2 = np.digitize(azNW, azimuthNW[isoline, : x_obs + 1])
        valid_k2 = (k2 != 0) & (k2 != (west_distance + 1))
        elevationNW[isoline, valid_k2] = elevation[isoline, k2[valid_k2] - 1]

    for isoline in range(east_distance):
        k = np.digitize(azE, azimuth[:, x_obs + isoline])
        valid_k = (k != 0) & (k != height)
        elevationE[isoline, valid_k] = elevation[k[valid_k], x_obs + isoline]

    for isoline in range(south_distance):
        k = np.digitize(azS, azimuth[y_obs + isoline, ::-1])
        valid_k = (k != 0) & (k != width)
        elevationS[isoline, valid_k] = elevation[
            y_obs + isoline, width - 1 - k[valid_k]
        ]

    for isoline in range(west_distance):
        k = np.digitize(azW, azimuth[::-1, isoline])
        valid_k = (k != 0) & (k != height)
        elevationW[isoline, valid_k] = elevation[height - 1 - k[valid_k], isoline]

    # max for each angle (2nd dimension) for each sunmask
    sun_maskNE = np.max(elevationNE, axis=0)
    sun_maskE = np.max(elevationE, axis=0)
    sun_maskS = np.max(elevationS, axis=0)
    sun_maskW = np.max(elevationW, axis=0)
    sun_maskNW = np.max(elevationNW, axis=0)

    # Global azimuth (North to North) and sun mask (elevation angle)
    azNtoN = np.concatenate([azNE, azE, azS, azW, azNW])
    sun_mask = np.concatenate([sun_maskNE, sun_maskE, sun_maskS, sun_maskW, sun_maskNW])

    total_length_elevation = precision * 360 + 1
    helevation = np.zeros(total_length_elevation)

    # Corresponding azimuth (from -180° to 180°)
    az = np.linspace(-180, 180, total_length_elevation)
    for i in range(len(az)):
        # max of angle r across all the sunmasks
        helevation[i] = np.max(sun_mask[azNtoN == az[i]])

    return helevation


def _extract_data(latlon, raster, distance, ellipsoid):
    lat, lon = latlon
    lat_resolution, lon_resolution = _resolution(raster.transform)

    row, col = raster.index(lon, lat)

    if not 0 <= row < raster.height or not 0 <= col < raster.width:
        # TODO specific exception
        raise Exception("LatLon not covered by DEM")

    # deg corresponding to distance

    lat_distance, lon_distance = _dist2deg(latlon, distance, ellipsoid)
    topo_lat_distance = round(lat_distance / lat_resolution)
    topo_lon_distance = round(lon_distance / lon_resolution)

    window_base = rasterio.windows.Window(0, 0, raster.width, raster.height)
    window = rasterio.windows.Window(
        # clamp to the raster area
        col - topo_lon_distance,
        row - topo_lat_distance,
        topo_lon_distance * 2,
        topo_lat_distance * 2,
    )

    window = window.intersection(window_base)
    study_area = raster.read(1, window=window)
    study_transform = raster.window_transform(window)

    return study_area, study_transform


def _elevation_angle(
    z_obs,
    study_area,
    latlon,
    lat_grid,
    lon_grid,
    ellipsoid,
):
    latlon = np.deg2rad(latlon)
    lat_grid = np.deg2rad(lat_grid)
    lon_grid = np.deg2rad(lon_grid)

    # Compute cartesian coordinates of point A and B located at altitude
    # z_obs and study_area from the ellipsoid surface (ellipsoidal heights)
    x_A, y_A, z_A = _geographic2cartesian(*latlon, z_obs, ellipsoid)
    x_B, y_B, z_B = _geographic2cartesian(lat_grid, lon_grid, study_area, ellipsoid)

    # Scalar product between AB and normal to the point A
    inner_product = (
        (x_B - x_A) * np.cos(latlon[1]) * np.cos(latlon[0])
        + (y_B - y_A) * np.sin(latlon[1]) * np.cos(latlon[0])
        + (z_B - z_A) * np.sin(latlon[0])
    )

    # Angular elevation computation
    norm = np.sqrt((x_B - x_A) ** 2 + (y_B - y_A) ** 2 + (z_B - z_A) ** 2)
    alpha = np.rad2deg(np.arcsin(inner_product / norm))
    return alpha


def _geographic2cartesian(lat, lon, h, ellipsoid):
    a = ellipsoid.semi_major_metre
    e = _eccentricity(ellipsoid)

    sinlat = np.sin(lat)
    coslat = np.cos(lat)
    # ellipsoid normal
    # TODO optimize for vector case : only 1D (instead of meshgrid) + repeat
    # test if not scalar
    N = a / np.sqrt(1 - (e ** 2) * (sinlat ** 2))

    x = (N + h) * np.cos(lon) * coslat
    y = (N + h) * np.sin(lon) * coslat
    z = (N * (1 - e ** 2) + h) * sinlat

    return x, y, z


# angle of the rhumb line between the points
def _azimuth(latlon, lats, lons, ellipsoid):
    latlon = np.deg2rad(latlon)
    lats = np.deg2rad(lats)
    lons = np.deg2rad(lons)

    L1 = _isometric_latitude(latlon[0], ellipsoid)
    L2 = _isometric_latitude(lats, ellipsoid)
    L2_grid = np.repeat(L2[:, np.newaxis], len(lons), axis=1)

    dlons = latlon[1] - lons
    dlons_grid = np.repeat(dlons[np.newaxis, :], len(lats), axis=0)

    az = np.arctan2(dlons_grid, (L1 - L2_grid))
    return np.rad2deg(az)


def _isometric_latitude(lat, ellipsoid):
    e = _eccentricity(ellipsoid)
    term1 = np.tan((np.pi / 4) + (lat / 2))
    sinlat = np.sin(lat)
    num = 1 - e * sinlat
    denom = 1 + e * sinlat
    term2 = (num / denom) ** (e / 2)
    isometric_latitude = np.log(term1 * term2)
    return isometric_latitude


def _pixel_positions(dem, transform):
    height, width = dem.shape
    dlat, dlon = _resolution(transform)
    w, s, e, n = rasterio.transform.array_bounds(height, width, transform)
    # centers of pixels
    lats = np.linspace(n - dlat / 2, s + dlat / 2, height)
    lons = np.linspace(w + dlon / 2, e - dlon / 2, width)
    return lats, lons


def _dist2deg(latlon, distance, ellipsoid):
    lat, _ = latlon

    a = ellipsoid.semi_major_metre
    e = _eccentricity(ellipsoid)
    # in meter / radian
    rad_lat = np.deg2rad(lat)
    dlat = a * (1.0 - e ** 2) / (1.0 - e ** 2 * math.sin(rad_lat) ** 2) ** (3 / 2)
    dlon = a * math.cos(rad_lat) / math.sqrt(1.0 - e ** 2 * math.sin(rad_lat) ** 2)

    distance_eps = 0.01
    latMin = 0
    latMax = np.deg2rad(10)
    lonMin = 0
    lonMax = np.deg2rad(10)

    while True:
        delta_lat = (latMin + latMax) / 2
        delta_lon = (lonMin + lonMax) / 2

        dist_var_lat = dlat * delta_lat
        dist_var_lon = dlon * delta_lon

        if (
            abs(dist_var_lat - distance) < distance_eps
            and abs(dist_var_lon - distance) < distance_eps
        ):
            break
        if dist_var_lat < distance:
            latMin = delta_lat
        else:
            latMax = delta_lat
        if dist_var_lon < distance:
            lonMin = delta_lon
        else:
            lonMax = delta_lon

    return (np.rad2deg(delta_lat), np.rad2deg(delta_lon))


def _eccentricity(ellipsoid):
    f = 1 / ellipsoid.inverse_flattening
    e = math.sqrt(2 * f - f ** 2)
    return e


def _resolution(transform):
    return -transform[4], transform[0]
