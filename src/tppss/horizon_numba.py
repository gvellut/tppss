import logging
import math

from numba import njit, prange
import numpy as np
from pyproj import CRS
import rasterio

logger = logging.getLogger(__name__)

KM = 1000.0
MILES = 1609.34


def horizon(latlon, raster, distance=25 * KM, precision=1, height=0):
    """
    Compute horizon profile from a given location using DEM data.
    
    This is the main entry point - keeps the same interface as the original.
    Internal computations are accelerated with numba.
    """
    if raster.crs.is_projected:
        raise Exception("Only geographic CRS are supported")

    crs = CRS.from_wkt(raster.crs.to_wkt())
    ellipsoid = crs.ellipsoid
    
    # Extract ellipsoid parameters once for numba functions
    semi_major = ellipsoid.semi_major_metre
    inverse_flattening = ellipsoid.inverse_flattening

    logger.debug("Extracting data...")
    study_area, study_transform = _extract_data(latlon, raster, distance, ellipsoid)

    lats, lons = _pixel_positions(study_area, study_transform)

    y_obs, x_obs = rasterio.transform.rowcol(study_transform, latlon[1], latlon[0])

    logger.debug("Computing azimuths...")
    # Azimuth of all pixel corners seen from observer point
    # bottom right corners
    lat_res, lon_res = _resolution(study_transform)
    
    # Pre-compute eccentricity
    e = _eccentricity_from_params(inverse_flattening)
    
    azimuth = _azimuth_numba(
        latlon[0], latlon[1], 
        lats - lat_res / 2, 
        lons + lon_res / 2, 
        e
    )

    # Corresponding elevation angle
    # add observer height
    z_obs = study_area[y_obs, x_obs] + height

    logger.debug("Computing elevations...")
    elevation = _elevation_angle_numba(
        z_obs,
        study_area,
        latlon[0],
        latlon[1],
        lats,
        lons,
        semi_major,
        e,
    )

    logger.debug(f"Computing sun mask... precision={precision}")
    helevation = _compute_mask_numba(
        x_obs, y_obs, study_area, precision, azimuth, elevation
    )

    # Finalization
    # Set horizon
    # TODO negative values should be allowed (when located at peak of mountain)
    helevation[helevation < 0] = 0
    hzenith = 90 - helevation
    # Store azimuth from 0 to 360 (for using with solar azimuth functions)
    hazimuth = np.linspace(0, 360, len(helevation))

    return helevation, hzenith, hazimuth


@njit(cache=True)
def _eccentricity_from_params(inverse_flattening):
    """Compute eccentricity from inverse flattening."""
    f = 1.0 / inverse_flattening
    e = math.sqrt(2.0 * f - f * f)
    return e


@njit(cache=True)
def _isometric_latitude_numba(lat, e):
    """Compute isometric latitude for a single value."""
    term1 = math.tan((math.pi / 4.0) + (lat / 2.0))
    sinlat = math.sin(lat)
    num = 1.0 - e * sinlat
    denom = 1.0 + e * sinlat
    term2 = (num / denom) ** (e / 2.0)
    return math.log(term1 * term2)


@njit(cache=True, parallel=True)
def _azimuth_numba(lat_obs, lon_obs, lats, lons, e):
    """
    Compute azimuth of all pixel corners seen from observer point.
    Returns azimuth in degrees.
    """
    lat_obs_rad = math.radians(lat_obs)
    lon_obs_rad = math.radians(lon_obs)
    
    nlats = len(lats)
    nlons = len(lons)
    
    # Pre-compute isometric latitudes for all lats
    L1 = _isometric_latitude_numba(lat_obs_rad, e)
    
    lats_rad = np.empty(nlats)
    L2 = np.empty(nlats)
    for i in prange(nlats):
        lats_rad[i] = math.radians(lats[i])
        L2[i] = _isometric_latitude_numba(lats_rad[i], e)
    
    lons_rad = np.empty(nlons)
    dlons = np.empty(nlons)
    for j in prange(nlons):
        lons_rad[j] = math.radians(lons[j])
        dlons[j] = lon_obs_rad - lons_rad[j]
    
    az = np.empty((nlats, nlons))
    for i in prange(nlats):
        dL = L1 - L2[i]
        for j in range(nlons):
            az[i, j] = math.degrees(math.atan2(dlons[j], dL))
    
    return az


@njit(cache=True)
def _geographic2cartesian_single(lat, lon, h, semi_major, e):
    """Convert geographic to cartesian for a single point."""
    sinlat = math.sin(lat)
    coslat = math.cos(lat)
    N = semi_major / math.sqrt(1.0 - (e * e) * (sinlat * sinlat))
    
    x = (N + h) * math.cos(lon) * coslat
    y = (N + h) * math.sin(lon) * coslat
    z = (N * (1.0 - e * e) + h) * sinlat
    
    return x, y, z


@njit(cache=True, parallel=True)
def _elevation_angle_numba(
    z_obs, study_area, lat_obs, lon_obs, lats, lons, semi_major, e
):
    """
    Compute elevation angle for all pixels.
    """
    lat_obs_rad = math.radians(lat_obs)
    lon_obs_rad = math.radians(lon_obs)
    
    height, width = study_area.shape
    
    # Observer cartesian coordinates
    x_A, y_A, z_A = _geographic2cartesian_single(
        lat_obs_rad, lon_obs_rad, z_obs, semi_major, e
    )

    # Pre-compute cos/sin for observer
    cos_lon_obs = math.cos(lon_obs_rad)
    sin_lon_obs = math.sin(lon_obs_rad)
    cos_lat_obs = math.cos(lat_obs_rad)
    sin_lat_obs = math.sin(lat_obs_rad)
    
    # Pre-compute lat/lon in radians
    lats_rad = np.empty(height)
    lons_rad = np.empty(width)
    for i in range(height):
        lats_rad[i] = math.radians(lats[i])
    for j in range(width):
        lons_rad[j] = math.radians(lons[j])
    
    alpha = np.empty((height, width))
    
    for i in prange(height):
        lat_rad = lats_rad[i]
        sinlat = math.sin(lat_rad)
        coslat = math.cos(lat_rad)
        N = semi_major / math.sqrt(1.0 - (e * e) * (sinlat * sinlat))
        
        for j in range(width):
            lon_rad = lons_rad[j]
            h = study_area[i, j]
            
            # Cartesian coordinates of point B
            x_B = (N + h) * math.cos(lon_rad) * coslat
            y_B = (N + h) * math.sin(lon_rad) * coslat
            z_B = (N * (1.0 - e * e) + h) * sinlat
            
            # Vector AB
            dx = x_B - x_A
            dy = y_B - y_A
            dz = z_B - z_A
            
            # Scalar product with normal at A
            inner_product = (
                dx * cos_lon_obs * cos_lat_obs
                + dy * sin_lon_obs * cos_lat_obs
                + dz * sin_lat_obs
            )
            
            # Norm of AB
            norm = math.sqrt(dx * dx + dy * dy + dz * dz)
            
            # Elevation angle in degrees
            if norm > 0:
                alpha[i, j] = math.degrees(math.asin(inner_product / norm))
            else:
                alpha[i, j] = 0.0
    
    return alpha


@njit(cache=True)
def _digitize_scalar(value, bins):
    """
    Find the index where value would be inserted in bins (sorted ascending).
    Equivalent to np.digitize for a single value.
    """
    n = len(bins)
    if value < bins[0]:
        return 0
    if value >= bins[n - 1]:
        return n
    
    # Binary search
    lo = 0
    hi = n
    while lo < hi:
        mid = (lo + hi) // 2
        if bins[mid] <= value:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True)
def _fill_elevation_NE(
    elevationNE, azNE, azimuthNE, elevation, x_obs,
    north_distance, east_distance, length_elevation
):
    """Fill elevation matrix for NE sector."""
    for isoline in range(north_distance):
        for az_idx in range(length_elevation):
            az_val = azNE[az_idx]
            k = _digitize_scalar(az_val, azimuthNE[isoline, x_obs - 1:])
            if k != 0 and k != east_distance + 1:
                elevationNE[isoline, az_idx] = elevation[isoline, x_obs - 1 + k]


@njit(cache=True)
def _fill_elevation_NW(
    elevationNW, azNW, azimuthNW, elevation, x_obs,
    north_distance, west_distance, length_elevation
):
    """Fill elevation matrix for NW sector."""
    for isoline in range(north_distance):
        for az_idx in range(length_elevation):
            az_val = azNW[az_idx]
            k = _digitize_scalar(az_val, azimuthNW[isoline, :x_obs + 1])
            if k != 0 and k != west_distance + 1:
                elevationNW[isoline, az_idx] = elevation[isoline, k - 1]


@njit(cache=True)
def _fill_elevation_E(
    elevationE, azE, azimuth, elevation, x_obs, east_distance, height, length_elevation
):
    """Fill elevation matrix for E sector."""
    for isoline in range(east_distance):
        col_idx = x_obs + isoline
        for az_idx in range(length_elevation):
            az_val = azE[az_idx]
            k = _digitize_scalar(az_val, azimuth[:, col_idx])
            if k != 0 and k != height:
                elevationE[isoline, az_idx] = elevation[k, col_idx]


@njit(cache=True)
def _fill_elevation_S(
    elevationS, azS, azimuth, elevation, y_obs, south_distance, width, length_elevation
):
    """Fill elevation matrix for S sector."""
    for isoline in range(south_distance):
        row_idx = y_obs + isoline
        for az_idx in range(length_elevation):
            az_val = azS[az_idx]
            # Reverse the row for digitize
            k = _digitize_scalar(az_val, azimuth[row_idx, ::-1])
            if k != 0 and k != width:
                elevationS[isoline, az_idx] = elevation[row_idx, width - 1 - k]


@njit(cache=True)
def _fill_elevation_W(
    elevationW, azW, azimuth, elevation, west_distance, height, length_elevation
):
    """Fill elevation matrix for W sector."""
    for isoline in range(west_distance):
        for az_idx in range(length_elevation):
            az_val = azW[az_idx]
            # Reverse the column for digitize
            k = _digitize_scalar(az_val, azimuth[::-1, isoline])
            if k != 0 and k != height:
                elevationW[isoline, az_idx] = elevation[height - 1 - k, isoline]


@njit(cache=True)
def _max_along_axis0(arr):
    """Compute max along axis 0."""
    nrows, ncols = arr.shape
    result = np.empty(ncols)
    for j in range(ncols):
        max_val = arr[0, j]
        for i in range(1, nrows):
            if arr[i, j] > max_val:
                max_val = arr[i, j]
        result[j] = max_val
    return result


@njit(cache=True)
def _compute_helevation(az, azNtoN, sun_mask, total_length_elevation):
    """Compute final horizon elevation array."""
    helevation = np.zeros(total_length_elevation)
    
    for i in range(total_length_elevation):
        az_val = az[i]
        max_val = -np.inf
        for j in range(len(azNtoN)):
            if abs(azNtoN[j] - az_val) < 1e-10:  # floating point comparison
                if sun_mask[j] > max_val:
                    max_val = sun_mask[j]
        if max_val > -np.inf:
            helevation[i] = max_val
    
    return helevation


def _compute_mask_numba(x_obs, y_obs, study_area, precision, azimuth, elevation):
    """
    Compute the horizon mask using numba-accelerated inner loops.
    """
    # Elevation vector length (degrees)
    length_elevation = precision * 90 + 1
    height, width = study_area.shape

    # Specific azimuth values for NE (-180 to -90) and NW (90 to 180) areas
    azimuthNE = azimuth.copy()
    azimuthNE[:y_obs, x_obs - 1] = azimuthNE[:y_obs, x_obs - 1] - 360
    azimuthNW = azimuth.copy()
    azimuthNW[:y_obs, x_obs] = azimuthNW[:y_obs, x_obs] + 360

    # Initialization
    north_distance = y_obs
    east_distance = width - x_obs
    south_distance = height - y_obs
    west_distance = x_obs

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

    # Main computation using numba functions
    if north_distance > 0 and east_distance > 0:
        _fill_elevation_NE(elevationNE, azNE, azimuthNE, elevation, x_obs, 
                          north_distance, east_distance, length_elevation)
    
    if north_distance > 0 and west_distance > 0:
        _fill_elevation_NW(elevationNW, azNW, azimuthNW, elevation, x_obs,
                          north_distance, west_distance, length_elevation)
    
    if east_distance > 0:
        _fill_elevation_E(elevationE, azE, azimuth, elevation, x_obs,
                         east_distance, height, 2 * length_elevation - 1)
    
    if south_distance > 0:
        _fill_elevation_S(elevationS, azS, azimuth, elevation, y_obs,
                         south_distance, width, 2 * length_elevation - 1)
    
    if west_distance > 0:
        _fill_elevation_W(elevationW, azW, azimuth, elevation,
                         west_distance, height, 2 * length_elevation - 1)

    # Max for each angle (2nd dimension) for each sunmask
    if north_distance > 0:
        sun_maskNE = _max_along_axis0(elevationNE)
        sun_maskNW = _max_along_axis0(elevationNW)
    else:
        sun_maskNE = np.zeros(length_elevation)
        sun_maskNW = np.zeros(length_elevation)
    
    if east_distance > 0:
        sun_maskE = _max_along_axis0(elevationE)
    else:
        sun_maskE = np.zeros(2 * length_elevation - 1)
    
    if south_distance > 0:
        sun_maskS = _max_along_axis0(elevationS)
    else:
        sun_maskS = np.zeros(2 * length_elevation - 1)
    
    if west_distance > 0:
        sun_maskW = _max_along_axis0(elevationW)
    else:
        sun_maskW = np.zeros(2 * length_elevation - 1)

    # Global azimuth (North to North) and sun mask (elevation angle)
    azNtoN = np.concatenate((azNE, azE, azS, azW, azNW))
    sun_mask = np.concatenate((sun_maskNE, sun_maskE, sun_maskS, sun_maskW, sun_maskNW))

    total_length_elevation = precision * 360 + 1
    
    # Corresponding azimuth (from -180° to 180°)
    az = np.linspace(-180, 180, total_length_elevation)
    helevation = _compute_helevation(az, azNtoN, sun_mask, total_length_elevation)

    return helevation


def _extract_data(latlon, raster, distance, ellipsoid):
    """Extract relevant DEM data around the observation point."""
    lat, lon = latlon
    lat_resolution, lon_resolution = _resolution(raster.transform)

    row, col = raster.index(lon, lat)

    if not 0 <= row < raster.height or not 0 <= col < raster.width:
        raise Exception("LatLon not covered by DEM")

    # deg corresponding to distance
    lat_distance, lon_distance = _dist2deg(latlon, distance, ellipsoid)
    topo_lat_distance = round(lat_distance / lat_resolution)
    topo_lon_distance = round(lon_distance / lon_resolution)

    window_base = rasterio.windows.Window(0, 0, raster.width, raster.height)
    window = rasterio.windows.Window(
        col - topo_lon_distance,
        row - topo_lat_distance,
        topo_lon_distance * 2,
        topo_lat_distance * 2,
    )

    window = window.intersection(window_base)
    study_area = raster.read(1, window=window)
    study_transform = raster.window_transform(window)

    return study_area, study_transform


def _pixel_positions(dem, transform):
    """Get lat/lon positions for pixel centers."""
    height, width = dem.shape
    dlat, dlon = _resolution(transform)
    w, s, e, n = rasterio.transform.array_bounds(height, width, transform)
    # centers of pixels
    lats = np.linspace(n - dlat / 2, s + dlat / 2, height)
    lons = np.linspace(w + dlon / 2, e - dlon / 2, width)
    return lats, lons


@njit(cache=True)
def _dist2deg_numba(lat, distance, semi_major, e):
    """
    Convert distance to degrees for lat and lon.
    """
    rad_lat = math.radians(lat)
    sin_lat = math.sin(rad_lat)
    cos_lat = math.cos(rad_lat)
    
    e2 = e * e
    
    dlat = semi_major * (1.0 - e2) / ((1.0 - e2 * sin_lat * sin_lat) ** 1.5)
    dlon = semi_major * cos_lat / math.sqrt(1.0 - e2 * sin_lat * sin_lat)

    distance_eps = 0.01
    latMin = 0.0
    latMax = math.radians(10.0)
    lonMin = 0.0
    lonMax = math.radians(10.0)

    while True:
        delta_lat = (latMin + latMax) / 2.0
        delta_lon = (lonMin + lonMax) / 2.0

        dist_var_lat = dlat * delta_lat
        dist_var_lon = dlon * delta_lon

        if (abs(dist_var_lat - distance) < distance_eps and 
            abs(dist_var_lon - distance) < distance_eps):
            break
        
        if dist_var_lat < distance:
            latMin = delta_lat
        else:
            latMax = delta_lat
        
        if dist_var_lon < distance:
            lonMin = delta_lon
        else:
            lonMax = delta_lon

    return math.degrees(delta_lat), math.degrees(delta_lon)


def _dist2deg(latlon, distance, ellipsoid):
    """Convert distance to degrees."""
    lat, _ = latlon
    semi_major = ellipsoid.semi_major_metre
    e = _eccentricity_from_params(ellipsoid.inverse_flattening)
    return _dist2deg_numba(lat, distance, semi_major, e)


def _resolution(transform):
    """Get resolution from transform."""
    return -transform[4], transform[0]
