import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from rasterio.windows import from_bounds


def render_atmospheric_horizon(
    dem_path, observer_latlon, observer_height=1.7, max_dist_km=50
):
    obs_lat, obs_lon = observer_latlon
    max_dist_m = max_dist_km * 1000

    # Constants for Earth
    R_EARTH = 6371000  # meters
    # Standard Refraction coefficient (makes the earth look slightly flatter)
    # usually 0.13, effectively increasing R by ~1.15x.
    # Using 1.0 for geometric or 1.15 for optical approximation.
    REFRACTION_FACTOR = 1.15
    R_EFFECTIVE = R_EARTH * REFRACTION_FACTOR

    with rasterio.open(dem_path) as src:
        # 1. Calculate a Window (Bounding Box) to avoid reading the whole file
        # Approx: 1 degree lat ~ 111km.
        deg_buffer = (max_dist_km / 111.0) * 1.2  # 20% safety margin

        window = from_bounds(
            obs_lon - deg_buffer,
            obs_lat - deg_buffer,
            obs_lon + deg_buffer,
            obs_lat + deg_buffer,
            src.transform,
        )

        # Read data and transform
        dem_data = src.read(1, window=window)
        dem_transform = src.window_transform(window)
        nodata = src.nodata if src.nodata is not None else -9999

    # 2. Create Coordinate Grids
    height, width = dem_data.shape
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))

    # Convert pixels to Lat/Lon
    xs, ys = rasterio.transform.xy(dem_transform, rows, cols, offset="center")
    lons = np.array(xs).reshape(dem_data.shape)
    lats = np.array(ys).reshape(dem_data.shape)

    # 3. Convert to Meters (Local approximation)
    # We use the Cosine of the Observer Latitude to scale Longitude
    lat_scale = 111132.0  # Meters per degree latitude
    lon_scale = 111132.0 * np.cos(np.radians(obs_lat))  # Meters per degree longitude

    dy = (lats - obs_lat) * lat_scale
    dx = (lons - obs_lon) * lon_scale

    dist_m = np.sqrt(dx**2 + dy**2)

    # 4. Filter Data
    # Mask out: Nodata, Distance > 50km, and Distance == 0 (Observer)
    mask = (dem_data != nodata) & (dist_m <= max_dist_m) & (dist_m > 1)

    valid_dist = dist_m[mask]
    valid_dem = dem_data[mask]
    valid_dx = dx[mask]
    valid_dy = dy[mask]

    # 5. Calculate Angles

    # A. Azimuth (0 to 360)
    azimuths_rad = np.arctan2(valid_dx, valid_dy)  # dy is North
    azimuths_deg = (np.degrees(azimuths_rad) + 360) % 360

    # B. Earth Curvature Drop
    # Drop = Distance^2 / (2 * R)
    curvature_drop = (valid_dist**2) / (2 * R_EFFECTIVE)

    # C. Relative Height
    # Get observer elevation from the pixel at his location (or close to it)
    try:
        row_obs, col_obs = src.index(obs_lon, obs_lat)
        # Handle case where observer is outside the loaded window (rare due to buffer)
        # Simple fallback: average of local window if exact pixel is tricky inside this
        # logic
        obs_base_elev = np.median(dem_data[dem_data != nodata])  # Fallback
        # Ideally: read the specific pixel value for observer Z
    except Exception:
        obs_base_elev = 0

    # Adjust Target Height: Real Z - Curvature Drop - (Obs Z + Eye Level)
    rel_height = valid_dem - curvature_drop - (obs_base_elev + observer_height)

    # D. Elevation Angle
    elev_angles = np.degrees(np.arctan(rel_height / valid_dist))

    # 6. Extract the Skyline (Argmax)
    # We want the MAX angle for every azimuth bin, BUT we need the DISTANCE of that max
    # angle.

    # Binning: 3600 bins (0.1 degree precision)
    bins = np.linspace(0, 360, 3601)
    bin_indices = np.digitize(azimuths_deg, bins)

    horizon_angles = []
    horizon_dists = []
    horizon_azimuths = []

    # Iterate bins (This loop is fast enough for 3600 items)
    for i in range(1, len(bins)):
        # Indices of pixels in this azimuth slice
        in_bin = bin_indices == i

        if np.any(in_bin):
            # Get angles in this bin
            bin_angles = elev_angles[in_bin]
            bin_dists = valid_dist[in_bin]

            # Find the index of the max angle
            max_idx = np.argmax(bin_angles)

            horizon_angles.append(bin_angles[max_idx])
            horizon_dists.append(bin_dists[max_idx])
            horizon_azimuths.append(bins[i - 1])
        else:
            # No terrain in this direction (or masked out)
            horizon_angles.append(-90)  # Below ground
            horizon_dists.append(max_dist_m)
            horizon_azimuths.append(bins[i - 1])

    return np.array(horizon_azimuths), np.array(horizon_angles), np.array(horizon_dists)


# ==========================================
# RUNNING THE TOOL
# ==========================================

# 1. Configuration
tif_file = "/Users/guilhem/Documents/projects/dtm/dem_wgs84_b.tif"  # Or local path
observer = (45.902437, 6.144651)  # (Lat, Lon) - Chamonix area example
eye_level = 1.7  # meters

# 2. Compute
print("Computing visibility...")
az, el, dists = render_atmospheric_horizon(
    tif_file, observer, observer_height=eye_level
)

# 3. Plotting with Atmospheric Shading
print("Plotting...")
plt.figure(figsize=(15, 6), facecolor="white")
ax = plt.gca()

# Set Sky Color
sky_color = "#87CEEB"  # Light Sky Blue
ax.set_facecolor(sky_color)

# Normalize distance for the colormap
# Close = 0 (Dark), Far = 1 (Light/Blue)
norm = plt.Normalize(vmin=0, vmax=50000)  # 0 to 50km

# Choose a colormap that goes from Dark Grey/Green to Sky Blue/White
# 'copper' is nice for silhouettes, 'Blues_r' is good for fog, 'gist_earth' is realistic
cmap = cm.get_cmap("gray")

# SCATTER PLOT
# We use scatter to color every degree individually based on distance
# s=15 ensures the dots merge into a line
plt.scatter(az, el, c=dists, cmap=cmap, norm=norm, s=15, edgecolors="none", zorder=10)

# FILL BELOW
# To make it look solid, we can fill simply with black or dark gray,
# or repeat the scatter logic vertically (computationally expensive).
# Simple solution: Fill black, but keep the colored scatter on top as the "edge".
plt.fill_between(az, el, -90, color="#202020", zorder=5)

# Styling
plt.colorbar(label="Distance to Horizon (meters)")
plt.xlim(0, 360)
plt.ylim(min(max(el), -5), max(el) + 2)  # Crop to interesting area
plt.xlabel("Azimuth (Degrees)")
plt.ylabel("Elevation Angle (Degrees)")
plt.title("Horizon View with Atmospheric Depth (Max 50km)")
plt.grid(True, alpha=0.2, linestyle="--")

plt.tight_layout()
plt.show()
