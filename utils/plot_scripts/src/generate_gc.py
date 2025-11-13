import numpy as np 
import sys 

def generate_gc(lon1, lat1, lon2, lat2, num_points=100, radius=6371.0):
    """
    Generate a great-circle path between two points using only NumPy.

    Parameters:
    - lon1, lat1: coordinates of the first point (in degrees)
    - lon2, lat2: coordinates of the second point (in degrees)
    - num_points: number of interpolated points
    - radius: sphere radius in meters (default: Earth's mean radius)

    Returns:
    - lons: longitudes along the path (degrees)
    - lats: latitudes along the path (degrees)
    - dists: distances from the first point (meters)
    """
    # Convert to radians
    lon1, lat1 = np.radians([lon1, lat1])
    lon2, lat2 = np.radians([lon2, lat2])

    # Convert to Cartesian coordinates
    def sph2cart(lon, lat):
        return np.array([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat)
        ])

    A = sph2cart(lon1, lat1)
    B = sph2cart(lon2, lat2)

    # Angle between A and B
    omega = np.arccos(np.clip(np.dot(A, B), -1.0, 1.0))

    if omega == 0:
        # Points are the same
        lats = np.full(num_points, np.degrees(lat1))
        lons = np.full(num_points, np.degrees(lon1))
        dists = np.zeros(num_points)
        return lons, lats, dists

    # Interpolation fractions
    f = np.linspace(0, 1, num_points)

    # Slerp (spherical linear interpolation)
    sin_omega = np.sin(omega)
    points = (np.sin((1 - f) * omega)[:, None] * A + np.sin(f * omega)[:, None] * B) / sin_omega

    # Convert back to lat/lon
    lats = np.degrees(np.arcsin(points[:, 2]))
    lons = np.degrees(np.arctan2(points[:, 1], points[:, 0]))

    # Distance from first point along the arc
    dists = f * omega * radius

    return lons, lats, dists

def main():
    if len(sys.argv) !=9:
        print("Usage: python generate_gc.py lon1 lat1 lon2 lat2 zmin zmax n outfile")
        print("Example: python generate_gc.py 100.0 30.0 110.0 35.0 0.0 80.0 128 gc.dat")
        print(len(sys.argv))
        sys.exit(1)

    lon1 = float(sys.argv[1])
    lat1 = float(sys.argv[2])
    lon2 = float(sys.argv[3])
    lat2 = float(sys.argv[4])
    zmin = float(sys.argv[5])
    zmax = float(sys.argv[6])
    n = int(sys.argv[7])
    outfile = sys.argv[8]

    lons, lats, dists = generate_gc(lon1, lat1, lon2, lat2, n)
    z = np.linspace(zmin, zmax, n)
    data = np.zeros((n*n, 4))
    for i in range(n):
        for j in range(n):
            idx = i * n + j
            data[idx, 0] = lons[j]
            data[idx, 1] = lats[j]
            data[idx, 3] = dists[j]
            data[idx, 2] = z[i]

    np.savetxt(outfile, data, fmt="%f")

if __name__ == "__main__":
    main()