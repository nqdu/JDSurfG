import numpy as np 

def main():
    import sys 
    if len(sys.argv) != 8:
        print("Usage: python generate_plane.py lonmin lonmax latmin latmax depth n outfile")
        print("Example: python generate_plane.py 99.7 109.2 26.5 35.3  10 128 plane.dat")
        sys.exit(1)

    lonmin = float(sys.argv[1])
    lonmax = float(sys.argv[2])
    latmin = float(sys.argv[3])
    latmax = float(sys.argv[4])
    depth = float(sys.argv[5])
    n = int(sys.argv[6])
    outfile = sys.argv[7]

    # Generate grid points
    lons = np.linspace(lonmin, lonmax, n)
    lats = np.linspace(latmin, latmax, n) 
    data = np.zeros((n*n, 3))
    for i in range(n):
        for j in range(n):
            idx = i * n + j
            data[idx, 0] = lons[j]
            data[idx, 1] = lats[i]
            data[idx, 2] = depth
    np.savetxt(outfile, data, fmt="%f")

if __name__ == "__main__":
    main()