import numpy as np 
from scipy.interpolate import griddata
import sys 
from numba import jit 

@jit(nopython=True)
def bsspline_coefs(w):
    cs4 = w*w*w/6.0
    cs1 = 1.0/6.0 + w*(w-1.0)/2.0 - cs4
    cs3 = w + cs1 - 2.0*cs4
    cs2 = 1.0 - cs1 - cs3 - cs4

    return [cs1, cs2, cs3, cs4]

@jit(nopython=True)
def bspline_interp_3d(x, y, z, x0, y0, z0, dx, dy, values):
    """
    B-spline interpolation in x/y direction and linear in z direction
    x, y: coordinates with homogeneous intervals
    z: 3D depth array z[nz,ny,nx] (can be different at each x,y location)
    x0, y0, z0: target interpolation point
    dx, dy: grid spacing in x and y
    values: 3D array of values to interpolate
    """
    nx = len(x)
    ny = len(y)
    nz = z.shape[0]
    
    # Find grid indices
    ix = -int((x0 - x[0]) / dx)
    iy = int((y0 - y[0]) / dy)
    
    # Handle edge cases in x/y
    if ix < 1:
        ix = 1
    if ix > nx - 3:
        ix = nx - 3
    if iy < 1:
        iy = 1
    if iy > ny - 3:
        iy = ny - 3
    
    # B-spline coefficients in x/y
    wx = -(x0 - x[ix]) / dx
    wy = (y0 - y[iy]) / dy
    coefx = bsspline_coefs(wx)
    coefy = bsspline_coefs(wy)
    
    # Interpolate using B-spline in x/y and linear in z
    result = 0.0
    
    for j in range(4):
        for i in range(4):
            idx_x = ix + i - 1
            idx_y = iy + j - 1
            
            # B-spline weight for this x,y position
            weight_xy = coefx[i] * coefy[j]
            
            # Find depth bracket at this specific (x,y) location
            iz = 0
            wz = 0.0
            
            if z0 <= z[0, idx_y, idx_x]:
                # Below minimum depth - use first layer
                iz = 0
                wz = 0.0
            elif z0 >= z[nz - 1, idx_y, idx_x]:
                # Above maximum depth - use last layer
                iz = nz - 2
                wz = 1.0
            else:
                # Find bracket at this (x,y) location
                for k in range(nz - 1):
                    if z0 >= z[k, idx_y, idx_x] and z0 <= z[k + 1, idx_y, idx_x]:
                        iz = k
                        dz = z[k + 1, idx_y, idx_x] - z[k, idx_y, idx_x]
                        if dz > 0:
                            wz = (z0 - z[k, idx_y, idx_x]) / dz
                        else:
                            wz = 0.0
                        break
            
            # Linear interpolation in z at this (x,y) location
            val_lower = values[iz, idx_y, idx_x]
            val_upper = values[iz + 1, idx_y, idx_x] if iz < nz - 1 else val_lower
            val_z = val_lower * (1.0 - wz) + val_upper * wz
            
            # Accumulate with B-spline weight
            result += weight_xy * val_z
    
    return result

@jit(nopython=True)
def interpolate(grid,orig):
    n = grid.shape[0]
    interp_values = np.zeros(n)

    # get lon/lat in original file
    lon = np.unique(orig[:, 0])
    lat = np.unique(orig[:, 1])[::-1]
    dlon = lon[1] - lon[0]
    dlat = abs(lat[1] - lat[0])

    # get depth
    nx = len(lat)
    ny = len(lon)
    nz = orig[:,0].size // (nx * ny)
    #print(nx,ny,nz)
    depth = np.zeros((nz, ny, nx), dtype=float)
    values = np.zeros((nz, ny, nx), dtype=float)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                idx = k * ny * nx + j * nx + i
                depth[k, j, i] = orig[idx, 2]
                values[k, j, i] = orig[idx, 3]

    # loop every point
    for i in range(n):
        lon0 = grid[i,0]
        lat0 = grid[i,1]
        dep0 = grid[i,2]

        interp_values[i] = bspline_interp_3d(lat, lon, depth, lat0, lon0, dep0, dlat, dlon, values)

        # # find nearest grid point
        # iy = int((lon0 - lon[0]) / dlon)
        # ix = int((lat0 - lat[0]) / dlat)

        # # check if ix/iy is out of bounds
        # if ix < 1 or ix > nx - 3 or iy < 1 or iy > ny - 3:
        #     interp_values[i] = np.nan
        #     continue

        # # get bspline coefficients
        # wy = ((lon0 - lon[ix]) / dlon)
        # wx = ((lat0 - lat[iy]) / dlat)
        # coefx = bsspline_coefs(wx)
        # coefy = bsspline_coefs(wy)

        # # find closest depth 
        # iz_loc = np.zeros((4,4), dtype=np.int32)
        # for jj in range(-1,3):
        #     for ii in range(-1,3):
        #         flag = False
        #         iz = -1
        #         for k in range(nz-1):
        #             if dep0 >= depth[k, iy+jj, ix+ii] and dep0 <= depth[k+1, iy+jj, ix+ii]:
        #                 iz = k
        #                 flag = True
        #                 break
        #         if not flag:
        #             if dep0 < depth[0, iy+jj, ix+ii]:
        #                 iz = -1
        #             else:
        #                 iz = nz - 1
        #         iz_loc[jj+1, ii+1] = iz
        
        # sums = 0.
        # for jj in range(4):
        #     for kk in range(4):
        #         iz = iz_loc[jj, kk]
        #         if iz == -1:
        #             val = values[0, iy+jj-1, ix+kk-1]
        #         elif iz == nz - 1:
        #             val = values[nz - 1, iy+jj-1, ix+kk-1]
        #         else:
        #             z1 = depth[iz+1, iy+jj-1, ix+kk-1]
        #             z0 = depth[iz, iy+jj-1, ix+kk-1]
        #             d1 = z1 - z0
        #             val1 = values[iz+1, iy+jj-1, ix+kk-1]
        #             val0 = values[iz, iy+jj-1, ix+kk-1]
        #             val = val0 + (dep0 - z0) * (val1 - val0) / d1
        #         sums += coefx[kk] * coefy[jj] * val
        
        # # assign value
        # interp_values[i] = sums

    return interp_values

def main():
    if len(sys.argv) != 4:
        print("Usage: python interp3d.py input_grid original_file output_file")
        print("Example: python interp3d.py grid.txt model.txt interp_data.txt")
        sys.exit(1)

    grid_file = sys.argv[1]
    original_file = sys.argv[2]
    output_file = sys.argv[3]

    # load file
    grid = np.loadtxt(grid_file)
    orig = np.loadtxt(original_file)

    # get grdid points and values
    points = orig[:, :3]
    values = orig[:, -1]

    if False:
        # griddata interpolation
        interp_values = griddata(points, values, grid[:, :3], method='linear')
    else:
        interp_values = interpolate(grid, orig)

    # handle NaN values (outside convex hull)
    nan_mask = np.isnan(interp_values)
    if np.any(nan_mask):
        interp_values[nan_mask] = griddata(points, values, grid[nan_mask, :3], method='nearest')

    # save output
    np.savetxt(output_file, np.hstack((grid, interp_values.reshape(-1, 1))), fmt="%f")

if __name__ == "__main__":
    main()