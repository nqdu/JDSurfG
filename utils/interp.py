import numpy as np
from scipy.interpolate import griddata 
from numba import jit 
from numpy import sin,cos,arctan
import sys 

def gcpoints(lon,lat,lon1,lat1,n):
    '''
    get great circle points with end points (lon,lat),(lon1,lat1)
    return n points in degree
    
    Parameters:
        lon,lat, lon1,lat1 : end points of a great circle
        n                  : number of points
        
    Returns:
        x,y : lon and lat for great circle, in degree
    '''
    theta0,phi0,theta1,phi1=list(map(lambda x:np.deg2rad(x),[lat,lon,lat1,lon1]))

    if phi0 > phi1:
        phi0,phi1,theta0,theta1=phi1,phi0,theta1,theta0

    if phi0 == phi1:
        theta = np.linspace(theta0,theta1,n)
        phi = np.ones(n,)*phi0
    else:
        # construct a plane equation ax+by+z=0,solve a and b
        d = np.array([-sin(theta0),-sin(theta1)]).reshape(2,1)
        A = np.array([[cos(theta0)*cos(phi0),cos(theta0)*sin(phi0)],[cos(theta1)*cos(phi1),cos(theta1)*sin(phi1)]])
        m = np.linalg.inv(A)@d
        a = m[0,0]
        b = m[1,0]

        phi=np.linspace(phi0,phi1,n)
        theta= arctan(-a*cos(phi)-b*sin(phi))

    
    phi,theta=list(map(lambda x: np.rad2deg(x),[phi,theta]))

    return phi,theta

def locations2degrees(lat1, long1, lat2, long2):
    lat1, lat2, long1, long2 = np.broadcast_arrays(lat1, lat2, long1, long2)
    lat1 = np.radians(np.asarray(lat1))
    lat2 = np.radians(np.asarray(lat2))
    long1 = np.radians(np.asarray(long1))
    long2 = np.radians(np.asarray(long2))
    long_diff = long2 - long1

    gd = np.degrees(
        np.arctan2(
        np.sqrt((
        np.cos(lat2) * np.sin(long_diff)) ** 2 +
        (np.cos(lat1) * np.sin(lat2) - np.sin(lat1) *
        np.cos(lat2) * np.cos(long_diff)) ** 2),
        np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) *
        np.cos(long_diff)))

    return gd

def lonlat2xyz(lon,lat,dep):
    r = 6371.
    r2d = np.pi / 180
    x = (r - dep) * np.cos(lat*r2d) * np.cos(lon*r2d)
    y = (r - dep) * np.cos(lat*r2d) * np.sin(lon*r2d)
    z = (r - dep) * np.sin(lat*r2d)

    return x,y,z

def get_xyz(data,is_spherical):
    if not is_spherical:
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]

        return x,y,z 
    else:
        x,y,z = lonlat2xyz(data[:,0],data[:,1],data[:,2])

        return x,y,z

def get_distance(x1,y1,x2,y2,is_spher):
    if is_spher:
        dist = 6371 * np.deg2rad(locations2degrees(y1,x1,y2,x2))
    else:
        dist = np.hypot(x1-x2,y1-y2)
    
    return dist 

def load_topo(topo_file):
    topo = np.loadtxt(topo_file,skiprows=2)
    f = open(topo_file,"r")
    line = f.readline()
    nz,nx = map(lambda x:int(x),line.split())
    line = f.readline()
    xmin,xmax,zmin,zmax = map(lambda x:float(x),line.split())
    n = nz * nx

    xtopo = np.zeros((n))
    ztopo = np.zeros((n))

    for i in range(nz):
        for j in range(nx):
            id = i * nx + j 
            xtopo[id] = xmin + (xmax - xmin) / (nx-1) * j
            ztopo[id] = zmin + (zmax - zmin) / (nz-1) * i
    
    topo *= -0.001
    return xtopo,ztopo,topo

def main():
    if len(sys.argv) != 7:
        print("Usage: ./this x1 y1 x2 y2 zmax model_name")
        exit(1)
    x1,y1,x2,y2,zmax = map(lambda x:float(x),sys.argv[1:6])
    model_name = sys.argv[6]
    print(model_name)
    
    # get spherical flag
    is_spherical = True

    # load model
    model = np.loadtxt(model_name)
    model_x,model_y,model_z = get_xyz(model[:,:3],is_spherical)

    # profile
    n = 64
    profx = np.zeros((n))
    profy = np.zeros((n))
    dist = profx.copy()
    if not is_spherical:
        profx = np.linspace(x1,x2,n)
        profy = np.linspace(y1,y2,n)
    else:
        profx,profy = gcpoints(x1,y1,x2,y2,n)

    for i in range(0,n):
        dist[i] = get_distance(profx[i],profy[i],profx[0],profy[0],is_spherical)

    # get coordinate    
    cord = np.zeros((n*n,3))
    for i in range(n):
        for j in range(n):
            id = i * n + j
            cord[id,0] = profx[i]
            cord[id,1] = profy[i]
            cord[id,2] = zmax / (n-1) * j

    # get velocity
    cordx,cordy,cordz = get_xyz(cord,is_spherical)
    v = griddata((model_x,model_y,model_z),model[:,-1],(cordx,cordy,cordz))
    idx = np.where(np.isnan(v))[0]
    v[idx] = griddata((model_x,model_y,model_z),model[:,-1],(cordx[idx],cordy[idx],cordz[idx]),'nearest')

    f = open("out_v.txt","w")
    for i in range(n):
        for j in range(n):
            id = i * n + j
            f.write("%f %f %f\n"%(dist[i],cord[id,2],v[id]))
    f.close()

main()