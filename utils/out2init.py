import numpy as np
import sys 

def main():
    if len(sys.argv) !=4 :
        print("python out2init.py outfile MOD.init MOD.out")
        print("example: python out2init.py results/mod_iter10.dat MOD MOD.out")
        exit(1)
    
    infile = sys.argv[1]
    refile = sys.argv[2]
    outfile = sys.argv[3]

    # load ref model
    f = open(refile,"r")
    line = f.readline()
    info = line.split()
    nx,ny,nz = int(info[0]),int(info[1]),int(info[2])
    line = f.readline()
    info = line.split()
    lat0,lon0 = float(info[0]),float(info[1])
    line = f.readline()
    info = line.split()
    dx,dy = float(info[0]),float(info[1])
    line = f.readline()
    info = line.split()
    z = np.zeros((nz))
    for i in range(nz):
        z[i] = float(info[i])
    f.close()

    # load output model
    d = np.loadtxt(infile)[:,3]
    d = d.reshape(nz,ny,nx)

    # write model
    f = open(outfile,"w")
    f.write("%d %d %d           #nlat nlon nz\n"%(nx,ny,nz))
    f.write("%g %g          # lat0 lon0\n" %(lat0,lon0))
    f.write("%g %g          # dlat dlon\n"%(dx,dy))
    for iz in range(nz):
        f.write("%g "%(z[iz]))
    f.write("\n")
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                f.write("%g "%(d[iz,iy,ix]))
            f.write("\n")
    f.close()

if __name__ == "__main__":
    main()