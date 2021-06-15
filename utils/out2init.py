import sys 
import numpy as np 
def main():
    if len(sys.argv) != 3:
        print("please run like this:")
        print("this.py outfile inputfile")
        exit()
    infile = sys.argv[1]
    outname = sys.argv[2]

    # load files and get z components
    d = np.loadtxt(infile)
    nlon = np.unique(d[:,0]).size
    nlat = np.unique(d[:,1]).size
    z:np.ndarray = np.unique(d[:,2])
    nz = z.size

    # loop to save files
    f = open(outname,'w')
    for i in range(nz):
        f.write("%f "%(z[i]))
    f.write('\n')
    for i in range(nz):
        for j in range(nlon):
            for k in range(nlat):
                idx = i * nlat * nlon + j * nlat + k
                f.write('%f '%(d[idx,3]))
            f.write('\n')
    f.close()

if __name__ == "__main__":
    main()
