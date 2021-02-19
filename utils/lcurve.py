import numpy as np
from numba import jit 
from scipy.sparse import csr_matrix
from glob import glob
import matplotlib.pyplot as plt

@jit
def add_regularization(nx,ny,nz,indices,indptr,data):
    # get parameters required
    n = (nx-2) * (ny -2) * (nz  - 1) # model dimension
    nonzeros = len(data)
    m = 0

    count = 0
    nar = 0
    for k in range(nz-1):
        for j in range(ny-2):
            for i in range(nx-2):
                if i==0 or i == nx-3 or j==0 or j==ny-3 or k==0 or k==nz-2:
                    # and more restrictions to boundary points
                    #if nar + 1 > nonzeros:
                    #    exit()
                        #print("please increase sparse ratio!")
                        #exit()

                    rwc = count + m
                    data[nar] = 2.0
                    indices[nar] = k * (ny -2) * (nx -2) + j * (nx -2) + i
                    indptr[rwc + 1] = indptr[rwc]
                    nar += 1
                    count += 1
                    continue
                else:
                    #if nar  + 7 > nonzeros:
                        #print("please increase sparse ratio!")
                    #    exit()
                    rwc = count + m;  # current row
                    indptr[rwc +1] = indptr[rwc] + 7
                    clc = k * (ny -2) * (nx -2) + j * (nx -2) + i # current column
                    data[nar] = 6.0
                    indices[nar] = clc

                    # x direction
                    data[nar + 1] = -1.0
                    indices[nar + 1] = clc - 1
                    data[nar + 2] = -1.0
                    indices[nar + 2] = clc + 1

                    # y direction
                    data[nar + 3] = -1.
                    indices[nar + 3] = clc - (nx - 2)
                    data[nar + 4] = -1.
                    indices[nar + 4] = clc + (nx - 2)

                    # z direction
                    data[nar + 5] = -1.
                    indices[nar + 5] = clc - (nx - 2) * (ny - 2)
                    data[nar + 6] = -1.
                    indices[nar + 6] = clc + (nx - 2) * (ny - 2)

                    nar += 7
                    count += 1
    nonzeros = nar

    return nonzeros

@jit
def convert(md,x,nx,ny,nz):
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                idx = i * (ny-2) * (nx-2) + (j-1) * (nx-2) + k-1 
                if i<nz-1 and j>=1 and j<=ny-2 and  k>=1 and k<=nx-2:
                    #print(idx)
                    x[idx] = md[i*ny*nx+j*nx+k,3]
    

def main():
    nx = 33 
    ny = 36
    nz = 25
    n = (nx-2) * (ny-2) * (nz - 1)

    # generate csr_matrix
    nonzeros = n * 7
    indptr = np.zeros((n+1),dtype=np.int)
    indices = np.zeros((nonzeros),dtype=np.int )
    data = np.zeros((nonzeros),dtype = np.float32)
    nonzeros = add_regularization(nx,ny,nz,indices,indptr,data)
    indices = indices[:nonzeros]
    data = data[:nonzeros]
    smat = csr_matrix((data, indices, indptr),shape=(n,n))

    # load every folder
    filenames = glob("storage/*")
    d = np.zeros((len(filenames),3))

    for i,f in enumerate(filenames):
        # get smooth param
        smooth = float(f.split('results')[-1])
        d[i,0] = smooth

        # compute rms
        ds = np.loadtxt(f + "/res_surf1.dat")
        dg = np.loadtxt(f + "/res_grav1.dat")
        ns = ds.shape[0]
        ng = dg.shape[0]
        rms = np.sum((ds[:,1] - ds[:,2])**2 / 1.4**2)
        rmg = np.sum((dg[:,2] - dg[:,3])**2 * ns/ng / 10**2)
        rms = np.sqrt(rms + rmg)

        # compute Lm
        md = np.loadtxt(f + "/joint_mod_iter1.dat")
        md -= np.loadtxt(f + "/joint_mod_iter0.dat")
        x = np.zeros((n))
        convert(md,x,nx,ny,nz)
        y = np.sqrt(np.sum((smat * x)**2))

        d[i,1] = rms 
        d[i,2] = y 

        print(smooth)
    
    idx = np.argsort(d[:,0])
    d = d[idx,:]

    # save lcurve
    np.savetxt("lcurve.dat",d,fmt='%f')

    # draw
    plt.plot(d[:,1],d[:,2])
    plt.scatter(d[:,1],d[:,2])
    plt.show()
    

main()
