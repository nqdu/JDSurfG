import numpy as np
from numba import jit 
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import sys

@jit(nopython=True)
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

@jit(nopython=True)
def convert(md,nx,ny,nz):
    n = (nx-2) * (ny-2) * (nz - 1)
    x = np.zeros((n),dtype=float)
    for i in range(nz):
        for j in range(ny):
            for k in range(nx):
                idx = i * (ny-2) * (nx-2) + (j-1) * (nx-2) + k-1 
                if i<nz-1 and j>=1 and j<=ny-2 and  k>=1 and k<=nx-2:
                    #print(idx)
                    x[idx] = md[i*ny*nx+j*nx+k]
    return x

def get_laplacian(nx,ny,nz):
    # generate csr type laplacian matrix 

    n = (nx-2) * (ny-2) * (nz - 1)
    nonzeros = n * 7
    indptr = np.zeros((n+1),dtype=int)
    indices = np.zeros((nonzeros),dtype=int )
    data = np.zeros((nonzeros),dtype = float)
    nonzeros = add_regularization(nx,ny,nz,indices,indptr,data)
    indices = indices[:nonzeros]
    data = data[:nonzeros]
    smat = csr_matrix((data, indices, indptr),shape=(n,n))

    return smat


def compute_roughness(nx,ny,nz,smat,result_dir:str):
    # load model and compute m - m0
    md = np.loadtxt(result_dir + '/mod_iter1.dat') [:,-1]
    md -= np.loadtxt(result_dir + '/mod_iter0.dat') [:,-1]
    x = convert(md,nx,ny,nz)
    roughness_smooth = np.sqrt(np.sum((smat * x)**2))
    roughness_damp = np.sqrt(np.sum((x)**2))

    return roughness_smooth,roughness_damp

def main():
    if len(sys.argv) == 1:
        print("compute |m-m0|, |L(m-m0)| for mod_iter1.dat")
        print("Usage: python lcurve.py init_model result_dir")
        print("example: python lcurve.py MOD.init results/")

        exit(1)
    
    # get input args
    init_model_str = sys.argv[1]
    result_dir = sys.argv[2]

    # read model dimension/
    # read corner 
    f = open(init_model_str,"r")
    line = f.readline()
    nx,ny,nz = map(lambda x: int(x), line.split()[:3])
    f.close()

    # get laplacian
    smat = get_laplacian(nx,ny,nz)

    # get roughness
    r_smooth,r_damp = compute_roughness(nx,ny,nz,smat,result_dir)

    print("roughness smooth |L(m-m0)| = %g" %(r_smooth))
    print("roughness damp |(m-m0)| = %g" %(r_damp))
    
main()