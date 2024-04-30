import numpy as np
import matplotlib.pyplot as plt 

# input
x0,y0 = 35.3,99.7
nlon,nlat = 36,33
dx,dy = 0.3,0.3
nlevelx = 0.6
nlevely = 0.6
nlevelz = 0.4
zstr = " 0.0000  2.0000  4.0000  6.0000  8.0000 10.0000 12.0000 14.0000 16.0000 18.0000 20.0000 24.0000 28.0000 32.0000 36.0000 40.0000 45.0000 50.0000 55.0000 60.0000 70.0000 80.0000 90.0000 100.0000 110.0000"
use_sph = True
shift_depth = False

# write velocity
z = list(map(lambda x: float(x),zstr.split()))
nz = len(z)
v = np.zeros((nz,nlon,nlat))
for i in range(nz):
    v[i,:,:] = 2.5 + 0.02 * z[i]

f = open("MOD.in","w")
f.write("%d %d %d\n" %(nlat,nlon,nz))
f.write("%f %f\n" %(x0,y0))
f.write("%f %f\n" %(dx,dy))
for i in range(nz):
    f.write("%f " %(z[i]))
f.write("\n")
for i in range(nz):
    for j in range(nlon):
        for k in range(nlat):
            f.write("%f "%(v[i,j,k]))
        f.write("\n")
f.close()

xx = np.arange(nlat-2) + 1
yy = np.arange(nlon-2) + 1
zz = np.arange(nz-2) + 1

vtrue = v.copy()
for i in range(nz-2):
    for j in range(nlon-2):
        for k in range(nlat-2):
            v0 = 0.1 * np.sin(xx[k] * nlevelx) * np.sin(yy[j] * nlevely) * np.sin(zz[i] * nlevelz)
            vtrue[i+1,j+1,k+1] = v[i+1,j+1,k+1] * (1.0 + v0)

# plot 
plt.figure(1,figsize=(14,7))
plt.subplot(1,2,1)
plt.contourf(x0 - np.arange(nlat) * dx,-np.array(z),(vtrue[:,3,:] - v[:,3,:]) / v[:,3,:] * 100)
plt.colorbar(orientation='horizontal')
plt.subplot(1,2,2)
plt.contourf(x0 - np.arange(nlat) * dx,y0 + np.arange(nlon) * dx,(vtrue[3,:,:] - v[3,:,:]) / v[3,:,:] * 100)
plt.colorbar(orientation='horizontal')
plt.savefig("veloc3d.jpg")

f = open("MOD.true","w")
f.write("%d %d %d\n" %(nlat,nlon,nz))
f.write("%f %f\n" %(x0,y0))
f.write("%f %f\n" %(0.3,0.3))
for i in range(nz):
    f.write("%f " %(z[i]))
f.write("\n")
for i in range(nz):
    for j in range(nlon):
        for k in range(nlat):
            f.write("%f "%(vtrue[i,j,k]))
        f.write("\n")
f.close()

f = open("mod_true.dat","w")
for i in range(nz):
    for j in range(nlon):
        for k in range(nlat):
            f.write("%f %f %f %f\n"%(y0 + dy * j, x0 - dx * k,z[i],vtrue[i,j,k]))
        f.write("\n")
f.close()