import numpy as np
import matplotlib.pyplot as plt 
import sys 


def main():
    if len(sys.argv) != 3:
        print("Usage: ./plot.py wavetype(Rc,Rg,Lc,Lg) Tid")
        print("Usage: python plot_kernel.py Rc 6")
        exit(1)

    wavetype = sys.argv[1]
    nt = int(sys.argv[2])
    
    data = np.loadtxt(f"output/disper_{wavetype}_{nt}.txt")
    data = data[:data.shape[0]-1,:]

    plt.figure(1,figsize=(14,7))
    plt.subplot(131)
    plt.title("vs")
    plt.plot(data[:,1],-data[:,0],label='analytical')
    plt.plot(data[:,2],-data[:,0],label='fd')
    plt.legend()

    plt.subplot(132)
    plt.title("rho")
    plt.plot(data[:,3],-data[:,0],label='analytical')
    plt.plot(data[:,4],-data[:,0],label='fd')
    plt.legend()

    plt.subplot(133)
    plt.title("vp")
    plt.plot(data[:,5],-data[:,0],label='analytical')
    plt.plot(data[:,6],-data[:,0],label='fd')

    plt.legend()

    if wavetype == "Rc":
        plt.suptitle("Rayleigh Phase")
    elif wavetype == "Rg":
        plt.suptitle("Rayleigh Group")
    elif wavetype == "Lc":
        plt.suptitle("Love Phase")
    else:
        plt.suptitle("Love Group")

    plt.savefig("kernel.jpg")

main()

