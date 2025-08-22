Manual of JDSurfG
==================
<a id="toc"></a>

## Table of Contents
- [1. Introduction](#intro)
- [2. Preliminaries](#preliminaries)
- [3. Installation](#install)
- [4. Direct Surface Wave Tomography Module](#dsurf)
  - [4.1 Parameter file](#dsurf-param)
  - [4.2 Initial Model and True Model Files](#dsurf-model-file)
  - [4.3 Dispersion Data File](#dsurf-data-file)
  - [4.4 Run This Module](#dsurf-run)
  - [4.5 Output Files](#dsurf-output)
- [5. Joint Inversion Module](#joint)
  - [5.1 Parameter File](#joint-param)
  - [5.2 Gravity Data File](#joint-data-file)
  - [5.3 Gravity Matrix File](#joint-gravity-matrix)
  - [5.4 Gravity Reference Model File](#joint-ref-model)
  - [5.5 Run This Module](#joint-run)
  - [5.6 Output Files](#joint-output)
- [6. Gravity Module](#gravity-module)
  - [6.1 Run This Module](#gravity-run)
  - [6.2 Output Files](#gravity-output)
- [7. Advanced Options](#advanced-options)
  - [7.1 OpenMP Mode](#openmp-mode)
  - [7.2 Empirical Relations](#empirical-relations)
  - [7.3 Random Seed](#random-seed)
- [8. Tips and Tools](#tips-and-tools)
  - [8.1 Warnings](#warnings)
  - [8.2 L-curve Analysis](#l-curve-analysis)
  - [8.3 Checkpoint Restart](#checkpoint-restart)

 <a id="intro"></a>
## 1. Introduction

This document provides an overview of how to use the **JDSurfG** package to perform 3-D joint inversion of shear wave velocity in the crust and upper mantle using surface wave dispersion and gravity anomaly data.

The package is composed of three independent modules:

1. Generation of the gravity matrix in spherical coordinates  
2. Direct surface wave tomographic inversion  
3. Joint inversion of dispersion and gravity anomaly data  

In the joint inversion process, we apply the direct surface wave tomography method ([Fang et al., 2015](https://academic.oup.com/gji/article/201/3/1251/759155)) to compute the 3-D (pseudo) sensitivity kernel of surface wave travel times. For the gravity forward matrix, we use the adaptive Gauss–Legendre integration method ([Li et al., 2011](https://www.sciencedirect.com/science/article/pii/S0926985111000206)).

Using empirical relationships between seismic velocity and rock density ([Brocher, 2005](https://pubs.geoscienceworld.org/ssa/bssa/article/95/6/2081/146858/Empirical-Relations-between-Elastic-Wavespeeds-and)), this package serves as a valuable tool for investigating the 3-D shear wave structure of the crust and upper mantle. For further technical details, please refer to [Du et al., 2021](https://academic.oup.com/gji/article/227/3/1961/6333361).

This package is written in **C++** and **Fortran**.  
The C++ components handle core logic tasks, including I/O, data structures, and interfaces.  
The Fortran modules—adapted from the [DSurfTomo](https://github.com/HongjianFang/DSurfTomo/tree/stable/src) with significant modifications (e.g., parallelization, group velocity ray tracing, and analytical derivatives)—are used to compute surface wave travel times and the associated 1D/3D sensitivity kernels.

> **Note:**  
> This package has only been tested by a limited number of users under similar computational environments. As such, there may still be bugs. I welcome bug reports and suggestions, and will make improvements accordingly. Over time, I plan to introduce additional features to enhance the package.

[Back to Table of Contents](#toc)

 <a id="preliminaries"></a>

## 2. Preliminaries
This package uses the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ library to handle multi-dimensional arrays. You need to install Eigen on your system before compiling this package.

Make sure your **C++ compiler supports C++11** standards (e.g., **GCC ≥ 4.8**).

The package is built using the [CMake](https://cmake.org/) build system. Ensure that your `cmake` version is **≥ 3.1.0**.

[Back to Table of Contents](#toc)

 <a id="install"></a>

## 3. Installation

After downloading the code from GitHub, you can compile it using the following commands:

```bash
cd JDSurfG
mkdir build
cd build
cmake .. -DCXX=g++ -DFC=gfortran -DEIGEN_INC=/path/to/eigen
make -j4
```
This will generate four executable files in the bin/ directory:
```bash 
mkmat DSurfTomo JointSG syngrav
```
If these files appear, the installation was successful. To check some functionalities of the package, you can build with testing enabled:
```bash
cmake .. -DCXX=g++ -DFC=gfortran -DEIGEN_INC=/path/to/eigen \
         -DBUILD_TEST=TRUE
make -j4
```
[Back to Table of Contents](#toc)

 <a id="dsurf"></a>

## 4. Direct Surface Wave Tomography Module

This module requires 3 input files (or 4 if you're running a checkerboard or other synthetic test):

- `DSurfTomo.in`: Contains information about the input model, dispersion data, and inversion parameters.
- `surfdataSC.dat`: Dispersion data for each station pair.
- `MOD`: Initial model.
- `MOD.true`: True model (required only for synthetic tests).

Before proceeding to the file formats, please consider the following:

1. The input model is assumed to be evenly spacing in latitude and longitude, but grid spacing may vary in depth.
2. The input model should be slightly larger than the target region to allow for **B-spline smoothing** during forward computation. For a target region `[lat0,lat1],[lon0,lon1],[0,depth]`,   
the input model should at least cover  
`[lat0-dlat,lat1+dlat]`,`[lon0-dlon,lon1+dlon]`,
`[0,depth+dz]`
3. Stations should be positioned sufficiently far from the model boundaries to avoid wavefronts traveling along edges, which may lead to warnings when running the program.

---

 <a id="dsurf-param"></a>

### 4.1 Parameter File (`DSurfTomo.in`)

This is a self-explanatory text file. Here's a sample template:

```ini
# DSurfTomo inversion Parameters

# Number of iterations
NITERS = 16
ITER_CURRENT = 0  # Current iteration index

# Velocity range in km/s
MIN_VELOC = 2.0
MAX_VELOC = 5.0

# Synthetic test settings
SYN_TEST = 1           # 0: real data, 1: synthetic
NOISE_LEVEL = 1.5      # Std. dev. of Gaussian noise

# Inversion method: 0 = LSMR, 1 = CG, 2 = LBFGS
INV_METHOD = 2

# LSMR options
SMOOTH = 15.0          # 2nd-order Tikhonov regularization
DAMP = 0.01            # Damping for LSMR
NTHREADS = 2           # Number of threads

# CG/LBFGS smoothing
SMOOTH_IN_KM = 0       # If 1, smoothing in km; else in grid units
SIGMA_H = 0.25         # Horizontal smoothing parameter
SIGMA_V = 0.25         # Vertical smoothing parameter
ITER_START = 0         # Start iteration for using previous model info

# Line search settings
MAX_REL_STEP = 0.04    # Max relative change allowed per iteration
```
Lines beginning with # are treated as comments and ignored. You may freely add your own comments. Here is the parameters description below:

- **`NITERS`**  
  Number of iterations to run.

- **`ITER_CURRENT`**  
  Current iteration of the initial model. If set to `10`, the initial model will be indexed as `mod_iter10` in the result directory. Useful for continuing from a previous run.

- **`MIN_VELOC`, `MAX_VELOC`**  
  Minimum and maximum permitted shear wave velocities (in km/s). These serve as prior constraints on the velocity model.

- **`SYN_TEST`**  
  Whether to perform a synthetic test. If enabled (`1`), the `MOD.true` file should be provided.

- **`NOISE_LEVEL`**  
  Gaussian noise level added to synthetic data. Defined as:  
  ```math 
  d_{obs}^i = d_{syn}^i + \text{noiselevel} \times \mathcal{N}(0, 1) 
  ```

- **`INV_METHOD`**  
  Select the inversion method:
  - `0`: LSMR (recommended for small-scale problems)
  - `1`: CG (only for test)
  - `2`: L-BFGS (recommended for large-scale problems)

- **`SMOOTH`, `DAMP`**  
  Second and first-order Tikhonov regularization coefficients, used only in LSMR mode.

- **`NTHREADS`**  
  Number of threads used in the LSMR solver. Typically 2–4 is sufficient.

- **`SMOOTH_IN_KM`**  
  If `1`, smoothing parameters are in kilometers; otherwise, they are index-based. Used for smoothing in L-BFGS/CG mode (not LSMR).

- **`SIGMA_H`, `SIGMA_V`**  
  Gaussian smoothing parameters for horizontal and vertical directions. The search direction is smoothed using:  
  ```math
  \bar{g}(\mathbf{x}) = \frac{1}{W(\mathbf{x})} \int_V g(\mathbf{x}) e^{-\frac{1}{2} \mathbf{x}^T \Sigma^{-1} \mathbf{x}} \, dV
  ```
  where the normalization factor is:  
  ```math
  W(\mathbf{x}) = \int_V e^{-\frac{1}{2} \mathbf{x}^T \Sigma^{-1} \mathbf{x}} \, dV
  ```
  and  
  ```math
  \Sigma = \text{Diag}[\sigma_h^2, \sigma_h^2, \sigma_v^2]
  ```
  If `SMOOTH_IN_KM = 1`, then ( $\sigma_h$, $\sigma_v$ ) are in kilometers and $\mathbf{x}$ represents spatial coordinates; otherwise, they are in grid indices.

- **`ITER_START`**  
  Use only the previous models starting from this iteration in L-BFGS/CG mode. This must match the current iteration when testing new smoothing parameters.

- **`MAX_REL_STEP`**  
  In L-BFGS/CG mode, this sets the maximum relative change allowed between iterations during line search. A value of `0.03` or `0.04` is typically sufficient.

 <a id="dsurf-model-file"></a>

### 4.2 Initial Model and True Model Files

> ⚠️ **Note:** The format of `MOD.true` in the current version is **different** from previous versions!

The **initial model** and **true model** files are named `MOD` and `MOD.true`, respectively.  
- `MOD` is used as the initial input model.  
- `MOD.true` is required for synthetic tests (e.g., checkerboard test). The format is below:

- **Line 1:**  
   Three integers specifying the number of grid points in the **latitude**, **longitude**, and **depth** directions:


- **Line 2:**  
Coordinates of the origin point (Northwest corner), in degrees:

- **Line 3:**  
Grid intervals in the latitudinal and longitudinal directions (both in degrees):


- **Line 4:**  
Depth values (in km) of each grid point. There should be exactly `nz` depth values.

- **Remaining lines:**  
Shear wave velocity values stored in **row-major** order as `V[nz, nlon, nlat]`.  
- `V[0, 0, 0]` corresponds to the **top Northwest** point at depth `dep[0]`
- `V[nz-1, nlon-1, nlat-1]` corresponds to the **bottom Southeast** point at depth `dep[nz-1]`, it should be written into this file by:
```python
# Python script to write velocity model to 'MOD' file
with open("MOD", "w") as f:
 for i in range(nz):
     for j in range(nlon):
         for k in range(nlat):
             f.write(f"{V[i, j, k]} ")
         f.write("\n")
```

 <a id="dsurf-data-file"></a>
### 4.3 Dispersion Data File

The dispersion data file (e.g., `example/surfdataSC.dat`) contains surface wave dispersion data (distance/travel-time). Below is an excerpt of this file:
```ini
1 19 # Rayleigh Phase
4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40
1 19 # Rayleigh Group
4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40
0 # Love Phase
0 # Love Group
# 31.962600 108.646000 0 1 2 0
33.732800 105.764100 3.0840
34.342500 106.020600 3.1154
33.357400 104.991700 3.0484
33.356800 106.139500 3.0820
34.128300 107.817000 3.1250
33.229300 106.800200 3.0575
# 29.905000 107.232500 0 1 2 0
34.020000 102.060100 2.7532
```

This file contains **three types of information**:
1. **Data Header**: Lists wave type and period information.
2. **Master Station (Source)**: Lines beginning with `#` that define the source location and mode details.
3. **Receiver Stations**: Following lines (not starting with `#`) list receiver locations and dispersion values.

---
The data header defines the number of modes and their corresponding period points.  
Each wave type has its own section. Example: Using fundamental and first modes of Rayleigh phase velocity with periods [4, 7.5, 10] and [3, 4] respectively:
```ini
2 3 2 # Rayleigh phase (2 modes: 3 and 2 periods)
4 7.5 10
3 4
```

You may also define other wave and velocity types:
If **no data** is available for a wave type, use `0` in its line.

---

Each source station is denoted by a comment line starting with `#`. The format is: 
```ini
latitude longitude mode period-index wavetype velotype
```

- `mode`: Mode index (0 = fundamental, 1 = first mode, etc.)
- `wavetype`: 2 = Rayleigh, 1 = Love
- `velotype`: 0 = phase velocity, 1 = group velocity
- `period-index`: Index of the period (starting from 1)

Example:  
If your dataset has fundamental mode Rayleigh phase velocity for periods 4s, 5s, 6s, 7s and this entry corresponds to 5s, then the values would be:
```ini
lat lon 0 2 2 0
```
---
Receiver stations follow the master station entry. Each line contains:
```ini
latitude longitude dispersion_value
```
<a id="dsurf-run"></a>

### 4.4 Run This Module
After you have prepared all the required files in a folder,
then you could enter into this folder,and run this command
in your shell:
```bash
./DSurfTomo -h
```
Then you could type in all required files according to the 
order on the screen.

<a id="dsurf-output"></a>

### 4.5 Output Files
This module outputs one folder: `results/`. It contains all intermediate models (`mod_iter*`) and data files (`res*.dat`) for each iteration. There may also be temporary binary files (`*.bin`) that store the optimization history of NLCG/L-BFGS optimization.

The format of the model file is:
```
longitude latitude depth Vs
```
The format of the traveltime file is (9 columns):
```
distance observation synthetic lat1 lon1 lat2 lon2 period wavetype

```

[Back to Table of Contents](#toc)

<a id="joint"></a>

## 5. Joint Inversion Module

This module requires 4–6 input files:

- `JointSG.in`: Contains information about the input model, dispersion data, and inversion parameters.
- `surfdataSC.dat`: Surface wave dispersion data.
- `obsgrav.dat`: Gravity anomaly dataset.
- `gravmat.dat`: Sensitivity kernel of gravity to density.
- `MOD`: Initial model.
- `MOD.true`: Checkerboard model.
- `MOD.ref`: Gravity reference model, used to compute gravity anomaly.

The formats of each file are described below.

<a id="joint-param"></a>

## 5.1 Parameter File

The `JointSG.in` file is similar to `DSurfTomo.in` with only minor differences. Below is an example of its content:

```ini
# synthetic test 
NOISE_LEVEL_SWD  =  1.5  # swd noise level
NOISE_LEVEL_GRAV =  5.   # gravity noise level 
WEIGHT_SWD       =  1.5
WEIGHT_GRAV      =  5.
RELATIVE_P       =  0.5  # relative factor
```
The weighted cost function for the joint inversion problem, as described in [Julia et al., 2000](https://academic.oup.com/gji/article/143/1/99/894713), is:
```math
 L = \frac{p}{N_1 \sigma_1^2}\sum_{i=1}^{N_1} (d_{1,i}- 
        d_{1,i}^o)^2 +   
        \frac{1-p}{N_2 \sigma_2^2} \sum_{i=1}^{N_2}
        (d_{2,i}- d_{2,i}^o)^2
```
where $\mathbf{d}_{1,2}$ is the data vector for each dataset. $N_{1,2}$
is the size of each dataset and $\sigma_{1,2}$ (The `WEIGTH_SWD` and `WEIGHT_GRAV` in `JointSG.in`) is the corresponding
uncertainties. `p` (`RELATIVE_P`) is a parameter in the range [0,1] to control the contribution
of each dataset manually.

The additional parameter is `REMOVE_AVG`. When this integer is equal to 1, the program will automatically 
remove the average value of observed and synthetic gravity data. This option is applied to remove large-scale density anomalies in the deep mantle.

<a id="joint-data-file"></a>

### 5.2 Gravity Data File
This file uses a simple format to store gravity data. Each line records the longitude, latitude, and gravity anomaly for a specific point.

Below is an example showing the first nine lines of the file:

```text
100.000000 35.000000 -224.206921
100.000000 34.950000 -224.110699
100.000000 34.900000 -241.774731
100.000000 34.850000 -240.245968
100.000000 34.800000 -234.187542
100.000000 34.750000 -246.670605
100.000000 34.700000 -238.575939
100.000000 34.650000 -241.357250
100.000000 34.600000 -260.059429
```
<a id="joint-gravity-matrix"></a>

### 5.3 Gravity Matrix File

This file is generated by the gravity module. For more details, refer to [Gravity module](#6-gravity-module)

<a id="joint-ref-model"></a>

### 5.4 Gravity Reference Model File

A gravity anomaly reflects density perturbations relative to a normal ellipsoid. To focus on local anomalies in the research area, absolute anomalies must be converted to relative anomalies before inversion by removing the average value of observed anomalies. For real data inversion, the following three steps are used to compute relative anomalies:

- Calculate the density distribution of the current S-wave model and the reference model using empirical relationships.
- Compute the density anomaly and apply the gravity matrix to obtain gravity anomalies.
- Subtract the average value from the synthetic anomaly data.

The reference model can be user-defined or set to the default, which is the average of the initial model at each depth. When processing real data, ensure the `REMOVE_AVG` parameter is set to `1` in the `JointSG.in` file.

<a id="joint-run"></a>

### 5.5 Run This Module

Once all required files are prepared in a directory, navigate to that directory and execute the following command in your terminal:

```bash
./JointSG -h
```
Then, provide the required files in the order prompted on the screen.

<a id="joint-output"></a>
### 5.6 Output Files

This module generates output files in the `resultsJ/` directory.

The format of files in the `kernel/` directory is identical to those described in [Surface Wave Output Files](#surface-wave-output-files). The `results/` directory includes the following for each iteration:
- Models (`joint_mod_iter*`)
- Surface wave travel times (`res_surf*.dat`)
- Gravity anomalies (`res_grav*.dat`)

The gravity anomaly files use the following format:
```
lon lat observed-anomalies synthetic-anomalies
```

[Back to Table of Contents](#toc)

<a id="gravity-module"></a>

## 6. Gravity Module

This module requires two input files:

- `MOD`: Only the header and depth information are needed.
- `obsgrav.dat`: Contains the coordinates of each observation point.

These file formats are described in previous chapters.

<a id="gravity-run"></a>

### 6.1 Run This Module

Once all required files are prepared in a directory, navigate to that directory and execute the following command in your terminal:

```bash
./mkmat -h
```
Then, input the required files in the order prompted on the screen.

<a id="gravity-output"></a>

### 6.2 Output Files 

This module generates a single output file: `gravmat.bin`. This file contains the derivative matrix for gravity anomalies, stored in Compressed Sparse Row (CSR) format. To understand how to read this file, refer to the code in `src/shared/csr_matrix.cpp`.


[Back to Table of Contents](#toc)


<a id="advanced-options"></a>

## 7. Advanced Options

<a id="openmp-mode"></a>

### 7.1 OpenMP Mode

All four programs support execution on a single node with `OpenMP`. To specify the number of threads, set the following environment variable before running any program:

```bash
export OMP_NUM_THREADS=8
```
The maximum number of threads depends on your platform. To determine this, use the following command:
```bash
cat /proc/cpuinfo | grep cores | wc -l
```
The maximum number of threads is `half` the output of this command.

<a id="empirical-relations"></a>

### 7.2 Empirical Relations
The program uses empirical relations to link rock density and velocity. To customize these relations for your study area, modify the function in `src/SWD/exmpirical.cpp`:
```cpp
void empirical_relation(float vsz, float &vpz, float &rhoz)
{
    vpz = 0.9409 + 2.0947*vsz - 0.8206*pow(vsz,2) +
          0.2683*pow(vsz,3) - 0.0251*pow(vsz,4);
    rhoz = 1.6612*vpz - 0.4721*pow(vpz,2) +
           0.0671*pow(vpz,3) - 0.0043*pow(vpz,4) +
           0.000106*pow(vpz,5);
}
```
If you have defined your own function, ensure you also update the derivative function in the same file.

<a id="random-seed"></a>

### 7.3 Random Seed
To ensure consistent random numbers across independent numerical experiments, the random seed is fixed in `src/utils/gaussian.cpp`:
```cpp
static default_random_engine e(11);
```
You can change 11 to another integer if needed.

<a id="tips-and-tools"></a>


[Back to Table of Contents](#toc)

## 8. Tips and Tools

<a id="warnings"></a>

### 8.1 Warnings
Both DSurfTomo and JointSG rely on two modules: (1) forward computation of dispersion curves from a 1-D model and (2) ray-tracing on a 2-D spherical surface. Errors may arise due to incompatible input model formats, poor initial models, or inappropriate station locations. The modules provide on-screen messages or output files to assist with debugging. The first warnings is like:
```bash
Warning: "WARNING: improper initial value in disper-no zero found"
```
This warning appears on the screen, and a fort.66 file is generated in the running directory. It indicates that the dispersion module cannot compute dispersion curves for a problematic 1-D model. Possible causes include:

- First iteration issues: The initial model format may not match the required format. Compare the model in fort.66 with your input model.
- Inappropriate regularization: Unsuitable Tikhonov regularization parameters (e.g., damping or smoothing factors) may cause abrupt model changes.

Another warning is like :
```bash
Warning: "Source lies outside bounds of model (lat,long)="
```
This warning occurs when receivers are too close to the model’s boundary, causing ray paths to align with it. To resolve this, enlarge the model to ensure receivers are sufficiently within the bounds.

<a id="l-curve-analysis"></a>

### 8.2 L-curve Analysis

L-curve analysis is commonly used to identify optimal smoothing and damping factors for inverse problems. In our numerical experiments, the damping factor has a smaller impact on the final result compared to the smoothing factor. Therefore, we recommend setting the damping factor to a small value (e.g., `0.001`) during inversion and focusing on tuning the smoothing factor.

A Python script, `utils/lcurve.py`, is provided in the `utils/` directory to assist with tuning the smoothing factor. This script outputs the roughness of the inverted model, calculated as $||L(m - m_0)||_2$, and the model difference, $||m - m_0||_2$.

<a id="checkpoint-restart"></a>

## 8.3 Checkpoint Restart

To start an inversion from an intermediate result (e.g., `mod_iter5.dat`) instead of the initial model, use the Python script `utils/out2init.py`. Run it as follows:

```bash
python out2init.py mod_iter5.dat MOD05
```
This script converts the output model to the initial model format. It is particularly useful for joint inversion and for restarting from a checkpoint if issues occur.


[Back to Table of Contents](#toc)
