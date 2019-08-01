Thank you for downloading the Uncertainty framework. 

Created by: Michal Polic

-----------------------------------------------------
Requirements: 
 - Gflags, Eigen, Blas and Lapack, Ceres, Cuda, Magma dense

Optional:
 - Glog, Matlab, AliceVision, SuiteSparse
-----------------------------------------------------

## Installation on Windows

This tutorial shows how to create the project in Visual Studio by Cmake. At first, you need to download the repository. You can use any GIT client for that. We recommend the GitExtensions which is simple and intuitive.
1) Install libraries dependencies:
- Download the precompiled libraries
You can download the prebuild libraries for win64 at: 'cmp.felk.cvut.cz/~policmic/files/uncer_libs_win64.zip' and unpack it into the root project directory.
- Or install them through VCPKG
- Or build them manually
2) Generate the Visual Studio project and compile it
3) Try it by calling
 uncertatinty.exe -in=./in/01_Cube.jacob


## Installation on Linux

1) Install Cmake
2) Install all the required libraries

 sudo apt-get install cmake
 sudo apt-get install build-essential
 sudo apt-get install libboost-all-dev
 sudo apt-get install libeigen3-dev
 sudo apt-get install libsuitesparse-dev
 sudo apt-get install libfreeimage-dev
 sudo apt-get install libgoogle-glog-dev  # may cause problems in MEX ( unresolved symbol _Unwind.. )
 sudo apt-get install libgflags-dev
 sudo apt-get install libglew-dev
 sudo apt-get install qtbase5-dev
 sudo apt-get install libqt5opengl5-dev
 sudo apt-get install libatlas-base-dev
 sudo apt-get install libsuitesparse-dev
 git clone https://ceres-solver.googlesource.com/ceres-solver
 cd ceres-solver
 mkdir build
 cd build
 cmake ..
 sudo make install

Download and compile Magma (http://icl.cs.utk.edu/magma/software/index.html)
Set up the root of the source code and the directory with compiled libraries. It is not necessary to install magma.

3) Generate the project (Makefiles)
 mkdir build
 cd build
 cmake .. -DMAGMA_ROOT=/path/to/magma/sources
4) Try the framework by calling
 uncertainty -in=./Release/in/01_Cube.jacob


-----------------------------------------------------

Output:
The covariance matrices are symmetric. Therfore, our output is just the upper triangle of their covariance matrices.
So, the camera which is in our model represented by 9 parameters has upper triangle composed of 45 values and the 
covariances of points are represented by 6 values:
-----------
| 1  2  3 |
| 2  4  5 |
| 3  5  6 |
----------- 
