Thank you for downloading the Uncertainty framework. 

Created by: Michal Polic

-----------------------------------------------------
Requirements: 
 - Gflags, Eigen, Blas and Lapack, Ceres, Cuda, Magma dense

Optional:
 - Glog, Matlab, OpenMVG, SuiteSparse

-----------------------------------------------------

Installation on Windows: 
This tutorial shows how to create the project in Visual Studio by Cmake. At first, you need to download the repository. You can use any GIT client for that. We recommend the GitExtensions which is simple and intuitive. 
1) Download the precompiled libraries or compile the libs described in the Cmake
2) Build the project, compile it and try it by calling 'uncertatinty.exe 2 in/01_Cube.jacob' in command line


Installation on Linux:
1) Install Cmake
2) Create build directory in the project root and call there "cmake .."
3) Install all the required libraries (described in cmake output)
4) Download and compile Magma (http://icl.cs.utk.edu/magma/software/index.html)
   Set up the root of the source code and the directory with compiled libraries. It is not necessary to install magma.
5) Try the framework by calling 'uncertatinty 2 Release/in/01_Cube.jacob'


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