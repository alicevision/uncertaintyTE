Thank you for downloading the Uncertainty framework. 

Created by: Michal Polic
License agreements: you can do whatever you want.

-----------------------------------------------------
Requirements: 
 - Magma -> Magma is a mathematical library for computation on GPUs. 
         To install Magma, you will have to install: Fortran compiler, Cuda libraries, BLAS, LPACK 
 - Eigen 

-----------------------------------------------------

Installation on Windows: 
This tutorial shows how to create the project in Visual Studio by Cmake. At first, you need to download the repository. You can use any GIT client for that. We recommend the GitExtensions which is simple and intuitive. 
1) Install GitExtensions (https://gitextensions.github.io) 
      Download the project from internet: right click to the directory -> select GitExt Clone
      Select repository to clone: git@bitbucket.org:policmic/uncertatinty.git and confirm Clone.
2) Download and install the last version of Cmake ( https://cmake.org/download )
3) Open Cmake (cmake-gui)
   Select the source code, 'Browse Source...', as the root directory of the project
   Select the path to the binaries, 'Browse Build...', as a subdirectory 
   (e.g. <project_dir>/build), run 'Configure' and run 'Generate'
4) Open the project '<project_dir>/build/Uncertatinty.sln'
   Right click at 'uncertainty' and select 'Build'
5) After build will be executable in './build/[Release,Debug]/uncertatinty.exe'
   Try it by calling 'uncertatinty.exe 2 in/01_Cube.jacob' in command line


Installation on Linux:
1) Install Cmake from command line: 
  > sudo apt-get install cmake
  > sudo apt-get install ccmake (NOT REQUIRED BUT RECOMMENDED)
  > sudo apt-get install cmake-qt-gui (NOT REQUIRED BUT RECOMMENDED)
2) Install BLAS and LAPACK:
  > sudo apt-get install libblas-dev checkinstall
  > sudo apt-get install libblas-doc checkinstall
  > sudo apt-get install liblapacke-dev checkinstall
  > sudo apt-get install liblapack-doc checkinstall
3) Install Cuda development toolkit. 
  Download from: 
    https://developer.nvidia.com/cuda-downloads
  Follow the instructions:
    https://askubuntu.com/questions/799184/how-can-i-install-cuda-on-ubuntu-16-04
4?) The library may require MKL library. At least it requires the variable MKLROOT.  Skip this step and try return here if building of the project fail.
(The MKL lib can be download here https://software.intel.com/en-us/mkl. Follow provided steps in MKL README. After the installation add environment variables by running script '/opt/intel/compilers_and_libraries_xxxx/linux/bin/compilervars.sh' calling 'sh compilervars.sh intel64'.)
5) Few codes should be added to '~/.bashrc'. You can use any text editor, e.g. nano 'nano ~/.bashrc' and add at the end: 
  source /opt/intel/compilers_and_libraries_2017.4.196/linux/bin/compilervars.sh intel64
  export CUDADIR=/usr/local/cuda-8.0
  export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64:/usr/local/magma/lib:$LD_LIBRARY_PATH
  export MKLROOT=/opt/intel/compilers_and_libraries_2017.4.196/linux/mkl
  export PATH=/usr/local/cuda-8.0/bin:${MKLROOT}/lib/intel64:/usr/local/magma/lib:$PATH
Adjust the paths and according your system.
6) Create project using Cmake.
  Run 'cmake ..' in subdirectory build (you should create it) in this project
7) Try the framework by calling 'uncertatinty 2 Release/in/01_Cube.jacob'
