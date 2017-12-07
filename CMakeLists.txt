CMAKE_MINIMUM_REQUIRED (VERSION 3.5)
SET(CMAKE_VERBOSE_MAKEFILE ON)

# About the project
PROJECT (UncertaintyTE
  LANGUAGES C CXX
  VERSION 1.2.0)
#SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CONFIGURATION_TYPES "Debug;Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


MESSAGE(STATUS "---------------------------------------------------------")
MESSAGE(STATUS "TE-inverse uncertainty framework v.${PROJECT_VERSION}")
MESSAGE(STATUS "CMAKE_CONFIGURATION_TYPES: ${CMAKE_CONFIGURATION_TYPES}")
MESSAGE(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
MESSAGE(STATUS "---------------------------------------------------------")

# Guard against in-source builds
IF (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  MESSAGE(FATAL_ERROR "In-source builds not allowed.")
ENDIF()

# Add path to find_package modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules")



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# REQUIRED AND OPTIONAL LIBRARIES
# -----------------------------------------------------------------------------
MESSAGE(STATUS "REQUIRED LIBRARIES")
MESSAGE(STATUS " Gflags")
MESSAGE(STATUS " Eigen")
MESSAGE(STATUS " Blas and Lapack")
MESSAGE(STATUS " Cuda")
MESSAGE(STATUS " Ceres")
MESSAGE(STATUS " Magma dense")
MESSAGE(STATUS "")
MESSAGE(STATUS "OPTIONAL LIBRARIES")
MESSAGE(STATUS " SuiteSparse")
MESSAGE(STATUS " Matlab mex")
MESSAGE(STATUS " OpenMVG (AliceVision)")
MESSAGE(STATUS "---------------------------------------------------------")

# -----------------------------------------------------------------------------
# SYSTEM LIBRARY - FOR LINUX
# -----------------------------------------------------------------------------
IF(UNIX)
	FIND_PACKAGE (Threads)
	SET(LIBS_FILES ${LIBS_FILES} ${CMAKE_THREAD_LIBS_INIT})
	
	FIND_PACKAGE(OpenMP)
	IF (OPENMP_FOUND)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	ENDIF(OPENMP_FOUND)
ENDIF(UNIX)

# -----------------------------------------------------------------------------
# PACKAGES FOR LINUX
# -----------------------------------------------------------------------------
IF(UNIX)
    MESSAGE(STATUS "You can download most of the repositories by calling: ")
    MESSAGE(STATUS " sudo apt-get install cmake")
    MESSAGE(STATUS " sudo apt-get install build-essential")
    MESSAGE(STATUS " sudo apt-get install libboost-all-dev")
    MESSAGE(STATUS " sudo apt-get install libeigen3-dev")
    MESSAGE(STATUS " sudo apt-get install libsuitesparse-dev")
    MESSAGE(STATUS " sudo apt-get install libfreeimage-dev")
    MESSAGE(STATUS " sudo apt-get install libgoogle-glog-dev")  # may cause problems in MEX ( unresolved symbol _Unwind.. )
    MESSAGE(STATUS " sudo apt-get install libgflags-dev")
    MESSAGE(STATUS " sudo apt-get install libglew-dev")
    MESSAGE(STATUS " sudo apt-get install qtbase5-dev")
    MESSAGE(STATUS " sudo apt-get install libqt5opengl5-dev")
    MESSAGE(STATUS " sudo apt-get install libatlas-base-dev")
    MESSAGE(STATUS " sudo apt-get install libsuitesparse-dev")
    MESSAGE(STATUS " git clone https://ceres-solver.googlesource.com/ceres-solver")
    MESSAGE(STATUS " cd ceres-solver")
    MESSAGE(STATUS " mkdir build")
    MESSAGE(STATUS " cd build")
    MESSAGE(STATUS " cmake ..")
    MESSAGE(STATUS " sudo make install")
    MESSAGE(STATUS "---------------------------------------------------------")
ENDIF(UNIX)


# -----------------------------------------------------------------------------
# PREBUILDED LIBRARIES FOR WINDOWS - SETUP PATHS
# -----------------------------------------------------------------------------
IF(WIN32)
    IF (NOT EXISTS "${PROJECT_SOURCE_DIR}/libs")
            MESSAGE(STATUS "You can download the prebuild libraries for win64 at: 'cmp.felk.cvut.cz/~policmic/files/uncer_libs_win64.zip'")
            MESSAGE(STATUS "Please unpack it into the root project directory.")
            MESSAGE(STATUS "---------------------------------------------------------")
    ELSE (NOT EXIST "${PROJECT_SOURCE_DIR}/libs")
            SET(USE_PREBUILD_WIN64_LIBS "ON")
            MESSAGE(STATUS "There were foud precompiled libraries. They will be used in case of not detected the instaled ones.")
            MESSAGE(STATUS "---------------------------------------------------------")

            # Third party paths
            SET(INTERNAL_LIBRARIES "${PROJECT_SOURCE_DIR}/libs/compiled_libraries/windows")
            SET(CERES_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/ceres-solver")
            SET(CUDA8_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/cuda-8.0")
            SET(EIGEN3_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/eigen-3.3.4")
            SET(GFLAGS_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/gflags/BUILD_")
            SET(MKL_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/IntelMKL")
            SET(MAGMA_ROOT_DIR_PRE "${PROJECT_SOURCE_DIR}/libs/third_party/magma-2.2.0")
            SET(SUITESPARSE_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/suitesparse-metis-for-windows")
    ENDIF(NOT EXISTS "${PROJECT_SOURCE_DIR}/libs")
ENDIF(WIN32)


# -----------------------------------------------------------------------------
# Gflags
# -----------------------------------------------------------------------------
find_package ( gflags )

IF(DEFINED gflags_VERSION_STRING OR DEFINED GFLAGS_VERSION_STRING)
        # different cmake variable with lower or upper case exist
        SET(GFLAGS_INCLUDE_DIR ${GFLAGS_INCLUDE_DIR} ${gflags_INCLUDE_DIR})
        SET(GFLAGS_LIBRARIES ${GFLAGS_LIBRARIES} ${gflags_LIBRARIES})
	MESSAGE( "Found Gflags ${GFLAGS_VERSION_STRING}" )
	MESSAGE( STATUS "> GFLAGS_INCLUDE_DIR:   ${GFLAGS_INCLUDE_DIR}" )
	MESSAGE( STATUS "> GFLAGS_LIBRARIES:   ${GFLAGS_LIBRARIES}" )

ELSE()
	IF(UNIX)
		MESSAGE(WARNING "Gflags was not found...")
	ELSE(UNIX)
		# Precompiled win64
		IF(${USE_PREBUILD_WIN64_LIBS})
			MESSAGE( "Found Gflags in precompiled libraries." )
			SET(GFLAGS_INCLUDE_DIR "${GFLAGS_ROOT_DIR}/include")
			SET(GFLAGS_LIBRARIES 
				optimized "${INTERNAL_LIBRARIES}/Release/gflags_static.lib"			
				debug "${INTERNAL_LIBRARIES}/Debug/gflags_static.lib" )	
		ELSE(${USE_PREBUILD_WIN64_LIBS})
		
			# Set up gflags install directory
			SET(GFLAGS_DIR CACHE PATH FORCE)
			IF(NOT GFLAGS_DIR STREQUAL "")

				MESSAGE( "Gflags should be on manualy selected path ..." )
				MESSAGE(STATUS "GFLAGS_DIR: ${GFLAGS_DIR}")
				IF(EXISTS "${GFLAGS_DIR}/include/gflags/gflags.h" AND EXISTS "${GFLAGS_DIR}/include/gflags/gflags.h")
				    SET(GFLAGS_INCLUDE_DIR "${GFLAGS_DIR}/include") # "C:/Program Files/gflags/include"
				    SET(GFLAGS_LIBRARIES 
						optimized "${INTERNAL_LIBRARIES}/Release/gflags_static.lib"			# install generate only Release or Debug, not both libs
						debug "${INTERNAL_LIBRARIES}/Debug/gflags_static.lib")				# so, it is generated in advance
				ELSE()
					MESSAGE(WARNING "Setup the directory where is installed gflags to variable GFLAGS_DIR. It is usually \"C:/Program Files/gflags\" and contains \"./include/gflags/gflags.h\" etc.")
				ENDIF()
			ELSE()
				MESSAGE(WARNING "Setup the directory where is installed gflags to variable GFLAGS_DIR. It is usually \"C:/Program Files/gflags\" and contains \"./include/gflags/gflags.h\" etc.")
			ENDIF()
		ENDIF(${USE_PREBUILD_WIN64_LIBS})
	ENDIF(UNIX)
ENDIF()

# SHLWAPI - Used by gflags on WIN32
if(WIN32)
find_library(SHLWAPI Shlwapi.lib) 
SET(GFLAGS_LIBRARIES ${GFLAGS_LIBRARIES} ${SHLWAPI})
else()
set(SHLWAPI)
endif()

# -----------------------------------------------------------------------------
# EIGEN
# -----------------------------------------------------------------------------
FIND_PACKAGE (Eigen3 3.3 NO_MODULE)

IF (EIGEN3_FOUND)
	MESSAGE( "Found Eigen${EIGEN3_VERSION_STRING}" )
	MESSAGE( STATUS "> EIGEN3_INCLUDE_DIRS:   ${EIGEN3_INCLUDE_DIRS}"   )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
  add_definitions("-DEIGEN_DONT_ALIGN_STATICALLY=1")
  add_definitions("-DEIGEN_DONT_VECTORIZE=1")

ELSE(EIGEN3_FOUND)
	IF(UNIX)
		IF(NOT EXISTS "/usr/local/include/eigen3")
			SET(INCLUDE_DIRS ${INCLUDE_DIRS} "/usr/local/include/eigen3")
			MESSAGE( "Eigen was included manualy, used path: /usr/local/include/eigen3" )
		ELSE(NOT EXISTS "/usr/local/include/eigen3")
			MESSAGE(WARNING "Eigen was not found...")
		ENDIF(NOT EXISTS "/usr/local/include/eigen3")	
	ELSE(UNIX)
		IF(${USE_PREBUILD_WIN64_LIBS})
			SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${EIGEN3_ROOT_DIR}")	# use the internal headers
			MESSAGE( "Eigen is included from precompiled libraries, used path: ${EIGEN3_ROOT_DIR}" )
		ELSE(${USE_PREBUILD_WIN64_LIBS})
			SET(EIGEN3_DIR CACHE PATH FORCE)
			IF(NOT EIGEN3_DIR STREQUAL "")
				MESSAGE(STATUS "EIGEN3_DIR: ${EIGEN3_DIR}")
				IF(EXISTS "${EIGEN3_DIR}")
					SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${EIGEN3_DIR}")			
				ELSE()
					MESSAGE(WARNING "Setup the directory where is installed eigen to variable EIGEN3_DIR. It is usually 'C:/Program Files/eigen'")
				ENDIF()
			ELSE()
				MESSAGE(WARNING "Setup the directory where is installed eigen to variable EIGEN3_DIR. It is usually 'C:/Program Files/eigen'")
			ENDIF()
		ENDIF(${USE_PREBUILD_WIN64_LIBS})	
	ENDIF(UNIX)
ENDIF(EIGEN3_FOUND)



# -----------------------------------------------------------------------------
# SuiteSparse 
# -----------------------------------------------------------------------------
FIND_PACKAGE(SuiteSparse QUIET)

IF( SuiteSparse_FOUND)
	MESSAGE( "Found SuiteSparse" )
	MESSAGE( STATUS "> SUITESPARSE_INCLUDE_DIRS:   ${SUITESPARSE_INCLUDE_DIRS}" )
	MESSAGE( STATUS "> SUITESPARSE_LIBRARIES:   ${SUITESPARSE_LIBRARIES}" )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${SUITESPARSE_INCLUDE_DIRS})
	SET(LIBS_FILES ${LIBS_FILES} ${SUITESPARSE_LIBRARIES}) 
	
ELSE (SuiteSparse_FOUND)
	IF(UNIX)
		MESSAGE(WARNING "SuiteSparse was not found...")
	ELSE(UNIX)
		IF(${USE_PREBUILD_WIN64_LIBS})
			MESSAGE("Found SuiteSparse in precompiled libraries.")
			SET(INCLUDE_DIRS ${INCLUDE_DIRS} 
				"${SUITESPARSE_ROOT_DIR}/install/include/suitesparse"        
				"${SUITESPARSE_ROOT_DIR}/install/include"
				"${SUITESPARSE_ROOT_DIR}/SuiteSparse/CXSparse/Include" )      
			SET(LIBS_FILES ${LIBS_FILES} 
  			optimized "${INTERNAL_LIBRARIES}/Release/suitesparseconfig.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libcholmod.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libccolamd.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libcamd.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libcolamd.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libamd.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/liblapack.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libblas.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/metis.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libcxsparse.lib"
  			optimized "${INTERNAL_LIBRARIES}/Release/libspqr.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/suitesparseconfigd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libcholmodd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libccolamdd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libcamdd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libcolamdd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libamdd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/lapack.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/blas.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/metisd.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libcxsparsed.lib"
  			debug "${INTERNAL_LIBRARIES}/Debug/libspqrd.lib" )
		ELSE(${USE_PREBUILD_WIN64_LIBS})
			MESSAGE( WARNING "SuiteSparse was not found... please install it first or download the precompiled libs." )
		ENDIF(${USE_PREBUILD_WIN64_LIBS}) 
	ENDIF(UNIX) 
ENDIF(SuiteSparse_FOUND)


# -----------------------------------------------------------------------------
# CERES
# -----------------------------------------------------------------------------
FIND_PACKAGE( Ceres )
IF ( CERES_FOUND )
	MESSAGE( "Found CERES ${CERES_VERSION}" )
	MESSAGE( STATUS "> CERES_INCLUDE_DIRS:   ${CERES_INCLUDE_DIRS}"   )
	MESSAGE( STATUS "> CERES_LIBRARIES: ${CERES_LIBRARIES}" )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})
	IF(WIN32)
		SET(LIBS_FILES ${LIBS_FILES} optimized ${CERES_LIBRARIES})
		IF(${USE_PREBUILD_WIN64_LIBS})
			SET(LIBS_FILES ${LIBS_FILES} debug "${INTERNAL_LIBRARIES}/Debug/ceres-debug.lib")
			MESSAGE( STATUS "Ceres debud library was found in precompiled libraries and added to the project.")
		
		ELSE(${USE_PREBUILD_WIN64_LIBS})
			SET(CERES_DEBUG_LIBRARY CACHE PATH FORCE)
			IF(NOT CERES_DEBUG_LIBRARY STREQUAL "")
				SET(LIBS_FILES ${LIBS_FILES} debug "${CERES_DEBUG_LIBRARY}")
			ELSE(NOT CERES_DEBUG_LIBRARY STREQUAL "")
				MESSAGE( STATUS "Ceres found but the debud compilation was not found... please add path to ceres-debug.lib to variable CERES_DEBUG_LIBRARY to enable debuging of the code." )
			ENDIF(NOT CERES_DEBUG_LIBRARY STREQUAL "")
		ENDIF(${USE_PREBUILD_WIN64_LIBS})	
	ELSE(WIN32)
		SET(LIBS_FILES ${LIBS_FILES} ${CERES_LIBRARIES})
	ENDIF(WIN32)
	
ELSE( CERES_FOUND )
	IF(UNIX)
		MESSAGE( WARNING "Ceres was not found ... please install Ceres libraries. You can follow the guideline: http://ceres-solver.org/installation.html")
	ELSE(UNIX)
		IF(${USE_PREBUILD_WIN64_LIBS})
			MESSAGE( "Ceres was included from precompiled libraries." )
			MESSAGE( STATUS "> Include paths: ./config, ./include, ./internal, ./internal/ceres and ./ceres/miniglog" )
			MESSAGE( STATUS "> Release library: ${INTERNAL_LIBRARIES}/Release/ceres.lib" )
			MESSAGE( STATUS "> Debug library: ${INTERNAL_LIBRARIES}/Debug/ceres-debug.lib" )
			SET(INCLUDE_DIRS ${INCLUDE_DIRS} 
				"${CERES_ROOT_DIR}/config"        
				"${CERES_ROOT_DIR}/include"
				"${CERES_ROOT_DIR}/internal"
				"${CERES_ROOT_DIR}/internal/ceres"
				"${CERES_ROOT_DIR}/internal/ceres/miniglog") 
			SET(LIBS_FILES ${LIBS_FILES} 
				optimized "${INTERNAL_LIBRARIES}/Release/ceres.lib"
				debug "${INTERNAL_LIBRARIES}/Debug/ceres-debug.lib")
			
		ELSE(${USE_PREBUILD_WIN64_LIBS})
			MESSAGE( WARNING "Ceres was not found ... please install Ceres libraries. You can follow the guideline: http://ceres-solver.org/installation.html")
		ENDIF(${USE_PREBUILD_WIN64_LIBS})
	ENDIF(UNIX)
ENDIF( CERES_FOUND )




# -----------------------------------------------------------------------------
# IntelMKL -> BLAS & LAPACK
# -----------------------------------------------------------------------------
option(UNC_MKL_SUPPORT "Include Intel MKL libraries (include BLAS and LAPACK)" OFF)
IF (UNC_MKL_SUPPORT)
	find_package(MKL QUIET)
	IF(MKL_FOUND)
		MESSAGE( "Found Intel MKL" )
		MESSAGE( STATUS "> MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}"   )
		MESSAGE( STATUS "> MKL_LIBRARIES: ${MKL_LIBRARIES}" )
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${MKL_INCLUDE_DIR})
		SET(LIBS_FILES ${LIBS_FILES} ${MKL_LIBRARIES})
		
	ELSE(MKL_FOUND)
		IF(UNIX)
			MESSAGE( WARNING "Intel MKL was not found ... please install Intel MKL libraries or shutdown UNC_MKL_SUPPORT and add another version of BLAS and LAPACK.")
		ELSE(UNIX)
			IF(${USE_PREBUILD_WIN64_LIBS})
				MESSAGE("Found Intel MKL in precompiled libraries." )
				SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${MKL_ROOT_DIR}/include")
				SET(LIBS_FILES ${LIBS_FILES}
  					"${INTERNAL_LIBRARIES}/Release/mkl_intel_lp64_dll.lib"
  					"${INTERNAL_LIBRARIES}/Release/mkl_intel_thread_dll.lib"
  					"${INTERNAL_LIBRARIES}/Release/mkl_core_dll.lib"
  					"${INTERNAL_LIBRARIES}/Release/libiomp5md.lib")
			
			ELSE(${USE_PREBUILD_WIN64_LIBS})
				SET(MKL_ALL_LIBRARIES CACHE PATH FORCE)
				SET(MKL_INCLUDE CACHE PATH FORCE)
				
				IF(NOT MKL_ALL_LIBRARIES STREQUAL "")
					MESSAGE(STATUS "MKL_ALL_LIBRARIES: ${MKL_ALL_LIBRARIES}")
					SET(LIBS_FILES ${LIBS_FILES} ${MKL_LIBRARIES})
				ELSE(NOT MKL_ALL_LIBRARIES STREQUAL "")
					MESSAGE("Intel MKL was not found ... please add path to required libraries (use Intel MKL link advisor)")
					MESSAGE("  e.g. MKL_ALL_LIBRARIES = <PATH_TO_LIBS>/mkl_intel_lp64_dll.lib,<PATH_TO_LIBS/mkl_intel_thread_dll.lib...")
				ENDIF(NOT MKL_ALL_LIBRARIES STREQUAL "")	
				
				IF(NOT MKL_INCLUDE STREQUAL "")
					MESSAGE(STATUS "MKL_INCLUDE: ${MKL_INCLUDE}")
					SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${MKL_INCLUDE)")
				ELSE(NOT MKL_INCLUDE STREQUAL "")
					MESSAGE( WARNING "Intel MKL was not found ... please add path to include directory (e.g. MKL_INCLUDE = <MKL_ROOT_DIR>/include)")
				ENDIF(NOT MKL_INCLUDE STREQUAL "")	
			ENDIF(${USE_PREBUILD_WIN64_LIBS})
		ENDIF(UNIX)
	ENDIF(MKL_FOUND)
ENDIF(UNC_MKL_SUPPORT)


# -----------------------------------------------------------------------------
# LAPACK AND BLAS LIBS
# -----------------------------------------------------------------------------
IF (NOT UNC_MKL_SUPPORT AND UNIX)
	find_package( LAPACK )
	IF(LAPACK_FOUND)
		SET(LIBS_FILES ${LIBS_FILES} ${LAPACK_LIBRARIES})
		MESSAGE("Found lapack and blas.")
	ELSE(LAPACK_FOUND)
		MESSAGE(STATUS "Lapack not found ... please make sure that you have instaled at least one LAPACK and BLAS library.")
		# TODO check suitesparse etc.
	ENDIF(LAPACK_FOUND)
ENDIF()


# -----------------------------------------------------------------------------
# MATLAB interface - MEX libs ( C:/Program Files/MATLAB/R2017a/extern/include; C:/Program Files/MATLAB/R2017a/extern/lib/win64/microsoft )
# -----------------------------------------------------------------------------

OPTION(UNC_MATLAB_SUPPORT "Include MATLAB libraries (enable input from MATLAB)" OFF)
IF (UNC_MATLAB_SUPPORT)

	FIND_PACKAGE(Matlab)
	IF(Matlab_FOUND)
		ADD_DEFINITIONS(-DUSE_MATLAB)
		MESSAGE("Found Matlab libraries")
		MESSAGE(STATUS "Matlab_INCLUDE_DIRS: ${Matlab_INCLUDE_DIRS}")
		MESSAGE(STATUS "Matlab_LIBRARIES: ${Matlab_LIBRARIES}")
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${Matlab_INCLUDE_DIRS})
		GET_FILENAME_COMPONENT(MATLAB_LIBS_ROOT ${Matlab_LIBRARIES} DIRECTORY)  # ${Matlab_LIBRARIES} doesn't contain all required libs
		IF(UNIX)
			SET(LIBS_FILES ${LIBS_FILES} 
				"${MATLAB_LIBS_ROOT}/libmx.so"
				"${MATLAB_LIBS_ROOT}/libmat.so"
				"${MATLAB_LIBS_ROOT}/libeng.so"
				"${MATLAB_LIBS_ROOT}/libmex.so")
		ELSE(UNIX)
			SET(LIBS_FILES ${LIBS_FILES} 
				"${MATLAB_LIBS_ROOT}/libmx.lib"
				"${MATLAB_LIBS_ROOT}/libmat.lib"
				"${MATLAB_LIBS_ROOT}/libeng.lib"
				"${MATLAB_LIBS_ROOT}/libmex.lib")
		ENDIF(UNIX)
	ELSE(Matlab_FOUND)
		SET(MEX_INCLUDE_DIR CACHE PATH FORCE)
		IF(NOT MEX_INCLUDE_DIR STREQUAL "")
  			MESSAGE(STATUS "MEX_INCLUDE_DIR: ${MEX_INCLUDE_DIR}")
  			SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${MEX_INCLUDE_DIR})     
		ELSE(NOT MEX_INCLUDE_DIR STREQUAL "")
  			MESSAGE(WARNING "Setup the include directory for MATLAB MEX interface (contatins mex.h).")
		ENDIF(NOT MEX_INCLUDE_DIR STREQUAL "")

		SET(MEX_LIBRARIES_DIR CACHE PATH FORCE)
		IF(NOT MEX_LIBRARIES_DIR STREQUAL "")
  			MESSAGE(STATUS "MEX_LIBRARIES_DIR: ${MEX_LIBRARIES_DIR}")
  			IF(WIN32)
  				SET(LIBS_FILES ${LIBS_FILES} 
  				"${MEX_LIBRARIES_DIR}/libmx.lib"
  				"${MEX_LIBRARIES_DIR}/libmat.lib"
  				"${MEX_LIBRARIES_DIR}/libeng.lib"
  				"${MEX_LIBRARIES_DIR}/libmex.lib")     
  			ELSE(WIN32)
  				SET(LIBS_FILES ${LIBS_FILES} 
  				"${MEX_LIBRARIES_DIR}/libmx.so"
  				"${MEX_LIBRARIES_DIR}/libmat.so"
  				"${MEX_LIBRARIES_DIR}/libeng.so"
  				"${MEX_LIBRARIES_DIR}/libmex.so")
  			ENDIF(WIN32)
		ELSE(NOT MEX_LIBRARIES_DIR STREQUAL "")
			MESSAGE(WARNING "Please setup the library directory for MATLAB MEX interface (contatins libmx.lib, libmat.lib, libeng.lib, libmex.lib).")
		ENDIF(NOT MEX_LIBRARIES_DIR STREQUAL "") 
	ENDIF(Matlab_FOUND)	
ENDIF(UNC_MATLAB_SUPPORT)


# -----------------------------------------------------------------------------
# MAGMA
# -----------------------------------------------------------------------------
FIND_PACKAGE(MAGMA QUIET)

IF(MAGMA_FOUND)
	MESSAGE("Found Magma")
	MESSAGE(STATUS "MAGMA_INCLUDE_DIRS: ${MAGMA_INCLUDE_DIRS}")
	MESSAGE(STATUS "MAGMA_LIBRARIES: ${MAGMA_LIBRARIES}")
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${MAGMA_INCLUDE_DIRS})
	SET(LIBS_FILES ${LIBS_FILES} ${MAGMA_LIBRARIES})
	
ELSE(MAGMA_FOUND)
	IF(UNIX)
		SET(MAGMA_ROOT_DIR CACHE PATH FORCE)
		IF(NOT MAGMA_ROOT_DIR STREQUAL "")
			MESSAGE(STATUS "MAGMA_ROOT_DIR: ${MAGMA_ROOT_DIR}")
			SET(INCLUDE_DIRS ${INCLUDE_DIRS}
    		"${MAGMA_ROOT_DIR}/include"
    		"${MAGMA_ROOT_DIR}/control"
    		"${MAGMA_ROOT_DIR}/sparse/include"
			"${MAGMA_ROOT_DIR}/sparse/control"
    		"${MAGMA_ROOT_DIR}/testing")   
		ELSE(NOT MAGMA_ROOT_DIR STREQUAL "")
			MESSAGE(WARNING "Magma not found ... please setup the MAGMA_ROOT_DIR")
		ENDIF(NOT MAGMA_ROOT_DIR STREQUAL "")  
 
		SET(MAGMA_LIBRARIES_DIR CACHE PATH FORCE)
		IF(NOT MAGMA_LIBRARIES_DIR STREQUAL "")
			MESSAGE(STATUS "MAGMA_LIBRARIES_DIR: ${MAGMA_LIBRARIES_DIR}")
			SET(LIBS_FILES ${LIBS_FILES}
    		"${MAGMA_LIBRARIES_DIR}/libmagma.a"
		  # "${MAGMA_LIBRARIES_DIR}/libmagma_sparse.a"
		    "${MAGMA_LIBRARIES_DIR}/libtester.a"
    		"${MAGMA_LIBRARIES_DIR}/liblapacktest.a") 
		ELSE(NOT MAGMA_LIBRARIES_DIR STREQUAL "")
			MESSAGE(WARNING "Magma not found ... please setup the MAGMA_LIBRARIES_DIR")
		ENDIF(NOT MAGMA_LIBRARIES_DIR STREQUAL "")  
		
	ELSE(UNIX)
		IF(${USE_PREBUILD_WIN64_LIBS})
			MESSAGE("Found Magma in precompiled libraries.")
			SET(INCLUDE_DIRS ${INCLUDE_DIRS}
  				"${MAGMA_ROOT_DIR_PRE}/include"
  				"${MAGMA_ROOT_DIR_PRE}/control"
  				"${MAGMA_ROOT_DIR_PRE}/sparse/include"
  				"${MAGMA_ROOT_DIR_PRE}/testing")
			SET(LIBS_FILES ${LIBS_FILES}
  				optimized "${INTERNAL_LIBRARIES}/Release/tester.lib"
  				optimized "${INTERNAL_LIBRARIES}/Release/lapacktest.lib"
  				optimized "${INTERNAL_LIBRARIES}/Release/magma.lib"
  				optimized "${INTERNAL_LIBRARIES}/Release/magma_sparse.lib"
  				debug "${INTERNAL_LIBRARIES}/Debug/tester.lib"
  				debug "${INTERNAL_LIBRARIES}/Debug/lapacktest.lib"
  				debug "${INTERNAL_LIBRARIES}/Debug/magma.lib"
  				debug "${INTERNAL_LIBRARIES}/Debug/magma_sparse.lib")
			
		ELSE(${USE_PREBUILD_WIN64_LIBS})
			
			SET(MAGMA_ROOT_DIR CACHE PATH FORCE)
			IF(NOT MAGMA_ROOT_DIR STREQUAL "")
				MESSAGE(STATUS "MAGMA_ROOT_DIR: ${MAGMA_ROOT_DIR}")
				SET(INCLUDE_DIRS ${INCLUDE_DIRS}
    			"${MAGMA_ROOT_DIR}/include"
    			"${MAGMA_ROOT_DIR}/control"
    			"${MAGMA_ROOT_DIR}/sparse/include"
				"${MAGMA_ROOT_DIR}/sparse/control"
    			"${MAGMA_ROOT_DIR}/testing")   
			ELSE(NOT MAGMA_ROOT_DIR STREQUAL "")
				MESSAGE(WARNING "Magma not found ... please setup the MAGMA_ROOT_DIR")
			ENDIF(NOT MAGMA_ROOT_DIR STREQUAL "")  
 
			SET(MAGMA_LIBRARIES_DIR CACHE PATH FORCE)
			IF(NOT MAGMA_LIBRARIES_DIR STREQUAL "")
				MESSAGE(STATUS "MAGMA_LIBRARIES_DIR: ${MAGMA_LIBRARIES_DIR}")
				SET(LIBS_FILES ${LIBS_FILES}
    			"${MAGMA_LIBRARIES_DIR}/libmagma.lib"
			    "${MAGMA_LIBRARIES_DIR}/libmagma_sparse.lib"
				"${MAGMA_LIBRARIES_DIR}/libtester.lib"
    			"${MAGMA_LIBRARIES_DIR}/liblapacktest.lib") 
			ELSE(NOT MAGMA_LIBRARIES_DIR STREQUAL "")
				MESSAGE(WARNING "Magma not found ... please setup the MAGMA_LIBRARIES_DIR")
			ENDIF(NOT MAGMA_LIBRARIES_DIR STREQUAL "")  
					
		ENDIF(${USE_PREBUILD_WIN64_LIBS})
	ENDIF(UNIX)
ENDIF(MAGMA_FOUND)


# -----------------------------------------------------------------------------
# CUDA
# -----------------------------------------------------------------------------
FIND_PACKAGE( CUDA )
IF ( CUDA_FOUND )
	MESSAGE( "Found CUDA ${CUDA_VERSION}" )
	MESSAGE( STATUS "> CUDA_INCLUDE_DIRS:   ${CUDA_INCLUDE_DIRS}"   )
	MESSAGE( STATUS "> CUDA_LIBRARIES: ${CUDA_LIBRARIES}" )
	MESSAGE( STATUS "> CUDA_CUBLAS_LIBRARIES: ${CUDA_CUBLAS_LIBRARIES}" )
	MESSAGE( STATUS "> CUDA_CUDART_LIBRARY: ${CUDA_cusparse_LIBRARY}" )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${CUDA_INCLUDE_DIRS})
	SET(LIBS_FILES ${LIBS_FILES} ${CUDA_LIBRARIES};${CUDA_CUBLAS_LIBRARIES};${CUDA_cusparse_LIBRARY})

ELSE( CUDA_FOUND )
	IF(WIN32)
		MESSAGE("Found CUDA 8.0 in precompiled libraries, it will be included manualy which may cause some trouble.")
		MESSAGE("Please install 'Cuda Develeopment Toolkit' on you computer.")
		MESSAGE( STATUS "> include path: ${CUDA8_ROOT_DIR}/include" )
		MESSAGE( STATUS "> libraries directory: ${INTERNAL_LIBRARIES}/Release" )
		MESSAGE( STATUS "> added libraries: cublas.lib, cusparse.lib, cublas_device.lib, cudadevrt.lib, cudart_static.lib" )
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${CUDA8_ROOT_DIR}/include") 
		SET(LIBS_FILES ${LIBS_FILES}
			"${INTERNAL_LIBRARIES}/Release/cublas.lib"
			"${INTERNAL_LIBRARIES}/Release/cusparse.lib"
			"${INTERNAL_LIBRARIES}/Release/cublas_device.lib"
			"${INTERNAL_LIBRARIES}/Release/cudadevrt.lib"
			"${INTERNAL_LIBRARIES}/Release/cudart_static.lib" )
	ELSE(WIN32)
		MESSAGE( WARNING "Install Cuda libraries first. Follow this guideline: https://askubuntu.com/questions/799184/how-can-i-install-cuda-on-ubuntu-16-04")
	ENDIF(WIN32)
ENDIF( CUDA_FOUND )



# -----------------------------------------------------------------------------
# OpenMVG
# -----------------------------------------------------------------------------
# TODO add find package ...
option(UNC_OPENMVG_SUPPORT "Include OpenMVG libraries (enable input from OpenMVG)" OFF)
IF (UNC_OPENMVG_SUPPORT)
	add_definitions(-DUSE_OPENMVG)
	
	SET(OPENMVG_SRC_DIR CACHE PATH FORCE)
	IF(NOT OPENMVG_SRC_DIR STREQUAL "")
		MESSAGE(STATUS "OPENMVG_SRC_DIR: ${OPENMVG_SRC_DIR}")
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} 
			"${OPENMVG_SRC_DIR}"
			"${OPENMVG_SRC_DIR}/third_party/jpeg"
			"${OPENMVG_SRC_DIR}/third_party/png"
			"${OPENMVG_SRC_DIR}/third_party/zlib"
			"${OPENMVG_SRC_DIR}/third_party/tiff"
			"${OPENMVG_SRC_DIR}/third_party"
			"${OPENMVG_SRC_DIR}/dependencies"
			"${OPENMVG_SRC_DIR}/dependencies/cereal/include"
			"${OPENMVG_SRC_DIR}/third_party/lemon"
			"${OPENMVG_SRC_DIR}/third_party/flann/src/cpp"
			"${OPENMVG_SRC_DIR}/dependencies/osi_clp/CoinUtils/src"
			"${OPENMVG_SRC_DIR}/dependencies/osi_clp/Clp/src"
			"${OPENMVG_SRC_DIR}/dependencies/osi_clp/Clp/src/OsiClp"
			"${OPENMVG_SRC_DIR}/dependencies/osi_clp/Osi/src/Osi")     
	ELSE()
		MESSAGE(WARNING "Setup the OpenMVG root directory (contatins ./dependencies, ./third_party).")
	ENDIF()
	
	SET(OPENMVG_BUILD_DIR CACHE PATH FORCE)
	IF(NOT OPENMVG_BUILD_DIR STREQUAL "")
		MESSAGE(STATUS "OPENMVG_BUILD_DIR: ${OPENMVG_BUILD_DIR}")
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} 
			"${OPENMVG_BUILD_DIR}/generated"
			"${OPENMVG_BUILD_DIR}/third_party/jpeg/config"
			"${OPENMVG_BUILD_DIR}/third_party/png/config"
			"${OPENMVG_BUILD_DIR}/third_party/tiff"
			"${OPENMVG_BUILD_DIR}/third_party/lemon") 

		IF(WIN32)  
  			SET(LIBS_FILES ${LIBS_FILES} 
  				optimized "${OPENMVG_BUILD_DIR}/third_party/flann/lib/Release/flann_cpp_s.lib"
  				debug "${OPENMVG_BUILD_DIR}/third_party/flann/lib/Debug/flann_cpp_s.lib")   
		ENDIF()
	ELSE()
		MESSAGE(WARNING "Setup the OpenMVG root directory (contatins ./src).")
	ENDIF()
	
	SET(OPENMVG_LIBRARIES_DIR CACHE PATH FORCE)
	IF(NOT OPENMVG_LIBRARIES_DIR STREQUAL "")
		MESSAGE(STATUS "OPENMVG_LIBRARIES_DIR: ${OPENMVG_LIBRARIES_DIR}")
    IF(WIN32)
      SET(LIBS_FILES ${LIBS_FILES} 
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_system.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_features.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_sfm.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/stlplus.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_lInftyComputerVision.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/lib_clp.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/lib_OsiClpSolver.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/lib_CoinUtils.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/lib_Osi.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_matching.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_image.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/png.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/jpeg.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/tiff.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/zlib.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_multiview.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_numeric.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/lemon.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/vlsift.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/openMVG_exif.lib"
  			optimized "${OPENMVG_LIBRARIES_DIR}/Release/easyexif.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_system.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_features.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_sfm.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/stlplus.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_lInftyComputerVision.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/lib_clp.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/lib_OsiClpSolver.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/lib_CoinUtils.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/lib_Osi.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_matching.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_image.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/png.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/jpeg.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/tiff.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/zlib.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_multiview.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_numeric.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/lemon.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/vlsift.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/openMVG_exif.lib"
  			debug "${OPENMVG_LIBRARIES_DIR}/Debug/easyexif.lib")   
      ELSE()
        SET(LIBS_FILES ${LIBS_FILES} 
        "${OPENMVG_LIBRARIES_DIR}/libopenMVG_system.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_features.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_sfm.a"
  			"${OPENMVG_LIBRARIES_DIR}/libstlplus.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_lInftyComputerVision.a"
  			"${OPENMVG_LIBRARIES_DIR}/liblib_clp.a"
  			"${OPENMVG_LIBRARIES_DIR}/liblib_OsiClpSolver.a"
  			"${OPENMVG_LIBRARIES_DIR}/liblib_CoinUtils.a"
  			"${OPENMVG_LIBRARIES_DIR}/liblib_Osi.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_matching.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_image.a"
  			#"${OPENMVG_LIBRARIES_DIR}/libpng.a"
  			#"${OPENMVG_LIBRARIES_DIR}/libjpeg.a"
  			#"${OPENMVG_LIBRARIES_DIR}/libtiff.a"
  			#"${OPENMVG_LIBRARIES_DIR}/libzlib.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_multiview.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_numeric.a"
  			"${OPENMVG_LIBRARIES_DIR}/liblemon.a"
  			"${OPENMVG_LIBRARIES_DIR}/libvlsift.a"
  			"${OPENMVG_LIBRARIES_DIR}/libopenMVG_exif.a"
  			"${OPENMVG_LIBRARIES_DIR}/libeasyexif.a")
      ENDIF()   
	ELSE()
		MESSAGE(WARNING "Setup the library directories for OpenMVG (contians the subdirecries ./Release, ./Debug.")
	ENDIF()
ENDIF()

# -----------------------------------------------------------------------------
# Add project source code
# -----------------------------------------------------------------------------
add_subdirectory(src)


# -----------------------------------------------------------------------------
# Compile MEX
# -----------------------------------------------------------------------------
OPTION(UNC_BUILD_MEX "Build uncertainty_mex file" OFF)

# Check if the matlab libs are setup
IF(UNC_BUILD_MEX AND NOT UNC_MATLAB_SUPPORT)
	MESSAGE( FATAL_ERROR "Please, include matlab libraries and headers (i.e. UNC_MATLAB_SUPPORT = ON) to build mex file (i.e. UNC_BUILD_MEX = ON)." )
ENDIF(UNC_BUILD_MEX AND NOT UNC_MATLAB_SUPPORT)

IF(UNC_BUILD_MEX) 
        IF(WIN32)
            add_library(uncertainty_mex SHARED ${HEADERS} ${SOURCES} ${MODUL_DEFINITION_FILES})
            target_include_directories(uncertainty_mex PUBLIC ${INCLUDE_DIRS})
            target_link_libraries (uncertainty_mex ${LIBS_FILES})
            add_definitions(-DMATLAB_MEX_FILE)
            SET_TARGET_PROPERTIES(uncertainty_mex PROPERTIES LINK_FLAGS "/export:mexFunction /dll")
            SET_TARGET_PROPERTIES(uncertainty_mex PROPERTIES SUFFIX ".mexw64")
			
	ELSE(WIN32)
            # Build shared library, create another project, generate make, execute it
            ADD_LIBRARY(unc STATIC ${HEADERS} ${SOURCES})
            target_include_directories(unc PUBLIC ${INCLUDE_DIRS})
            target_link_libraries (unc ${LIBS_FILES})
            set_property(TARGET unc PROPERTY POSITION_INDEPENDENT_CODE ON)

            # Generate shell script for compilation
            SET(CM_DIR "${PROJECT_SOURCE_DIR}/mex")
            SET(MEX_SH "${CM_DIR}/mex_link.sh")
            file(WRITE ${MEX_SH}
                "g++ -g -Wall -fPIC -std=c++11 -pthread -fPIC -DMX_COMPAT_32 -DMATLAB_MEX_FILE -DUSE_MATLAB "
                "-I${Matlab_INCLUDE_DIRS} -I../include -c ../mex/uncertainty_mex.cpp -o ../mex/uncertainty_mex.o \n")
            file(APPEND ${MEX_SH} 
                "mex -g -cxx CC='/usr/bin/g++' CXX='/usr/bin/g++' LD='/usr/bin/g++' GCC='/usr/bin/g++' "
                "-DUSE_MATLAB -L./ -L/media/policmic/DATA/MichalPolic/libs/magma-2.2.0/build/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib -L/usr/local/lib "
                "-lpthread -lunc -lmagma -llapacktest -ltester -lcudart_static -lrt -lcublas -lcusparse "
                "-llapack -lf77blas -latlas -lsuitesparseconfig -lspqr -lcholmod -lccolamd -lcamd -lcolamd "
                "-lcolamd -lamd -lgflags -lceres -ldl -lrt -largeArrayDims "
                "LDFLAGS='$LDFLAGS -fopenmp' -output ../mex/uncertainty_mex ../mex/uncertainty_mex.o\n") 
            file(APPEND ${MEX_SH} 
                "mv ../mex/uncertainty_mex.mexa64 ${CMAKE_BINARY_DIR}")   

            # Execute the script after the build
            ADD_CUSTOM_COMMAND(TARGET unc POST_BUILD COMMAND "sh" ${MEX_SH})
        ENDIF(WIN32)
ENDIF(UNC_BUILD_MEX) 




# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# PLATFORM SPECIFIC ERRORS
# -----------------------------------------------------------------------------

# The debud on windows may need additional libraries
SET(INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY "c:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017.4.210/windows/compiler/lib/intel64_win")
IF (WIN32)
	MESSAGE("For the debuging on windows platform may required an additional libraries (e.g. ifconsole.lib).")
	MESSAGE("Please add the library path to your project (e.g. 'c:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017.4.210/windows/compiler/lib/intel64_win')")
	
	#IF (EXISTS ${INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY})
	#	#get_directory_property(OUT_VAR LINK_DIRECTORIES)	
	#	MESSAGE("There was added path to system libs for debuging: ${INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY}")
	#	SET(LINK_DIRECTORIES ${LINK_DIRECTORIES} ${INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY})
	#	#link_directories(${INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY})
	#	
	#ELSE(EXISTS ${INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY})
	#	
	#ENDIF(EXISTS ${INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY})
ENDIF(WIN32)


# -----------------------------------------------------------------------------
# LINUX COMPILER FLAGS
# -----------------------------------------------------------------------------
IF( UNIX )
	SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 ")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
	SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/documents/libs/ceres-solver-install/lib/cmake/Ceres")
ENDIF()


# -----------------------------------------------------------------------------
# COPY OF THE EXAMPLES
# -----------------------------------------------------------------------------
SET(TEST_JACOBIANS_DIR "${PROJECT_SOURCE_DIR}/test/input_jacobians")
file(GLOB TEST_JACOBIANS_FILES "${TEST_JACOBIANS_DIR}/*")
IF(WIN32)
	SET(DLL_LIBRARIES_DIR "${INTERNAL_LIBRARIES}/Release/DLL")
	file(GLOB DLL_FILES "${DLL_LIBRARIES_DIR}/*")
ENDIF()

# Release
SET(RELEASE_DIR "${CMAKE_BINARY_DIR}/Release")
SET(RELEASE_IN_DIR "${CMAKE_BINARY_DIR}/Release/in")
IF(NOT EXISTS ${RELEASE_IN_DIR})
	file(MAKE_DIRECTORY ${RELEASE_IN_DIR})
	file(COPY ${TEST_JACOBIANS_FILES} DESTINATION ${RELEASE_IN_DIR})
	IF(WIN32)
		file(COPY ${DLL_FILES} DESTINATION ${RELEASE_DIR})
	ELSE()
		file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Release/out")
	ENDIF()
ENDIF()

# Debug
SET(DEBUG_DIR "${CMAKE_BINARY_DIR}/Debug")
SET(DEBUG_IN_DIR "${CMAKE_BINARY_DIR}/Debug/in")
IF(NOT EXISTS ${DEBUG_IN_DIR})
	file(MAKE_DIRECTORY ${DEBUG_IN_DIR})
	file(COPY ${TEST_JACOBIANS_FILES} DESTINATION ${DEBUG_IN_DIR})
	IF(WIN32)
		file(COPY ${DLL_FILES} DESTINATION ${DEBUG_DIR})
	ELSE()
		file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/Debug/out")
	ENDIF()
ENDIF()


# -----------------------------------------------------------------------------
# INSTALL LIBRARIES
# -----------------------------------------------------------------------------

# TODO ...



# -----------------------------------------------------------------------------
# STATUS
# -----------------------------------------------------------------------------
MESSAGE(STATUS " ")
MESSAGE(STATUS "---------------------------------------------------------")
MESSAGE(STATUS "SUMMARY")
MESSAGE(STATUS "---------------------------------------------------------")
MESSAGE(STATUS "UNC_BUILD_STATIC_LIB: ${UNC_BUILD_STATIC_LIB}")
MESSAGE(STATUS "UNC_BUILD_SHARED_LIB: ${UNC_BUILD_SHARED_LIB}")

MESSAGE(STATUS "UNC_OPENMVG_SUPPORT: ${UNC_OPENMVG_SUPPORT}")
MESSAGE(STATUS "UNC_MKL_SUPPORT: ${UNC_MKL_SUPPORT}")
MESSAGE(STATUS "UNC_MATLAB_SUPPORT: ${UNC_MATLAB_SUPPORT}")
MESSAGE(STATUS "UNC_BUILD_MEX: ${UNC_BUILD_MEX}")

MESSAGE(STATUS " ")
MESSAGE(STATUS "INCLUDE_DIRS: ${INCLUDE_DIRS}")
MESSAGE(STATUS "LIBS_FILES: ${LIBS_FILES}")
MESSAGE(STATUS " ")
MESSAGE(STATUS " ")


