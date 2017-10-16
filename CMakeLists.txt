CMAKE_MINIMUM_REQUIRED (VERSION 3.5)

# About the project
PROJECT (Uncertatinty)
SET (VERSION_MAJOR 1)
SET (VERSION_MINOR 1)
SET (VERSION_PATCH 0)
SET (VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_CONFIGURATION_TYPES "Debug;Release")
MESSAGE(STATUS "CMAKE_CONFIGURATION_TYPES: ${CMAKE_CONFIGURATION_TYPES}")
SET(CMAKE_BUILD_TYPE "Release")
MESSAGE(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# guard against in-source builds
if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed.")
endif()


# -----------------------------------------------------------------------------
# Linux compilation require c++11 flag
IF( NOT WIN32 )
	SET( GCC_COVERAGE_COMPILE_FLAGS "-lpthread")
	SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 ${GCC_COVERAGE_LINK_FLAGS}" )
	SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )
	SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "~/documents/libs/ceres-solver-install/lib/cmake/Ceres")
ENDIF()

# -----------------------------------------------------------------------------
# Set sources and create executable target
file(GLOB HEADERS "${PROJECT_SOURCE_DIR}/include/*.h")
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
file(GLOB MODUL_DEFINITION_FILES "${PROJECT_SOURCE_DIR}/src/*.def")

message("HEADERS: ${HEADERS}")


# Executable - binary file with txt input
add_executable(uncertainty ${HEADERS} ${SOURCES})

# Create also the static library
#ADD_LIBRARY(uncertainty_lib STATIC ${HEADERS} ${SOURCES})

# Compile a MEX lib
#add_library(uncertainty_mex SHARED ${HEADERS} ${SOURCES} ${MODUL_DEFINITION_FILES})
#add_definitions(-DMATLAB_MEX_FILE)
#SET_TARGET_PROPERTIES(uncertainty_mex PROPERTIES LINK_FLAGS "/export:mexFunction /dll")
#SET_TARGET_PROPERTIES(uncertainty_mex PROPERTIES SUFFIX ".mexw64")

# Include headers for both (executable and library)
SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${PROJECT_SOURCE_DIR}/include")

# -----------------------------------------------------------------------------
# Internal prebuilded libraries path
IF(WIN32)	# Windows
	SET(INTERNAL_LIBRARIES "${PROJECT_SOURCE_DIR}/libs/compiled_libraries/windows")
ELSE()		# Linux
	SET(INTERNAL_LIBRARIES "${PROJECT_SOURCE_DIR}/libs/compiled_libraries/linux")
ENDIF()

# Third party paths
SET(CERES_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/ceres-solver")
SET(CUDA8_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/cuda-8.0")
SET(EIGEN3_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/eigen-3.3.4")
SET(GFLAGS_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/gflags")
SET(MKL_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/IntelMKL")
IF(WIN32)
  SET(MAGMA_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/magma-2.2.0")
ENDIF()
SET(SUITESPARSE_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/suitesparse-metis-for-windows")
# -----------------------------------------------------------------------------

# SYSTEM LIBRARY - FOR LINUX
IF( NOT WIN32 )
	find_package (Threads)
	target_link_libraries (uncertainty ${CMAKE_THREAD_LIBS_INIT})
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif()
ENDIF()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# LINUX - INCLUDE ALL LIBRARIES IN DIRECTORY
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# TODO: REPLACE THIS WITH APROPRIATE LINKING OF LIBRARIES
#IF(NOT WIN32)
#	FILE(GLOB LINUX_LIB_SEARCH "${INTERNAL_LIBRARIES}/*")
#	SET(LIBS_FILES ${LIBS_FILES} ${LINUX_LIB_SEARCH} "-lpthread -lm")
#ENDIF()


# -----------------------------------------------------------------------------
# EIGEN
# -----------------------------------------------------------------------------
find_package (Eigen3 NO_MODULE)

IF (EIGEN3_FOUND)
	message( STATUS "Found Eigen${EIGEN3_VERSION_STRING}" )
	message( STATUS "> EIGEN3_INCLUDE_DIRS:   ${EIGEN3_INCLUDE_DIRS}"   )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})

	#SET(LIBS_FILES ${LIBS_FILES} "Eigen3::Eigen")
	
ELSE(EIGEN3_FOUND)
	message( STATUS "Eigen wasn't found.../n" )
	IF (UNIX)
		IF(NOT EXISTS "/usr/local/include/eigen3")
			MESSAGE( WARNING "First, you should install Eigen! You can use Cmake at directory 'libs/third_party/eigen3.3.4'. /n
							 It will add necessary paths and disable functions which are not compatible with system libraries.")
		ELSE()
			message( STATUS "Eigen was included manualy, used path: /usr/local/include/eigen3" )
			SET(INCLUDE_DIRS ${INCLUDE_DIRS} "/usr/local/include/eigen3")			# internal headers do not work, it must be instaled
		ENDIF()
	ELSEIF(WIN32)
		message( STATUS "Eigen will be included manualy, used path: ${EIGEN3_ROOT_DIR}" )
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${EIGEN3_ROOT_DIR}")	# use the internal headers
	ENDIF()
ENDIF(EIGEN3_FOUND)


# -----------------------------------------------------------------------------
# Gflags
# -----------------------------------------------------------------------------
find_package ( gflags )

IF(NOT GFLAGS_VERSION_STRING STREQUAL "")
	message( STATUS "Found Gflags ${GFLAGS_VERSION_STRING}" )
	message( STATUS "> GFLAGS_INCLUDE_DIR:   ${GFLAGS_INCLUDE_DIR}" )
	message( STATUS "> GFLAGS_LIBRARIES:   ${GFLAGS_LIBRARIES}" )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${GFLAGS_INCLUDE_DIR})
	SET(LIBS_FILES ${LIBS_FILES} ${GFLAGS_LIBRARIES}) 
	
ELSE()
	SET(GFLAGS_DIR CACHE PATH FORCE)
	IF(NOT GFLAGS_DIR STREQUAL "")
		MESSAGE(STATUS "GFLAGS_DIR: ${GFLAGS_DIR}")
		IF(EXISTS "${GFLAGS_DIR}/include/gflags/gflags.h" AND EXISTS "${GFLAGS_DIR}/include/gflags/gflags.h")
			SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${GFLAGS_DIR}/include")			# "C:/Program Files/gflags/include"
			SET(LIBS_FILES ${LIBS_FILES} 
				optimized "${INTERNAL_LIBRARIES}/Release/gflags_static.lib"			# install generate only Release or Debug, not both libs
				debug "${INTERNAL_LIBRARIES}/Debug/gflags_static.lib" )				# so, it is generated in advance
		ELSE()
			MESSAGE(WARNING "Setup the directory where is installed GFLAGS. It is usually \"C:/Program Files/gflags\" and contains \"./include/gflags/gflags.h\" etc.")
		ENDIF()
	ELSE()
		MESSAGE(WARNING "Setup the directory where is installed GFLAGS. It is usually \"C:/Program Files/gflags\" and contains \"./include/gflags/gflags.h\" etc.")
	ENDIF()
ENDIF()


# -----------------------------------------------------------------------------
# SuiteSparse 
# -----------------------------------------------------------------------------
# For this version of Suitesparse do not exist the configuration file

SET(INCLUDE_DIRS ${INCLUDE_DIRS} 
	"${SUITESPARSE_ROOT_DIR}/install/include/suitesparse"        
	"${SUITESPARSE_ROOT_DIR}/install/include"
	"${SUITESPARSE_ROOT_DIR}/SuiteSparse/CXSparse/Include" )      

IF(WIN32)
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
ELSE()
  
ENDIF()


# -----------------------------------------------------------------------------
# CERES
# -----------------------------------------------------------------------------
find_package( Ceres )
IF ( CERES_FOUND )
	message( STATUS "Found CERES ${CERES_VERSION}" )
	message( STATUS "> CERES_INCLUDE_DIRS:   ${CERES_INCLUDE_DIRS}"   )
	message( STATUS "> CERES_LIBRARIES: ${CERES_LIBRARIES}" )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${CERES_INCLUDE_DIRS})
	SET(LIBS_FILES ${LIBS_FILES} ${CERES_LIBRARIES})

ELSE( CERES_FOUND )
	IF(WIN32)
		MESSAGE( STATUS "Ceres will be included manualy...")
		MESSAGE( STATUS "> include paths: ./config, ./include, ./internal, ./internal/ceres and ./ceres/miniglog" )
		MESSAGE( STATUS "> Release added library: ${INTERNAL_LIBRARIES}/Release/ceres.lib" )
		MESSAGE( STATUS "> Debug added library: ${INTERNAL_LIBRARIES}/Debug/ceres-debug.lib" )
		SET(INCLUDE_DIRS ${INCLUDE_DIRS} 
			"${CERES_ROOT_DIR}/config"        
			"${CERES_ROOT_DIR}/include"
			"${CERES_ROOT_DIR}/internal"
			"${CERES_ROOT_DIR}/internal/ceres"
			"${CERES_ROOT_DIR}/ceres/miniglog") 
		SET(LIBS_FILES ${LIBS_FILES} 
			optimized "${INTERNAL_LIBRARIES}/Release/ceres.lib"
			debug "${INTERNAL_LIBRARIES}/Debug/ceres-debug.lib")
	ELSEIF( UNIX )
		MESSAGE( WARNING "First, you should install Ceres libraries. Follow this guideline: http://ceres-solver.org/installation.html")
	ENDIF()
ENDIF( CERES_FOUND )



# -----------------------------------------------------------------------------
# MATLAB interface - MEX libs ( C:/Program Files/MATLAB/R2017a/extern/include; C:/Program Files/MATLAB/R2017a/extern/lib/win64/microsoft )
# -----------------------------------------------------------------------------
# MANUALLY REWRITTE THE CERES -> CERES DEBUG FOR DEBUG MODE  ..\libs\compiled_libraries\windows\Debug\ceres-debug.lib

option(MATLAB_SUPPORT "Include MATLAB libraries (enable input from MATLAB)" OFF)
IF (MATLAB_SUPPORT)
	add_definitions(-DUSE_MATLAB)

  SET(MEX_INCLUDE_DIR CACHE PATH FORCE)
  IF(NOT MEX_INCLUDE_DIR STREQUAL "")
  	MESSAGE(STATUS "MEX_INCLUDE_DIR: ${MEX_INCLUDE_DIR}")
  	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${MEX_INCLUDE_DIR})     
  ELSE()
  	MESSAGE(WARNING "Setup the include directory for MATLAB MEX interface (contatins mex.h).")
  ENDIF()
  SET(MEX_LIBRARIES_DIR CACHE PATH FORCE)
  IF(NOT MEX_LIBRARIES_DIR STREQUAL "")
  	MESSAGE(STATUS "MEX_LIBRARIES_DIR: ${MEX_LIBRARIES_DIR}")
  	IF(WIN32)
  		SET(LIBS_FILES ${LIBS_FILES} 
  		"${MEX_LIBRARIES_DIR}/libmx.lib"
  		"${MEX_LIBRARIES_DIR}/libmat.lib"
  		"${MEX_LIBRARIES_DIR}/libeng.lib"
  		"${MEX_LIBRARIES_DIR}/libmex.lib")     
  	ELSE()
  		SET(LIBS_FILES ${LIBS_FILES} 
  		"${MEX_LIBRARIES_DIR}/libmx.so"
  		"${MEX_LIBRARIES_DIR}/libmat.so"
  		"${MEX_LIBRARIES_DIR}/libeng.so"
  		"${MEX_LIBRARIES_DIR}/libmex.so")
  	ENDIF()
  ELSE()
  	MESSAGE(WARNING "Setup the library directory for MATLAB MEX interface (contatins libmx.lib, libmat.lib, libeng.lib, libmex.lib).")
  ENDIF() 
  
ENDIF()


# -----------------------------------------------------------------------------
# MAGMA
# -----------------------------------------------------------------------------
# For MAGMA do not exist the configuration file
IF(WIN32)
  SET(INCLUDE_DIRS ${INCLUDE_DIRS}
  	"${MAGMA_ROOT_DIR}/include"
  	"${MAGMA_ROOT_DIR}/control"
  	"${MAGMA_ROOT_DIR}/sparse/include"
  	"${MAGMA_ROOT_DIR}/testing")
 
	SET(LIBS_FILES ${LIBS_FILES}
  	optimized "${INTERNAL_LIBRARIES}/Release/tester.lib"
  	optimized "${INTERNAL_LIBRARIES}/Release/lapacktest.lib"
  	optimized "${INTERNAL_LIBRARIES}/Release/magma.lib"
  	optimized "${INTERNAL_LIBRARIES}/Release/magma_sparse.lib"
  	debug "${INTERNAL_LIBRARIES}/Debug/tester.lib"
  	debug "${INTERNAL_LIBRARIES}/Debug/lapacktest.lib"
  	debug "${INTERNAL_LIBRARIES}/Debug/magma.lib"
  	debug "${INTERNAL_LIBRARIES}/Debug/magma_sparse.lib")
   
ELSE() 
  SET(MAGMA_ROOT_DIR CACHE PATH FORCE)
	IF(NOT MAGMA_ROOT_DIR STREQUAL "")
		MESSAGE(STATUS "MAGMA_ROOT_DIR: ${MAGMA_ROOT_DIR}")
		SET(INCLUDE_DIRS ${INCLUDE_DIRS}
    	"${MAGMA_ROOT_DIR}/include"
    	"${MAGMA_ROOT_DIR}/control"
    	"${MAGMA_ROOT_DIR}/sparse/include"
      "${MAGMA_ROOT_DIR}/sparse/control"
    	"${MAGMA_ROOT_DIR}/testing")   
	ELSE()
		MESSAGE(WARNING "Setup the MAGMA_ROOT_DIR")
	ENDIF()  
 
  SET(MAGMA_INCLUDE_DIR CACHE PATH FORCE)
	IF(NOT MAGMA_INCLUDE_DIR STREQUAL "")
		MESSAGE(STATUS "MAGMA_INCLUDE_DIR: ${MAGMA_INCLUDE_DIR}")
		SET(LIBS_FILES ${LIBS_FILES}
    	"${MAGMA_INCLUDE_DIR}/libtester.a"
    	"${MAGMA_INCLUDE_DIR}/liblapacktest.a"
    	"${MAGMA_INCLUDE_DIR}/libmagma.a"
    	"${MAGMA_INCLUDE_DIR}/libmagma_sparse.a") 
	ELSE()
		MESSAGE(WARNING "Setup the MAGMA_INCLUDE_DIR")
	ENDIF()  
ENDIF()


# -----------------------------------------------------------------------------
# CUDA
# -----------------------------------------------------------------------------
find_package( CUDA )
IF ( CUDA_FOUND )
	message( STATUS "Found CUDA ${CUDA_VERSION}" )
	message( STATUS "> CUDA_INCLUDE_DIRS:   ${CUDA_INCLUDE_DIRS}"   )
	message( STATUS "> CUDA_LIBRARIES: ${CUDA_LIBRARIES}" )
	message( STATUS "> CUDA_CUBLAS_LIBRARIES: ${CUDA_CUBLAS_LIBRARIES}" )
	message( STATUS "> CUDA_CUDART_LIBRARY: ${CUDA_cusparse_LIBRARY}" )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${CUDA_INCLUDE_DIRS})
	SET(LIBS_FILES ${LIBS_FILES} ${CUDA_LIBRARIES};${CUDA_CUBLAS_LIBRARIES};${CUDA_cusparse_LIBRARY})

  message("LIBS_FILES: ${LIBS_FILES}")

ELSE( CUDA_FOUND )
	IF(WIN32)
		message( STATUS "CUDA 8.0 will be included manualy...")
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
	ELSEIF( UNIX )
		MESSAGE( WARNING "Install Cuda libraries first. Follow this guideline: https://askubuntu.com/questions/799184/how-can-i-install-cuda-on-ubuntu-16-04")
	ENDIF()
ENDIF( CUDA_FOUND )


# -----------------------------------------------------------------------------
# IntelMKL -> BLAS & LAPACK ( MAGMA was tested with this implementation of BLAS and LAPACK )
# -----------------------------------------------------------------------------
# For IntelMKL do not exist the configuration file
SET(INCLUDE_DIRS ${INCLUDE_DIRS} "${MKL_ROOT_DIR}/include")
IF(WIN32)
	SET(LIBS_FILES ${LIBS_FILES}
		"${INTERNAL_LIBRARIES}/Release/mkl_intel_lp64_dll.lib"
		"${INTERNAL_LIBRARIES}/Release/mkl_intel_thread_dll.lib"
		"${INTERNAL_LIBRARIES}/Release/mkl_core_dll.lib"
		"${INTERNAL_LIBRARIES}/Release/libiomp5md.lib")
ENDIF()



# -----------------------------------------------------------------------------
# OpenMVG
# -----------------------------------------------------------------------------
option(OPENMVG_SUPPORT "Include OpenMVG libraries (enable input from OpenMVG)" OFF)
IF (OPENMVG_SUPPORT)
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
    ELSE()
      #SET(LIBS_FILES ${LIBS_FILES} 
      #  "${OPENMVG_BUILD_DIR}/third_party/flann/lib/flann_cpp_s.a")  
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
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# INCLUDE ALL DIRECTORIES
# -----------------------------------------------------------------------------
target_include_directories(uncertainty PUBLIC ${INCLUDE_DIRS})
#target_include_directories(uncertainty_lib PUBLIC ${INCLUDE_DIRS})
#target_include_directories(uncertainty_mex PUBLIC ${INCLUDE_DIRS})

# -----------------------------------------------------------------------------
# ADD ALL LIBRARIES
# -----------------------------------------------------------------------------
target_link_libraries (uncertainty ${LIBS_FILES})
#target_link_libraries (uncertainty_lib ${LIBS_FILES})
#target_link_libraries (uncertainty_mex ${LIBS_FILES})

# -----------------------------------------------------------------------------
# Copy the examples and dll files if not exist
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


# TODO: install library and generate UncertaintyConfig.cmake
