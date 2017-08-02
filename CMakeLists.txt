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
ENDIF()

# -----------------------------------------------------------------------------
# Set sources and create executable target
file(GLOB HEADERS "${PROJECT_SOURCE_DIR}/include/*.h")
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
add_executable(uncertainty ${HEADERS} ${SOURCES})

# Create also the static library
ADD_LIBRARY(uncertainty_lib STATIC ${HEADERS} ${SOURCES})

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
SET(MAGMA_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/magma-2.2.0")
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
# EIGEN
# -----------------------------------------------------------------------------
find_package (Eigen3 NO_MODULE)

IF (EIGEN3_FOUND)
	message( STATUS "Found Eigen${EIGEN3_VERSION_STRING}" )
	message( STATUS "> EIGEN3_INCLUDE_DIRS:   ${EIGEN3_INCLUDE_DIRS}"   )
	SET(INCLUDE_DIRS ${INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS})
	SET(LIBS_FILES ${LIBS_FILES} "Eigen3::Eigen")
	
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


# -----------------------------------------------------------------------------
# CERES
# -----------------------------------------------------------------------------
find_package( CERES )
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
# For MATLAB interface do not exist the configuration file

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
	SET(LIBS_FILES ${LIBS_FILES} 
		"${MEX_LIBRARIES_DIR}/libmx.lib"
		"${MEX_LIBRARIES_DIR}/libmat.lib"
		"${MEX_LIBRARIES_DIR}/libeng.lib"
		"${MEX_LIBRARIES_DIR}/libmex.lib")     
ELSE()
	MESSAGE(WARNING "Setup the library directory for MATLAB MEX interface (contatins libmx.lib, libmat.lib, libeng.lib, libmex.lib).")
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
	SET(LIBS_FILES ${LIBS_FILES} ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_cusparse_LIBRARY})

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
# MAGMA
# -----------------------------------------------------------------------------
# For MAGMA do not exist the configuration file
SET(INCLUDE_DIRS ${INCLUDE_DIRS}
	"${MAGMA_ROOT_DIR}/include"
	"${MAGMA_ROOT_DIR}/control"
	"${MAGMA_ROOT_DIR}/sparse/include"
	"${MAGMA_ROOT_DIR}/testing")
IF(WIN32)
	SET(LIBS_FILES ${LIBS_FILES}
	optimized "${INTERNAL_LIBRARIES}/Release/tester.lib"
	optimized "${INTERNAL_LIBRARIES}/Release/lapacktest.lib"
	optimized "${INTERNAL_LIBRARIES}/Release/magma.lib"
	optimized "${INTERNAL_LIBRARIES}/Release/magma_sparse.lib"
	debug "${INTERNAL_LIBRARIES}/Debug/tester.lib"
	debug "${INTERNAL_LIBRARIES}/Debug/lapacktest.lib"
	debug "${INTERNAL_LIBRARIES}/Debug/magma.lib"
	debug "${INTERNAL_LIBRARIES}/Debug/magma_sparse.lib")
ENDIF()


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
# -----------------------------------------------------------------------------
# LINUX - INCLUDE ALL LIBRARIES IN DIRECTORY
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# TODO: REPLACE THIS WITH APROPRIATE LINKING OF LIBRARIES
IF(NOT WIN32)
	FILE(GLOB LINUX_LIB_SEARCH "${INTERNAL_LIBRARIES}/*")
	SET(LIBS_FILES ${LIBS_FILES} ${LINUX_LIB_SEARCH} "-lpthread -lm")
ENDIF()




# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# INCLUDE ALL DIRECTORIES
# -----------------------------------------------------------------------------

target_include_directories(uncertainty PUBLIC ${INCLUDE_DIRS})
target_include_directories(uncertainty_lib PUBLIC ${INCLUDE_DIRS})

# -----------------------------------------------------------------------------
# ADD ALL LIBRARIES
# -----------------------------------------------------------------------------
target_link_libraries (uncertainty ${LIBS_FILES})
target_link_libraries (uncertainty_lib ${LIBS_FILES})


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