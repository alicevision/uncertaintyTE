CMAKE_MINIMUM_REQUIRED (VERSION 3.5)
set(CMAKE_VERBOSE_MAKEFILE ON)

# About the project
PROJECT (UncertaintyTE
  LANGUAGES C CXX
  VERSION 1.2.0)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "---------------------------------------------------------")
message(STATUS "UncertaintyTE framework v.${PROJECT_VERSION}")
message(STATUS "CMAKE_CONFIGURATION_TYPES: ${CMAKE_CONFIGURATION_TYPES}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "---------------------------------------------------------")

# Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed.")
endif()

# Add path to find_package modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")


option(UNC_MKL_SUPPORT "Build with Intel MKL (include BLAS and LAPACK)" OFF)
option(UNC_MATLAB_SUPPORT "Build with MATLAB support (enable input from MATLAB)" OFF)
option(UNC_ALICEVISION_SUPPORT "Build with AliceVision (enable input from AliceVision)" OFF)
option(UNC_BUILD_MEX "Build MEX extension for Matlab" OFF)
option(UNC_BUILD_TESTS "Build unit tests" OFF)


# -----------------------------------------------------------------------------
# REQUIRED AND OPTIONAL LIBRARIES
# -----------------------------------------------------------------------------
message(STATUS "REQUIRED LIBRARIES")
message(STATUS " Gflags")
message(STATUS " Eigen")
message(STATUS " Blas and Lapack")
message(STATUS " Cuda")
message(STATUS " Ceres")
message(STATUS " Magma dense")
message(STATUS "")
message(STATUS "OPTIONAL LIBRARIES")
message(STATUS " SuiteSparse")
message(STATUS " Matlab mex")
message(STATUS " AliceVision")
message(STATUS " Boost")
message(STATUS "---------------------------------------------------------")


# -----------------------------------------------------------------------------
# SYSTEM LIBRARY - FOR LINUX
# -----------------------------------------------------------------------------
if(UNIX)
  find_package (Threads)
  set(LIBS_FILES ${LIBS_FILES} ${CMAKE_THREAD_LIBS_INIT})
	
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif(OPENMP_FOUND)
endif(UNIX)


# -----------------------------------------------------------------------------
# PREBUILDED LIBRARIES FOR WINDOWS - SETUP PATHS
# -----------------------------------------------------------------------------
if(WIN32)
    if(EXISTS "${PROJECT_SOURCE_DIR}/libs")
      set(USE_PREBUILD_WIN64_LIBS "ON")
      message(STATUS "Use precompiled libraries.")

      # Third party paths
      set(INTERNAL_LIBRARIES "${PROJECT_SOURCE_DIR}/libs/compiled_libraries/windows")
      set(CERES_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/ceres-solver")
      set(CUDA8_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/cuda-8.0")
      set(EIGEN3_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/eigen-3.3.4")
      set(GFLAGS_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/gflags/BUILD_")
      set(MKL_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/IntelMKL")
      set(MAGMA_ROOT_DIR_PRE "${PROJECT_SOURCE_DIR}/libs/third_party/magma-2.2.0")
      set(SUITESPARSE_ROOT_DIR "${PROJECT_SOURCE_DIR}/libs/third_party/suitesparse-metis-for-windows")
    endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/libs")
endif(WIN32)


# -----------------------------------------------------------------------------
# Gflags
# -----------------------------------------------------------------------------
find_package(gflags REQUIRED)
if(DEFINED gflags_VERSION_STRING OR DEFINED GFLAGS_VERSION_STRING)
        # different cmake variable with lower or upper case exist
        set(GFLAGS_INCLUDE_DIR ${GFLAGS_INCLUDE_DIR} ${gflags_INCLUDE_DIR})
        set(GFLAGS_LIBRARIES ${GFLAGS_LIBRARIES} ${gflags_LIBRARIES})
  message( "Found Gflags ${GFLAGS_VERSION_STRING}" )
  message( STATUS "> GFLAGS_INCLUDE_DIR:   ${GFLAGS_INCLUDE_DIR}" )
  message( STATUS "> GFLAGS_LIBRARIES:   ${GFLAGS_LIBRARIES}" )
endif()

# SHLWAPI - Used by gflags on WIN32
#if(WIN32)
#find_library(SHLWAPI Shlwapi.lib) 
#set(GFLAGS_LIBRARIES ${GFLAGS_LIBRARIES} ${SHLWAPI})
#else()
#set(SHLWAPI)
#endif()

# -----------------------------------------------------------------------------
# EIGEN
# -----------------------------------------------------------------------------
find_package(Eigen3 3.3 NO_MODULE REQUIRED)
if(EIGEN3_FOUND)
  message( "Found Eigen${EIGEN3_VERSION_STRING}" )
  message( STATUS "> EIGEN3_INCLUDE_DIRS:   ${EIGEN3_INCLUDE_DIRS}"   )
  set(EIGEN_FLAGS "-DEIGEN_DONT_ALIGN_STATICALLY=1 -DEIGEN_DONT_VECTORIZE=1")
endif(EIGEN3_FOUND)

# -----------------------------------------------------------------------------
# SuiteSparse 
# -----------------------------------------------------------------------------
find_package(SuiteSparse REQUIRED)

if( SuiteSparse_FOUND)
  message( "Found SuiteSparse" )
  message( STATUS "> SUITESPARSE_INCLUDE_DIRS:   ${SUITESPARSE_INCLUDE_DIRS}" )
  message( STATUS "> SUITESPARSE_LIBRARIES:   ${SUITESPARSE_LIBRARIES}" )
endif(SuiteSparse_FOUND)


# -----------------------------------------------------------------------------
# CERES
# -----------------------------------------------------------------------------
find_package(Ceres REQUIRED)
if(CERES_FOUND)
  message( "Found CERES ${CERES_VERSION}" )
  message( STATUS "> CERES_INCLUDE_DIRS:   ${CERES_INCLUDE_DIRS}"   )
  message( STATUS "> CERES_LIBRARIES: ${CERES_LIBRARIES}" )
endif()


# -----------------------------------------------------------------------------
# IntelMKL -> BLAS & LAPACK
# -----------------------------------------------------------------------------
if(UNC_MKL_SUPPORT)
	find_package(MKL QUIET)
  if(MKL_FOUND)
    message( "Found Intel MKL" )
    message( STATUS "> MKL_INCLUDE_DIR: ${MKL_INCLUDE_DIR}"   )
    message( STATUS "> MKL_LIBRARIES: ${MKL_LIBRARIES}" )
  endif()
endif()


# -----------------------------------------------------------------------------
# LAPACK AND BLAS LIBS
# -----------------------------------------------------------------------------
if(NOT UNC_MKL_SUPPORT AND UNIX)
  find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND)
    message("Found lapack and blas: ${LAPACK_LIBRARIES}")
  endif()
endif()


# -----------------------------------------------------------------------------
# MATLAB interface - MEX libs ( C:/Program Files/MATLAB/R2017a/extern/include; C:/Program Files/MATLAB/R2017a/extern/lib/win64/microsoft )
# -----------------------------------------------------------------------------
if(UNC_MATLAB_SUPPORT)
  find_package(Matlab)
  if(Matlab_FOUND)
    add_definitions(-DUSE_MATLAB)
    message("Found Matlab libraries")
    message(STATUS "Matlab_INCLUDE_DIRS: ${Matlab_INCLUDE_DIRS}")
    message(STATUS "Matlab_LIBRARIES: ${Matlab_LIBRARIES}")
    set(INCLUDE_DIRS ${INCLUDE_DIRS} ${Matlab_INCLUDE_DIRS})
		GET_FILENAME_COMPONENT(MATLAB_LIBS_ROOT ${Matlab_LIBRARIES} DIRECTORY)  # ${Matlab_LIBRARIES} doesn't contain all required libs
    if(UNIX)
      set(Matlab_LIBRARIES ${Matlab_LIBRARIES}
				"${MATLAB_LIBS_ROOT}/libmx.so"
				"${MATLAB_LIBS_ROOT}/libmat.so"
				"${MATLAB_LIBS_ROOT}/libeng.so"
				"${MATLAB_LIBS_ROOT}/libmex.so")
    else()
      set(Matlab_LIBRARIES ${Matlab_LIBRARIES}
				"${MATLAB_LIBS_ROOT}/libmx.lib"
				"${MATLAB_LIBS_ROOT}/libmat.lib"
				"${MATLAB_LIBS_ROOT}/libeng.lib"
				"${MATLAB_LIBS_ROOT}/libmex.lib")
    endif()
  endif()
endif(UNC_MATLAB_SUPPORT)


# -----------------------------------------------------------------------------
# MAGMA
# -----------------------------------------------------------------------------
find_package(MAGMA REQUIRED)
if(MAGMA_FOUND)
  message("Found Magma")
  message(STATUS "MAGMA_INCLUDE_DIRS: ${MAGMA_INCLUDE_DIRS}")
  message(STATUS "MAGMA_LIBRARIES: ${MAGMA_LIBRARIES}")
endif(MAGMA_FOUND)


# -----------------------------------------------------------------------------
# CUDA
# -----------------------------------------------------------------------------
find_package(CUDA 7.0 REQUIRED)
if(CUDA_FOUND)
  message( "Found CUDA ${CUDA_VERSION}" )
  message( STATUS "> CUDA_INCLUDE_DIRS:   ${CUDA_INCLUDE_DIRS}"   )
  message( STATUS "> CUDA_LIBRARIES: ${CUDA_LIBRARIES}" )
  message( STATUS "> CUDA_CUBLAS_LIBRARIES: ${CUDA_CUBLAS_LIBRARIES}" )
  message( STATUS "> CUDA_CUDART_LIBRARY: ${CUDA_cusparse_LIBRARY}" )
endif()


# -----------------------------------------------------------------------------
# AliceVision
# -----------------------------------------------------------------------------
if(UNC_ALICEVISION_SUPPORT)
  find_package(AliceVision)
  if(AliceVision_FOUND)
    add_definitions(-DUSE_ALICEVISION)
  endif()
endif()


# -----------------------------------------------------------------------------
# Boost
# -----------------------------------------------------------------------------
find_package(Boost 1.55 REQUIRED)

if(Boost_FOUND)
  message("Found Boost")
  message(STATUS "> Boost_INCLUDE_DIRS:   ${Boost_INCLUDE_DIRS}")
  message(STATUS "> Boost_LIBRARIES:   ${Boost_LIBRARIES}")
else()
  message(WARNING "Failed to find the Boost library.")
endif()

# -----------------------------------------------------------------------------
# Add project source code
# -----------------------------------------------------------------------------
add_subdirectory(src)


# -----------------------------------------------------------------------------
# PLATFORM SPECIFIC ERRORS
# The debug on windows may need additional libraries
# -----------------------------------------------------------------------------
if(WIN32)
  set(INTEL_COMPILER_DEBUD_LIBRARIES_DIRECTORY "c:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017.4.210/windows/compiler/lib/intel64_win")
  message("For the debuging on windows platform may required an additional libraries (e.g. ifconsole.lib).")
  message("Please add the library path to your project (e.g. 'c:/Program Files (x86)/IntelSWTools/compilers_and_libraries_2017.4.210/windows/compiler/lib/intel64_win')")
endif()


# -----------------------------------------------------------------------------
# UNIT TEST macro
# -----------------------------------------------------------------------------
macro(UNIT_TEST NAMESPACE NAME EXTRA_LIBS)
  if(UNC_BUILD_TESTS)
    add_executable(${NAMESPACE}_test_${NAME} ${NAME}_test.cpp)

    set_property(TARGET ${NAMESPACE}_test_${NAME} PROPERTY FOLDER uncertaintyTE/test)

    target_link_libraries(${NAMESPACE}_test_${NAME}
                          ${EXTRA_LIBS} # Extra libs MUST be first.
                          ${BOOST_LIBRARIES} uncertainty_shared)

  target_include_directories(${NAMESPACE}_test_${NAME}
        PUBLIC
        .
        ${INCLUDE_DIRS})

    add_test(NAME ${NAMESPACE}_test_${NAME}
       WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
       COMMAND $<TARGET_file:${NAMESPACE}_test_${NAME}> --log_level=all)
  endif()
endmacro()

# -----------------------------------------------------------------------------
# STATUS
# -----------------------------------------------------------------------------
message(STATUS " ")
message(STATUS "---------------------------------------------------------")
message(STATUS "UNC_ALICEVISION_SUPPORT: ${UNC_ALICEVISION_SUPPORT}")
message(STATUS "UNC_MKL_SUPPORT: ${UNC_MKL_SUPPORT}")
message(STATUS "UNC_MATLAB_SUPPORT: ${UNC_MATLAB_SUPPORT}")
message(STATUS "UNC_BUILD_MEX: ${UNC_BUILD_MEX}")
message(STATUS "UNC_BUILD_TESTS: ${UNC_BUILD_TESTS}")
message(STATUS "---------------------------------------------------------")

