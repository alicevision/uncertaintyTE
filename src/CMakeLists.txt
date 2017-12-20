# -----------------------------------------------------------------------------
# Set sources and create executable target
# -----------------------------------------------------------------------------

file(GLOB MODUL_DEFINITION_FILES "*.def")
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# ==============================================================================
# Uncertainty Library
# ==============================================================================

file(GLOB UNCERTAINTYTE_HEADERS "uncertaintyTE/*.h")
file(GLOB UNCERTAINTYTE_SOURCES "uncertaintyTE/*.cpp")

add_library(uncertaintyTE SHARED ${UNCERTAINTYTE_HEADERS} ${UNCERTAINTYTE_SOURCES})
target_include_directories(uncertaintyTE
  PUBLIC
    ${EIGEN3_INCLUDE_DIRS}
    ${SUITESPARSE_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${MKL_INCLUDE_DIR} # optional
    ${Matlab_INCLUDE_DIRS} # optional
    ${MEX_INCLUDE_DIR} # optional
    ${MAGMA_INCLUDE_DIRS}
    ${CUDA_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
  )
target_link_libraries(uncertaintyTE
  ${EIGEN3_INCLUDE_DIRS}
  ${SUITESPARSE_LIBRARIES}
  ${CERES_LIBRARIES}
  ${MKL_LIBRARIES} # optional
  ${Matlab_LIBRARIES} # optional
  ${MAGMA_LIBRARIES}
  ${CUDA_LIBRARIES}
  ${Boost_LIBRARIES}
  )
target_compile_definitions(uncertaintyTE
  PUBLIC
    ${EIGEN_FLAGS}
  )
install(TARGETS uncertaintyTE
    DESTINATION lib
    EXPORT uncertaintyTE-targets
  )

# ==============================================================================
# Uncertainty command line application
# ==============================================================================
add_executable(uncertainty software/main_uncertainty.cpp)
target_include_directories(uncertainty
    PUBLIC
    .
    ${INCLUDE_DIRS}
    ${GFLAGS_INCLUDE_DIR})

target_link_libraries(uncertainty
  uncertaintyTE
  ${GFLAGS_LIBRARIES}
  )

install(DIRECTORY uncertaintyTE
  DESTINATION include
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
  )


# ==============================================================================
# MATLAB extension
# ==============================================================================
if(UNC_BUILD_MEX)
  find_package(Matlab)
  matlab_add_mex(
    NAME uncertainty_mex
    SRC matlab/mex_uncertainty.cpp
    LINK_TO uncertaintyTE)
endif(UNC_BUILD_MEX)

# ==============================================================================
# UNIT TEST framework
# ==============================================================================
if(UNC_BUILD_TESTS)
  enable_testing()
endif()


# Tests 
add_definitions(-DTESTS_ROOT_DIR="${CMAKE_SOURCE_DIR}/test")
add_subdirectory(test_eigValues)




