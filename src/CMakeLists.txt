# -----------------------------------------------------------------------------
# Set sources and create executable target
# -----------------------------------------------------------------------------
FILE(GLOB HEADERS "uncertaintyTE/*.h")
FILE(GLOB SOURCES "uncertaintyTE/*.cpp")

FILE(GLOB MODUL_DEFINITION_FILES "${PROJECT_SOURCE_DIR}/src/*.def")
SET(LIBS_FILES ${LIBS_FILES}) 

# -----------------------------------------------------------------------------
# Library output
# -----------------------------------------------------------------------------
OPTION(UNC_BUILD_STATIC_LIB "Build uncertainty as static library" OFF)
OPTION(UNC_BUILD_SHARED_LIB "Build uncertainty as shared library" OFF)

IF (UNC_BUILD_STATIC_LIB)
	ADD_LIBRARY(uncertainty_static STATIC ${HEADERS} ${SOURCES})
	target_include_directories(uncertainty_static PUBLIC ${INCLUDE_DIRS})
	target_link_libraries (uncertainty_static ${LIBS_FILES})
        install(TARGETS uncertainty_static
          DESTINATION lib
          EXPORT uncertainty-targets
        )
ENDIF (UNC_BUILD_STATIC_LIB)

IF (UNC_BUILD_SHARED_LIB)
	ADD_LIBRARY(uncertainty_shared SHARED ${HEADERS} ${SOURCES})
	target_include_directories(uncertainty_shared PUBLIC ${INCLUDE_DIRS})
	target_link_libraries (uncertainty_shared ${LIBS_FILES})
        install(TARGETS uncertainty_shared
          DESTINATION lib
          EXPORT uncertainty-targets
        )
ENDIF (UNC_BUILD_SHARED_LIB)

# Executable - binary file with txt input
IF(NOT UNC_BUILD_MEX)
    ADD_EXECUTABLE(uncertainty software/main_uncertainty.cpp)
    target_include_directories(uncertainty PUBLIC . ${INCLUDE_DIRS})
    target_link_libraries (uncertainty ${LIBS_FILES})
    if(UNC_BUILD_STATIC_LIB)
      target_link_libraries (uncertainty uncertainty_static)
    else()
      target_link_libraries (uncertainty uncertainty_shared)
    endif()
ENDIF(NOT UNC_BUILD_MEX)


install(DIRECTORY uncertaintyTE
  DESTINATION include
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
  )

