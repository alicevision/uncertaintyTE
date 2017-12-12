# -----------------------------------------------------------------------------
# Set sources and create executable target
# -----------------------------------------------------------------------------
FILE(GLOB HEADERS "uncertaintyTE/*.h")
FILE(GLOB SOURCES "uncertaintyTE/*.cpp")

FILE(GLOB MODUL_DEFINITION_FILES "${PROJECT_SOURCE_DIR}/src/*.def")
SET(LIBS_FILES ${LIBS_FILES}) 
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# -----------------------------------------------------------------------------
# Library output
# -----------------------------------------------------------------------------
OPTION(UNC_BUILD_STATIC_LIB "Build uncertainty as static library" OFF)
OPTION(UNC_BUILD_SHARED_LIB "Build uncertainty as shared library" ON)

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
ADD_EXECUTABLE(uncertainty software/main_uncertainty.cpp)
target_include_directories(uncertainty
    PUBLIC
    .
    ${INCLUDE_DIRS}
    ${GFLAGS_INCLUDE_DIR})

target_link_libraries(uncertainty
    ${LIBS_FILES}
    ${GFLAGS_LIBRARIES})

if(UNC_BUILD_STATIC_LIB)
    target_link_libraries (uncertainty uncertainty_static)
else()
    target_link_libraries (uncertainty uncertainty_shared)
endif()


install(DIRECTORY uncertaintyTE
  DESTINATION include
  COMPONENT headers
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
  )



# ==============================================================================
# UNIT TEST framework
# ==============================================================================
if(UNC_BUILD_TESTS)
  enable_testing()
endif()

# MACRO to ease UnitTesting
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
			 COMMAND $<TARGET_FILE:${NAMESPACE}_test_${NAME}> --log_level=all)
  endif()
endmacro()


# Tests 
add_definitions(-DTESTS_ROOT_DIR="${CMAKE_SOURCE_DIR}/test")
add_subdirectory(test_eigValues)




