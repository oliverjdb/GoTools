PROJECT(GoTesselation)

# find packages
# FIND_PACKAGE (Eigen3 3.3 REQUIRED NO_MODULE)


# Include directories

INCLUDE_DIRECTORIES(
  ${GoTesselation_SOURCE_DIR}/include
  ${GoToolsCore_SOURCE_DIR}/include
  ${GoTrivariateModel_SOURCE_DIR}/include
  ${GoCompositeModel_SOURCE_DIR}/include
  ${GoTopology_SOURCE_DIR}/include
  ${GoTrivariate_SOURCE_DIR}/include
  ${parametrization_SOURCE_DIR}/include
  ${GoIgeslib_SOURCE_DIR}/include  
  ${GoTools_COMMON_INCLUDE_DIRS}
  /usr/local/include/eigen3/
)


# Linked in libraries

SET(DEPLIBS
  GoToolsCore
  GoTrivariateModel
  GoCompositeModel
  GoTopology
  GoTrivariate
  GoIgeslib
  parametrization
  ttl
  )

# Make the GoTesselation library

FILE(GLOB_RECURSE GoTesselation_SRCS src/*.C include/*.h)
if (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTesselation SHARED ${GoTesselation_SRCS})
else (BUILD_AS_SHARED_LIBRARY)
    ADD_LIBRARY(GoTesselation ${GoTesselation_SRCS})
endif (BUILD_AS_SHARED_LIBRARY)
TARGET_LINK_LIBRARIES(GoTesselation ${DEPLIBS})
SET_PROPERTY(TARGET GoTesselation
  PROPERTY FOLDER "GoTesselation/Libs")
SET_TARGET_PROPERTIES(GoTesselation PROPERTIES SOVERSION ${GoTools_ABI_VERSION})

# Apps and tests
MACRO(ADD_APPS SUBDIR PROPERTY_FOLDER IS_TEST)
  FILE(GLOB_RECURSE GoTesselation_APPS ${SUBDIR}/*.C)
  FOREACH(app ${GoTesselation_APPS})
    GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
    ADD_EXECUTABLE(${appname} ${app})
    TARGET_LINK_LIBRARIES(${appname} GoTesselation ${DEPLIBS})
    SET_TARGET_PROPERTIES(${appname}
      PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${SUBDIR})
    SET_PROPERTY(TARGET ${appname}
      PROPERTY FOLDER "GoTesselation/${PROPERTY_FOLDER}")
    IF(${IS_TEST})
      ADD_TEST(${appname} ${SUBDIR}/${appname}
		--log_format=XML --log_level=all --log_sink=../Testing/${appname}.xml)
    ENDIF(${IS_TEST})

    # Get rid of -Wconversion.  This is because annoying warning from the CUDA code in Eigen
    # TARGET_COMPILE_OPTIONS(${appname} PUBLIC -std=c++11 -Wall -Wno-unused-but-set-variable -fPIC)

  ENDFOREACH(app)
ENDMACRO(ADD_APPS)

IF(GoTools_COMPILE_APPS)
  ADD_APPS(app "Apps" FALSE)
  ADD_APPS(examples "Examples" FALSE)
ENDIF(GoTools_COMPILE_APPS)

IF(GoTools_COMPILE_TESTS)
  SET(DEPLIBS ${DEPLIBS} ${Boost_LIBRARIES})
  ADD_APPS(unit_tests "Unit Tests" TRUE)
  ADD_APPS(regression_tests "Regression Tests" TRUE)
ENDIF(GoTools_COMPILE_TESTS)


# Copy data
FILE(COPY ${GoTesselation_SOURCE_DIR}/data
  DESTINATION ${GoTesselation_BINARY_DIR})
# FILE(COPY ${GoTesselation_SOURCE_DIR}/examples/data
#   DESTINATION ${GoTesselation_BINARY_DIR}/examples)

# 'install' target

IF(WIN32)
  # Windows
  # lib
  INSTALL(TARGETS GoTesselation DESTINATION ${GoTools_INSTALL_PREFIX}/lib)
  #include
  INSTALL(DIRECTORY include/GoTools/tesselation
    DESTINATION ${GoTools_INSTALL_PREFIX}/include/GoTools
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  #lib
  INSTALL(TARGETS GoTesselation DESTINATION lib COMPONENT tesselation)
  # include
  INSTALL(DIRECTORY include/GoTools/tesselation
    DESTINATION include/GoTools
    COMPONENT tesselation-dev
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)

SET(CPACK_STRIP_FILES ${CPACK_STRIP_FILES} libGoTesselation.so)

