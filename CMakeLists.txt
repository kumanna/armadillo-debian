#
# (C) 2009 NICTA

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_aux/cmake/Modules/")

INCLUDE(CheckIncludeFile)
INCLUDE(CheckLibraryExists)
INCLUDE(CheckFunctionExists) 
#INCLUDE(CheckPrototypeExists)   # seems to be only present in cmake >= 2.6

#SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

PROJECT( armadillo CXX )

#CMAKE_REQUIRED_FLAGS = string of compile command line flags
#CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#CMAKE_REQUIRED_LIBRARIES = list of libraries to link

IF(WIN32)
  IF(NOT CYGWIN)
    MESSAGE(FATAL_ERROR "Sorry, automatic installation is currently not supported on this platform. Please try the manual installation.")
  ENDIF(NOT CYGWIN)
ENDIF(WIN32)


IF(APPLE)
  SET(ARMA_OS macos)

  SET(ARMA_USE_LAPACK true)
  SET(ARMA_USE_BLAS true)

ELSE(APPLE)
  SET(ARMA_OS unix)

  INCLUDE(ARMA_FindLAPACK)
  INCLUDE(ARMA_FindBLAS)

  IF(LAPACK_FOUND)
    SET(ARMA_USE_LAPACK true)
  ENDIF(LAPACK_FOUND)
  
  IF(BLAS_FOUND)
    SET(ARMA_USE_BLAS true)
  ENDIF(BLAS_FOUND)
ENDIF(APPLE)


INCLUDE(ARMA_FindCLAPACKATLAS)
INCLUDE(ARMA_FindCBLAS)

IF(CLAPACKATLAS_FOUND AND CBLAS_FOUND)
  IF(${CLAPACKATLAS_INCLUDE_DIR} STREQUAL ${CBLAS_INCLUDE_DIR})
    SET(ARMA_USE_ATLAS true)
    SET(ARMA_ATLAS_INCLUDE_DIR ${CLAPACKATLAS_INCLUDE_DIR}/)  # the trailing / character is critical
  ENDIF(${CLAPACKATLAS_INCLUDE_DIR} STREQUAL ${CBLAS_INCLUDE_DIR})
ENDIF(CLAPACKATLAS_FOUND AND CBLAS_FOUND)


INCLUDE(FindBoost)
FIND_PACKAGE(Boost)

IF(Boost_FOUND)
  SET(ARMA_USE_BOOST      true)
  SET(ARMA_USE_BOOST_BOOL true)
ELSE(Boost_FOUND)
  SET(ARMA_USE_BOOST_BOOL false)
ENDIF(Boost_FOUND)


#CHECK_PROTOTYPE_EXISTS(isfinite "cmath" ARMA_SYS_ISFINITE)
#CHECK_PROTOTYPE_EXISTS(isinf    "cmath" ARMA_SYS_ISINF)
#CHECK_PROTOTYPE_EXISTS(isnan    "cmath" ARMA_SYS_ISNAN)

SET(ARMA_LIBS ${ARMA_LIBS} ${CLAPACKATLAS_LIBRARIES})
SET(ARMA_LIBS ${ARMA_LIBS} ${CBLAS_LIBRARIES})
SET(ARMA_LIBS ${ARMA_LIBS} ${BLAS_LIBRARIES})
SET(ARMA_LIBS ${ARMA_LIBS} ${LAPACK_LIBRARIES})

IF(APPLE)
  SET(ARMA_LIBS ${ARMA_LIBS} "-framework Accelerate")  # or "-framework accelerate" ?
  MESSAGE(STATUS "MacOS X detected. Added '-framework Accelerate' to compiler flags")
ENDIF(APPLE)

MESSAGE(STATUS "ARMA_LIBS = ${ARMA_LIBS}")

MESSAGE(STATUS "Generating 'include/armadillo_bits/config.hpp'")
CONFIGURE_FILE(include/armadillo_bits/config.hpp.cmake include/armadillo_bits/config.hpp)

MESSAGE(STATUS "Generating 'examples/Makefile'")
CONFIGURE_FILE(examples/Makefile.cmake examples/Makefile)


IF(ARMA_USE_BOOST)
  INCLUDE_DIRECTORIES(include ${Boost_INCLUDE_DIR})
ELSE(ARMA_USE_BOOST)
  INCLUDE_DIRECTORIES(include)
ENDIF(ARMA_USE_BOOST)

ADD_LIBRARY( armadillo SHARED src/pull_libs )
TARGET_LINK_LIBRARIES( armadillo ${ARMA_LIBS} )

SET_TARGET_PROPERTIES(
armadillo 
PROPERTIES
VERSION 0.6.10
SOVERSION 0
)

# As Red Hat Enterprise Linux (and related systems such as Fedora)
# does not search /usr/local/lib by default, we need to place the
# library in /usr

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "/usr")
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

MESSAGE(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")


# Allow for the "lib" directory to be specified on the command line

IF(NOT LIB_INSTALL_DIR)
  IF(UNIX AND NOT APPLE)
  
    # The variable "CMAKE_SIZEOF_VOID_P" seems to be set only in CMake >= 2.6.0
    
    IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
      MESSAGE(STATUS "Detected 64 bit system")
      SET(LIB_INSTALL_DIR "lib64")
    ELSE(CMAKE_SIZEOF_VOID_P EQUAL 8)
      SET(LIB_INSTALL_DIR "lib")
    ENDIF(CMAKE_SIZEOF_VOID_P EQUAL 8)
  
  ELSE(UNIX AND NOT APPLE)
    SET(LIB_INSTALL_DIR "lib")
  ENDIF(UNIX AND NOT APPLE)
ENDIF(NOT LIB_INSTALL_DIR)

MESSAGE(STATUS "LIB_INSTALL_DIR = ${LIB_INSTALL_DIR}")

# Allow for the "include" directory to be specified on the command line

IF(NOT INCLUDE_INSTALL_DIR)
  SET(INCLUDE_INSTALL_DIR "include")
ENDIF(NOT INCLUDE_INSTALL_DIR)

MESSAGE(STATUS "INCLUDE_INSTALL_DIR = ${INCLUDE_INSTALL_DIR}")


# Note that the trailing / character in "include/" is critical

INSTALL(DIRECTORY include/ DESTINATION ${INCLUDE_INSTALL_DIR}
PATTERN ".svn" EXCLUDE
PATTERN "*.cmake" EXCLUDE
PATTERN "*~" EXCLUDE
PATTERN "*orig" EXCLUDE
)


INSTALL(TARGETS armadillo
LIBRARY DESTINATION ${LIB_INSTALL_DIR}
)



