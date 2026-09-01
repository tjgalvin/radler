# SPDX-License-Identifier: LGPL-3.0-only

cmake_minimum_required(VERSION 3.18)
include(ExternalProject)
include(FetchContent)

# When Radler is compiled as an ExternalProject inside another project, set this
# option to On. See, e.g., the wsclean CMake file for an example.
option(COMPILE_AS_EXTERNAL_PROJECT OFF)
option(PORTABLE "Build portable binaries (with slightly decreased performance)"
       OFF)
option(BUILD_PYTHON_BINDINGS "Build python bindings" OFF)
option(BUILD_DOCUMENTATION "Build documentation" OFF)
option(BUILD_DOCSTRINGS "Build docstrings (implies BUILD_DOCUMENTATION)" OFF)

set(RADLER_VERSION 1.0.0) # Keep in sync with `pyproject.toml` file
project(
  radler
  VERSION ${RADLER_VERSION}
  LANGUAGES CXX)
if(RADLER_VERSION MATCHES "^([0-9]+)\\.([0-9]+)\\.([0-9]+)")
  set(RADLER_VERSION_MAJOR "${CMAKE_MATCH_1}")
  set(RADLER_VERSION_MINOR "${CMAKE_MATCH_2}")
  set(RADLER_VERSION_PATCH "${CMAKE_MATCH_3}")
else()
  message(FATAL_ERROR "Failed to parse RADLER_VERSION='${RADLER_VERSION}'")
endif()

if(POLICY CMP0076)
  cmake_policy(SET CMP0076 NEW)
endif()
if(POLICY CMP0115)
  # Prevent warning about extensions, specifically about Doxyfile
  cmake_policy(SET CMP0115 NEW)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(
  -Wall
  -Wnon-virtual-dtor
  -Wzero-as-null-pointer-constant
  -Wduplicated-branches
  -Wundef
  -Wvla
  -Wpointer-arith
  -Wextra
  -Wno-unused-parameter)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-DNDEBUG)
endif()

set(SCHAAPCOMMON_MODULES fitters math threading)
if(COMPILE_AS_EXTERNAL_PROJECT)
  message(
    STATUS "Radler is compiled as an external project within another project.")
  if(NOT DEFINED AOCOMMON_INCLUDE_DIR OR NOT DEFINED AOCOMMON_CMAKE_DIR)
    message(
      FATAL_ERROR "AOCOMMON_INCLUDE_DIR or AOCOMMON_CMAKE_DIR is not specified."
    )
  endif()
  if(NOT DEFINED SCHAAPCOMMON_SOURCE_DIR)
    message(
      FATAL_ERROR
        "SCHAAPCOMMON_SOURCE_DIR not specified. Please add -DSCHAAPCOMMON_SOURCE_DIR to the CMAKE_ARGS."
    )
  endif()
else()
  FetchContent_Declare(
    aocommon
    GIT_REPOSITORY https://gitlab.com/aroffringa/aocommon.git
    GIT_TAG e1b7980e02211ade95eef176a3ac7f9ba3dcc55c)
  FetchContent_MakeAvailable(aocommon)
  set(AOCOMMON_INCLUDE_DIR "${aocommon_SOURCE_DIR}/include")
  set(AOCOMMON_CMAKE_DIR "${aocommon_SOURCE_DIR}/CMake")
  # While we include aocommon for the Radler target later on, Schaapcommon needs
  # the path included before fetching.
  include_directories(${AOCOMMON_INCLUDE_DIR})

  FetchContent_Declare(
    schaapcommon
    GIT_REPOSITORY https://git.astron.nl/RD/schaapcommon.git
    GIT_TAG fa27b39d4371f4575bac2fcfe17075a41776c679)

  FetchContent_MakeAvailable(schaapcommon)
  set(SCHAAPCOMMON_SOURCE_DIR "${schaapcommon_SOURCE_DIR}")
endif()

# User may optionally set `TARGET_CPU` if `PORTABLE=OFF`
include(${AOCOMMON_CMAKE_DIR}/SetTargetCPU.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Boost date_time is needed in aocommon. date_time is header-only since Boost
# version 1.73. (Using "COMPONENTS date_time" no longer works from Boost 1.89.)
find_package(Boost 1.73 REQUIRED)

# Threads
find_package(Threads REQUIRED)

# GSL
find_package(GSL REQUIRED)
if(NOT GSL_LIBRARY
   OR NOT GSL_INCLUDE_DIRS
   OR NOT GSL_CBLAS_LIBRARY)
  message(FATAL_ERROR "GSL not found, but required to build Radler!")
endif()

# Find and include HDF5
find_package(
  HDF5
  COMPONENTS CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS} -DH5_USE_110_API)

# Casacore
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore 3.6.0 REQUIRED COMPONENTS casa tables measures)

# CFitsio
find_package(CFITSIO REQUIRED)

# Python
if(DEFINED ENV{CIBUILDWHEEL})
  find_package(
    Python
    COMPONENTS Interpreter Development.Module
    REQUIRED)
else()
  find_package(
    Python
    COMPONENTS Interpreter Development
    REQUIRED)
endif()
message(STATUS "Using python version ${Python_VERSION}")

find_package(pybind11 CONFIG REQUIRED HINTS ${Python_SITELIB})

# When building the python bindings, set the rpath to the install directory, so
# that importing the python module doesn't require to set the LD_LIBRARY_PATH
# manually.
if(NOT ${COMPILE_AS_EXTERNAL_PROJECT} AND ${BUILD_PYTHON_BINDINGS})
  # Include GNUInstallDirs for CMAKE_INSTALL_FULL_LIBDIR include(GNUInstallDirs)
  include(GNUInstallDirs)
  # Use, i.e. don't skip the full RPATH for the build tree.
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  # When building, don't use the install RPATH already (but later on when
  # installing).
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  # Add the automatically determined parts of the RPATH which point to
  # directories outside the build tree to the install RPATH.
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # The RPATH to be used when installing, but only if it's not a system
  # directory.
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
       "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
  endif("${isSystemDir}" STREQUAL "-1")
endif()

set(RADLER_TARGET_INCLUDE_DIRS
    ${AOCOMMON_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${CASACORE_INCLUDE_DIRS}
    ${CFITSIO_INCLUDE_DIR}
    ${GSL_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
    ${SCHAAPCOMMON_SOURCE_DIR}/include)

set(RADLER_TARGET_LIBS
    ${CMAKE_THREAD_LIBS_INIT}
    ${CASACORE_LIBRARIES}
    ${CFITSIO_LIBRARY}
    pybind11::embed
    schaapcommon
    schaapthreading
    GSL::gsl
    GSL::gslcblas)

# Source directories
add_subdirectory(cpp)

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

# Compile tests
if(NOT ${COMPILE_AS_EXTERNAL_PROJECT} AND ${BUILD_TESTING})
  include(CTest)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework)

  # Set up a test_data directory in the cpp/test directory.
  set(TEST_DATA_DIR ${CMAKE_BINARY_DIR}/cpp/test/test_data)
  file(MAKE_DIRECTORY ${TEST_DATA_DIR})

  add_subdirectory(cpp/test)
  add_subdirectory(cpp/algorithms/test)
  add_subdirectory(cpp/math/test)
  add_subdirectory(cpp/utils/test)

  if(BUILD_PYTHON_BINDINGS)
    add_subdirectory(python/test)
  endif()

  # demo targets are excluded from a default build
  add_subdirectory(cpp/demo)
endif()

if(BUILD_DOCSTRINGS)
  set(BUILD_DOCUMENTATION "ON")
endif()
if(BUILD_DOCUMENTATION)
  add_subdirectory(doc)
endif()
