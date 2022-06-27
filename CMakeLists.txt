# SPDX-License-Identifier: LGPL-3.0-only

cmake_minimum_required(VERSION 3.8)
include(ExternalProject)

# When Radler is compiled as an ExternalProject inside another project, set this
# option to On. See, e.g., the wsclean CMake file for an example.
option(COMPILE_AS_EXTERNAL_PROJECT OFF)
option(PORTABLE "Generate portable code" OFF)
option(BUILD_PYTHON_BINDINGS "Build python bindings" OFF)
option(BUILD_DOCUMENTATION "Build documentation" OFF)
option(BUILD_DOCSTRINGS "Build docstrings (implies BUILD_DOCUMENTATION)" OFF)

set(RADLER_VERSION 0.0.0)
project(radler VERSION ${RADLER_VERSION})
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

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckCXXCompilerFlag)
if(PORTABLE)
  if(DEFINED TARGET_CPU)
    message(WARNING "You have selected to build PORTABLE binaries. "
                    "TARGET_CPU settings will be ignored.")
    unset(TARGET_CPU CACHE)
  endif()
endif()

if(NOT COMPILE_AS_EXTERNAL_PROJECT)
  # Include submodules
  set(ExternalSubmoduleDirectories aocommon pybind11 schaapcommon)
  foreach(ExternalSubmodule ${ExternalSubmoduleDirectories})
    if(NOT EXISTS ${CMAKE_SOURCE_DIR}/external/${ExternalSubmodule})
      message(
        FATAL_ERROR
          "The external submodule '${ExternalSubmodule}' is missing in the external/ subdirectory. "
          "This is likely the result of downloading a git tarball without submodules. "
          "This is not supported: git tarballs do not provide the required versioning "
          "information for the submodules. Please perform a git clone of this repository."
      )
    endif()
  endforeach()

  # Find and include git submodules
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    # Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
      message(STATUS "Submodule update")
      execute_process(
        COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --checkout
                --depth 1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GIT_SUBMOD_RESULT)
      if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(
          FATAL_ERROR
            "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules"
        )
      endif()
    endif()
  endif()
endif()

if(COMPILE_AS_EXTERNAL_PROJECT)
  message(
    STATUS "Radler is compiled as an external project within another project.")
  if(NOT DEFINED AOCOMMON_INCLUDE_DIR)
    message(
      FATAL_ERROR
        "AOCOMMON_INCLUDE_DIR not specified. Please add -DAOCOMMON_INCLUDE_DIR to the CMAKE_ARGS."
    )
  endif()
  if(NOT DEFINED SCHAAPCOMMON_SOURCE_DIR)
    message(
      FATAL_ERROR
        "SCHAAPCOMMON_SOURCE_DIR not specified. Please add -DSCHAAPCOMMON_SOURCE_DIR to the CMAKE_ARGS."
    )
  endif()
  if(NOT DEFINED PYBIND11_SOURCE_DIR)
    message(
      FATAL_ERROR
        "PYBIND11_SOURCE_DIR not specified. Please add -DPYBIND11_SOURCE_DIR to the CMAKE_ARGS."
    )
  endif()
else()
  set(AOCOMMON_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/aocommon/include)
  set(SCHAAPCOMMON_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/schaapcommon)
  set(PYBIND11_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/pybind11)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Boost date_time is needed in aocommon
find_package(
  Boost
  COMPONENTS date_time
  REQUIRED)

# Threads
find_package(Threads REQUIRED)

# Find and include HDF5
find_package(
  HDF5
  COMPONENTS C CXX
  REQUIRED)
add_definitions(${HDF5_DEFINITIONS} -DH5_USE_110_API)

# Casacore
set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
find_package(Casacore REQUIRED COMPONENTS fits casa ms tables measures)

# CFitsio
find_package(CFITSIO REQUIRED)

# Python3
find_package(PythonInterp REQUIRED) # Call before finding PythonLibs!
find_package(PythonLibs 3 REQUIRED)
message(STATUS "Using python version ${PYTHON_VERSION_STRING}")

if(COMPILE_AS_EXTERNAL_PROJECT)
  add_subdirectory(${SCHAAPCOMMON_SOURCE_DIR}
                   ${PROJECT_BINARY_DIR}/schaapcommon)
  add_subdirectory(${PYBIND11_SOURCE_DIR} ${PROJECT_BINARY_DIR}/pybind11)
else()
  add_subdirectory(${SCHAAPCOMMON_SOURCE_DIR})
  add_subdirectory(${PYBIND11_SOURCE_DIR})
endif()
target_include_directories(schaapcommon SYSTEM PUBLIC ${AOCOMMON_INCLUDE_DIR})

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
    ${HDF5_INCLUDE_DIRS}
    ${pybind11_INCLUDE_DIRS}
    ${SCHAAPCOMMON_SOURCE_DIR}/include)

set(RADLER_TARGET_LIBS
    ${Boost_DATE_TIME_LIBRARY} ${CMAKE_THREAD_LIBS_INIT} ${CASACORE_LIBRARIES}
    ${CFITSIO_LIBRARY} pybind11::embed schaapcommon)

# Source directories
add_subdirectory(cpp)

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

# Compile tests
if(NOT ${COMPILE_AS_EXTERNAL_PROJECT} AND ${BUILD_TESTING})
  include(CTest)
  find_package(
    Boost
    COMPONENTS unit_test_framework
    REQUIRED)

  # Set up a test_data directory in the cpp/test directory.
  set(DATA_DIR ${CMAKE_BINARY_DIR}/cpp/test/test_data)
  file(MAKE_DIRECTORY ${DATA_DIR})

  add_subdirectory(cpp/test)
  add_subdirectory(cpp/algorithms/test)
  add_subdirectory(cpp/math/test)
  # ------------------------------------------------------------------------------
  # Generate test_config.h header
  configure_file(${CMAKE_SOURCE_DIR}/cmake/config/test_config.h.in
                 ${CMAKE_BINARY_DIR}/test_config.h)

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
