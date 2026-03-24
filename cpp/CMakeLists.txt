# SPDX-License-Identifier: LGPL-3.0-only

# CMake function to keep directory structure when installing headers.
function(install_headers_with_directory HEADER_LIST)
  foreach(HEADER ${HEADER_LIST})
    string(REGEX MATCH ".*\/" DIR ${HEADER})
    install(FILES ${HEADER}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/include/radler/${DIR})
  endforeach(HEADER)
endfunction(install_headers_with_directory)

add_library(radler STATIC "")
set_target_properties(radler PROPERTIES POSITION_INDEPENDENT_CODE True)

set(RADLER_FILES
    radler.cc
    component_list.cc
    image_set.cc
    algorithms/asp_algorithm.cc
    algorithms/deconvolution_algorithm.cc
    algorithms/generic_clean.cc
    algorithms/iuwt_deconvolution_algorithm.cc
    algorithms/ls_deconvolution.cc
    algorithms/more_sane.cc
    algorithms/multiscale_algorithm.cc
    algorithms/parallel_deconvolution.cc
    algorithms/python_deconvolution.cc
    algorithms/simple_clean.cc
    algorithms/subminor_loop.cc
    algorithms/threaded_deconvolution_tools.cc
    algorithms/iuwt/image_analysis.cc
    algorithms/iuwt/iuwt_decomposition.cc
    algorithms/iuwt/iuwt_mask.cc
    algorithms/multiscale/multiscale_transforms.cc
    math/dijkstra_splitter.cc
    math/peak_finder.cc
    math/component_optimization.cc
    math/rms_image.cc
    utils/casa_mask_reader.cc
    utils/compressed_mask.cc
    work_table.cc)

# A number of files perform the 'core' high-performance floating point
# operations. In these files, NaNs are avoided and thus -ffast-math is allowed.
# Note that visibilities can be NaN hence this can not be turned on for all
# files.
set_source_files_properties(
  image_set.cpp
  algorithms/generic_clean.cpp
  algorithms/multiscale_algorithm.cpp
  algorithms/threaded_deconvolution_tools.cpp
  algorithms/simple_clean.cpp
  algorithms/subminor_loop.cpp
  algorithms/multiscale/multiscale_transforms.cpp
  PROPERTIES COMPILE_FLAGS -ffast-math)

# Using pybind11 requires using -fvisibility=hidden. See
# https://pybind11.readthedocs.io/en/stable/faq.html
set_source_files_properties(algorithms/python_deconvolution.cpp
                            PROPERTIES COMPILE_FLAGS -fvisibility=hidden)

target_sources(radler PRIVATE ${RADLER_FILES})
target_link_libraries(radler ${RADLER_TARGET_LIBS})
target_include_directories(radler SYSTEM PUBLIC ${RADLER_TARGET_INCLUDE_DIRS})

# Allows including the paths relative to base algorithm, in line with Google
# style
# https://google.github.io/styleguide/cppguide.html#Names_and_Order_of_Includes
target_include_directories(
  radler PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

target_compile_options(radler PRIVATE -O3 -Wall -Wzero-as-null-pointer-constant)

if(NOT COMPILE_AS_EXTERNAL_PROJECT)
  include(GNUInstallDirs)
  install(TARGETS radler ARCHIVE DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR})
else()
  install(TARGETS radler ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()

set(PUBLIC_HEADERS
    component_list.h
    image_set.h
    psf_offset.h
    radler.h
    settings.h
    work_table.h
    work_table_entry.h)

install_headers_with_directory("${PUBLIC_HEADERS}")

if(NOT COMPILE_AS_EXTERNAL_PROJECT)
  configure_file(${PROJECT_SOURCE_DIR}/cmake/config/radler-config.cmake.in
                 ${PROJECT_BINARY_DIR}/CMakeFiles/radler-config.cmake @ONLY)
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/config/radler-config-version.cmake.in
    ${PROJECT_BINARY_DIR}/CMakeFiles/radler-config-version.cmake @ONLY)

  # Install configuration files
  install(FILES ${PROJECT_BINARY_DIR}/CMakeFiles/radler-config.cmake
                ${PROJECT_BINARY_DIR}/CMakeFiles/radler-config-version.cmake
          DESTINATION ${CMAKE_INSTALL_PREFIX}/share/radler)
endif()
