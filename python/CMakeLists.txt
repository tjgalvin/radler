# SPDX-License-Identifier: LGPL-3.0-only
project(pyradler)

# Create the binding library pyradler is a temporary alias for the radler target
pybind11_add_module(pyradler pyradler.cc pysettings.cc pywork_table.cc
                    pycomponent_list.cc pywrappers.cc)

target_include_directories(${PROJECT_NAME}
                           PRIVATE "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/cpp>")

# Required to hide symols in target, see
# https://pybind11.readthedocs.io/en/stable/faq.html#someclass-declared-with-greater-visibility-than-the-type-of-its-field-someclass-member-wattributes
target_compile_options(
  ${PROJECT_NAME} PRIVATE -O3 -Wall -Wzero-as-null-pointer-constant
                          -fvisibility=hidden)
target_link_libraries(${PROJECT_NAME} PRIVATE radler)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME radler)

# If the PYTHON_LIBRARY_DIR is not specified, install it in proper lib dir
if(NOT DEFINED ${PYTHON_LIBRARY_DIR})
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c
            "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_DIST_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(PYTHON_DIST_PATH
     MATCHES
     "(\\/lib.*\\/)(python${PYTHON_VERSION_MAJOR}\\.${PYTHON_VERSION_MINOR}\\/.*)"
  )

    set(PYTHON_SITE_PACKAGES_DIR ${CMAKE_MATCH_2})

    if(NOT COMPILE_AS_EXTERNAL_PROJECT)
      include(GNUInstallDirs)
      set(PYTHON_LIBRARY_DIR
          "${CMAKE_INSTALL_FULL_LIBDIR}/${PYTHON_SITE_PACKAGES_DIR}")
    else()
      set(PYTHON_LIBRARY_DIR
          "${CMAKE_INSTALL_PREFIX}/lib/${PYTHON_SITE_PACKAGES_DIR}")
    endif()
  else()
    message(
      FATAL_ERROR "Failed to parse PYTHON_DIST_PATH='${PYTHON_DIST_PATH}'")
  endif()
endif()

# Install pyradler in PYTHON_LIBRARY_DIR directory
install(
  TARGETS pyradler
  COMPONENT python
  LIBRARY DESTINATION ${PYTHON_LIBRARY_DIR})
