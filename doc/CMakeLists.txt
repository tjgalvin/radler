find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

set(RADLER_SO_PATH ${CMAKE_BINARY_DIR}/python)

set(DOXYGEN_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/doxygen")
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_GENERATED ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_GENERATED} @ONLY)

add_custom_target(
  doxygen_xml
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_GENERATED}
  WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIR}
  COMMENT "Generating documentation with Doxygen")

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/html)

add_custom_target(
  doc
  ${CMAKE_COMMAND}
  -E
  env
  RADLER_SO_PATH=${RADLER_SO_PATH}
  ${SPHINX_EXECUTABLE}
  -b
  html
  ${SPHINX_SOURCE}
  ${SPHINX_BUILD}
  # Tell Breathe where to find the Doxygen output
  -Dbreathe_projects.Radler=${DOXYGEN_OUTPUT_DIR}/xml
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with Sphinx"
  DEPENDS doxygen_xml pyradler)

if(BUILD_DOCSTRINGS)
  set(FILES_TO_DOCUMENT ${CMAKE_CURRENT_SOURCE_DIR}/../cpp/settings.h)

  set(DOCSTRINGS_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/docstrings")
  file(MAKE_DIRECTORY ${DOCSTRINGS_OUTPUT_DIR})

  # Find libclang for cases where pybind11_mkdoc itself fails (like on DAS-6)
  find_library(LIBCLANG_PATH libclang${CMAKE_SHARED_LIBRARY_SUFFIX})
  find_program(PYBIND11_MKDOC "pybind11-mkdoc" REQUIRED)
  set(PYBIND11_MKDOC_BASE_COMMAND "${PYBIND11_MKDOC}" -std=c++17)
  if(LIBCLANG_PATH)
    get_filename_component(LLVM_DIR_PATH ${LIBCLANG_PATH}/../.. ABSOLUTE)
    set(PYBIND11_MKDOC_BASE_COMMAND
        ${CMAKE_COMMAND} -E env LIBCLANG_PATH=${LIBCLANG_PATH}
        LLVM_DIR_PATH=${LLVM_DIR_PATH} ${PYBIND11_MKDOC_BASE_COMMAND})
  endif()

  foreach(RADLER_TARGET_INCLUDE_DIR ${RADLER_TARGET_INCLUDE_DIRS})
    set(PYBIND11_MKDOC_BASE_COMMAND ${PYBIND11_MKDOC_BASE_COMMAND} -I
                                    ${RADLER_TARGET_INCLUDE_DIR})
  endforeach()
  foreach(FILE_TO_DOCUMENT ${FILES_TO_DOCUMENT})
    get_filename_component(DOCSTRINGS_NAME ${FILE_TO_DOCUMENT} NAME_WLE)
    list(
      APPEND
      PYBIND11_MKDOC_COMMAND
      COMMAND
      ${PYBIND11_MKDOC_BASE_COMMAND}
      -o
      ${DOCSTRINGS_NAME}_docstrings.h
      ${FILE_TO_DOCUMENT}
      # Add a comma to the resulting file and make it ISO C++ compliant.
      COMMAND
      sed
      -i
      "'s=\\(__EXPAND(__COUNT(__VA_ARGS__, 7, 6, 5, 4, 3, 2, 1\\)))=\\1, 0))='"
      ${DOCSTRINGS_NAME}_docstrings.h)
  endforeach()

  add_custom_target(
    docstrings
    ${PYBIND11_MKDOC_COMMAND}
    WORKING_DIRECTORY ${DOCSTRINGS_OUTPUT_DIR}
    COMMENT "Creating docstrings from doxygen comments")
endif()
