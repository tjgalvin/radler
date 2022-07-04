# SPDX-License-Identifier: LGPL-3.0-only

# List of tests
set(TEST unittests_radler_pybindings)

# Full path to PYTESTFILES should be specified
set(PYTESTFILES
    ${CMAKE_CURRENT_SOURCE_DIR}/test_component_list.py
    ${CMAKE_CURRENT_SOURCE_DIR}/test_psf.py
    ${CMAKE_CURRENT_SOURCE_DIR}/test_radler.py
    ${CMAKE_CURRENT_SOURCE_DIR}/test_settings.py
    ${CMAKE_CURRENT_SOURCE_DIR}/test_work_table.py)

# Add pybindings tests, can be invoked with "ctest -L pybindings"
add_test(NAME ${TEST} COMMAND python3 -m pytest -v -s --exitfirst
                              --junit-xml=${TEST}.xml ${PYTESTFILES})
set_tests_properties(
  ${TEST} PROPERTIES LABELS "pybindings" ENVIRONMENT
                     "PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH}")
