# SPDX-License-Identifier: LGPL-3.0-only

add_executable(multiscale_example multiscale_example.cc)
target_link_libraries(multiscale_example PRIVATE radler)
target_include_directories(multiscale_example
                           PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>)
