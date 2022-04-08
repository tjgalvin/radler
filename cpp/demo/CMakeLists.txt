# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(multiscale_example multiscale_example.cc)
target_link_libraries(multiscale_example PRIVATE radler)
target_include_directories(multiscale_example
                           PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/cpp>)
