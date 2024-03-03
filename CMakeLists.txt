# spaceseries: CMakeLists.txt

# Copyright (C) 2024  Johnathan K Burchill

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


PROJECT(spaceseries)

CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

ADD_DEFINITIONS(-DHAVE_GETOPT_H)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Debugging enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -std=gnu99 ")
else (CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Optimization enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=gnu99 ")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)

INCLUDE_DIRECTORIES(include)

# GSL
FIND_PACKAGE(GSL REQUIRED)

# Quasi dipole library
FIND_LIBRARY(QUASIDIPOLE quasidipole)

# requires -lm on linux
FIND_LIBRARY(MATH m)

SET(LIBS ${LIBS} ${MATH} ${QUASIDIPOLE} ${GSL_LIBRARY} cdf )

INCLUDE_DIRECTORIES(${INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} "include")

ADD_LIBRARY(spaceseries spaceseries.c interpolation.c util.c programstate.c)
target_link_libraries(spaceseries ${MATH})

ADD_LIBRARY(binnedstatistics statistics.c)
target_link_libraries(binnedstatistics ${GSL_LIBRARY} ${MATH})

ADD_EXECUTABLE(swarm_space_series main.c data.c parseargs.c)
TARGET_INCLUDE_DIRECTORIES(swarm_space_series PRIVATE)
TARGET_LINK_LIBRARIES(swarm_space_series spaceseries ${LIBS})

ADD_EXECUTABLE(analyze_space_series analyze_space_series.c)
TARGET_INCLUDE_DIRECTORIES(analyze_space_series PRIVATE)
TARGET_LINK_LIBRARIES(analyze_space_series ${MATH} spaceseries binnedstatistics cdf)

install(TARGETS swarm_space_series DESTINATION $ENV{HOME}/bin)
install(TARGETS analyze_space_series DESTINATION $ENV{HOME}/bin)
