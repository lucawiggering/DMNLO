# CMake file for DM@NLO.
#
# Licensed under the terms of the EUPL version 1.1 or later.
# See the LICENCE file for more information.

cmake_minimum_required(VERSION 3.0)
project(dmnlo LANGUAGES C CXX Fortran)

# for debugging
# set(CMAKE_VERBOSE_MAKEFILE on )

# set the install directory to the dmnlo folder
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

# specification of the build type. 'Release' enforces the -O3 optimization
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(ExternalProject)

option(BUILD_LOOPTOOLS "Build LoopTools" ON)
set(BUILD_SLHALIB TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/external/cmake)


# LoopTools-2.16
if(BUILD_LOOPTOOLS)
	ExternalProject_Add(lt
		SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/LoopTools-2.16
		CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FFLAGS=-std=legacy CCFLAGS=-fPIC CPPFLAGS=-fPIC ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/LoopTools-2.16/configure --prefix=<INSTALL_DIR>
		BUILD_COMMAND ${MAKE})
	ExternalProject_Get_Property(lt INSTALL_DIR)
	set(LT_DIR ${INSTALL_DIR})
	include_directories(${LT_DIR}/include)
	set(LOOPTOOLS_LIBRARIES ${LT_DIR}/lib/libooptools.a)
else()
	find_package(LoopTools REQUIRED)
	include_directories(${LOOPTOOLS_INCLUDE_DIRS})
endif()

# SLHALib-2.2
if(BUILD_SLHALIB)
	ExternalProject_Add(slha
		SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/SLHALib-2.2
		CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env FFLAGS=-std=legacy CCFLAGS=-fPIC CPPFLAGS=-fPIC ${CMAKE_CURRENT_SOURCE_DIR}/external/lib/SLHALib-2.2/configure --prefix=<INSTALL_DIR>
		BUILD_COMMAND ${MAKE})
	ExternalProject_Get_Property(slha INSTALL_DIR)
	set(SLHA_DIR ${INSTALL_DIR})
	include_directories(${SLHA_DIR}/include)
	set(SLHA_LIBRARIES ${SLHA_DIR}/lib/libSLHA.a)
else()
#	find_package(SLHALib REQUIRED)
#	include_directories(${SLHALIB_INCLUDE_DIRS})
endif()


add_subdirectory(run_dmnlo)
add_subdirectory(src)


