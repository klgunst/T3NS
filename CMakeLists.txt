cmake_minimum_required(VERSION 3.9.0)
project(T3NS 
    LANGUAGES C 
    DESCRIPTION "An implementation of the Three-Legged Tree Tensor Network algorithm."
    VERSION 0.1.0)
set(T3NS_MAIL "Klaas.Gunst@UGent.be")
string(TIMESTAMP THISTIME "Compiled %b %d %Y %H:%M:%S")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

option(DEBUG 		"Debug symbols used" 			  OFF)
option(DEBUGFLAG        "Compile with DDEBUG flag" 		  OFF)
option(MKL 		"Compile using MKL" 			  OFF)
option(ENABLE_XHOST     "Enable processor-specific optimizations" ON)
option(DAVID_INFO     	"Print intermediate results for the Davidson algorithm" OFF)
option(BUILD_TESTING 	"Compile the tests" 			  ON)
option(PERFORMANCETEST  "Compile the performance tests" 	  OFF)
option(BUILD_DOXYGEN    "Use Doxygen to create a HTML/PDF manual" OFF)
set(MAX_SYMMETRIES "5" CACHE STRING "")

# Find Lapack or MKL
if(MKL)
find_package(MKL)
if(MKL_FOUND)
    message("MKL found")
    include_directories(${MKL_INCLUDE_DIRS})
    set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
    add_definitions(-DT3NS_MKL)
else()
    message("MKL not found")
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    include_directories(${LAPACKE_INCLUDE_DIRS})
    set(LAPACK_LIBRARIES "cblas ${LAPACKE_LIBRARIES}")
endif()
else()
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    include_directories(${LAPACKE_INCLUDE_DIRS})
    set(LAPACK_LIBRARIES "cblas ${LAPACKE_LIBRARIES}")
endif()

# Find OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# Find HDF5
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

# show all warnings
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wremarks -Wchecks -w3 -wd2547 -wd10382 -wd11074 -wd11076 -wd279")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wpedantic")
endif()

# for debug
if(DEBUG)
    if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    else()
	#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer")
    endif()
endif(DEBUG)

add_definitions(-DMAX_SYMMETRIES=${MAX_SYMMETRIES})
# set DDEBUG flag
if(DEBUGFLAG)
    add_definitions(-DDEBUG)
endif(DEBUGFLAG)
if(DAVID_INFO)
    add_definitions(-DDAVID_INFO)
endif(DAVID_INFO)

if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()
message("-- The CMake build type is ${CMAKE_BUILD_TYPE}")

# for aggressive optimization or not
if(ENABLE_XHOST)
    if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHOST")
    else()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    endif()
endif(ENABLE_XHOST)

# Cmake includes
include(GNUInstallDirs)

# Include files
configure_file(
    "${PROJECT_SOURCE_DIR}/include/T3NSConfig.h.in"
    "${PROJECT_BINARY_DIR}/T3NSConfig.h")
include_directories(include ${CMAKE_BINARY_DIR})
add_subdirectory(src)

if(BUILD_DOXYGEN)
    find_package(Doxygen)
    if (NOT DOXYGEN_FOUND)
	message (FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly.")
    endif(NOT DOXYGEN_FOUND)
    configure_file (${PROJECT_SOURCE_DIR}/src/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
    add_custom_target (doc COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
endif(BUILD_DOXYGEN)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
