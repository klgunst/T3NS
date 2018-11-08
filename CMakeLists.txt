cmake_minimum_required(VERSION 3.9.0)
project(T3NS 
    LANGUAGES C 
    DESCRIPTION "An implementation of the Three-legged Tree Tensor Network algorithm."
    VERSION 0.1.0)
set(T3NS_MAIL "Klaas.Gunst@UGent.be")
string(TIMESTAMP THISTIME "Compiled %b %d %Y %H:%M:%S")

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")

option(DEBUG 		"Debug symbols used" 			  OFF)
option(DEBUGFLAG        "Compile with DDEBUG flag" 		  OFF)
option(LOWOPT 		"Optimization level of compilation O/O3"  OFF)
option(MKL 		"Compile using MKL" 			  OFF)

option(BUILD_TESTING 	"Compile the tests" 			  ON)
option(PERFORMANCETEST  "Compile the performance tests" 	  OFF)

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
    set(LAPACK_LIBRARIES "${LAPACKE_LIBRARIES}")
endif()
else()
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    include_directories(${LAPACKE_INCLUDE_DIRS})
    set(LAPACK_LIBRARIES "${LAPACKE_LIBRARIES}")
endif()

# Find OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# Find GSL
find_package(GSL REQUIRED)

# Find HDF5
find_package(HDF5 REQUIRED)

# show all warnings
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wremarks -Wchecks -w3 -wd2547 -wd10382 -wd11074 -wd11076")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif()

# for debug
if(DEBUG)
    if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
    else()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
	#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=address")
    endif()
endif(DEBUG)

# set DDEBUG flag
if(DEBUGFLAG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDEBUG")
endif(DEBUGFLAG)

# for aggressive optimization or not
if(LOWOPT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
else()
    if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHOST -O3")
    else()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -O3")
    endif()
endif(LOWOPT)

# Cmake includes
include(GNUInstallDirs)

# Include files
configure_file(
    "${PROJECT_SOURCE_DIR}/include/T3NSConfig.h.in"
    "${PROJECT_BINARY_DIR}/T3NSConfig.h")
include_directories(include ${CMAKE_BINARY_DIR})
add_subdirectory(src)

include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
