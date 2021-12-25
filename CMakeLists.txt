#########
# Setup #
#########

# Specify the minimum version for CMake
cmake_minimum_required(VERSION 2.8)

# Message
message(STATUS "CRP2")


# Find includes in the build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Project's name
project(crp2)

# Set the output folder where your program will be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_STANDARD 14)

set(PATH_CRP2 ${CMAKE_CURRENT_SOURCE_DIR}/src)


###################
# Warnings/Errors #
###################

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Werror=return-type")



###############
# Source Code #
###############

file(GLOB_RECURSE SRC_FILES_CRP2
    ${PATH_CRP2}/*.cpp
    ${PATH_CRP2}/*.c)

set(SRC_FILES ${SRC_FILES_CRP2})

#################
# Library Paths #
#################
# Need to come before target is created.


###########
# Targets #
###########

add_executable(crp2 ${SRC_FILES})

#######################
# Include Directories #
#######################

include_directories(${PATH_CRP2})

#################
# Linker Flags  #
#################

