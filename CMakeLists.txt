cmake_minimum_required(VERSION 3.10)

project(nerdss)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -fsanitize=address -fno-omit-frame-pointer")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(default_build_type "Release")
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to default, '${default_build_type}'.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
endif()

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please use a different directory (e.g. 'build').")
endif()

# Get system information
cmake_host_system_information(RESULT NUM_CPUS QUERY NUMBER_OF_PHYSICAL_CORES)
cmake_host_system_information(RESULT NUM_THREADS QUERY NUMBER_OF_LOGICAL_CORES)

# Print system information
message(STATUS "We are on a ${CMAKE_SYSTEM_NAME} system")
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR} with ${NUM_CPUS} physical cores and ${NUM_THREADS} logical cores.")

include_directories(include)

#add_definitions(-D_DEBUG) # comment this line if no need _DEBUG

# ---- Option: CUDA GPU acceleration ----
option(NERDSS_CUDA "Enable CUDA GPU acceleration" OFF)

# Collect sources from specific subdirectories (matching Makefile behavior).
# The debug, mpi, and io_mpi directories require MPI and are excluded from
# the serial and CUDA builds.
set(SERIAL_SOURCE_DIRS
    boundary_conditions classes error math parser reactions
    system_setup trajectory_functions io)

set(SOURCES "")
foreach(dir ${SERIAL_SOURCE_DIRS})
    file(GLOB dir_sources "src/${dir}/*.cpp")
    list(APPEND SOURCES ${dir_sources})
endforeach()

# ---- Serial (CPU-only) target ----
add_executable(nerdss ${SOURCES} EXEs/nerdss.cpp)

# Set up external libraries
find_package(GSL REQUIRED)
target_link_libraries(nerdss GSL::gsl GSL::gslcblas)

# Set up header directories
include_directories(include $(GSL_INCLUDE_DIR))

# ---- CUDA GPU target ----
if(NERDSS_CUDA)
    # Require CMake 3.18+ for better CUDA support; 3.10 is the minimum
    include(CheckLanguage)
    check_language(CUDA)

    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 11)
        set(CMAKE_CUDA_STANDARD_REQUIRED YES)

        message(STATUS "CUDA compiler found: ${CMAKE_CUDA_COMPILER}")
        message(STATUS "CUDA version: ${CMAKE_CUDA_COMPILER_VERSION}")

        # Collect CUDA source files
        file(GLOB CUDA_SOURCES "src/gpu/*.cu")

        # Build the GPU-accelerated executable
        add_executable(nerdss_cuda ${SOURCES} ${CUDA_SOURCES} EXEs/nerdss.cpp)
        target_compile_definitions(nerdss_cuda PRIVATE NERDSS_CUDA)
        target_link_libraries(nerdss_cuda GSL::gsl GSL::gslcblas)
        target_include_directories(nerdss_cuda PRIVATE include ${GSL_INCLUDE_DIR})

        # Set CUDA architectures (common GPUs: 60=Pascal, 70=Volta, 75=Turing, 80=Ampere, 86=Ampere, 89=Ada, 90=Hopper)
        if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
            set_target_properties(nerdss_cuda PROPERTIES CUDA_ARCHITECTURES "60;70;75;80;86")
        endif()

        # Link cuRAND
        find_library(CURAND_LIBRARY curand HINTS ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}/../lib64)
        if(CURAND_LIBRARY)
            target_link_libraries(nerdss_cuda ${CURAND_LIBRARY})
        else()
            target_link_libraries(nerdss_cuda curand)
        endif()

        # Allow .cpp files to include CUDA headers
        set_target_properties(nerdss_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

        message(STATUS "NERDSS CUDA target enabled: nerdss_cuda")
    else()
        message(WARNING "NERDSS_CUDA requested but no CUDA compiler found. GPU target disabled.")
    endif()
endif()