cmake_minimum_required(VERSION 3.10)
project(MAGRATHEA-PATHFINDER)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler & Flags
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-unused-parameter -Wuninitialized -Winit-self -Wno-shift-count-overflow -Wno-shift-count-negative -pedantic -march=native")

# (Optional) Paths - Set to paths if installed manually
set(HDF5_PATH "/obs/mbreton/hdf5_src")
set(HEALPIX_PATH "/data/home/mbreton/Healpix_3.81")
set(CFITSTIO_PATH "/data/home/mbreton/cfitsio")

# Options
add_definitions(
    -DIOGROUPSIZE=1
    -DEXTENT=4
    -DORDER=2
    -DVERBOSE
    -DSIZEOFOCTREE=48
)

# Find HDF5 Library
find_path(HDF5_INCLUDE_DIR hdf5.h 
    PATHS 
        /usr/include 
        /usr/local/include 
        ${HDF5_PATH}/include
)
find_library(HDF5_LIBRARY NAMES hdf5 
    PATHS
        /usr/include 
        /usr/local/include  
        ${HDF5_PATH}/lib
)
message(STATUS "HDF5_INCLUDE_DIR: ${HDF5_INCLUDE_DIR}")
message(STATUS "HDF5_LIBRARY: ${HDF5_LIBRARY}")

# Check if HDF5 is found
if(NOT HDF5_INCLUDE_DIR OR NOT HDF5_LIBRARY)
    message(FATAL_ERROR "HDF5 library and/or include files not found at ${HDF5_PATH}. Please ensure HDF5 is installed and the path is correct.")
endif()

# Find Healpix and CFITSTIO Libraries
find_path(CHEALPIX_INCLUDE_DIR chealpix.h 
    PATHS 
        /usr/include 
        /usr/local/include 
        ${HEALPIX_PATH}/include
)
find_path(HEALPIX_CXX_INCLUDE_DIR healpix_map.h PATHS ${CHEALPIX_INCLUDE_DIR}/healpix_cxx)
find_library(CHEALPIX_LIBRARY NAMES chealpix PATHS ${HEALPIX_PATH}/lib)
find_library(HEALPIX_CXX_LIBRARY NAMES healpix_cxx PATHS ${HEALPIX_PATH}/lib)
find_library(CFITSTIO_LIBRARY NAMES cfitsio 
    PATHS   
    /usr/include 
    /usr/local/include
    ${CFITSTIO_PATH}/lib
)


message(STATUS "CHEALPIX_INCLUDE_DIR: ${CHEALPIX_INCLUDE_DIR}")
message(STATUS "CHEALPIX_LIBRARY: ${CHEALPIX_LIBRARY}")
message(STATUS "HEALPIX_CXX_INCLUDE_DIR: ${HEALPIX_CXX_INCLUDE_DIR}")
message(STATUS "HEALPIX_CXX_LIBRARY: ${HEALPIX_CXX_LIBRARY}")
message(STATUS "CFITSTIO_LIBRARY: ${CFITSTIO_LIBRARY}")


# Include Directories
include_directories(${HDF5_INCLUDE_DIR})

# Set librairies
set(LIBS pthread z tbb ${HDF5_LIBRARY})

# Executables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(generate_cones src/generate_cones.cpp)
target_link_libraries(generate_cones ${LIBS})

add_executable(create_octree src/create_octree.cpp)
target_link_libraries(create_octree ${LIBS})

add_executable(observer_velocity src/observer_velocity.cpp)
target_link_libraries(observer_velocity ${LIBS})

add_executable(rays src/rays.cpp)
target_link_libraries(rays ${LIBS})

add_executable(catalogues src/catalogues.cpp)
target_link_libraries(catalogues ${LIBS})

# Check if both Healpix and CFITSTIO are found
if(HEALPIX_CXX_INCLUDE_DIR AND HEALPIX_CXX_LIBRARY AND CHEALPIX_INCLUDE_DIR AND CHEALPIX_LIBRARY AND CFITSTIO_LIBRARY)

    message(STATUS "Healpix and CFITSIO libraries found. Proceeding with hmaps compilation.")

    # Enable OpenMP
    find_package(OpenMP REQUIRED)

    # Include Directories    
    include_directories(${CHEALPIX_INCLUDE_DIR})
    include_directories(${HEALPIX_CXX_INCLUDE_DIR})

    list(APPEND LIBS ${CFITSTIO_LIBRARY} ${CHEALPIX_LIBRARY} ${HEALPIX_CXX_LIBRARY})

    add_executable(hmaps src/hmaps.cpp)
    target_link_libraries(hmaps ${LIBS} OpenMP::OpenMP_CXX)
    
    add_executable(hmaps_velocityfield src/hmaps.cpp)
    target_link_libraries(hmaps_velocityfield ${LIBS} OpenMP::OpenMP_CXX)
    target_compile_definitions(hmaps_velocityfield PRIVATE -DVELOCITYFIELD)
    
    else()
    message(WARNING "Healpix or CFITSIO library not found. Skipping hmaps compilation.")
endif()    




# Clean target
add_custom_target(clean_all
    COMMAND rm -rf CMakeFiles CMakeCache.txt Makefile cmake_install.cmake bin/*
    COMMENT "Cleaning up build artifacts"
)
