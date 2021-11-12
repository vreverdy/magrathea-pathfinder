# Makefile for compiling the MAGRATHEA code on Linux systems

#-- Compiler --#

CC = mpic++

#-- Flags --#

CF = -O3 -Wall -Wextra -Wno-unused-parameter  -Wuninitialized -Winit-self -Wno-shift-count-overflow -Wno-shift-count-negative -pedantic -std=c++17

#-- Paths --#

# Path to the HDF5, HEALPIX C++ libraries - must be set by user if they are not installed in the default system path

# The HDF5 library is mandatory to compile the code, even if we do not necessarily use it the get the data
HDF5_PATH = /opt/hdf5/1.10.1/openmpi/3.0.0/gcc/7.3.0/
# Need the Healpix library to compute Healpix maps. Not needed for other purposes
HEALPIX_PATH = /data/home/mbreton/Healpix_3.70
# Needed when using the Healpix library
CFITSTIO_PATH = /data/home/mbreton/cfitsio

#-- Options --#

OPTIONS =
# I/O Node group. IOGROUPSIZE = 1 Means no ticket system
OPTIONS += -DIOGROUPSIZE=1
# If memory issues (only used in octree.update(), so that we do not use shrink_to_fit which cause a copy of an already huge array)
OPTIONS += -DMEMORY_SAVING
# EXTENT (In Ramses Units). Octree will have cells between +-EXTENT/2 along each dimension.
# Photons do not propagate outside of the octree.
# EXTENT = 1 is the size of the simulation box, might need to put a higher value for narrow cones. 
# Must be equal to 2^n with n an integer. 
# As long as the octree volume is larger than that of the data, the value for EXTENT does not matter
OPTIONS += -DEXTENT=4
# Interpolation order for ray-tracing. NGP = 0, CIC = 1, TSC = 2
OPTIONS += -DORDER=2
# VERBOSE
OPTIONS += -DVERBOSE
# Is GCC version below 7.x 
#OPTIONS += -DGCCBELOW7

# DO NOT TOUCH UNLESS YOU KNOW WHAT YOU ARE DOING
OPTIONS += -DSIZEOFOCTREE=48 # size given by sizeof(*(octree.data())), for gravity.h (standard octree)


#########-----------------   DO NOT TOUCH BELOW   ------------------########

#-- Libraries and Includes --#

LIBRARIES = -lpthread
# HDF5
INCLUDES = -I/$(strip $(HDF5_PATH))/include
LIBRARIES += -L/$(strip $(HDF5_PATH))/lib  -lhdf5 -lz
# Healpix + Cfitsio
CF_HEALPIX = -fopenmp 
LIBRARIES_HEALPIX =  -L/$(strip $(CFITSTIO_PATH))
INCLUDES_HEALPIX = -I/$(strip $(HEALPIX_PATH))/include -I/$(strip $(HEALPIX_PATH))/src/cxx/Healpix_cxx -I/$(strip $(HEALPIX_PATH))/src/cxx/cxxsupport -I/$(strip $(HEALPIX_PATH))/src/C/autotools/
LIBRARIES_HEALPIX += -L/$(strip $(HEALPIX_PATH))/lib -L/$(strip $(HEALPIX_PATH))/src/cxx/optimized_gcc/lib  -lchealpix -lcfitsio -lhealpix_cxx # Check the compiling option!

#-- Finalise --#

CF += $(OPTIONS)

#-- Compile --#

EXEC = bin/create_octree bin/observer_velocity bin/rays bin/catalogues bin/hmaps bin/hmaps_velocityfield

all: $(EXEC)

bin/create_octree: src/create_octree.cpp src/create_octree.h 
	$(CC) src/create_octree.cpp -o bin/create_octree $(CF) $(INCLUDES) $(LIBRARIES)
bin/observer_velocity: src/observer_velocity.cpp src/observer_velocity.h 
	$(CC) src/observer_velocity.cpp -o bin/observer_velocity $(CF) $(INCLUDES) $(LIBRARIES)
bin/rays: src/rays.cpp src/rays.h
	$(CC) src/rays.cpp -o bin/rays $(CF) $(INCLUDES) $(LIBRARIES)
bin/catalogues: src/catalogues.cpp src/catalogues.h  src/lensing.h
	$(CC) src/catalogues.cpp -o bin/catalogues $(CF) $(INCLUDES) $(LIBRARIES)
bin/hmaps: src/hmaps.cpp src/hmaps.h src/lensing.h
	$(CC) src/hmaps.cpp -o bin/hmaps $(CF) $(CF_HEALPIX) $(INCLUDES) $(INCLUDES_HEALPIX) $(LIBRARIES) $(LIBRARIES_HEALPIX)
# -DVELOCITYFIELD enables the use of velocity fields for the octree (more memory needed). Useful to compute redshift with peculiar velocities at any location
bin/hmaps_velocityfield: src/hmaps.cpp src/hmaps.h 
	$(CC) src/hmaps.cpp -o bin/hmaps_velocityfield $(CF) $(CF_HEALPIX) -DVELOCITYFIELD $(INCLUDES) $(INCLUDES_HEALPIX) $(LIBRARIES) $(LIBRARIES_HEALPIX)

clean:
	rm -f $(EXEC) *~ */*~ */*/*~ 
clena:
	$(info You type too fast, but I understood what you mean ^_^)
	rm -f $(EXEC) *~ */*~ */*/*~ 
	

