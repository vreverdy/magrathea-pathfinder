#!/bin/bash
#SBATCH -J job_catalogs_preparation
#SBATCH --nodes=128
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=2
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --constraint=HSW24
#SBATCH --mem=115G
#SBATCH --output error_output/job_catalogs_preparation_%I.o
#SBATCH --error error_output/job_catalogs_preparation_%I.e

#export LD_LIBRARY_PATH=/home/mbreton/hdf5-1.8.16/lib:$LD_LIBRARY_PATH


NCPU=1
exec_line='mpirun -np '$NCPU
home=/data/home/mbreton/magrathea-pathfinder/bin/ # executable directory
EXEC=$home/create_octree
base=/data/mbreton/test_data/simulation/   # simulation directory
sim=boxlen82.03125_n128_lcdmw7 # simulation name
ncoarse=7
ncones=8

inputconegrav='cone_grav_fullsky_pastnfut_00001'
inputconepart='cone_part_fullsky_past_00001'
inputminicone='cone_grav_fullsky_pastnfut_00001'  #CENTRAL BUFFER ZONE FOR CONE NARROW
conenameout='cone_lensing_fullsky_past_00001'
conedir='conedir'

paramfile='/data/mbreton/test_data/simulation/boxlen82.03125_n128_lcdmw7/param_lcdmw7v2.txt'
evolfile='/data/mbreton/test_data/simulation/boxlen82.03125_n128_lcdmw7/ramses_input_lcdmw7v2.dat'

DIR=$base/$sim # The directory must contain the keyword 'boxlen' followed by the box length in Mpc/h

cd $DIR
cd post/slicer_ncoarse/

mkdir $conenameout
cd $conenameout

mkdir $conedir

cd $conedir
date

# PARAMETER FILE
cat > 'raytracer_'$inputconegrav'_prepa.txt' << EOF

# Constants
seed = 42 # Seed for pseudo-random generation
allocation = 402653184 # Pre-allocate the number of cells per MPI task for each octree (can be set to zero but will cause lots of reallocations). Should be close (but larger) to the number of cells in the octree
microcoeff = $((2**${ncoarse}/8)) # Size on inner sphere shared by all process.
mpc = 3.08568E22 # (Mpc value in SI)
rhoch2 = 1.8783467E-26 # (Value of rhoc*h^2)
#
paramfile = $paramfile # Cosmological parameter file.
evolfile = $evolfile # Ramses evolution file. 

# Input files
typefile = 1 # Type of input data. 0 = binary, 1 = hdf5 and 2 = ascii
isfullsky = 1 # Type of lightcone. 1 = fullsky, 0 = narrow
inputtype = cells # "cells" or "particles". Input type to create the octree
coarseonly = 0 # In preparation mode, only accounts for coarse cells.
celldir = $DIR/post/slicer_ncoarse/$inputconegrav/ 	 # Directory containing gravity cells
cellfmt = fof_${inputconegrav}_cube_%05d 		 # Format of binary cubes containing cells	
conedir = $DIR/post/slicer_ncoarse/$conenameout/$conedir/ # Directory in which the binary octree is written
conefmt = fof_${inputconegrav}_cone_%05d # Cone file format 				  
partdir = $DIR/post/slicer_ncoarse/$inputconepart/# Directory containing DM particles, if inputtype = 'particles
minicone = $DIR/post/slicer_ncoarse/$inputminicone/ # Small fullsky light-cone used as buffer zone for narrow light-cones


####### GENERAL #######
ncoarse = $ncoarse # Coarse level of refinement in the simulation
ncones = $ncones # Total number of cones (can be equal to the number of MPI tasks)
buffer = 0.0 # Buffer zones for the square narrow cones, in radians (used only if isfullsky = 0)
firstcone = 0  # First cone for which we compute the octree
lastcone = 7 # Last cone for which we compute the octree
correction = 1 # If equal to 1, allow to correct the octree structure (8 children cells per parent cell)
coarsecorrection = 0 # Allow to correct the density of coarse cells. If rho = 0, the value if computed using neighbours
acorrection = 0 # Allow to correct the value of the scale factor in cells (Specific Full Universe Run)
#########################


EOF

$exec_line $EXEC 'raytracer_'$inputconegrav'_prepa.txt' > 'raytracer_'$inputconegrav'_prepa.log'

ls -lrts

date
