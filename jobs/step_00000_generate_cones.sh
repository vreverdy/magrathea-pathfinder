#!/bin/bash
#SBATCH -J job_generate_cones
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --constraint=HSW24
#SBATCH --mem=115G
#SBATCH --output error_output/job_generate_cones_%I.o
#SBATCH --error error_output/job_generate_cones_%I.e

#export LD_LIBRARY_PATH=/home/mbreton/hdf5-1.8.16/lib:$LD_LIBRARY_PATH


NCPU=1

exec_line='mpirun -np '$NCPU
home=/data/home/mbreton/magrathea-pathfinder/bin/ # executable directory
EXEC=$home/generate_cones
base=/data/mbreton/test_data/simulation/   # simulation directory
sim=boxlen82.03125_n128_lcdmw7 # simulation name
ncones=8

inputconegrav='cone_grav_fullsky_pastnfut_00001'
conenameout='cone_lensing_fullsky_past_00001'
conedir='conedir'

DIR=$base/$sim # The directory must contain the keyword 'boxlen' followed by the box length in Mpc/h

cd $DIR
cd post/slicer_ncoarse/

mkdir $conenameout
cd $conenameout

mkdir $conedir


cd $conedir
date
# PARAMETER FILE
cat > 'raytracer_'$inputconegrav'_generate_cones.txt' << EOF

# Input files
typefile = 1 # Type of input data. 0 = binary, 1 = hdf5 and 2 = ascii
isfullsky = 1 # Type of lightcone. 1 = fullsky, 0 = narrow
celldir = $DIR/post/slicer_ncoarse/$inputconegrav/ 	 # Directory containing gravity cells
conedir = $DIR/post/slicer_ncoarse/$conenameout/$conedir/ # Directory in which the binary octree is written


####### GENERAL #######
ncones = $ncones # Total number of cones (can be equal to the number of MPI tasks)
buffer = 0.0 # Buffer zones for the square narrow cones, in radians (used only if isfullsky = 0)
#########################


EOF
#

$exec_line $EXEC 'raytracer_'$inputconegrav'_generate_cones.txt' > 'raytracer_'$inputconegrav'_generate_cones.log'

ls -lrts

date
