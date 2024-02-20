#!/bin/bash
#SBATCH -J job_catalogs_preparation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=2
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --constraint=HSW24
#SBATCH --mem=115G
#SBATCH --output error_output/job_catalogs_preparation_%I.o
#SBATCH --error error_output/job_catalogs_preparation_%I.e

#export LD_LIBRARY_PATH=/home/mbreton/hdf5-1.8.16/lib:$LD_LIBRARY_PATH



exec_line='mpirun -np 1'
home=/data/home/mbreton/magrathea-pathfinder/bin/ # executable directory
EXEC=$home/observer_velocity
base=/data/mbreton/test_data/simulation/   # simulation directory
sim=boxlen82.03125_n128_lcdmw7 # simulation name
ncoarse=7    

inputconepart='cone_part_fullsky_past_00001'
conenameout='cone_lensing_fullsky_past_00001'
conedir='conedir'

DIR=$base/$sim # The directory must contain the keyword 'boxlen' followed by the box length in Mpc/h

cd $DIR
cd post/slicer_ncoarse/

cd $conenameout
cd $conedir

date
cat > 'raytracer_observer_velocity.txt' << EOF
		  
partdir = $DIR/post/slicer_ncoarse/$inputconepart/# Directory containing DM particles 
ncoarse = $ncoarse # Coarse level of refinement in the simulation
velocity_field_v0 = tsc # interpolation for the velocity field:  "cic" (Cloud-In-Cell), "tsc" (Triangular Shaped Cloud)


EOF


$exec_line $SLURM_NTASKS $EXEC 'raytracer_observer_velocity.txt' > 'raytracer_observer_velocity.log'

ls -lrts

date
