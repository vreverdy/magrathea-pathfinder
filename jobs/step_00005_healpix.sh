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

module purge
module load gcc/7.3.0
module load openmpi/3.0.0-gcc-7.3.0 
module load hdf5/1.10.1-gcc-7.3.0

#export LD_LIBRARY_PATH=/home/mbreton/hdf5-1.8.16/lib:$LD_LIBRARY_PATH


NCPU=8
exec_line='mpirun -np '$NCPU
home=/data/home/mbreton/magrathea-pathfinder/bin/ # executable directory
EXEC=$home/hmaps
base=/data/mbreton/test_data/simulation/   # simulation directory
sim=boxlen82.03125_n128_lcdmw7 # simulation name
ncoarse=7
ncones=8

inputconegrav='cone_grav_narrow_pastnfut_00002'
inputminicone='cone_grav_fullsky_pastnfut_00001'  #CENTRAL BUFFER ZONE FOR CONE NARROW
conenameout='cone_lensing_narrow_past_00002'
healpix_outputdir='healpix_maps'
conedir='conedir'

paramfile='/data/mbreton/test_data/simulation/boxlen82.03125_n128_lcdmw7/param_lcdmw7v2.txt'
evolfile='/data/mbreton/test_data/simulation/boxlen82.03125_n128_lcdmw7/ramses_input_lcdmw7v2.dat'

DIR=$base/$sim # The directory must contain the keyword 'boxlen' followed by the box length in Mpc/h

cd $DIR
cd post/slicer_ncoarse/

mkdir $conenameout
cd $conenameout

mkdir $healpix_outputdir


cd $healpix_outputdir
date
cat > 'raytracer_'$inputconegrav'_healpix.txt' << EOF
# Constants
seed = 42 # Seed for pseudo-random generation
mpc = 3.08568E22 # (Mpc value in SI)
rhoch2 = 1.8783467E-26 # (Value of rhoc*h^2)
microcoeff = $((2**$ncoarse/8)) # Size on inner sphere shared by all process.

# Input files
typefile = 1 # Type of input data. 0 = binary, 1 = hdf5 and 2 = ascii
isfullsky = 0 # Type of lightcone. 1 = fullsky, 0 = narrow
paramfile = $paramfile # Cosmological parameter file.
evolfile = $evolfile # Ramses evolution file. 

celldir = $DIR/post/slicer_ncoarse/$inputconegrav/ 	 # Directory containing gravity cells
conedir = $DIR/post/slicer_ncoarse/$conenameout/$conedir/ # Directory in which the binary octree is written
conefmt = fof_${inputconegrav}_cone_%05d # Cone file format 				  
# Simulation parameters

####### GENERAL #######
ncoarse = $ncoarse # Coarse level of refinement in the simulation
ncones = $ncones # Total number of cones (can be equal to the number of MPI tasks)
buffer = 0.0 # Buffer zones for the square narrow cones, in radians (ued only if isfullsky = 0)
nsteps = 4 # Number of integration steps per cell (not only coarse cells but also refined ones)
stop_ray = t # Important, corresponds to the stop criterion of bundles.    Possible cases : redshift, lambda, a, t, r
#########################

###### HEALPIX ########
nside = 256 # total of 12*nside*nside pixels
map_components = lensing, lensing_born, deflection, steps
nb_z_maps = 2 # Number of maps within redshifts
z_stop_min = 0.01 # Minimum redshift (and the only redshift if nb_z_maps = 1)
z_stop_max = 0.025 # Maximum redshift (cannot go beyond the maximum redshift of the lightcone)
#######################

## If you choose 'lensing' in map_components:
beam = bundle # Computation of the jacobian matrix : "infinitesimal" or "bundle"
## if beam = 'bundle' :
stop_bundle = lambda # Important, corresponds to the stop criterion of bundles (if beam = 'bundle'). Possible cases : redshift, lambda, a, t, r, plane
plane = normal # On which plane we compute the jacobian matrix with bundles or computation of iterations for catalogues. "normal" (normal to the observer line of sight) or "sachs" (normal to the photon) 
openingmin = 0.0001 # Bundle opening angle (at the observer location, in radians)

# output
outputdir = $DIR/post/slicer_ncoarse/$conenameout/$healpix_outputdir/
outputprefix = map_$sim # (Prefix for catalogs name)

nbundlemin = 4 # Number of photons per bundle (ununsed)

EOF


$exec_line $SLURM_NTASKS $EXEC 'raytracer_'$inputconegrav'_healpix.txt' > 'raytracer_'$inputconegrav'_healpix.log'

ls -lrts

date
