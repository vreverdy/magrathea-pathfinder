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
EXEC=$home/rays
base=/data/mbreton/test_data/simulation/   # simulation directory
sim=boxlen82.03125_n128_lcdmw7 # simulation name
ncoarse=7
ncones=8

inputconegrav='cone_grav_narrow_pastnfut_00002'
conenameout='cone_lensing_narrow_past_00002'
conedir='conedir'
rays='rays'

paramfile='/data/mbreton/test_data/simulation/boxlen82.03125_n128_lcdmw7/param_lcdmw7v2.txt'
evolfile='/data/mbreton/test_data/simulation/boxlen82.03125_n128_lcdmw7/ramses_input_lcdmw7v2.dat'


DIR=$base/$sim # The directory must contain the keyword 'boxlen' followed by the box length in Mpc/h

cd $DIR
cd post/slicer_ncoarse/

mkdir $conenameout
cd $conenameout

mkdir $rays


cd $rays
date
cat > 'raytracer_'$inputconegrav'_rays.txt' << EOF

# Constants
seed = 42 # Seed for pseudo-random generation
microcoeff = $((2**$ncoarse/8)) # Size on inner sphere shared by all process.
mpc = 3.08568E22 # (Mpc value in SI)
rhoch2 = 1.8783467E-26 # (Value of rhoc*h^2)

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
buffer = 0.0 # Buffer zones for the square narrow cones, in radians (used only if isfullsky = 0)
stop_ray = t # Important, corresponds to the stop criterion of bundles.    Possible cases : redshift, lambda, a, t, r
openingmin = 0.0001 # Bundle opening angle (at the observer location, in radians)
nsteps = 4 # Number of integration steps per cell (not only coarse cells but also refined ones)
#########################


base=catalog_boxlen300_n512_lcdmw7v2_00000_fullsky_single_halos_0.000100 # Put the base of the catalog name, if we use it to compute other catalogues or to compute rays (see ray_targets).
halos = 1 # 1 if the targets are halos, 0 if targets are particles (to read in ray-targets)


###### RAYS ########
ray_targets = random # use 'random' or 'catalogue' (in that case, will look for catalogues at outputdir+base) 
nstat = 128 # Compared to the homogeneous case, save statistics every nstats steps, so ref.size()/nstat lines of statistics. Which gives roughly trajectory.size()*nstat/ref.size() lines of statis for a perturbed case. 
nbundlemin = 4 # Number of photons per bundle
nbundlecnt = 1 # Number of bundle parametrisations. Set to 1 by default 
openingcnt = 1 # Number of iterations on opening angle. Set to 1 by default
statistic = distance # Variable on which the average is computed. Possible cases: distance, distance2, homogeneous, inhomogeneous, invdistance, invdistance2
makestat = 0 # Compute statistics
savemode = 1 # -1, 0 ou 1 respecctively allow to save only statistics, each bundle or each photon.
#### IF ray_targets = 'random' !!
ntrajectories = 8 # Number of bundles per cone (= per MPI task) in random directions
#### IF ray_targets = 'catalogue'
massmin = 20000
massmax = 200000
#########################

outputdir = $DIR/post/slicer_ncoarse/$conenameout/$rays/ # Output directory for rays
outputprefix = raytracing_$sim # (Prefix for output files)


EOF


$exec_line $SLURM_NTASKS $EXEC 'raytracer_'$inputconegrav'_rays.txt' > 'raytracer_'$inputconegrav'_rays.log'

ls -lrts

date