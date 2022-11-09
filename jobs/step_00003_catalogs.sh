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

#module purge
#module load gcc/6.3.0
#module load openmpi/2.1.0-gcc-6.3.0
#module load hdf5/1.10.0-gcc-6.3.0
#export LD_LIBRARY_PATH=/home/mbreton/hdf5-1.8.16/lib:$LD_LIBRARY_PATH


NCPU=1
exec_line='mpirun -np '$NCPU
home=/data/home/mbreton/magrathea-pathfinder/bin/ # executable directory
EXEC=$home/catalogues
base=/data/mbreton/   # simulation directory
sim=boxlen110_n1024_lcdmplanck18_00000 # simulation name
ncoarse=10
ncones=64

inputconegrav='cone_grav_narrow_pastnfut_00001'
inputconepart='cone_part_narrow_past_00001/'
inputminicone='cone_grav_fullsky_pastnfut_00001'  #CENTRAL BUFFER ZONE FOR CONE NARROW
conenameout='cone_lensing_narrow_past_00001'
halodir='/fof_b02000m/' # For haloes
#halodir='' # For particles
conedir='conedir_3930mini'
catalogs='catalogs_3930mini'

paramfile='/data/mbreton/boxlen110_n1024_lcdmplanck18_00000/param_lcdmplanck18.txt'
evolfile='/data/mbreton/boxlen110_n1024_lcdmplanck18_00000/ramses_input_lcdmplanck18.dat'

DIR=$base/$sim # The directory must contain the keyword 'boxlen' followed by the box length in Mpc/h

cd $DIR
cd post/slicer_ncoarse/

mkdir $conenameout
cd $conenameout

mkdir $catalogs
mkdir $catalogs/rejected


cd $catalogs
date
cat > 'raytracer_'$inputconegrav'_catalogs.txt' << EOF
# Constants
seed = 42 # Seed for pseudo-random generation
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
sourcedir = $DIR/post/slicer_ncoarse/$inputconepart/$halodir/ # Directory containing sources (DM particles or halos, depending on 'halos = 0 or 1')

# Simulation parameters

####### GENERAL #######
ncoarse = $ncoarse # Coarse level of refinement in the simulation
ncones = $ncones # Total number of cones (can be equal to the number of MPI tasks)
nsteps = 4 # Number of integration steps per cell (not only coarse cells but also refined ones)
#########################
base=catalog_boxlen300_n512_lcdmw7v2_00000_fullsky_single_halos_0.000100 # Put the base of the catalog name, if we use it to compute other catalogues or to compute rays (see ray_targets).

# Observer peculiar velocity
v0x = 0 
v0y = 0 
v0z = 0

####### RELATIVISTIC CATALOGS #######
use_previous_catalogues = 0 # use '0' if need to iterate on source position, or '1' if we already have computed a catalog (in that case, will look for catalogues at outputdir+base), or '2' to rerun rejected sources
halos = 1 # 1 if the targets are halos, 0 if targets are particles
npart = 30000 # Roughly the maximum number of particle wanted (with halos = 0)
zmin = 0.0   # If halos = 0, minimum particle redshift
zmax = 0.025 # If halos = 0, maximum particle redshift. If equal to 0, compute on the full volume
firstcone = 0  # First cone for which we compute the catalogs
lastcone = 0 # Last cone for which we compute the catalogs
cat_accuracy = 1e-8 # Convergence threshold for the newton method at the source (in radians)
#######################

# Regarding the jacobian matrix
jacobiantype = infinitesimal # Computation of the jacobian matrix : "infinitesimal" or "infinitesimal_born" for the use of an infinitesimal beam or "bundle"
stop_bundle = plane # Important, corresponds to the stop criterion of bundles. Possible cases : redshift, lambda, a, t, r, plane
plane = sachs # On which plane we compute the jacobian matrix with bundles or computation of iterations for catalogues. "normal" (normal to the observer line of sight) or "sachs" (normal to the photon)
openingmin = 0.0001 # Bundle opening angle (at the observer location, in radians)


outputdir = $DIR/post/slicer_ncoarse/$conenameout/$catalogs/ # Output directory for catalogs
outputprefix = catalogue_$sim # (Prefix for catalogs name)

EOF


$exec_line $SLURM_NTASKS $EXEC 'raytracer_'$inputconegrav'_catalogs.txt' > 'raytracer_'$inputconegrav'_catalogs.log'

ls -lrts

date
