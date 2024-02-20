<a id="top"></a>

<h3 align="center"> MAGRATHEA-PATHFINDER  </h3>

  <p align="center">
    A code for 3D AMR ray-tracing in simulations
  </p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li> <a href="#about-the-project">About The Project</a> </li>
    <li>
      <a href="#installation">Installation</a>
      <ul>
        <li><a href="#requirements">Requirements</a></li>
      </ul>
      <ul>
        <li><a href="#compilation">Compilation</a></li>
      </ul>
    </li>
    <li> <a href="#running-the-code">Running the code</a> </li>
    <li>
       <a href="#jobs">Jobs</a>
       <ul>
        <li><a href="#step-0---generate-cones">Step 0 - Generate cones</a></li>
        <li><a href="#step-1---preparation">Step 1 - Preparation</a></li>
        <li><a href="#step-2---observer-velocity">Step 2 - Observer velocity</a></li>
        <li><a href="#step-3---rays">Step 3 - Rays</a></li>
        <li><a href="#step-4---catalogues">Step 4 - Catalogues</a></li>
        <li><a href="#step-5---maps">Step 5 - Maps</a></li>
        <li><a href="#other-files">Other files</a></li>
       </ul>
    </li>
    <li>
      <a href="#testing-the-code">Testing the code</a>
    </li>
    <li><a href="#further-documentation">Further documentation</a></li>
    <li><a href="#citations">Citations</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

The goal of the present framework is to propagate photons 
within cosmological simulations to construct observables.
 
Authors: Vincent Reverdy, Michel-Andrès Breton

# About the project

Magrathea-Pathfinder is a high-performance computational framework for relativistic raytracing at cosmological scales.
It is built on top of the Magrathea (Multi-processor Adaptive Grid Refinement Analysis for THEoretical Astrophysics) metalibrary (Reverdy, 2014) available at https://github.com/vreverdy/magrathea

<p align="right">(<a href="#top">back to top</a>)</p>


# Installation

## Requirements 

In any case, you need to link paths to the HDF5 library.


To produce maps, the Healpix library is required, 
as well as the cfitsio library for the outputs.
 

<p align="right">(<a href="#top">back to top</a>)</p>


## Compilation

Clone the repository and compile the code

Compile the code
```
git clone https://github.com/vreverdy/magrathea-pathfinder
cd magrathea-pathfinder/
make
```
This will produce all the executables in `./bin/`

The compiler must be GCC version must be at least 4.8.1 (Does not work for GCC 8.X, but does for all other versions up to GCC 11.X, this will be fixed soon)
If the version of GCC is inferior to 7.x.x, it must be indicated in the Makefile
(uncomment the `OPTIONS += -DGCCBELOW7` line)

<p align="right">(<a href="#top">back to top</a>)</p>


# Running the code


The command line

`
	mpirun -np NTASKS ./bin/executable parameter.txt
`

runs the code in parallel using MPI. `NTASKS` is the number of tasks used (can be equal to 1).
When enough cpus are available, the code uses C++11 std::thread multithreading for each task.
The executable `./bin/executable` depends on the use, and so is the parameter file. 

See Section [Jobs](#jobs) for more informations on the executables.

<p align="right">(<a href="#top">back to top</a>)</p>


# Jobs


Instead of running the code with command lines, we recommand to use the (self-explanatory) shell scripts at your disposal in ./jobs/, which are 

- step_00000_generate_cones.sh
- step_00001_preparation.sh
- step_00002_observer_velocity.sh  
- step_00003_rays.sh  
- step_00004_catalogs.sh  
- step_00005_healpix.sh

To run the jobs, just type (for example for the first job)

`sh step_00001_preparation.sh`
or
`./step_00001_preparation.sh`


This will produce a parameter file and a log file in the output directory.

<p align="right">(<a href="#top">back to top</a>)</p>


## Step 0 - Generate cone directions

Magrathea-Pathfinder uses a cone-shaped domain decomposition. The first step is to compute (on a single CPU) the cone directions and aperture.

The associated executable and job are respectively `./bin/generate_cones` and `step_00000_generate_cones.sh`, which produces a text file containing for each cone: `vertex [x, y, z], base [x,y,z], angle (in radians)` 

<p align="right">(<a href="#top">back to top</a>)</p>


## Step 1 - Preparation

Produce a N-dimensional gravity octree with Adaptive-Mesh Refinement (AMR) from simulation data 
in Magrathea format. For most uses the octree is 3-dimensional.

The associated executable and job are respectively `./bin/create_octree` and `step_00001_preparation.sh`, which will produce a binary file per cone.


<p align="right">(<a href="#top">back to top</a>)</p>

### Domain decomposition


Magrathea-Pathfinder divides the simulation in cone-shaped sub-domains. These cones share a common spherical region
at their vertex (which coincides with the observer's position). These cones contain some overlaping regions, 
so that photons do not leave the cone due to gravitational lensing. In Magrathea-Pathfinder, photons are propagated in the cones that are assigned to them. 
Because these cones cover different areas, they do not communicate (except for some specific applications as in Section [Step 3 - Rays](#step-3---rays)).

If the number of MPI tasks is the same as the number of cones, then each task takes care of one cone.
If the number of MPI tasks is smaller than the number of cones, then individual tasks will produce multiple cones.

<p align="right">(<a href="#top">back to top</a>)</p>

### Information needed


To produce a gravity octree which will be used during ray-tracing, for each cell we need to have access to:

- 3D position
- The scale factor
- The gravitational potential
- The force (cartesian derivatives of the gravitational potential)
- The density (optional)
- The gravitational potential time derivative  (optional)

<p align="right">(<a href="#top">back to top</a>)</p>


### Input format


The input format can be either binary, HDF5 or ascii.

For binary and HDF5, the format is specific to the RayGal simulations (Rasera et al. 2022, see also Section [Testing the code](#testing-the-code)) 
which come from the codes RAMSES (Teyssier 2002) and pFoF (Roy et al. 2014).

ASCII: 
The input data are ascii files which must end with suffix '.dat'
These are 10 column files with 

level (AMR level), x_pos, y_pos, z_pos, density, phi (grav. potential), dphidx, dphidy, dphidz, a (scale factor)
with all the informations in Ramses Units (comoving, normalised so that the box length is unity) 

To produce a narrow light-cone with ascii inputs we also need a 'info_narrow_cone.txt' file which contains 4 rows.
Example of 'info_narrow_cone.txt' file:

``` ini
# info_narrow_cone.txt
phi = 20
theta = 20
thetay = 10
thetaz = 10
```

(phi, theta) gives the orientation of the cone w.r.t the x-axis (in degrees).
In this example the cone is rotated by 20 degrees along the phi and theta angles in
spherical coordinates. (thetay, thetaz) gives the semi-aperture of the lightcone at the observer. 
The cone has a square basis and in this example, a 10x10 deg^2 aperture.

<p align="right">(<a href="#top">back to top</a>)</p>

### Input type


To produce a gravity octree, Magrathea-Pathfinder needs as inputs either a gravity grid or a distribution of particles (specific to HDF5 format) 
For the latter informations are interpolated as some cell locations with refined mesh if we have more than 8 particles per cell.
In any case, the input data needs to contain all the information in 4.1.

<p align="right">(<a href="#top">back to top</a>)</p>

### Outputs


Magrathea-Pathfinder produces binary files which contain the gravity information on an AMR grid. Each file is associated to one cone.
The output directory is given by 'conedir' in the parameter file.

<p align="right">(<a href="#top">back to top</a>)</p>

## Step 2 - Observer velocity


Mono-cpu job which compute the velocity field at the observer. This should run with a single task.

The associated executable and job are respectively `./bin/observer_velocity` and `step_00002_observer_velocity.sh`  

Magrathea-Pathfinder reads the particle file at z = 0, compute the velocity field in the eight cells around the observer, and average over them. Since the result depends on the size of the neighbouring cells, we choose to 
compute the velocity field using 5 different levels, starting from the coarse level of the simulation to
finer levels. To perform the interpolation, the user can choose to use the CIC or TSC algorithms.

The output file is a *_observer_velocity.log file in the same directory as the cones.
The value of the velocity field at the observer is written for each level, in Ramses Units and SI.

<p align="right">(<a href="#top">back to top</a>)</p>

## Step 3 - Rays


This jobs propagates photons in the simulation, and compute various statistics.

The associated executable and job are respectively ./bin/rays and step_00002_rays.sh  

This jobs allows the user to propagate photons in the simulation, either in random directions 
or toward sources that have already been ray-traced.

We must run the code with as many MPI tasks as there are cones.

<p align="right">(<a href="#top">back to top</a>)</p>

### Random directions


It is possible to launch an arbitrary number of photons in a beam around a central one. It is also possible 
to choose the number of beams per cone. The user can choose to output all the individual trajectories,
or statistics regarding functions or the angular diameter distance as a function of redshift in each cone or averaged over the full lightcone.

<p align="right">(<a href="#top">back to top</a>)</p>

### Catalogue directions


Alternatively, it is possible to output the full photon trajectories when launched in the observed directions of sources.
To do so, one needs to produce source catalogues (see Section [Step 4 - Catalogues](#step-4---catalogues)). 

<p align="right">(<a href="#top">back to top</a>)</p>

### Outputs


The outputs are ascii files in the 'rays' directory.

The format for the full trajectory of individual rays is:

step, a, t, x, y, z, k^0, k^1, k^2, k^3, level, a*, rho, phi, dphi/dx, dphi/dy, dphi/dz, dphi/dl, dphi/dl*, redshift, (ds/dl)^2, error, angular diameter distance, isw, isw*, comoving distance, l,    dphi/dt, s

where the wavevector k^\mu = dx^\mu/dl.

All the informations are in SI, except for the density rho which is normalised by the mean matter density and the position (x,y,z, comoving distance) which are in Ramses Units.
t is the conformal time in seconds.
The column 'angular diameter distance' is non zero only if we use a light bundle.
The column 'l' is the affine parameter
The last column 's' is the total distance travelled.
The columns which contain '*' yields a different (less accurate) method to compute a quantity that we also provide differently.

The format of statistics file is:

redshift, distance (FLRW), distance, standard deviation

<p align="right">(<a href="#top">back to top</a>)</p>


## Step 4 - Catalogues


This job produces relativistic source catalogues using ray-tracing.

The associated executable and job are respectively ./bin/catalogues and step_00003_catalogues.sh  

<p align="right">(<a href="#top">back to top</a>)</p>

### Catalogue production


Produce catalogue from sources. This option is enabled when we have `use_previous_catalogues = 0`

The input source files must be in HDF5 RayGal format.

Magrathea-Pathfinder uses an iterative method to find the null geodesic that connects the source to the observer, 
up to a precision given by 'cat_accuracy'. 

The output catalogues for each cone are written as ascii files in the 'catalogue' directory. The columns are:

Id, beta1, beta2, theta1, theta2, error, z0, z1, z2, z3, z4, z5, a11, a12, a21, a22, npart

beta and theta are the comoving and observed angles respectively. For the angles, the indexes '1' and '2' refer to
phi and theta, the angles in spherical coordinates. The angles are given in radians.

error gives the angle at the observer between the photon at the source and the source.

z0 = FLRW redshift \
z1 = z0 + potential \
z2 = z1 + RSD \
z3 = z2 + T. Doppler \
z4 = z3 + ISW/RS \
z5 = exact calculation

a11, a12, a21 and a22 are the components of the lensing distortion matrix.

npart gives the number of particles within the halo (equal to 1 if the target are DM particles).

Additionally, there are `*_err.txt` and `*_rejected.txt` files.

The former are sources which cannot be accessed by ray-tracing (generally because it is outside, or at the limits of our numerical light-cone). The latter are sources which did not converge during the root-finding algorithm
(see Section [Re-run rejected sources](#re-run-rejected-sources) in order to re-run them in order to obtain a full sample).

<p align="right">(<a href="#top">back to top</a>)</p>

### Re-run the lensing matrix


If one wants a different definition for the lensing distortion matrix, it is possible to
read the catalogues in order to launch the photons directly in the observed direction (hence we do not need to
re-run the root-finding algorithm the connect the observer to sources).

This option is enabled when we have `use_previous_catalogues = 1`, together with the token `base` which
gives the name format of the catalogues.

<p align="right">(<a href="#top">back to top</a>)</p>


### Re-run rejected sources


Sometimes, our root-finding method does not converge quickly enough. In this case the sources are put in `*_rejected.txt` files

By using `use_previous_catalogues = 2`, these sources are re-run, starting from the last tentative observed angle. Photons are launch in a fine grid around this direction and we take the pixel for which the distance to the source is the closest at the same comoving distance.
This iterative method is much slower than the previous ones.

<p align="right">(<a href="#top">back to top</a>)</p>


## Step 5 - Maps


This job propagates photons in the light-cone in the direction of Healpix pixel and produces Healpix maps

The associated executable and job are respectively `./bin/hmaps` and `step_00004_healpix.sh`  

The photons can be stopped at some given surface, which roughly coincides with the redshifts requested (with `z_stop_min`, `z_stop_max`, and `nb_z_maps` ) in inputs.

<p align="right">(<a href="#top">back to top</a>)</p>

### Velocity field


Alternatively, it is possible to use `./bin/hmaps_velocityfield` and `step_00004_healpix_velocityfield.sh`  

In this case, the octree is modified to contain additional slots for the velocity field. On the fly, Magrathea-Pathfinder computes the velocity field for each cell of the octree from particles.

The main difference with respect to the 'normal' version is that here the redshift contains the Doppler terms.

WARNING: computing the velocity field on the fly with ALL the DM particles in the vicinity of the redshifts of interest will considerably slow down the code.

<p align="right">(<a href="#top">back to top</a>)</p>

### Outputs


For each pixel at these surfaces we can output different kind of maps with the following keywords:

- lensing	(Output the lensing statistics for infinitesimal or bundle methods)
- lensing_born  (same but use the infinitesimal method around a FLRW trajectory)
- dr		(relative difference on the comoving distance)
- dl		(relative difference on the affine parameter)
- dt		(relative difference on the conformal time)
- da		(relative difference on the scale factor)
- dz		(relative difference on the redshift)
- ds		(relative difference on the distance travelled)
- isw		(ISW/RS effect)
- dens		(local density)
- dens_max	(maximum overdensity of the trajectory)
- steps		(number of integration steps)
- phi		(gravitational potential)

In it possible to give a list of arguments containing the above keywords, the output maps will contain some index
which refer to the position of the argument in the list. 

For example, when putting `map_components = lensing, lensing_born`, in the .log file we have:

Index: 0, component: kappa ('lensing') \
Index: 1, component: gamma1 ('lensing') \
Index: 2, component: gamma2 ('lensing') \
Index: 3, component: inverse magnification ('lensing') \
Index: 4, component: kappa ('lensing_born') WARNING: with lensing_born we compute the jacobian matrix with a single photon (no bundle method available) \
Index: 5, component: gamma1 ('lensing_born') \
Index: 6, component: gamma2 ('lensing_born') \
Index: 7, component: inverse magnification ('lensing_born') 

This gives the relation between index in the output names and nature of the map.

<p align="right">(<a href="#top">back to top</a>)</p>

## Other files


For all the steps, we need two additional files

### Evolution file

`data/ramses_input_xxxxx.dat`

Beginning a is 0.001 roughly, final a is 1. exactly and  number of elements is  1000 roughly
a, H/H0, tau_supercomoving*H0, tlookback*H0, tproper*H0
a = 1 at z = 0 (today), t and tau = 0 at z = 0,
a is the expansion factor, h the Hubble constant defined as hexp=1/a*da/dtau
 tau the supercomoving conformal time (integral of dt/a^2), t the look-back time, tproper is integral of dt/a

<p align="right">(<a href="#top">back to top</a>)</p>


### Cosmological parameters


`data/param_xxxxx.txt`

Ramses parameter file read by Magrathea-Pathfinder to get the hubble parameter.

Example of parameter file:
``` ini
h=0.72000000                
Omega_m=0.25733000        
Omega_b=0.043557099
omega_r=0.000080763524
n_s=0.96300000
sigma_8=0.80100775  
w=-1.0000000
```
<p align="right">(<a href="#top">back to top</a>)</p>


# Testing the code


Download the test data using the command line

`wget http://cosmo.obspm.fr//efiler2_raygal/misc/test_data.tar.gz`

In `test_data/simulation/` you will find the gravity and particle data.

You can directly run the jobs in `magrathea_pathfinder/jobs`,
only changing the paths to the one of the data.

You can compare your results to the ones in `test_data/expected_results/`

<p align="right">(<a href="#top">back to top</a>)</p>


# Further documentation


An exhaustive documentation generated with Doxygen
for every function is available in `doc/doxygen/latex/refman.pdf`

For more information on the underlying Magrathea metalibrary, please see:
https://github.com/vreverdy/magrathea

<p align="right">(<a href="#top">back to top</a>)</p>


# Citations

Vincent Reverdy
Michel-Andres Breton

If you use Magrathea-Pathfinder for your research, please cite:
- Reverdy (2014), PhD thesis (https://hal.archives-ouvertes.fr/tel-02095297/document)
- Breton & Reverdy (2021), arXiv:2111.08744

<p align="right">(<a href="#top">back to top</a>)</p>

# Contact

For any question, recommandation or contribution regarding Magrathea-Pathfinder,
do not hesitate to send an email to 

michel-andres.breton@obspm.fr \
vince.rev@gmail.com

<p align="right">(<a href="#top">back to top</a>)</p>
