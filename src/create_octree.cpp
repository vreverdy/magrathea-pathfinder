/* ******************************* CREATE_OCTREE ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Create_octree
// DESCRIPTION :    Main function of the raytracer
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           raytracer.cpp
/// \brief          Main function of the raytracer
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
#include <ctime>
// Include C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>
#include <random>
#include <string>
#include <vector>
#include <deque>
#include <array>
#include <tuple>
#include <atomic>
#include <mutex>
// Include libs
#include <mpi.h>
// Include project
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/hypersphere.h"
#include "magrathea/hypercube.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/timer.h"
#include "cone.h"
#include "utility.h"
#include "input.h"
#include "output.h"
#include "integrator.h"
#include "miscellaneous.h"
#include "create_octree.h"
// Octree
#ifdef VELOCITYFIELD
#include "gravity2.h" // with velocity slots
#else
#include "gravity.h"
#endif
// Include C++ HDF5 project
#include "TReadHDF5.h"


// Misc
using namespace magrathea;
// -------------------------------------------------------------------------- //



// ---------------------------------- MAIN ---------------------------------- //
// Main function
/// \brief          Main function.
/// \details        Main function of the program. It takes a namelist as a 
///                 parameter and run the raytracing using MPI parallelization.
///                 To compile it, use make
///                 To run it, use:
///                 mpirun -np NTASKS ./raytracer raytracer.txt
/// \param[in]      argc Number of arguments.
/// \param[in]      argv List of arguments.
/// \return         Zero on success, error code otherwise.
int main(int argc, char* argv[])
{
    // Constants

    using integer = int;
    using uint = unsigned int;
    using real = double;
    using floating = float;
    using point = std::array<real, 3>;
    using position = std::ratio<0>;
    using extent = std::ratio<EXTENT>;
    using indexing = __uint128_t;
    static constexpr uint zero = 0;
    static constexpr uint one = 1;
    static constexpr real rone = 1.; 
    static constexpr uint two = 2;
    static constexpr uint dimension = 3;
    static constexpr uint nreference = 5; // Used to set homogeneous octree
    static constexpr real rposition = static_cast<real>(position::num)/static_cast<real>(position::den);
    static constexpr point center({{rposition, rposition, rposition}});
    static constexpr real diameter = static_cast<real>(extent::num)/static_cast<real>(extent::den);
    static const std::string namelist = argc > 1 ? std::string(argv[1]) : std::string("raytracer.txt"); 

    // Parameters
    std::map<std::string, std::string> parameter;
    integer nthreads = std::thread::hardware_concurrency();
    integer ntasks = nthreads*zero;
    integer rank = zero;

    // Message passing interface
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Read parameter file
    Miscellaneous::TicketizeFunction(rank, ntasks, [=, &parameter]{parameter = Input::parse(namelist);});
    // Convert strings and put it in struct
    Create_octree::ReadParamFile(parameters, parameter);
    // Initialization
    FileList conefile(parameters.conefmt, zero, parameters.ncones, zero, parameters.conedir);
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<uint, dimension>, std::string, dimension, position, extent> filetree;
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>, Gravity<floating, dimension>, dimension, position, extent> octree;
    HyperSphere<dimension, point> sphere(center, diameter/two);
    HyperSphere<dimension, point> microsphere(center, rone/parameters.microcoeff);
    std::vector<Cone<point> > cone(parameters.ncones);
    std::vector<Cone<point> > coneIfRot(parameters.ncones);
    std::array<std::vector<real>, two+two> cosmology;
    Evolution<Photon<real, dimension> > reference;
    real h = zero;
    real omegam = zero;
    real lboxmpch = zero;
    real amin = zero;
    real thetay(0), thetaz(0);
    std::array< std::array< double, 3 >, 3 > rotm1 = {{zero}};
    const point vobs0 = {0,0,0}; // No peculiar velocity for homogeneous quantities
 
    if(rank == 0) 
	std::cout<<"#### MAGRATHEA_PATHFINDER "<<std::endl; 
    // Generate cones
    Miscellaneous::TicketizeFunction(rank, ntasks, [=, &cone, &coneIfRot, &parameter]{
	Miscellaneous::read_cone_orientation(cone, coneIfRot, parameters);
    });
    if(!parameters.isfullsky){
        Miscellaneous::TicketizeFunction(rank, ntasks, [=, &parameter, &rotm1, &thetay, &thetaz]{
            Miscellaneous::get_narrow_specs(parameters, rotm1, thetay, thetaz);
        });
    }
    // Read cosmology
    cosmology = Input::acquire(parameters, h, omegam, lboxmpch);
    if(rank == 0){
        std::cout<<"## Parameter file "<<parameters.paramfile<<std::endl;
        std::cout<<"## Evolution file "<<parameters.evolfile<<std::endl;
        std::cout<<"## Cosmology : h = "<<h<<", omega_m = "<<omegam<<", rhoch2 = "<<parameters.rhoch2<<std::endl;
        std::cout<<"## Box length : "<<lboxmpch<<" Mpc/h, 1 Mpc = "<<parameters.mpc<<" meters "<<std::endl;
        std::cout<<"## Coarse level : "<<parameters.ncoarse<<", number of cones = "<<parameters.ncones<<std::endl;
    }
    const double length = lboxmpch*parameters.mpc/h;

    // Construct homogeneous tree 

    Input::homogenize(octree.assign(nreference, zero));
    reference.append(Integrator::launch(center[zero], center[one], center[two], center[zero]+diameter/two, center[one], center[two]));

    // Propagate a photon in a homogeneous cosmology

    Integrator::integrate<-1>(reference, cosmology, octree, vobs0, length, EXTENT*std::pow(two, static_cast<uint>(std::log2(std::get<0>(cosmology).size()/std::pow(two, nreference)+one)+one)+one));
    cosmology = Input::correct(cosmology, reference);
    reference.fullclear();
    octree.fullclear();
   
    // Construct octree

    if (parameters.typefile == 2) { // ASCII
#ifdef VERBOSE
        if(rank == 0) 
	    std::cout<<"# Preparation mode : read ASCII files"<<std::endl;
#endif
	Create_octree::PreparationASCII(octree, parameters, ntasks, rank, cone, coneIfRot, rotm1, conefile, microsphere, h,  omegam, lboxmpch, amin, cosmology);
    } else if (parameters.typefile == 1) { // HDF5
#ifdef VERBOSE
        if(rank == 0) 
	    std::cout<<"# Preparation mode : read HDF5 files"<<std::endl;
#endif
	if(parameters.inputtype == "cells"){
	    Create_octree::PreparationHDF5_from_cells(octree, parameters, ntasks, rank, cone, coneIfRot, rotm1, thetay, thetaz, conefile, microsphere, h,  omegam, lboxmpch, amin, cosmology);
	} else if (parameters.inputtype == "particles"){ 
	    Create_octree::PreparationHDF5_from_particles(octree, parameters, ntasks, rank, cone, coneIfRot, thetay, thetaz, conefile, microsphere);
	} else{
	    std::cout<<"# Please choose 'cells' or 'particles' for inputtype"<<std::endl;
	    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	    std::terminate();
	}
    } else if (parameters.typefile == 0) { // Binary
#ifdef VERBOSE
        if(rank == 0) 
	    std::cout<<"# Preparation mode : read Binary files"<<std::endl;
#endif
	Create_octree::PreparationBinary(octree, parameters, ntasks, rank, cone, conefile, microsphere, filetree, h,  omegam, lboxmpch, amin, cosmology);
    } else{
        if(rank == 0) {
	    std::cout<<"# Please choose Preparation or Propagation mode"<<std::endl;
	    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	    std::terminate();
	}
    }

    // Finalization
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    if(rank == 0)
        std::cout<<"# Run completed !"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
