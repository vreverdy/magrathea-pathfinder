/* ******************************* OBSERVER_VELOCITY ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Observer_velocity
// DESCRIPTION :    Main function of the raytracer to compute the observer's velocity
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           observer_velocity.cpp
/// \brief          Main function of the raytracer to compute the observer's velocity
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
#include "observer_velocity.h"
// Octree
#ifdef VELOCITYFIELD
#include "gravity2.h" // with velocity slots
#else
#include "gravity.h"
#endif
// Include C++ HDF5 project
#include "TReadHDF5.h"

// In miscellaneous.h we use std::sample. For versions of gcc lower than 7.x, please use std::experimental::sample. Check the Makefiles

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
    using element = std::pair<SimpleHyperOctreeIndex<indexing, 3>, Gravity<floating, 3> >;
    static constexpr uint zero = 0;
    static constexpr uint dimension = 3;
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
    Observer_velocity::ReadParamFile(parameters, parameter);
    // Initialization
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>, Gravity<floating, dimension>, dimension, position, extent> octree;

    if(rank == 0) 
	std::cout<<"#### MAGRATHEA_PATHFINDER "<<std::endl; 

    // Execution
    if(rank == 0)
	std::cout<<"# Compute the observer's peculiar velocity"<<std::endl;
    point observer; 
    observer[0] = 0;
    observer[1] = 0;
    observer[2] = 0;
    double aexp(0), unit_t(0), unit_l(0);
    static const double c = magrathea::Constants<double>::c();
    std::size_t found;
    std::vector<std::string> filelistprior;
    std::string partfile;

    // Get filenames in directory
    Miscellaneous::getFilesinDir(parameters.partdir, filelistprior); 
    for(uint ifiling = 0; ifiling < filelistprior.size(); ifiling++){
	// Get file which contain the keyword 'shell' in its name
	found = filelistprior[ifiling].find("shell");
	if (found!=std::string::npos){
	    TReadHDF5::getAttribute(parameters.partdir + filelistprior[ifiling], "metadata/ramses_info", "aexp", aexp);
	    // Get only the file at a = 1 (today). Due to the simulation time step it is generally not exactly 1
	    if(aexp >= 0.9995){
	        partfile = parameters.partdir + filelistprior[ifiling];
		// Get conversion factors
       		TReadHDF5::getAttribute(partfile, "metadata/ramses_info", "unit_l", unit_l);
       		TReadHDF5::getAttribute(partfile, "metadata/ramses_info", "unit_t", unit_t);
		break;
	    }
	}
    }
#ifdef VERBOSE
    std::cout<<"# Particle file : "<<partfile<<std::endl;
#endif
    uint lvlmin = parameters.ncoarse + log2(EXTENT);
    double factor = unit_l*1e-2/(c*unit_t);
    for(uint ilvl = lvlmin; ilvl < lvlmin+5; ilvl++){
	octree.clear();
	// If for some reason the file is not at a = 1, throw error
	if(aexp < 0.9995){
	    std::cout<<"# Bad file : need particles at aexp = 1."<<std::endl;
	    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	    std::terminate();
	} else{
	    std::vector<float> pos_part, vel_part;
	    // Take all the particles in file
	    const double fraction = 1;
	    // Get particles from file
	    TReadHDF5::fillVectors_part(fraction, partfile, "data", "position_part", pos_part, "velocity_part", vel_part);			    
	    SimpleHyperOctreeIndex<indexing, dimension> id;
	    double half = 0.5*EXTENT*pow(2, -static_cast<int>(ilvl));
	    // Create an octree with 8 cells around the observer at the centre
	    for(int iz = -1; iz <= 1; iz += 2){
		for(int iy = -1; iy <= 1; iy += 2){
		    for(int ix = -1; ix <= 1; ix += 2){
		        id = id.template compute<double,position,extent>(ilvl, half*ix, half*iy, half*iz);
			octree.append(element(id, Gravity<floating, 3>()));
		    }
		}
	    }
	    // Sort octree
	    octree.update();
	    // Compute velocity field in neighbour cells
	    if(parameters.velocity_field_v0 == "cic"){
		Observer_velocity::CreateOctreeVelocityWithCIC(octree, pos_part, vel_part);
	    } else if(parameters.velocity_field_v0 == "tsc"){
		Observer_velocity::CreateOctreeVelocityWithTSC(octree, pos_part, vel_part);
	    } else{
		std::cout<<"# Wrong velocity field name. Please type cic or tsc for velocity_field_v0"<<std::endl;
		std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		std::terminate();
	    }
	    // Finalize
	    const unsigned int lvlmax = (std::get<0>(*std::max_element(std::begin(octree), std::end(octree), [](const element& x, const element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level()); 
	    const unsigned int lvlmin = (std::get<0>(*std::min_element(std::begin(octree), std::end(octree), [](const element& x, const element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level());
	    // Normalise velocity field by the mass
	    Utility::parallelize(octree.size(), [=, &octree](const uint i){
		double mass = std::get<1>(octree[i]).rho();   
		mass = mass + (mass == 0);
		std::get<1>(octree[i]).dphidx() /=  mass;   
		std::get<1>(octree[i]).dphidy() /=  mass;   
		std::get<1>(octree[i]).dphidz() /=  mass; 
	    });
	    // Correct rho and v = 0. If rho = 0, get from parent 
	    for (unsigned int ilvl = lvlmin+1; ilvl <= lvlmax; ilvl++) {
		Utility::parallelize(octree.size(), [=, &octree](const uint i){
		    if(std::get<0>(octree[i]).level() == ilvl){
		        Gravity<floating, 3> data;
		        if (!std::isnormal(std::get<1>(octree[i]).rho())) {
		            data = std::get<1>(*octree.find(std::get<0>(octree[i]).parent()));
		            std::get<1>(octree[i]).dphidxyz() = {data.dphidx(), data.dphidy(), data.dphidz()};
		        }
		    }
	        });
	    }
	}
	double vx(0), vy(0), vz(0);
	// Sum contribution from neighbour cells
	for(uint i = 0; i < octree.size(); i++){
	    vx += std::get<1>(octree[i]).dphidx();
	    vy += std::get<1>(octree[i]).dphidy();
	    vz += std::get<1>(octree[i]).dphidz();
	}
	// Output the velocity field at the observer
	std::cout<<"# level = "<<ilvl<<" | scheme = "<<parameters.velocity_field_v0<<" | Observers peculiar velocity : v0x = "<<vx/octree.size()<<" v0y = "<<vy/octree.size()<<" v0z = "<<vz/octree.size()<<", in RAMSES UNITS"<<std::endl;
	std::cout<<"# level = "<<ilvl<<" | scheme = "<<parameters.velocity_field_v0<<" | Observers peculiar velocity : v0x = "<<factor*vx/octree.size()<<" v0y = "<<factor*vy/octree.size()<<" v0z = "<<factor*vz/octree.size()<<", divided by c, in SI (to use in parameter files)"<<std::endl;
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
