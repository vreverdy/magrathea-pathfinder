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
#include "generate_cones.h"
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
    using point = std::array<real, 3>;
    using extent = std::ratio<EXTENT>;

    static constexpr uint zero = 0;
    static constexpr uint two = 2;
    static constexpr point center({{0, 0, 0}});
    static constexpr real diameter = static_cast<real>(extent::num)/static_cast<real>(extent::den);
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
    Generate_cones::ReadParamFile(parameters, parameter);
    // Initialization
    HyperSphere<dimension, point> sphere(center, diameter/two);
    std::vector<Cone<point> > cone(parameters.ncones);
    std::vector<Cone<point> > coneIfRot(parameters.ncones);
    real thetay(0), thetaz(0);
    std::array< std::array< double, 3 >, 3 > rotm1 = {{zero}};
 
    if(rank == 0) 
	std::cout<<"#### MAGRATHEA_PATHFINDER "<<std::endl; 
    // Generate cones
    if(parameters.isfullsky){
	if(rank == 0){
	    std::cout<<"## Fullsky cone"<<std::endl;
	}
        Generate_cones::GenerateFullskyCones(parameters.ncones, cone, coneIfRot, sphere);
    } else{
	if(rank == 0){
	    std::cout<<"## Narrow cone"<<std::endl;
	}
        Generate_cones::GenerateNarrowCones(parameters, cone, coneIfRot, sphere, rotm1, thetay, thetaz);
    }

    // Write cone orientations in txt file
    Miscellaneous::write_cone_orientation(cone, coneIfRot, parameters);

    // Finalization
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    if(rank == 0)
        std::cout<<"# Run completed !"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
