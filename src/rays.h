/* ********************************** RAYS ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Rays
// DESCRIPTION :    Some rays function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           rays.h
/// \brief          Some rays function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef RAYS_H_INCLUDED
#define RAYS_H_INCLUDED


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

#ifdef GCCBELOW7
#include <experimental/algorithm>
#endif
// Include libs
#include <mpi.h>
// Include project
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/hypersphere.h"
#include "magrathea/abstracthypersphere.h"
#include "magrathea/hypercube.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/timer.h"
#include "cone.h"
#include "utility.h"
#include "input.h"
#include "output.h"
#include "TReadHDF5.h"


using namespace magrathea;

struct parameters_t {
    // Common
    real buffer;
    std::string celldir;
    std::string conedir;
    std::string conefmt;
    std::string evolfile;
    uint isfullsky;
    real mpc; 
    uint ncoarse;
    uint ncones;
    std::string paramfile;
    real rhoch2;
    uint seed;
    uint typefile;

    // Specific
    std::string base;
    uint halos;
    uint makestat;
    uint massmin;
    uint massmax;
    uint microcoeff;
    uint nbundlecnt;
    uint nbundlemin;
    uint nstat;
    uint nsteps;
    uint ntrajectories;
    uint openingcnt;
    real openingmin;
    std::string outputdir;
    std::string outputprefix;
    std::string ray_targets;
    integer savemode;
    std::string statistic;
    std::string stop_ray;

} parameters;

class Rays {

	//Methodes
	public:
	
	// Read parameter file
	template <class Parameters, class Map> static void ReadParamFile(Parameters& parameters, Map& parameter);

};

// Read parameter file
/// \brief          Read parameter file.
/// \details        Read and put in a structure the parameters.
/// \tparam         Parameters structure type
/// \tparam         Map map type
/// \param[in,out]  parameters Structure containing the parameters.
/// \param[in]      parameter Contains parameters to be rewritten
template <class Parameters, class Map> 
void Rays::ReadParamFile(Parameters& parameters, Map& parameter){

    parameters.statistic = parameter["statistic"];
    parameters.conefmt = parameter["conefmt"];
    parameters.ncones = std::stoul(parameter["ncones"]);
    parameters.buffer = std::stod(parameter["buffer"]);
    parameters.conedir = parameter["conedir"];
    parameters.seed = std::stoul(parameter["seed"]);
    parameters.microcoeff = std::stoul(parameter["microcoeff"]);
    parameters.ntrajectories = std::stoul(parameter["ntrajectories"]);
    parameters.isfullsky = std::stoul(parameter["isfullsky"]);
    parameters.paramfile = parameter["paramfile"];
    parameters.evolfile = parameter["evolfile"];
    parameters.rhoch2 = std::stod(parameter["rhoch2"]);
    parameters.mpc = std::stod(parameter["mpc"]);
    parameters.ncoarse = std::stoul(parameter["ncoarse"]);
    parameters.ray_targets = parameter["ray_targets"];
    parameters.massmin = std::stoul(parameter["massmin"]);
    parameters.massmax = std::stoul(parameter["massmax"]);
    parameters.nbundlemin = std::stoul(parameter["nbundlemin"]);
    parameters.nbundlecnt = std::stoul(parameter["nbundlecnt"]);
    parameters.openingmin = std::stod(parameter["openingmin"]);
    parameters.openingcnt = std::stoul(parameter["openingcnt"]);
    parameters.stop_ray = parameter["stop_ray"];
    parameters.outputdir = parameter["outputdir"];
    parameters.outputprefix = parameter["outputprefix"];
    parameters.nsteps = std::stoul(parameter["nsteps"]);
    parameters.savemode = std::stol(parameter["savemode"]);
    parameters.makestat = std::stoul(parameter["makestat"]);
    parameters.nstat = std::stoul(parameter["nstat"]);
    parameters.celldir = parameter["celldir"];
    parameters.typefile = std::stoul(parameter["typefile"]);
    parameters.base = parameter["base"];
    parameters.halos = std::stoul(parameter["halos"]);

}




#endif // RAYS_H_INCLUDED
