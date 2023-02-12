/* ******************************* RAYS ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Rays
// DESCRIPTION :    Main function of the raytracer
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           rays.cpp
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
#include "rays.h"
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
    using evolution = Evolution<Photon<real, 3> >;
    static constexpr uint zero = 0;
    static constexpr uint one = 1;
    static constexpr real rone = 1.; 
    static constexpr uint two = 2;
    static constexpr uint three = 3;
    static constexpr uint four = 4;
    static constexpr uint five = 5;
    static constexpr uint dimension = 3;
    static constexpr uint nreference = 5; // Used to set homogeneous octree
    static constexpr real pi = Constants<real>::pi();
    static constexpr real rposition = static_cast<real>(position::num)/static_cast<real>(position::den);
    static constexpr point center({{rposition, rposition, rposition}});
    static constexpr real diameter = static_cast<real>(extent::num)/static_cast<real>(extent::den);
    static constexpr uint digits = std::numeric_limits<real>::max_digits10;
    static constexpr char dot = '.';
    static constexpr char dotc = 'd';
    static const std::string outputsep = "_"; // (Separator used in file names)
    static const std::string outputint = "%05d"; // (Integer format in file names)
    static const std::string outputopening = "%8.6f"; // (Opening angle format in file names)
    static const std::string outputsuffix = ".txt"; // (Suffix for text result files)
    static const std::string outputrk4 = "rk4"; // (Keyword used in file names)
    static const std::string outputtree = "octree"; // (Keyword used in file names)
    static const std::string outputstat = "stat"; // (Keyword used in file names)
    static const std::string all = "all";
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
    Rays::ReadParamFile(parameters, parameter);
    // Initialization
    std::vector<std::string> statistics = (parameters.statistic == all) ? (std::vector<std::string>({"distance", "distance2", "homogeneous", "inhomogeneous"})) : (std::vector<std::string>({parameters.statistic}));
    uint statcnt = statistics.size();
    uint interpcase = zero;
    uint statcase = zero;
    std::string interp;
    std::string stat;
    FileList conefile(parameters.conefmt, zero, parameters.ncones, zero, parameters.conedir);
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>, Gravity<floating, dimension>, dimension, position, extent> octree;
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>, Gravity<floating, dimension>, dimension, position, extent> homotree;
    std::uniform_real_distribution<real> distribution(zero, one);
    HyperSphere<dimension, point> sphere(center, diameter/two);
    HyperSphere<dimension, point> microsphere(center, rone/parameters.microcoeff);
    std::vector<std::string> filelist;
    std::vector<point> tiling(parameters.ncones);
    std::vector<point> tilingbis(parameters.ncones);
    std::vector<Cone<point> > cone(parameters.ncones);
    std::vector<Cone<point> > coneIfRot(parameters.ncones);
    std::vector<evolution> trajectory(parameters.ntrajectories);
    std::vector<Photon<real, dimension> > photons(parameters.ntrajectories);
    std::vector<real> random(parameters.ntrajectories);
    std::array<std::vector<real>, two+two> cosmology;
    Photon<real, dimension> photon;
    Evolution<Photon<real, dimension> > reference;
    real h = zero;
    real omegam = zero;
    real lboxmpch = zero;
    real amin = zero;
    uint nbundle = zero;
    real opening = zero;
    std::ofstream stream;
    std::string filename;
    std::mutex mutex;
    uint statmod = zero;
    uint statlength = zero;
    uint statsize = zero;
    uint statgsize = zero;
    std::vector<real> statrefx;
    std::vector<real> statrefy;
    std::vector<real> statrefz;
    std::vector<std::vector<real> > statx;
    std::vector<std::vector<real> > staty;
    std::vector<real> statmean;
    std::vector<real> statstd;
    std::vector<real> statgmean;
    std::vector<real> statgstd;
    const point vobs0 = {0,0,0}; // No peculiar velocity for homogeneous quantities
    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());

    if(rank == 0) 
	std::cout<<"#### MAGRATHEA_PATHFINDER "<<std::endl; 
    // Generate cones
    Miscellaneous::TicketizeFunction(rank, ntasks, [=, &cone, &coneIfRot, &parameter]{
	Miscellaneous::read_cone_orientation(cone, coneIfRot, parameters);
    });
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
    octree.fullclear();

    // Construct octree
    Miscellaneous::loadOctree(rank, octree, conefile);


    // Execution
    if(rank == 0)
        std::cout<<"# Compute rays and distance statistics"<<std::endl;
    // Initialization
    reference.fullclear();
    const double invlength = rone/length;
    // Create FLRW octree
    homotree.assign(parameters.ncoarse/two, zero);
    // Launch a photon toward the x-axis (will be used as reference FLRW ray)
    photon = Integrator::launch(center[zero], center[one], center[two], center[zero]+diameter/two, center[one], center[two]);
    uint ntrajectoriesMax(0);
    // Launch photons toward random directions
    if(parameters.ray_targets == "random"){
        ntrajectoriesMax = parameters.ntrajectories;
        for (uint itrajectory = zero; itrajectory < parameters.ntrajectories; ++itrajectory) {
            photons[itrajectory] = Integrator::launch(microsphere, coneIfRot[rank], coneIfRot, engine1, distribution);
            random[itrajectory] = distribution(engine1)*two*pi;
        }
    // Launch photons toward observed sources (with mass threshold for haloes)
    } else if (parameters.ray_targets == "catalogue"){
	std::vector < std::array< double, 18 > > previous_catalogue;
	// Read catalogue
        Miscellaneous::ReadFromCat(rank, parameters, previous_catalogue);
	// Select sources within mass bin
	previous_catalogue.erase(std::remove_if(previous_catalogue.begin(), previous_catalogue.end(), [](const std::array< double, 18 >& elem) {return (elem[17] >= parameters.massmax) || (elem[17] < parameters.massmin);}), previous_catalogue.end());
	ntrajectoriesMax = previous_catalogue.size();
#ifdef VERBOSE
	std::cout<<"Rank : "<<rank<<" number of sources "<<ntrajectoriesMax<<std::endl;
#endif
	photons.resize(ntrajectoriesMax);
	random.resize(ntrajectoriesMax);
	// Initialise photons
        for (uint itrajectory = zero; itrajectory < ntrajectoriesMax; ++itrajectory) {
            photons[itrajectory] = Integrator::launch(center[zero], center[one], center[two], previous_catalogue[itrajectory][3], previous_catalogue[itrajectory][4]); // 3 and 4 are the observed angle in catalogues
            random[itrajectory] = distribution(engine1)*two*pi;
        }
    } else{
        std::cout<<"# WARNING : ray_target must be 'random' or 'catalog'"<<std::endl;
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	std::terminate();
    }
    // Loop over configurations
    for (uint ibundle = zero; ibundle < parameters.nbundlecnt; ++ibundle) {
        nbundle = parameters.nbundlemin+ibundle;
        for (uint iopening = zero; iopening < parameters.openingcnt; ++iopening) {
            opening = (iopening+one)*parameters.openingmin;
            interp = parameters.stop_ray;
            for (uint istat = zero; istat < statcnt; ++istat) {
                stat = statistics[istat];
                // Reference
		filename = Output::name(parameters.outputprefix, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interp, outputsep, std::make_pair(outputint, rank));
                std::replace(filename.begin(), filename.end(), dot, dotc);
		filename = Output::name(parameters.outputdir, filename);
                reference = Integrator::propagate<-1>(photon, nbundle, opening, real(), interp, cosmology, homotree, vobs0, length, EXTENT*parameters.nsteps*(one << (parameters.ncoarse-parameters.ncoarse/two))*two, real(), std::signbit(parameters.savemode) ? Output::name() : Output::name(filename, outputsuffix));
                // Integration without statistics
                if (parameters.makestat == zero) {
                    Utility::parallelize(ntrajectoriesMax, [=, &photons, &nbundle, &opening, &random, &interp, &cosmology, &octree, &vobs0, &length, &amin, &filename, &reference](const uint i){Integrator::propagate<1>(photons[i], nbundle, opening, random[i], interp, cosmology, octree, vobs0, length, parameters.nsteps, amin, std::signbit(parameters.savemode) ? Output::name() : Output::name(parameters.savemode ? Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, outputint) : Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, std::make_pair(outputint, zero)), outputsuffix), reference);});
                } else {
                    // Integration with statistics
                    // Clear statistics arrays
                    statrefx.clear();
                    statrefy.clear();
                    statrefz.clear();
                    statx.clear();
                    staty.clear();
                    statmean.clear();
                    statstd.clear();
                    statgmean.clear();
                    statgstd.clear();
                    // Interpolation case
                    if (interp == "redshift") {
                        interpcase = zero;
                    } else if (interp == "a") {
                        interpcase = one;
                    } else if (interp == "t") {
                        interpcase = two;
                    } else if (interp == "r") {
                        interpcase = three;
		    } else if (interp == "lambda"){
			interpcase = four;
                    } else {
                        interpcase = zero;
                    }
                    // Statistics case
                    if (stat == "distance") {
                        statcase = zero;
                    } else if (stat == "distance2") {
                        statcase = one;
                    } else if (stat == "homogeneous") {
                        statcase = two;
                    } else if (stat == "inhomogeneous") {
                        statcase = three;
                    } else if (stat == "invdistance") {
                        statcase = four;
                    } else if (stat == "invdistance2") {
                        statcase = five;
                    } else {
                        statcase = zero;
                    }
                    // Integration
                    Utility::parallelize(ntrajectoriesMax, [=, &photons, &nbundle, &opening, &random, &interp, &cosmology, &octree, &vobs0, &length, &amin, &filename, &reference, &mutex, &interpcase, &statcase, &statx, &staty](const uint i){
                        evolution result = Integrator::propagate(photons[i], nbundle, opening, random[i], interp, cosmology, octree, vobs0, length, parameters.nsteps, amin, std::signbit(parameters.savemode) ? Output::name() : Output::name(parameters.savemode ? Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, outputint) : Output::name(filename, outputsep, std::make_pair(outputint, i), outputsep, std::make_pair(outputint, zero)), outputsuffix), reference);
                        std::vector<std::vector<real> > tmp(two, std::vector<real>(result.size()));
                        for (uint j = zero; j < result.size(); ++j) { // Loop on all central rays
                            if (interpcase == zero) {
                                tmp[zero][j] = result[j].redshift();
                            } else if (interpcase == one) {
                                tmp[zero][j] = result[j].t();
                            } else if (interpcase == two) {
                                tmp[zero][j] = one/result[j].a() - one;
                            } else if (interpcase == three) {
                                tmp[zero][j] = std::sqrt(std::pow(result[j].x()-result[zero].x(), two)+std::pow(result[j].y()-result[zero].y(), two)+std::pow(result[j].z()-result[zero].z(), two));
                            } else if (interpcase == four){
				tmp[zero][j] = result[j].lambda();
			    }
                            if (statcase == zero) {
                                tmp[one][j] = result[j].distance()*invlength;
                            } else if (statcase == one) {
                                tmp[one][j] = result[j].distance()*invlength*result[j].distance()*invlength;
                            } else if (statcase == two) {
                                tmp[one][j] = result[zero].a()/result[j].a()-one;
                            } else if (statcase == three) {
                                tmp[one][j] = result[j].redshift();
                            } else if (statcase == four) {
                                tmp[one][j] = rone/(result[j].distance()*invlength);
                            } else if (statcase == five) {
                                tmp[one][j] = rone/(result[j].distance()*invlength*result[j].distance()*invlength);
                            } 
                        }
                        if (result.size() > zero) {
                            mutex.lock();
                            statx.emplace_back(std::move(tmp[zero]));
                            staty.emplace_back(std::move(tmp[one]));
                            mutex.unlock();
                        }
                    }); 
		    octree.fullclear();
                    // Transfer reference
                    reference.container().erase(std::remove_if(reference.container().begin(), reference.container().end(), [=, &amin](const Photon<real, dimension>& p){return std::isnormal(amin) && (p.a() < amin);}), reference.container().end());
                    statmod = std::max(one, static_cast<uint>(reference.size())/std::max(parameters.nstat, one));
                    for (uint j = zero; j < reference.size(); ++j) {
                        if (j%statmod == zero) {
                            if (interpcase == zero) {
                                statrefx.emplace_back(reference[j].redshift());
                            } else if (interpcase == one) {
                                statrefx.emplace_back(reference[j].t());
                            } else if (interpcase == two) {
                                statrefx.emplace_back(one/reference[j].a() - one );
                            } else if (interpcase == three) {
                                statrefx.emplace_back(std::sqrt(std::pow(reference[j].x()-reference[zero].x(), two)+std::pow(reference[j].y()-reference[zero].y(), two)+std::pow(reference[j].z()-reference[zero].z(), two)));
                            } else if (interpcase == four) {
				statrefx.emplace_back(reference[j].lambda());
			    }
                            if (statcase == zero) {
                                statrefy.emplace_back(reference[j].distance()*invlength);
                            } else if (statcase == one) {
                                statrefy.emplace_back(reference[j].distance()*invlength*reference[j].distance()*invlength);
                            } else if (statcase == two) {
                                statrefy.emplace_back(reference[zero].a()/reference[j].a()-one);
                            } else if (statcase == three) {
                                statrefy.emplace_back(reference[j].redshift());
                            } else if (statcase == four) {
                                statrefy.emplace_back(rone/(reference[j].distance()*invlength));
                            } else if (statcase == five) {
                                statrefy.emplace_back(rone/(reference[j].distance()*invlength*reference[j].distance()*invlength));
                            }
                            statrefz.emplace_back(reference[j].redshift());
                        }
                    }
                    // Resize
                    statlength = statrefx.size();
                    statsize = statx.size();
                    statmean.resize(statlength, real());
                    statstd.resize(statlength, real());
                    statgmean.resize(statlength, real());
                    statgstd.resize(statlength, real());
                    statmean.shrink_to_fit();
                    statstd.shrink_to_fit();
                    statgmean.shrink_to_fit();
                    statgstd.shrink_to_fit();
                    // Reinterpolate
                    Utility::parallelize(statsize, [=, &statrefx, &statx, &staty](const uint j){staty[j] = Utility::reinterpolate(statrefx, statx[j], staty[j]);}); // WARNING ! Interpolates outside of range !
                    // Reduction
		    // Compte mean and standard deviation
                    MPI_Allreduce(&statsize, &statgsize, one, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
                    std::fill(statmean.begin(), statmean.end(), real());
                    for (uint j = zero; j < statsize; ++j) {
                        for (uint k = zero; k < statlength; ++k) {
                            statmean[k] += staty[j][k];
                        }
                    }
                    MPI_Allreduce(statmean.data(), statgmean.data(), statlength, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    for (uint k = zero; k < statlength; ++k) {
                        statmean[k] /= real(statsize);
                        statgmean[k] /= real(statgsize);
                    }
                    std::fill(statstd.begin(), statstd.end(), real());
                    for (uint j = zero; j < statsize; ++j) {
                        for (uint k = zero; k < statlength; ++k) {
                            statstd[k] += (staty[j][k]-statgmean[k])*(staty[j][k]-statgmean[k]);
                        }
                    }
                    MPI_Allreduce(statstd.data(), statgstd.data(), statlength, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    std::fill(statstd.begin(), statstd.end(), real());
                    for (uint j = zero; j < statsize; ++j) {
                        for (uint k = zero; k < statlength; ++k) {
                            statstd[k] += (staty[j][k]-statmean[k])*(staty[j][k]-statmean[k]);
                        }
                    }
                    for (uint k = zero; k < statlength; ++k) {
                        statstd[k] = std::sqrt(statstd[k]/real(statsize-one));
                        statgstd[k] = std::sqrt(statgstd[k]/real(statgsize-one));
                    }
                    // Output statistics
		    filename = Output::name(parameters.outputprefix, outputsep, outputstat, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interp, outputsep, stat, outputsep, std::make_pair(outputint, rank));
                    std::replace(filename.begin(), filename.end(), dot, dotc); 
		    filename = Output::name(parameters.outputdir, filename);
                    stream.open(Output::name(filename, outputsuffix));
                    Output::save(stream, statrefz, statrefy, statmean, statstd, digits, statsize);
                    stream.close();
                    filename = Output::name(parameters.outputprefix, outputsep, outputstat, outputsep, std::make_pair(outputint, nbundle), outputsep, std::make_pair(outputopening, opening), outputsep, interp, outputsep, stat);
                    std::replace(filename.begin(), filename.end(), dot, dotc);
		    filename = Output::name(parameters.outputdir, filename);
                    if (rank == zero) {
#ifdef VERBOSE
			std::cout<<"#"<<Output::name(filename, outputsuffix)<<std::endl;
#endif
                        stream.open(Output::name(filename, outputsuffix));
                        Output::save(stream, statrefz, statrefy, statgmean, statgstd, digits, statgsize);
                        stream.close();
                    }
                } 
            }
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
