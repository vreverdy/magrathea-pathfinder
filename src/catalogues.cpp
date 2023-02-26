/* ******************************* CATALOGUES ********************************
 */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Catalogues
// DESCRIPTION :    Main function of the raytracer to create catalogues
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           catalogues.cpp
/// \brief          Main function of the raytracer to create catalogues
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
/// (michel-andres.breton@obspm.fr) \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

// ------------------------------ PREPROCESSOR ------------------------------ //
#include <ctime>
// Include C++
#include <algorithm>
#include <array>
#include <atomic>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
// Include libs
#include <mpi.h>
// Include project
#include "catalogues.h"
#include "cone.h"
#include "input.h"
#include "integrator.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/hypercube.h"
#include "magrathea/hypersphere.h"
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/timer.h"
#include "miscellaneous.h"
#include "output.h"
#include "utility.h"
// Octree
#ifdef VELOCITYFIELD
#include "gravity2.h" // with velocity slots
#else
#include "gravity.h"
#endif
// Include C++ HDF5 project
#include "TReadHDF5.h"

// In miscellaneous.h we use std::sample. For versions of gcc lower than 7.x,
// please use std::experimental::sample. This can be done in the Makefile

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
int main(int argc, char *argv[]) {
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
    static constexpr uint two = 2;
    static constexpr uint dimension = 3;
    static constexpr uint nreference = 5; // Used to set homogeneous octree
    static constexpr real rposition =
        static_cast<real>(position::num) / static_cast<real>(position::den);
    static constexpr point center({{rposition, rposition, rposition}});
    static constexpr real diameter =
        static_cast<real>(extent::num) / static_cast<real>(extent::den);
    static const std::string all = "all";
    static const std::string outputsep = "_";    // (Separator used in file names)
    static const std::string outputint = "%05d"; // (Integer format in file names)
    static const std::string outputopening =
        "%8.6f"; // (Opening angle format in file names)
    static const std::string namelist =
        argc > 1 ? std::string(argv[1]) : std::string("raytracer.txt");

    // Parameters
    std::map<std::string, std::string> parameter;
    integer nthreads = std::thread::hardware_concurrency();
    integer ntasks = nthreads * zero;
    integer rank = zero;

    // Message passing interface
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Read parameter file
    Miscellaneous::TicketizeFunction(
        rank, ntasks, [=, &parameter] { parameter = Input::parse(namelist); });
    // Convert strings and put it in struct
    Catalogues::ReadParamFile(parameters, parameter);
    // Initialization
    FileList conefile(parameters.conefmt, zero, parameters.ncones, zero,
                      parameters.conedir);
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>,
                      Gravity<floating, dimension>, dimension, position, extent>
        octree;
    HyperSphere<dimension, point> sphere(center, diameter / two);
    std::vector<std::string> filelist;
    std::vector<Cone<point>> cone(parameters.ncones);
    std::vector<Cone<point>> coneIfRot(parameters.ncones);
    std::array<std::vector<real>, two + two> cosmology;
    Evolution<Photon<real, dimension>> reference;
    real h = zero;
    real omegam = zero;
    real lboxmpch = zero;
    std::ofstream stream;
    std::string filename;
    std::array<std::array<double, 3>, 3> rotm1 = {{zero}};
    const point vobs = {parameters.v0x, parameters.v0y, parameters.v0z};
    const point vobs0 = {0, 0,
                         0}; // No peculiar velocity for homogeneous quantities

    if (rank == 0)
        std::cout << "#### MAGRATHEA_PATHFINDER " << std::endl;
    // Generate cones
    if (parameters.use_previous_catalogues ==
        0) { // Only useful to assign sources to cone
        Miscellaneous::TicketizeFunction(
            rank, ntasks, [=, &cone, &coneIfRot, &parameter] {
                Miscellaneous::read_cone_orientation(cone, coneIfRot, parameters);
            });
    }

    if (!parameters.isfullsky) {
        real thetay, thetaz;
        Miscellaneous::TicketizeFunction(
            rank, ntasks, [=, &parameter, &rotm1, &thetay, &thetaz] {
                Miscellaneous::get_narrow_specs(parameters, rotm1, thetay, thetaz);
            });
    }
    // Read cosmology
    cosmology = Input::acquire(parameters, h, omegam, lboxmpch);
    if (rank == 0) {
        std::cout << "## Parameter file " << parameters.paramfile << std::endl;
        std::cout << "## Evolution file " << parameters.evolfile << std::endl;
        std::cout << "## Cosmology : h = " << h << ", omega_m = " << omegam
                  << ", rhoch2 = " << parameters.rhoch2 << std::endl;
        std::cout << "## Box length : " << lboxmpch
                  << " Mpc/h, 1 Mpc = " << parameters.mpc << " meters "
                  << std::endl;
        std::cout << "## Coarse level : " << parameters.ncoarse
                  << ", number of cones = " << parameters.ncones << std::endl;
    }
    const double length = lboxmpch * parameters.mpc / h;

    // Construct homogeneous tree
    Input::homogenize(octree.assign(nreference, zero));
    reference.append(Integrator::launch(center[zero], center[one], center[two],
                                        center[zero] + diameter / two,
                                        center[one], center[two]));
    // Propagate a photon in a homogeneous cosmology
    Integrator::integrate<-1>(
        reference, cosmology, octree, vobs0, length,
        EXTENT * std::pow(two, static_cast<uint>(
                                   std::log2(std::get<0>(cosmology).size() /
                                                 std::pow(two, nreference) +
                                             one) +
                                   one) +
                                   one));
    cosmology = Input::correct(cosmology, reference);
    octree.fullclear();

    // Execution
#ifndef VELOCITYFIELD
    if (rank == 0)
        std::cout << "# Compute relativistic catalogs" << std::endl;
    reference.fullclear();
    for (uint icone = parameters.firstcone; icone < parameters.lastcone + one;
         ++icone) {
        if (icone % static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
            std::vector<std::array<double, 18>> previous_catalogue;
            if (parameters.use_previous_catalogues) {
                Miscellaneous::ReadFromCat(
                    icone, parameters,
                    previous_catalogue); // For 'rejected', read the id here, and choose
                                         // only those that are inside the cone
#ifdef VERBOSE
                std::cout << "# Rank " << rank << " cone " << icone
                          << " number of sources in the cone : "
                          << previous_catalogue.size() << std::endl;
#endif
                if (previous_catalogue.size() == 0) {
                    std::cout << "# Rank " << rank << " cone " << icone
                              << " no data, skip this cone" << std::endl;
                    continue;
                }
            }
            // Load octree
            Miscellaneous::loadOctree(icone, octree, conefile);
            // Set observer position
            point observer;
            observer[0] = 0;
            observer[1] = 0;
            observer[2] = 0;

            if (parameters.use_previous_catalogues == 1 ||
                parameters.use_previous_catalogues ==
                    3) { // Rerun previously computed catalogs
                // Set filename
                std::string conetype = parameters.isfullsky ? "fullsky" : "narrow";
                std::string sourcetype = parameters.halos ? "halos" : "part";
                std::string jacobinfo =
                    (parameters.beam == "bundle")
                        ? Output::name(
                              parameters.stop_bundle, outputsep, parameters.plane,
                              outputsep,
                              std::make_pair(outputopening, parameters.openingmin))
                        : parameters.beam;
                filename = Output::name(parameters.outputdir, parameters.outputprefix,
                                        outputsep, conetype, outputsep, jacobinfo,
                                        outputsep, sourcetype);
                if (parameters.use_previous_catalogues == 1) {
                    filename = Output::name(filename, outputsep,
                                            std::make_pair(outputint, icone));
                    // Re-run sources in already computed catalogue
                    Catalogues::relCat_with_previous_cat(vobs, filename, observer,
                                                         previous_catalogue, parameters,
                                                         cosmology, octree, length, h);
                } else if (parameters.use_previous_catalogues ==
                           3) { // Rerun previously computed catalogs
                    // Set filename
                    filename = Output::name(filename, outputsep, "flexion", outputsep,
                                            std::make_pair(outputint, icone));
                    // Re-run sources with flexion raytracing
                    Catalogues::relCat_with_previous_cat_flexion(
                        vobs, filename, observer, previous_catalogue, parameters,
                        cosmology, octree, length, h);
                }
            } else if (parameters.use_previous_catalogues == 0 ||
                       parameters.use_previous_catalogues ==
                           2) { // If need to compute catalogue or re-run
                                // non-converged sources

                std::vector<std::array<double, 8>> caractVect_source;
                // Read source files (usually depends on some convention in the
                // filenames)
                if (parameters.typefile == 1) {
                    Catalogues::ReadParticlesHDF5(rank, parameters, caractVect_source);
                } else if (parameters.typefile == 2) {
                    Catalogues::ReadParticlesASCII(rank, parameters, caractVect_source);
                } else {
                    if (rank == 0) {
                        std::cout
                            << "Targets are only available with the HDF5 or ASCII format"
                            << std::endl;
                        std::cout << "# Error at file " << __FILE__
                                  << ", line : " << __LINE__ << std::endl;
                        std::terminate();
                    }
                }

                std::vector<std::array<double, 8>> targets_position;
                // Assign sources to cone
                targets_position =
                    Miscellaneous::getTargets(caractVect_source, cone[icone], cone);
#ifdef VERBOSE
                std::cout << "# Rank " << rank << " cone " << icone
                          << " HALOS : number of sources in the cone : "
                          << targets_position.size() << " of "
                          << caractVect_source.size() << std::endl;
#endif
                Miscellaneous::fullclear_vector(caractVect_source);
                // If re-run rejected sources
                if (parameters.use_previous_catalogues == 2) {
                    std::vector<std::array<double, 8>> targets_position_tmp;
                    // Sort sources from previously computed catalogue with index
                    std::sort(
                        previous_catalogue.begin(), previous_catalogue.end(),
                        [](const std::array<double, 18> &a,
                           const std::array<double, 18> &b) { return a[0] < b[0]; });
                    // Sort sources from full dataset with index
                    std::sort(targets_position.begin(), targets_position.end(),
                              [](const std::array<double, 8> &a,
                                 const std::array<double, 8> &b) { return a[6] < b[6]; });

                    for (uint i = 0; i < previous_catalogue.size(); i++) {
                        std::array<double, 8> point_tmp;
                        point_tmp[6] = previous_catalogue[i][0];
                        const unsigned long int marked = std::distance(
                            std::begin(targets_position),
                            std::upper_bound(std::begin(targets_position),
                                             std::end(targets_position), point_tmp,
                                             [](const std::array<double, 8> &first,
                                                const std::array<double, 8> &second) {
                                                 return first[6] < second[6];
                                             }));
                        // If we find a source that was rejected, then append the vector to
                        // re-run it from the beginning
                        targets_position_tmp.push_back(
                            targets_position[marked - (marked > 0)]);
                    }
                    Miscellaneous::fullclear_vector(targets_position);
                    targets_position = targets_position_tmp;
                    Miscellaneous::fullclear_vector(targets_position_tmp);
                }
                // Set filename
                std::string conetype = parameters.isfullsky ? "fullsky" : "narrow";
                std::string sourcetype = parameters.halos ? "halos" : "part";
                std::string jacobinfo =
                    (parameters.beam == "bundle")
                        ? Output::name(
                              parameters.stop_bundle, outputsep, parameters.plane,
                              outputsep,
                              std::make_pair(outputopening, parameters.openingmin))
                        : parameters.beam;
                if (parameters.use_previous_catalogues == 2) {
                    filename = Output::name(parameters.outputdir, "/rejected/",
                                            parameters.outputprefix, outputsep, conetype,
                                            outputsep, jacobinfo, outputsep, sourcetype,
                                            outputsep, std::make_pair(outputint, icone));
                } else {
                    filename = Output::name(parameters.outputdir, parameters.outputprefix,
                                            outputsep, conetype, outputsep, jacobinfo,
                                            outputsep, sourcetype, outputsep,
                                            std::make_pair(outputint, icone));
                }
                // Run the root-finding method to connect the observer and sources
                Catalogues::relCat(vobs, rotm1, filename, observer, targets_position,
                                   previous_catalogue, parameters, cosmology, octree,
                                   length, h);
            } else {
                std::cout << "# WARNING : If 'use_previous_catalogues' is non-zero, "
                             "must be equal to 1 or 2"
                          << std::endl;
                std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                          << std::endl;
                std::terminate();
            }
        }
    } // for icone
#else
    std::cout << "Do no forget to disable the -DVELOCITYFIELD option in the "
                 "Makefile if you want to compute relativistic catalogues"
              << std::endl;
    std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
              << std::endl;
    std::terminate();
#endif

    // Finalization
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    if (rank == 0)
        std::cout << "# Run completed !" << std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
