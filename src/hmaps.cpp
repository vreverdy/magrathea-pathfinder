/* ******************************* HMAPS ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Hmaps
// DESCRIPTION :    Main function of the raytracer to create Healpix maps
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           hmaps.cpp
/// \brief          Main function of the raytracer to create Healpix maps
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
#include "cone.h"
#include "hmaps.h"
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

// Include HEALPIX
#include <chealpix.h>
#include <healpix_map.h>
#include <pointing.h>

// In miscellaneous.h we use std::sample. For versions of gcc lower than 7.x,
// please use std::experimental::sample. Check the Makefiles

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
#ifdef VELOCITYFIELD
    using element =
        std::pair<SimpleHyperOctreeIndex<indexing, 3>, Gravity<floating, 3>>;
#endif
    static constexpr uint INDEX_LENSING = 0;
    static constexpr uint INDEX_LENSING_BORN = 1;
    static constexpr uint INDEX_DR = 2;
    static constexpr uint INDEX_DL = 3;
    static constexpr uint INDEX_DT = 4;
    static constexpr uint INDEX_DA = 5;
    static constexpr uint INDEX_DZ = 6;
    static constexpr uint INDEX_DS = 7;
    static constexpr uint INDEX_ISW = 8;
    static constexpr uint INDEX_DENSITY = 9;
    static constexpr uint INDEX_DENSITY_MAX = 10;
    static constexpr uint INDEX_NSTEPS = 11;
    static constexpr uint INDEX_POTENTIAL = 12;
    static constexpr uint INDEX_DEFLECTION = 13;
    static constexpr uint INDEX_FLEXION = 14;
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
    Hmaps::ReadParamFile(parameters, parameter);
    // Initialization
    FileList conefile(parameters.conefmt, zero, parameters.ncones, zero,
                      parameters.conedir);
    SimpleHyperOctree<real, SimpleHyperOctreeIndex<indexing, dimension>,
                      Gravity<floating, dimension>, dimension, position, extent>
        octree;
    HyperSphere<dimension, point> sphere(center, diameter / two);
    std::vector<Cone<point>> cone(parameters.ncones);
    std::vector<Cone<point>> coneIfRot(parameters.ncones);
    std::array<std::vector<real>, two + two> cosmology;
    Photon<real, dimension> photon;
    Evolution<Photon<real, dimension>> reference;
    real h = zero;
    real omegam = zero;
    real lboxmpch = zero;
    const point vobs0 = {0, 0,
                         0}; // No peculiar velocity for homogeneous quantities
    std::mt19937 engine1(parameters.seed > zero ? parameters.seed + rank
                                                : std::random_device()());

    if (rank == 0)
        std::cout << "#### MAGRATHEA_PATHFINDER " << std::endl;
    // Generate cones
    Miscellaneous::TicketizeFunction(
        rank, ntasks, [=, &cone, &coneIfRot, &parameter] {
            Miscellaneous::read_cone_orientation(cone, coneIfRot, parameters);
        });
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
    // Propagate a photon in a very refined homogeneous FLRW metric
    Integrator::integrate<-1>(
        reference, cosmology, octree, vobs0, length,
        EXTENT * std::pow(two, static_cast<uint>(
                                   std::log2(std::get<0>(cosmology).size() /
                                                 std::pow(two, nreference) +
                                             one) +
                                   one) +
                                   one));
    // Correct cosmology with photon
    cosmology = Input::correct(cosmology, reference);
    octree.fullclear();

    // Execution
    if (rank == 0)
        std::cout << "# Compute Healpix maps, nside = " << parameters.nside
                  << std::endl;
    std::vector<real> z_stop_vec(parameters.nb_z_maps);
    // Get redshift values from redshift range and number of redhsifts
    for (uint i = 0; i < parameters.nb_z_maps; i++) {
        z_stop_vec[i] = parameters.z_stop_min +
                        i * (parameters.z_stop_max - parameters.z_stop_min) /
                            (parameters.nb_z_maps - (parameters.nb_z_maps > 1));
    }
    // Get number of pixels in a fullsky map from with nside
    long npix = nside2npix(parameters.nside);

    std::vector<std::string> map_components;
    // Tokenize map types to put in vector
    Miscellaneous::Tokenize(parameters.map_components, map_components, ", ");
    std::vector<uint> index_components(map_components.size());
    uint nmaps(0);

    if (rank == 0)
        std::cout << "Components of Healpix maps:" << std::endl;
    // For each map type, assign a unique index. 'Lensing' types produce multiple
    // maps
    for (uint i = 0; i < map_components.size(); i++) {
        if (map_components[i] == "lensing") {
            index_components[i] = INDEX_LENSING;
            if (parameters.beam == "bundle") {
                if (rank == 0) {
                    std::cout << "Index: " << nmaps << ", component: kappa ('lensing')"
                              << std::endl;
                    std::cout << "Index: " << nmaps + 1
                              << ", component: gamma1 ('lensing')" << std::endl;
                    std::cout << "Index: " << nmaps + 2
                              << ", component: gamma2 ('lensing')" << std::endl;
                    std::cout << "Index: " << nmaps + 3 << ", component: w ('lensing')"
                              << std::endl;
                    std::cout << "Index: " << nmaps + 4
                              << ", component: inverse magnification ('lensing')"
                              << std::endl;
                }
                nmaps += 4;
            } else if (parameters.beam == "infinitesimal") {
                if (rank == 0) {
                    std::cout << "Index: " << nmaps << ", component: kappa ('lensing')"
                              << std::endl;
                    std::cout << "Index: " << nmaps + 1
                              << ", component: gamma1 ('lensing')" << std::endl;
                    std::cout << "Index: " << nmaps + 2
                              << ", component: gamma2 ('lensing')" << std::endl;
                    std::cout << "Index: " << nmaps + 3
                              << ", component: inverse magnification ('lensing')"
                              << std::endl;
                }
                nmaps += 3;
            } else {
                std::cout << "# WARNING: With 'lensing', please choose beam 'bundle' "
                             "or 'infinitesimal'"
                          << std::endl;
                std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                          << std::endl;
                std::terminate();
            }
        } else if (map_components[i] == "lensing_born") {
            index_components[i] = INDEX_LENSING_BORN;
            if (rank == 0) {
                std::cout << "Index: " << nmaps
                          << ", component: kappa ('lensing_born') WARNING: with "
                             "lensing_born we compute the jacobian matrix with an "
                             "infinitesimal beam (no bundle method available)"
                          << std::endl;
                std::cout << "Index: " << nmaps + 1
                          << ", component: gamma1 ('lensing_born')" << std::endl;
                std::cout << "Index: " << nmaps + 2
                          << ", component: gamma2 ('lensing_born')" << std::endl;
                std::cout << "Index: " << nmaps + 3
                          << ", component: inverse magnification ('lensing_born')"
                          << std::endl;
            }
            nmaps += 3;
        } else if (map_components[i] == "deflection") {
            index_components[i] = INDEX_DEFLECTION;
            if (rank == 0) {
                std::cout << "Index: " << nmaps
                          << ", component: alpha1 ('deflection') [dphi in spherical "
                             "coordinates]"
                          << std::endl;
                std::cout << "Index: " << nmaps + 1
                          << ", component: alpha2 ('deflection') [dtheta in spherical "
                             "coordinates]"
                          << std::endl;
            }
            nmaps++;
        } else if (map_components[i] == "flexion") {
            index_components[i] = INDEX_FLEXION; // "d111", "d112", "d122", "d211", "d212", "d222"
            if (rank == 0) {
                std::cout << "Index: " << nmaps
                          << ", component: D111 ('flexion') [in radians^-1]"
                          << std::endl;
                std::cout << "Index: " << nmaps + 1
                          << ", component: D112 ('flexion') [in radians^-1]"
                          << std::endl;
                std::cout << "Index: " << nmaps + 2
                          << ", component: D122 ('flexion') [in radians^-1]"
                          << std::endl;
                std::cout << "Index: " << nmaps + 3
                          << ", component: D211 ('flexion') [in radians^-1]"
                          << std::endl;
                std::cout << "Index: " << nmaps + 4
                          << ", component: D212 ('flexion') [in radians^-1]"
                          << std::endl;
                std::cout << "Index: " << nmaps + 5
                          << ", component: D222 ('flexion') [in radians^-1]"
                          << std::endl;
            }
            nmaps += 5;
        } else {
            if (rank == 0)
                std::cout << "Index: " << nmaps << ", component: '" << map_components[i]
                          << "'" << std::endl;
            if (map_components[i] == "dr")
                index_components[i] = INDEX_DR;
            else if (map_components[i] == "dl")
                index_components[i] = INDEX_DL;
            else if (map_components[i] == "dt")
                index_components[i] = INDEX_DT;
            else if (map_components[i] == "da")
                index_components[i] = INDEX_DA;
            else if (map_components[i] == "dz")
                index_components[i] = INDEX_DZ;
            else if (map_components[i] == "ds")
                index_components[i] = INDEX_DS;
            else if (map_components[i] == "isw")
                index_components[i] = INDEX_ISW;
            else if (map_components[i] == "dens")
                index_components[i] = INDEX_DENSITY;
            else if (map_components[i] == "dens_max")
                index_components[i] = INDEX_DENSITY_MAX;
            else if (map_components[i] == "steps")
                index_components[i] = INDEX_NSTEPS;
            else if (map_components[i] == "phi")
                index_components[i] = INDEX_POTENTIAL;

            else {
                std::cout << "# WARNING: Need good name for map_components"
                          << std::endl;
                std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                          << std::endl;
                std::terminate();
            }
        }
        nmaps++;
    }
    // Create vector of Healpix maps
    std::vector<Healpix_Map<float>> map(parameters.nb_z_maps * nmaps);
    // Initialise maps
    for (uint i = 0; i < parameters.nb_z_maps * nmaps; i++) {
        map[i].SetNside((int)parameters.nside, (Healpix_Ordering_Scheme)0);
        map[i].fill(0.);
    }

#ifdef VERBOSE
    if (rank == 0)
        std::cout << "# Pixels attribution" << std::endl;
#endif
    std::vector<long> pixel;
    std::vector<uint> ntrajectories;
    // Assign pixels to cones
    for (uint iconerank = zero; iconerank < parameters.ncones; ++iconerank) {
        if (iconerank % static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
            Hmaps::getPixels_per_cone2(parameters, npix, pixel, ntrajectories,
                                       iconerank, coneIfRot);
        }
    }
    pixel.shrink_to_fit();
#ifdef VERBOSE
    std::cout << "# Pixels distributed : rank " << rank << " takes care of "
              << pixel.size() << " pixels" << std::endl;
#endif
    point observer;
    observer[0] = 0;
    observer[1] = 0;
    observer[2] = 0;
    std::vector<real> interpRefvec(parameters.nb_z_maps);
    std::vector<double> thomo(parameters.nb_z_maps), rhomo(parameters.nb_z_maps),
        lambdahomo(parameters.nb_z_maps), redshifthomo(parameters.nb_z_maps),
        ahomo(parameters.nb_z_maps);

    // For each redshift, estimate the homogeneous quantities [scale factor,
    // redshift, comoving distance, conformal time, affine parameter]
    for (unsigned int i = 0; i < parameters.nb_z_maps; i++) {
        photon.redshift() = z_stop_vec[i];
        const unsigned long int marked = std::distance(
            std::begin(reference),
            std::upper_bound(std::begin(reference), std::end(reference), photon,
                             [](const Photon<double, 3> &first,
                                const Photon<double, 3> &second) {
                                 return first.redshift() < second.redshift();
                             }));
        const unsigned long int firstid = marked - (marked > 0);
        const double f =
            (reference[firstid + 1].redshift() - photon.redshift()) /
            (reference[firstid + 1].redshift() - reference[firstid].redshift());
        thomo[i] =
            reference[firstid].t() * f + reference[firstid + 1].t() * (1 - f);
        rhomo[i] =
            reference[firstid].chi() * f + reference[firstid + 1].chi() * (1 - f);
        lambdahomo[i] = reference[firstid].lambda() * f +
                        reference[firstid + 1].lambda() * (1 - f);
        ahomo[i] = 1. / (1. + z_stop_vec[i]);
        redshifthomo[i] = z_stop_vec[i];
        if (parameters.stop_ray == "t") {
            interpRefvec = thomo;
        } else if (parameters.stop_ray == "r") {
            interpRefvec = rhomo;
        } else if (parameters.stop_ray == "lambda") {
            interpRefvec = lambdahomo;
        } else if (parameters.stop_ray == "a") {
            interpRefvec = ahomo;
        } else if (parameters.stop_ray == "redshift") {
            interpRefvec = z_stop_vec;
        } else if (parameters.stop_ray == "s") {
            interpRefvec = rhomo;
        } else {
            std::cout << "Need good interpolation 't, r, lambda, a, redshift or ell'"
                      << std::endl;
            std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                      << std::endl;
            std::terminate();
        }

        if (rank == 0) {
            std::cout << std::setprecision(17)
                      << "# Homogeneous quantities : z = " << redshifthomo[i]
                      << " a = " << ahomo[i] << " t = " << thomo[i]
                      << " chi = " << rhomo[i] << " lambda = " << lambdahomo[i]
                      << std::endl;
        }
    }

    reference.fullclear();

    // Fill map
#ifndef VELOCITYFIELD
    // No velocity field calculation
    uint firsttrajectory(0);
    long itraj = 0;
    for (uint icone = zero; icone < parameters.ncones; ++icone) {
        if (icone % static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
#ifdef VERBOSE
            std::cout << "# Rank " << rank << " using cone " << icone << std::endl;
#endif
            // Load octree from binary file
            Miscellaneous::loadOctree(icone, octree, conefile);
#ifdef VERBOSE
            std::cout << "# Rank " << rank << " using cone " << icone
                      << ", now computing Healpix maps" << std::endl;
#endif
            // Fill Healpix maps with ray-tracing routines
            Hmaps::FillMap(parameters, map_components, index_components,
                           ntrajectories[itraj], firsttrajectory, octree, vobs0, map,
                           nmaps, pixel, cosmology, observer, length, interpRefvec,
                           rhomo, thomo, lambdahomo, redshifthomo, ahomo);
            firsttrajectory += ntrajectories[itraj];
            itraj++;
        }
    }
#else
    // Account for velocity field
#ifdef VERBOSE
    if (rank == 0)
        std::cout
            << "# Maps at inhomogeneous redshift(s). Velocity field computed with "
            << parameters.velocity_field << std::endl;
#endif
    const point vobs = {parameters.v0x, parameters.v0y, parameters.v0z};
    std::size_t found;
    std::vector<std::string> filelistprior, shellList;

    // Get all filenames in particle directory
    Miscellaneous::getFilesinDir(parameters.partdir, filelistprior);
    for (uint ifiling = 0; ifiling < filelistprior.size(); ++ifiling) {
        // Only keep files with 'shell' in their name
        found = filelistprior[ifiling].find("shell");
        if (found != std::string::npos) {
            shellList.push_back(parameters.partdir + filelistprior[ifiling]);
        }
    }
    Miscellaneous::fullclear_vector(filelistprior);
    // Shuffle filenames so that different MPI tasks read different files
    std::shuffle(std::begin(shellList), std::end(shellList), engine1);

    uint firsttrajectory(0);
    long itraj = 0;
    real thetay(0), thetaz(0);
    std::array<std::array<double, 3>, 3> rotm1 = {{zero}};
    // Get thetay, thetaz and rotm1
    if (!parameters.isfullsky) {
        Miscellaneous::TicketizeFunction(
            rank, ntasks, [=, &parameter, &rotm1, &thetay, &thetaz] {
                Miscellaneous::get_narrow_specs(parameters, rotm1, thetay, thetaz);
            });
    }
    // Create octree from particles
    for (uint icone = zero; icone < parameters.ncones; ++icone) {
        if (icone % static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
#ifdef VERBOSE
            std::cout << "# Rank " << rank << " using cone " << icone
                      << ", getting particles..." << std::endl;
#endif
            // Load octree from binary file
            Miscellaneous::loadOctree(icone, octree, conefile);
            // WARNING ! Set rho (density) to zero in each cell of the octree;
            Utility::parallelize(octree.size(), [=, &octree](const uint i) {
                std::get<1>(octree[i]).rho() = 0;
            });
            // Read Particle files to compute the velocity field in each AMR cell
            for (uint ifile = 0; ifile < shellList.size(); ifile++) {
                std::vector<float> pos_part, vel_part;
                // Miscellaneous::TicketizeFunction(rank, ntasks, [=, &shellList,
                // &pos_part, &vel_part]{ // Uncomment in case of IO problems
                //  Get particles from files
                Hmaps::fill_particles_vectors(parameters, rotm1, cone[icone],
                                              shellList[ifile], pos_part, vel_part,
                                              z_stop_vec, thetay, thetaz);
                //});
                if (pos_part.size() > 0) {
#ifdef VERBOSE
                    std::cout << "# Rank " << rank << " using cone " << icone
                              << ", file : " << ifile << "/" << shellList.size()
                              << " npart " << pos_part.size() / 3 << " : "
                              << shellList[ifile] << std::endl;
#endif
                    // From particles, compute the velocity field with CIC or TSC scheme
                    if (parameters.velocity_field == "cic") {
                        Hmaps::CreateOctreeVelocityWithCIC(
                            octree, pos_part,
                            vel_part); // octree with velocity field using CIC interpolation
#ifdef VERBOSE
                        std::cout << "# Rank " << rank << " using cone " << icone
                                  << ", Octree created with CIC" << std::endl;
#endif
                    } else if (parameters.velocity_field == "tsc") {
                        Hmaps::CreateOctreeVelocityWithTSC(
                            octree, pos_part,
                            vel_part); // octree with velocity field using TSC interpolation
#ifdef VERBOSE
                        std::cout << "# Rank " << rank << " using cone " << icone
                                  << ", Octree created with TSC" << std::endl;
#endif
                    } else {
                        std::cout << "No velocity field computed, please type cic or tsc"
                                  << std::endl;
                        std::cout << "# Error at file " << __FILE__
                                  << ", line : " << __LINE__ << std::endl;
                        std::terminate();
                    }
                }
            }
            // Finalize
            const unsigned int lvlmax =
                (std::get<0>(
                     *std::max_element(std::begin(octree), std::end(octree),
                                       [](const element &x, const element &y) {
                                           return std::get<0>(x).level() <
                                                  std::get<0>(y).level();
                                       }))
                     .level());
            const unsigned int lvlmin =
                (std::get<0>(
                     *std::min_element(std::begin(octree), std::end(octree),
                                       [](const element &x, const element &y) {
                                           return std::get<0>(x).level() <
                                                  std::get<0>(y).level();
                                       }))
                     .level());
            // Normalise velocity field with mass
            Utility::parallelize(octree.size(), [=, &octree](const uint i) {
                double mass = std::get<1>(octree[i]).rho();
                mass = mass + (mass == 0);
                std::get<1>(octree[i]).vx() /= mass;
                std::get<1>(octree[i]).vy() /= mass;
                std::get<1>(octree[i]).vz() /= mass;
            });
            // If the density is zero in one cell, get value from parent
            for (unsigned int ilvl = lvlmin + 1; ilvl <= lvlmax; ilvl++) {
                Utility::parallelize(octree.size(), [=, &octree](const uint i) {
                    if (std::get<0>(octree[i]).level() == ilvl) {
                        Gravity<floating, 3> data;
                        if (!std::isnormal(std::get<1>(octree[i]).rho())) {
                            data = std::get<1>(*octree.find(std::get<0>(octree[i]).parent()));
                            std::get<1>(octree[i]).vxyz() = {data.vx(), data.vy(), data.vz()};
                        }
                    }
                });
            }

            // Fill map
#ifdef VERBOSE
            std::cout << "# Rank " << rank << " using cone " << icone
                      << ", now computing Healpix maps, ntrajectory = "
                      << ntrajectories[itraj] << std::endl;
#endif
            // Fill Healpix maps with ray-tracing routines
            Hmaps::FillMap(parameters, map_components, index_components,
                           ntrajectories[itraj], firsttrajectory, octree, vobs, map,
                           nmaps, pixel, cosmology, observer, length, interpRefvec,
                           rhomo, thomo, lambdahomo, redshifthomo, ahomo);
#ifdef VERBOSE
            std::cout << "# Rank " << rank << " using cone " << icone
                      << ", Healpix maps computed" << std::endl;
#endif
            firsttrajectory += ntrajectories[itraj];
            itraj++;
        }
    }
#endif

#ifdef VERBOSE
    std::cout
        << "# Rank " << rank
        << " : Healpix maps filled, now waiting for other procs to finalize "
        << std::endl;
#endif
    Miscellaneous::fullclear_vector(pixel);
    octree.fullclear();
    // Wait for all tasks
    MPI_Barrier(MPI_COMM_WORLD);

    // Finalisation
    float *mapdata1;
    float *mapdata2;
    mapdata1 = (float *)malloc(npix * sizeof(float));
    mapdata2 = (float *)malloc(npix * sizeof(float));
    const char coordsys = 'G';
    // Create vector of Healpix maps
    Healpix_Map<float> map_tmp;
    map_tmp.SetNside((int)parameters.nside, (Healpix_Ordering_Scheme)0);
    int my_iz(-1), my_i(-1);

    // Loop over all the maps
    for (uint iz = 0; iz < parameters.nb_z_maps; iz++) {
        for (uint i = 0; i < nmaps; i++) {
            // Copy map in temporary array
            map[nmaps * iz + i].Map().copyToPtr<float>(mapdata1);
            MPI_Barrier(MPI_COMM_WORLD);
            // Add contribution from every MPI task mapdata1 into the array mapdata2
            const int icone = (nmaps * iz + i) % static_cast<uint>(ntasks);
            MPI_Reduce((void *)mapdata1, (void *)mapdata2, (int)npix, MPI_FLOAT,
                       MPI_SUM, icone, MPI_COMM_WORLD);
            if (rank == icone) {
                my_iz = iz;
                my_i = i;
            }
            // Write the outputs if all the tasks contain a full map, or if the total
            // number of maps are contained std::cout<<"rank = "<<rank<<" map
            // "<<nmaps*iz + i<<" ntasks "<<ntasks<<" ratio "<<(nmaps*iz + i) %
            // static_cast<uint>(ntasks)<<" max
            // "<<static_cast<uint>(nmaps*parameters.nb_z_maps - 1)<<" my "<<my_iz<<"
            // "<<my_i<<std::endl;
            if (((nmaps * iz + i) % static_cast<uint>(ntasks) ==
                 static_cast<uint>(ntasks - 1)) |
                ((nmaps * iz + i) ==
                 static_cast<uint>(nmaps * parameters.nb_z_maps - 1))) {
                // Select tasks
                if (my_iz != -1) {
                    // Clean the map
                    for (long pix = 0; pix < npix; pix++) {
                        // If pixel is not normal (i.e. zero, infinity or NaN), then put
                        // undefined value
                        if (!std::isnormal(mapdata2[pix])) {
                            mapdata2[pix] = Healpix_undef;
                        }
                    }
#ifdef VERBOSE
                    std::cout << " #-- Writing output with rank = " << rank
                              << " redshift = " << z_stop_vec[my_iz] << " id = " << my_i
                              << std::endl;
#endif
                    std::string outputfilename;
                    outputfilename = parameters.outputdir + parameters.outputprefix +
                                     "_nside" + std::to_string(parameters.nside) + "_" +
                                     std::to_string(my_i) + "_z_" +
                                     std::to_string(z_stop_vec[my_iz]) + "_" +
                                     parameters.stop_ray + ".fits";
                    // Write Healpix map
                    write_healpix_map(mapdata2, parameters.nside, outputfilename.data(),
                                      0, &coordsys);
                    my_iz = -1;
                }
            }
        }
    }

    free(mapdata1);
    free(mapdata2);

    // Finalization
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    if (rank == 0)
        std::cout << "# Run completed !" << std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
