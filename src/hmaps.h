/* ********************************** HMAPS **********************************
 */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Hmaps
// DESCRIPTION :    Some hmaps function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           hmaps.h
/// \brief          Some hmaps function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
/// (michel-andres.breton@obspm.fr) \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef HMAPS_H_INCLUDED
#define HMAPS_H_INCLUDED

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

#ifdef GCCBELOW7
#include <experimental/algorithm>
#endif
// Include libs
#include <mpi.h>
// Include project
#include "TReadHDF5.h"
#include "cone.h"
#include "input.h"
#include "lensing.h"
#include "magrathea/abstracthypersphere.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/hypercube.h"
#include "magrathea/hypersphere.h"
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/timer.h"
#include "output.h"
#include "utility.h"
// Include HEALPIX
#include <chealpix.h>
#include <healpix_map.h>
#include <pointing.h>

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
#ifdef VELOCITYFIELD
    real v0x;
    real v0y;
    real v0z;
    std::string velocity_field;
    std::string partdir;
#endif

    std::string beam;
    std::string map_components;
    uint microcoeff;
    uint nb_z_maps;
    uint nbundlemin;
    long nside;
    uint nsteps;
    real openingmin;
    std::string outputdir;
    std::string outputprefix;
    std::string plane;
    std::string stop_ray;
    std::string stop_bundle;
    real z_stop_min;
    real z_stop_max;

} parameters;

class Hmaps {

    // Methodes
public:
    // Read parameter file
    template <class Parameters, class Map>
    static void ReadParamFile(Parameters &parameters, Map &parameter);

    // Target
    template <class Parameter, typename Integer, class Cones>
    static void
    getPixels_per_cone(const Parameter &parameters, const Integer npix,
                       std::vector<Integer> &pixel,
                       std::vector<unsigned int> &ntrajectories,
                       const unsigned int iconerank, const Cones &cones);
    template <class Parameter, typename Integer, class Cones>
    static void
    getPixels_per_cone2(const Parameter &parameters, const Integer npix,
                        std::vector<Integer> &pixel,
                        std::vector<unsigned int> &ntrajectories,
                        const unsigned int iconerank, const Cones &cones);

    // Fill particles
#ifdef VELOCITYFIELD
    template <class Parameter, class Cone, typename Type1, typename Type2>
    static void fill_particles_vectors(
        const Parameter &parameters,
        const std::array<std::array<double, 3>, 3> &rotm1, const Cone &cone,
        const std::vector<std::string> &shellList, std::vector<Type1> &pos_part,
        std::vector<Type1> &vel_part, const std::vector<Type2> &z_stop_vec,
        const double thetay, const double thetaz);
    template <class Parameter, class Cone, typename Type1, typename Type2>
    static void fill_particles_vectors(
        const Parameter &parameters,
        const std::array<std::array<double, 3>, 3> &rotm1, const Cone &cone,
        const std::string shellList, std::vector<Type1> &pos_part,
        std::vector<Type1> &vel_part, const std::vector<Type2> &z_stop_vec,
        const double thetay, const double thetaz);
    // Velocity field octree
    template <
        typename Type1,
        template <typename Type, class Index, class Data, unsigned int Dimension,
                  class Position, class Extent, class Element, class Container>
        class Octree,
        typename Type, class Index, class Data, unsigned int Dimension,
        class Position, class Extent, class Element, class Container>
    static void CreateOctreeVelocityWithCIC(
        Octree<Type, Index, Data, Dimension, Position, Extent, Element, Container>
            &octree,
        const std::vector<Type1> &pos_part, const std::vector<Type1> &vel_part);
    template <
        typename Type1,
        template <typename Type, class Index, class Data, unsigned int Dimension,
                  class Position, class Extent, class Element, class Container>
        class Octree,
        typename Type, class Index, class Data, unsigned int Dimension,
        class Position, class Extent, class Element, class Container>
    static void CreateOctreeVelocityWithTSC(
        Octree<Type, Index, Data, Dimension, Position, Extent, Element, Container>
            &octree,
        const std::vector<Type1> &pos_part, const std::vector<Type1> &vel_part);
#endif
    // Fill maps
    template <class Parameter, class Octree, class Map, class Pixel,
              class Cosmology, class Point, typename Integer, typename Real>
    static void
    FillMapPropagate(const Parameter &parameters, const Integer ntrajectory,
                     const Integer firsttrajectory, const Octree &octree,
                     const Point &vobs, Map &map, const Integer nmaps,
                     const Pixel &pixel, const Cosmology &cosmology,
                     const point &observer, const Real length,
                     const std::vector<Real> &z_stop_vec);

    template <class Parameter, class Octree, class Map, class Pixel,
              class Cosmology, class Point, typename Integer, typename Real>
    static void
    FillMap(const Parameter &parameters,
            const std::vector<std::string> &map_components,
            const std::vector<unsigned int> &index_components,
            const Integer ntrajectory, const Integer firsttrajectory,
            const Octree &octree, const Point &vobs, Map &map,
            const unsigned int nmaps, const Pixel &pixel,
            const Cosmology &cosmology, const point &observer, const Real length,
            const std::vector<Real> &interpRefvec, const std::vector<Real> &rhomo,
            const std::vector<Real> &thomo, const std::vector<Real> &lambdahomo,
            const std::vector<Real> &redshifthomo,
            const std::vector<Real> &ahomo);
};

// Read parameter file
/// \brief          Read parameter file.
/// \details        Read and put in a structure the parameters.
/// \tparam         Parameters structure type
/// \tparam         Map map type
/// \param[in,out]  parameters Structure containing the parameters.
/// \param[in]      parameter Contains parameters to be rewritten
/// \return         Filled parameters structure.
template <class Parameters, class Map>
void Hmaps::ReadParamFile(Parameters &parameters, Map &parameter) {

    parameters.seed = std::stoul(parameter["seed"]);
    parameters.mpc = std::stod(parameter["mpc"]);
    parameters.rhoch2 = std::stod(parameter["rhoch2"]);
    parameters.typefile = std::stoul(parameter["typefile"]);
    parameters.isfullsky = std::stoul(parameter["isfullsky"]);
    parameters.paramfile = parameter["paramfile"];
    parameters.evolfile = parameter["evolfile"];
    parameters.celldir = parameter["celldir"];
    parameters.conedir = parameter["conedir"];
    parameters.conefmt = parameter["conefmt"];
    parameters.openingmin = std::stod(parameter["openingmin"]);
    parameters.nside = (long)std::stoul(parameter["nside"]);
    parameters.outputdir = parameter["outputdir"];
    parameters.z_stop_min = std::stod(parameter["z_stop_min"]);
    parameters.z_stop_max = std::stod(parameter["z_stop_max"]);
    parameters.nb_z_maps = std::stoul(parameter["nb_z_maps"]);
    parameters.ncoarse = std::stoul(parameter["ncoarse"]);
    parameters.ncones = std::stoul(parameter["ncones"]);
    parameters.buffer = std::stod(parameter["buffer"]);
    parameters.stop_ray = parameter["stop_ray"];
    parameters.stop_bundle = parameter["stop_bundle"];
    parameters.plane = parameter["plane"];
    parameters.beam = parameter["beam"];
    parameters.outputprefix = parameter["outputprefix"];
    parameters.nbundlemin = std::stoul(parameter["nbundlemin"]);
    parameters.nsteps = std::stoul(parameter["nsteps"]);
    parameters.microcoeff = std::stoul(parameter["microcoeff"]);
    parameters.map_components = parameter["map_components"];
#ifdef VELOCITYFIELD
    parameters.v0x = std::stod(parameter["v0x"]);
    parameters.v0y = std::stod(parameter["v0y"]);
    parameters.v0z = std::stod(parameter["v0z"]);
    parameters.velocity_field = parameter["velocity_field"];
    parameters.partdir = parameter["partdir"];
#endif
}

/// \brief          Get targets.
/// \details        Get targets inside cone.
/// \tparam         Parameter Parameter type
/// \tparam         Integer Integer type
/// \tparam         Cones Cones type
/// \param[in]      parameters Parameter structure.
/// \param[in]      npix Number of pixels.
/// \param[in,out]  pixel Vector of pixels inside cones that are treated by the
/// same MPI task. \param[in]	    ntrajectories Number of pixels per cone.
/// \param[in]	    iconerank Number of the cone of interest.
/// \param[in]	    Cones Geometry of all the cones.
/// return  	    Filled vector of targets inside cones for a given MPI task
template <class Parameter, typename Integer, class Cones>
void Hmaps::getPixels_per_cone(const Parameter &parameters, const Integer npix,
                               std::vector<Integer> &pixel,
                               std::vector<unsigned int> &ntrajectories,
                               const unsigned int iconerank,
                               const Cones &cones) {
    // Fill temporary pixel vector with -1
    std::vector<long> pixeltmp(npix, -1);
    ntrajectories.push_back(0);
    // Loop over all the pixels
    Utility::parallelize(npix, [=, &parameters, &iconerank, &cones,
                                &ntrajectories, &pixeltmp](long pix) {
        double ang1(0), ang2(0), vec[3];
        bool closest = false;
        // Convert pixel to 3D vector
        pix2vec_ring(parameters.nside, pix, vec);
        // Get the angle between cone base and pixel
        ang1 = cones[iconerank].direction(0) * vec[0];
        ang1 += cones[iconerank].direction(1) * vec[1];
        ang1 += cones[iconerank].direction(2) * vec[2];
        ang1 = std::acos(ang1 / cones[iconerank].length());
        // If the angle is smaller than the cone width
        if (ang1 <= cones[iconerank].angle()) {
            closest = true;
            // Loop over all the other cones to find if there is another one that is
            // closer
            for (uint icone = zero; icone < parameters.ncones; icone++) {
                if (icone != static_cast<uint>(iconerank)) {
                    // Compute angle for other cone
                    ang2 = cones[icone].direction(0) * vec[0];
                    ang2 += cones[icone].direction(1) * vec[1];
                    ang2 += cones[icone].direction(2) * vec[2];
                    ang2 = std::acos(ang2 / cones[icone].length());
                    // If other cone is closer, then the pixel is not attributed
                    if (ang2 < ang1 && ang2 < cones[icone].angle()) {
                        closest = false;
                        break;
                    }
                }
            }
        }
        if (closest) {
            pixeltmp[pix] = pix;
        }
    });
    // Erase pixels which are not inside the cone
    pixeltmp.erase(std::remove(std::begin(pixeltmp), std::end(pixeltmp), -1),
                   std::end(pixeltmp));
    // Count number of pixels inside the cone
    ntrajectories.back() = pixeltmp.size();
    // Insert pixels in the vectors for a single MPI task which will take care of
    // several cones
    pixel.insert(std::end(pixel), std::begin(pixeltmp), std::end(pixeltmp));
}

/// \brief          Get targets.
/// \details        Get targets inside cone.
/// \tparam         Parameter Parameter type
/// \tparam         Integer Integer type
/// \tparam         Cones Cones type
/// \param[in]      parameters Parameter structure.
/// \param[in]      npix Number of pixels.
/// \param[in,out]  pixel Vector of pixels inside cones that are treated by the
/// same MPI task. \param[in]	    ntrajectories Number of pixels per cone.
/// \param[in]	    iconerank Number of the cone of interest.
/// \param[in]	    Cones Geometry of all the cones.
/// return  	    Filled vector of targets inside cones for a given MPI task
template <class Parameter, typename Integer, class Cones>
void Hmaps::getPixels_per_cone2(const Parameter &parameters, const Integer npix,
                                std::vector<Integer> &pixel,
                                std::vector<unsigned int> &ntrajectories,
                                const unsigned int iconerank,
                                const Cones &cones) {
    // Fill temporary pixel vector with -1
    std::vector<long> pixeltmp(npix, -1);
    ntrajectories.push_back(0);
    const double cone_length = cones[0].length();
    // Loop over all the pixels
    Utility::parallelize(npix, [=, &parameters, &iconerank, &cones,
                                &ntrajectories, &pixeltmp,
                                &cone_length](long pix) {
        double vec[3];
        // Convert pixel to 3D vector
        pix2vec_ring(parameters.nside, pix, vec);
        vec[0] *= 0.9 * cone_length;
        vec[1] *= 0.9 * cone_length;
        vec[2] *= 0.9 * cone_length;
        double reference(0), length(0), distance(0);
        bool ok = false;
        reference = 0;
        length = 0;
        // Check if pixel is inside the cone
        if (cones[iconerank].inside(vec)) {
            ok = true;
            // Compute scalar product of cone base and pixel direction
            for (unsigned int idim = 0; idim < 3; ++idim) {
                length +=
                    (cones[iconerank].base(idim) - cones[iconerank].vertex(idim)) *
                    (vec[idim] - cones[iconerank].vertex(idim));
            }
            length /= cones[iconerank].template pow<2>(cones[iconerank].length());
            // Compute the ditance between cone base and pixel
            for (unsigned int idim = 0; idim < 3; ++idim) {
                reference += cones[iconerank].template pow<2>(
                    vec[idim] -
                    (cones[iconerank].vertex(idim) +
                     (cones[iconerank].base(idim) - cones[iconerank].vertex(idim)) *
                         length));
            }
            // Loop over all the other cones
            for (unsigned int icone = 0; icone < cones.size(); ++icone) {
                length = 0;
                distance = 0;
                if (icone != iconerank) {
                    // Compute scalar product of other cone base and pixel direction
                    for (unsigned int idim = 0; idim < 3; ++idim) {
                        length +=
                            (cones[icone].base(idim) - cones[iconerank].vertex(idim)) *
                            (vec[idim] - cones[icone].vertex(idim));
                    }
                    if (!(length < 0)) {
                        length /= cones[icone].template pow<2>(cones[icone].length());
                        // Compute the ditance between other cone base and pixel
                        for (unsigned int idim = 0; idim < 3; ++idim) {
                            distance += cones[icone].template pow<2>(
                                vec[idim] -
                                (cones[icone].vertex(idim) +
                                 (cones[icone].base(idim) - cones[icone].vertex(idim)) *
                                     length));
                        }
                        // Check if other cone is closer to pixel
                        if (distance < reference) {
                            ok = false;
                            icone = cones.size();
                        }
                    } //  if length
                }     // if
            }         // for icone
            if (ok) {
                pixeltmp[pix] = pix;
            }
        } //  if
    });

    // Erase pixels which are not inside the cone
    pixeltmp.erase(std::remove(std::begin(pixeltmp), std::end(pixeltmp), -1),
                   std::end(pixeltmp));
    // Count number of pixels inside the cone
    ntrajectories.back() = pixeltmp.size();
    // Insert pixels in the vectors for a single MPI task which will take care of
    // several cones
    pixel.insert(std::end(pixel), std::begin(pixeltmp), std::end(pixeltmp));
}

#ifdef VELOCITYFIELD
// Fill vectors of particles
/// \brief          Rewrite particle in vector.
/// \details        Rewrite particle position/velocity from HDF5 in a vector.
/// \tparam         Parameter Parameter type
/// \tparam         Cone cone type
/// \tparam         Type1 type1 type
/// \tparam         Type2 type2 type
/// \param[in]      parameters Parameter structure.
/// \param[in]      rotm1 Rotation matrix (used for narrow cones).
/// \param[in]      cone Cone parameters
/// \param[in]      shellList List of HDF5 files containing particles
/// \param[in,out]  pos_part Position of particles
/// \param[in,out]  vel_part Velocity of particles
/// \param[in]      z_stop_vec Redshift at which we want to compute Healpix maps
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
template <class Parameter, class Cone, typename Type1, typename Type2>
void Hmaps::fill_particles_vectors(
    const Parameter &parameters,
    const std::array<std::array<double, 3>, 3> &rotm1, const Cone &cone,
    const std::vector<std::string> &shellList, std::vector<Type1> &pos_part,
    std::vector<Type1> &vel_part, const std::vector<Type2> &z_stop_vec,
    const double thetay, const double thetaz) {

    unsigned long long int marker(0);
    const double c = magrathea::Constants<double>::c();

#ifdef VERBOSE
    std::cout << "# Now loading the position and velocity of particles from the "
                 "following files : "
              << std::endl;
#endif
    Miscellaneous::fullclear_vector(pos_part);
    Miscellaneous::fullclear_vector(vel_part);

    // Loop over all the particle files
    for (uint ifiling = 0; ifiling < shellList.size(); ifiling++) {
        double amaxing(0), amining(0);
        // Get the limits of particle shell
        TReadHDF5::getAttribute(shellList[ifiling], "metadata/cone_info", "amax",
                                amaxing);
        TReadHDF5::getAttribute(shellList[ifiling], "metadata/cone_info", "amin",
                                amining);
        double da = amaxing - amining;
        // Need to define a buffer zone around the reference redshift.
        // Must be chosen so the maximum observed redshift is far enough (take dz =
        // 5e-3. Compare to redshift step between shells and put
        // int(dz/redshiftstep) + 1)
        amaxing += 4 * da;
        amining -= 4 * da;
#ifdef VERBOSE
        std::cout << shellList[ifiling] << std::endl;
#endif
        // Loop over all the reference redshifts
        for (uint iz = 0; iz < z_stop_vec.size(); iz++) {
            double scale_factor = 1. / (z_stop_vec[iz] + 1);
            // If file is inside the buffer zone of any reference redshift, then get
            // the particles
            if (scale_factor <= amining && scale_factor >= amaxing) {
                TReadHDF5::fillVectors_part(parameters, shellList[ifiling], thetay,
                                            thetaz, cone, "data", "position_part",
                                            pos_part, "velocity_part", vel_part);
                float unit_l(0), unit_t(0);
                TReadHDF5::getAttribute(shellList[ifiling], "metadata/ramses_info",
                                        "unit_l", unit_l); // comobile
                TReadHDF5::getAttribute(shellList[ifiling], "metadata/ramses_info",
                                        "unit_t", unit_t); // superconformal unit
                double factor = unit_l * 1e-2 / (c * unit_t);
                std::transform(
                    std::begin(vel_part) + marker, std::end(vel_part),
                    std::begin(vel_part) + marker,
                    std::bind1st(std::multiplies<double>(), factor)); // SI units
                marker = vel_part.size();
                break;
            }
        }
    }
    // For narrow cones, need to rotate position and velocity of particles
    if (!parameters.isfullsky) {
        std::vector<Type1> pos_part_tmp = pos_part;
        std::vector<Type1> vel_part_tmp = vel_part;
        Utility::parallelize(pos_part.size() / 3, [=, &pos_part, &vel_part,
                                                   &pos_part_tmp, &vel_part_tmp,
                                                   &rotm1](const uint i) {
            pos_part[3 * i] = pos_part_tmp[3 * i] * rotm1[0][0] +
                              pos_part_tmp[3 * i + 1] * rotm1[0][1] +
                              pos_part_tmp[3 * i + 2] * rotm1[0][2];
            pos_part[3 * i + 1] = pos_part_tmp[3 * i] * rotm1[1][0] +
                                  pos_part_tmp[3 * i + 1] * rotm1[1][1] +
                                  pos_part_tmp[3 * i + 2] * rotm1[1][2];
            pos_part[3 * i + 2] = pos_part_tmp[3 * i] * rotm1[2][0] +
                                  pos_part_tmp[3 * i + 1] * rotm1[2][1] +
                                  pos_part_tmp[3 * i + 2] * rotm1[2][2];
            vel_part[3 * i] = vel_part_tmp[3 * i] * rotm1[0][0] +
                              vel_part_tmp[3 * i + 1] * rotm1[0][1] +
                              vel_part_tmp[3 * i + 2] * rotm1[0][2];
            vel_part[3 * i + 1] = vel_part_tmp[3 * i] * rotm1[1][0] +
                                  vel_part_tmp[3 * i + 1] * rotm1[1][1] +
                                  vel_part_tmp[3 * i + 2] * rotm1[1][2];
            vel_part[3 * i + 2] = vel_part_tmp[3 * i] * rotm1[2][0] +
                                  vel_part_tmp[3 * i + 1] * rotm1[2][1] +
                                  vel_part_tmp[3 * i + 2] * rotm1[2][2];
        });
    }
}
// Fill vectors of particles
/// \brief          Rewrite particle in vector.
/// \details        Rewrite particle position/velocity from HDF5 in a vector.
/// Single file version \tparam         Parametr Parameter type \tparam Cone
/// cone type \tparam         Type1 type1 type \tparam         Type2 type2 type
/// \param[in]      parameters Parameter structure.
/// \param[in]      rotm1 Rotation matrix (used for narrow cones).
/// \param[in]      cone Cone parameters
/// \param[in]      shellList List of HDF5 files containing particles
/// \param[in,out]  pos_part Position of particles
/// \param[in,out]  vel_part Velocity of particles
/// \param[in]      z_stop_vec Redshift at which we want to compute Healpix maps
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
template <class Parameter, class Cone, typename Type1, typename Type2>
void Hmaps::fill_particles_vectors(
    const Parameter &parameters,
    const std::array<std::array<double, 3>, 3> &rotm1, const Cone &cone,
    const std::string shellList, std::vector<Type1> &pos_part,
    std::vector<Type1> &vel_part, const std::vector<Type2> &z_stop_vec,
    const double thetay, const double thetaz) {
    std::vector<std::string> shellLists{shellList};
    Hmaps::fill_particles_vectors(parameters, rotm1, cone, shellLists, pos_part,
                                  vel_part, z_stop_vec, thetay, thetaz);
}

// Compute CIC velocity field
/// \brief          Add velocity to Octree using CIC
/// \details        Add velocity to Octree using CIC
/// \tparam         Type1 type type
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in,out]  octree Octree be to filled with velocity field
/// \param[in]      pos_part Position of particles
/// \param[in]      vel_part Velocity of particles
template <
    typename Type1,
    template <typename Type, class Index, class Data, unsigned int Dimension,
              class Position, class Extent, class Element, class Container>
    class Octree,
    typename Type, class Index, class Data, unsigned int Dimension,
    class Position, class Extent, class Element, class Container>
void Hmaps::CreateOctreeVelocityWithCIC(
    Octree<Type, Index, Data, Dimension, Position, Extent, Element, Container>
        &octree,
    const std::vector<Type1> &pos_part, const std::vector<Type1> &vel_part) {

    const unsigned int lvlmax =
        (std::get<0>(*std::max_element(std::begin(octree), std::end(octree),
                                       [](const Element &x, const Element &y) {
                                           return std::get<0>(x).level() <
                                                  std::get<0>(y).level();
                                       }))
             .level());
    const unsigned int lvlmin =
        (std::get<0>(*std::min_element(std::begin(octree), std::end(octree),
                                       [](const Element &x, const Element &y) {
                                           return std::get<0>(x).level() <
                                                  std::get<0>(y).level();
                                       }))
             .level());

    // Loop over all the levels
    for (uint ilvl = lvlmin; ilvl <= lvlmax; ilvl++) {
        const double invextension = 1 / (EXTENT / pow(2, ilvl));
        const double half = 0.5 * EXTENT / pow(2, ilvl);
#ifdef VERBOSE
        std::cout << "# Compute velocity field at level " << ilvl << std::endl;
#endif
        // Loop over all the particles
        Utility::parallelize(pos_part.size() / 3, [=, &octree, &pos_part, &vel_part,
                                                   &invextension,
                                                   &half](const uint i) {
            Index idxvertex;
            Data data;
            unsigned long long int marker(0);
            double vratio(0);
            // Loop over the 8 neighboring cells of a given particle
            for (int iz = -1; iz <= 1; iz += 2) {
                for (int iy = -1; iy <= 1; iy += 2) {
                    for (int ix = -1; ix <= 1; ix += 2) {
                        // Create an index at the level of interest for the neighboring cell
                        idxvertex = idxvertex.template compute<Type, Position, Extent>(
                            ilvl, pos_part[3 * i] + half * ix,
                            pos_part[3 * i + 1] + half * iy,
                            pos_part[3 * i + 2] + half * iz);
                        // Given an index in the octree that is consistent with the created
                        // index
                        marker = std::distance(
                            std::begin(octree),
                            std::upper_bound(
                                std::begin(octree), std::end(octree),
                                Element(idxvertex, data),
                                [](const Element &first, const Element &second) {
                                    return std::get<0>(first) < std::get<0>(second);
                                }));
                        // If the index exists in the octree, compute CIC
                        if (std::get<0>(*(std::begin(octree) + marker - (marker > 0))) ==
                            idxvertex) {
                            vratio =
                                (1 -
                                 std::abs(
                                     pos_part[3 * i] -
                                     idxvertex.template center<Type, Position, Extent>(0)) *
                                     invextension) *
                                (1 -
                                 std::abs(
                                     pos_part[3 * i + 1] -
                                     idxvertex.template center<Type, Position, Extent>(1)) *
                                     invextension) *
                                (1 -
                                 std::abs(
                                     pos_part[3 * i + 2] -
                                     idxvertex.template center<Type, Position, Extent>(2)) *
                                     invextension);
                            std::get<1>(octree[marker - (marker > 0)]).rho() += vratio;
                            std::get<1>(octree[marker - (marker > 0)]).vx() +=
                                vel_part[3 * i] * vratio;
                            std::get<1>(octree[marker - (marker > 0)]).vy() +=
                                vel_part[3 * i + 1] * vratio;
                            std::get<1>(octree[marker - (marker > 0)]).vz() +=
                                vel_part[3 * i + 2] * vratio;
                        }
                    }
                }
            }
        });
    }
}

// Compute TSC velocity field
/// \brief          Add velocity to Octree using TSC
/// \details        Add velocity to Octree using TSC
/// \tparam         Type1 Type1 type
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in,out]  octree Octree be to filled with velocity field
/// \param[in]      pos_part Position of particles
/// \param[in]      vel_part Velocity of particles
template <
    typename Type1,
    template <typename Type, class Index, class Data, unsigned int Dimension,
              class Position, class Extent, class Element, class Container>
    class Octree,
    typename Type, class Index, class Data, unsigned int Dimension,
    class Position, class Extent, class Element, class Container>
void Hmaps::CreateOctreeVelocityWithTSC(
    Octree<Type, Index, Data, Dimension, Position, Extent, Element, Container>
        &octree,
    const std::vector<Type1> &pos_part, const std::vector<Type1> &vel_part) {

    const unsigned int lvlmax =
        (std::get<0>(*std::max_element(std::begin(octree), std::end(octree),
                                       [](const Element &x, const Element &y) {
                                           return std::get<0>(x).level() <
                                                  std::get<0>(y).level();
                                       }))
             .level());
    const unsigned int lvlmin =
        (std::get<0>(*std::min_element(std::begin(octree), std::end(octree),
                                       [](const Element &x, const Element &y) {
                                           return std::get<0>(x).level() <
                                                  std::get<0>(y).level();
                                       }))
             .level());
    // Loop over all the levels
    for (uint ilvl = lvlmin; ilvl <= lvlmax; ilvl++) {
        const double half = 0.5 * EXTENT / pow(2, ilvl);
        const double twohalves = 2 * half;
#ifdef VERBOSE
        std::cout << "# Compute velocity field at level " << ilvl << std::endl;
#endif
        // Loop over all the particles
        Utility::parallelize(pos_part.size() / 3, [=, &octree, &pos_part, &vel_part,
                                                   &half,
                                                   &twohalves](const uint i) {
            Index idxvertex;
            Data data;
            std::array<double, Index::dimension()> dist;
            unsigned long long int marker(0);
            double vratio(0), weightx(0), weighty(0), weightz(0);
            // Loop over the 27 neighboring cells of a given particle
            for (int iz = -1; iz <= 1; iz += 1) {
                const int aiz = abs(iz);
                for (int iy = -1; iy <= 1; iy += 1) {
                    const int aiy = abs(iy);
                    for (int ix = -1; ix <= 1; ix += 1) {
                        const int aix = abs(ix);
                        // Create an index at the level of interest for the neighboring cell
                        idxvertex = idxvertex.template compute<Type, Position, Extent>(
                            ilvl, pos_part[3 * i] + twohalves * ix,
                            pos_part[3 * i + 1] + twohalves * iy,
                            pos_part[3 * i + 2] + twohalves * iz);
                        // Given an index in the octree that is consistent with the created
                        // index
                        marker = std::distance(
                            std::begin(octree),
                            std::upper_bound(
                                std::begin(octree), std::end(octree),
                                Element(idxvertex, data),
                                [](const Element &first, const Element &second) {
                                    return std::get<0>(first) < std::get<0>(second);
                                }));
                        // If the index exists in the octree, compute TSC
                        if (std::get<0>(*(std::begin(octree) + marker - (marker > 0))) ==
                            idxvertex) {
                            for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                                dist[idim] = std::abs(
                                    (std::get<0>(*(std::begin(octree) + marker - (marker > 0)))
                                         .template center<Type, Position, Extent>(idim) -
                                     pos_part[3 * i + idim]) /
                                    (twohalves));
                            }
                            weightx = aix * 0.5 * (1.5 - dist[0]) * (1.5 - dist[0]) +
                                      (1 - aix) * (0.75 - dist[0] * dist[0]);
                            weighty = aiy * 0.5 * (1.5 - dist[1]) * (1.5 - dist[1]) +
                                      (1 - aiy) * (0.75 - dist[1] * dist[1]);
                            weightz = aiz * 0.5 * (1.5 - dist[2]) * (1.5 - dist[2]) +
                                      (1 - aiz) * (0.75 - dist[2] * dist[2]);
                            vratio = weightx * weighty * weightz;
                            std::get<1>(octree[marker - (marker > 0)]).rho() += vratio;
                            std::get<1>(octree[marker - (marker > 0)]).vx() +=
                                vel_part[3 * i] * vratio;
                            std::get<1>(octree[marker - (marker > 0)]).vy() +=
                                vel_part[3 * i + 1] * vratio;
                            std::get<1>(octree[marker - (marker > 0)]).vz() +=
                                vel_part[3 * i + 2] * vratio;
                        }
                    }
                }
            }
        });
    }
}
#endif

// Fill Healpix maps
/// \brief          Produce Healpix maps
/// \details        Produce Healpix maps
/// \tparam         Parameter Parameter type
/// \tparam         Octree octree type
/// \tparam         Map map type
/// \tparam         Pixel pixel type
/// \tparam         Cosmology cosmology type
/// \tparam         Point Point type
/// \tparam         Integer integer type
/// \tparam         Real real/float type
/// \param[in]      parameters Parameter structure.
/// \param[in]      map_components Vector of strings containing the map types.
/// \param[in]      index_components Vector containing the index relative to the
/// map types. \param[in]      ntrajectory Number of pixels to be filled
/// \param[in]      firsttrajectory first index for pixel for a given cone
/// \param[in]      octree Octree
/// \param[in]      vobs Observer velocity
/// \param[in,out]  map Map
/// \param[in]      nmaps Number of Healpix maps
/// \param[in]      pixel Pixel vector
/// \param[in]      cosmology Cosmological table
/// \param[in]      observer Observer position
/// \param[in]      length R.U to S.I units for length
/// \param[in]      interpRefvec vector of interpolations at a given quantity
/// for a given surface \param[in]      rhomo vector of FLRW distance at the
/// redshifts of interest \param[in]      thomo vector of FLRW time at the
/// redshifts of interest \param[in]      lambdahomo vector of FLRW affine
/// parameter at the redshifts of interest \param[in]      redshifthomo vector
/// of FLRW redshift at the redshifts of interest \param[in]      ahomo vector
/// of FLRW scale factor at the redshifts of interest
template <class Parameter, class Octree, class Map, class Pixel,
          class Cosmology, class Point, typename Integer, typename Real>
void Hmaps::FillMap(
    const Parameter &parameters, const std::vector<std::string> &map_components,
    const std::vector<unsigned int> &index_components,
    const Integer ntrajectory, const Integer firsttrajectory,
    const Octree &octree, const Point &vobs, Map &map, const unsigned int nmaps,
    const Pixel &pixel, const Cosmology &cosmology, const point &observer,
    const Real length, const std::vector<Real> &interpRefvec,
    const std::vector<Real> &rhomo, const std::vector<Real> &thomo,
    const std::vector<Real> &lambdahomo, const std::vector<Real> &redshifthomo,
    const std::vector<Real> &ahomo) {

    // Loop over all the pixels in the cone
    Utility::parallelize(ntrajectory, [=, &parameters, &map_components,
                                       &index_components, &ntrajectory,
                                       &firsttrajectory, &octree, &vobs, &map,
                                       &nmaps, &pixel, &cosmology, &observer,
                                       &length, &interpRefvec, &rhomo, &thomo,
                                       &lambdahomo, &redshifthomo,
                                       &ahomo](const unsigned int itrajectory) {
        const uint itrajectorys = itrajectory + firsttrajectory;
        magrathea::Evolution<Photon<double, 3>> trajectorycenter,
            trajectorycenter_born;
        Photon<double, 3> photoncenter;
        std::vector<std::array<std::array<double, 2>, 2>> jacobian(
            parameters.nb_z_maps),
            jacobian_born(parameters.nb_z_maps);
        std::vector<uint> firstid(parameters.nb_z_maps);
        std::vector<Real> interpRefvecBundle(parameters.nb_z_maps),
            f(parameters.nb_z_maps), distance(parameters.nb_z_maps);
        double phi(0), theta(0);
        unsigned int marked(0);
        double vec1[3];
        const bool only_born =
            (index_components.size() == 1 && index_components[0] == 1);
        const bool compute_born =
            std::find(map_components.begin(), map_components.end(),
                      "lensing_born") != map_components.end();
        const bool compute_lensing =
            std::find(map_components.begin(), map_components.end(), "lensing") !=
            map_components.end();
        std::vector<point> kiTargets(parameters.nb_z_maps),
            posTargets(parameters.nb_z_maps);
        // Convert pixel to 3D vector
        pix2vec_ring(parameters.nside, pixel[itrajectorys], vec1);
        // Convert to angles
        phi = std::atan2(vec1[1], vec1[0]);
        theta = std::acos(vec1[2]);

        // If we only use 'lensing_born', then no need for real ray-tracing
        if (!only_born) {
            // Integrate central ray
            photoncenter = Integrator::launch(observer[0], observer[1], observer[2],
                                              vec1[0], vec1[1], vec1[2]);
            trajectorycenter.append(photoncenter);
            // Integrate photon trajectory
            Integrator::integrate(trajectorycenter, parameters.stop_ray,
                                  interpRefvec.back(), cosmology, octree, vobs,
                                  length, parameters.nsteps);
            // Integrator::integrate(trajectorycenter, cosmology, octree, vobs,
            // length, parameters.nsteps); // Uncomment to compute upper surface (with
            // stop_rays == redshift)

            // Get index and factors for constant-parameter surfaces
            // Loop over all the redshifts of interest
            for (uint i = 0; i < parameters.nb_z_maps; i++) {
                // Stop central ray at some parameter
                if (parameters.stop_ray == "r") {
                    photoncenter.chi() = interpRefvec[i];
                    marked = std::distance(
                        std::begin(trajectorycenter),
                        std::upper_bound(std::begin(trajectorycenter),
                                         std::end(trajectorycenter), photoncenter,
                                         [](const Photon<double, 3> &first,
                                            const Photon<double, 3> &second) {
                                             return first.chi() < second.chi();
                                         }));
                    firstid[i] = marked - (marked > 0);
                    f[i] = (trajectorycenter[firstid[i] + 1].chi() - interpRefvec[i]) /
                           (trajectorycenter[firstid[i] + 1].chi() -
                            trajectorycenter[firstid[i]].chi());
                } else if (parameters.stop_ray == "lambda") {
                    photoncenter.lambda() = interpRefvec[i];
                    marked = std::distance(
                        std::begin(trajectorycenter),
                        std::upper_bound(std::begin(trajectorycenter),
                                         std::end(trajectorycenter), photoncenter,
                                         [](const Photon<double, 3> &first,
                                            const Photon<double, 3> &second) {
                                             return first.lambda() < second.lambda();
                                         }));
                    firstid[i] = marked - (marked > 0);
                    f[i] = (trajectorycenter[firstid[i] + 1].lambda() - interpRefvec[i]) /
                           (trajectorycenter[firstid[i] + 1].lambda() -
                            trajectorycenter[firstid[i]].lambda());
                } else if (parameters.stop_ray == "t") {
                    photoncenter.t() = interpRefvec[i];
                    marked = std::distance(
                        std::begin(trajectorycenter),
                        std::upper_bound(std::begin(trajectorycenter),
                                         std::end(trajectorycenter), photoncenter,
                                         [](const Photon<double, 3> &first,
                                            const Photon<double, 3> &second) {
                                             return first.t() < second.t();
                                         }));
                    firstid[i] = marked - (marked > 0);
                    f[i] = (trajectorycenter[firstid[i] + 1].t() - interpRefvec[i]) /
                           (trajectorycenter[firstid[i] + 1].t() -
                            trajectorycenter[firstid[i]].t());
                } else if (parameters.stop_ray == "redshift") {
                    photoncenter.redshift() = interpRefvec[i];
                    marked = std::distance(
                        std::begin(trajectorycenter),
                        std::find_if(std::begin(trajectorycenter),
                                     std::end(trajectorycenter),
                                     [=, &photoncenter](const Photon<double, 3> &first) {
                                         return first.redshift() > photoncenter.redshift();
                                     })); // lower surface
                    // marked = std::distance(std::begin(trajectorycenter),
                    // std::find_if(std::rbegin(trajectorycenter),
                    // std::rend(trajectorycenter), [=, &photoncenter](const
                    // Photon<double, 3>& first){return first.redshift() <
                    // photoncenter.redshift();}).base()); // Uncomment to compute upper
                    // surface (with stop_rays == redshift)
                    firstid[i] = marked - (marked > 0);
                    f[i] =
                        (trajectorycenter[firstid[i] + 1].redshift() - interpRefvec[i]) /
                        (trajectorycenter[firstid[i] + 1].redshift() -
                         trajectorycenter[firstid[i]].redshift());
                } else if (parameters.stop_ray == "a") {
                    photoncenter.a() = interpRefvec[i];
                    marked = std::distance(
                        std::begin(trajectorycenter),
                        std::upper_bound(std::begin(trajectorycenter),
                                         std::end(trajectorycenter), photoncenter,
                                         [](const Photon<double, 3> &first,
                                            const Photon<double, 3> &second) {
                                             return first.a() > second.a();
                                         }));
                    firstid[i] = marked - (marked > 0);
                    f[i] = (trajectorycenter[firstid[i] + 1].a() - interpRefvec[i]) /
                           (trajectorycenter[firstid[i] + 1].a() -
                            trajectorycenter[firstid[i]].a());
                } else if (parameters.stop_ray == "s") {
                    photoncenter.s() = interpRefvec[i];
                    marked = std::distance(
                        std::begin(trajectorycenter),
                        std::upper_bound(std::begin(trajectorycenter),
                                         std::end(trajectorycenter), photoncenter,
                                         [](const Photon<double, 3> &first,
                                            const Photon<double, 3> &second) {
                                             return first.s() < second.s();
                                         }));
                    firstid[i] = marked - (marked > 0);
                    f[i] = (trajectorycenter[firstid[i] + 1].s() - interpRefvec[i]) /
                           (trajectorycenter[firstid[i] + 1].s() -
                            trajectorycenter[firstid[i]].s());
                }
            }

            // Check if we need to compute the lensing jacobian matrix
            if (compute_lensing) {
                for (uint i = 0; i < parameters.nb_z_maps; i++) { // different maps
                    // Get position of photon at the surface and the plane (screen for
                    // bundle method)
                    posTargets[i][0] = trajectorycenter[firstid[i]].x() * f[i] +
                                       trajectorycenter[firstid[i] + 1].x() * (1 - f[i]);
                    posTargets[i][1] = trajectorycenter[firstid[i]].y() * f[i] +
                                       trajectorycenter[firstid[i] + 1].y() * (1 - f[i]);
                    posTargets[i][2] = trajectorycenter[firstid[i]].z() * f[i] +
                                       trajectorycenter[firstid[i] + 1].z() * (1 - f[i]);
                    distance[i] = std::sqrt(posTargets[i][0] * posTargets[i][0] +
                                            posTargets[i][1] * posTargets[i][1] +
                                            posTargets[i][2] * posTargets[i][2]);
                    if (parameters.beam == "bundle") {
                        if (parameters.plane == "sachs") {
                            kiTargets[i][0] =
                                trajectorycenter[firstid[i]].dxdl() * f[i] +
                                trajectorycenter[firstid[i] + 1].dxdl() * (1 - f[i]);
                            kiTargets[i][1] =
                                trajectorycenter[firstid[i]].dydl() * f[i] +
                                trajectorycenter[firstid[i] + 1].dydl() * (1 - f[i]);
                            kiTargets[i][2] =
                                trajectorycenter[firstid[i]].dzdl() * f[i] +
                                trajectorycenter[firstid[i] + 1].dzdl() * (1 - f[i]);
                        } else if (parameters.plane == "normal") {
                            kiTargets = posTargets;
                        } else {
                            std::cout << "# WARNING: with beam = 'bundle', please choose "
                                         "plane = 'normal' or 'sachs'"
                                      << std::endl;
                            std::cout << "# Error at file " << __FILE__
                                      << ", line : " << __LINE__ << std::endl;
                            std::terminate();
                        }
                    }
                    // Stop criterion for bundle
                    if (parameters.beam == "bundle") {
                        if (parameters.stop_bundle == "lambda") {
                            interpRefvecBundle[i] =
                                trajectorycenter[firstid[i]].lambda() * f[i] +
                                trajectorycenter[firstid[i] + 1].lambda() * (1 - f[i]);
                        } else if (parameters.stop_bundle == "r") {
                            interpRefvecBundle[i] =
                                trajectorycenter[firstid[i]].chi() * f[i] +
                                trajectorycenter[firstid[i] + 1].chi() * (1 - f[i]);
                            ;
                        } else if (parameters.stop_bundle == "redshift") {
                            interpRefvecBundle[i] =
                                trajectorycenter[firstid[i]].redshift() * f[i] +
                                trajectorycenter[firstid[i] + 1].redshift() * (1 - f[i]);
                            ;
                        } else if (parameters.stop_bundle == "a") {
                            interpRefvecBundle[i] =
                                trajectorycenter[firstid[i]].a() * f[i] +
                                trajectorycenter[firstid[i] + 1].a() * (1 - f[i]);
                        } else if (parameters.stop_bundle == "t") {
                            interpRefvecBundle[i] =
                                trajectorycenter[firstid[i]].t() * f[i] +
                                trajectorycenter[firstid[i] + 1].t() * (1 - f[i]);
                        } else if (parameters.stop_bundle == "plane") {
                            interpRefvecBundle[i] = kiTargets[i][0] * posTargets[i][0] +
                                                    kiTargets[i][1] * posTargets[i][1] +
                                                    kiTargets[i][2] * posTargets[i][2];
                        } else {
                            std::cout << "# WARNING: Please choose an existing stop_bundle "
                                         "parameter"
                                      << std::endl;
                            std::cout << "# Error at file " << __FILE__
                                      << ", line : " << __LINE__ << std::endl;
                            std::terminate();
                        }
                    }
                }
                // Compute Lensing Jacobian matrix
                if (parameters.beam == "bundle") {
                    jacobian = Lensing::dbetadtheta(
                        parameters, kiTargets, interpRefvecBundle, observer, phi, theta,
                        distance, cosmology, octree, vobs, length);
                } else if (parameters.beam == "infinitesimal") {
                    jacobian = Lensing::dbetadtheta_infinitesimal(
                        distance, trajectorycenter, octree, length);
                } else {
                    std::cout << "# WARNING: beam must be 'bundle' or 'infinitesimal'"
                              << std::endl;
                    std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                              << std::endl;
                    std::terminate();
                }
            }
        }
        // Check if we need to compute the lensing jacobian matrix with Born
        // approximation
        if (compute_born) {
            photoncenter = Integrator::launch(observer[0], observer[1], observer[2],
                                              vec1[0], vec1[1], vec1[2]);
            trajectorycenter_born.append(photoncenter);
            // Integrate photon on a FLRW trajectory
            Integrator::integrate<-1>(trajectorycenter_born, parameters.stop_ray,
                                      interpRefvec.back(), cosmology, octree, vobs,
                                      length, parameters.nsteps);
            // Compute distortions around the FLRW trajectory
            jacobian_born = Lensing::dbetadtheta_infinitesimal(
                rhomo, trajectorycenter_born, octree, length);
        }

        // Fill maps
        for (uint iz = 0; iz < parameters.nb_z_maps; iz++) { // different maps
            // Check in trajectory is fine
            if ((compute_born && trajectorycenter_born.size() == 1) ||
                (!only_born &&
                 (firstid[iz] == trajectorycenter.size() - 1 || firstid[iz] == 0))) {
                for (uint i = 0; i < nmaps; i++) {
                    map[nmaps * iz + i][pixel[itrajectorys]] = 0;
                }
            } else {
                uint icomp(0);
                // Loops over map types
                for (uint j = 0; j < map_components.size(); j++) {
                    if (index_components[j] == 0) {
                        // For 'bundle' need to compute the image rotation, while for
                        // 'infinitesimal' we assume that there isn't
                        if (parameters.beam == "bundle") {
                            const double a11(jacobian[iz][0][0]), a12(jacobian[iz][0][1]),
                                a21(jacobian[iz][1][0]), a22(jacobian[iz][1][1]);
                            // Check if jacobian is good
                            if (a11 == 42) {
                                for (uint i = 0; i < nmaps; i++) {
                                    map[nmaps * iz + i][pixel[itrajectorys]] = 0;
                                }
                                break;
                            } else {
                                const double w = -std::atan2(a12 - a21, a11 + a22);
                                const double kappa = 1 - (a11 + a22) * 0.5 / std::cos(w);
                                const double gamma1 = -0.5 * (std::cos(w) * (a11 - a22) +
                                                              std::sin(w) * (a12 + a21));
                                const double gamma2 = -0.5 * (std::cos(w) * (a12 + a21) +
                                                              std::sin(w) * (a22 - a11));
                                const double invmagnification = a11 * a22 - a12 * a21;
                                map[nmaps * iz + icomp][pixel[itrajectorys]] = kappa;
                                map[nmaps * iz + icomp + 1][pixel[itrajectorys]] = gamma1;
                                map[nmaps * iz + icomp + 2][pixel[itrajectorys]] = gamma2;
                                map[nmaps * iz + icomp + 3][pixel[itrajectorys]] = w;
                                map[nmaps * iz + icomp + 4][pixel[itrajectorys]] =
                                    invmagnification;
                                icomp += 4;
                            }
                        } else if (parameters.beam == "infinitesimal") {
                            const double a11(jacobian[iz][0][0]), a12(jacobian[iz][0][1]),
                                a22(jacobian[iz][1][1]);
                            // Check if jacobian is good
                            if (a11 == 42) {
                                for (uint i = 0; i < nmaps; i++) {
                                    map[nmaps * iz + i][pixel[itrajectorys]] = 0;
                                }
                                break;
                            } else {
                                const double kappa = 1 - (a11 + a22) * 0.5;
                                const double gamma1 = -0.5 * (a11 - a22);
                                const double gamma2 = -a12;
                                const double invmagnification = a11 * a22 - a12 * a12;
                                map[nmaps * iz + icomp][pixel[itrajectorys]] = kappa;
                                map[nmaps * iz + icomp + 1][pixel[itrajectorys]] = gamma1;
                                map[nmaps * iz + icomp + 2][pixel[itrajectorys]] = gamma2;
                                map[nmaps * iz + icomp + 3][pixel[itrajectorys]] =
                                    invmagnification;
                                icomp += 3;
                            }
                        }
                    } else if (index_components[j] == 1) {
                        const double a11(jacobian_born[iz][0][0]),
                            a12(jacobian_born[iz][0][1]), a22(jacobian_born[iz][1][1]);
                        // Check if jacobian is good
                        if (a11 == 42) {
                            for (uint i = 0; i < nmaps; i++) {
                                map[nmaps * iz + i][pixel[itrajectorys]] = 0;
                            }
                            break;
                        } else {
                            const double kappa = 1 - (a11 + a22) * 0.5;
                            const double gamma1 = -0.5 * (a11 - a22);
                            const double gamma2 = -a12;
                            const double invmagnification = a11 * a22 - a12 * a12;
                            map[nmaps * iz + icomp][pixel[itrajectorys]] = kappa;
                            map[nmaps * iz + icomp + 1][pixel[itrajectorys]] = gamma1;
                            map[nmaps * iz + icomp + 2][pixel[itrajectorys]] = gamma2;
                            map[nmaps * iz + icomp + 3][pixel[itrajectorys]] =
                                invmagnification;
                            icomp += 3;
                        }
                    } else if (index_components[j] == 2) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].chi() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].chi() * (1 - f[iz])) /
                                rhomo[iz] -
                            1;
                    } else if (index_components[j] == 3) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].lambda() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].lambda() * (1 - f[iz])) /
                                lambdahomo[iz] -
                            1;
                    } else if (index_components[j] == 4) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].t() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].t() * (1 - f[iz])) /
                                thomo[iz] -
                            1;
                    } else if (index_components[j] == 5) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].a() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].a() * (1 - f[iz])) /
                                ahomo[iz] -
                            1;
                    } else if (index_components[j] == 6) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].redshift() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].redshift() * (1 - f[iz])) /
                                redshifthomo[iz] -
                            1;
                    } else if (index_components[j] == 7) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].s() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].s() * (1 - f[iz])) /
                                rhomo[iz] -
                            1;
                    } else if (index_components[j] == 8) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].isw() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].isw() * (1 - f[iz]));
                    } else if (index_components[j] == 9) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].rho() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].rho() * (1 - f[iz]));
                    } else if (index_components[j] == 10) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (*std::max_element(std::begin(trajectorycenter),
                                               std::begin(trajectorycenter) + firstid[iz],
                                               [](const Photon<double, 3> &first,
                                                  const Photon<double, 3> &second) {
                                                   return first.rho() < second.rho();
                                               }))
                                .rho();
                    } else if (index_components[j] == 11) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] = firstid[iz];
                    } else if (index_components[j] == 12) {
                        map[nmaps * iz + icomp][pixel[itrajectorys]] =
                            (trajectorycenter[firstid[iz]].phi() * f[iz] +
                             trajectorycenter[firstid[iz] + 1].phi() * (1 - f[iz]));
                    } else if (index_components[j] == 13) {
                        double beta1(0), beta2(0);
                        if (compute_lensing) {
                            beta1 = std::atan2(posTargets[iz][1], posTargets[iz][0]);
                            beta2 = std::acos(posTargets[iz][2] / distance[iz]);
                        } else {
                            double x = trajectorycenter[firstid[iz]].x() * f[iz] +
                                       trajectorycenter[firstid[iz] + 1].x() * (1 - f[iz]);
                            double y = trajectorycenter[firstid[iz]].y() * f[iz] +
                                       trajectorycenter[firstid[iz] + 1].y() * (1 - f[iz]);
                            double z = trajectorycenter[firstid[iz]].z() * f[iz] +
                                       trajectorycenter[firstid[iz] + 1].z() * (1 - f[iz]);
                            double r = std::sqrt(x * x + y * y + z * z);
                            beta1 = std::atan2(y, x);
                            beta2 = std::acos(z / r);
                        }
                        map[nmaps * iz + icomp][pixel[itrajectorys]] = phi - beta1;
                        map[nmaps * iz + icomp + 1][pixel[itrajectorys]] = theta - beta2;
                        icomp++;
                    }
                    icomp++;
                }
            }
        }
    });
}

// Fill Healpix maps with spherical bundles [UNUSED]
/// \brief          Produce Healpix maps using bundles
/// \details        Produce Healpix maps using bundles at constant FLRW redshift
/// \tparam         Parameter Parameter type
/// \tparam         Octree octree type
/// \tparam         Map map type
/// \tparam         Pixel pixel type
/// \tparam         Cosmology cosmology type
/// \tparam         Integer integer type
/// \tparam         Real real/float type
/// \param[in]      parameters Parameter structure
/// \param[in]      ntrajectory Number of pixels to be filled
/// \param[in]      firsttrajectory first index for pixel for a given cone
/// \param[in]      octree Octree
/// \param[in]      vobs Observer velocity
/// \param[in,out]  map Map
/// \param[in]      nmaps Number of Healpix maps
/// \param[in]      pixel Pixel vector
/// \param[in]      cosmology Cosmological table
/// \param[in]      observer Observer position
/// \param[in]      length R.U to S.I units for length
/// \param[in]      z_stop_vec Vector containing redshifts at which we compute
/// scalar quantities
template <class Parameter, class Octree, class Map, class Pixel,
          class Cosmology, class Point, typename Integer, typename Real>
void Hmaps::FillMapPropagate(const Parameter &parameters,
                             const Integer ntrajectory,
                             const Integer firsttrajectory,
                             const Octree &octree, const Point &vobs, Map &map,
                             const Integer nmaps, const Pixel &pixel,
                             const Cosmology &cosmology, const point &observer,
                             const Real length,
                             const std::vector<Real> &z_stop_vec) {

    const double amin = one / (one + z_stop_vec.back());
    magrathea::Evolution<Photon<double, 3>> reference;
    Photon<double, 3> photonref;
    // Launch a ray toward the x-axis
    photonref = Integrator::launch(0., 0., 0., 1., 0., 0.);
    const point vobs0 = {0, 0,
                         0}; // No peculiar velocity for homogeneous quantities
    std::vector<unsigned long int> firstid_ref(parameters.nb_z_maps);
    std::vector<double> f_ref(parameters.nb_z_maps),
        distance_ref(parameters.nb_z_maps);
    Octree homotree;
    // Create homogeneous octree
    homotree.assign(parameters.ncoarse / 2., zero);
    // Compute FLRW ray in the homogeneous octree
    reference = Integrator::propagate<-1>(
        photonref, parameters.nbundlemin, parameters.openingmin, real(),
        parameters.stop_ray, cosmology, homotree, vobs0, length,
        parameters.nsteps * (1 << (parameters.ncoarse - parameters.ncoarse / 2)) *
            2,
        double());

    // Loop over all the redshifts to get angular distance and interpolation
    // factors
    for (uint iz = 0; iz < parameters.nb_z_maps; iz++) {
        photonref.a() = one / (one + z_stop_vec[iz]);
        const unsigned long int marked = std::distance(
            std::begin(reference),
            std::upper_bound(std::begin(reference), std::end(reference), photonref,
                             [](const Photon<double, 3> &first,
                                const Photon<double, 3> &second) {
                                 return first.a() > second.a();
                             }));
        firstid_ref[iz] = marked - (marked > 0);
        const double previous = reference[firstid_ref[iz]].a();
        const double next = reference[firstid_ref[iz] + 1].a();
        f_ref[iz] = (next - photonref.a()) / (next - previous);
        distance_ref[iz] =
            reference[firstid_ref[iz]].distance() * f_ref[iz] +
            reference[firstid_ref[iz] + 1].distance() * (1 - f_ref[iz]);
    }
    // Loop over the pixels
    Utility::parallelize(ntrajectory, [=, &map, &octree, &distance_ref,
                                       &parameters, &cosmology,
                                       &pixel](const uint itrajectory) {
        const uint itrajectorys = itrajectory + firsttrajectory;
        magrathea::Evolution<Photon<double, 3>> trajectorycenter;
        Photon<double, 3> photoncenter;
        std::vector<std::array<std::array<double, 2>, 2>> jacobian(
            parameters.nb_z_maps);
        std::vector<uint> firstid(parameters.nb_z_maps);
        double vec1[3];
        std::vector<double> f(parameters.nb_z_maps), distance(parameters.nb_z_maps);
        std::vector<point> kiTargets(parameters.nb_z_maps),
            posTargets(parameters.nb_z_maps);
        // Convert pixel to 3D vector
        pix2vec_ring(parameters.nside, pixel[itrajectorys], vec1);
        // Launch photon toward pixel
        photoncenter = Integrator::launch(observer[0], observer[1], observer[2],
                                          vec1[0], vec1[1], vec1[2]);
        // Propagate spherical bundle in the direction of pixel
        trajectorycenter = Integrator::propagate(
            photoncenter, parameters.nbundlemin, parameters.openingmin, double(),
            parameters.stop_ray, cosmology, octree, vobs0, length,
            parameters.nsteps, amin);

        // Check if trajectory is good
        if (trajectorycenter.size()) {
            // Loop over all the redshifts of reference to get the angular diameter
            // distance
            for (uint iz = 0; iz < parameters.nb_z_maps; iz++) {
                photoncenter.a() = one / (one + z_stop_vec[iz]);
                const unsigned long int marked = std::distance(
                    std::begin(trajectorycenter),
                    std::upper_bound(std::begin(trajectorycenter),
                                     std::end(trajectorycenter), photoncenter,
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.a() > second.a();
                                     }));
                firstid[iz] = marked - (marked > 0);
                const double previous = trajectorycenter[firstid[iz]].a();
                const double next = trajectorycenter[firstid[iz] + 1].a();
                f[iz] = (next - photoncenter.a()) / (next - previous);
                distance[iz] =
                    trajectorycenter[firstid[iz]].distance() * f[iz] +
                    trajectorycenter[firstid[iz] + 1].distance() * (1 - f[iz]);
            }
            // Loop over all the redshifts of reference
            for (uint iz = 0; iz < parameters.nb_z_maps; iz++) {
                // Check if everything is fine
                if (firstid[iz] == trajectorycenter.size() - 1 || firstid[iz] == 0) {
                    for (uint i = 0; i < nmaps; i++) {
                        map[nmaps * iz + i][pixel[itrajectorys]] = 0;
                    }
                } else {
                    const double dd0 = distance[iz] / distance_ref[iz];
                    const double sf = trajectorycenter[firstid[iz]].a() * f[iz] +
                                      trajectorycenter[firstid[iz] + 1].a() * (1 - f[iz]);
                    const double chi =
                        trajectorycenter[firstid[iz]].chi() * f[iz] +
                        trajectorycenter[firstid[iz] + 1].chi() * (1 - f[iz]);
                    map[nmaps * iz + 0][pixel[itrajectorys]] = dd0;
                    map[nmaps * iz + 1][pixel[itrajectorys]] = 1. / sf - 1;
                    map[nmaps * iz + 2][pixel[itrajectorys]] =
                        chi * sf / (distance_ref[iz] / length);
                    map[nmaps * iz + 3][pixel[itrajectorys]] = chi;
                }
            }
        } else {
            for (uint iz = 0; iz < parameters.nb_z_maps; iz++) {
                for (uint i = 0; i < nmaps; i++) {
                    map[nmaps * iz + i][pixel[itrajectorys]] = 0;
                }
            }
        }
    });
}

#endif // HMAPS_H_INCLUDED
