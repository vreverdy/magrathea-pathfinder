/* ********************************** GENERATE_CONES
 * ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Generate_cones
// DESCRIPTION :    Some cone generation function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           generate_cones.h
/// \brief          Some cone generation function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
/// (michel-andres.breton@obspm.fr) \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef GENERATE_CONES_H_INCLUDED
#define GENERATE_CONES_H_INCLUDED

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

using namespace magrathea;

struct parameters_t {
    // Common
    std::string celldir;
    std::string conedir;
    uint typefile;
    uint isfullsky;
    uint ncones;
    real buffer;

} parameters;

class Generate_cones {
    // Methodes
public:
    // Read parameter file
    template <class Parameters, class Map>
    static void ReadParamFile(Parameters &parameters, Map &parameter);

    // Cones generation
    template <class Cone, template <unsigned int, class, typename> class Sphere,
              unsigned int Dimension, class Vector, typename Scalar,
              typename Integer>
    static void GenerateFullskyCones(const Integer ncones,
                                     std::vector<Cone> &cone,
                                     std::vector<Cone> &coneIfRot,
                                     Sphere<Dimension, Vector, Scalar> &sphere);
    template <class Parameter, class Cone,
              template <unsigned int, class, typename> class Sphere,
              unsigned int Dimension, class Vector, typename Scalar>
    static void GenerateNarrowCones(const Parameter &parameters,
                                    std::vector<Cone> &cone,
                                    std::vector<Cone> &coneIfRot,
                                    Sphere<Dimension, Vector, Scalar> &sphere,
                                    std::array<std::array<double, 3>, 3> &rotm1,
                                    Scalar &thetay, Scalar &thetaz);
};

// Read parameter file
/// \brief          Read parameter file.
/// \details        Read and put in a structure the parameters.
/// \tparam         Parameters structure type
/// \tparam         Map map type
/// \param[in,out]  parameters Structure containing the parameters.
/// \param[in]      parameter Contains parameters to be rewritten
template <class Parameters, class Map>
void Generate_cones::ReadParamFile(Parameters &parameters, Map &parameter) {
    parameters.celldir = parameter["celldir"];
    parameters.conedir = parameter["conedir"];
    parameters.typefile = std::stoul(parameter["typefile"]);
    parameters.isfullsky = std::stoul(parameter["isfullsky"]);
    parameters.ncones = std::stoul(parameter["ncones"]);
    parameters.buffer = std::stod(parameter["buffer"]);
}

// Fullsky cones generation
/// \brief          Generate fullsky cones.
/// \details        Generate fullsky cones.
/// \tparam         Cone cone type
/// \tparam         Sphere sphere type
/// \tparam         Dimension dimension
/// \tparam         Vector vector type
/// \tparam         Scalar scalar type
/// \tparam         Integer integer type
/// \param[in]      ncones Number of cones .
/// \param[in,out]  cone Cones generated
/// \param[in]      coneIfRot Cones generated
/// \param[in]      sphere Central sphere
/// \return         Generate the shape of fullsky cones.
template <class Cone, template <unsigned int, class, typename> class Sphere,
          unsigned int Dimension, class Vector, typename Scalar,
          typename Integer>
void Generate_cones::GenerateFullskyCones(
    const Integer ncones, std::vector<Cone> &cone, std::vector<Cone> &coneIfRot,
    Sphere<Dimension, Vector, Scalar> &sphere) {
    std::vector<std::array<double, 3>> tiling(ncones);
    // Cone angle from the maximum distance between points generated on the
    // sphere. We multiply by an arbitrary factor which seems ideal to produce
    // wide enough cones
    double alpha = 1.8 * std::asin(sphere
                                       .template uniform<Dimension - 1>(
                                           std::begin(tiling), std::end(tiling))
                                       .first /
                                   sphere.diameter());
    // Assign properties to each cone
    Utility::parallelize(ncones,
                         [=, &tiling, &cone, &sphere, &alpha](const uint i) {
                             cone[i].assign(sphere.position(), tiling[i], alpha);
                         });
    // No rotation for fullsky cones
    for (uint i = 0; i < ncones; i++) {
        coneIfRot[i] = cone[i];
    }
}

// Narrow cones generation
/// \brief          Generate narrow cones.
/// \details        Generate narrow cones.
/// \tparam         Parameter Parameter type.
/// \tparam         Cone cone type
/// \tparam         Sphere sphere type
/// \tparam         Dimension dimension
/// \tparam         Vector vector type
/// \tparam         Scalar scalar type
/// \param[in]      parameters Parameter structure
/// \param[in,out]  cone Cones generated
/// \param[in,out]      coneIfRot Cones generated
/// \param[in]      sphere Central sphere
/// \param[in,out]  rotm1 Rotation matrix for narrow cone cells
/// \param[in,out]  thetay Semi-angle for solid angle in direction y
/// \param[in,out]  thetaz Semi-angle for solid angle in direction z
template <class Parameter, class Cone,
          template <unsigned int, class, typename> class Sphere,
          unsigned int Dimension, class Vector, typename Scalar>
void Generate_cones::GenerateNarrowCones(
    const Parameter &parameters, std::vector<Cone> &cone,
    std::vector<Cone> &coneIfRot, Sphere<Dimension, Vector, Scalar> &sphere,
    std::array<std::array<double, 3>, 3> &rotm1, Scalar &thetay,
    Scalar &thetaz) {
    static constexpr double pi = Constants<double>::pi();
    std::size_t found;
    std::vector<std::array<double, 3>> tiling(parameters.ncones),
        tilingbis(parameters.ncones);
    double theta_rot(0), phi_rot(0);
    std::array<std::array<double, 3>, 3> rotation = {{0}};
    std::vector<std::string> filelistingprior;
    std::string filelisting;

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.celldir, filelistingprior);
    for (uint ifiling = 0; ifiling < filelistingprior.size(); ++ifiling) {
        if (parameters.typefile == 1) {
            // Find any hdf5 file
            found = filelistingprior[ifiling].find(".h5");
            if (found != std::string::npos) {
                filelisting = parameters.celldir + filelistingprior[ifiling];
                break;
            }
        } else if (parameters.typefile == 2) {
            // For ASCII, search for the 'info_narrow_cone.txt' file
            if (filelistingprior[ifiling] == "info_narrow_cone.txt") {
                filelisting = parameters.celldir + filelistingprior[ifiling];
                break;
            }
        } else {
            std::cout << "# WARNING : Narrow cones can only be computed using HDF5 "
                         "or ASCII input files"
                      << std::endl;
            std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                      << std::endl;
            std::terminate();
        }
    }
    // Get informations from HDF5 files
    if (parameters.typefile == 1) {
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "phi", phi_rot);
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "theta",
                                theta_rot);
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "thetay",
                                thetay);
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "thetaz",
                                thetaz);
        // Get informations from ASCII files
    } else if (parameters.typefile == 2) {
        std::map<std::string, std::string> parameterASCII;
        parameterASCII = Input::parse(filelisting);
        phi_rot = std::stoul(parameterASCII["phi"]);
        theta_rot = std::stod(parameterASCII["theta"]);
        thetay = std::stod(parameterASCII["thetay"]);
        thetaz = std::stod(parameterASCII["thetaz"]);
    } else {
        std::cout << "# WARNING : Narrow cones can only be computed using HDF5 or "
                     "ASCII input files"
                  << std::endl;
        std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                  << std::endl;
        std::terminate();
    }
    thetaz *= magrathea::Constants<double>::deg();
    thetay *= magrathea::Constants<double>::deg();
    phi_rot *= magrathea::Constants<double>::deg();
    theta_rot *= magrathea::Constants<double>::deg();

    // Compute the rotation matrix to rotate the cone
    rotation[0][0] = std::cos(theta_rot) * std::cos(phi_rot);
    rotation[0][1] = std::cos(theta_rot) * std::sin(phi_rot);
    rotation[0][2] = -std::sin(theta_rot);
    rotation[1][0] = -std::sin(phi_rot);
    rotation[1][1] = std::cos(phi_rot);
    rotation[1][2] = 0;
    rotation[2][0] = std::cos(phi_rot) * std::sin(theta_rot);
    rotation[2][1] = std::sin(phi_rot) * std::sin(theta_rot);
    rotation[2][2] = std::cos(theta_rot);
    // Inverse matrix
    rotm1 = Utility::invMatrix3d(rotation);
    uint iloop(0);
    std::array<double, 3> ciblage;
    ciblage[0] = 1;
    ciblage[1] = 0;
    ciblage[2] = 0;
    double fullsky = 4 * pi;
    double portion = 2 * thetay * (std::sin(thetaz) - std::cos(pi / 2. + thetaz));
    // Inverse fraction of the sky
    uint fsp = fullsky / portion;
    tiling.resize(parameters.ncones * fsp);

    // Generate random points on the full sky. Then check is we have a number of
    // points equal to 'ncones' in the area of interest
    sphere.template uniform<Dimension - 1>(std::begin(tiling), std::end(tiling));
    for (uint i = 0; i < tiling.size(); ++i) {
        double phi = std::atan2(tiling[i][1], tiling[i][0]);
        double theta = std::acos(tiling[i][2] / sphere.radius());
        // If point is inside the region of interest
        if (std::abs(phi) < thetay - parameters.buffer &&
            theta > pi / 2 - thetaz + parameters.buffer &&
            theta < pi / 2 + thetaz - parameters.buffer) {
            iloop++;
        }
    }
    float irecuploop = 0;
    // Generally, we will not exactly have the good number of points, then need to
    // iterate on the initial number of points on the full sky
    while (iloop != parameters.ncones) {
        if (irecuploop < 10)
            tiling.resize(tiling.size() +
                          static_cast<int>((static_cast<float>(parameters.ncones) -
                                            static_cast<float>(iloop)) *
                                           static_cast<float>(fsp) *
                                           (1 / (1 + irecuploop))));
        else
            tiling.resize(tiling.size() + parameters.ncones - iloop);
        sphere.template uniform<Dimension - 1>(std::begin(tiling),
                                               std::end(tiling));
        iloop = 0;
        for (uint i = 0; i < tiling.size(); ++i) {
            double phi = std::atan2(tiling[i][1], tiling[i][0]);
            double theta = std::acos(tiling[i][2] / sphere.radius());
            if (std::abs(phi) < thetay - parameters.buffer &&
                theta > pi / 2 - thetaz + parameters.buffer &&
                theta < pi / 2 + thetaz - parameters.buffer) {
                iloop++;
            }
        }
        irecuploop++;
    }
    iloop = 0;
    // Put final result in tilingbis
    for (uint i = 0; i < tiling.size(); ++i) {
        double phi = std::atan2(tiling[i][1], tiling[i][0]);
        double theta = std::acos(tiling[i][2] / sphere.radius());
        if (std::abs(phi) < thetay - parameters.buffer &&
            theta > pi / 2 - thetaz + parameters.buffer &&
            theta < pi / 2 + thetaz - parameters.buffer) {
            tilingbis[iloop][0] = tiling[i][0];
            tilingbis[iloop][1] = tiling[i][1];
            tilingbis[iloop][2] = tiling[i][2];
            iloop++;
        }
    }

    double resulting(0);
    resulting = sphere.diameter();
    // Check the maximum distance between two cones
    for (unsigned int i = 0; i < parameters.ncones; ++i) {
        for (unsigned int j = 0; j < parameters.ncones; ++j) {
            if (i != j) {
                double tmp = 0;
                for (unsigned int idim = 0; idim < 3; ++idim) {
                    tmp += std::pow(tilingbis[i][idim] - tilingbis[j][idim], 2);
                }
                tmp = std::sqrt(tmp);
                resulting = (resulting < tmp) ? (resulting) : (tmp);
            }
        }
    }

    // Assign an angle depending on the maximum distance between two cones. We
    // multiply by an arbitrary factor which seems ideal to produce wide enough
    // cones
    double alpha = 1.8 * std::asin(resulting / sphere.diameter());
    // Need a minimum angle to avoid thin cones. 0.1 = 6 degrees
    const double anglemin = 0.001;
    alpha = (anglemin > alpha) ? anglemin : alpha;
    std::cout << "# Angle proposed : "
              << 1.8 * std::asin(resulting / sphere.diameter())
              << " angle chosen : " << alpha << std::endl;
    // Assign properties to cones
    Utility::parallelize(parameters.ncones, [=, &tilingbis, &cone](const uint i) {
        cone[i].assign(sphere.position(), tilingbis[i], alpha);
    });
    // For narrow cones, need rotation
    for (uint i = 0; i < parameters.ncones; i++) {
        coneIfRot[i] = cone[i];
        coneIfRot[i].base(0) = cone[i].base(0) * rotm1[0][0] +
                               cone[i].base(1) * rotm1[0][1] +
                               cone[i].base(2) * rotm1[0][2];
        coneIfRot[i].base(1) = cone[i].base(0) * rotm1[1][0] +
                               cone[i].base(1) * rotm1[1][1] +
                               cone[i].base(2) * rotm1[1][2];
        coneIfRot[i].base(2) = cone[i].base(0) * rotm1[2][0] +
                               cone[i].base(1) * rotm1[2][1] +
                               cone[i].base(2) * rotm1[2][2];
    }
}

#endif // GENERATE_CONES_H_INCLUDED
