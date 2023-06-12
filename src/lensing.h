/* ******************************* LENSING ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Lensing
// DESCRIPTION :    Integration utilities for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           lensing.h
/// \brief          Integration utilities for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
/// (michel-andres.breton@obspm.fr) \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef LENSING_H_INCLUDED
#define LENSING_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
// Include libs
// Include project
#include "integrator.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/hypersphere.h"
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "utility.h"
// Octree
#ifdef VELOCITYFIELD
#include "gravity2.h" // with velocity slots
#else
#include "gravity.h"
#endif
#include "cone.h"
#include "output.h"
#include "photon.h"
// Misc
// -------------------------------------------------------------------------- //

// --------------------------------- CLASS ---------------------------------- //
// Integration utilities for raytracing
/// \brief          Integration utilities for raytracing.
/// \details        Provides a list of integration routines for geodesics
///                 integration.
class Lensing {
    // Initialization
    /// \name           Initialization
    //@{
public:
    // Jacobian matrices
    template <int Order = ORDER, bool RK4 = true, bool Verbose = false,
              class Parameter, class Point, class Cosmology, class Octree,
              class Type>
    static std::array<std::array<double, 2>, 2>
    dbetadtheta(const Parameter &parameters, const Point &kiTarget,
                const Type interpRef, const Point &observer, const Type phi,
                const Type theta, const Type dist, const Cosmology &cosmology,
                const Octree &octree, const Point &vobs, const Type length);
    template <int Order = ORDER, bool RK4 = true, bool Verbose = false,
              class Parameter, class Point, class Cosmology, class Octree,
              class Type>
    static std::vector<std::array<std::array<double, 2>, 2>>
    dbetadtheta(const Parameter &parameters, const std::vector<Point> &kiTargets,
                const std::vector<Type> &interpRef, const Point &observer,
                const Type phi, const Type theta, const std::vector<Type> &dist,
                const Cosmology &cosmology, const Octree &octree,
                const Point &vobs, const Type length);

    template <
        int Order = ORDER, bool RK4 = true, bool Verbose = false, class Type,
        class Trajectory,
        template <typename Kind, class Index, class Data, unsigned int Dimension,
                  class Position, class Extent, class Element, class Container>
        class Octree,
        typename Kind, class Index, class Data, unsigned int Dimension,
        class Position, class Extent, class Element, class Container>
    static std::array<std::array<double, 2>, 2>
    dbetadtheta_infinitesimal(const Type distance, const Trajectory &trajectory,
                              const Octree<Kind, Index, Data, Dimension, Position,
                                           Extent, Element, Container> &octree,
                              const Type length);
    template <
        int Order = ORDER, bool RK4 = true, bool Verbose = false, class Type,
        class Trajectory,
        template <typename Kind, class Index, class Data, unsigned int Dimension,
                  class Position, class Extent, class Element, class Container>
        class Octree,
        typename Kind, class Index, class Data, unsigned int Dimension,
        class Position, class Extent, class Element, class Container>
    static std::vector<std::array<std::array<double, 2>, 2>>
    dbetadtheta_infinitesimal(const std::vector<Type> &dist,
                              const Trajectory &trajectory,
                              const Octree<Kind, Index, Data, Dimension, Position,
                                           Extent, Element, Container> &octree,
                              const Type length);

    // Flexion
    template <int Order = ORDER, bool RK4 = true, bool Verbose = false,
              class Parameter, class Point, class Cosmology, class Octree,
              class Type>
    static std::array<double, 6>
    flexion(const Parameter &parameters, const Point &central_position,
            const Point &kiTarget, const Type interpRef, const Point &observer,
            const Type phi, const Type theta, const Type dist,
            const Cosmology &cosmology, const Octree &octree, const Point &vobs,
            const Type length);
    template <int Order = ORDER, bool RK4 = true, bool Verbose = false,
              class Parameter, class Point, class Cosmology, class Octree,
              class Type>
    static std::vector<std::array<double, 6>>
    flexion(const Parameter &parameters, const std::vector<Point> &central_positions,
            const std::vector<Point> &kiTargets,
            const std::vector<Type> &interpRef, const Point &observer,
            const Type phi, const Type theta, const std::vector<Type> &dist,
            const Cosmology &cosmology, const Octree &octree, const Point &vobs,
            const Type length);

    //@}

    // Test
    /// \name           Test
    //@{
public:
    static int example();
    //@}
};
// -------------------------------------------------------------------------- //

//  Jacobian calculation for angles deformations
/// \brief          Compute the jacobian for lensing
/// \details        Compute the jacobian for lensing
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2
/// for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Parameter Parameter type.
/// \tparam	    Point point type.
/// \tparam	    Cosmology cosmology type
/// \tparam         Octree octree type
/// \tparam         Type Scalar type.
/// \param[in]      parameters Parameter structure.
/// \param[in]      kiTarget Vector normal to the plane needed to compute the
/// distortion matrix. \param[in]      interpRef value for interpolation of the
/// bundle. \param[in]      observer Point observer position \param[in] phiInit
/// initial trajectory angle \param[in]      thetaInit initial trajectory angle
/// \param[in]      dist distance at different evaluations of the matrix
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \return         2x2 array Aij giving the jacobian from seen angles to true
/// angles
template <int Order, bool RK4, bool Verbose, class Parameter, class Point,
          class Cosmology, class Octree, class Type>
std::array<std::array<double, 2>, 2>
Lensing::dbetadtheta(const Parameter &parameters, const Point &kiTarget,
                     const Type interpRef, const Point &observer,
                     const Type phiInit, const Type thetaInit, const Type dist,
                     const Cosmology &cosmology, const Octree &octree,
                     const Point &vobs, const Type length) {

    std::vector<Point> kiTargets{kiTarget};
    std::vector<Type> interpRefvec{interpRef};
    std::vector<Type> distance{dist};
    return Lensing::dbetadtheta<Order, RK4, Verbose>(
        parameters, kiTargets, interpRefvec, observer, phiInit, thetaInit,
        distance, cosmology, octree, vobs, length)[0];
}

//  Jacobian calculation for angles deformations for multiple redshifts
/// \brief          Compute the jacobian for lensing
/// \details        Compute the jacobian for lensing
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2
/// for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Parameter Parameter type.
/// \tparam	    Point point type.
/// \tparam	    Cosmology cosmology type
/// \tparam         Octree octree type
/// \tparam         Type Scalar type.
/// \param[in]      parameters Parameter structure.
/// \param[in]      kiTargets Vector of vectors normal to the plane needed
///		    to compute the distortion matrix at some given surface.
/// \param[in]      interpRefvec values for interpolation of the bundle.
/// \param[in]      observer Point observer position
/// \param[in]      phiInit initial trajectory angle
/// \param[in]      thetaInit initial trajectory angle
/// \param[in]      dist distance at different evaluations of the matrix
/// \param[in]      interpolation String of a given interpolation
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \return         2x2 array Aij giving the jacobian from seen angles to true
/// angles
template <int Order, bool RK4, bool Verbose, class Parameter, class Point,
          class Cosmology, class Octree, class Type>
std::vector<std::array<std::array<double, 2>, 2>> Lensing::dbetadtheta(
    const Parameter &parameters, const std::vector<Point> &kiTargets,
    const std::vector<Type> &interpRefvec, const Point &observer,
    const Type phiInit, const Type thetaInit, const std::vector<Type> &dist,
    const Cosmology &cosmology, const Octree &octree, const Point &vobs,
    const Type length) {

    const unsigned int size = interpRefvec.size();
    std::vector<std::array<std::array<double, 2>, 2>> jacobian(size);
    std::array<Photon<double, 3>, 4> photons;
    std::array<magrathea::Evolution<Photon<double, 3>>, 4> trajectories;
    std::vector<std::array<Point, 4>> bundle_position(size);
    const double cp(std::cos(phiInit)), sp(std::sin(phiInit)),
        ct(std::cos(thetaInit)), st(std::sin(thetaInit));
    const uint n_bundle = 4;

    // Initialization
    std::vector<Point> rini(n_bundle);
    Point target, e1, e2;
    // Direction of target
    target = {cp * st, sp * st, ct};
    // screen perpendicular to the target direction
    e1 = {-sp, cp, 0};
    e2 = {-cp * ct, -sp * ct, st};
    const double to = std::tan(parameters.openingmin);

    // Initialise direction of bundle photons
    /*  2

    0   C   1

        3 */
    rini[0] = {target[0] + to * e1[0], target[1] + to * e1[1],
               target[2]};
    rini[1] = {target[0] - to * e1[0],
               target[1] - to * e1[1], target[2]};
    rini[2] = {target[0] - to * e2[0], target[1] - to * e2[1],
               target[2] - to * e2[2]};
    rini[3] = {target[0] + to * e2[0],
               target[1] + to * e2[1],
               target[2] + to * e2[2]};

    // Launch photons
    for (unsigned int i = 0; i < n_bundle; i++) {
        photons[i] = Integrator::launch(observer[0], observer[1], observer[2],
                                        rini[i][0], rini[i][1], rini[i][2]);
    }

    // Integration
    for (unsigned int itrajectory = 0; itrajectory < n_bundle; itrajectory++) {
        trajectories[itrajectory].append(photons[itrajectory]);
        Integrator::integrate<Order, RK4, Verbose>(
            trajectories[itrajectory], parameters.stop_bundle, interpRefvec.back(),
            cosmology, octree, vobs, length, parameters.nsteps, kiTargets.back());

        for (uint iref = 0; iref < size; iref++) {
            unsigned int marked(0), firstid(0);
            double next(0), previous(0);
            // Interpolate bundle photons at some parameter
            if (parameters.stop_bundle == "lambda") {
                if (trajectories[itrajectory].back().lambda() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        jacobian[iref2][0][0] = 42;
                        jacobian[iref2][0][1] = 42;
                        jacobian[iref2][1][0] = 1;
                        jacobian[iref2][1][1] = 0;
                    }
                    return jacobian;
                }
                photons[itrajectory].lambda() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.lambda() < second.lambda();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].lambda();
                previous = trajectories[itrajectory][firstid].lambda();
            } else if ((parameters.stop_bundle == "r") ||
                       (parameters.stop_bundle == "radius")) {
                if (trajectories[itrajectory].back().chi() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        jacobian[iref2][0][0] = 42;
                        jacobian[iref2][0][1] = 42;
                        jacobian[iref2][1][0] = 1;
                        jacobian[iref2][1][1] = 0;
                    }
                    return jacobian;
                }
                photons[itrajectory].chi() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.chi() < second.chi();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].chi();
                previous = trajectories[itrajectory][firstid].chi();
            } else if (parameters.stop_bundle == "redshift") {
                if (trajectories[itrajectory].back().redshift() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        jacobian[iref2][0][0] = 42;
                        jacobian[iref2][0][1] = 42;
                        jacobian[iref2][1][0] = 1;
                        jacobian[iref2][1][1] = 0;
                    }
                    return jacobian;
                }
                photons[itrajectory].redshift() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.redshift() < second.redshift();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].redshift();
                previous = trajectories[itrajectory][firstid].redshift();
            } else if (parameters.stop_bundle == "a") {
                if (trajectories[itrajectory].back().a() > interpRefvec[iref]) {
                    jacobian[iref][0][0] = 42;
                    jacobian[iref][0][1] = 42;
                    jacobian[iref][1][0] = 1;
                    jacobian[iref][1][1] = 3;
                    return jacobian;
                }
                photons[itrajectory].a() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.a() > second.a();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].a();
                previous = trajectories[itrajectory][firstid].a();
            } else if ((parameters.stop_bundle == "t") ||
                       (parameters.stop_bundle == "eta")) {
                if (trajectories[itrajectory].back().t() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        jacobian[iref2][0][0] = 42;
                        jacobian[iref2][0][1] = 42;
                        jacobian[iref2][1][0] = 1;
                        jacobian[iref2][1][1] = 0;
                    }
                    return jacobian;
                }
                photons[itrajectory].t() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.t() < second.t();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].t();
                previous = trajectories[itrajectory][firstid].t();
            } else if (parameters.stop_bundle == "plane") {
                if ((trajectories[itrajectory].back().x() * kiTargets[iref][0] +
                     trajectories[itrajectory].back().y() * kiTargets[iref][1] +
                     trajectories[itrajectory].back().z() * kiTargets[iref][2]) <
                    interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        jacobian[iref2][0][0] = 42;
                        jacobian[iref2][0][1] = 42;
                        jacobian[iref2][1][0] = 1;
                        jacobian[iref2][1][1] = 0;
                    }
                    return jacobian;
                }
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     interpRefvec[iref],
                                     [=, &kiTargets](const double &val,
                                                     const Photon<double, 3> &first) {
                                         return val < (first.x() * kiTargets[iref][0] +
                                                       first.y() * kiTargets[iref][1] +
                                                       first.z() * kiTargets[iref][2]);
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].x() * kiTargets[iref][0] +
                       trajectories[itrajectory][firstid + 1].y() * kiTargets[iref][1] +
                       trajectories[itrajectory][firstid + 1].z() * kiTargets[iref][2];
                previous = trajectories[itrajectory][firstid].x() * kiTargets[iref][0] +
                           trajectories[itrajectory][firstid].y() * kiTargets[iref][1] +
                           trajectories[itrajectory][firstid].z() * kiTargets[iref][2];
            } else {
                std::cout
                    << "# WARNING: Please choose an existing stop_bundle parameter"
                    << std::endl;
                std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                          << std::endl;
                std::terminate();
            }
            const double f = (next - interpRefvec[iref]) / (next - previous);
            // Get position of bundle photons at some surface
            bundle_position[iref][itrajectory][0] =
                trajectories[itrajectory][firstid].x() * f +
                trajectories[itrajectory][firstid + 1].x() * (1 - f);
            bundle_position[iref][itrajectory][1] =
                trajectories[itrajectory][firstid].y() * f +
                trajectories[itrajectory][firstid + 1].y() * (1 - f);
            bundle_position[iref][itrajectory][2] =
                trajectories[itrajectory][firstid].z() * f +
                trajectories[itrajectory][firstid + 1].z() * (1 - f);
        } // iref for different interp
    }     // itrajectory for

    // Compute jacobian
    for (uint iref = 0; iref < size; iref++) {
        Point rh, rv;
        // Compute finite differences
        for (unsigned int idim = 0; idim < 3; idim++) {
            rh[idim] = bundle_position[iref][0][idim] - bundle_position[iref][1][idim];
            rv[idim] = bundle_position[iref][2][idim] - bundle_position[iref][3][idim];
        }
        // Need normal to screen (either normal to the photon or normal to the
        // comoving direction of the source)
        double kx = kiTargets[iref][0];
        double ky = kiTargets[iref][1];
        double kz = kiTargets[iref][2];
        const double phik = std::atan2(ky, kx);
        const double thetak =
            std::acos(kz / std::sqrt(kx * kx + ky * ky + kz * kz));
        const double cpk(std::cos(phik)), spk(std::sin(phik)),
            ctk(std::cos(thetak)), stk(std::sin(thetak));
        // screen perpendicular to the target direction
        e1 = {-spk, cpk, 0};
        e2 = {-cpk * ctk, -spk * ctk, stk};
        // Compute lensing jacobian matrix
        const double inv_denominator = 0.5 / (dist[iref] * to);
        jacobian[iref][0][0] = std::inner_product(std::begin(rh), std::end(rh), std::begin(e1), 0.) * inv_denominator;
        jacobian[iref][0][1] = std::inner_product(std::begin(rv), std::end(rv), std::begin(e1), 0.) * inv_denominator;
        jacobian[iref][1][0] = -std::inner_product(std::begin(rh), std::end(rh), std::begin(e2), 0.) * inv_denominator;
        jacobian[iref][1][1] = -std::inner_product(std::begin(rv), std::end(rv), std::begin(e2), 0.) * inv_denominator;
    }

    return jacobian;
}

//  Jacobian calculation for angles deformations for one trajectory
/// \brief          Compute the jacobian for lensing
/// \details        Compute the jacobian for lensing
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2
/// for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Type Scalar type.
/// \tparam	    Trajectory trajectory type.
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Kind Kind type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in]      distance scalar comoving distance to stop lensing
/// integration \param[in]      phi initial trajectory angle \param[in] theta
/// initial trajectory angle \param[in]      trajectory Photon trajectory used
/// to compute the jacobian matrix. \param[in]      octree Octree. \param[in]
/// length Spatial length in SI units. \return         2x2 array Aij giving the
/// jacobian from seen angles to true angles
template <
    int Order, bool RK4, bool Verbose, class Type, class Trajectory,
    template <typename Kind, class Index, class Data, unsigned int Dimension,
              class Position, class Extent, class Element, class Container>
    class Octree,
    typename Kind, class Index, class Data, unsigned int Dimension,
    class Position, class Extent, class Element, class Container>
std::array<std::array<double, 2>, 2> Lensing::dbetadtheta_infinitesimal(
    const Type distance, const Trajectory &trajectory,
    const Octree<Kind, Index, Data, Dimension, Position, Extent, Element,
                 Container> &octree,
    const Type length) {

    std::vector<Type> dist{distance};
    return Lensing::dbetadtheta_infinitesimal<Order, RK4, Verbose>(
        dist, trajectory, octree, length)[0];
}
//  Jacobian calculation for angles deformations for one trajectory
/// \brief          Compute the jacobian for lensing
/// \details        Compute the jacobian for lensing
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2
/// for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Type Scalar type.
/// \tparam	    Trajectory trajectory type.
/// \tparam         Octree octree type
/// \tparam         Kind Kind type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in]      dist vector of comoving distances to stop lensing
/// integrations \param[in]      phi initial trajectory angle \param[in] theta
/// initial trajectory angle \param[in]      trajectory Photon trajectory used
/// to compute the jacobian matrix. \param[in]      octree Octree. \param[in]
/// length Spatial length in SI units. \return         2x2 array Aij giving the
/// jacobian from seen angles to true angles
template <
    int Order, bool RK4, bool Verbose, class Type, class Trajectory,
    template <typename Kind, class Index, class Data, unsigned int Dimension,
              class Position, class Extent, class Element, class Container>
    class Octree,
    typename Kind, class Index, class Data, unsigned int Dimension,
    class Position, class Extent, class Element, class Container>
std::vector<std::array<std::array<double, 2>, 2>>
Lensing::dbetadtheta_infinitesimal(
    const std::vector<Type> &dist, const Trajectory &trajectory,
    const Octree<Kind, Index, Data, Dimension, Position, Extent, Element,
                 Container> &octree,
    const Type length) {

    unsigned int size = dist.size();
    std::vector<std::array<std::array<double, 2>, 2>> result(size);
    static const double c2 = magrathea::Constants<double>::c2();
    double dxxp(0), dxyp(0), dxzp(0), dyyp(0), dyzp(0), dzzp(0);
    double d1d1p(0), d2d2p(0), d1d2p(0);
    Data dataxm1, dataxp1, dataym1, datayp1, datazm1, datazp1;
    unsigned int levelshift = std::log2(EXTENT);
    double dl(0), x(0), y(0), z(0), r(0), dr(0), gxsx(0);
    std::vector<Element> elemsTsc(27);

    std::vector<std::vector<double>> a11(size), a12(size), a22(size);

    for (unsigned int j = 0; j < size; j++) {
        a11[j].resize(trajectory.size());
        a12[j].resize(trajectory.size());
        a22[j].resize(trajectory.size());
    }
    for (unsigned int i = 1; i < trajectory.size(); i++) {
        // Get position and half-cell distance
        x = trajectory[i].x();
        y = trajectory[i].y();
        z = trajectory[i].z();
        r = trajectory[i].chi();
        const double phi = std::atan2(y, x);
        const double theta = std::acos(z / r);
        const double cp(std::cos(phi)), sp(std::sin(phi)), ct(std::cos(theta)),
            st(std::sin(theta));
        const double s2p(std::sin(2 * phi)), s2t(std::sin(2 * theta));
        // Derivation step (arbitrary)
        dl = std::pow(2, -(trajectory[i].level() - 1 - levelshift));
        // Compute force and derivative
        dataxm1 = octree.tsc(elemsTsc, x - 0.5 * dl, y, z);
        dataxp1 = octree.tsc(elemsTsc, x + 0.5 * dl, y, z);
        dataym1 = octree.tsc(elemsTsc, x, y - 0.5 * dl, z);
        datayp1 = octree.tsc(elemsTsc, x, y + 0.5 * dl, z);
        datazm1 = octree.tsc(elemsTsc, x, y, z - 0.5 * dl);
        datazp1 = octree.tsc(elemsTsc, x, y, z + 0.5 * dl);
        dxxp = (dataxp1.dphidx() - dataxm1.dphidx()) / dl;
        dxyp = (datayp1.dphidx() - dataym1.dphidx()) / dl;
        dxzp = (datazp1.dphidx() - datazm1.dphidx()) / dl;
        dyyp = (datayp1.dphidy() - dataym1.dphidy()) / dl;
        dyzp = (datazp1.dphidy() - datazm1.dphidy()) / dl;
        dzzp = (datazp1.dphidz() - datazm1.dphidz()) / dl;
        // Convert to angular derivative (see Barreira et al. 2016, eq 39-41). d1d2p
        // corrected : need to go from orthonormal vector frame to coordinate frame.
        d1d1p = sp * sp * dxxp + cp * cp * dyyp - s2p * dxyp;
        d2d2p = cp * cp * ct * ct * dxxp + sp * sp * ct * ct * dyyp +
                st * st * dzzp + s2p * ct * ct * dxyp - sp * s2t * dyzp -
                cp * s2t * dxzp;
        d1d2p = ct * cp * sp * (dyyp - dxxp) + (cp * cp - sp * sp) * ct * dxyp +
                sp * st * dxzp - cp * st * dyzp;
        // Finalize computation
        dr = trajectory[i].chi() - trajectory[i - 1].chi();
        for (unsigned int j = 0; j < size; j++) {
            gxsx = (dist[j] - trajectory[i].chi()) * trajectory[i].chi() / dist[j];
            a11[j][i] = a11[j][i - 1] + gxsx * d1d1p * dr;
            a12[j][i] = a12[j][i - 1] + gxsx * d1d2p * dr;
            a22[j][i] = a22[j][i - 1] + gxsx * d2d2p * dr;
        }
    }

    // Sum and finalize
    unsigned int marked(0), firstid(0);
    Trajectory photon;
    photon.resize(1);
    for (unsigned int j = 0; j < size; j++) {
        photon[0].chi() = dist[j];
        marked =
            std::distance(std::begin(trajectory),
                          std::upper_bound(std::begin(trajectory),
                                           std::end(trajectory), photon[0],
                                           [](const Photon<double, 3> &first,
                                              const Photon<double, 3> &second) {
                                               return first.chi() < second.chi();
                                           }));
        firstid = marked - (marked > 0);
        result[j][0][0] = 1 - 2 * a11[j][firstid] / c2 * length;
        // We assume a12 = a21
        result[j][0][1] = result[j][1][0] = -2 * a12[j][firstid] / c2 * length;
        result[j][1][1] = 1 - 2 * a22[j][firstid] / c2 * length;
    }
    return result;
}

//  Hessian calculation for angles deformations
/// \brief          Compute the flexion for lensing
/// \details        Compute the flexion for lensing
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2
/// for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Parameter Parameter type.
/// \tparam	    Point point type.
/// \tparam	    Cosmology cosmology type
/// \tparam         Octree octree type
/// \tparam         Type Scalar type.
/// \param[in]      parameters Parameter structure.
/// \param[in]      kiTarget Vector normal to the plane needed to compute the
/// distortion matrix. \param[in]      interpRef value for interpolation of the
/// bundle. \param[in]      observer Point observer position \param[in] phiInit
/// initial trajectory angle \param[in]      thetaInit initial trajectory angle
/// \param[in]      dist distance at different evaluations of the matrix
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \return         1d array Dijk giving the Hessian from seen angles to true
/// angles

template <int Order, bool RK4, bool Verbose, class Parameter, class Point,
          class Cosmology, class Octree, class Type>
std::array<double, 6> Lensing::flexion(
    const Parameter &parameters, const Point &central_position, const Point &kiTarget,
    const Type interpRef, const Point &observer, const Type phiInit,
    const Type thetaInit, const Type dist, const Cosmology &cosmology,
    const Octree &octree, const Point &vobs, const Type length) {

    std::vector<Point> central_positions{central_position};
    std::vector<Point> kiTargets{kiTarget};
    std::vector<Type> interpRefvec{interpRef};
    std::vector<Type> distance{dist};
    return Lensing::flexion<Order, RK4, Verbose>(
        parameters, central_positions, kiTargets, interpRefvec, observer, phiInit,
        thetaInit, distance, cosmology, octree, vobs, length)[0];
}

//  Hessian calculation for angles deformations for multiple redshifts
/// \brief          Compute the flexion for lensing
/// \details        Compute the flexion for lensing
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2
/// for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Parameter Parameter type.
/// \tparam	    Point point type.
/// \tparam	    Cosmology cosmology type
/// \tparam         Octree octree type
/// \tparam         Type Scalar type.
/// \param[in]      parameters Parameter structure.
/// \param[in]      kiTargets Vector of vectors normal to the plane needed
///		    to compute the distortion matrix at some given surface.
/// \param[in]      interpRefvec values for interpolation of the bundle.
/// \param[in]      observer Point observer position
/// \param[in]      phiInit initial trajectory angle
/// \param[in]      thetaInit initial trajectory angle
/// \param[in]      dist distance at different evaluations of the matrix
/// \param[in]      interpolation String of a given interpolation
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \return         1d array Dijk giving the Hessian from seen angles to true
/// angles
template <int Order, bool RK4, bool Verbose, class Parameter, class Point,
          class Cosmology, class Octree, class Type>
std::vector<std::array<double, 6>> Lensing::flexion(
    const Parameter &parameters, const std::vector<Point> &central_positions,
    const std::vector<Point> &kiTargets, const std::vector<Type> &interpRefvec,
    const Point &observer, const Type phiInit, const Type thetaInit,
    const std::vector<Type> &dist, const Cosmology &cosmology,
    const Octree &octree, const Point &vobs, const Type length) {

    const unsigned int size = interpRefvec.size();
    std::vector<std::array<double, 6>> hessian(size);
    std::array<Photon<double, 3>, 8> photons;
    std::array<magrathea::Evolution<Photon<double, 3>>, 8> trajectories;
    std::vector<std::array<Point, 8>> bundle_position(size);
    const uint n_bundle = 8;
    const double cp(std::cos(phiInit)), sp(std::sin(phiInit)),
        ct(std::cos(thetaInit)), st(std::sin(thetaInit));

    // Initialization
    std::vector<Point> rini(n_bundle);
    Point target, e1, e2;
    target = {cp * st, sp * st, ct};
    // screen perpendicular to the target direction
    e1 = {-sp, cp, 0};
    e2 = {-cp * ct, -sp * ct, st};
    const double to = std::tan(parameters.openingmin);

    // Initialise direction of bundle photons

    /*  4   2   6

        0   C   1

        5   3   7 */

    rini[0] = {target[0] + to * e1[0], target[1] + to * e1[1], target[2]};
    rini[1] = {target[0] - to * e1[0], target[1] - to * e1[1], target[2]};
    rini[2] = {target[0] - to * e2[0], target[1] - to * e2[1],
               target[2] - to * e2[2]};
    rini[3] = {target[0] + to * e2[0], target[1] + to * e2[1],
               target[2] + to * e2[2]};
    rini[4] = {target[0] + to * e1[0] - to * e2[0],
               target[1] + to * e1[1] - to * e2[1], target[2] - to * e2[2]};
    rini[5] = {target[0] + to * e1[0] + to * e2[0],
               target[1] + to * e1[1] + to * e2[1], target[2] + to * e2[2]};
    rini[6] = {target[0] - to * e1[0] - to * e2[0],
               target[1] - to * e1[1] - to * e2[1], target[2] - to * e2[2]};
    rini[7] = {target[0] - to * e1[0] + to * e2[0],
               target[1] - to * e1[1] + to * e2[1], target[2] + to * e2[2]};
    // Launch photons
    for (unsigned int i = 0; i < n_bundle; i++) {
        photons[i] = Integrator::launch(observer[0], observer[1], observer[2],
                                        rini[i][0], rini[i][1], rini[i][2]);
    }

    // Integration
    for (unsigned int itrajectory = 0; itrajectory < n_bundle; itrajectory++) {
        trajectories[itrajectory].append(photons[itrajectory]);
        Integrator::integrate<Order, RK4, Verbose>(
            trajectories[itrajectory], parameters.stop_bundle, interpRefvec.back(),
            cosmology, octree, vobs, length, parameters.nsteps, kiTargets.back());

        for (uint iref = 0; iref < size; iref++) {
            unsigned int marked(0), firstid(0);
            double next(0), previous(0);
            // Interpolate bundle photons at some parameter
            if (parameters.stop_bundle == "lambda") {
                if (trajectories[itrajectory].back().lambda() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        hessian[iref2][0] = 42;
                    }
                    return hessian;
                }
                photons[itrajectory].lambda() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.lambda() < second.lambda();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].lambda();
                previous = trajectories[itrajectory][firstid].lambda();
            } else if ((parameters.stop_bundle == "r") ||
                       (parameters.stop_bundle == "radius")) {
                if (trajectories[itrajectory].back().chi() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        hessian[iref2][0] = 42;
                    }
                    return hessian;
                }
                photons[itrajectory].chi() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.chi() < second.chi();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].chi();
                previous = trajectories[itrajectory][firstid].chi();
            } else if (parameters.stop_bundle == "redshift") {
                if (trajectories[itrajectory].back().redshift() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        hessian[iref2][0] = 42;
                    }
                    return hessian;
                }
                photons[itrajectory].redshift() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.redshift() < second.redshift();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].redshift();
                previous = trajectories[itrajectory][firstid].redshift();
            } else if (parameters.stop_bundle == "a") {
                if (trajectories[itrajectory].back().a() > interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        hessian[iref2][0] = 42;
                    }
                    return hessian;
                }
                photons[itrajectory].a() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.a() > second.a();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].a();
                previous = trajectories[itrajectory][firstid].a();
            } else if ((parameters.stop_bundle == "t") ||
                       (parameters.stop_bundle == "eta")) {
                if (trajectories[itrajectory].back().t() < interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        hessian[iref2][0] = 42;
                    }
                    return hessian;
                }
                photons[itrajectory].t() = interpRefvec[iref];
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     photons[itrajectory],
                                     [](const Photon<double, 3> &first,
                                        const Photon<double, 3> &second) {
                                         return first.t() < second.t();
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].t();
                previous = trajectories[itrajectory][firstid].t();
            } else if (parameters.stop_bundle == "plane") {
                if ((trajectories[itrajectory].back().x() * kiTargets[iref][0] +
                     trajectories[itrajectory].back().y() * kiTargets[iref][1] +
                     trajectories[itrajectory].back().z() * kiTargets[iref][2]) <
                    interpRefvec[iref]) {
                    for (uint iref2 = 0; iref2 < size; iref2++) {
                        hessian[iref2][0] = 42;
                    }
                    return hessian;
                }
                marked = std::distance(
                    std::begin(trajectories[itrajectory]),
                    std::upper_bound(std::begin(trajectories[itrajectory]),
                                     std::end(trajectories[itrajectory]),
                                     interpRefvec[iref],
                                     [=, &kiTargets](const double &val,
                                                     const Photon<double, 3> &first) {
                                         return val < (first.x() * kiTargets[iref][0] +
                                                       first.y() * kiTargets[iref][1] +
                                                       first.z() * kiTargets[iref][2]);
                                     }));
                firstid = marked - (marked > 0);
                next = trajectories[itrajectory][firstid + 1].x() * kiTargets[iref][0] +
                       trajectories[itrajectory][firstid + 1].y() * kiTargets[iref][1] +
                       trajectories[itrajectory][firstid + 1].z() * kiTargets[iref][2];
                previous = trajectories[itrajectory][firstid].x() * kiTargets[iref][0] +
                           trajectories[itrajectory][firstid].y() * kiTargets[iref][1] +
                           trajectories[itrajectory][firstid].z() * kiTargets[iref][2];
            } else {
                std::cout
                    << "# WARNING: Please choose an existing stop_bundle parameter"
                    << std::endl;
                std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                          << std::endl;
                std::terminate();
            }
            const double f = (next - interpRefvec[iref]) / (next - previous);
            // Get position of bundle photons at some surface
            bundle_position[iref][itrajectory][0] =
                trajectories[itrajectory][firstid].x() * f +
                trajectories[itrajectory][firstid + 1].x() * (1 - f);
            bundle_position[iref][itrajectory][1] =
                trajectories[itrajectory][firstid].y() * f +
                trajectories[itrajectory][firstid + 1].y() * (1 - f);
            bundle_position[iref][itrajectory][2] =
                trajectories[itrajectory][firstid].z() * f +
                trajectories[itrajectory][firstid + 1].z() * (1 - f);
        } // iref for different interp
    }     // itrajectory for

    // Compute Hessian
    for (uint iref = 0; iref < size; iref++) {
        Point rhh, rvv, rhv;
        // Compute finite differences
        for (unsigned int idim = 0; idim < 3; idim++) {
            rhh[idim] = bundle_position[iref][0][idim] + bundle_position[iref][1][idim] - 2 * central_positions[iref][idim];
            rvv[idim] = bundle_position[iref][2][idim] + bundle_position[iref][3][idim] - 2 * central_positions[iref][idim];
            rhv[idim] = (bundle_position[iref][4][idim] + bundle_position[iref][7][idim] - bundle_position[iref][5][idim] - bundle_position[iref][6][idim]);
        }

        // Need normal to screen (either normal to the photon or normal to the
        // comoving direction of the source)
        double kx = kiTargets[iref][0];
        double ky = kiTargets[iref][1];
        double kz = kiTargets[iref][2];
        const double phik = std::atan2(ky, kx);
        const double thetak =
            std::acos(kz / std::sqrt(kx * kx + ky * ky + kz * kz));
        const double cpk(std::cos(phik)), spk(std::sin(phik)),
            ctk(std::cos(thetak)), stk(std::sin(thetak));
        // screen perpendicular to the target direction
        e1 = {-spk, cpk, 0};
        e2 = {-cpk * ctk, -spk * ctk, stk};
        // Compute lensing jacobian matrix
        const double inv_denominator = 1. / (dist[iref] * to * to); // Only 1 dist[ref] so that dimension of flexion is rad^-1
        hessian[iref][0] = std::inner_product(std::begin(rhh), std::end(rhh), std::begin(e1), 0.) * inv_denominator;
        hessian[iref][1] = std::inner_product(std::begin(rhv), std::end(rhv), std::begin(e1), 0.) * inv_denominator;
        hessian[iref][2] = std::inner_product(std::begin(rvv), std::end(rvv), std::begin(e1), 0.) * inv_denominator;
        hessian[iref][3] = std::inner_product(std::begin(rhh), std::end(rhh), std::begin(e2), 0.) * inv_denominator;
        hessian[iref][4] = std::inner_product(std::begin(rhv), std::end(rhv), std::begin(e2), 0.) * inv_denominator;
        hessian[iref][5] = std::inner_product(std::begin(rvv), std::end(rvv), std::begin(e2), 0.) * inv_denominator;
    }

    return hessian;
}

// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
#endif // LENSING_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
