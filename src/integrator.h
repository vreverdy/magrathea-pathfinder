/* ******************************* INTEGRATOR ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Integrator
// DESCRIPTION :    Integration utilities for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           integrator.h
/// \brief          Integration utilities for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef INTEGRATOR_H_INCLUDED
#define INTEGRATOR_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <type_traits>
#include <limits>
#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <random>
// Include libs
// Include project
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/constants.h"
#include "magrathea/hypersphere.h"
#include "magrathea/evolution.h"
#include "utility.h"
// Octree
#ifdef VELOCITYFIELD
#include "gravity2.h" // with velocity slots
#else
#include "gravity.h"
#endif
#include "photon.h"
#include "cone.h"
#include "output.h"
// Misc
// -------------------------------------------------------------------------- //



// --------------------------------- CLASS ---------------------------------- //
// Integration utilities for raytracing
/// \brief          Integration utilities for raytracing.
/// \details        Provides a list of integration routines for geodesics 
///                 integration.
class Integrator final
{  
    // Initialization
    /// \name           Initialization
    //@{
    public:
        template <typename Type = double, class Sphere, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension = Sphere::dimension(), class = typename std::enable_if<Dimension == 3>::type> static Photon<Type, Dimension> launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, Engine& engine, Distribution& distribution);
        template <typename Type = double, class Sphere, class Container, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension = Sphere::dimension(), class = typename std::enable_if<(Dimension == 3) && (std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, Cone<Vector, Scalar> >::value)>::type> static Photon<Type, Dimension> launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, const Container& cones, Engine& engine, Distribution& distribution);
        template <typename Type, unsigned int Dimension = 3, class = typename std::enable_if<Dimension == 3>::type> static Photon<Type, Dimension> launch(const Type xbegin, const Type ybegin, const Type zbegin, const Type xend, const Type yend, const Type zend);
	//MODIF
        template <typename Type, unsigned int Dimension = 3, class = typename std::enable_if<Dimension == 3>::type> static Photon<Type, Dimension> launch(const Type xbegin, const Type ybegin, const Type zbegin, const Type phi, const Type theta);
	// FIN MODIF
        template <bool Center = false, typename Type, unsigned int Dimension, class Container = std::vector<Photon<Type, Dimension> >, class = typename std::enable_if<Dimension == 3>::type> static Container launch(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation = Type());
    //@}
    
    // Computation
    /// \name           Computation
    //@{
    public:
        template <int Order = ORDER, class Array, class Element, class Cosmology, class Octree, class Type, unsigned int Dimension = Octree::dimension(), class Data = typename std::tuple_element<1, decltype(Octree::element())>::type, class Position = decltype(Octree::position()), class Extent = decltype(Octree::extent())> static Array& dphotondl(Array& output, const Array& input, const Cosmology& cosmology, const Octree& octree, const Type length, std::vector<Element>& elemsTsc);
    //@}

    // Evolution
    /// \name           Evolution
    //@{
    public:
        template <int Order = ORDER, bool RK4 = true, bool Verbose = false, class Cosmology, class Octree, class Type, class Trajectory, unsigned int Dimension = Octree::dimension(), class Element = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Trajectory>()[0])>::type>::type, class Data = typename std::tuple_element<1, decltype(Octree::element())>::type, class Core = decltype(Element::template type<1>()), unsigned int Size = std::tuple_size<Core>::value, class Position = decltype(Octree::position()), class Extent = decltype(Octree::extent()), class Point> static Trajectory& integrate(Trajectory& trajectory, const Cosmology& cosmology, const Octree& octree, const Point& vobs, const Type length, const unsigned int nsteps = 1);
        template <int Order = ORDER, bool RK4 = true, bool Verbose = false, class Cosmology, class Octree, class Type, class Trajectory, unsigned int Dimension = Octree::dimension(), class Element = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Trajectory>()[0])>::type>::type, class Data = typename std::tuple_element<1, decltype(Octree::element())>::type, class Core = decltype(Element::template type<1>()), unsigned int Size = std::tuple_size<Core>::value, class Position = decltype(Octree::position()), class Extent = decltype(Octree::extent()), class Point = std::array<Type, 3>> static Trajectory& integrate(Trajectory& trajectory, const std::string interpolation, const Type interpRef, const Cosmology& cosmology, const Octree& octree, const Point& vobs, const Type length, const unsigned int nsteps = 1, const Point& kiTarget = Point());

        template <int Order = ORDER, bool RK4 = true, bool Verbose = false, class Cosmology, class Octree, class Type, unsigned int Dimension, class Homogeneous = std::vector<Photon<Type, Dimension> >, class Point, class = typename std::enable_if<(Dimension == 3) && (Dimension == Octree::dimension())>::type> static magrathea::Evolution<Photon<Type, Dimension> > propagate(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation, const std::string& interpolation, const Cosmology& cosmology, const Octree& octree, const Point& vobs, const Type length, const unsigned int nsteps = 1, const Type amin = Type(), const std::string& filenames = std::string(), const Homogeneous& homogeneous = Homogeneous());
        

    //@}



    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ----------------------------- INITIALIZATION ----------------------------- //
// Launch a photon in a cone
/// \brief          Launch a photon in a cone.
/// \details        Launches a photon in the provided cone.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \tparam         Engine Random engine type.
/// \tparam         Distribution Random distribution type.
/// \param[in]      sphere Sphere to pick up the center and the surface.
/// \param[in]      cone Current cone.
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         Initial photon in the cone. 
template <typename Type, class Sphere, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension, class> 
Photon<Type, Dimension> Integrator::launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, Engine& engine, Distribution& distribution)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    static const unsigned int x = 0;
    static const unsigned int y = 1;
    static const unsigned int z = 2;
    Vector position = Vector();
    Photon<Type, Dimension> result;
    
    // Compute the photon direction
    do {
        position = sphere.template random<Dimension-1>(engine, distribution);
    } while (!cone.inside(position));
    
    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = sphere.center(x);
    result.y() = sphere.center(y);
    result.z() = sphere.center(z);
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = (position[x]-sphere.center(x))/sphere.radius();
    result.dydl() = (position[y]-sphere.center(y))/sphere.radius();
    result.dzdl() = (position[z]-sphere.center(z))/sphere.radius();

    // Finalization
    return result;
}

// Launch a photon in a serie of cones
/// \brief          Launch a photon in a serie of cones.
/// \details        Launches a photon in the provided cone with the guarantee
///                 that this cone is the closest one comparatively to a list
///                 of cones.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Sphere Sphere type.
/// \tparam         Container Container of cones type.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \tparam         Engine Random engine type.
/// \tparam         Distribution Random distribution type.
/// \param[in]      sphere Sphere to pick up the center and the surface.
/// \param[in]      cone Current cone.
/// \param[in]      cones List of cones to compare distance.
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         Initial photon in the cone. 
template <typename Type, class Sphere, class Container, class Vector, typename Scalar, class Engine, class Distribution, unsigned int Dimension, class> 
Photon<Type, Dimension> Integrator::launch(const Sphere& sphere, const Cone<Vector, Scalar>& cone, const Container& cones, Engine& engine, Distribution& distribution)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    static const unsigned int x = 0;
    static const unsigned int y = 1;
    static const unsigned int z = 2;
    const unsigned int size = cones.size();
    Vector position = Vector();
    Scalar length = Scalar();
    Scalar reference = Scalar();
    Scalar distance = Scalar();  
    bool ok = false;
    Photon<Type, Dimension> result;
    
    // Compute the photon direction
    do {
	length = Scalar();
	reference = Scalar();
        position = sphere.template random<Dimension-1>(engine, distribution);
        if (cone.inside(position)) {
            ok = true;
            for (unsigned int idim = 0; idim < Dimension; ++idim) {
                length += (cone.base(idim)-cone.vertex(idim))*(position[idim]-cone.vertex(idim));
            }
            length /= cone.template pow<2>(cone.length());
            for (unsigned int idim = 0; idim < Dimension; ++idim) {
                reference += cone.template pow<2>(position[idim]-(cone.vertex(idim)+(cone.base(idim)-cone.vertex(idim))*length));
            }
            for (unsigned int icone = 0; icone < size; ++icone) {
		length = Scalar();
		distance = Scalar();
                if (&cones[icone] != &cone) {
                    for (unsigned int idim = 0; idim < Dimension; ++idim) {
                        length += (cones[icone].base(idim)-cone.vertex(idim))*(position[idim]-cones[icone].vertex(idim));
                    }
		    if( !(length <  0)){
                        length /= cones[icone].template pow<2>(cones[icone].length());
                        for (unsigned int idim = 0; idim < Dimension; ++idim) {
                            distance += cones[icone].template pow<2>(position[idim]-(cones[icone].vertex(idim)+(cones[icone].base(idim)-cones[icone].vertex(idim))*length));
                        }
                        if (distance < reference) {
                            ok = false;
                            icone = size;
                        }
		    }
                }
            }
        }
    } while (!ok);
    
    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = sphere.center(x);
    result.y() = sphere.center(y);
    result.z() = sphere.center(z);
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = (position[x]-sphere.center(x))/sphere.radius();
    result.dydl() = (position[y]-sphere.center(y))/sphere.radius();
    result.dzdl() = (position[z]-sphere.center(z))/sphere.radius();

    // Finalization
    return result;
}

// Launch a photon going from a position to another
/// \brief          Launch a photon going from a position to another.
/// \details        Launches a photon starting from a point and going to another
///                 one.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \param[in]      xbegin Starting x coordinate.
/// \param[in]      ybegin Starting y coordinate.
/// \param[in]      zbegin Starting z coordinate.
/// \param[in]      xend Ending x coordinate.
/// \param[in]      yend Ending y coordinate.
/// \param[in]      zend Ending z coordinate.
/// \return         Initial photon.
template <typename Type, unsigned int Dimension, class>
Photon<Type, Dimension> Integrator::launch(const Type xbegin, const Type ybegin, const Type zbegin, const Type xend, const Type yend, const Type zend)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    const Type norm = std::sqrt((xend-xbegin)*(xend-xbegin)+(yend-ybegin)*(yend-ybegin)+(zend-zbegin)*(zend-zbegin)); 
    Photon<Type, Dimension> result;

    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = xbegin;
    result.y() = ybegin;
    result.z() = zbegin;
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = (xend-xbegin)/norm;
    result.dydl() = (yend-ybegin)/norm;
    result.dzdl() = (zend-zbegin)/norm;
    
    // Finalization
    return result;
}

// Launch a photon going from a position to a specified direction
/// \brief          Launch a photon going from a position to a specified direction.
/// \details        Launches a photon starting from a point and going to specified 
///                 direction.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \param[in]      xbegin Starting x coordinate.
/// \param[in]      ybegin Starting y coordinate.
/// \param[in]      zbegin Starting z coordinate.
/// \param[in]      phi Phi angular coordinate.
/// \param[in]      theta Theta angular coordinate.
/// \return         Initial photon.
template <typename Type, unsigned int Dimension, class>
Photon<Type, Dimension> Integrator::launch(const Type xbegin, const Type ybegin, const Type zbegin, const Type phi, const Type theta)
{
    // Initialization
    static const Type zero = Type();
    static const Type one = Type(1);
    Photon<Type, Dimension> result;

    // Set the photon cosmology
    result.a() = one;
    
    // Set the photon position
    result.t() = zero;
    result.x() = xbegin;
    result.y() = ybegin;
    result.z() = zbegin;
    
    // Set the photon direction
    result.dtdl() = one;
    result.dxdl() = std::cos(phi)*std::sin(theta);
    result.dydl() = std::sin(phi)*std::sin(theta);
    result.dzdl() = std::cos(theta);

    // Finalization
    return result;
}

// Launch photons around one photon
/// \brief          Launch photons around one photon.
/// \details        Launches a group of photon on a cone around one photon.
/// \tparam         Center Adds the center photon in the resulting vector if
///                 true.
/// \tparam         Type Photon type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Container Container of photons type.
/// \param[in]      photon Central photon.
/// \param[in]      count Number of photons to return.
/// \param[in]      angle Half-angle at the cone vertex.
/// \param[in]      rotation Arbitrary rotation to optionally apply on the
///                 resulting circle of photons.
/// \return         Circle of photons.
template <bool Center, typename Type, unsigned int Dimension, class Container, class>
Container Integrator::launch(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation)
{
    // Initialization
    static const unsigned int zero = 0;
    static const Type one = 1;
    static const Type four = 4;
    static const Type pi = four*std::atan(one);
    const Type step = (pi+pi)/Type(count+(!count));
    const Type r = std::sqrt((photon.dxdl()*photon.dxdl())+(photon.dydl()*photon.dydl())+(photon.dzdl()*photon.dzdl()));
    const Type rcos = r*std::cos(angle);
    const Type rsin = r*std::sin(angle);
    const Type costheta = std::cos(std::acos(photon.dzdl()/r));
    const Type sintheta = std::sin(std::acos(photon.dzdl()/r));
    const Type cosphi = std::cos(std::atan2(photon.dydl(), photon.dxdl()));
    const Type sinphi = std::sin(std::atan2(photon.dydl(), photon.dxdl()));
    const Type cospsi = std::cos(rotation);
    const Type sinpsi = std::sin(rotation);
    Type x = zero;
    Type y = zero;
    Type z = zero;
    Container result(count+Center, photon);
    
    // Loop over photons
    if (r > zero) {
        for (unsigned int istep = zero; istep < count; ++istep) {
            x = rsin*std::cos(istep*step);
            y = rsin*std::sin(istep*step);
            z = rcos;
            result[istep+Center].dxdl() = -cosphi*sinpsi*costheta*x-cosphi*cospsi*costheta*y+cosphi*sintheta*z+sinphi*sinpsi*y-sinphi*cospsi*x;
            result[istep+Center].dydl() = -sinphi*sinpsi*costheta*x-sinphi*cospsi*costheta*y-cosphi*sinpsi*y+cosphi*cospsi*x+sintheta*sinphi*z;
            result[istep+Center].dzdl() = sintheta*sinpsi*x+sintheta*cospsi*y+costheta*z;
        }
    }
    
    // Finalization
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- COMPUTATION ------------------------------ //
// Derivative of a photon
/// \brief          Derivative of a photon.
/// \details        Computes the derivative of the core components of a photon.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2 for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         Array Core array type.
/// \tparam	    Element data element for neighbouring cells
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Data Data type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \param[in,out]  output Output data.
/// \param[in]      input Input data.
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length.
/// \param[in]      elemsTsc Indexes of neighbouring cells
/// \return         Reference to the output data.
template <int Order, class Array, class Element, class Cosmology, class Octree, class Type, unsigned int Dimension, class Data, class Position, class Extent> 
Array& Integrator::dphotondl(Array& output, const Array& input, const Cosmology& cosmology, const Octree& octree, const Type length, std::vector<Element>& elemsTsc)
{
    // Initialization
    static const unsigned int a = 0;
    static const unsigned int t = 1;
    static const unsigned int x = 2;
    static const unsigned int y = 3;
    static const unsigned int z = 4;
    static const unsigned int dtdl = 5;
    static const unsigned int dxdl = 6;
    static const unsigned int dydl = 7;
    static const unsigned int dzdl = 8;
    static const Type two = 2;
    static const Type c2 = magrathea::Constants<Type>::c2();
    // Interpolate at photon position
    Data data = (Order == 0) ? (octree.ngp(input[x], input[y], input[z])) : ((Order == 1) ? (octree.cic(input[x], input[y], input[z])) : (Order == 2) ? (octree.tsc(elemsTsc, input[x], input[y], input[z])) : (Data()));
    // Estimate the total derivative of the potential over the affine parameter
    const Type dphidl = data.dphidt()*input[dtdl] + (data.dphidx()*input[dxdl] + data.dphidy()*input[dydl] + data.dphidz()*input[dzdl]);
    const Type dadt = Utility::interpolate(input[t], std::get<0>(cosmology), std::get<2>(cosmology));
    const Type scale = length;
    // Computation
    output[a] = input[dtdl]*dadt;
    output[t] = input[dtdl];
    output[x] = input[dxdl]/scale;
    output[y] = input[dydl]/scale;
    output[z] = input[dzdl]/scale;
    // Geodesic equations
    //output[dtdl] = -(two*dadt/input[a]*input[dtdl]*input[dtdl])-(two/c2*input[dtdl])*(data.dphidx()*input[dxdl]+data.dphidy()*input[dydl]+data.dphidz()*input[dzdl]);
    output[dtdl] = -(two*dadt/input[a]*input[dtdl]*input[dtdl])-(two/c2*input[dtdl])*(dphidl - data.dphidt()*input[dtdl]);
    output[dxdl] = -(two*dadt/input[a]*input[dtdl]*input[dxdl])+(two/c2*dphidl*input[dxdl])-(two*data.dphidx()*input[dtdl]*input[dtdl]);
    output[dydl] = -(two*dadt/input[a]*input[dtdl]*input[dydl])+(two/c2*dphidl*input[dydl])-(two*data.dphidy()*input[dtdl]*input[dtdl]);
    output[dzdl] = -(two*dadt/input[a]*input[dtdl]*input[dzdl])+(two/c2*dphidl*input[dzdl])-(two*data.dphidz()*input[dtdl]*input[dtdl]);
    
    // Finalization
    return output;
}
// -------------------------------------------------------------------------- //



// -------------------------------- EVOLUTION ------------------------------- //
// Geodesics integration
/// \brief          Geodesics integration.
/// \details        Integrates the geodesics equation of a photon.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2 for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Trajectory Trajectory type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Element Photon type.
/// \tparam         Data Data type.
/// \tparam         Core Core data type.
/// \tparam         Size Core size.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         Point Point type.
/// \param[in,out]  trajectory Trajectory.
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \return         Reference to the trajectory data.
template <int Order, bool RK4, bool Verbose, class Cosmology, class Octree, class Type, class Trajectory, unsigned int Dimension, class Element, class Data, class Core, unsigned int Size, class Position, class Extent, class Point> 
Trajectory& Integrator::integrate(Trajectory& trajectory, const Cosmology& cosmology, const Octree& octree, const Point& vobs, const Type length, const unsigned int nsteps)
{
    // Initialization
    static const Type zero = 0;
    static const Type one = 1;
    static const Type two = 2;
    static const Type six = 6;
    static const Type c = magrathea::Constants<Type>::c();
    static const Type c2 = magrathea::Constants<Type>::c2();
    static const Type position = Type(Position::num)/Type(Position::den);
    static const Type extent = Type(Extent::num)/Type(Extent::den);
    static const Type min = position-(extent/two);
    static const Type max = position+(extent/two);
    static const Data empty = Data();
    static const Data homogeneous = empty.copy().a(one);
    const Type scale = length;
    std::vector<std::pair<magrathea::SimpleHyperOctreeIndex<__uint128_t, 3>, Gravity<float, 3> >> elemsTsc(27);
    Type norm = Type();
    Data data = Data();
    Element photon = Element();
    Type ratio = Type();
    Type dl = Type();
    std::array<Core, 4> dcoredl = std::array<Core, 4>();

    // Integrate 
    if (!trajectory.empty()) {
        // Get initial data at the observer
        data = (Order == 0) ? (octree.ngp(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())) : ((Order == 1) ? (octree.cic(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())) : (Order == 2) ? (octree.tsc(elemsTsc, trajectory.back().x(), trajectory.back().y(), trajectory.back().z())):  (homogeneous));
	// Normalise with k^µ k_µ = 0
        norm = std::sqrt((c2*(one+two*(data.phi()/c2))*trajectory.back().dtdl()*trajectory.back().dtdl())/((one-two*(data.phi()/c2))*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl())));
        trajectory.back().dxdl() *= norm;
        trajectory.back().dydl() *= norm;
        trajectory.back().dzdl() *= norm;
	// Initialise photon
        trajectory.back().a() = Utility::interpolate(trajectory.back().t(), std::get<0>(cosmology), std::get<1>(cosmology));
        trajectory.back().level() = std::get<0>(*octree.locate(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())).level();
        trajectory.back().ah() = data.a();
        trajectory.back().ah() = zero;
        trajectory.back().rho() = data.rho();
        trajectory.back().phi() = data.phi();
        trajectory.back().dphidx() = data.dphidx();
        trajectory.back().dphidy() = data.dphidy();
        trajectory.back().dphidz() = data.dphidz();
        trajectory.back().dphidt() = data.dphidt();
        trajectory.back().dphidl() = zero;
        trajectory.back().laplacian() = zero;
        trajectory.back().redshift() = zero;
        trajectory.back().dsdl2() = (trajectory.back().a()*trajectory.back().a())*(-(c2*(one+two*(trajectory.back().phi()/c2))*trajectory.back().dtdl()*trajectory.back().dtdl())+((one-two*(trajectory.back().phi()/c2))*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl())));
        trajectory.back().error() = one-((one-two*(trajectory.back().phi()/c2))*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl()))/(c2*(one+two*(trajectory.back().phi()/c2))*trajectory.back().dtdl()*trajectory.back().dtdl());
        trajectory.back().distance() = zero;
        trajectory.back().isw() = zero;
        trajectory.back().chi() = zero;
        trajectory.back().iswold() = zero;
	trajectory.back().lambda() = zero;
	trajectory.back().s() = zero;
	ratio = trajectory.back().a()*trajectory.back().a()*(scale/c)/nsteps;
        dl = std::get<0>(*octree.locate(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())).template extent<Type, Position, Extent>()*ratio;
	// a0 given by the scale factor at the observer
	const Type a0 = trajectory.back().a();
#ifdef VELOCITYFIELD
        const Type v0n = (vobs[0]*trajectory[0].dxdl() + vobs[1]*trajectory[0].dydl() + vobs[2]*trajectory[0].dzdl())/(c*trajectory[0].dtdl());
#endif
	const Type phi0c2 = trajectory.back().phi()/c2;
        // Integrate
        while (data != empty) {
            // Photon index
            photon.index() = trajectory.back().index()+one;
            // Photon core
	    // RK4 integration
            if (RK4) {
                dphotondl<Order>(dcoredl[0], trajectory.back().core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl/two*dcoredl[0][i];
                }
                dphotondl<Order>(dcoredl[1], photon.core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl/two*dcoredl[1][i];
                }
                dphotondl<Order>(dcoredl[2], photon.core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < std::tuple_size<Core>::value; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl*dcoredl[2][i];
                }
                dphotondl<Order>(dcoredl[3], photon.core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+(dl/six)*(dcoredl[0][i]+two*dcoredl[1][i]+two*dcoredl[2][i]+dcoredl[3][i]);
                }
	    // Euler integration
            } else {
                dphotondl<Order>(photon.core(), trajectory.back().core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl*photon.core(i);
                }
            }
            // Photon extra
	    // Get data at new photon position
            data = (Order == 0) ? (octree.ngp(photon.x(), photon.y(), photon.z())) : ((Order == 1) ? (octree.cic(photon.x(), photon.y(), photon.z())) : (Order == 2) ? (octree.tsc(elemsTsc, photon.x(), photon.y(), photon.z())) : (homogeneous));
	    // If photon outside of the Extent, put empty data
            data = (!(photon.a() < zero) && ((photon.x() > min) && (photon.x() < max) && (photon.y() > min) && (photon.y() < max) && (photon.z() > min) && (photon.z() < max))) ? (data) : (empty);
            photon.level() = std::get<0>(*octree.locate(photon.x(), photon.y(), photon.z())).level();
            photon.ah() = data.a();
            photon.rho() = data.rho();
            photon.phi() = data.phi();
            photon.dphidx() = data.dphidx();
            photon.dphidy() = data.dphidy();
            photon.dphidz() = data.dphidz();
            photon.dphidt() = -data.dphidt(); // forward dphidt
	    photon.dphidl() = data.dphidt()*photon.dtdl() + (photon.dphidx()*photon.dxdl() + photon.dphidy()*photon.dydl() + photon.dphidz()*photon.dzdl());
            photon.laplacian() = (data.phi()-trajectory.back().phi())/dl;
            photon.dsdl2() = (photon.a()*photon.a())*((-c2*(one+two*(photon.phi()/c2))*photon.dtdl()*photon.dtdl())+((one-two*(photon.phi()/c2))*(photon.dxdl()*photon.dxdl()+photon.dydl()*photon.dydl()+photon.dzdl()*photon.dzdl())));
            photon.error() = one-((one-two*(photon.phi()/c2))*(photon.dxdl()*photon.dxdl()+photon.dydl()*photon.dydl()+photon.dzdl()*photon.dzdl()))/(c2*(one+two*(photon.phi()/c2))*photon.dtdl()*photon.dtdl());
            photon.distance() = zero;
            photon.isw() += (trajectory.back().dphidt()+photon.dphidt())*(photon.t()-trajectory.back().t())/c2;
            photon.iswold() -= 2*(photon.phi()-trajectory.back().phi())/c2 - ((trajectory.back().dphidx()+photon.dphidx())*(photon.x()-trajectory.back().x())*scale + (trajectory.back().dphidy()+photon.dphidy())*(photon.y()-trajectory.back().y())*scale + (trajectory.back().dphidz()+photon.dphidz())*(photon.z()-trajectory.back().z())*scale)/c2;
            photon.chi() = std::sqrt(photon.x()*photon.x()+photon.y()*photon.y()+photon.z()*photon.z());
	    photon.lambda() += dl; 
	    photon.s() += std::sqrt((photon.x()-trajectory.back().x())*(photon.x()-trajectory.back().x()) + (photon.y()-trajectory.back().y())*(photon.y()-trajectory.back().y()) + (photon.z()-trajectory.back().z())*(photon.z()-trajectory.back().z()));
#ifdef VELOCITYFIELD
            photon.redshift() = a0/photon.a()*(one + phi0c2 - photon.phi()/c2 + photon.isw() + (data.vx()*photon.dxdl() + data.vy()*photon.dydl() + data.vz()*photon.dzdl())/(c*photon.dtdl()) - v0n)-one;
            //photon.redshift() = a0/photon.a()*(one + (data.vx()*photon.dxdl() + data.vy()*photon.dydl() + data.vz()*photon.dzdl())/(c*photon.dtdl()) - v0n)-one;
#else
            photon.redshift() = a0/photon.a()*(one + phi0c2 - photon.phi()/c2)-one;
#endif
            // Next step
	    // Continue integrating if the photon in still in the Octree
            if (data != empty) {
		ratio = photon.a()*photon.a()*(scale/c)/nsteps;
                dl = std::get<0>(*(octree.locate(photon.x(), photon.y(), photon.z()))).template extent<Type, Position, Extent>()*ratio;
                trajectory.append(photon);
            }
        }        
        // Erase last element if non compatible
        if (!trajectory.empty()) {
            if (std::signbit(trajectory.back().redshift()) || (std::signbit(trajectory.back().a()))) {
                trajectory.pop();
            }
        }
    }

    // Finalization
    return trajectory;
}
// Geodesics integration
/// \brief          Geodesics integration.
/// \details        Integrates the geodesics equation of a photon until some condition.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2 for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Trajectory Trajectory type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Element Photon type.
/// \tparam         Data Data type.
/// \tparam         Core Core data type.
/// \tparam         Size Core size.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         Point Point type.
/// \param[in,out]  trajectory Trajectory.
/// \param[in]      interpolation Which type of stop criterion for the integration
/// \param[in]      interpRef Where to stop the integration
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \param[in]      kiTarget Normal to the plane needed for a given photon.
/// \return         Reference to the trajectory data.
template <int Order, bool RK4, bool Verbose, class Cosmology, class Octree, class Type, class Trajectory, unsigned int Dimension, class Element, class Data, class Core, unsigned int Size, class Position, class Extent, class Point> 
Trajectory& Integrator::integrate(Trajectory& trajectory, const std::string interpolation, const Type interpRef, const Cosmology& cosmology, const Octree& octree, const Point& vobs, const Type length, const unsigned int nsteps, const Point& kiTarget)
{

    // Initialization
    static const Type zero = 0;
    static const Type one = 1;
    static const Type two = 2;
    static const Type six = 6;
    static const Type c = magrathea::Constants<Type>::c();
    static const Type c2 = magrathea::Constants<Type>::c2();
    static const Type position = Type(Position::num)/Type(Position::den);
    static const Type extent = Type(Extent::num)/Type(Extent::den);
    static const Type min = position-(extent/two);
    static const Type max = position+(extent/two);
    static const Data empty = Data();
    static const Data homogeneous = empty.copy().a(one);
    const Type scale = length;
    std::vector<std::pair<magrathea::SimpleHyperOctreeIndex<__uint128_t, 3>, Gravity<float, 3> >> elemsTsc(27);
    Type norm = Type();
    Data data = Data();
    Element photon = Element();
    Type ratio = Type();
    Type dl = Type();
    std::array<Core, 4> dcoredl = std::array<Core, 4>();
    // Integrate 
    if (!trajectory.empty()) {
        // Get initial data at the observer
        data = (Order == 0) ? (octree.ngp(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())) : ((Order == 1) ? (octree.cic(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())) : (Order == 2) ? (octree.tsc(elemsTsc, trajectory.back().x(), trajectory.back().y(), trajectory.back().z())):  (homogeneous));
	// Normalise with k^µ k_µ = 0
        norm = std::sqrt((c2*(one+two*(data.phi()/c2))*trajectory.back().dtdl()*trajectory.back().dtdl())/((one-two*(data.phi()/c2))*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl())));
        trajectory.back().dxdl() *= norm;
        trajectory.back().dydl() *= norm;
        trajectory.back().dzdl() *= norm;
        trajectory.back().a() = Utility::interpolate(trajectory.back().t(), std::get<0>(cosmology), std::get<1>(cosmology));
        trajectory.back().level() = std::get<0>(*octree.locate(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())).level();
        trajectory.back().ah() = data.a();
        trajectory.back().ah() = zero;
        trajectory.back().rho() = data.rho();
        trajectory.back().phi() = data.phi();
        trajectory.back().dphidx() = data.dphidx();
        trajectory.back().dphidy() = data.dphidy();
        trajectory.back().dphidz() = data.dphidz();
        trajectory.back().dphidt() = data.dphidt();
        trajectory.back().dphidl() = zero;
        trajectory.back().laplacian() = zero;
        trajectory.back().redshift() = zero;
        trajectory.back().dsdl2() = (trajectory.back().a()*trajectory.back().a())*(-(c2*(one+two*(trajectory.back().phi()/c2))*trajectory.back().dtdl()*trajectory.back().dtdl())+((one-two*(trajectory.back().phi()/c2))*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl())));
        trajectory.back().error() = one-((one-two*(trajectory.back().phi()/c2))*(trajectory.back().dxdl()*trajectory.back().dxdl()+trajectory.back().dydl()*trajectory.back().dydl()+trajectory.back().dzdl()*trajectory.back().dzdl()))/(c2*(one+two*(trajectory.back().phi()/c2))*trajectory.back().dtdl()*trajectory.back().dtdl());
        trajectory.back().distance() = zero;
        trajectory.back().isw() = zero;
        trajectory.back().chi() = zero;
        trajectory.back().iswold() = zero;
	trajectory.back().lambda() = zero;
	trajectory.back().s() = zero;
	ratio = trajectory.back().a()*trajectory.back().a()*(scale/c)/nsteps;
        dl = std::get<0>(*octree.locate(trajectory.back().x(), trajectory.back().y(), trajectory.back().z())).template extent<Type, Position, Extent>()*ratio;
	// a0 given by the scale factor at the observer
	const Type a0 = trajectory.back().a();
#ifdef VELOCITYFIELD
        const Type v0n = (vobs[0]*trajectory[0].dxdl() + vobs[1]*trajectory[0].dydl() + vobs[2]*trajectory[0].dzdl())/(c*trajectory[0].dtdl());
#endif
	const Type phi0c2 = trajectory.back().phi()/c2;

        // Advance
        while (data != empty) {
            // Photon index
            photon.index() = trajectory.back().index()+one;
            // Photon core
	    // RK4 integration
            if (RK4) {
                dphotondl<Order>(dcoredl[0], trajectory.back().core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl/two*dcoredl[0][i];
                }
                dphotondl<Order>(dcoredl[1], photon.core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl/two*dcoredl[1][i];
                }
                dphotondl<Order>(dcoredl[2], photon.core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < std::tuple_size<Core>::value; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl*dcoredl[2][i];
                }
                dphotondl<Order>(dcoredl[3], photon.core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+(dl/six)*(dcoredl[0][i]+two*dcoredl[1][i]+two*dcoredl[2][i]+dcoredl[3][i]);
                }
	    // Euler integration
            } else {
                dphotondl<Order>(photon.core(), trajectory.back().core(), cosmology, octree, length, elemsTsc);
                for (unsigned int i = 0; i < Size; ++i) {
                    photon.core(i) = trajectory.back().core(i)+dl*photon.core(i);
                }
            }
            // Photon extra
	    // Get data at new photon position
            data = (Order == 0) ? (octree.ngp(photon.x(), photon.y(), photon.z())) : ((Order == 1) ? (octree.cic(photon.x(), photon.y(), photon.z())) : (Order == 2) ? (octree.tsc(elemsTsc, photon.x(), photon.y(), photon.z())) : (homogeneous));
	    // If the photon is outside of the Extent, put empty data
            data = (!(photon.a() < zero) && ((photon.x() > min) && (photon.x() < max) && (photon.y() > min) && (photon.y() < max) && (photon.z() > min) && (photon.z() < max))) ? (data) : (empty);
            photon.level() = std::get<0>(*octree.locate(photon.x(), photon.y(), photon.z())).level();
            photon.ah() = data.a();
            photon.rho() = data.rho();
            photon.phi() = data.phi();
            photon.dphidx() = data.dphidx();
            photon.dphidy() = data.dphidy();
            photon.dphidz() = data.dphidz();
            photon.dphidt() = -data.dphidt(); // forward dphidt
	    photon.dphidl() = data.dphidt()*photon.dtdl() + (photon.dphidx()*photon.dxdl() + photon.dphidy()*photon.dydl() + photon.dphidz()*photon.dzdl());
            photon.laplacian() = (data.phi()-trajectory.back().phi())/dl;
            photon.dsdl2() = (photon.a()*photon.a())*((-c2*(one+two*(photon.phi()/c2))*photon.dtdl()*photon.dtdl())+((one-two*(photon.phi()/c2))*(photon.dxdl()*photon.dxdl()+photon.dydl()*photon.dydl()+photon.dzdl()*photon.dzdl())));
            photon.error() = one-((one-two*(photon.phi()/c2))*(photon.dxdl()*photon.dxdl()+photon.dydl()*photon.dydl()+photon.dzdl()*photon.dzdl()))/(c2*(one+two*(photon.phi()/c2))*photon.dtdl()*photon.dtdl());
            photon.distance() = zero;
            photon.isw() += (trajectory.back().dphidt()+photon.dphidt())*(photon.t()-trajectory.back().t())/c2;
            photon.iswold() -= 2*(photon.phi()-trajectory.back().phi())/c2 - ((trajectory.back().dphidx()+photon.dphidx())*(photon.x()-trajectory.back().x())*scale + (trajectory.back().dphidy()+photon.dphidy())*(photon.y()-trajectory.back().y())*scale + (trajectory.back().dphidz()+photon.dphidz())*(photon.z()-trajectory.back().z())*scale)/c2;
            photon.chi() = std::sqrt(photon.x()*photon.x()+photon.y()*photon.y()+photon.z()*photon.z());
	    photon.lambda() += dl; 
	    photon.s() += std::sqrt((photon.x()-trajectory.back().x())*(photon.x()-trajectory.back().x()) + (photon.y()-trajectory.back().y())*(photon.y()-trajectory.back().y()) + (photon.z()-trajectory.back().z())*(photon.z()-trajectory.back().z()));
#ifdef VELOCITYFIELD
            photon.redshift() = a0/photon.a()*(one + phi0c2 - photon.phi()/c2 + photon.isw() + (data.vx()*photon.dxdl() + data.vy()*photon.dydl() + data.vz()*photon.dzdl())/(c*photon.dtdl()) - v0n)-one;
            //photon.redshift() = a0/photon.a()*(one + (data.vx()*photon.dxdl() + data.vy()*photon.dydl() + data.vz()*photon.dzdl())/(c*photon.dtdl()) - v0n)-one;
#else
            photon.redshift() = a0/photon.a()*(one + phi0c2 - photon.phi()/c2)-one;
#endif
            // Next step
            if (data != empty) {
		ratio = photon.a()*photon.a()*(scale/c)/nsteps;
                dl = std::get<0>(*(octree.locate(photon.x(), photon.y(), photon.z()))).template extent<Type, Position, Extent>()*ratio;
                trajectory.append(photon);
            }
	    // Stop the integration if the photon reaches some value for a given parameter
	    if((interpolation == "r") || (interpolation == "radius")){
		if(trajectory.back().chi() > interpRef){
		    break;
		}
	    } else if (interpolation == "lambda"){
		if(trajectory.back().lambda() > interpRef){
		    break;
		}
	    } else if (interpolation == "redshift"){
		if(trajectory.back().redshift() > interpRef){
		    break;
		}
	    } else if ((interpolation == "t") || (interpolation == "eta")){
		if(trajectory.back().t() > interpRef){
		    break;
		}
	    } else if (interpolation == "a"){
		if(trajectory.back().a() < interpRef){
		    break;
		}
	    } else if (interpolation == "plane"){
		if( (photon.x()*kiTarget[0] + photon.y()*kiTarget[1] + photon.z()*kiTarget[2]) > interpRef){
		    break;
		}
	    } else if (interpolation == "s"){
		if(trajectory.back().s() > interpRef){
		    break;
		}
	    } else{
		std::cout<<"# WARNING : Wrong stop criterion for integration"<<std::endl;
		std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		std::terminate();
	    }
        }
        
        // Erase last element if non compatible
        if (!trajectory.empty()) {
            if (std::signbit(trajectory.back().redshift()) || (std::signbit(trajectory.back().a()))) {
                trajectory.pop();
            }
        }
    }

    // Finalization
    return trajectory;
}

// Propagation of a ray bundle
/// \brief          Propagation of a ray bundle.
/// \details        Propagates a ray bundle calling the integrator for each
///                 photon.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2 for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Homogeneous Homogeneous reference.
/// \tparam         Point point type.
/// \param[in]      photon Central photon initial data.
/// \param[in]      count Number of other photons to use.
/// \param[in]      angle Half-angle at the cone vertex.
/// \param[in]      rotation Arbitrary rotation to optionally apply on the
///                 resulting circle of photons.
/// \param[in]      interpolation Stop condition : redshift, a, t, r.
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      length Spatial length in SI units.
/// \param[in]      nsteps Number of lambda steps per grid.
/// \param[in]      amin If different from zero, all photons should end by this
///                 value of a.
/// \param[in]      filenames File names of the output. If empty, no output. 
///                 If at least on percent sign, all trajectories are saved.
///                 Otherwise, only the central one is saved.
/// \param[in]      Homogeneous Optional homogeneous trajectory. If not provided
///                 the angular diameter distance use the inhomogeneous value of
///                 a. If provided, the homogeneous value of a for the given 
///                 radius is used.
/// \return         Central photon trajectory.
template <int Order, bool RK4, bool Verbose, class Cosmology, class Octree, class Type, unsigned int Dimension, class Homogeneous, class Point, class>
magrathea::Evolution<Photon<Type, Dimension> > Integrator::propagate(const Photon<Type, Dimension>& photon, const unsigned int count, const Type angle, const Type rotation, const std::string& interpolation, const Cosmology& cosmology, const Octree& octree, const Point& vobs, const Type length, const unsigned int nsteps, const Type amin, const std::string& filenames, const Homogeneous& homogeneous)
{
    // Initialization
    static const Type zero = 0;
    static const Type one = 1;
    static const Type two = 2;
    static const unsigned int first = 0;
    static const unsigned int center = 0;
    static const unsigned int x = 0;
    static const unsigned int y = 1;
    static const unsigned int z = 2;
    static const char percent = '%';
    static const unsigned int digits = std::numeric_limits<Type>::max_digits10;
    static const Type quarter = (Type(octree.extent().num)/Type(two*octree.extent().den))/two;
    static const Type limit = one/(two*two*two);
    std::vector<Photon<Type, Dimension> > initial = launch<true>(photon, count, angle, rotation); 
    std::vector<magrathea::Evolution<Photon<Type, Dimension> > > trajectories(initial.size());
    unsigned int ntrajectories = trajectories.size();
    unsigned int size = zero;
    std::vector<Type> last(ntrajectories);
    std::vector<std::vector<Type> > ref(ntrajectories);
    std::vector<std::array<std::vector<Type>, Dimension> > xyz(ntrajectories);
    std::array<Type, Dimension> coord = std::array<Type, Dimension>();
    std::pair<std::vector<Type>, std::vector<Type> > flrw;
    std::ofstream stream;

    // Integration
    // Loop over all the photons of a given bundle
    for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
	// Launch photon
        trajectories[itrajectory].append(initial[itrajectory]);
	// Integrate 
        integrate<Order, RK4, Verbose>(trajectories[itrajectory], cosmology, octree, vobs, length, nsteps);
        size = trajectories[itrajectory].size();
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            xyz[itrajectory][idim].resize(size);
        }
        for (unsigned int istep = 0; istep < size; ++istep) {
            xyz[itrajectory][x][istep] = trajectories[itrajectory][istep].x();
            xyz[itrajectory][y][istep] = trajectories[itrajectory][istep].y();
            xyz[itrajectory][z][istep] = trajectories[itrajectory][istep].z();
        }
	// Check if the integration is well done
        if (!size) {
            ntrajectories = size;
        } else {
            coord[x] = trajectories[itrajectory].back().x()-trajectories[itrajectory].front().x();
            coord[y] = trajectories[itrajectory].back().y()-trajectories[itrajectory].front().y();
            coord[z] = trajectories[itrajectory].back().z()-trajectories[itrajectory].front().z();
            last[itrajectory] = std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
	    // Check if photons could go beyond 1/4th of the Extent
            ntrajectories *= (last[itrajectory] > quarter);
	    // Check if values for the scale factor from data are normal and within the expected range
            ntrajectories *= (!(std::isnormal(amin) && std::isnormal(trajectories[itrajectory].back().ah()) && (trajectories[itrajectory].back().ah() < one))) || (!(trajectories[itrajectory].back().ah() > amin));
        }
    }
    ntrajectories *= (std::abs((*std::max_element(last.begin(), last.end()))-(*std::min_element(last.begin(), last.end())))/(*std::max_element(last.begin(), last.end())) < limit);

    // Interpolation for the photon from bundle at some given parameter
    if (ntrajectories > 1) {
        if (interpolation == "redshift") {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].redshift();
                }
            }
        } else if (interpolation == "a") {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = 1/trajectories[itrajectory][istep].a() - 1;
                }
            }
        } else if ((interpolation == "t") || (interpolation == "eta")) {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].t();
                }
            }
	} else if (interpolation == "lambda") {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].lambda();
                }
            }
        } else if ((interpolation == "r") || (interpolation == "radius")) {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    coord[x] = trajectories[itrajectory][istep].x()-trajectories[itrajectory][first].x();
                    coord[y] = trajectories[itrajectory][istep].y()-trajectories[itrajectory][first].y();
                    coord[z] = trajectories[itrajectory][istep].z()-trajectories[itrajectory][first].z();
                    ref[itrajectory][istep] = std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]);
                }
            }
        } else {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                size = trajectories[itrajectory].size();
                ref[itrajectory].resize(size);
                for (unsigned int istep = 0; istep < size; ++istep) {
                    ref[itrajectory][istep] = trajectories[itrajectory][istep].redshift();
                }
            }
        }
    }
    // Computation of the angular diameter distance
    if (ntrajectories > 1) {
        size = trajectories[center].size();
        for (unsigned int istep = 0; istep < size; ++istep) {
            for (unsigned int itrajectory = 1; itrajectory < ntrajectories; ++itrajectory) {
		// interpolate x position at some reference parameter given by central ray
                coord[x] = Utility::interpolate(ref[center][istep], ref[itrajectory], xyz[itrajectory][x])-xyz[center][x][istep];
		// interpolate y position at some reference parameter given by central ray 
                coord[y] = Utility::interpolate(ref[center][istep], ref[itrajectory], xyz[itrajectory][y])-xyz[center][y][istep];
		// interpolate z position at some reference parameter given by central ray
                coord[z] = Utility::interpolate(ref[center][istep], ref[itrajectory], xyz[itrajectory][z])-xyz[center][z][istep]; 
		// Compute absolute distance from the central ray to each bundle ray
                trajectories[center][istep].distance() += std::sqrt(coord[x]*coord[x]+coord[y]*coord[y]+coord[z]*coord[z]); 
            }
	    // Mean distance from central ray to bundle rays
            trajectories[center][istep].distance() /= (ntrajectories-1); 
        }
	// Convert comoving to angular diameter distance
        for (unsigned int istep = 0; istep < size; ++istep) {
            trajectories[center][istep].distance() *= (length*trajectories[center][istep].a())/angle;
        }
    }

    // Output
    if ((ntrajectories > 0) && (!filenames.empty())) {
        if (filenames.find(percent) == std::string::npos) {
            stream.open(filenames);
            Output::save(stream, trajectories[center], digits); 
            stream.close();
        } else {
            for (unsigned int itrajectory = 0; itrajectory < ntrajectories; ++itrajectory) {
                stream.open(Output::name(std::make_pair(filenames, itrajectory)));
                Output::save(stream, trajectories[itrajectory], digits); 
                stream.close();
            }
        }
    }

    // Finalization
    return ntrajectories > 0 ? trajectories[0] : magrathea::Evolution<Photon<Type, Dimension> >();
}
// -------------------------------------------------------------------------- //



// -------------------------------------------------------------------------- //


// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Integrator.
/// \return         0 if no error.
/*int Integrator::example()
{
    // Initialize
    std::cout<<"BEGIN = Integrator::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::random_device device;
    std::mt19937 engine(device());
    std::uniform_real_distribution<double> distribution(0, 1);
    std::array<double, 3> beg = std::array<double, 3>({{0., 0., 0.}});
    std::array<double, 3> end = std::array<double, 3>({{16., 23., 42.}});
    std::array<double, 10> array = std::array<double, 10>();
    std::array<std::vector<double>, 4> cosmology = std::array<std::vector<double>, 4>();
    magrathea::SimpleHyperOctree<double, magrathea::SimpleHyperOctreeIndex<unsigned long long int, 3>, Gravity<float, 3> > octree(0, 2);
    magrathea::HyperSphere<3> sphere = magrathea::HyperSphere<3>::unit();
    magrathea::Evolution<Photon<double, 3> > trajectory;
    Photon<double, 3> photon;
    Cone<> cone(beg, end, 0.42);
    std::vector<Cone<> > cones(3, cone);
    double one = 1;
    
    // Construction
    Integrator integrator;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"Integrator() : "                            ; Integrator(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"integrator = Integrator() : "               ; integrator = Integrator(); std::cout<<std::endl;

    // Initialization
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Initialization : "                                                                  <<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(sphere, cones[0], engine, distribution) : "                       <<integrator.launch(sphere, cones[0], engine, distribution)<<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(sphere, cones[0], cones, engine, distribution) : "                <<integrator.launch(sphere, cones[0], cones, engine, distribution)<<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(beg[0], beg[1], beg[2], end[0], end[1], end[2]) : "               <<integrator.launch(beg[0], beg[1], beg[2], end[0], end[1], end[2])<<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.launch(photon, 3, 0.42, 0.1).size() : "                                  <<integrator.launch(photon, 3, 0.42, 0.1).size()<<std::endl;
    
    // Computation
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Computation : "                                                                     <<std::endl;
    std::cout<<std::setw(width*2)<<"integrator.dphotondl(array, array, cosmology, octree, one, one, one)[0] : "         <<integrator.dphotondl(array, array, cosmology, octree, one, one, one)[0]<<std::endl;

    // Evolution
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Evolution : "                                                                                                               <<std::endl;
    std::cout<<std::setw(width*3)<<"integrator.integrate(trajectory, cosmology, octree, one, one).size() : "                                                    <<integrator.integrate(trajectory, cosmology, octree, one, one).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"integrator.propagate(photon, 3, 0.42, 0.1, \"a\", cosmology, octree, one, one).size()"                                      <<integrator.propagate(photon, 3, 0.42, 0.1, "a", cosmology, octree, one, one).size()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Integrator::example()"<<std::endl;
    return 0;
}*/
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
#endif // INTEGRATOR_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
