/* ******************************** UTILITY ********************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Utility
// DESCRIPTION :    List of utilities for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           utility.h
/// \brief          List of utilities for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <chrono>
#include <array>
#include <tuple>
#include <utility>
#include <thread>
#include <vector>
#include <algorithm>
// Include libs
// Include project
#include "magrathea/hypercube.h"
#include "magrathea/hypersphere.h"
#include "cone.h"
// Misc
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// List of utilities for raytracing
/// \brief          List of utilities for raytracing.
/// \details        Provides a list of general purpose utilities for raytracing
///                 like parallelization, collision detection or interpolation.
class Utility final
{
    // Parallelization
    /// \name           Parallelization
    //@{
    public:
        template <int Default = 0, typename Type, class Function, class = typename std::enable_if<(std::is_convertible<decltype(std::declval<Type>()+std::declval<Type>()), int>::value) && (!std::is_function<typename std::result_of<Function(Type)>::type>::value)>::type> static double parallelize(const Type nsteps, Function&& function, const int nthreads = (Default != 0) ? (Default) : (std::thread::hardware_concurrency()));
        template <int Default = 0, typename Type, class Function, class = typename std::enable_if<(!std::is_void<decltype(std::declval<Type>()/std::declval<Type>())>::value) && (!std::is_function<typename std::result_of<Function(Type)>::type>::value)>::type> static double parallelize(const Type& first, const Type& last, const Type& increment, Function&& function, const int nthreads = (Default != 0) ? (Default) : (std::thread::hardware_concurrency()));
        template <int Default = 0, typename Iterator, class Function, class = typename std::enable_if<(!std::is_void<decltype(*std::declval<Iterator>())>::value) && (!std::is_function<typename std::result_of<Function(decltype(*std::declval<Iterator>()))>::type>::value)>::type> static double parallelize(const Iterator& first, const Iterator& last, Function&& function, const int nthreads = (Default != 0) ? (Default) : (std::thread::hardware_concurrency()));
    //@}
    
    // Geometry
    /// \name           Geometry
    //@{
    public:
        template <unsigned int Dimension, unsigned int Index = 0, class Vector, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type> static constexpr Scalar distance(const Vector& first, const Vector& second); 
        template <unsigned int Dimension, unsigned int Index = 0, class Vector, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type> static constexpr Scalar dot(const Vector& first, const Vector& second);
        template <unsigned int Dimension, unsigned int Index = 0, class Vector, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type> static inline Vector cross(const Vector& first, const Vector& second);
        template <unsigned int Dimension, unsigned int Index = 0, class Vector, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type> static inline Vector join(const Vector& first, const Vector& second);
        template <class Operator, unsigned int Dimension, unsigned int Index = 0, class Vector, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type> static inline Vector apply(const Vector& first, const Vector& second);
        template <unsigned int Dimension, class Vector, typename Scalar> static constexpr Scalar radius(const magrathea::HyperSphere<Dimension, Vector, Scalar>& hypersphere);
        template <unsigned int Dimension, class Vector, typename Scalar> static constexpr Scalar radius(const magrathea::HyperCube<Dimension, Vector, Scalar>& hypercube);
        template <unsigned int Dimension, class Vector, typename Scalar, typename Type = Scalar> static inline magrathea::HyperCube<Dimension, Vector, Scalar> cubify(const magrathea::HyperSphere<Dimension, Vector, Scalar>& hypersphere, Type factor = Type(1));
        template <unsigned int Dimension, class Vector, typename Scalar, typename Type = Scalar> static inline magrathea::HyperSphere<Dimension, Vector, Scalar> spherify(const magrathea::HyperCube<Dimension, Vector, Scalar>& hypercube, Type factor = Type(1));
        template <template <unsigned int, class, typename> class First, template <unsigned int, class, typename> class Second, unsigned int Dimension, class Vector, typename Scalar> static inline bool collide(const First<Dimension, Vector, Scalar>& first, const Second<Dimension, Vector, Scalar>& second);
        template <template <unsigned int, class, typename> class Object, unsigned int Dimension, class Vector, typename Scalar, class = typename std::enable_if<Dimension == 3>::type> static inline bool collide(const Object<Dimension, Vector, Scalar>& object, const Cone<Vector, Scalar>& cone, const double& inspheresize);
    //@}
    
    // Interpolation
    /// \name           Interpolation
    //@{
    public:
        template <typename Type, class Container, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static inline Type interpolate(const Type x0, const Container& x, const Container& y);
        template <typename Type, class Container, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static inline Type rinterpolate(const Type x0, const Container& x, const Container& y);
        template <typename Type, class Container, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static inline Type interpolate(const Type x0, const Container& x, const Container& y, const Container& dydx);
        template <class Container, class... Containers, class = typename std::enable_if<sizeof...(Containers) == 2 || sizeof...(Containers) == 3>::type> static inline Container reinterpolate(const Container& x0, Containers&&... containers);
        template <int Direction = 0, typename Type, class Container, typename T = Type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static inline Type differentiate(const Type x0, const Container& x, const Container& y, const unsigned int neighbourhood = 1);
        template <int Derivative = 0, typename Type, class Container, typename T = Type, class = typename std::enable_if<(Derivative >= 0) && (Derivative <= 3) && (std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value)>::type> static inline Type filter(const Type x0, const Container& x, const Container& y, const unsigned int neighbourhood = 1);
    //@}
    
    // Evolution
    /// \name           Evolution
    //@{
    public:
        template <class Container, typename Type, class = typename std::enable_if<std::is_same<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static Container reverse(const Container& container, const Type value);
        template <class Container, typename Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, class Function, class = typename std::enable_if<!std::is_function<typename std::result_of<Function(Type, Type)>::type>::value>::type> static Container smooth(const Container& x, const Container& y, Function&& kernel, const unsigned int window = 0);
        template <class Container, typename Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static Container integrate(const Container& x, const Container& y, const Type value = Type());
        template <int Direction = 0, class Container, typename Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type> static Container derive(const Container& x, const Container& y, const unsigned int neighbourhood = 1);
    //@}

	template < class Type > static inline std::array< std::array<Type, 2>, 2 > invMatrix2d(const std::array< std::array<Type, 2>, 2 >& A );
	template < class Type > static inline std::array< std::array<Type, 3>, 3 > invMatrix3d(const std::array< std::array<Type, 3>, 3 >& A );

    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ----------------------------- PARALLELIZATION ---------------------------- //
// Parallelize a loop
/// \brief          Parallelize a loop.
/// \details        Executes the provided function on each index of the loop
///                 using the specified number of threads. 
/// \tparam         Default Default concurrency where zero means hardware
///                 concurrency.
/// \tparam         Type Loop index type.
/// \tparam         Function Function type taking a loop index as argument.
/// \param[in]      nsteps Total number of steps.
/// \param[in]      function Function.
/// \param[in]      nthreads Number of threads.
/// \return         Elapsed time in seconds.
template <int Default, typename Type, class Function, class>
double Utility::parallelize(const Type nsteps, Function&& function, const int nthreads)
{
    static const Type zero = Type();
    const std::chrono::high_resolution_clock::time_point tbegin = std::chrono::high_resolution_clock::now();
    const Type ntasks = std::max(static_cast<int>(1), nthreads);
    const Type group = std::max(Type(nsteps > zero), nsteps/ntasks);
    std::vector<std::thread> threads;
    Type istep = zero;
    threads.reserve(ntasks);
    // Each thread takes care of some tasks
    for (istep = zero; istep < floor(nsteps/ntasks)*ntasks; istep += group) {
        threads.push_back(std::thread([=, &nsteps, &group, &function](){for (Type i = istep; i < std::min(istep+group, nsteps); ++i) function(i);}));
    }
    // Run the last few steps
    for (istep = istep; istep < nsteps; ++istep) {
        function(istep);
    }
    std::for_each(threads.begin(), threads.end(), [](std::thread& current){current.join();});

    return std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::high_resolution_clock::now()-tbegin).count();
}

// Parallelize an iteration over a range of values
/// \brief          Parallelize an iteration over a range of values.
/// \details        Executes the provided function on each values of the range
///                 using the provided increment and the specified number of 
///                 threads.
/// \tparam         Default Default concurrency where zero means hardware
///                 concurrency.
/// \tparam         Type Value type.
/// \tparam         Function Function type taking a value as argument.
/// \param[in]      first First value.
/// \param[in]      last Last value.
/// \param[in]      increment Increment.
/// \param[in]      function Function.
/// \param[in]      nthreads Number of threads.
/// \return         Elapsed time in seconds.
template <int Default, typename Type, class Function, class>
double Utility::parallelize(const Type& first, const Type& last, const Type& increment, Function&& function, const int nthreads)
{
    const std::chrono::high_resolution_clock::time_point tbegin = std::chrono::high_resolution_clock::now();
    const long long int ntasks = std::max(static_cast<int>(1), nthreads);
    const long long int nsteps = (last-first)/increment;
    const long long int size = nsteps+((first < last) ? ((first+Type(nsteps)*increment) < last) : (((first+Type(nsteps))*increment) > last));
    const long long int group = std::max(static_cast<long long int>(size > 0), size/ntasks);
    std::vector<std::thread> threads;
    long long int ivalue = 0;
    threads.reserve(ntasks);
    for (ivalue = 0; ivalue < size-group; ivalue += group) {
        threads.push_back(std::thread([=, &first, &increment, &size, &group, &function](){for (long long int i = ivalue; i < std::min(ivalue+group, size); ++i) function(first+Type(i)*increment);}));
    }
    for (; ivalue < size; ++ivalue) {
        function(first+Type(ivalue)*increment);
    }
    std::for_each(threads.begin(), threads.end(), [](std::thread& current){current.join();});
    return std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::high_resolution_clock::now()-tbegin).count();
}

// Parallelize an iteration over a range of iterators
/// \brief          Parallelize an iteration over a range of iterators.
/// \details        Executes the provided function on each element of the 
///                 iterator range using the specified number of threads. 
/// \tparam         Default Default concurrency where zero means hardware
///                 concurrency.
/// \tparam         Iterator Iterator type.
/// \tparam         Function Function type taking an element as argument.
/// \param[in]      first First iterator.
/// \param[in]      last Last iterator.
/// \param[in]      function Function.
/// \param[in]      nthreads Number of threads.
/// \return         Elapsed time in seconds.
template <int Default, typename Iterator, class Function, class>
double Utility::parallelize(const Iterator& first, const Iterator& last, Function&& function, const int nthreads)
{
    const std::chrono::high_resolution_clock::time_point tbegin = std::chrono::high_resolution_clock::now();
    const long long int ntasks = std::max(static_cast<int>(1), nthreads);
    const long long int group = std::max(static_cast<long long int>(first < last), static_cast<long long int>((last-first)/ntasks));
    std::vector<std::thread> threads;
    Iterator it = first;
    threads.reserve(ntasks);
    for (it = first; it < last-group; it += group) {
        threads.push_back(std::thread([=, &last, &group, &function](){std::for_each(it, std::min(it+group, last), function);}));
    }
    std::for_each(it, last, function);
    std::for_each(threads.begin(), threads.end(), [](std::thread& current){current.join();});
    return std::chrono::duration_cast<std::chrono::duration<double> >(std::chrono::high_resolution_clock::now()-tbegin).count();
}
// -------------------------------------------------------------------------- //



// -------------------------------- GEOMETRY -------------------------------- //
// Euclidian distance
/// \brief          Euclidian distance.
/// \details        Computes the euclidian distance between two points.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Index Inner computation index.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      first First point.
/// \param[in]      second Second point.
/// \return         Value of the distance.
template <unsigned int Dimension, unsigned int Index, class Vector, typename Scalar> 
constexpr Scalar Utility::distance(const Vector& first, const Vector& second)
{
    return (Index < Dimension) ? ((Index > 0) ? (((first[Index]-second[Index])*(first[Index]-second[Index]))+distance<Dimension, Index+(Index < Dimension)>(first, second)) : (std::sqrt(((first[Index]-second[Index])*(first[Index]-second[Index]))+distance<Dimension, Index+(Index < Dimension)>(first, second)))) : (Scalar());
}

// Dot product
/// \brief          Dot product.
/// \details        Computes the dot product of two vectors.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Index Inner computation index.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      first First point.
/// \param[in]      second Second point.
/// \return         Value of the dot product.
template <unsigned int Dimension, unsigned int Index, class Vector, typename Scalar> 
constexpr Scalar Utility::dot(const Vector& first, const Vector& second)
{
    return (Index < Dimension) ? ((first[Index]*second[Index])+dot<Dimension, Index+(Index < Dimension)>(first, second)) : (Scalar());
}

// Three dimensional cross product
/// \brief          Three dimensional cross product.
/// \details        Computes the cross product of two vectors in three
///                 dimensions.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Index Inner computation index.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      first First point.
/// \param[in]      second Second point.
/// \return         Result of the cross product.
template <unsigned int Dimension, unsigned int Index, class Vector, typename Scalar> 
inline Vector Utility::cross(const Vector& first, const Vector& second)
{
    Vector result = Vector();
    if ((Dimension == 3) && (!std::is_void<Scalar>::value)) {
        result[Index+0] = first[Index+1]*second[Index+2]-first[Index+2]*second[Index+1];
        result[Index+1] = first[Index+2]*second[Index+0]-first[Index+0]*second[Index+2];
        result[Index+2] = first[Index+0]*second[Index+1]-first[Index+1]*second[Index+0];
    }
    return result;
}

// Join two points
/// \brief          Join two points.
/// \details        Computes the vector resulting from going from the first
///                 point to the second point.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Index Inner computation index.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      first First point.
/// \param[in]      second Second point.
/// \return         Resulting joining vector.
template <unsigned int Dimension, unsigned int Index, class Vector, typename Scalar>
inline Vector Utility::join(const Vector& first, const Vector& second)
{
    Vector result = Vector();
    for (unsigned int idim = Index; idim < Index+Dimension; ++idim) {
        result[idim] = second[idim]-first[idim];
    }
    return result;
}

// Apply operator to all elements
/// \brief          Apply operator to all elements.
/// \details        Computes the vector resulting from the application of an
///                 operator to all elements of two vectors.
/// \tparam         Operator Operator type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Index Inner computation index.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      first First point.
/// \param[in]      second Second point.
/// \return         Resulting vector.
template <class Operator, unsigned int Dimension, unsigned int Index, class Vector, typename Scalar>
inline Vector Utility::apply(const Vector& first, const Vector& second)
{
    static const Operator op;
    Vector result = Vector();
    for (unsigned int idim = Index; idim < Index+Dimension; ++idim) {
        result[idim] = op(first[idim], second[idim]);
    }
    return result;
}

// Radius of an hypersphere
/// \brief          Radius of an hypersphere.
/// \details        Computes the radius of the provided hypersphere.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      hypersphere Hypersphere.
/// \return         Value of the radius.
template <unsigned int Dimension, class Vector, typename Scalar> 
constexpr Scalar Utility::radius(const magrathea::HyperSphere<Dimension, Vector, Scalar>& hypersphere)
{
    return hypersphere.radius();
}

// Radius of an hypercube
/// \brief          Radius of an hypercube.
/// \details        Computes the radius corresponding to half the diagonal of
///                 the provided hypercube.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      hypercube Hypercube.
/// \return         Value of the radius.
template <unsigned int Dimension, class Vector, typename Scalar> 
constexpr Scalar Utility::radius(const magrathea::HyperCube<Dimension, Vector, Scalar>& hypercube)
{
    return hypercube.diagonal()/Scalar(2);
}

// Hypersphere to hypercube conversion
/// \brief          Hypersphere to hypercube conversion.
/// \details        Converts a hypersphere to a hypercubecube with its diagonal 
///                 equals to the diameter with the provided scaling.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      hypersphere Input hypersphere.
/// \param[in]      scale Scaling factor.
/// \return         Hypercube resulting from the hypersphere conversion.
template <unsigned int Dimension, class Vector, typename Scalar, typename Type> 
inline magrathea::HyperCube<Dimension, Vector, Scalar> Utility::cubify(const magrathea::HyperSphere<Dimension, Vector, Scalar>& hypersphere, Type scale)
{
    return magrathea::HyperCube<Dimension, Vector, Scalar>(hypersphere.position(), (radius(hypersphere)/std::sqrt(Scalar(Dimension)))*Scalar(2)*scale);
}

// Hypercube to hypersphere conversion
/// \brief          Hypercube to hypersphere conversion.
/// \details        Converts a hypercube to a hypersphere with its diameter 
///                 equals to the diagonal with the provided scaling.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      hypercube Input hypercube.
/// \param[in]      scale Scaling factor.
/// \return         Hypersphere resulting from the hypercube conversion.
template <unsigned int Dimension, class Vector, typename Scalar, typename Type> 
inline magrathea::HyperSphere<Dimension, Vector, Scalar> Utility::spherify(const magrathea::HyperCube<Dimension, Vector, Scalar>& hypercube, Type scale)
{
    return magrathea::HyperSphere<Dimension, Vector, Scalar>(hypercube.position(), radius(hypercube)*scale);
}

// Collision between two hyperobjects
/// \brief          Collision between two hyperobjects.
/// \details        Detects collision between two geometrical objects in
///                 arbitrary dimension treated as hyperspheres.
/// \tparam         First First object type.
/// \tparam         Second Second object type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      first First object.
/// \param[in]      second Second object.
/// \return         True if collision, false otherwise.
template <template <unsigned int, class, typename> class First, template <unsigned int, class, typename> class Second, unsigned int Dimension, class Vector, typename Scalar>
inline bool Utility::collide(const First<Dimension, Vector, Scalar>& first, const Second<Dimension, Vector, Scalar>& second)
{
    return (distance<Dimension>(first.position(), second.position()) < (radius(first)+radius(second)));
}

// Collision between an hyperobject and a cone
/// \brief          Collision between an hyperobject and a cone.
/// \details        Detects collision between a geometrical object and a three
///                 dimensional cone.
/// \tparam         Object Object type.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
/// \param[in]      object Geometrical object.
/// \param[in]      cone Three dimensional cone.
/// \return         True if collision, false otherwise.
template <template <unsigned int, class, typename> class Object, unsigned int Dimension, class Vector, typename Scalar, class>
inline bool Utility::collide(const Object<Dimension, Vector, Scalar>& object, const Cone<Vector, Scalar>& cone, const double& inspheresize)
{
    const Scalar norm = cone.length();
    Scalar length = Scalar();
    Scalar distance = Scalar();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        length += (cone.base(idim)-cone.vertex(idim))*(object.position(idim)-cone.vertex(idim));
    }

    length /= norm;
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        distance += cone.template pow<2>(object.position(idim)-(cone.vertex(idim)+(cone.base(idim)-cone.vertex(idim))*(length/norm))); 
    }

    return ((std::sqrt(distance) < (5*radius(object)+length*std::tan(cone.angle()) + 0.4*inspheresize*std::exp(-0.7*(length - inspheresize*std::cos(cone.angle()))*(length - inspheresize*std::cos(cone.angle()))/(inspheresize*inspheresize)))) && !(length < Scalar()) && (length < (norm+radius(object)*Scalar(2))));
}


// -------------------------------------------------------------------------- //



// ------------------------------ INTERPOLATION ----------------------------- //
// Linear interpolation
/// \brief          Linear interpolation.
/// \details        Interpolates the value at the given position linearly.
/// \tparam         Type Data type.
/// \tparam         Container Container type.
/// \param[in]      x0 Interpolation abscissa.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \return         Interpolated ordinate.
template <typename Type, class Container, class> 
inline Type Utility::interpolate(const Type x0, const Container& x, const Container& y)
{
    static const Type zero = Type();
    const long long int n = std::min(x.size(), y.size());
    long long int i = std::distance(std::begin(x), std::upper_bound(std::begin(x), std::begin(x)+n, x0, [](const Type x1, const Type x2){return x1 < x2;}));
    i += (n > 1)*((i <= 0)-(i >= n));
    return (n > 1) ? (y[i-1]+((y[i]-y[i-1])*(x0-x[i-1])/(x[i]-x[i-1]))) : ((n > 0) ? (y[0]) : (zero));
}

// Linear interpolation
/// \brief          Linear interpolation.
/// \details        Interpolates the value at the given position linearly, with decreasing x.
/// \tparam         Type Data type.
/// \tparam         Container Container type.
/// \param[in]      x0 Interpolation abscissa.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \return         Interpolated ordinate.
template <typename Type, class Container, class> 
inline Type Utility::rinterpolate(const Type x0, const Container& x, const Container& y)
{
    static const Type zero = Type();
    const long long int n = std::min(x.size(), y.size());
    long long int i = std::distance(std::begin(x), std::upper_bound(std::begin(x), std::begin(x)+n, x0, [](const Type x1, const Type x2){return x1 > x2;}));
    i += (n > 1)*((i <= 0)-(i >= n));
    return (n > 1) ? (y[i-1]+((y[i]-y[i-1])*(x0-x[i-1])/(x[i]-x[i-1]))) : ((n > 0) ? (y[0]) : (zero));
}

// Cubic spline interpolation
/// \brief          Cubic spline interpolation.
/// \details        Interpolates the value at the given position using a cubic
///                 spline with the specified derivative.
/// \tparam         Type Data type.
/// \tparam         Container Container type.
/// \param[in]      x0 Interpolation abscissa.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \param[in]      dydx Derivatives.
/// \return         Interpolated ordinate.
template <typename Type, class Container, class> 
inline Type Utility::interpolate(const Type x0, const Container& x, const Container& y, const Container& dydx)
{
    static const Type one = Type(1);
    const long long int n = std::min(x.size(), y.size());
    long long int i = std::distance(std::begin(x), std::upper_bound(std::begin(x), std::begin(x)+n, x0, [](const Type x1, const Type x2){return x1 < x2;}));
    Type t = one;
    i += (n > 1)*((i <= 0)-(i >= n));
    t = (n > 1) ? ((x0-x[i-1])/(x[i]-x[i-1])) : (one-one);
    return (n > 1) ? (one-t)*y[i-1]+t*y[i]+t*(one-t)*((dydx[i-1]*(x[i]-x[i-1])-(y[i]-y[i-1]))*(one-t)+(-dydx[i]*(x[i]-x[i-1])+(y[i]-y[i-1]))*t) : ((n > 0) ? (y[0]) : (one-one));
}

// Container reinterpolation
/// \brief          Container reinterpolation.
/// \details        Interpolates each value of the container.
/// \tparam         Container Container type.
/// \tparam         Containers Containers types.
/// \param[in]      x0 Interpolation abscissae.
/// \param[in]      containers Input containers.
/// \return         Interpolated ordinates.
template <class Container, class... Containers, class> 
inline Container Utility::reinterpolate(const Container& x0, Containers&&... containers)
{
    const long long int n = std::distance(std::begin(x0), std::end(x0));
    Container result = x0;
    for (long long int i = 0; i < n; ++i) {
        result[i] = interpolate(x0[i], std::forward<Containers>(containers)...);
    }
    return result;
}

// Fornberg differentiation
/// \brief          Fornberg differentiation.
/// \details        Computes the value of the derivative at the given position
///                 using a fourth order Fornberg algorithm. The greater the 
///                 neighbourhood, the greater the computing distance is.
/// \tparam         Direction Sign for centered, backward or forward 
///                 differentiation.
/// \tparam         Type Data type.
/// \tparam         Container Container type.
/// \tparam         T Conversion type.
/// \param[in]      x0 Differentiation abscissa.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \param[in]      neighbourhood Computing distance.
/// \return         Derivative at the provided abscissa.
template <int Direction, typename Type, class Container, typename T, class>
inline Type Utility::differentiate(const Type x0, const Container& x, const Container& y, const unsigned int neighbourhood)
{
    static const unsigned int one = 1;
    static const unsigned int size = 9;
    static const std::array<Type, size> centered({{T(1)/T(280), T(-4)/T(105), T(1)/T(5), T(-4)/T(5), T(0), T(4)/T(5), T(-1)/T(5), T(4)/T(105), T(-1)/T(280)}});
    static const std::array<Type, size> forward({{T(-761)/T(280), T(8), T(-14), T(56)/T(3), T(-35)/T(2), T(56)/T(5), T(-14)/T(3), T(8)/T(7), T(-1)/T(8)}});
    static const std::array<Type, size> backward({{T(761)/T(280), T(-8), T(14), T(-56)/T(3), T(35)/T(2), T(-56)/T(5), T(14)/T(3), T(-8)/T(7), T(1)/T(8)}});
    const long long int n = std::min(x.size(), y.size());
    const long long int i = std::distance(std::begin(x), std::upper_bound(std::begin(x), std::begin(x)+n, x0, [](const Type x1, const Type x2){return x1 < x2;}));
    const long long int j = i+((n > 1)*((i <= 0)-(i >= n)));
    const Type h = (n > 1) ? ((x[j]-x[j-1])*std::max(one, neighbourhood)) : (Type());
    Type result = Type();
    if (n > 1) {
        if (Direction == 0) {
            for (unsigned int k = 0; k < size/2; ++k) {
                result += centered[k]*interpolate(x0-(size/2-k)*h, x, y)+centered[size-(k+1)]*interpolate(x0+(size/2-k)*h, x, y);
            }
            result += centered[size/2]*interpolate(x0, x, y);
        } else {
            for (unsigned int k = 0; k < size; ++k) {
                result += (Direction > 0) ? (forward[size-(k+1)]*interpolate(x0+(size-(k+1))*h, x, y)) : (backward[size-(k+1)]*interpolate(x0-(size-(k+1))*h, x, y));
            }
        }
    }
    return result/h;
}

// Savitzky-Golay filter
/// \brief          Savitzky-Golay filter.
/// \details        Computes the smoothed value of the n-th derivative at the 
///                 given position using a fourth order Savitzky-Golay algorithm. 
///                 The greater the neighbourhood, the greater the computing 
///                 distance is.
/// \tparam         Derivative Order of the derivative.
/// \tparam         Type Data type.
/// \tparam         Container Container type.
/// \tparam         T Conversion type.
/// \param[in]      x0 Differentiation abscissa.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \param[in]      neighbourhood Computing distance.
/// \return         Filtered derivative at the provided abscissa.
template <int Derivative, typename Type, class Container, typename T, class> 
inline Type Utility::filter(const Type x0, const Container& x, const Container& y, const unsigned int neighbourhood)
{
    static const unsigned int one = 1;
    static const unsigned int order = 4;
    static const std::array<Type, 4> normalization({{T(231), T(1188), T(462), T(198)}});
    static const std::array<Type, order*2+1> zeroth({{T(-21), T(14), T(39), T(54), T(59), T(54), T(39), T(14), T(-21)}});
    static const std::array<Type, order*2+1> first({{T(86), T(-142), T(-193), T(-126), T(0), T(126), T(193), T(142), T(-86)}});
    static const std::array<Type, order*2+1> second({{T(28), T(7), T(-8), T(-17), T(-20), T(-17), T(-8), T(7), T(28)}});
    static const std::array<Type, order*2+1> third({{T(-14), T(7), T(13), T(9), T(0), T(-9), T(-13), T(-7), T(14)}});
    static const std::array<Type, order*2+1> coeff = (Derivative == 0) ? (zeroth) : ((Derivative == 1) ? (first) : ((Derivative == 2) ? (second) : (third)));
    const long long int n = std::min(x.size(), y.size());
    const long long int i = std::distance(std::begin(x), std::upper_bound(std::begin(x), std::begin(x)+n, x0, [](const Type x1, const Type x2){return x1 < x2;}));
    const long long int j = i+((n > 1)*((i <= 0)-(i >= n)));
    const Type h = (n > 1) ? ((x[j]-x[j-1])*std::max(one, neighbourhood)) : (Type());
    Type result = Type();
    if (n > 1) {
        for (unsigned int k = 0; k < order; ++k) {
            result += coeff[k]*interpolate(x0-(order-k)*h, x, y)+coeff[order+order-k]*interpolate(x0+(order-k)*h, x, y);
        }
        result += coeff[order]*interpolate(x0, x, y);
    }
    return result/(normalization[Derivative]*std::pow(h, Type(Derivative)));
}
// -------------------------------------------------------------------------- //



// -------------------------------- EVOLUTION ------------------------------- //
// Reverse
/// \brief          Reverse.
/// \details        Reverses a container by substracting a value.
/// \tparam         Container Container type.
/// \tparam         Type Data type.
/// \param[in]      container Container.
/// \param[in]      value Value to be substracted.
/// \return         Reversed vector.
template <class Container, typename Type, class> 
Container Utility::reverse(const Container& container, const Type value)
{
    Container result(container);
    parallelize(result.begin(), result.end(), [=, &value](Type& element){element -= value;});
    return result;
}

// Kernel smoother
/// \brief          Kernel smoother.
/// \details        Smoothes a serie of data using a kernel smoother and a 
///                 window function on the first neigbours.
/// \tparam         Container Container type.
/// \tparam         Type Data type.
/// \tparam         Function Function type taking two data as argument.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \param[in]      kernel Kernel function.
/// \param[in]      window Cut on the provided number of neighbours.
/// \return         Smoothed ordinates.
template <class Container, typename Type, class Function, class> 
Container Utility::smooth(const Container& x, const Container& y, Function&& kernel, const unsigned int window)
{
    static const Type zero = Type();
    const unsigned int size = std::min(x.size(), y.size());
    const unsigned int win = (window > 0) ? (window) : (size);
    Container result(size);
    parallelize(size, [=, &x, &y, &win, &result](const unsigned int i){const unsigned int n = std::min(win, std::min(i, size-1-i)); Type w = zero; Type num = zero; Type den = zero; for (unsigned int j = i-n; j <= i+n; ++j) {w = kernel(x[i], x[j]); num += w*y[j]; den += w;} result[i] = num/den;});
    return result;
}

// Integration
/// \brief          Integration.
/// \details        Computes the integral of the ordinates regarding to the
///                 abscissae and starting from the specified value.
/// \tparam         Container Container type.
/// \tparam         Type Data type.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \return         Integral of ordinates regarding to the abscissae.
template <class Container, typename Type, class>
Container Utility::integrate(const Container& x, const Container& y, const Type value)
{
    const unsigned int size = std::min(x.size(), y.size());
    Container result(size, value);
    for (unsigned int i = 1; i < size; ++i) {
        result[i] = result[i-1]+y[i]*(x[i]-x[i-1]);
    }
    return result;
}

// Derivative
/// \brief          Derivative.
/// \details        Computes the derivative of the ordinates regarding to the
///                 abscissae.
/// \tparam         Direction Sign for centered, backward or forward 
///                 differentiation.
/// \tparam         Container Container type.
/// \tparam         Type Data type.
/// \param[in]      x Abscissae.
/// \param[in]      y Ordinates.
/// \param[in]      neighbourhood Computing distance.
/// \return         Derivative of ordinates regarding to the abscissae.
template <int Direction, class Container, typename Type, class>
Container Utility::derive(const Container& x, const Container& y, const unsigned int neighbourhood)
{
    const unsigned int size = std::min(x.size(), y.size());
    Container result(size);
    parallelize(size, [=, &x, &y, &result, &neighbourhood](const unsigned int i){result[i] = differentiate<Direction>(x[i], x, y, neighbourhood);});
    return result;  
}


// Inversion de matrices
/// \brief          Inverse matrix.
/// \details        Inverse 2D matrix.
/// \tparam         Type type of matrix components.
/// \param[in]      A array 2D Matrix.
/// return  	    Array 2D inverted matrix
template < class Type >  
inline std::array< std::array<Type, 2>, 2 > Utility::invMatrix2d(const std::array< std::array<Type, 2>, 2 >& A ) 
{
    std::array< std::array<Type, 2>, 2 > result;
    const double invdet = 1./(A[0][0]*A[1][1]-A[0][1]*A[1][0]);
    result[0][0] = A[1][1]*invdet;
    result[0][1] = -A[0][1]*invdet;
    result[1][0] = -A[1][0]*invdet;
    result[1][1] = A[0][0]*invdet;
   return result;
}

/// \brief          Inverse matrix.
/// \details        Inverse 3D matrix.
/// \tparam         Type type of matrix components.
/// \param[in]      A array 3D Matrix.
/// return  	    Array 3D inverted matrix
template < class Type >  
inline std::array< std::array<Type, 3>, 3 > Utility::invMatrix3d(const std::array< std::array<Type, 3>, 3 >& A ) 
{
    std::array< std::array<Type, 3>, 3 > result;
    const double invdet = 1./(A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]));
    result[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
    result[0][1] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
    result[0][2] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
    result[1][0] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
    result[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
    result[1][2] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
    result[2][0] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
    result[2][1] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
    result[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;

   return result;
}

// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Utility.
/// \return         0 if no error.
int Utility::example()
{
    // Initialize
    std::cout<<"BEGIN = Utility::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<double, 3> first = std::array<double, 3>({{4., 8., 15.}});
    std::array<double, 3> second = std::array<double, 3>({{16., 23., 42.}});
    std::vector<int> vector(42);
    std::vector<double> x(42);
    std::vector<double> y(42);
    magrathea::HyperCube<3> cube = magrathea::HyperCube<3>::unit();
    magrathea::HyperSphere<3> sphere = magrathea::HyperSphere<3>::unit();
    Cone<> cone(first, second, 0.42);
    std::iota(x.begin(), x.end(), 0);
    std::iota(y.begin(), y.end(), 42);

    // Construction
    Utility utility;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"Utility() : "                                 ; Utility(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"utility = Utility() : "                       ; utility = Utility(); std::cout<<std::endl;
    
    // Parallelization
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Parallelization : "                                                                 <<std::endl;
    std::cout<<std::setw(width*2)<<"utility.parallelize<1>(42, [](unsigned int){;}) : "                                 <<utility.parallelize<1>(42, [](unsigned int){;})<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.parallelize<1>(0., 42., 0.5, [](double){;}) : "                             <<utility.parallelize<1>(0., 42., 0.5, [](double){;})<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.parallelize<1>(vector.begin(), vector.end(), [](int& d){d += 42;}) : "      <<utility.parallelize<1>(vector.begin(), vector.end(), [](int& d){d += 42;})<<std::endl;

    // Geometry
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Geometry : "                                                                        <<std::endl;
    std::cout<<std::setw(width*2)<<"utility.distance<3>(first, second) : "                                              <<utility.distance<3>(first, second)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.dot<3>(first, second) : "                                                   <<utility.dot<3>(first, second)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.cross<3>(first, second)[0] : "                                              <<utility.cross<3>(first, second)[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.join<3>(first, second)[0] : "                                               <<utility.join<3>(first, second)[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.apply<std::plus<double>, 3>(first, second)[0] : "                           <<utility.apply<std::plus<double>, 3>(first, second)[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.radius(sphere) : "                                                          <<utility.radius(sphere)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.radius(cube) : "                                                            <<utility.radius(cube)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.cubify(sphere, 2) : "                                                       <<utility.cubify(sphere, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.spherify(cube, 2) : "                                                       <<utility.spherify(cube, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.collide(cube, sphere) : "                                                   <<utility.collide(cube, sphere)<<std::endl;
    std::cout<<std::setw(width*2)<<"utility.collide(cube, cone, 0.05) : "                                                     <<utility.collide(cube, cone, 0.05)<<std::endl;

    // Interpolation
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Interpolation : "                             <<std::endl;
    std::cout<<std::setw(width)<<"utility.interpolate(0.42, x, y) : "           <<utility.interpolate(0.42, x, y)<<std::endl;
    std::cout<<std::setw(width)<<"utility.interpolate(0.42, x, y, y) : "        <<utility.interpolate(0.42, x, y, y)<<std::endl;
    std::cout<<std::setw(width)<<"utility.reinterpolate(x, x, y).size() : "     <<utility.reinterpolate(x, x, y).size()<<std::endl;
    std::cout<<std::setw(width)<<"utility.differentiate(0.42, x, y, 10) : "     <<utility.differentiate(0.42, x, y, 10)<<std::endl;
    std::cout<<std::setw(width)<<"utility.filter(0.42, x, y, 10) : "            <<utility.filter(0.42, x, y, 10)<<std::endl;

    // Evolution
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Evolution : "                                                                                                               <<std::endl;
    std::cout<<std::setw(width*3)<<"utility.reverse(x, 42.).size() : "                                                                                          <<utility.reverse(x, 42.).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"utility.smooth(x, y, [](double xi, double xj){return std::exp(-((xi-xj)*(xi-xj))/2*std::pow(0.01, 2));}, 10).size() : "     <<utility.smooth(x, y, [](double xi, double xj){return std::exp(-((xi-xj)*(xi-xj))/2*std::pow(0.01, 2));}, 10).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"utility.integrate(x, y).size() : "                                                                                          <<utility.integrate(x, y).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"utility.derive(x, y, 10).size() : "                                                                                         <<utility.derive(x, y, 10).size()<<std::endl;
        
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Utility::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
#endif // UTILITY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
