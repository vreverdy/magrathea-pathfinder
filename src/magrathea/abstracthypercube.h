/* *************************** ABSTRACTHYPERCUBE **************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          AbstractHyperCube
// DESCRIPTION :    Abstract function provider for n-dimensional cubes
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           abstracthypercube.h
/// \brief          Abstract function provider for n-dimensional cubes
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABSTRACTHYPERCUBE_H_INCLUDED
#define ABSTRACTHYPERCUBE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <random>
#include <algorithm>
// Include libs
// Include project
#include "abstractshape.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Abstract function provider for n-dimensional cubes
/// \brief          Abstract function provider for n-dimensional cubes.
/// \details        Provides a common base for n-dimensional cubes thanks to 
///                 CRTP. To use it, one has to derive from this class and pass 
///                 the derived class itself as the CRTP parameter. The derived 
///                 classes should provide two immutable functions : 
///                 <ul><li>position()</li>
///                 <li>extent()</li></ul>
///                 in order to get the position of the center along one 
///                 coordinate and the edge length of the hypercube.
/// \tparam         Crtp Derived CRTP class.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
class AbstractHyperCube
: public AbstractShape
{
    // Setup
    static_assert(std::is_same<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type, Scalar>::value, "ERROR = AbstractHyperCube : non matching vector and scalar types");

    // Protected lifecycle
    /// \name           Protected lifecycle
    //@{
    protected: 
        inline ~AbstractHyperCube();
    //@}
    
    // Position
    /// \name           Position
    //@{
    public:
        constexpr Scalar center(const unsigned int idim) const;
        inline Vector center() const;
        constexpr Scalar minimum(const unsigned int idim) const;
        inline Vector minimum() const;
        constexpr Scalar maximum(const unsigned int idim) const;
        inline Vector maximum() const;
    //@}

    // Measures
    /// \name           Measures
    //@{
    public:
        constexpr Scalar length() const;
        constexpr Scalar volume() const;
        template <unsigned int Subdimension = Dimension-(Dimension != 0)> constexpr Scalar surface() const;
        template <unsigned int Subdimension = Dimension-(Dimension != 0)> constexpr Scalar area() const;
        template <unsigned int Subdimension = Dimension> constexpr Scalar diagonal() const;
    //@}
    
    // Distribution
    /// \name           Distribution
    //@{
    public:
        template <unsigned int Subdimension = Dimension, class = typename std::enable_if<(Subdimension != 0) || (Subdimension == 0)>::type> Vector random() const;
        template <unsigned int Subdimension = Dimension, class Engine, class Distribution, class = typename std::enable_if<((Subdimension != 0) || (Subdimension == 0)) && (std::decay<Engine>::type::min() != std::decay<Engine>::type::max()) && (!std::is_void<typename std::decay<Distribution>::type::result_type>::value)>::type> Vector random(Engine&& engine, Distribution&& distribution) const;
    //@}
    
    // Collision
    /// \name           Collision
    //@{
    public:
        template <class OtherVector, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<OtherVector>()[0])>::type>::type, Scalar>::value>::type> inline bool inside(const OtherVector& point) const;
        template <class OtherVector, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<OtherVector>()[0])>::type>::type, Scalar>::value>::type> inline bool outside(const OtherVector& point) const;
    //@}
    
    // Constants
    /// \name           Constants
    //@{
    public:
        static constexpr unsigned int dimension();
        template <unsigned int Subdimension = Dimension> static constexpr unsigned int elements();
        template <unsigned int FirstSubdimension = 0, unsigned int LastSubdimension = Dimension-(Dimension != 0)> static constexpr unsigned int subelements();
    //@}

    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// -------------------------- PROTECTED LIFECYCLE --------------------------- //
// Protected destructor
/// \brief          Protected destructor.
/// \details        Avoids direct instantiation of the class, and only allows 
///                 it through its derived children.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::~AbstractHyperCube()
= default;
// -------------------------------------------------------------------------- //



// -------------------------------- POSITION -------------------------------- //
// Center coordinate
/// \brief          Center coordinate.
/// \details        Computes the specified coordinate of the center of the 
///                 hypercube.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate \f$x_{i}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::center(const unsigned int idim) const
{
    return static_cast<const Crtp&>(*this).position(idim);
}

// Center vector
/// \brief          Center vector.
/// \details        Computes the position vector of the center of the 
///                 hypercube.
/// \return         The position vector \f$\vec{x}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline Vector AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::center() const
{
    Vector result = Vector();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result[idim] = static_cast<const Crtp&>(*this).position(idim);
    }
    return result;
}

// Minimum coordinate
/// \brief          Minimum coordinate.
/// \details        Computes the specified coordinate of the minimum boundary 
///                 of the hypercube.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate \f$x_{i}-\frac{l}{2}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::minimum(const unsigned int idim) const
{
    return static_cast<const Crtp&>(*this).position(idim)-static_cast<const Crtp&>(*this).extent()/2;
}

// Minimum vector
/// \brief          Minimum vector.
/// \details        Computes the position vector of the minimum boundary of the 
///                 hypercube.
/// \return         The position vector \f$\vec{x}-\frac{\vec{l}}{2}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline Vector AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::minimum() const
{
    Vector result = Vector();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result[idim] = static_cast<const Crtp&>(*this).position(idim)-static_cast<const Crtp&>(*this).extent()/2;
    }
    return result;
}

// Maximum coordinate
/// \brief          Maximum coordinate.
/// \details        Computes the specified coordinate of the maximum boundary 
///                 of the hypercube.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate \f$x_{i}+\frac{l}{2}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::maximum(const unsigned int idim) const
{
    return static_cast<const Crtp&>(*this).position(idim)+static_cast<const Crtp&>(*this).extent()/2;
}

// Maximum vector
/// \brief          Maximum vector.
/// \details        Computes the position vector of the maximum boundary of the 
///                 hypercube.
/// \return         The position vector \f$\vec{x}+\frac{\vec{l}}{2}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline Vector AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::maximum() const
{
    Vector result = Vector();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result[idim] = static_cast<const Crtp&>(*this).position(idim)+static_cast<const Crtp&>(*this).extent()/2;
    }
    return result;
}
// -------------------------------------------------------------------------- //



// -------------------------------- MEASURES -------------------------------- //
// Length
/// \brief          Length.
/// \details        Computes the edge length of the hypercube.
/// \return         \f$l\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::length() const
{
    return static_cast<const Crtp&>(*this).extent()/2;
}

// Volume
/// \brief          Volume.
/// \details        Computes the volume of the hypercube.
/// \return         \f$l^{D}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::volume() const
{
    return pow<Dimension>(static_cast<const Crtp&>(*this).extent());
}

// Outer surface
/// \brief          Outer surface.
/// \details        Computes the total outer surface of the hypercube.
/// \tparam         Subdimension Dimension of elements.
/// \return         \f$2^{D-d}\binom{D}{d} \times l^{D}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension> 
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::surface() const
{
    return elements<Subdimension>()*pow<Subdimension>(static_cast<const Crtp&>(*this).extent());
}

// Element area
/// \brief          Element area.
/// \details        Computes the surface of a single element of the hypercube.
/// \tparam         Subdimension Dimension of elements.
/// \return         \f$l^{D}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension> 
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::area() const
{
    return (Subdimension <= Dimension) ? (pow<Subdimension>(static_cast<const Crtp&>(*this).extent())) : (0);
}

// Diagonal
/// \brief          Diagonal.
/// \details        Computes the length of the diagonal in the specified 
///                 subdimension.
/// \tparam         Subdimension Dimension of elements.
/// \return         \f$\sqrt{D} \times l\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension>
constexpr Scalar AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::diagonal() const
{
    return (Subdimension <= Dimension) ? (rt<2>(Scalar(Subdimension))*static_cast<const Crtp&>(*this).extent()) : (0);
}
// -------------------------------------------------------------------------- //



// ------------------------------ DISTRIBUTION ------------------------------ //
// Basic random location in the hypercube
/// \brief          Basic random location in the hypercube.
/// \details        Generates a random location, located on the subdimensional 
///                 elements of the hypercube. For example for a subdimension of
///                 2 of a 3-dimensional hypercube, the function will generates
///                 a random point located on the surface of the cube. 
/// \tparam         Subdimension Dimension space.
/// \return         Random position vector.
/// \warning        As the internal engine is a static one, do not use this 
///                 function in parallel.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension, class> 
Vector AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::random() const
{
    static const long long int ilim = std::numeric_limits<long long int>::max();
    static const Scalar slim = ilim;
    static const Scalar stwo = 2;
    static std::random_device device;
    static std::mt19937 engine(device());
    std::uniform_int_distribution<long long int> distribution(-ilim, ilim);
    Scalar array[Dimension];
    Vector result = Vector();
    if (Subdimension <= Dimension) {
        for (unsigned int idim = 0; idim < Subdimension; ++idim) {
            array[idim] = (Scalar(distribution(engine))/slim)*(static_cast<const Crtp&>(*this).extent()/stwo);
        }
        for (unsigned int idim = Subdimension; idim < Dimension; ++idim) {
            array[idim] = ((Scalar(distribution(engine) >= 0)) ? (1) : (-1))*(static_cast<const Crtp&>(*this).extent()/stwo);
        }
        std::shuffle(&array[0], &array[0]+Dimension, engine);
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            result[idim] = static_cast<const Crtp&>(*this).position(idim)+array[idim];
        }
    }
    return result;
}

// Generic random location in the hypercube
/// \brief          Generic random location in the hypercube.
/// \details        Generates a random location, located on the subdimensional 
///                 elements of the hypercube. For example for a subdimension of
///                 2 of a 3-dimensional hypercube, the function will generates
///                 a random point located on the surface of the cube. As this
///                 function uses the passed random engine and distribution it
///                 is completely thread safe.
/// \tparam         Subdimension Dimension space.
/// \tparam         Engine (Random engine type.)
/// \tparam         Distribution (Random distribution type.)
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         Random position vector.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension, class Engine, class Distribution, class> 
Vector AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::random(Engine&& engine, Distribution&& distribution) const
{
    static const Scalar two = 2;
    const Scalar a = two/(distribution.max()-distribution.min());
    const Scalar b = (-distribution.max()-distribution.min())/(distribution.max()-distribution.min());
    Scalar array[Dimension];
    Vector result = Vector();
    if (Subdimension <= Dimension) {
        for (unsigned int idim = 0; idim < Subdimension; ++idim) {
            array[idim] = (a*distribution(engine)+b)*(static_cast<const Crtp&>(*this).extent()/two);
        }
        for (unsigned int idim = Subdimension; idim < Dimension; ++idim) {
            array[idim] = ((a*distribution(engine)+b >= 0) ? (1) : (-1))*(static_cast<const Crtp&>(*this).extent()/two);
        }
        std::shuffle(&array[0], &array[0]+Dimension, engine);
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            result[idim] = static_cast<const Crtp&>(*this).position(idim)+array[idim];
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// -------------------------------- COLLISION ------------------------------- //
// Point inside
/// \brief          Point inside.
/// \details        Checks whether a point is inside the hypercube.
/// \tparam         OtherVector Other position vector type.
/// \param[in]      point Position of the point.
/// \return         True if the point is inside the hypercube, false otherwise.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <class OtherVector, class> 
inline bool AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::inside(const OtherVector& point) const
{
    const Scalar half = static_cast<const Crtp&>(*this).extent()/2;
    bool result = true;
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result = ((result) && !(static_cast<const Crtp&>(*this).position(idim)-half > point[idim]) && !(point[idim] > static_cast<const Crtp&>(*this).position(idim)+half));
    }
    return result;
}

// Point outside
/// \brief          Point outside.
/// \details        Checks whether a point is outside the hypercube.
/// \tparam         OtherVector Other position vector type.
/// \param[in]      point Position of the point.
/// \return         True if the point is outside the hypercube, false otherwise.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <class OtherVector, class> 
inline bool AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::outside(const OtherVector& point) const
{
    const Scalar half = static_cast<const Crtp&>(*this).extent()/2;
    bool result = false;
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result = ((result) || !(point[idim] > static_cast<const Crtp&>(*this).position(idim)-half) || !(static_cast<const Crtp&>(*this).position(idim)+half > point[idim]));
    }
    return result;
}
// -------------------------------------------------------------------------- //



// -------------------------------- CONSTANTS ------------------------------- //
// Number of space dimension
/// \brief          Number of space dimension.
/// \details        Computes the number of space dimension of the hypercube.
/// \return         Dimension. 
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr unsigned int AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::dimension()
{
    return Dimension;
}

// Number of elements
/// \brief          Number of elements.
/// \details        Computes the number of elements of the specified 
///                 subdimension in the hypercube. For example, for a 
///                 specified subdimension 0, it computes the number of points.
/// \tparam         Subdimension Dimension of elements.
/// \return         \f$2^{D-d}\binom{D}{d}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension> 
constexpr unsigned int AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::elements()
{
    return (Subdimension <= Dimension) ? (pow<(Dimension-Subdimension)*(Subdimension <= Dimension)>(2)*combinations<Dimension*(Subdimension <= Dimension), Subdimension*(Subdimension <= Dimension), false, unsigned int>()) : (0);
}

// Sum of elements
/// \brief          Sum of elements.
/// \details        Computes the sum of the hypercube elements between the two
///                 given dimensions.
/// \tparam         First First dimension of elements.
/// \tparam         Last Last dimension of elements.
/// \return         \f$\sum_{d}2^{D-d}\binom{D}{d}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int First, unsigned int Last>
constexpr unsigned int AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::subelements()
{
    return (First != Last) ? (elements<First>()+subelements<First+((First <= Last) ? (1) : (-1))*(First != Last), Last>()) : (elements<First*(First == Last)>());
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of AbstractHyperCube.
/// \return         0 if no error.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
int AbstractHyperCube<Crtp, Dimension, Vector, Scalar>::example()
{
    std::cout<<"BEGIN = AbstractHyperCube::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = AbstractHyperCube::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = AbstractHyperCube::example()"<<std::endl;
    return 1;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABSTRACTHYPERCUBE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
