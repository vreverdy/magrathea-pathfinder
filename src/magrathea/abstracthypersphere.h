/* ************************** ABSTRACTHYPERSPHERE *************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          AbstractHyperSphere
// DESCRIPTION :    Abstract function provider for n-dimensional spheres
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           abstracthypersphere.h
/// \brief          Abstract function provider for n-dimensional spheres
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABSTRACTHYPERSPHERE_H_INCLUDED
#define ABSTRACTHYPERSPHERE_H_INCLUDED
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
// Abstract function provider for n-dimensional spheres
/// \brief          Abstract function provider for n-dimensional spheres.
/// \details        Provides a common base for n-dimensional spheres thanks to 
///                 CRTP. To use it, one has to derive from this class and pass 
///                 the derived class itself as the CRTP parameter. The derived 
///                 classes should provide two immutable functions : 
///                 <ul><li>position()</li>
///                 <li>extent()</li></ul>
///                 in order to get the position of the center along one 
///                 coordinate and the radius length of the hypersphere.
/// \tparam         Crtp Derived CRTP class.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
class AbstractHyperSphere
: public AbstractShape
{
    // Setup
    static_assert(std::is_same<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type, Scalar>::value, "ERROR = AbstractHyperSphere : non matching vector and scalar types");
    
    // Protected lifecycle
    /// \name           Protected lifecycle
    //@{
    protected: 
        inline ~AbstractHyperSphere();
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
        constexpr Scalar radius() const;
        constexpr Scalar diameter() const;
        constexpr Scalar volume() const;
        constexpr Scalar surface() const;
    //@}

    // Distribution
    /// \name           Distribution
    //@{
    public:
        template <unsigned int Subdimension = Dimension, class = typename std::enable_if<(Subdimension+1 == Dimension) || (Subdimension == Dimension)>::type> Vector random() const;
        template <unsigned int Subdimension = Dimension, class Engine, class Distribution, class = typename std::enable_if<((Subdimension+1 == Dimension) || (Subdimension == Dimension)) && (std::decay<Engine>::type::min() != std::decay<Engine>::type::max()) && (!std::is_void<typename std::decay<Distribution>::type::result_type>::value)>::type> Vector random(Engine&& engine, Distribution&& distribution) const;
        template <unsigned int Subdimension = Dimension, typename Iterator, class = typename std::enable_if<(Subdimension+1 == Dimension) && (Subdimension >= 1) && (Subdimension <= 2) && (std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Vector>::value)>::type> std::pair<Scalar, Scalar> uniform(const Iterator& first, const Iterator& last) const;
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
    //@}

    // Helpers
    /// \name           Helpers
    //@{
    public:
        template <unsigned int OtherDimension = Dimension, typename OtherScalar = Scalar> static constexpr OtherScalar sn();
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
inline AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::~AbstractHyperSphere()
= default;
// -------------------------------------------------------------------------- //



// -------------------------------- POSITION -------------------------------- //
// Center coordinate
/// \brief          Center coordinate.
/// \details        Computes the specified coordinate of the center of the 
///                 hypersphere.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate \f$x_{i}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::center(const unsigned int idim) const
{
    return static_cast<const Crtp&>(*this).position(idim);
}

// Center vector
/// \brief          Center vector.
/// \details        Computes the position vector of the center of the 
///                 hypersphere.
/// \return         The position vector \f$\vec{x}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline Vector AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::center() const
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
///                 of the hypersphere.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate \f$x_{i}-r\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::minimum(const unsigned int idim) const
{
    return static_cast<const Crtp&>(*this).position(idim)-static_cast<const Crtp&>(*this).extent();
}

// Minimum vector
/// \brief          Minimum vector.
/// \details        Computes the position vector of the minimum boundary of the 
///                 hypersphere.
/// \return         The position vector \f$\vec{x}-\vec{r}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline Vector AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::minimum() const
{
    Vector result = Vector();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result[idim] = static_cast<const Crtp&>(*this).position(idim)-static_cast<const Crtp&>(*this).extent();
    }
    return result;
}

// Maximum coordinate
/// \brief          Maximum coordinate.
/// \details        Computes the specified coordinate of the maximum boundary 
///                 of the hypersphere.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate \f$x_{i}+r\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::maximum(const unsigned int idim) const
{
    return static_cast<const Crtp&>(*this).position(idim)+static_cast<const Crtp&>(*this).extent();
}

// Maximum vector
/// \brief          Maximum vector.
/// \details        Computes the position vector of the maximum boundary of the 
///                 hypersphere.
/// \return         The position vector \f$\vec{x}+\vec{r}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
inline Vector AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::maximum() const
{
    Vector result = Vector();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result[idim] = static_cast<const Crtp&>(*this).position(idim)+static_cast<const Crtp&>(*this).extent();
    }
    return result;
}
// -------------------------------------------------------------------------- //



// -------------------------------- MEASURES -------------------------------- //
// Radius
/// \brief          Radius.
/// \details        Computes the radius of the hypersphere.
/// \return         \f$r\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::radius() const
{
    return static_cast<const Crtp&>(*this).extent();
}

// Diameter
/// \brief          Diameter.
/// \details        Computes the diameter of the hypersphere.
/// \return         \f$2r\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::diameter() const
{
    return 2*static_cast<const Crtp&>(*this).extent();
}

// Volume
/// \brief          Volume.
/// \details        Computes the volume of the hypersphere.
/// \return         \f$\frac{s_{n}r^{n}}{n}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::volume() const
{
    return (sn<Dimension>()*pow<Dimension>(static_cast<const Crtp&>(*this).extent()))/Scalar(Dimension);
}

// Outer surface
/// \brief          Outer surface.
/// \details        Computes the total outer surface of the hypercube.
/// \return         \f$s_{n}r^{n-1}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr Scalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::surface() const
{
    return sn<Dimension>()*pow<Dimension-(Dimension > 0)>(static_cast<const Crtp&>(*this).extent())*(Dimension > 0);
}
// -------------------------------------------------------------------------- //



// ------------------------------ DISTRIBUTION ------------------------------ //
// Basic random location in the hypersphere
/// \brief          Basic random location in the hypersphere.
/// \details        Generates a random location, located in the volume or on the
///                 surface of the hypersphere. For example for a subdimension 
///                 of 2 of a 3-dimensional hypersphere, the function will 
///                 generates a random point located on the surface of the 
///                 sphere. 
/// \tparam         Subdimension Dimension space.
/// \return         Random position vector.
/// \warning        As the internal engine is a static one, do not use this 
///                 function in parallel.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension, class> 
Vector AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::random() const
{
    static const unsigned long long int ilim = std::numeric_limits<unsigned long long int>::max();
    static const Scalar slim = ilim;
    static std::random_device device;
    static std::mt19937 engine(device());
    std::uniform_int_distribution<unsigned long long int> distribution(0, ilim);
    std::normal_distribution<typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type> normal(0, 1);
    typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type array[Dimension];
    typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type coefficient = 0;
    Vector result = Vector();
    if ((Subdimension+1 == Dimension) || (Subdimension == Dimension)) {
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            array[idim] = normal(engine);
            coefficient += pow<2>(array[idim]);
        }
        coefficient = (static_cast<const Crtp&>(*this).extent()/(rt<2>(coefficient)*slim))*((Subdimension == Dimension) ? (distribution(engine)) : (ilim));
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            result[idim] = static_cast<const Crtp&>(*this).position(idim)+coefficient*array[idim]; 
        }
    }
    return result;
}

// Generic random location in the hypersphere
/// \brief          Generic random location in the hypersphere.
/// \details        Generates a random location, located in the volume or on the
///                 surface of the hypersphere. For example for a subdimension 
///                 of 2 of a 3-dimensional hypersphere, the function will 
///                 generates a random point located on the surface of the 
///                 sphere. As this function uses the passed random engine and 
///                 distribution it is completely thread safe.
/// \tparam         Subdimension Dimension space.
/// \tparam         Engine (Random engine type.)
/// \tparam         Distribution (Random distribution type.)
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         Random position vector.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension, class Engine, class Distribution, class> 
Vector AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::random(Engine&& engine, Distribution&& distribution) const
{
    static const Scalar one = 1;
    const Scalar a = one/(distribution.max()-distribution.min());
    const Scalar b = -a*distribution.min();
    std::normal_distribution<typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type> normal(0, 1);
    typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type array[Dimension];
    typename std::conditional<std::is_floating_point<Scalar>::value, Scalar, double>::type coefficient = 0;
    Vector result = Vector();
    if ((Subdimension+1 == Dimension) || (Subdimension == Dimension)) {
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            array[idim] = normal(engine);
            coefficient += pow<2>(array[idim]);
        }
        coefficient = (static_cast<const Crtp&>(*this).extent()/(rt<2>(coefficient)))*((Subdimension == Dimension) ? (a*distribution(engine)+b) : (one));
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            result[idim] = static_cast<const Crtp&>(*this).position(idim)+coefficient*array[idim]; 
        }
    }
    return result;
}

// Uniform distribution
/// \brief          Uniform distribution.
/// \details        Generates a random distribution of points located on the 
///                 surface of the hypersphere. Currently the function only 
///                 works in two or three dimension corresponding to a value 
///                 of one or two for the subdimension. In three dimensions,
///                 it uses a spiral to approximately distribute the points. 
///                 The function eventually returns a pair of distance 
///                 corresponding to the minimum and maximum distances between
///                 two generated points.
/// \tparam         Subdimension Dimension space.
/// \tparam         Iterator (Pointer or iterator to vector types.)
/// \param[in]      first Beginning of the interval.
/// \param[in]      end End of the interval.
/// \return         Minimum and maximum distance between two points.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int Subdimension, typename Iterator, class> 
std::pair<Scalar, Scalar> AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::uniform(const Iterator& first, const Iterator& last) const
{
    static const Scalar zero = 0;
    static const Scalar one = 1;
    static const Scalar two = one+one;
    static const Scalar twopi = two*pi<1, Scalar>();
    static const Scalar spiral = twopi*(two-golden<1, Scalar>());
    const unsigned int n = (first < last) ? (last-first) : (0);
    const Scalar extent = static_cast<const Crtp&>(*this).extent();
    Vector position = Vector();
    Scalar tmp = zero;
    Vector tmpv = Vector();
    std::pair<Scalar, Scalar> result(zero, zero);
    if (n > 0) {
        for (unsigned int idim = 0; idim < Dimension; ++idim) {
            position[idim] = static_cast<const Crtp&>(*this).position(idim);
        }
        if (Subdimension+1 == Dimension) {
            if (Dimension == 2) {
                result.first = two*extent;
                tmp = twopi/Scalar(n);
                for (unsigned int i = 0; i < n; ++i) {
                    (*(first+i))[0] = position[0]+(extent*std::cos(tmp*Scalar(i))); 
                    (*(first+i))[1] = position[1]+(extent*std::sin(tmp*Scalar(i)));
                }
            } else if (Dimension == 3) {
                result.first = two*extent;
                tmp = two/Scalar(n);
                for (unsigned int i = 0; i < n; ++i) {
                    tmpv[0] = (Scalar(i)*tmp)-one+(tmp/two);
                    tmpv[1] = rt<2>(one-pow<2>(tmpv[0]));
                    tmpv[2] = Scalar(i)*spiral;
                    (*(first+i))[0] = position[0]+(extent*tmpv[1]*std::cos(tmpv[2])); 
                    (*(first+i))[1] = position[1]+(extent*tmpv[0]);
                    (*(first+i))[2] = position[2]+(extent*tmpv[1]*std::sin(tmpv[2]));
                }
            }
        }
        if ((n > 1) && ((result.first > zero) || (result.first < zero))) {
            for (unsigned int i = 0; i < n; ++i) {
                for (unsigned int j = 0; j < n; ++j) {
                    if (i != j) {
                        tmp = zero;
                        for (unsigned int idim = 0; idim < Dimension; ++idim) {
                            tmp += pow<2>((*(first+i))[idim]-(*(first+j))[idim]);
                        }
                        tmp = rt<2>(tmp);
                        result.first = (result.first < tmp) ? (result.first) : (tmp);
                        result.second = (result.second > tmp) ? (result.second) : (tmp);
                    }
                }
            }
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// -------------------------------- COLLISION ------------------------------- //
// Point inside
/// \brief          Point inside.
/// \details        Checks whether a point is inside the hypersphere.
/// \tparam         OtherVector Other position vector type.
/// \param[in]      point Position of the point.
/// \return         True if the point is inside the hypersphere, false 
///                 otherwise.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <class OtherVector, class> 
inline bool AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::inside(const OtherVector& point) const
{
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result += pow<2>(static_cast<const Crtp&>(*this).position(idim)-point[idim]);
    }
    return !(rt<2>(result) > static_cast<const Crtp&>(*this).extent());
}

// Point outside
/// \brief          Point outside.
/// \details        Checks whether a point is outside the hypersphere.
/// \tparam         OtherVector Other position vector type.
/// \param[in]      point Position of the point.
/// \return         True if the point is outside the hypersphere, false 
///                 otherwise.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <class OtherVector, class> 
inline bool AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::outside(const OtherVector& point) const
{
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        result += pow<2>(static_cast<const Crtp&>(*this).position(idim)-point[idim]);
    }
    return !(rt<2>(result) < static_cast<const Crtp&>(*this).extent());
}
// -------------------------------------------------------------------------- //



// -------------------------------- CONSTANTS ------------------------------- //
// Number of space dimension
/// \brief          Number of space dimension.
/// \details        Computes the number of space dimension of the hypercube.
/// \return         Dimension. 
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
constexpr unsigned int AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::dimension()
{
    return Dimension;
}
// -------------------------------------------------------------------------- //



// --------------------------------- HELPERS -------------------------------- //
// Surface of a n-dimensional unit sphere
/// \brief          Surface of a n-dimensional unit sphere.
/// \details        Computes the surface of a n-dimensional sphere with a unit
///                 radius.
/// \tparam         OtherDimension Other number of space dimension.
/// \tparam         OtherScalar (Other scalar data type.)
/// \return         \f$s_{n}\f$.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int OtherDimension, typename OtherScalar>
constexpr OtherScalar AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::sn()
{
    return (OtherDimension%2 == 0) ? (((OtherDimension > 0)*2*pow<OtherDimension/2>(pi<1, OtherScalar>()))/factorial<(OtherDimension/2)-(OtherDimension > 1), OtherScalar>()) : ((pow<(OtherDimension+(OtherDimension > 0))/2>(2)*pow<(OtherDimension-(OtherDimension > 0))/2>(pi<1, OtherScalar>()))/factorial<OtherDimension-2*(OtherDimension > 1), true, OtherScalar>());
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of AbstractHyperSphere.
/// \return         0 if no error.
template <class Crtp, unsigned int Dimension, class Vector, typename Scalar>
int AbstractHyperSphere<Crtp, Dimension, Vector, Scalar>::example()
{
    std::cout<<"BEGIN = AbstractHyperSphere::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = AbstractHyperSphere::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = AbstractHyperSphere::example()"<<std::endl;
    return 1;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABSTRACTHYPERSPHERE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

