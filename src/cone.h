/* ********************************** CONE ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Cone
// DESCRIPTION :    Three-dimensional cone
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           cone.h
/// \brief          Three-dimensional cone
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
/// (michel-andres.breton@obspm.fr) \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef CONE_H_INCLUDED
#define CONE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <type_traits>
#include <utility>
// Include libs
// Include project
#include "magrathea/abstractshape.h"
#include "magrathea/abstractsubstance.h"
// Misc
// -------------------------------------------------------------------------- //

// ---------------------------------- CLASS --------------------------------- //
// Three-dimensional cone
/// \brief          Three-dimensional cone.
/// \details        Implementation of a basic cone in three dimensions.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
template <class Vector = std::array<double, 3>,
          typename Scalar =
              typename std::remove_cv<typename std::remove_reference<
                  decltype(std::declval<Vector>()[0])>::type>::type>
class Cone final : public magrathea::AbstractShape,
                   public magrathea::AbstractSubstance<Cone<Vector, Scalar>,
                                                       Vector, Vector, Scalar> {
    // Setup
public:
    using magrathea::AbstractSubstance<Cone<Vector, Scalar>, Vector, Vector,
                                       Scalar>::operator=;

    // Lifecycle
    /// \name           Lifecycle
    //@{
public:
    template <class... Misc>
    explicit inline Cone(Misc &&...misc);
    //@}

    // Data
    /// \name           Data
    //@{
public:
    template <
        unsigned int... Values, class... Misc,
        class Template =
            decltype(std::declval<magrathea::AbstractSubstance<
                         Cone<Vector, Scalar>, Vector, Vector, Scalar>>()
                         .template data<0, Values...>(std::declval<Misc>()...)),
        class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template vertex(Misc &&...misc);
    template <
        unsigned int... Values, class... Misc,
        class Template =
            decltype(std::declval<const magrathea::AbstractSubstance<
                         Cone<Vector, Scalar>, Vector, Vector, Scalar>>()
                         .template data<0, Values...>(std::declval<Misc>()...)),
        class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template vertex(Misc &&...misc) const;
    template <
        unsigned int... Values, class... Misc,
        class Template =
            decltype(std::declval<magrathea::AbstractSubstance<
                         Cone<Vector, Scalar>, Vector, Vector, Scalar>>()
                         .template data<1, Values...>(std::declval<Misc>()...)),
        class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template base(Misc &&...misc);
    template <
        unsigned int... Values, class... Misc,
        class Template =
            decltype(std::declval<const magrathea::AbstractSubstance<
                         Cone<Vector, Scalar>, Vector, Vector, Scalar>>()
                         .template data<1, Values...>(std::declval<Misc>()...)),
        class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template base(Misc &&...misc) const;
    template <
        unsigned int... Values, class... Misc,
        class Template =
            decltype(std::declval<magrathea::AbstractSubstance<
                         Cone<Vector, Scalar>, Vector, Vector, Scalar>>()
                         .template data<2, Values...>(std::declval<Misc>()...)),
        class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template angle(Misc &&...misc);
    template <
        unsigned int... Values, class... Misc,
        class Template =
            decltype(std::declval<const magrathea::AbstractSubstance<
                         Cone<Vector, Scalar>, Vector, Vector, Scalar>>()
                         .template data<2, Values...>(std::declval<Misc>()...)),
        class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template angle(Misc &&...misc) const;
    //@}

    // Position
    /// \name           Position
    //@{
public:
    inline Scalar direction(const unsigned int idim) const;
    inline Vector direction() const;
    //@}

    // Measures
    /// \name           Measures
    //@{
public:
    inline Scalar length() const;
    inline Scalar radius() const;
    inline Scalar diameter() const;
    inline Scalar circle() const;
    inline Scalar surface() const;
    inline Scalar volume() const;
    //@}

    // Collision
    /// \name           Collision
    //@{
public:
    template <class OtherVector,
              class = typename std::enable_if<std::is_convertible<
                  typename std::remove_cv<typename std::remove_reference<
                      decltype(std::declval<OtherVector>()[0])>::type>::type,
                  Scalar>::value>::type>
    inline bool inside(const OtherVector &point) const;
    template <class OtherVector,
              class = typename std::enable_if<std::is_convertible<
                  typename std::remove_cv<typename std::remove_reference<
                      decltype(std::declval<OtherVector>()[0])>::type>::type,
                  Scalar>::value>::type>
    inline bool outside(const OtherVector &point) const;
    //@}

    // Test
    /// \name           Test
    //@{
public:
    static int example();
    //@}
};
// -------------------------------------------------------------------------- //

// ------------------------------- LIFECYCLE -------------------------------- //
// Explicit generic constructor
/// \brief          Explicit generic constructor.
/// \details        Provides a generic interface to all constructors of the
///                 base class.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
template <class Vector, typename Scalar>
template <class... Misc>
inline Cone<Vector, Scalar>::Cone(Misc &&...misc)
    : magrathea::AbstractShape(),
      magrathea::AbstractSubstance<Cone<Vector, Scalar>, Vector, Vector,
                                   Scalar>(std::forward<Misc>(misc)...) {
    ;
}
// -------------------------------------------------------------------------- //

// ---------------------------------- DATA ---------------------------------- //
// Access to the vertex data
/// \brief          Access to the vertex data.
/// \details        Provides an access to the vertex data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Cone<Vector, Scalar>::vertex(Misc &&...misc) {
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the vertex data
/// \brief          Immutable access to the vertex data.
/// \details        Provides an immutable access to the vertex data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Cone<Vector, Scalar>::vertex(Misc &&...misc) const {
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Access to the base data
/// \brief          Access to the base data.
/// \details        Provides an access to the base data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Cone<Vector, Scalar>::base(Misc &&...misc) {
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the base data
/// \brief          Immutable access to the base data.
/// \details        Provides an immutable access to the base data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Cone<Vector, Scalar>::base(Misc &&...misc) const {
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Access to the angle data
/// \brief          Access to the angle data.
/// \details        Provides an access to the angle data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Cone<Vector, Scalar>::angle(Misc &&...misc) {
    return this->template data<2, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the angle data
/// \brief          Immutable access to the angle data.
/// \details        Provides an immutable access to the angle data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Cone<Vector, Scalar>::angle(Misc &&...misc) const {
    return this->template data<2, Values...>(std::forward<Misc>(misc)...);
}
// -------------------------------------------------------------------------- //

// -------------------------------- POSITION -------------------------------- //
// Direction coordinate
/// \brief          Direction coordinate.
/// \details        Computes the specified coordinate of the vector going from
///                 the vertex to the base center.
/// \param[in]      idim Index of the dimension.
/// \return         The coordinate of \f$b_{i}-v_{i}\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::direction(const unsigned int idim) const {
    return base(idim) - vertex(idim);
}

// Direction vector
/// \brief          Direction vector.
/// \details        Computes the vector going from the vertex to the base
///                 center.
/// \return         The vector \f$\vec{b}-\vec{v}\f$.
template <class Vector, typename Scalar>
inline Vector Cone<Vector, Scalar>::direction() const {
    Vector result = Vector();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result[idim] = base(idim) - vertex(idim);
    }
    return result;
}
// -------------------------------------------------------------------------- //

// -------------------------------- MEASURES -------------------------------- //
// Length
/// \brief          Length.
/// \details        Computes the height length of the cone.
/// \return         \f$h\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::length() const {
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result += pow<2>(base(idim) - vertex(idim));
    }
    return rt<2>(result);
}

// Radius
/// \brief          Radius.
/// \details        Computes the base radius of the cone.
/// \return         \f$r\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::radius() const {
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result += pow<2>(base(idim) - vertex(idim));
    }
    return std::tan(angle()) * rt<2>(result);
}

// Diameter
/// \brief          Diameter.
/// \details        Computes the base diameter of the cone.
/// \return         \f$d\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::diameter() const {
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result += pow<2>(base(idim) - vertex(idim));
    }
    return 2 * std::tan(angle()) * rt<2>(result);
}

// Circle
/// \brief          Circle.
/// \details        Computes the base area of the cone.
/// \return         \f$B\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::circle() const {
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result += pow<2>(base(idim) - vertex(idim));
    }
    return pi<1, Scalar>() * pow<2>(std::tan(angle()) * rt<2>(result));
}

// Surface
/// \brief          Surface.
/// \details        Computes the outer surface of the cone.
/// \return         \f$\pi\times r^{2}+\pi\times r \times h\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::surface() const {
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result += pow<2>(base(idim) - vertex(idim));
    }
    return (pi<1, Scalar>() * std::tan(angle())) * (std::tan(angle()) + 1) *
           result;
}

// Volume
/// \brief          Volume.
/// \details        Computes the volume of the cone.
/// \return         \f$\frac{\pi\times r^{2}\times h}{3}\f$.
template <class Vector, typename Scalar>
inline Scalar Cone<Vector, Scalar>::volume() const {
    Scalar result = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        result += pow<2>(base(idim) - vertex(idim));
    }
    result = rt<2>(result);
    return (pi<1, Scalar>() * pow<2>(std::tan(angle()) * result) * result) / 3;
}
// -------------------------------------------------------------------------- //

// -------------------------------- COLLISION ------------------------------- //
// Point inside
/// \brief          Point inside.
/// \details        Checks whether a point is inside the cone.
/// \tparam         OtherVector Other position vector type.
/// \param[in]      point Position of the point.
/// \return         True if the point is inside the cone, false otherwise.
template <class Vector, typename Scalar>
template <class OtherVector, class>
inline bool Cone<Vector, Scalar>::inside(const OtherVector &point) const {
    Scalar dot = Scalar();
    Scalar dist = Scalar();
    Scalar ref = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        dot += (base(idim) - vertex(idim)) * (point[idim] - vertex(idim));
        dist += pow<2>(point[idim] - vertex(idim));
        ref += pow<2>(base(idim) - vertex(idim));
    }
    return ((dot > Scalar()) &&
            !(std::acos(dot / (rt<2>(dist * ref))) > angle()) && !(dot > ref));
}

// Point outside
/// \brief          Point outside.
/// \details        Checks whether a point is outside the cone.
/// \tparam         OtherVector Other position vector type.
/// \param[in]      point Position of the point.
/// \return         True if the point is outside the cone, false otherwise.
template <class Vector, typename Scalar>
template <class OtherVector, class>
inline bool Cone<Vector, Scalar>::outside(const OtherVector &point) const {
    Scalar dot = Scalar();
    Scalar dist = Scalar();
    Scalar ref = Scalar();
    for (unsigned int idim = 0; idim < 3; ++idim) {
        dot += (base(idim) - vertex(idim)) * (point[idim] - vertex(idim));
        dist += pow<2>(point[idim] - vertex(idim));
        ref += pow<2>(base(idim) - vertex(idim));
    }
    return (!(dot > Scalar()) ||
            !(std::acos(dot / rt<2>(dist * ref)) < angle()) || !(dot < ref));
}
// -------------------------------------------------------------------------- //

// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Cone.
/// \return         0 if no error.
template <class Vector, typename Scalar>
int Cone<Vector, Scalar>::example() {
    // Initialize
    std::cout << "BEGIN = Cone::example()" << std::endl;
    std::cout << std::boolalpha << std::left;
    const unsigned int width = 40;
    std::array<double, 3> array({{4, 8, 15}});
    std::array<double, 3> array2({{16, 23, 42}});

    // Construction
    Cone<std::array<double, 3>, double> cone(array, array2, 0.16);

    // Lifecycle
    std::cout << std::endl;
    std::cout << std::setw(width * 2) << "Lifecycle : " << std::endl;
    std::cout << std::setw(width * 2) << "Cone<>() : " << Cone<>() << std::endl;
    std::cout << std::setw(width * 2)
              << "Cone<std::array<double, 3>, double>() : "
              << Cone<std::array<double, 3>, double>() << std::endl;
    std::cout << std::setw(width * 2)
              << "Cone<std::array<double, 3>, double>(array, array2, 0.16) : "
              << Cone<std::array<double, 3>, double>(array, array2, 0.16)
              << std::endl;
    std::cout << std::setw(width * 2)
              << "Cone<std::array<double, 3>, double>(cone) : "
              << Cone<std::array<double, 3>, double>(cone) << std::endl;

    // Data
    std::cout << std::endl;
    std::cout << std::setw(width) << "Data : " << std::endl;
    std::cout << std::setw(width)
              << "cone.vertex()[0] = 23 : " << (cone.vertex()[0] = 23)
              << std::endl;
    std::cout << std::setw(width) << "cone.vertex()[0] : " << cone.vertex()[0]
              << std::endl;
    std::cout << std::setw(width)
              << "cone.base()[0] = 23 : " << (cone.base()[0] = 23) << std::endl;
    std::cout << std::setw(width) << "cone.base()[0] : " << cone.base()[0]
              << std::endl;
    std::cout << std::setw(width)
              << "cone.angle() = 0.23 : " << (cone.angle() = 0.23) << std::endl;
    std::cout << std::setw(width) << "cone.angle() : " << cone.angle()
              << std::endl;

    // Position
    std::cout << std::endl;
    std::cout << std::setw(width) << "Position : " << std::endl;
    std::cout << std::setw(width) << "cone.direction(0) : " << cone.direction(0)
              << std::endl;
    std::cout << std::setw(width)
              << "cone.direction()[1] : " << cone.direction()[1] << std::endl;

    // Measures
    std::cout << std::endl;
    std::cout << std::setw(width) << "Measures : " << std::endl;
    std::cout << std::setw(width) << "cone.length() : " << cone.length()
              << std::endl;
    std::cout << std::setw(width) << "cone.radius() : " << cone.radius()
              << std::endl;
    std::cout << std::setw(width) << "cone.diameter() : " << cone.diameter()
              << std::endl;
    std::cout << std::setw(width) << "cone.circle() : " << cone.circle()
              << std::endl;
    std::cout << std::setw(width) << "cone.surface() : " << cone.surface()
              << std::endl;
    std::cout << std::setw(width) << "cone.volume() : " << cone.volume()
              << std::endl;

    // Collision
    std::cout << std::endl;
    std::cout << std::setw(width) << "Collision : " << std::endl;
    std::cout << std::setw(width) << "cone.inside(array) : " << cone.inside(array)
              << std::endl;
    std::cout << std::setw(width)
              << "cone.outside(array) : " << cone.outside(array) << std::endl;

    // Finalize
    std::cout << std::noboolalpha << std::right << std::endl;
    std::cout << "END = Cone::example()" << std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
#endif // CONE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
