/* ****************************** HYPERSPHERE ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          HyperSphere
// DESCRIPTION :    N-dimensional sphere
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           hypersphere.h
/// \brief          N-dimensional sphere
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef HYPERSPHERE_H_INCLUDED
#define HYPERSPHERE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <array>
#include <utility>
// Include libs
// Include project
#include "abstractsubstance.h"
#include "abstracthypersphere.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// N-dimensional sphere
/// \brief          N-dimensional sphere.
/// \details        Implementation of a basic hypersphere of arbitrary dimension.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
template <unsigned int Dimension = 3, class Vector = std::array<double, Dimension>, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type>
class HyperSphere final
: public AbstractHyperSphere<HyperSphere<Dimension, Vector, Scalar>, Dimension, Vector, Scalar>
, public AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar>
{
    // Setup 
    public: using AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar>::operator=; 
    
    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        template <class... Misc> explicit inline HyperSphere(Misc&&... misc); 
    //@}
    
    // Data 
    /// \name           Data 
    //@{
    public: 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template position(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template position(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template extent(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template extent(Misc&&... misc) const; 
    //@}

    // Predefined
    /// \name           Predefined
    //@{
    public:
        static constexpr HyperSphere<Dimension, Vector, Scalar> unit();
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
template <unsigned int Dimension, class Vector, typename Scalar>
template <class... Misc>
inline HyperSphere<Dimension, Vector, Scalar>::HyperSphere(Misc&&... misc)
: AbstractHyperSphere<HyperSphere<Dimension, Vector, Scalar>, Dimension, Vector, Scalar>()
, AbstractSubstance<HyperSphere<Dimension, Vector, Scalar>, Vector, Scalar>(std::forward<Misc>(misc)...)
{
    ;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- DATA ---------------------------------- //
// Access to the position data
/// \brief          Access to the position data. 
/// \details        Provides an access to the position data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template HyperSphere<Dimension, Vector, Scalar>::position(Misc&&... misc)
{
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the position data
/// \brief          Immutable access to the position data. 
/// \details        Provides an immutable access to the position data by 
///                 forwarding parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template HyperSphere<Dimension, Vector, Scalar>::position(Misc&&... misc) const
{
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Access to the extent data
/// \brief          Access to the extent data. 
/// \details        Provides an access to the extent data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template HyperSphere<Dimension, Vector, Scalar>::extent(Misc&&... misc)
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the extent data
/// \brief          Immutable access to the extent data. 
/// \details        Provides an immutable access to the extent data by 
///                 forwarding parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result. 
template <unsigned int Dimension, class Vector, typename Scalar>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template HyperSphere<Dimension, Vector, Scalar>::extent(Misc&&... misc) const
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Unit hypersphere
/// \brief          Unit hypersphere.
/// \details        Creates an hypersphere with a position of zero and with an 
///                 extent of one.
/// \return         Copy of a unit hypersphere.
template <unsigned int Dimension, class Vector, typename Scalar>
constexpr HyperSphere<Dimension, Vector, Scalar> HyperSphere<Dimension, Vector, Scalar>::unit()
{
    return HyperSphere<Dimension, Vector, Scalar>(Vector(), Scalar(1));
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of HyperSphere.
/// \return         0 if no error.
template <unsigned int Dimension, class Vector, typename Scalar>
int HyperSphere<Dimension, Vector, Scalar>::example()
{
    // Initialize
    std::cout<<"BEGIN = HyperSphere::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<std::array<double, 3>, 100> points = std::array<std::array<double, 3>, 100>();
    std::array<double, 3> array({{4, 8, 15}});
    std::mt19937 e;
    std::uniform_real_distribution<double> d(0, 100);

    // Construction
    HyperSphere<3, std::array<double, 3>, double> hypersphere;
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"HyperSphere<>() : "                                                                 <<HyperSphere<>()<<std::endl;
    std::cout<<std::setw(width*2)<<"HyperSphere<3, std::array<double, 3>, double>() : "                                 <<HyperSphere<3, std::array<double, 3>, double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"HyperSphere<3, std::array<double, 3>, double>(array, 16) : "                        <<HyperSphere<3, std::array<double, 3>, double>(array, 16)<<std::endl;
    std::cout<<std::setw(width*2)<<"HyperSphere<3, std::array<double, 3>, double>(hypersphere) : "                      <<HyperSphere<3, std::array<double, 3>, double>(hypersphere)<<std::endl;

    // Data
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Data : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.position()[0] = 23 : "            <<(hypersphere.position()[0] = 23)<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.position()[0] : "                 <<hypersphere.position()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.extent() = 42 : "                 <<(hypersphere.extent() = 42)<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.extent() : "                      <<hypersphere.extent()<<std::endl;
    
    // Position
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Position : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.center(0) : "                     <<hypersphere.center(0)<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.center()[1] : "                   <<hypersphere.center()[1]<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.minimum(0) : "                    <<hypersphere.minimum(0)<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.minimum()[1] : "                  <<hypersphere.minimum()[1]<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.maximum(0) : "                    <<hypersphere.maximum(0)<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.maximum()[1] : "                  <<hypersphere.maximum()[1]<<std::endl;

    // Measures
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Measures : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.radius() : "                      <<hypersphere.radius()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.diameter() : "                    <<hypersphere.diameter()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.volume() : "                      <<hypersphere.volume()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.surface() : "                     <<hypersphere.surface()<<std::endl;

    // Distribution
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Distribution : "                                                                    <<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.random()[0] : "                                                         <<hypersphere.random()[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.random<3>()[0] : "                                                      <<hypersphere.random<3>()[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.random<2>()[0] : "                                                      <<hypersphere.random<2>()[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.random(e, d)[0] : "                                                     <<hypersphere.random(e, d)[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.random<3>(e, d)[0] : "                                                  <<hypersphere.random<3>(e, d)[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.random<2>(e, d)[0] : "                                                  <<hypersphere.random<2>(e, d)[0]<<std::endl;
    std::cout<<std::setw(width*2)<<"hypersphere.uniform<2>(points.begin(), points.end()).first : "                      <<hypersphere.uniform<2>(points.begin(), points.end()).first<<std::endl;

    // Collision
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Collision : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.inside(array) : "                 <<hypersphere.inside(array)<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.outside(array) : "                <<hypersphere.outside(array)<<std::endl;

    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Predefined : "                                <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.unit() : "                        <<hypersphere.unit()<<std::endl;
    
    // Constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Constants : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.dimension() : "                   <<hypersphere.dimension()<<std::endl;

    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Helpers : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn() : "                          <<hypersphere.sn()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn<5>() : "                       <<hypersphere.sn<5>()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn<4>() : "                       <<hypersphere.sn<4>()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn<3>() : "                       <<hypersphere.sn<3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn<2>() : "                       <<hypersphere.sn<2>()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn<1>() : "                       <<hypersphere.sn<1>()<<std::endl;
    std::cout<<std::setw(width)<<"hypersphere.sn<0>() : "                       <<hypersphere.sn<0>()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = HyperSphere::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // HYPERSPHERE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

