/* ******************************* HYPERCUBE ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          HyperCube
// DESCRIPTION :    N-dimensional cube
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           hypercube.h
/// \brief          N-dimensional cube
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef HYPERCUBE_H_INCLUDED
#define HYPERCUBE_H_INCLUDED
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
#include "abstracthypercube.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// N-dimensional cube
/// \brief          N-dimensional cube.
/// \details        Implementation of a basic hypercube of arbitrary dimension.
/// \tparam         Dimension Number of space dimension.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type.
template <unsigned int Dimension = 3, class Vector = std::array<double, Dimension>, typename Scalar = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Vector>()[0])>::type>::type>
class HyperCube final
: public AbstractHyperCube<HyperCube<Dimension, Vector, Scalar>, Dimension, Vector, Scalar>
, public AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar>
{
    // Setup 
    public: using AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar>::operator=; 
    
    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        template <class... Misc> explicit inline HyperCube(Misc&&... misc); 
    //@}
    
    // Data 
    /// \name           Data 
    //@{
    public: 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template position(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template position(Misc&&... misc) const; 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template extent(Misc&&... misc); 
        template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar> >().template data<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type> inline Template extent(Misc&&... misc) const; 
    //@}

    // Predefined
    /// \name           Predefined
    //@{
    public:
        static constexpr HyperCube<Dimension, Vector, Scalar> unit();
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
inline HyperCube<Dimension, Vector, Scalar>::HyperCube(Misc&&... misc)
: AbstractHyperCube<HyperCube<Dimension, Vector, Scalar>, Dimension, Vector, Scalar>()
, AbstractSubstance<HyperCube<Dimension, Vector, Scalar>, Vector, Scalar>(std::forward<Misc>(misc)...)
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
inline Template HyperCube<Dimension, Vector, Scalar>::position(Misc&&... misc)
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
inline Template HyperCube<Dimension, Vector, Scalar>::position(Misc&&... misc) const
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
inline Template HyperCube<Dimension, Vector, Scalar>::extent(Misc&&... misc)
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
inline Template HyperCube<Dimension, Vector, Scalar>::extent(Misc&&... misc) const
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Unit hypercube
/// \brief          Unit hypercube.
/// \details        Creates an hypercube with a position of zero and with an 
///                 extent of one.
/// \return         Copy of a unit hypercube.
template <unsigned int Dimension, class Vector, typename Scalar>
constexpr HyperCube<Dimension, Vector, Scalar> HyperCube<Dimension, Vector, Scalar>::unit()
{
    return HyperCube<Dimension, Vector, Scalar>(Vector(), Scalar(1));
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of HyperCube.
/// \return         0 if no error.
template <unsigned int Dimension, class Vector, typename Scalar>
int HyperCube<Dimension, Vector, Scalar>::example()
{
    // Initialize
    std::cout<<"BEGIN = HyperCube::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<double, 3> array({{4, 8, 15}});
    std::mt19937 e;
    std::uniform_real_distribution<double> d(0, 100);

    // Construction
    HyperCube<3, std::array<double, 3>, double> hypercube;
  
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"HyperCube<>() : "                                                                   <<HyperCube<>()<<std::endl;
    std::cout<<std::setw(width*2)<<"HyperCube<3, std::array<double, 3>, double>() : "                                   <<HyperCube<3, std::array<double, 3>, double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"HyperCube<3, std::array<double, 3>, double>(array, 16) : "                          <<HyperCube<3, std::array<double, 3>, double>(array, 16)<<std::endl;
    std::cout<<std::setw(width*2)<<"HyperCube<3, std::array<double, 3>, double>(hypercube) : "                          <<HyperCube<3, std::array<double, 3>, double>(hypercube)<<std::endl;

    // Data
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Data : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.position()[0] = 23 : "              <<(hypercube.position()[0] = 23)<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.position()[0] : "                   <<hypercube.position()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.extent() = 42 : "                   <<(hypercube.extent() = 42)<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.extent() : "                        <<hypercube.extent()<<std::endl;

    // Position
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Position : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.center(0) : "                       <<hypercube.center(0)<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.center()[1] : "                     <<hypercube.center()[1]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.minimum(0) : "                      <<hypercube.minimum(0)<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.minimum()[1] : "                    <<hypercube.minimum()[1]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.maximum(0) : "                      <<hypercube.maximum(0)<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.maximum()[1] : "                    <<hypercube.maximum()[1]<<std::endl;

    // Measures
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Measures : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.length() : "                        <<hypercube.length()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.volume() : "                        <<hypercube.volume()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.surface() : "                       <<hypercube.surface()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.surface<3>() : "                    <<hypercube.surface<3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.surface<2>() : "                    <<hypercube.surface<2>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.surface<1>() : "                    <<hypercube.surface<1>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.surface<0>() : "                    <<hypercube.surface<0>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.area() : "                          <<hypercube.area()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.area<3>() : "                       <<hypercube.area<3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.area<2>() : "                       <<hypercube.area<2>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.area<1>() : "                       <<hypercube.area<1>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.area<0>() : "                       <<hypercube.area<0>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.diagonal() : "                      <<hypercube.diagonal()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.diagonal<3>() : "                   <<hypercube.diagonal<3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.diagonal<2>() : "                   <<hypercube.diagonal<2>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.diagonal<1>() : "                   <<hypercube.diagonal<1>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.diagonal<0>() : "                   <<hypercube.diagonal<0>()<<std::endl;

    // Distribution
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Distribution : "                              <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random()[0] : "                     <<hypercube.random()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<3>()[0] : "                  <<hypercube.random<3>()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<2>()[0] : "                  <<hypercube.random<2>()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<1>()[0] : "                  <<hypercube.random<1>()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<0>()[0] : "                  <<hypercube.random<0>()[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<3>(e, d)[0] : "              <<hypercube.random<3>(e, d)[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<2>(e, d)[0] : "              <<hypercube.random<2>(e, d)[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<1>(e, d)[0] : "              <<hypercube.random<1>(e, d)[0]<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.random<0>(e, d)[0] : "              <<hypercube.random<0>(e, d)[0]<<std::endl;

    // Collision
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Collision : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.inside(array) : "                   <<hypercube.inside(array)<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.outside(array) : "                  <<hypercube.outside(array)<<std::endl;

    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Predefined : "                                <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.unit() : "                          <<hypercube.unit()<<std::endl;
    
    // Constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Constants : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"hypercube.dimension() : "                     <<hypercube.dimension()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.elements() : "                      <<hypercube.elements()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.elements<4>() : "                   <<hypercube.elements<4>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.elements<3>() : "                   <<hypercube.elements<3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.elements<2>() : "                   <<hypercube.elements<2>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.elements<1>() : "                   <<hypercube.elements<1>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.elements<0>() : "                   <<hypercube.elements<0>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements() : "                   <<hypercube.subelements()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<4, 3>() : "             <<hypercube.subelements<4, 3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<3, 4>() : "             <<hypercube.subelements<3, 4>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<0>() : "                <<hypercube.subelements<0>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<0, 3>() : "             <<hypercube.subelements<0, 3>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<1, 2>() : "             <<hypercube.subelements<1, 2>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<2, 1>() : "             <<hypercube.subelements<2, 1>()<<std::endl;
    std::cout<<std::setw(width)<<"hypercube.subelements<0, 0>() : "             <<hypercube.subelements<0, 0>()<<std::endl;
   
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = HyperCube::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // HYPERCUBE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
