/* ********************************** STEP ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          Step 
// DESCRIPTION :    Basic implementation of an evolution step 
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           step.h 
/// \brief          Basic implementation of an evolution step 
/// \author         Vincent Reverdy (vince.rev@gmail.com) 
/// \date           2012-2013 
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef STEP_H_INCLUDED
#define STEP_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



//------------------------------- PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <tuple>
#include <array>
#include <ratio>
#include <string>
#include <sstream>
// Include libs
// Include project
#include "abstractstep.h"
// Misc
namespace magrathea {
//--------------------------------------------------------------------------- //



//---------------------------------- CLASS ---------------------------------- //
// Basic implementation of an evolution step 
/// \brief          Basic implementation of an evolution step. 
/// \details        This class is the direct derivation of AbstractStep. It 
///                 provides the most basic and generic contents object 
///                 without adding new functionalities to the abstract class. 
///                 It can be used in most cases as a generic container for 
///                 evolution data. 
/// \tparam         Scalar Scalar type of id.
/// \tparam         Array Array type of core quantities.
/// \tparam         Tuple Tuple type of extra quantities.
template <class Scalar = unsigned int, class Array = std::array<double, 0>, class Tuple = std::tuple<> >
class Step final
: public AbstractStep<Step<Scalar, Array, Tuple>, Scalar, Array, Tuple>
{
    // Setup 
    public: using AbstractStep<Step<Scalar, Array, Tuple>, Scalar, Array, Tuple>::operator=; 

    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        template <class... Misc> explicit inline Step(Misc&&... misc); 
    //@}

    // Test 
    /// \name           Test 
    //@{
    public: 
        static int example(); 
    //@}
};
//--------------------------------------------------------------------------- //



//-------------------------------- LIFECYCLE -------------------------------- //
// Explicit generic constructor 
/// \brief          Explicit generic constructor. 
/// \details        Provides a generic interface to all constructors of the 
///                 base class. 
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
template <class Scalar, class Array, class Tuple>
template <class... Misc>
inline Step<Scalar, Array, Tuple>::Step(Misc&&... misc)
: AbstractStep<Step<Scalar, Array, Tuple>, Scalar, Array, Tuple>(std::forward<Misc>(misc)...)
{
    ;
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of Step. 
/// \return         0 if no error. 
template <class Scalar, class Array, class Tuple>
int Step<Scalar, Array, Tuple>::example()
{
    // Initialize
    std::cout<<"BEGIN = Contents::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<double, 3> arr({{42., 42., 42.}});
    std::tuple<std::array<double, 3>, std::array<double, 0>, std::tuple<> > dat(arr, std::array<double, 0>(), std::make_tuple());
    std::array<std::array<double, 3>, 3> all({{arr, arr, arr}});
    std::ostringstream stream;

    // Construction
    Step<int, std::array<double, 0>, std::tuple<> > i(4);
    Step<int, std::array<double, 0>, std::tuple<> > j(8);
    Step<double, std::array<double, 0>, std::tuple<> > d(15.16);
    Step<std::array<double, 3>, std::array<double, 0>, std::tuple<> > a(std::array<double, 3>({{23, 42, 4}}));
    Step<std::string, std::array<double, 0>, std::tuple<> > s("The answer is 42");
    Step<std::array<std::array<double, 3>, 3>, std::array<std::array<double, 3>, 3>, std::array<std::array<double, 3>, 3> > aaa(all, all, all);

    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"Step<int, std::array<double, 0>, std::tuple<> >() : "                               <<Step<int, std::array<double, 0>, std::tuple<> >()<<std::endl;
    std::cout<<std::setw(width*2)<<"Step<int, std::array<double, 0>, std::tuple<> >(d) : "                              <<Step<int, std::array<double, 0>, std::tuple<> >(d)<<std::endl;
    std::cout<<std::setw(width*2)<<"Step<double, std::array<double, 0>, std::tuple<> >(i) : "                           <<Step<double, std::array<double, 0>, std::tuple<> >(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"Step<int, std::array<double, 0>, std::tuple<> >(42) : "                             <<Step<int, std::array<double, 0>, std::tuple<> >(42)<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"i = j : "                                     <<(i = j)<<std::endl;
    std::cout<<std::setw(width)<<"i = d : "                                     <<(i = d)<<std::endl;
    std::cout<<std::setw(width)<<"i == d : "                                    <<(i == d)<<std::endl;
    std::cout<<std::setw(width)<<"i != d : "                                    <<(i != d)<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"i.assign() : "                                <<i.assign()<<std::endl;
    std::cout<<std::setw(width)<<"i.assign(j) : "                               <<i.assign(j)<<std::endl;
    std::cout<<std::setw(width)<<"i.assign(d) : "                               <<i.assign(d)<<std::endl;
    std::cout<<std::setw(width)<<"d.assign(i) : "                               <<d.assign(i)<<std::endl;
    std::cout<<std::setw(width)<<"i.assign(42) : "                              <<i.assign(42)<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"i.nullify() : "                               <<i.nullify()<<std::endl;
    std::cout<<std::setw(width)<<"i.copy() : "                                  <<i.copy()<<std::endl;
    std::cout<<std::setw(width)<<"i.cast() : "                                  <<i.cast()<<std::endl;

    // Data
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Data : "                                                                            <<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<0>(a.data())[0] = 0 : "                                                    <<(std::get<0>(a.data())[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<0>(a.data())[0] : "                                                        <<(std::get<0>(a.data())[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0>()[0] = 0 : "                                                              <<(a.data<0>()[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0>()[0] : "                                                                  <<(a.data<0>()[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0>(std::get<0>(dat)) : "                                                     <<(a.data<0>(std::get<0>(dat)))<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0, 0>() = 0 : "                                                              <<(a.data<0, 0>() = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0, 0>() : "                                                                  <<(a.data<0, 0>())<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0, 0>(std::get<0>(std::get<0>(dat))) : "                                     <<(a.data<0, 0>(std::get<0>(std::get<0>(dat))))<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0>(0) = 0 : "                                                                <<(a.data<0>(0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0>(0) : "                                                                    <<(a.data<0>(0))<<std::endl;
    std::cout<<std::setw(width*2)<<"a.data<0>(0, std::get<0>(dat)[0]) : "                                               <<(a.data<0>(0, std::get<0>(dat)[0]))<<std::endl;

    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"std::get<0>(a.get())[0] : "                   <<std::get<0>(a.get())[0]<<std::endl;
    std::cout<<std::setw(width)<<"a.get<0>()[0] : "                             <<a.get<0>()[0]<<std::endl;
    std::cout<<std::setw(width)<<"a.get<0, 0>() : "                             <<a.get<0, 0>()<<std::endl;
    std::cout<<std::setw(width)<<"a.get<0>(0) : "                               <<a.get<0>(0)<<std::endl;

    // Setters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Setters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"a.set(dat) : "                                <<a.set(dat)<<std::endl;
    std::cout<<std::setw(width)<<"a.set<0>(arr) : "                             <<a.set<0>(arr)<<std::endl;
    std::cout<<std::setw(width)<<"a.set<0, 0>(15) : "                           <<a.set<0, 0>(15)<<std::endl;
    std::cout<<std::setw(width)<<"a.set<0>(0, 16) : "                           <<a.set<0>(0, 16)<<std::endl;

    // Id
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Id : "                                                                              <<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id()[0][0] = 0 : "                                                              <<(aaa.id()[0][0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id()[0][0] : "                                                                  <<(aaa.id()[0][0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(all) : "                                                                     <<(aaa.id(all))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0>()[0] = 0 : "                                                              <<(aaa.id<0>()[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0>()[0] : "                                                                  <<(aaa.id<0>()[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0>(std::get<0>(all)) : "                                                     <<(aaa.id<0>(std::get<0>(all)))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0, 0>() = 0 : "                                                              <<(aaa.id<0, 0>() = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0, 0>() : "                                                                  <<(aaa.id<0, 0>())<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0, 0>(std::get<0>(std::get<0>(all))) : "                                     <<(aaa.id<0, 0>(std::get<0>(std::get<0>(all))))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(0)[0] = 0 : "                                                                <<(aaa.id(0)[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(0)[0] : "                                                                    <<(aaa.id(0)[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(0, all[0]) : "                                                               <<(aaa.id(0, all[0]))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(0, 0) = 0 : "                                                                <<(aaa.id(0, 0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(0, 0) : "                                                                    <<(aaa.id(0, 0))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id(0, 0, all[0][0]) : "                                                         <<(aaa.id(0, 0, all[0][0]))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0>(0) = 0 : "                                                                <<(aaa.id<0>(0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0>(0) : "                                                                    <<(aaa.id<0>(0))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.id<0>(0, std::get<0>(all)[0]) : "                                               <<(aaa.id<0>(0, std::get<0>(all)[0]))<<std::endl;

    // Core
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Core : "                                                                            <<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core()[0][0] = 0 : "                                                            <<(aaa.core()[0][0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core()[0][0] : "                                                                <<(aaa.core()[0][0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(all) : "                                                                   <<(aaa.core(all))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0>()[0] = 0 : "                                                            <<(aaa.core<0>()[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0>()[0] : "                                                                <<(aaa.core<0>()[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0>(std::get<0>(all)) : "                                                   <<(aaa.core<0>(std::get<0>(all)))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0, 0>() = 0 : "                                                            <<(aaa.core<0, 0>() = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0, 0>() : "                                                                <<(aaa.core<0, 0>())<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0, 0>(std::get<0>(std::get<0>(all))) : "                                   <<(aaa.core<0, 0>(std::get<0>(std::get<0>(all))))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(0)[0] = 0 : "                                                              <<(aaa.core(0)[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(0)[0] : "                                                                  <<(aaa.core(0)[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(0, all[0]) : "                                                             <<(aaa.core(0, all[0]))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(0, 0) = 0 : "                                                              <<(aaa.core(0, 0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(0, 0) : "                                                                  <<(aaa.core(0, 0))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core(0, 0, all[0][0]) : "                                                       <<(aaa.core(0, 0, all[0][0]))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0>(0) = 0 : "                                                              <<(aaa.core<0>(0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0>(0) : "                                                                  <<(aaa.core<0>(0))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.core<0>(0, std::get<0>(all)[0]) : "                                             <<(aaa.core<0>(0, std::get<0>(all)[0]))<<std::endl;

    // Extra
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Extra : "                                                                           <<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra()[0][0] = 0 : "                                                           <<(aaa.extra()[0][0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra()[0][0] : "                                                               <<(aaa.extra()[0][0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(all) : "                                                                  <<(aaa.extra(all))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0>()[0] = 0 : "                                                           <<(aaa.extra<0>()[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0>()[0] : "                                                               <<(aaa.extra<0>()[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0>(std::get<0>(all)) : "                                                  <<(aaa.extra<0>(std::get<0>(all)))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0, 0>() = 0 : "                                                           <<(aaa.extra<0, 0>() = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0, 0>() : "                                                               <<(aaa.extra<0, 0>())<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0, 0>(std::get<0>(std::get<0>(all))) : "                                  <<(aaa.extra<0, 0>(std::get<0>(std::get<0>(all))))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(0)[0] = 0 : "                                                             <<(aaa.extra(0)[0] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(0)[0] : "                                                                 <<(aaa.extra(0)[0])<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(0, all[0]) : "                                                            <<(aaa.extra(0, all[0]))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(0, 0) = 0 : "                                                             <<(aaa.extra(0, 0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(0, 0) : "                                                                 <<(aaa.extra(0, 0))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra(0, 0, all[0][0]) : "                                                      <<(aaa.extra(0, 0, all[0][0]))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0>(0) = 0 : "                                                             <<(aaa.extra<0>(0) = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0>(0) : "                                                                 <<(aaa.extra<0>(0))<<std::endl;
    std::cout<<std::setw(width*2)<<"aaa.extra<0>(0, std::get<0>(all)[0]) : "                                            <<(aaa.extra<0>(0, std::get<0>(all)[0]))<<std::endl;

    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, i) : "                  <<i<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, d) : "                  <<d<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, a) : "                  <<a<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, s) : "                  <<s<<std::endl;

    // Types
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Types : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"typeid(a.type()).name() : "                   <<typeid(a.type()).name()<<std::endl;
    std::cout<<std::setw(width)<<"typeid(a.type<0>()).name() : "                <<typeid(a.type<0>()).name()<<std::endl;
    std::cout<<std::setw(width)<<"typeid(a.type<0, 0>()).name() : "             <<typeid(a.type<0, 0>()).name()<<std::endl;
    std::cout<<std::setw(width)<<"typeid(a.type<0>(0)).name() : "               <<typeid(a.type<0>(0)).name()<<std::endl;

    // Properties
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Properties : "                                <<std::endl;
    std::cout<<std::setw(width)<<"i.types() : "                                 <<i.types()<<std::endl;

    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Helpers : "                                                                         <<std::endl;
    std::cout<<std::setw(width*2)<<"i.monomial(42) : "                                                                  <<i.monomial(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"i.monomial<2>(42) : "                                                               <<i.monomial<2>(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"i.monomial<2, std::ratio<2> >(42) : "                                               <<i.monomial<2, std::ratio<2> >(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"i.monomial<std::ratio<2> >(42) : "                                                  <<i.monomial<std::ratio<2> >(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"i.monomial<std::ratio<2>, 2>(42) : "                                                <<i.monomial<std::ratio<2>, 2>(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"typeid(i.transmute<std::tuple<int, int> >()).name() : "                             <<typeid(i.transmute<std::tuple<int, int> >()).name()<<std::endl;
    std::cout<<std::setw(width*2)<<"typeid(i.transmute<std::tuple<int, int> >(std::tuple<int, int>())).name() : "       <<typeid(i.transmute<std::tuple<int, int> >(std::tuple<int, int>())).name()<<std::endl;
    std::cout<<std::setw(width*2)<<"typeid(i.transmute<std::tuple<int, int> >(42)).name() : "                           <<typeid(i.transmute<std::tuple<int, int> >(42)).name()<<std::endl;
    std::cout<<std::setw(width*2)<<"typeid(i.transmute<std::tuple<int, int> >(42, 42)).name() : "                       <<typeid(i.transmute<std::tuple<int, int> >(42, 42)).name()<<std::endl;
    std::cout<<std::setw(width*2)<<"i.printable<std::string>() : "                                                      <<i.printable<std::string>()<<std::endl;
    std::cout<<std::setw(width*2)<<"i.printable<std::tuple<> >() : "                                                    <<i.printable<std::tuple<> >()<<std::endl;
    std::cout<<std::setw(width*2)<<"i.print(stream, 42) : "                                                             <<i.print(stream, 42)<<std::endl;
    std::cout<<std::setw(width*2)<<"i.print(stream, std::array<int, 6>({{4, 8, 15, 16, 23, 42}})) : "                   <<i.print(stream, std::array<int, 6>({{4, 8, 15, 16, 23, 42}}))<<std::endl;
    std::cout<<std::setw(width*2)<<"i.print(stream, std::make_tuple(4, 8, 15, 16, 23, 42)) : "                          <<i.print(stream, std::make_tuple(4, 8, 15, 16, 23, 42))<<std::endl;
    std::cout<<std::setw(width*2)<<"i.print(stream) : "                                                                 <<i.print(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"i.print(stream, std::make_tuple()) : "                                              <<i.print(stream, std::make_tuple())<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Step::example()"<<std::endl;
    return 0;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // STEP_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

