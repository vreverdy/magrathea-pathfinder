/* ******************************* ABOUTOBJECT ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          AboutObject 
// DESCRIPTION :    Basic about object with information on something 
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           aboutobject.h 
/// \brief          Basic about object with information on something 
/// \author         Vincent Reverdy (vince.rev@gmail.com) 
/// \date           2012-2013 
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABOUTOBJECT_H_INCLUDED
#define ABOUTOBJECT_H_INCLUDED
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
#include "abstractaboutobject.h"
// Misc
namespace magrathea {
//--------------------------------------------------------------------------- //



//---------------------------------- CLASS ---------------------------------- //
// Basic about object with information on something 
/// \brief          Basic about object with information on something. 
/// \details        This class is the direct derivation of 
///                 AbstractAboutObject. It provides the most basic and 
///                 generic contents object without adding new functionalities 
///                 to the abstract class. It can be used in most cases as a 
///                 generic container of groups of information on something. 
/// \tparam         Types Variadic list of components types.
template <class... Types>
class AboutObject final
: public AbstractAboutObject<AboutObject<Types...>, Types...>
{
    // Setup 
    public: using AbstractAboutObject<AboutObject<Types...>, Types...>::operator=; 

    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        template <class... Misc> explicit inline AboutObject(Misc&&... misc); 
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
template <class... Types>
template <class... Misc>
inline AboutObject<Types...>::AboutObject(Misc&&... misc)
: AbstractAboutObject<AboutObject<Types...>, Types...>(std::forward<Misc>(misc)...)
{
    ;
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of AboutObject. 
/// \return         0 if no error. 
template <class... Types>
int AboutObject<Types...>::example()
{
    // Initialize
    std::cout<<"BEGIN = Contents::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<double, 3> arr({{42., 42., 42.}});
    std::tuple<std::array<double, 3> > dat(arr);
    std::ostringstream stream;

    // Construction
    AboutObject<int> i(4);
    AboutObject<int> j(8);
    AboutObject<double> d(15.16);
    AboutObject<std::array<double, 3> > a(std::array<double, 3>({{23, 42, 4}}));
    AboutObject<std::string> s("The answer is 42");

    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"AboutObject<int>() : "                                                              <<AboutObject<int>()<<std::endl;
    std::cout<<std::setw(width*2)<<"AboutObject<int>(d) : "                                                             <<AboutObject<int>(d)<<std::endl;
    std::cout<<std::setw(width*2)<<"AboutObject<double>(i) : "                                                          <<AboutObject<double>(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"AboutObject<int>(std::make_tuple(42)) : "                                           <<AboutObject<int>(std::make_tuple(42))<<std::endl;
    std::cout<<std::setw(width*2)<<"AboutObject<int>(42) : "                                                            <<AboutObject<int>(42)<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"i = j : "                                     <<(i = j)<<std::endl;
    std::cout<<std::setw(width)<<"i = d : "                                     <<(i = d)<<std::endl;
    std::cout<<std::setw(width)<<"i = std::make_tuple(42) : "                   <<(i = std::make_tuple(42))<<std::endl;
    std::cout<<std::setw(width)<<"i == d : "                                    <<(i == d)<<std::endl;
    std::cout<<std::setw(width)<<"i != d : "                                    <<(i != d)<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"i.assign() : "                                <<i.assign()<<std::endl;
    std::cout<<std::setw(width)<<"i.assign(j) : "                               <<i.assign(j)<<std::endl;
    std::cout<<std::setw(width)<<"i.assign(d) : "                               <<i.assign(d)<<std::endl;
    std::cout<<std::setw(width)<<"d.assign(i) : "                               <<d.assign(i)<<std::endl;
    std::cout<<std::setw(width)<<"i.assign(std::make_tuple(42)) : "             <<i.assign(std::make_tuple(42))<<std::endl;
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
    std::cout<<"END = AboutObject::example()"<<std::endl;
    return 0;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABOUTOBJECT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

