/* ********************************* NARRAY ********************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          NArray
// DESCRIPTION :    Basic n-dimensional mathematical array
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           narray.h
/// \brief          Basic n-dimensional mathematical array
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef NARRAY_H_INCLUDED
#define NARRAY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <initializer_list>
#include <algorithm>
#include <random>
#include <cmath>
// Include libs
// Include project
#include "abstractnarray.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Basic n-dimensional mathematical array
/// \brief          Basic n-dimensional mathematical array.
/// \details        This class is the direct derivation of AbstractNArray. It 
///                 provides the most basic n-dimensional mathematical array
///                 without adding new functionalities to the abstract class.
/// \tparam         Type Data type.
/// \tparam         Size Number of elements.
template <typename Type = double, unsigned int Size = 1>
class NArray final
: public AbstractNArray<unsigned int, Size, NArray, Type, Size>
{
    // Setup
    public: using AbstractNArray<unsigned int, Size, NArray, Type, Size>::operator=;
    
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        inline NArray();
        template <typename FundamentalType = Type, class = typename std::enable_if<std::is_fundamental<FundamentalType>::value>::type> inline NArray(const NArray<FundamentalType, Size>& source);
        template <typename OtherType = Type, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline NArray(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc, class = typename std::enable_if<sizeof...(Misc) != 0>::type> explicit inline NArray(const Misc&... misc);
    //@}

    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}

    // Data members
    /// \name           Data members
    //@{
    protected:
        Type _data[Size];                                                       ///< Data contents. 
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit empty constructor
/// \brief          Implicit empty constructor.
/// \details        Does nothing.
template <typename Type, unsigned int Size>
inline NArray<Type, Size>::NArray()
: AbstractNArray<unsigned int, Size, NArray, Type, Size>()
{
    ;
}

// Implicit conversion constructor
/// \brief          Implicit conversion constructor.
/// \details        Provides an implicit conversion from a fundamental type
///                 contents.
/// \tparam         FundamentalType (Fundamental data type.)
/// \param[in]      source Source of the copy.   
template <typename Type, unsigned int Size>
template <typename FundamentalType, class> 
inline NArray<Type, Size>::NArray(const NArray<FundamentalType, Size>& source)
: AbstractNArray<unsigned int, Size, NArray, Type, Size>(source)
{
    ;
}

// Implicit initializer list constructor
/// \brief          Implicit initializer list constructor.
/// \details        Provides an implicit conversion from an initializer list.
/// \tparam         OtherType (Other data type.)
/// \tparam         Misc (Miscellaneous types.)  
/// \param[in]      source Source of the copy. 
/// \param[in]      misc Miscellaneous arguments.  
template <typename Type, unsigned int Size>
template <typename OtherType, class... Misc, class>
inline NArray<Type, Size>::NArray(const std::initializer_list<OtherType>& source, const Misc&... misc)
: AbstractNArray<unsigned int, Size, NArray, Type, Size>(source, misc...)
{
    ;
} 

// Explicit generic constructor
/// \brief          Explicit generic constructor.
/// \details        Provides a generic interface to all constructors of 
///                 the base class. Before calling the associated constructor of 
///                 the base class, the contents is initialized.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.  
template <typename Type, unsigned int Size>
template <class... Misc, class> 
inline NArray<Type, Size>::NArray(const Misc&... misc)
: AbstractNArray<unsigned int, Size, NArray, Type, Size>(misc...)
{
    ;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of NArray.
/// \return         0 if no error.
template <typename Type, unsigned int Size>
int NArray<Type, Size>::example()
{
    // Initialize
    std::cout<<"BEGIN = NArray::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    std::uniform_int_distribution<int> distribution(0, 100);
    std::mt19937 engine;
    
    // Construction
    const NArray<int, 6> cnumbers({4, 8, 15, 16, 23, 42});
    NArray<int, 6> inumbers({4, 8, 15, 16, 23, 42});
    NArray<double, 6> dnumbers({42, 23, 16, 15, 8, 4});
    NArray<bool, 6> mask({true, true, false, true, false, true});
    NArray<int, 6> ten(10);
    NArray<double, 4> small(42);
    NArray<double, 9> big({42, 23, 16, 15, 8, 4, 2, 1, 0});
    NArray<unsigned int, 6> indexes({0, 1, 2, 3, 4, 5});
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"NArray<>() : "                                <<NArray<>()<<std::endl;
    std::cout<<std::setw(width)<<"NArray<double, 6>() : "                       <<NArray<double, 6>()<<std::endl;
    std::cout<<std::setw(width)<<"NArray<double, 6>(dnumbers) : "               <<NArray<double, 6>(dnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"NArray<double, 6>(cnumbers) : "               <<NArray<double, 6>(cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"NArray<double, 6>({42}) : "                   <<NArray<double, 6>({42})<<std::endl;
    std::cout<<std::setw(width)<<"NArray<int, 4>({0, 1, 2, 3, 4}) : "           <<NArray<int, 4>({0, 1, 2, 3, 4})<<std::endl;
    std::cout<<std::setw(width)<<"NArray<bool, 4>({0, 1, 0, 1}) : "             <<NArray<bool, 4>({0, 1, 0, 1})<<std::endl;
    std::cout<<std::setw(width)<<"NArray<double, 6>(42) : "                     <<NArray<double, 6>(42)<<std::endl;
    std::cout<<std::setw(width)<<"NArray<int, 6>() = inumbers : "               <<(NArray<int, 6>() = inumbers)<<std::endl;

    // Access
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Access : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.data() : "                           <<inumbers.data()<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.data() : "                           <<cnumbers.data()<<std::endl;
    std::cout<<std::setw(width)<<"*inumbers.data() = 5 : "                      <<(*inumbers.data() = 5)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers : "                                  <<inumbers<<std::endl;
    std::cout<<std::setw(width)<<"*inumbers.data() = 4 : "                      <<(*inumbers.data() = 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers : "                                  <<inumbers<<std::endl;
    
    // Iterators
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Iterators : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"std::for_each(big.begin(), big.end(), [](double& x){x *= 2;}) : "                   ; std::for_each(big.begin(), big.end(), [](double& x){x *= 2;}); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"std::for_each(big.begin(), big.end(), [](double x){std::cout<<x<<' ';}) : "         ; std::for_each(big.begin(), big.end(), [](double x){std::cout<<x<<' ';}); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"std::for_each(big.cbegin(), big.cend(), [](double x){std::cout<<x<<' ';}) : "       ; std::for_each(big.cbegin(), big.cend(), [](double x){std::cout<<x<<' ';}); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"std::for_each(big.rbegin(), big.rend(), [](double& x){x /= 2;}) : "                 ; std::for_each(big.rbegin(), big.rend(), [](double& x){x /= 2;}); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"std::for_each(big.rbegin(), big.rend(), [](double x){std::cout<<x<<' ';}) : "       ; std::for_each(big.rbegin(), big.rend(), [](double x){std::cout<<x<<' ';}); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"std::for_each(big.crbegin(), big.crend(), [](double x){std::cout<<x<<' ';}) : "     ; std::for_each(big.crbegin(), big.crend(), [](double x){std::cout<<x<<' ';}); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"big.index(big.begin()+2) : "                                                        <<big.index(big.begin()+2)<<std::endl;

    // Comparison
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Comparison : "                                <<std::endl;
    std::cout<<std::setw(width)<<"(dnumbers/1).null(1) : "                      <<(dnumbers/1.).null(1.)<<std::endl;
    std::cout<<std::setw(width)<<"(dnumbers/10).null(1) : "                     <<(dnumbers/10.).null(1.)<<std::endl;
    std::cout<<std::setw(width)<<"(dnumbers/100).null(1) : "                    <<(dnumbers/100).null(1)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.eq(dnumbers+0.001, 0.01) : "         <<dnumbers.eq(dnumbers+0.001, 0.01)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.eq(dnumbers+0.0099, 0.01) : "        <<dnumbers.eq(dnumbers+0.0099, 0.01)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.eq(dnumbers+0.1, 0.01) : "           <<dnumbers.eq(dnumbers+0.1, 0.01)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.ne(dnumbers+0.001, 0.01) : "         <<dnumbers.ne(dnumbers+0.001, 0.01)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.ne(dnumbers+0.0099, 0.01) : "        <<dnumbers.ne(dnumbers+0.0099, 0.01)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.ne(dnumbers+0.1, 0.01) : "           <<dnumbers.ne(dnumbers+0.1, 0.01)<<std::endl;

    // Statistics
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Statistics : "                                <<std::endl;
    std::cout<<std::setw(width)<<"(inumbers-8).amin() : "                       <<(inumbers-8).amin()<<std::endl;
    std::cout<<std::setw(width)<<"(inumbers-8).amax() : "                       <<(inumbers-8).amax()<<std::endl;
    std::cout<<std::setw(width)<<"(inumbers-42).amin() : "                      <<(inumbers-42).amin()<<std::endl;
    std::cout<<std::setw(width)<<"(inumbers-42).amax() : "                      <<(inumbers-42).amax()<<std::endl;
    std::cout<<std::setw(width)<<"(inumbers-15).amin(mask) : "                  <<(inumbers-15).amin(mask)<<std::endl;
    std::cout<<std::setw(width)<<"(inumbers-15).amax(mask) : "                  <<(inumbers-15).amax(mask)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.mean() : "                           <<inumbers.mean()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.mean(0) : "                          <<inumbers.mean(0)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.mean(inumbers) : "                   <<inumbers.mean(inumbers)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.mean(1, mask) : "                    <<inumbers.mean(1, mask)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.sigma() : "                          <<inumbers.sigma()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.sigma(-1) : "                        <<inumbers.sigma(-1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.sigma(0, inumbers) : "               <<inumbers.sigma(0, inumbers)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.sigma(0, 1, mask) : "                <<inumbers.sigma(0, 1, mask)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.sigma(-1, 1, mask) : "               <<inumbers.sigma(-1, 1, mask)<<std::endl;
    
    // Application
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Application : "                                                                     <<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.modify(std::sqrt) : "                                                      <<dnumbers.modify(std::sqrt)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify<double>(std::sqrt) : "                                              <<inumbers.modify<double>(std::sqrt)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.modify<double, double, double>(std::pow, 2) : "                            <<dnumbers.modify<double, double, double>(std::pow, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify<double, double, double>(std::pow, 2) : "                            <<inumbers.modify<double, double, double>(std::pow, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.apply(std::sin) : "                                                        <<dnumbers.apply(std::sin)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.apply(std::cos) : "                                                        <<dnumbers.apply(std::cos)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply<double>(std::sin) : "                                                <<inumbers.apply<double>(std::sin)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply<double>(std::cos) : "                                                <<inumbers.apply<double>(std::cos)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.apply<double, double, double>(std::pow, 2.) : "                            <<dnumbers.apply<double, double, double>(std::pow, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply<double, double, double>(std::pow, 2.) : "                            <<inumbers.apply<double, double, double>(std::pow, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.apply<double, double, double>(std::pow, dnumbers) : "                      <<dnumbers.apply<double, double, double>(std::pow, dnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply<double, double, double>(std::pow, dnumbers) : "                      <<inumbers.apply<double, double, double>(std::pow, dnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.apply(mask, std::sin) : "                                                  <<dnumbers.apply(mask, std::sin)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers = cnumbers : "                                                             <<(inumbers = cnumbers)<<std::endl;

    // Count
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Count : "                                                                           <<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.unicity() : "                                                              <<dnumbers.unicity()<<std::endl;
    std::cout<<std::setw(width*2)<<"mask.unicity() : "                                                                  <<mask.unicity()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.unicity(std::equal_to<double>(), mask) : "                                 <<dnumbers.unicity(std::equal_to<double>(), mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"mask.unicity(std::equal_to<bool>(), mask) : "                                       <<mask.unicity(std::equal_to<bool>(), mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.distinct() : "                                                             <<dnumbers.distinct()<<std::endl;
    std::cout<<std::setw(width*2)<<"mask.distinct() : "                                                                 <<mask.distinct()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.distinct(std::equal_to<double>(), mask) : "                                <<dnumbers.distinct(std::equal_to<double>(), mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"mask.distinct(std::equal_to<bool>(), mask) : "                                      <<mask.distinct(std::equal_to<bool>(), mask)<<std::endl;
    
    // Sort
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Sort : "                                                                            <<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers = {42, -23, 16, 15, -8, 4} : "                                             <<(dnumbers = {42, -23, 16, 15, -8, 4})<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.arrange() : "                                                              <<dnumbers.arrange()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.arrange(std::greater<double>()) : "                                        <<dnumbers.arrange(std::greater<double>())<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.arrange(std::less<double>(), 0, mask) : "                                  <<dnumbers.arrange(std::less<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.arrange(std::greater<double>(), 0, mask) : "                               <<dnumbers.arrange(std::greater<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.arrange(std::less<double>(), inumbers, mask) : "                           <<dnumbers.arrange(std::less<double>(), inumbers, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.arrange(std::greater<double>(), inumbers, mask) : "                        <<dnumbers.arrange(std::greater<double>(), inumbers, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.sort() : "                                                                 <<dnumbers.sort()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.sort(std::greater<double>()) : "                                           <<dnumbers.sort(std::greater<double>())<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.sort(std::less<double>(), 0, mask) : "                                     <<dnumbers.sort(std::less<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.sort(std::greater<double>(), 0, mask) : "                                  <<dnumbers.sort(std::greater<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.aarrange() : "                                                             <<dnumbers.aarrange()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.aarrange(std::greater<double>()) : "                                       <<dnumbers.aarrange(std::greater<double>())<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.aarrange(std::less<double>(), 0, mask) : "                                 <<dnumbers.aarrange(std::less<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.aarrange(std::greater<double>(), 0, mask) : "                              <<dnumbers.aarrange(std::greater<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.asort() : "                                                                <<dnumbers.asort()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.asort(std::greater<double>()) : "                                          <<dnumbers.asort(std::greater<double>())<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.asort(std::less<double>(), 0, mask) : "                                    <<dnumbers.asort(std::less<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.asort(std::greater<double>(), 0, mask) : "                                 <<dnumbers.asort(std::greater<double>(), 0, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.rearrange() : "                                                            <<dnumbers.rearrange()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.rearrange() : "                                                            <<dnumbers.rearrange()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.resort() : "                                                               <<dnumbers.resort()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.resort() : "                                                               <<dnumbers.resort()<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers = {42, -23, 16, 15, -8, 4} : "                                             <<(dnumbers = {42, -23, 16, 15, -8, 4})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers = cnumbers : "                                                             <<(inumbers = cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.sort(std::less<double>(), inumbers) : "                                    <<dnumbers.sort(std::less<double>(), inumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers = cnumbers : "                                                             <<(inumbers = cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers.sort(std::less<double>(), indexes) : "                                     <<dnumbers.sort(std::less<double>(), indexes)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.resort(indexes) : "                                                        <<inumbers.resort(indexes)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.rearrange(indexes) : "                                                     <<inumbers.rearrange(indexes)<<std::endl;
    std::cout<<std::setw(width*2)<<"dnumbers = {42, 23, 16, 15, 8, 4} : "                                               <<(dnumbers = {42, 23, 16, 15, 8, 4})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers = cnumbers : "                                                             <<(inumbers = cnumbers)<<std::endl;
    
    // Mathematical functions
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Mathematical functions : "                    <<std::endl;
    std::cout<<std::setw(width)<<"(dnumbers-8).abs() : "                        <<(dnumbers-8).abs()<<std::endl;
    std::cout<<std::setw(width)<<"(dnumbers-8).abs(mask) : "                    <<(dnumbers-8).abs(mask)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<-2>() : "                        <<dnumbers.pow<-2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<-1>() : "                        <<dnumbers.pow<-1>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<0>() : "                         <<dnumbers.pow<0>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<1>() : "                         <<dnumbers.pow<1>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<2>() : "                         <<dnumbers.pow<2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<3>() : "                         <<dnumbers.pow<3>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow<4>() : "                         <<dnumbers.pow<4>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow(2) : "                           <<dnumbers.pow(2)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow(dnumbers) : "                    <<dnumbers.pow(dnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.pow(dnumbers, mask) : "              <<dnumbers.pow(dnumbers, mask)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt<-2>() : "                         <<dnumbers.rt<-2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt<2>() : "                          <<dnumbers.rt<2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt<3>() : "                          <<dnumbers.rt<3>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt<4>() : "                          <<dnumbers.rt<4>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt(2) : "                            <<dnumbers.rt(2)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt(dnumbers) : "                     <<dnumbers.rt(dnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.rt(dnumbers, mask) : "               <<dnumbers.rt(dnumbers, mask)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.log() : "                            <<dnumbers.log()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.log<2>() : "                         <<dnumbers.log<2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.log<10>() : "                        <<dnumbers.log<10>()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.log() : "                            <<inumbers.log()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.log<2>() : "                         <<inumbers.log<2>()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.log<3>() : "                         <<inumbers.log<3>()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.log<10>() : "                        <<inumbers.log<10>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.log(2) : "                           <<dnumbers.log(2)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.log(dnumbers, mask) : "              <<dnumbers.log(dnumbers, mask)<<std::endl;

    // Norm
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Norm : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.norm()"                              <<inumbers.norm()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.norm()"                              <<dnumbers.norm()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.norm<0>()"                           <<dnumbers.norm<0>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.norm<1>()"                           <<dnumbers.norm<1>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.norm<2>(mask)"                       <<dnumbers.norm<2>(mask)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.renormalize<2>(mask)"                <<dnumbers.renormalize<2>(mask)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.renormalize<2>()"                    <<dnumbers.renormalize<2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.norm()"                              <<dnumbers.norm()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers = {42, 23, 16, 15, 8, 4} : "         <<(dnumbers = {42, 23, 16, 15, 8, 4})<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.normalize()"                         <<dnumbers.normalize()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.normalize<0>()"                      <<dnumbers.normalize<0>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.normalize<1>()"                      <<dnumbers.normalize<1>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.normalize<2>()"                      <<dnumbers.normalize<2>()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers.normalize<2>(mask)"                  <<dnumbers.normalize<2>(mask)<<std::endl;

    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Predefined : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::zero() : "                                                          <<NArray<int, 6>::zero()<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::one() : "                                                           <<NArray<int, 6>::one()<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::value() : "                                                         <<NArray<int, 6>::value()<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::value(42) : "                                                       <<NArray<int, 6>::value(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::value(42, mask) : "                                                 <<NArray<int, 6>::value(42, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<double, 6>::random() : "                                                     <<NArray<double, 6>::random()<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<double, 6>::random(-10, 10) : "                                              <<NArray<double, 6>::random(-10, 10)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<double, 6>::random(-10, 10, mask) : "                                        <<NArray<double, 6>::random(-10, 10, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::random(engine, distribution) : "                                    <<NArray<int, 6>::random(engine, distribution)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::random(engine, distribution) : "                                    <<NArray<int, 6>::random(engine, distribution)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::random(engine, distribution, mask) : "                              <<NArray<int, 6>::random(engine, distribution, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::indexed() : "                                                       <<NArray<int, 6>::indexed()<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::indexed([](unsigned int i){return 2*i;}) : "                        <<NArray<int, 6>::indexed([](unsigned int i){return 2*i;})<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::indexed([](unsigned int i){return i;}, mask) : "                    <<NArray<int, 6>::indexed([](unsigned int i){return i;}, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::progressive() : "                                                   <<NArray<int, 6>::progressive()<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::progressive(42) : "                                                 <<NArray<int, 6>::progressive(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::progressive(42, 2) : "                                              <<NArray<int, 6>::progressive(42, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::progressive(1, 2, std::multiplies<int>) : "                         <<NArray<int, 6>::progressive(1, 2, std::multiplies<int>())<<std::endl;
    std::cout<<std::setw(width*2)<<"NArray<int, 6>::progressive(1, 2, std::multiplies<int>, mask) : "                   <<NArray<int, 6>::progressive(1, 2, std::multiplies<int>(), mask)<<std::endl;    
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = NArray::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // NARRAY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/


