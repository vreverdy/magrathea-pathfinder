/* ****************************** STATICVECTOR ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          StaticVector
// DESCRIPTION :    Basic vectorized constant size container
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           staticvector.h
/// \brief          Basic vectorized constant size container
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef STATICVECTOR_H_INCLUDED
#define STATICVECTOR_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <initializer_list>
#include <array>
// Include libs
// Include project
#include "staticvectorizer.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Basic vectorized constant size container
/// \brief          Basic vectorized constant size container.
/// \details        This class is the direct derivation of StaticVectorizer. It 
///                 provides the most basic vectorized constant size container 
///                 without adding new functionalities to the abstract class.
/// \tparam         Type Data type.
/// \tparam         Size Number of elements.
template <typename Type = double, unsigned int Size = 1>
class StaticVector final
: public StaticVectorizer<unsigned int, Size, StaticVector, Type, Size>
{
    // Setup
    public: using StaticVectorizer<unsigned int, Size, magrathea::StaticVector, Type, Size>::operator=;
    
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        inline StaticVector();
        template <typename FundamentalType = Type, class = typename std::enable_if<std::is_fundamental<FundamentalType>::value>::type> inline StaticVector(const StaticVector<FundamentalType, Size>& source);
        template <typename OtherType = Type, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline StaticVector(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc, class = typename std::enable_if<sizeof...(Misc) != 0>::type> explicit inline StaticVector(const Misc&... misc);
    //@}
    
    // Vectorization
    /// \name           Vectorization
    //@{
    public:
        inline Type& operator[](const unsigned int i);
        inline const Type& operator[](const unsigned int i) const;
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
inline StaticVector<Type, Size>::StaticVector()
: StaticVectorizer<unsigned int, Size, StaticVector, Type, Size>()
, _data()
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
inline StaticVector<Type, Size>::StaticVector(const StaticVector<FundamentalType, Size>& source)
: StaticVectorizer<unsigned int, Size, StaticVector, Type, Size>()
, _data()
{
    Vectorizer::set(*this, source);
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
inline StaticVector<Type, Size>::StaticVector(const std::initializer_list<OtherType>& source, const Misc&... misc)
: StaticVectorizer<unsigned int, Size, StaticVector, Type, Size>()
, _data()
{
    Vectorizer::set(*this, source, misc...);
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
inline StaticVector<Type, Size>::StaticVector(const Misc&... misc)
: StaticVectorizer<unsigned int, Size, StaticVector, Type, Size>()
, _data()
{
    Vectorizer::set(*this, misc...);
}
// -------------------------------------------------------------------------- //



// ----------------------------- VECTORIZATION ------------------------------ //
// Direct access to the element
/// \brief          Direct access to the element.
/// \details        Provides a direct access to the specified element.
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <typename Type, unsigned int Size>
inline Type& StaticVector<Type, Size>::operator[](const unsigned int i)
{
    return _data[i];
}

// Immutable direct access to the element
/// \brief          Immutable direct access to the element.
/// \details        Provides a constant direct access to the specified element.
/// \param[in]      i Index of the element.
/// \return         Const reference to the element.
template <typename Type, unsigned int Size>
inline const Type& StaticVector<Type, Size>::operator[](const unsigned int i) const
{
    return _data[i];
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of StaticVector.
/// \return         0 if no error.
template <typename Type, unsigned int Size>
int StaticVector<Type, Size>::example()
{
    // Initialize
    std::cout<<"BEGIN = StaticVector::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    std::istringstream iss("1 2 3 4 5 6");
    
    // Construction
    const StaticVector<int, 6> cnumbers({4, 8, 15, 16, 23, 42});
    StaticVector<int, 6> inumbers({4, 8, 15, 16, 23, 42});
    StaticVector<int, 6> inumbers1({4, 8, 15, 16, 23, 42});
    StaticVector<int, 6> inumbers2({1, 2, 3, 4, 5, 6});
    StaticVector<double, 6> dnumbers1({42, 23, 16, 15, 8, 4});
    StaticVector<double, 6> dnumbers2({4, 8, 15, 16, 23, 42});
    StaticVector<bool, 6> mask({true, true, false, true, false, true});
    StaticVector<int, 6> ten(10);
    StaticVector<double, 4> small(42);
    StaticVector<double, 9> big({42, 23, 0, 15, 8, 4, 2, 1, 0});
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<>() : "                          <<StaticVector<>()<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<double, 6>() : "                 <<StaticVector<double, 6>()<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<double, 6>(dnumbers1) : "        <<StaticVector<double, 6>(dnumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<double, 6>(cnumbers) : "         <<StaticVector<double, 6>(cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<double, 6>({42}) : "             <<StaticVector<double, 6>({42})<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<int, 4>({0,1,2,3,4}) : "         <<StaticVector<int, 4>({0, 1, 2, 3, 4})<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<bool, 4>({0, 1, 0, 1}) : "       <<StaticVector<bool, 4>({0, 1, 0, 1})<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<double, 6>(42) : "               <<StaticVector<double, 6>(42)<<std::endl;
    std::cout<<std::setw(width)<<"StaticVector<int, 6>() = inumbers : "         <<(StaticVector<int, 6>() = inumbers)<<std::endl;
    
    // Vectorization
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Vectorization : "                             <<std::endl;
    std::cout<<std::setw(width)<<"inumbers[0] : "                               <<inumbers[0]<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers[0] : "                               <<cnumbers[0]<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.resize(6) : "                        <<inumbers.resize(6)<<std::endl;
    
    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"dnumbers2 = dnumbers1 : "                     <<(dnumbers2 = dnumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers2 = {4, 8, 15, 16, 23, 42} : "        <<(dnumbers2 = {4, 8, 15, 16, 23, 42})<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers2 = 42 : "                            <<(dnumbers2 = 42)<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers2 = cnumbers : "                      <<(dnumbers2 = cnumbers)<<std::endl;
    
    // Compound assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Compound assignment : "                       <<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 += inumbers2 : "                    <<(inumbers1 += inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 -= inumbers2 : "                    <<(inumbers1 -= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 *= inumbers2 : "                    <<(inumbers1 *= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 /= inumbers2 : "                    <<(inumbers1 /= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 %= inumbers2 : "                    <<(inumbers1 %= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 &= inumbers2 : "                    <<(inumbers1 &= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 |= inumbers2 : "                    <<(inumbers1 |= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 ^= inumbers2 : "                    <<(inumbers1 ^= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 += inumbers2 : "                    <<(inumbers1 += inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 <<= inumbers2 : "                   <<(inumbers1 <<= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 >>= inumbers2 : "                   <<(inumbers1 >>= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 += 4 : "                            <<(inumbers1 += 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 -= 4 : "                            <<(inumbers1 -= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 *= 4 : "                            <<(inumbers1 *= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 /= 4 : "                            <<(inumbers1 /= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 %= 4 : "                            <<(inumbers1 %= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 &= 4 : "                            <<(inumbers1 &= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 |= 4 : "                            <<(inumbers1 |= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 ^= 4 : "                            <<(inumbers1 ^= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 += 4 : "                            <<(inumbers1 += 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 <<= 4 : "                           <<(inumbers1 <<= 4)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 >>= 4 : "                           <<(inumbers1 >>= 4)<<std::endl;

    // Operators : main
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : main : "                          <<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 = cnumbers : "                      <<(inumbers1 = cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers2 : "                                 <<(inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 + inumbers2 : "                     <<(inumbers1 + inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 - inumbers2 : "                     <<(inumbers1 - inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 * inumbers2 : "                     <<(inumbers1 * inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 / inumbers2 : "                     <<(inumbers1 / inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 % inumbers2 : "                     <<(inumbers1 % inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 & inumbers2 : "                     <<(inumbers1 & inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 | inumbers2 : "                     <<(inumbers1 | inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 ^ inumbers2 : "                     <<(inumbers1 ^ inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 << inumbers2 : "                    <<(inumbers1 << inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 >> inumbers2 : "                    <<(inumbers1 >> inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 && inumbers2 : "                    <<(inumbers1 && inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 || inumbers2 : "                    <<(inumbers1 || inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 == inumbers2 : "                    <<(inumbers1 == inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 != inumbers2 : "                    <<(inumbers1 != inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 > inumbers2 : "                     <<(inumbers1 > inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 < inumbers2 : "                     <<(inumbers1 < inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 >= inumbers2 : "                    <<(inumbers1 >= inumbers2)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 <= inumbers2 : "                    <<(inumbers1 <= inumbers2)<<std::endl;

    // Operators : with rhs value
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : with rhs value : "                <<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 + 15 : "                            <<(inumbers1 + 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 - 15 : "                            <<(inumbers1 - 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 * 15 : "                            <<(inumbers1 * 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 / 15 : "                            <<(inumbers1 / 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 % 15 : "                            <<(inumbers1 % 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 & 15 : "                            <<(inumbers1 & 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 | 15 : "                            <<(inumbers1 | 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 ^ 15 : "                            <<(inumbers1 ^ 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 << 15 : "                           <<(inumbers1 << 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 >> 15 : "                           <<(inumbers1 >> 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 && 15 : "                           <<(inumbers1 && 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 || 15 : "                           <<(inumbers1 || 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 == 15 : "                           <<(inumbers1 == 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 != 15 : "                           <<(inumbers1 != 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 > 15 : "                            <<(inumbers1 > 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 < 15 : "                            <<(inumbers1 < 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 >= 15 : "                           <<(inumbers1 >= 15)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 <= 15 : "                           <<(inumbers1 <= 15)<<std::endl;

    // Operators : with lhs value
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : with lhs value : "                <<std::endl;
    std::cout<<std::setw(width)<<"15 + inumbers1 : "                            <<(15 + inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 - inumbers1 : "                            <<(15 - inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 * inumbers1 : "                            <<(15 * inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 / inumbers1 : "                            <<(15 / inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 % inumbers1 : "                            <<(15 % inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 & inumbers1 : "                            <<(15 & inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 | inumbers1 : "                            <<(15 | inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 ^ inumbers1 : "                            <<(15 ^ inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 << inumbers1 : "                           <<(15 << inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 >> inumbers1 : "                           <<(15 >> inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 && inumbers1 : "                           <<(15 && inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 || inumbers1 : "                           <<(15 || inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 == inumbers1 : "                           <<(15 == inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 != inumbers1 : "                           <<(15 != inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 > inumbers1 : "                            <<(15 > inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 < inumbers1 : "                            <<(15 < inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 >= inumbers1 : "                           <<(15 >= inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"15 <= inumbers1 : "                           <<(15 <= inumbers1)<<std::endl;

    // Operators : unary
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : unary : "                         <<std::endl;
    std::cout<<std::setw(width)<<"!inumbers1 : "                                <<(!inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"~inumbers1 : "                                <<(~inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"+inumbers1 : "                                <<(+inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"-inumbers1 : "                                <<(-inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"++inumbers1 : "                               <<(++inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"--inumbers1 : "                               <<(--inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1++ : "                               <<(inumbers1++)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1-- : "                               <<(inumbers1--)<<std::endl;

    // Access
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Access : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"inumbers() : "                                <<inumbers()<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers() : "                                <<cnumbers()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers(1) : "                               <<inumbers(1)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers(1) : "                               <<cnumbers(1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.at(1) : "                            <<inumbers.at(1)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.at(1) : "                            <<cnumbers.at(1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.front() : "                          <<inumbers.front()<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.front() : "                          <<cnumbers.front()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.front(1) : "                         <<inumbers.front(1)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.front(1) : "                         <<cnumbers.front(1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.back() : "                           <<inumbers.back()<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.back() : "                           <<cnumbers.back()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.back(1) : "                          <<inumbers.back(1)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.back(1) : "                          <<cnumbers.back(1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.cycle(1) : "                         <<inumbers.cycle(1)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.cycle(1) : "                         <<cnumbers.cycle(1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.cycle(-1) : "                        <<inumbers.cycle(-1)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.cycle(-1) : "                        <<cnumbers.cycle(-1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.cycle(6) : "                         <<inumbers.cycle(6)<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.cycle(6) : "                         <<cnumbers.cycle(6)<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Assignment : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign({0, 1, 2, 3}) : "                                                   <<inumbers.assign({0, 1, 2, 3})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign({0, 1, 2, 3}, 1, 2) : "                                             <<inumbers.assign({0, 1, 2, 3}, 1, 3)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(cnumbers) : "                                                       <<inumbers.assign(cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(cnumbers, 1, 2, 4, 2) : "                                           <<inumbers.assign(cnumbers, 1, 2, 4, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(cnumbers, 1, 2) : "                                                 <<inumbers.assign(cnumbers, 1, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(cnumbers, mask) : "                                                 <<inumbers.assign(cnumbers, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(42) : "                                                             <<inumbers.assign(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(42, 1, 2) : "                                                       <<inumbers.assign(42, 1, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(42) : "                                                             <<inumbers.assign(42, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.assign(&big[0],&big[0]+big.size()) : "                                     <<inumbers.assign(&big[0], &big[0]+big.size())<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill({0, 1, 2, 3}) : "                                                     <<inumbers.fill({0, 1, 2, 3})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill({0, 1, 2, 3}, 1, 2) : "                                               <<inumbers.fill({0, 1, 2, 3}, 1, 3)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(cnumbers) : "                                                         <<inumbers.fill(cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(cnumbers, 1, 2, 4, 2) : "                                             <<inumbers.fill(cnumbers, 1, 2, 4, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(cnumbers, 1, 2) : "                                                   <<inumbers.fill(cnumbers, 1, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(cnumbers, mask) : "                                                   <<inumbers.fill(cnumbers, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(42) : "                                                               <<inumbers.fill(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(42, 1, 2) : "                                                         <<inumbers.fill(42, 1, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(42) : "                                                               <<inumbers.fill(42, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.fill(&big[0],&big[0]+big.size()) :"                                        <<inumbers.fill(&big[0], &big[0]+big.size())<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace({0, 1, 2, 3}) : "                                                  <<inumbers.replace({0, 1, 2, 3})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace({0, 1, 2, 3}, 1, 2) : "                                            <<inumbers.replace({0, 1, 2, 3}, 1, 3)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(cnumbers) : "                                                      <<inumbers.replace(cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(cnumbers, 1, 2, 4, 2) : "                                          <<inumbers.replace(cnumbers, 1, 2, 4, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(cnumbers, 1, 2) : "                                                <<inumbers.replace(cnumbers, 1, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(cnumbers, mask) : "                                                <<inumbers.replace(cnumbers, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(42) : "                                                            <<inumbers.replace(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(42, 1, 2) : "                                                      <<inumbers.replace(42, 1, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(42) : "                                                            <<inumbers.replace(42, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.replace(&big[0],&big[0]+big.size()) : "                                    <<inumbers.replace(&big[0], &big[0]+big.size())<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.put(0, 1) : "                                                              <<inumbers.put(0, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.put(0, 2, 2) : "                                                           <<inumbers.put(0, 2, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.put(cnumbers, 1) : "                                                       <<inumbers.put(cnumbers, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.put(cnumbers, 2, 2) : "                                                    <<inumbers.put(cnumbers, 2, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.change(0, 1) : "                                                           <<inumbers.change(0, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.change(0, 2, 2) : "                                                        <<inumbers.change(0, 2, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.change(dnumbers1, 1) : "                                                   <<inumbers.change(dnumbers1, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.change(dnumbers1, 2, 2) : "                                                <<inumbers.change(dnumbers1, 2, 2)<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.reserve(6) : "                       <<inumbers.reserve(6)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.nullify() : "                        <<inumbers.nullify()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.swap(inumbers1) : "                  <<inumbers.swap(inumbers1)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.nullify(mask) : "                    <<inumbers.nullify(mask)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.swap(inumbers1, !mask) : "           <<inumbers.swap(inumbers1, !mask)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers2.copy() : "                          <<inumbers2.copy()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers2.cast<double>() : "                  <<inumbers2.cast<double>()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers = cnumbers : "                       <<(inumbers = cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1 = 0 : "                             <<(inumbers1 = 0)<<std::endl;

    // Comparison
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Comparison : "                                <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.null() : "                           <<inumbers.null()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers1.null() : "                          <<inumbers1.null()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.eq(cnumbers) : "                     <<inumbers.eq(cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.ne(cnumbers) : "                     <<inumbers.ne(cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.eq(42) : "                           <<inumbers.eq(42)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.ne(42) : "                           <<inumbers.ne(42)<<std::endl;
    
    // Statistics
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Statistics : "                                <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.min() : "                            <<inumbers.min()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.max() : "                            <<inumbers.max()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers1.min() : "                           <<dnumbers1.min()<<std::endl;
    std::cout<<std::setw(width)<<"dnumbers1.max() : "                           <<dnumbers1.max()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.min(mask) : "                        <<inumbers.min(mask)<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.max(mask) : "                        <<inumbers.max(mask)<<std::endl;

    // Application
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Application : "                                                                     <<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"mask : "                                                                            <<mask<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify([](int x){return x+x;}) : "                                         <<inumbers.modify([](int x){return x+x;})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify([](int x, int y){return x+y;}, 2) : "                               <<inumbers.modify([](int x, int y){return x+y;}, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify([](int x, int y){return x+y;}, cnumbers) : "                        <<inumbers.modify([](int x, int y){return x+y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify(mask, [](int x){return x+x;}) : "                                   <<inumbers.modify(mask, [](int x){return x+x;})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify(mask, [](int x, int y){return x+y;}, 2) : "                         <<inumbers.modify(mask, [](int x, int y){return x+y;}, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.modify(mask, [](int x, int y){return x+y;}, cnumbers) : "                  <<inumbers.modify(mask, [](int x, int y){return x+y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply([](int x){return x+x;}) : "                                          <<inumbers.apply([](int x){return x+x;})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply([](int x, int y){return x+y;}, 2) : "                                <<inumbers.apply([](int x, int y){return x+y;}, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply([](int x, int y){return x+y;}, cnumbers) : "                         <<inumbers.apply([](int x, int y){return x+y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply([](double x, double y){return x/y;}, cnumbers) : "                   <<inumbers.apply([](double x, double y){return x/y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply<double>([](double x, double y){return x/y;}, cnumbers) : "           <<inumbers.apply<double>([](double x, double y){return x/y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply(mask, [](int x){return x+x;}) : "                                    <<inumbers.apply(mask, [](int x){return x+x;})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply(mask, [](int x, int y){return x+y;}, 2) : "                          <<inumbers.apply(mask, [](int x, int y){return x+y;}, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply(mask, [](int x, int y){return x+y;}, cnumbers) : "                   <<inumbers.apply(mask, [](int x, int y){return x+y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply(mask, [](double x, double y){return x/y;}, cnumbers) : "             <<inumbers.apply(mask, [](double x, double y){return x/y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.apply<double>(mask, [](double x, double y){return x/y;}, cnumbers) : "     <<inumbers.apply<double>(mask, [](double x, double y){return x/y;}, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"mask : "                                                                            <<mask<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce() : "                                                               <<inumbers.reduce()<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce(std::plus<int>()) : "                                               <<inumbers.reduce(std::plus<int>())<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce(std::plus<int>(), 10) : "                                           <<inumbers.reduce(std::plus<int>(), 10)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce(mask) : "                                                           <<inumbers.reduce(mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce(mask, std::plus<int>()) : "                                         <<inumbers.reduce(mask, std::plus<int>())<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce([](int x, int y){return std::min(x, y);}) : "                       <<inumbers.reduce([](int x, int y){return std::min(x, y);})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.reduce([](int x, int y){return std::max(x, y);}) : "                       <<inumbers.reduce([](int x, int y){return std::max(x, y);})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.combine(std::multiplies<double>(), inumbers2, dnumbers2) : "               <<inumbers.combine(std::multiplies<double>(), inumbers2, dnumbers2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.combine(std::multiplies<double>()) : "                                     <<inumbers.combine(std::multiplies<double>())<<std::endl;
    
    // Count
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Count : "                                                                           <<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers[3] = 0 : "                                                                 <<(inumbers[3] = 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers : "                                                                        <<inumbers<<std::endl;
    std::cout<<std::setw(width*2)<<"mask : "                                                                            <<mask<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count() : "                                                                <<inumbers.count()<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count(34) : "                                                              <<inumbers.count(34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count(34.5) : "                                                            <<inumbers.count(34.5)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count(dnumbers1.change(0, 3)) : "                                          <<inumbers.count(dnumbers1.change(0, 3))<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count(true, mask) : "                                                      <<inumbers.count(true, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count([](int x){return x > 64;}) : "                                       <<inumbers.count([](int x){return x > 64;})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.count([](int x){return x > 64;}, mask) : "                                 <<inumbers.count([](int x){return x > 64;}, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.all(34) : "                                                                <<inumbers.all(34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.all(34, inumbers == 34) : "                                                <<inumbers.all(34, inumbers == 34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.all(34, inumbers == 256) : "                                               <<inumbers.all(34, inumbers == 256)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.any(34) : "                                                                <<inumbers.any(34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.any(34, inumbers == 34) : "                                                <<inumbers.any(34, inumbers == 34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.any(34, inumbers == 256) : "                                               <<inumbers.any(34, inumbers == 256)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.none(34) : "                                                               <<inumbers.none(34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.none(34, inumbers == 34) : "                                               <<inumbers.none(34, inumbers == 34)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.none(34, inumbers == 256) : "                                              <<inumbers.none(34, inumbers == 256)<<std::endl;

    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, inumbers1) : "          <<inumbers1<<std::endl;
    std::cout<<std::setw(width)<<"operator>>(iss, inumbers1) : "                ; iss>>inumbers1; std::cout<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, inumbers1) : "          <<inumbers1<<std::endl;
    
    // Static vectorization
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Static vectorization : "                      <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.size() : "                           <<inumbers.size()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.constant() : "                       <<inumbers.constant()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.boolean() : "                        <<inumbers.boolean()<<std::endl;
    std::cout<<std::setw(width)<<"mask.boolean() : "                            <<mask.boolean()<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.parameters()[0] : "                  <<cnumbers.parameters()[0]<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.type() : "                           <<inumbers.type()<<std::endl;
    
    // Size
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Size : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.empty() : "                          <<inumbers.empty()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.capacity() : "                       <<inumbers.capacity()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.tbytes() : "                         <<inumbers.tbytes()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.bytes() : "                          <<inumbers.bytes()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.space() : "                          <<inumbers.space()<<std::endl;

    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Predefined : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask() : "                                                                 <<inumbers.mask()<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask(false) : "                                                            <<inumbers.mask(false)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask(std::bitset<5>(30)) : "                                               <<inumbers.mask(std::bitset<5>(30))<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask(std::vector<bool>(6)) : "                                             <<inumbers.mask(std::vector<bool>(6))<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask({0, 1, 0, 1}) : "                                                     <<inumbers.mask({0, 1, 0, 1})<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask({0, 1, 0, 1}, 2) : "                                                  <<inumbers.mask({0, 1, 0, 1}, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask(42) : "                                                               <<inumbers.mask(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"inumbers.mask(&big[0], &big[0]+big.size()) : "                                      <<inumbers.mask(&big[0], &big[0]+big.size())<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = StaticVector::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // STATICVECTOR_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
