/* ******************************* VECTORIZED ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Vectorized
// DESCRIPTION :    Basic vectorized container
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           vectorized.h
/// \brief          Basic vectorized container
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013 
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef VECTORIZED_H_INCLUDED
#define VECTORIZED_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <initializer_list>
#include <type_traits>
#include <array>
// Include libs
// Include project
#include "vectorizer.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Basic vectorized container
/// \brief          Basic vectorized container.
/// \details        This class is the direct derivation of Vectorizer. It 
///                 provides the most basic vectorized container without adding 
///                 new functionalities to the abstract class.
/// \tparam         Type Data type.
/// \tparam         Size Number of elements.
template <typename Type = double, unsigned int Size = 1>
class Vectorized final
: public Vectorizer
{
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        inline Vectorized();
        template <typename FundamentalType = Type, class = typename std::enable_if<std::is_fundamental<FundamentalType>::value>::type> inline Vectorized(const Vectorized<FundamentalType, Size>& source);
        template <typename OtherType = Type, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Vectorized(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc, class = typename std::enable_if<sizeof...(Misc) != 0>::type> explicit inline Vectorized(const Misc&... misc);
    //@}
    
    // Vectorization
    /// \name           Vectorization
    //@{
    public:
        inline Type& operator[](const unsigned int i);
        inline const Type& operator[](const unsigned int i) const;
        inline Vectorized<Type, Size>& resize(const unsigned int n);
        inline Vectorized<Type, Size>& reserve(const unsigned int n);
    //@}

    // Stream
    /// \name           Stream
    //@{
    public:
        template <typename SelfType, unsigned int SelfSize> friend std::ostream& operator<<(std::ostream& lhs, const Vectorized<SelfType, SelfSize>& rhs);
    //@}

    // Static vectorization
    /// \name           Static vectorization
    //@{
    public:
        static constexpr unsigned int size();
        static constexpr unsigned int capacity();
        static constexpr bool constant();
        static constexpr bool boolean();
        static constexpr std::array<unsigned int, 1> parameters();
        static inline Type type();
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
inline Vectorized<Type, Size>::Vectorized()
: Vectorizer() 
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
inline Vectorized<Type, Size>::Vectorized(const Vectorized<FundamentalType, Size>& source)
: Vectorizer() 
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
inline Vectorized<Type, Size>::Vectorized(const std::initializer_list<OtherType>& source, const Misc&... misc)
: Vectorizer() 
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
inline Vectorized<Type, Size>::Vectorized(const Misc&... misc)
: Vectorizer() 
, _data()
{
    Vectorizer::set(*this, misc...);
}
// -------------------------------------------------------------------------- //



// ----------------------------- VECTORIZATION ------------------------------ //
// Direct access to the element
/// \brief          Direct access to the element.
/// \details        Provides a direct access to the specified element. This 
///                 function is required by the vectorization mechanism. 
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <typename Type, unsigned int Size>
inline Type& Vectorized<Type, Size>::operator[](const unsigned int i)
{
    return _data[i];
}

// Immutable direct access to the element
/// \brief          Immutable direct access to the element.
/// \details        Provides a constant direct access to the specified element.
///                 This function is required by the vectorization mechanism.
/// \param[in]      i Index of the element.
/// \return         Const reference to the element.
template <typename Type, unsigned int Size>
inline const Type& Vectorized<Type, Size>::operator[](const unsigned int i) const
{
    return _data[i];
}

// Resize the container
/// \brief          Resize the container.
/// \details        Resizes the container and returns a reference to it. This 
///                 function is required by the vectorization mechanism. 
/// \param[in]      n New size.
/// \return         Self reference.
/// \exception      std::length_error The container cannot be resized.
template <typename Type, unsigned int Size>
inline Vectorized<Type, Size>& Vectorized<Type, Size>::resize(const unsigned int n)
{
    if (n != Size) {
        throw std::length_error("ERROR = Vectorized::resize() : the container cannot be resized");
    }
    return *this;
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Adds each element to the stream using the fill character to
///                 separate the elements.
/// \tparam         SelfType Data type.
/// \tparam         SelfSize Number of elements.
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side container.
/// \return         Output stream.
template <typename SelfType, unsigned int SelfSize>
std::ostream& operator<<(std::ostream& lhs, const Vectorized<SelfType, SelfSize>& rhs)
{
    const int lim = rhs.size()-1;
    const char spacer = lhs.fill();
    for (int i = 0; i < lim; ++i) {
        lhs<<Vectorizer::get(rhs, i)<<spacer; 
    }
    return lhs<<Vectorizer::get(rhs, lim);
}
// -------------------------------------------------------------------------- //



// -------------------------- STATIC VECTORIZATION -------------------------- //
// Get the size of the container
/// \brief          Get the size of the container.
/// \details        Returns the current number of elements. This function is 
///                 required by the vectorization mechanism. 
/// \return         Copy of the size.
template <typename Type, unsigned int Size>
constexpr unsigned int Vectorized<Type, Size>::size()
{
    return Size;
}

// Get whether the container has a constant size
/// \brief          Get whether the container has a constant size.
/// \details        Returns true if the container has a constant size, false
///                 otherwise. This function is required by the vectorization 
///                 mechanism. 
/// \return         Copy of true.
template <typename Type, unsigned int Size>
constexpr bool Vectorized<Type, Size>::constant()
{
    return true;
}

// Get whether the container has a boolean type
/// \brief          Get whether the container has a boolean type.
/// \details        Returns true if the container has a boolean type, false
///                 otherwise. This function is required by the vectorization 
///                 mechanism. 
/// \return         Copy of true if the container has a boolean type.
template <typename Type, unsigned int Size>
constexpr bool Vectorized<Type, Size>::boolean()
{
    return std::is_same<bool, Type>::value;
}

// Get the template parameters
/// \brief          Get the template parameters.
/// \details        Returns an array containing the template parameters. This 
///                 function is required by the vectorization mechanism.
/// \return         Copy of an array of parameters.
template <typename Type, unsigned int Size>
constexpr std::array<unsigned int, 1> Vectorized<Type, Size>::parameters()
{
    return std::array<unsigned int, 1>{{Size}};
}

// Get the data type
/// \brief          Get the data type.
/// \details        Returns a copy of the default value of the data type.
/// \return         Copy of the default value of the data type.
template <typename Type, unsigned int Size>
inline Type Vectorized<Type, Size>::type()
{
    return Type();
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Vectorized.
/// \return         0 if no error.
template <typename Type, unsigned int Size>
int Vectorized<Type, Size>::example()
{
    // Initialize
    std::cout<<"BEGIN = Vectorized::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    std::integral_constant<int, 10> cst;
    
    // Construction
    Vectorized<> v;
    const Vectorized<int, 6> cnumbers({4, 8, 15, 16, 23, 42});
    Vectorized<int, 6> inumbers({4, 8, 15, 16, 23, 42});
    Vectorized<double, 6> dnumbers({42, 23, 16, 15, 8, 4});
    Vectorized<bool, 6> mask({true, true, false, true, false, true});
    Vectorized<int, 6> ten(10);
    Vectorized<double, 4> small(42);
    Vectorized<double, 9> big({42, 23, 16, 15, 8, 4, 2, 1, 0});
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<>() : "                            <<Vectorized<>()<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<double, 6>() : "                   <<Vectorized<double, 6>()<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<double, 6>(dnumbers) : "           <<Vectorized<double, 6>(dnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<double, 6>(cnumbers) : "           <<Vectorized<double, 6>(cnumbers)<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<double, 6>({42}) : "               <<Vectorized<double, 6>({42})<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<int, 4>({0, 1, 2, 3, 4}) : "       <<Vectorized<int, 4>({0, 1, 2, 3, 4})<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<bool, 4>({0, 1, 0, 1}) : "         <<Vectorized<bool, 4>({0, 1, 0, 1})<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<double, 6>(42) : "                 <<Vectorized<double, 6>(42)<<std::endl;
    std::cout<<std::setw(width)<<"Vectorized<int, 6>() = inumbers : "           <<(Vectorized<int, 6>() = inumbers)<<std::endl;

    // Vectorization
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Vectorization : "                             <<std::endl;
    std::cout<<std::setw(width)<<"inumbers[0] : "                               <<inumbers[0]<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers[0] : "                               <<cnumbers[0]<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.resize(6) : "                        <<inumbers.resize(6)<<std::endl;
    
    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, inumbers) : "           <<inumbers<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, dnumbers) : "           <<dnumbers<<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, ten) : "                <<ten<<std::endl;
    
    // Static vectorization
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Static vectorization : "                      <<std::endl;
    std::cout<<std::setw(width)<<"inumbers.size() : "                           <<inumbers.size()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.constant() : "                       <<inumbers.constant()<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.boolean() : "                        <<inumbers.boolean()<<std::endl;
    std::cout<<std::setw(width)<<"mask.boolean() : "                            <<mask.boolean()<<std::endl;
    std::cout<<std::setw(width)<<"cnumbers.parameters()[0] : "                  <<cnumbers.parameters()[0]<<std::endl;
    std::cout<<std::setw(width)<<"inumbers.type() : "                           <<inumbers.type()<<std::endl;

    // Check
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Check : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"v.check(inumbers, 10) : "                     <<v.check(inumbers, 10)<<std::endl;
    std::cout<<std::setw(width)<<"v.check(10, inumbers) : "                     <<v.check(10, inumbers)<<std::endl;
    std::cout<<std::setw(width)<<"v.check(inumbers, dnumbers) : "               <<v.check(inumbers, dnumbers)<<std::endl;

    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"v.get(cst) : "                                <<v.get(cst)<<std::endl;
    std::cout<<std::setw(width)<<"v.get(std::false_type()) : "                  <<v.get(std::false_type())<<std::endl;
    std::cout<<std::setw(width)<<"v.get(std::true_type()) : "                   <<v.get(std::true_type())<<std::endl;
    std::cout<<std::setw(width)<<"v.get(42) : "                                 <<v.get(42)<<std::endl;
    std::cout<<std::setw(width)<<"v.get(42, 1) : "                              <<v.get(42, 1)<<std::endl;
    std::cout<<std::setw(width)<<"v.get(width, 1) : "                           <<v.get(width, 1)<<std::endl;
    std::cout<<std::setw(width)<<"v.get(inumbers, 1) : "                        <<v.get(inumbers, 1)<<std::endl;
    std::cout<<std::setw(width)<<"v.get(cnumbers, 1) : "                        <<v.get(cnumbers, 1)<<std::endl;

    // Setters
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Setters : "                                                                         <<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers) : "                                                                 <<v.set(inumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 5) : "                                                              <<v.set(inumbers, 5)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, cnumbers) : "                                                       <<v.set(inumbers, cnumbers)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 0, 1, 2, 4, 1) : "                                                  <<v.set(inumbers, 0, 1, 2, 4, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 42, 2) : "                                                          <<v.set(inumbers, 42, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, cnumbers, 2, 2) : "                                                 <<v.set(inumbers, cnumbers, 2, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 42, 2, 10) : "                                                      <<v.set(inumbers, 42, 2, 10)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, {1., 2., 3.}) : "                                                   <<v.set(inumbers, {1., 2., 3.})<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, {4., 5., 6.}, 2, 2) : "                                             <<v.set(inumbers, {4., 5., 6.}, 2, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, {70, 80, 90}, 1, 10) : "                                            <<v.set(inumbers, {70, 80, 90}, 1, 10)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+4) : "                                             <<v.set(inumbers, &big[0], &big[0]+4)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+9) : "                                             <<v.set(inumbers, &big[0], &big[0]+9)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 0) : "                                                              <<v.set(inumbers, 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+3, 2) : "                                          <<v.set(inumbers, &big[0], &big[0]+3, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+20, 2) : "                                         <<v.set(inumbers, &big[0], &big[0]+20, 2)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 0) : "                                                              <<v.set(inumbers, 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+3, 1, 1) : "                                       <<v.set(inumbers, &big[0], &big[0]+3, 1, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+20, 2, 1) : "                                      <<v.set(inumbers, &big[0], &big[0]+20, 2, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 0) : "                                                              <<v.set(inumbers, 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+3, 1, 4) : "                                       <<v.set(inumbers, &big[0], &big[0]+3, 1, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set((inumbers, &big[0], &big[0]+20, 2, 4) : "                                     <<v.set(inumbers, &big[0], &big[0]+20, 2, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 0) : "                                                              <<v.set(inumbers, 0)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 42, mask) : "                                                       <<v.set(inumbers, 42, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, cnumbers, mask) : "                                                 <<v.set(inumbers, cnumbers, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 42) : "                                                             <<v.set(inumbers, 42)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, {1., 2., 3.}, mask) : "                                             <<v.set(inumbers, {1., 2., 3.}, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, 1) : "                                                              <<v.set(inumbers, 1)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+3, mask) : "                                       <<v.set(inumbers, &big[0], &big[0]+3, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, &big[0], &big[0]+20, mask) : "                                      <<v.set(inumbers, &big[0], &big[0]+20, mask)<<std::endl;
    std::cout<<std::setw(width*2)<<"v.set(inumbers, cnumbers) : "                                                       <<v.set(inumbers, cnumbers)<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Vectorized::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // VECTORIZED_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
