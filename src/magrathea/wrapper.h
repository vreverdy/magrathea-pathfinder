/* ********************************* WRAPPER ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Wrapper
// DESCRIPTION :    Basic value wrapper with getter and setter
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           wrapper.h
/// \brief          Basic value wrapper with getter and setter
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef WRAPPER_H_INCLUDED
#define WRAPPER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <string>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Basic value wrapper with getter and setter
/// \brief          Basic value wrapper with getter and setter.
/// \details        Provides a class that can wrap a value or an object and 
///                 allows access by <tt>operator()</tt>. It can be used as a 
///                 public class member to avoid the writing of trivial getter 
///                 and setter.
/// \tparam         Type Wrapped type.
template <typename Type = double>
class Wrapper final
{
    // Lifecycle
    /// \name           Lifecycle    
    //@{
    public:
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> constexpr Wrapper(const OtherType& source = Type());
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> constexpr Wrapper(const Wrapper<OtherType>& source);
    //@}

    // Assignment 
    /// \name           Assignment    
    //@{
    public:
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Wrapper<Type>& operator=(const OtherType& source);
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Wrapper<Type>& operator=(const Wrapper<OtherType>& source);
    //@}
    
    // Operators
    /// \name           Operators
    //@{
    public:
        constexpr operator Type() const;
        constexpr const Type& operator()() const;
        inline Type& operator()();
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Type& operator()(const OtherType& source);
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Type& operator()(const Wrapper<OtherType>& source);
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
    public:
        Type _data;                                                             ///< Wrapped object. 
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit value constructor
/// \brief          Implicit value constructor.
/// \details        Implicitely constructs the wrapper from a value of a 
///                 convertible type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Value.
template <typename Type>
template <typename OtherType, class> 
constexpr Wrapper<Type>::Wrapper(const OtherType& source)
: _data(source)
{
    ;
}

// Implicit wrapper constructor
/// \brief          Implicit wrapper constructor.
/// \details        Implicitely constructs the wrapper from a wrapper of a 
///                 convertible type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Other wrapper.
template <typename Type>
template <typename OtherType, class> 
constexpr Wrapper<Type>::Wrapper(const Wrapper<OtherType>& source)
: _data(source())
{
    ;
}
// -------------------------------------------------------------------------- //



// ------------------------------- ASSIGNMENT -------------------------------- //
// Value assignment operator
/// \brief          Value assignment operator.
/// \details        Assigns the contents from a value of a convertible type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Value.
/// \return         Self reference.
template <typename Type>
template <typename OtherType, class> 
inline Wrapper<Type>& Wrapper<Type>::operator=(const OtherType& source)
{
    _data = source;
    return *this;
}

// Wrapped assignment operator
/// \brief          Wrapped assignment operator.
/// \details        Assigns the contents from a wrapper of a convertible type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Other wrapper.
/// \return         Self reference.
template <typename Type>
template <typename OtherType, class> 
inline Wrapper<Type>& Wrapper<Type>::operator=(const Wrapper<OtherType>& source)
{
    _data = source();
    return *this;
}
// -------------------------------------------------------------------------- //



// -------------------------------- OPERATORS ------------------------------- //
// Cast operator
/// \brief          Cast operator.
/// \details        Implicitely casts the wrapper to the wrapped type.
/// \return         Copy.
template <typename Type>
constexpr Wrapper<Type>::operator Type() const
{
    return _data;
}

// Immutable getter operator
/// \brief          Immutable getter operator.
/// \details        Returns a const reference to the wrapped object.
/// \return         Const reference.
template <typename Type>
constexpr const Type& Wrapper<Type>::operator()() const
{
    return _data;
}

// Getter operator
/// \brief          Getter operator.
/// \details        Returns a reference to the wrapped object.
/// \return         Reference.
template <typename Type>
inline Type& Wrapper<Type>::operator()()
{
    return _data;
}

// Value setter operator
/// \brief          Value setter operator.
/// \details        Sets the contents from a value of a convertible type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Value.
/// \return         Reference.
template <typename Type>
template <typename OtherType, class> 
inline Type& Wrapper<Type>::operator()(const OtherType& source)
{
    return (_data = source);
}

// Wrapped setter operator
/// \brief          Wrapped setter operator.
/// \details        Sets the contents from a wrapper of a convertible type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Other wrapper.
/// \return         Reference.
template <typename Type>
template <typename OtherType, class>
inline Type& Wrapper<Type>::operator()(const Wrapper<OtherType>& source)
{
    return (_data = source());
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Wrapper.
/// \return         0 if no error.
template <typename Type> 
int Wrapper<Type>::example()
{
    // Initialize
    std::cout<<"BEGIN = Wrapper::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    static const Wrapper<int> cst(42);
    Wrapper<int> i(42);
    Wrapper<double> d(42.42);
    Wrapper<std::string> s("forty-two");
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"Wrapper<int>() : "                            <<Wrapper<int>()<<std::endl;
    std::cout<<std::setw(width)<<"Wrapper<int>(42) : "                          <<Wrapper<int>(42)<<std::endl;
    std::cout<<std::setw(width)<<"Wrapper<int>(42.42) : "                       <<Wrapper<int>(42.42)<<std::endl;
    std::cout<<std::setw(width)<<"Wrapper<int>(cst) : "                         <<Wrapper<int>(cst)<<std::endl;
    std::cout<<std::setw(width)<<"Wrapper<int>(i) : "                           <<Wrapper<int>(i)<<std::endl;
    std::cout<<std::setw(width)<<"Wrapper<int>(d) : "                           <<Wrapper<int>(d)<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"i = -42 : "                                   <<(i = -42)<<std::endl;
    std::cout<<std::setw(width)<<"i = 42.42 : "                                 <<(i = 42.42)<<std::endl;
    std::cout<<std::setw(width)<<"i = d : "                                     <<(i = d)<<std::endl;
    
    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"int(cst) : "                                  <<int(cst)<<std::endl;
    std::cout<<std::setw(width)<<"double(cst) : "                               <<double(cst)<<std::endl;
    std::cout<<std::setw(width)<<"cst() : "                                     <<cst()<<std::endl;
    std::cout<<std::setw(width)<<"i() : "                                       <<i()<<std::endl;
    std::cout<<std::setw(width)<<"d() : "                                       <<d()<<std::endl;
    std::cout<<std::setw(width)<<"s() : "                                       <<s()<<std::endl;
    std::cout<<std::setw(width)<<"s().size() : "                                <<s().size()<<std::endl;
    std::cout<<std::setw(width)<<"i(4) : "                                      <<i(4)<<std::endl;
    std::cout<<std::setw(width)<<"d(i) : "                                      <<d(i)<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Wrapper::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // WRAPPER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
