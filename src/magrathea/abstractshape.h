/* ***************************** ABSTRACTSHAPE ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          AbstractShape
// DESCRIPTION :    Common abstraction of n-dimensional shapes 
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           abstractshape.h
/// \brief          Common abstraction of n-dimensional shapes
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABSTRACTSHAPE_H_INCLUDED
#define ABSTRACTSHAPE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <limits>
#include <cmath>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Common abstraction of n-dimensional shapes
/// \brief          Common abstraction of n-dimensional shapes.
/// \details        Provides compile-time helper functions to deal with 
///                 n-dimensional geometrical objects like hypercubes or 
///                 hyperspheres.
class AbstractShape
{
    // Protected lifecycle
    /// \name           Protected lifecycle
    //@{
    protected: 
        inline ~AbstractShape();
    //@}

    // Lifecycle
    /// \name           Lifecycle
    //@{
    protected: 
        template <class... Dummy> explicit inline AbstractShape(const Dummy&... dummy);
    //@}
    
    // Mathematical functions
    /// \name           Mathematical functions
    //@{
    public:
        template <int Exponent, typename Type> static constexpr Type pow(const Type& value);
        template <int Root, typename Type> static constexpr Type rt(const Type& value);
        template <int Value, typename Type = long long int> static constexpr Type factorial();
        template <int Value, bool Odd, typename Type = long long int> static constexpr Type factorial(); 
        template <int Set, int Subset, typename Type = long long int> static constexpr Type binomial();
        template <int Set, int Subset, bool Repetition, typename Type = long long int> static constexpr Type combinations();
        template <int Set, int Subset, bool Repetition, typename Type = long long int> static constexpr Type variations();
        template <long long int Value, long long int Even = 1, long long int Odd = -1, typename Type = long long int> static constexpr Type alt();
    //@}

    // Constants
    /// \name           Constants
    //@{
    public:
        template <int Exponent = 1, typename Type = double> static constexpr Type golden();  
        template <int Exponent = 1, typename Type = double> static constexpr Type pi();  
        template <int Exponent = 1, typename Type = double> static constexpr Type sqrtpi(); 
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
inline AbstractShape::~AbstractShape()
= default;
// -------------------------------------------------------------------------- //



//-------------------------------- LIFECYCLE -------------------------------- //
// Explicit generic constructor 
/// \brief          Explicit generic constructor. 
/// \details        Provides a constructor with dummy parameters for hackery
///                 purposes.
/// \tparam         Dummy (Dummy types.)
/// \param[in]      dummy Dummy arguments.
template <class... Dummy>
inline AbstractShape::AbstractShape(const Dummy&... dummy)
{
    (void)(sizeof...(dummy));
}
//--------------------------------------------------------------------------- //



// ------------------------- MATHEMATICAL FUNCTIONS ------------------------- //
// Integral exponentiation
/// \brief          Integral exponentiation.
/// \details        Computes the integral exponentiation of the specified value.
/// \tparam         Exponent Value of the exponent \f$n\f$.
/// \tparam         Type (Result type.)
/// \param[in]      value Value of the argument \f$x\f$.
/// \return         \f$x^{n}\f$.
template <int Exponent, typename Type> 
constexpr Type AbstractShape::pow(const Type& value)
{
    return (Exponent > 1) ? (value*pow<(Exponent-1)*(Exponent > 1), Type>(value)) : ((Exponent < 0) ? (Type(1)/pow<(-Exponent)*(Exponent < 0), Type>(value)) : ((Exponent == 1) ? (value) : (Type(1))));
}

// Floating-point root
/// \brief          Floating-point root.
/// \details        Computes the floating-point root of the specified value.
/// \tparam         Degree Value of the degree \f$n\f$.
/// \tparam         Type (Result type.)
/// \param[in]      value Value of the argument \f$x\f$.
/// \return         \f$\sqrt[n]{x}\f$.
template <int Degree, typename Type>
constexpr Type AbstractShape::rt(const Type& value)
{
    return (Degree > 0) ? (((Degree == 0) || (Degree > 3)) ? (std::pow(value, Type(1)/Type(Degree))) : ((Degree == 3) ? (std::cbrt(value)) : ((Degree == 2) ? (std::sqrt(value)) : (value)))) : ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? ((Degree < 0) ? (Type(1)/rt<(-Degree)*(Degree < 0), Type>(value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value)) : (value));
}

// Factorial
/// \brief          Factorial.
/// \details        Computes the factorial of the specified integer.
/// \tparam         Value Value of the argument \f$n\f$.
/// \tparam         Type (Result type.)
/// \return         \f$n!\f$.
template <int Value, typename Type> 
constexpr Type AbstractShape::factorial()
{
    static_assert(Value >= 0, "ERROR = AbstractShape::factorial() : underflow");
    static_assert((Value <= 1) || (factorial<(Value-1)*(Value > 1), Type>() < std::numeric_limits<Type>::max()/Type((Value != 0) ? (Value) : (1))), "ERROR = AbstractShape::factorial() : overflow");
    return (Value > 1) ? (Type(Value)*factorial<(Value-1)*(Value > 1), Type>()) : (Type(Value >= 0));
}

// Factorial
/// \brief          Factorial.
/// \details        Computes the double factorial of the specified integer:
///                 multiplies the number by all inferior numbers that are even 
///                 or odd depending on the parameter.
/// \tparam         Value Value of the argument \f$n\f$.
/// \tparam         Odd Odd double factorial if true, even otherwise.
/// \tparam         Type (Result type.)
/// \return         \f$n!!\f$.
template <int Value, bool Odd, typename Type> 
constexpr Type AbstractShape::factorial()
{
    static_assert(Value >= 0, "ERROR = AbstractShape::factorial() : underflow");
    static_assert((Value <= 1) || (Value%2 != Odd) || (factorial<(Value-2)*(Value > 2), Odd, Type>() < std::numeric_limits<Type>::max()/Type((Value != 0) ? (Value) : (1))), "ERROR = AbstractShape::factorial() : overflow");
    return (Value > 1) ? ((Value%2 == Odd) ? (Type(Value)*factorial<(Value-2)*(Value > 2), Odd, Type>()) : (factorial<(Value-1)*(Value > 1), Odd, Type>())) : (Type(Value >= 0));
} 

// Binomial
/// \brief          Binomial.
/// \details        Computes the binomial coefficient of the specified integers,
///                 including extension to negative numbers.
/// \tparam         Set Value of the set \f$n\f$.
/// \tparam         Subset Value of the subset \f$k\f$.
/// \tparam         Type (Result type.)
/// \return         \f$\binom{n}{k}\f$.
template <int Set, int Subset, typename Type> 
constexpr Type AbstractShape::binomial()
{
    static_assert((!std::is_unsigned<Type>::value) || ((Set >= 0) && (Subset >= 0)), "ERROR = AbstractShape::binomial() : sign error");
    ////////////////////////////////////static_assert((Subset <= 0) || (Subset > Set) || (((binomial<Set-((Subset > 0) && (Subset <= Set)), Subset-((Subset > 0) && (Subset <= Set)), Type>()/2)+(binomial<Set-((Subset > 0) && (Subset <= Set)), Subset, Type>()/2)) < (std::numeric_limits<Type>::max()/2)), "ERROR = AbstractShape::binomial() : overflow");
    return (Set >= 0) ? (((Subset >= 0) && (Subset <= Set)) ? ((Subset > 0) ? (binomial<Set-((Subset > 0) && (Subset <= Set)), Subset-((Subset > 0) && (Subset <= Set)), Type>()+binomial<Set-((Subset > 0) && (Subset <= Set)), Subset, Type>()) : (Type(1))) : (Type(0))) : ((Subset >= 0) ? (Type((Subset%2 == 0) ? (1) : (-1))*binomial<(-Set+Subset-1)*((Set < 0) && (Subset >= 0)), Subset*((Set < 0) && (Subset >= 0)), Type>()) : ((Subset <= Set) ? (Type(((Set-Subset)%2 == 0) ? (1) : (-1))*binomial<(-Subset-1)*((Subset <= Set) && (Set < 0)), (Set-Subset)*((Subset <= Set) && (Set < 0)), Type>()) : (Type(0))));
}

// Combinations
/// \brief          Combinations.
/// \details        Computes the number of combinations with or without 
///                 repetitions for the specified integers.
/// \tparam         Set Value of the set \f$n\f$.
/// \tparam         Subset Value of the subset \f$k\f$.
/// \tparam         Repetition Repetition specifier.
/// \tparam         Type (Result type.)
/// \return         \f$C^{k}_{n}\f$ or \f${}^{R}C^{k}_{n}\f$.
template <int Set, int Subset, bool Repetition, typename Type>
constexpr Type AbstractShape::combinations()
{
    static_assert((Subset >= 0) && (Set >= Subset), "ERROR = AbstractShape::combinations() : domain error");
    return ((Subset >= 0) && (Set >= Subset)) ? ((Repetition) ? (binomial<Set+Subset-((Set != 0) || (Subset != 0)), Subset, Type>()) : (binomial<Set, Subset, Type>())) : (Type(0));
}

// Variations
/// \brief          Variations.
/// \details        Computes the number of variations with or without 
///                 repetitions for the specified integers.
/// \tparam         Set Value of the set \f$n\f$.
/// \tparam         Subset Value of the subset \f$k\f$.
/// \tparam         Repetition Repetition specifier.
/// \tparam         Type (Result type.)
/// \return         \f$V^{k}_{n}\f$ or \f${}^{R}V^{k}_{n}\f$.
template <int Set, int Subset, bool Repetition, typename Type> 
constexpr Type AbstractShape::variations()
{
    static_assert((Subset >= 0) && (Set >= Subset), "ERROR = AbstractShape::variations() : domain error");
    static_assert((Repetition) ? (pow<Subset, long double>(Set) <= static_cast<long double>(std::numeric_limits<Type>::max())) : (factorial<Subset, Type>() <= std::numeric_limits<Type>::max()/((binomial<Set, Subset, Type>() != 0) ? (binomial<Set, Subset, Type>()) : (Type(1)))), "ERROR = AbstractShape::variations() : overflow");
    return ((Subset >= 0) && (Set >= Subset)) ? ((Repetition) ? (pow<Subset, Type>(Set)) : (factorial<Subset, Type>()*binomial<Set, Subset, Type>())) : (Type(0));
}

// Alternating sequence
/// \brief          Alternating sequence.
/// \details        Computes the specified value of an alternating sequence.
/// \tparam         Value Value of the argument \f$n\f$.
/// \tparam         Even Value if even \f$n_{0}\f$.
/// \tparam         Odd Value if odd \f$n_{1}\f$.
/// \tparam         Type (Result type.)
/// \return         \f$n_{0}\f$ if odd, \f$n_{1}\f$ if even.
template <long long int Value, long long int Even, long long int Odd, typename Type> 
constexpr Type AbstractShape::alt()
{
    static_assert((!std::is_unsigned<Type>::value) || ((Even >= 0) && (Odd >= 0)), "ERROR = AbstractShape::alt() : sign error");
    return Type((Value >= 0) ? ((Value%2 == 0) ? (Even) : (Odd)) : (((-Value)%2 == 0) ? (Even) : (Odd)));
}
// -------------------------------------------------------------------------- //



// -------------------------------- CONSTANTS ------------------------------- //
// Golden ratio
/// \brief          Golden ratio.
/// \details        Computes the value of \f$\varphi = \frac{1+\sqrt{5}}{2}\f$ 
///                 at the given power.
/// \tparam         Exponent Value of the exponent \f$n\f$.
/// \tparam         Type (Result type.)
/// \return         \f$varphi^{n}\f$.
template <int Exponent, typename Type>
constexpr Type AbstractShape::golden()
{
    return Type(pow<Exponent>(1.6180339887498948482045868343656381177203091798057628621354486L));
}

// Pi
/// \brief          Pi.
/// \details        Computes the value of \f$pi\f$ at the given power.
/// \tparam         Exponent Value of the exponent \f$n\f$.
/// \tparam         Type (Result type.)
/// \return         \f$\pi^{n}\f$.
template <int Exponent, typename Type>
constexpr Type AbstractShape::pi()
{
    return Type(pow<Exponent>(3.1415926535897932384626433832795028841971693993751058209749445L));
}

// Square root of pi
/// \brief          Square root of pi.
/// \details        Computes the value of the square root of \f$pi\f$ at the 
///                 given power.
/// \tparam         Exponent Value of the exponent \f$n\f$.
/// \tparam         Type (Result type.)
/// \return         \f$\left(\sqrt{\pi}\right)^{n}\f$.
template <int Exponent, typename Type>
constexpr Type AbstractShape::sqrtpi()
{
    return Type(pow<Exponent>(1.7724538509055160272981674833411451827975494561223871282138077L));
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of AbstractShape.
/// \return         0 if no error.
int AbstractShape::example()
{
    std::cout<<"BEGIN = AbstractShape::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = AbstractShape::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = AbstractShape::example()"<<std::endl;
    return 1;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABSTRACTSHAPE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
