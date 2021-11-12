/* ******************************** CONSTANT ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Constant
// DESCRIPTION :    Numerical constant with constexpr constructor
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           constant.h
/// \brief          Numerical constant with constexpr constructor
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef CONSTANT_H_INCLUDED
#define CONSTANT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <array>
#include <ratio>
#include <cmath>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Numerical constant with constexpr constructor
/// \brief          Numerical constant with constexpr constructor.
/// \details        Provides a class for integral and numerical constants based
///                 on the structure of <tt>std::integral_constant</tt> with 
///                 some extra mathematical functions. It accepts single or a
///                 set of constants internally stored in an array. As all is 
///                 done at compile-time, all members works with copies of 
///                 arithmetic values instead of references. Finally, one should 
///                 define a constant to the maximum precision and calls the 
///                 cast function to convert to a smaller precision.
/// \tparam         Type Numerical type of the constant, that should have a
///                 constexpr constructor.
/// \tparam         Size Size of the set of constants, which is equal to one
///                 per default for a single value.
template <typename Type = double, unsigned int Size = 1>
class Constant final
{
    // Setup
    static_assert(std::is_arithmetic<Type>::value, "ERROR = Constant : non arithmetic type");
    static_assert(Size > 0, "ERROR = Constant : null size");
    
    // Lifecycle
    /// \name           Lifecycle    
    //@{
    public:
        constexpr Constant(const Type source);
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != 1>::type> constexpr Constant(const Types... source);
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> constexpr Constant(const OtherType source);
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> constexpr Constant(const Constant<OtherType, Size> source);
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> constexpr Constant(const std::array<OtherType, Size> source);
        constexpr Constant(const std::array<Type, Size> source);
    //@}
    
    // Operators
    /// \name           Operators
    //@{
    public:
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> operator()(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> operator()(const Types... values) const;
        constexpr Type operator()() const;
        constexpr Type operator[](const unsigned int i) const;
        constexpr operator Type() const;
    //@}
    
    // Management
    /// \name           Management
    //@{
    public:
        constexpr unsigned int size() const;
        template <unsigned int OtherSize> constexpr const Constant<Type, OtherSize> resize() const;
        constexpr const Constant<Type, Size> nullify() const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> set(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> set(const Types... values) const;
        template <typename OtherType = Type> constexpr OtherType get(const unsigned int i = 0) const;
        template <typename OtherType = Type> constexpr OtherType value(const unsigned int i = 0) const;
        template <typename OtherType = Type> constexpr const std::array<OtherType, Size> data() const;
        template <typename OtherType> constexpr const Constant<OtherType, Size> cast() const;
        constexpr const Constant<Type, Size> copy() const;
    //@}
    
    // Unary operations
    /// \name           Unary operations
    //@{
    public:
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> inv(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> inv(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> opp(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> opp(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> abs(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> abs(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> sgn(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> sgn(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> sq(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> sq(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> cb(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> cb(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> isqrt(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> isqrt(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> icbrt(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> icbrt(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> ilog2(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> ilog2(const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> ilog10(const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> ilog10(const Types... values) const;
        template <class Ratio = std::ratio<1>, typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> si(const Types&... values) const;
        template <class Ratio = std::ratio<1>, typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> si(const Types... values) const;
    //@}
    
    // Arithmetic operations
    /// \name           Arithmetic operations
    //@{
    public:
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> add(const Type rhs, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> add(const Type, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> sub(const Type rhs, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> sub(const Type, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> mul(const Type rhs, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> mul(const Type, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> div(const Type rhs, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> div(const Type, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> mod(const Type rhs, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> mod(const Type, const Types... values) const;
    //@}

    // Mathematical functions
    /// \name           Mathematical functions
    //@{
    public:
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> ratio(const Type num, const Type den, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> ratio(const Type, const Type, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> pow(const int n, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> pow(const int, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> irt(const int n, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> irt(const int, const Types... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) != Size>::type> constexpr const Constant<Type, Size> ilog(const int b, const Types&... values) const;
        template <typename... Types, class = typename std::enable_if<sizeof...(Types) == Size>::type> constexpr const Constant<Type, Size> ilog(const int, const Types... values) const;
    //@}
    
    // Stream
    /// \name           Stream
    //@{
    public:
        template <typename SelfType, unsigned int SelfSize> friend std::ostream& operator<<(std::ostream& lhs, const Constant<SelfType, SelfSize>& rhs);
    //@}

    // Helpers
    /// \name           Helpers
    //@{
    public:
        template <typename OtherType, unsigned int OtherSize, class Array, typename... Types, class = typename std::enable_if<sizeof...(Types) != OtherSize>::type> static constexpr const std::array<OtherType, OtherSize> transmute(const Array source, const Types&... values);
        template <typename OtherType, unsigned int OtherSize, class Array, typename... Types, class = typename std::enable_if<sizeof...(Types) == OtherSize>::type> static constexpr const std::array<OtherType, OtherSize> transmute(const Array, const Types... values);
        template <typename OtherType, class = typename std::enable_if<std::is_arithmetic<OtherType>::value>::type> static constexpr OtherType metapow(const OtherType x, const int n);
        template <typename OtherType, typename IntegralType = std::true_type, class = typename std::enable_if<(std::is_arithmetic<OtherType>::value) && ((std::is_same<IntegralType, std::true_type>::value) || (std::is_integral<IntegralType>::value))>::type> static constexpr OtherType metairt(const OtherType x, const int n, const IntegralType low = IntegralType(), const IntegralType high = IntegralType(), const IntegralType mid = IntegralType());
        template <typename OtherType, class = typename std::enable_if<std::is_arithmetic<OtherType>::value>::type> static constexpr OtherType metailog(const OtherType x, const int b);
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
        const std::array<Type, Size> _data;                                     ///< Set of numerical constants. 
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit constexpr value constructor
/// \brief          Implicit constexpr value constructor.
/// \details        Implicitely constructs the constant from a single constant
///                 arithmetic value of the same type.
/// \param[in]      source Single constant arithmetic value.
template <typename Type, unsigned int Size>
constexpr Constant<Type, Size>::Constant(const Type source)
: _data({{source}})
{
    static_assert(1 == Size, "ERROR = Constant::Constant() : not compatible size");
}

// Implicit constexpr values constructor
/// \brief          Implicit constexpr values constructor.
/// \details        Implicitely constructs the constant from a set of constant 
///                 arithmetic values.
/// \tparam         Types (Variadic types.)
/// \param[in]      source Set of constant arithmetic values.
template <typename Type, unsigned int Size>
template <typename... Types, class>
constexpr Constant<Type, Size>::Constant(const Types... source)
: _data({{source...}})
{
    static_assert(sizeof...(source) == Size, "ERROR = Constant::Constant() : not compatible size");
}

// Implicit constexpr value conversion constructor
/// \brief          Implicit constexpr value conversion constructor.
/// \details        Implicitely constructs the constant from a single constant 
///                 arithmetic value of another type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Single constant arithmetic value.
template <typename Type, unsigned int Size>
template <typename OtherType, class>
constexpr Constant<Type, Size>::Constant(const OtherType source)
: _data(transmute<Type, Size>(std::array<OtherType, Size>({{source}})))
{
    static_assert(1 == Size, "ERROR = Constant::Constant() : not compatible size");
}

// Implicit constexpr constant conversion constructor
/// \brief          Implicit constexpr constant conversion constructor.
/// \details        Implicitely constructs the constant from a constant of 
///                 another type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Other constant.
template <typename Type, unsigned int Size>
template <typename OtherType, class>
constexpr Constant<Type, Size>::Constant(const Constant<OtherType, Size> source)
: _data(transmute<Type, Size>(source.data()))
{
    ;
}

// Implicit constexpr array conversion constructor
/// \brief          Implicit constexpr array conversion constructor.
/// \details        Implicitely constructs the constant from a constant array of 
///                 another type.
/// \tparam         OtherType (Other data type.)
/// \param[in]      source Array.
template <typename Type, unsigned int Size>
template <typename OtherType, class> 
constexpr Constant<Type, Size>::Constant(const std::array<OtherType, Size> source)
: _data(transmute<Type, Size>(source))
{
    ;
}

// Implicit constexpr array constructor
/// \brief          Implicit constexpr array constructor.
/// \details        Implicitely constructs the constant from a constant array of 
///                 the same type.
/// \param[in]      source Array.
template <typename Type, unsigned int Size>
constexpr Constant<Type, Size>::Constant(const std::array<Type, Size> source)
: _data(source)
{
    ;
}
// -------------------------------------------------------------------------- //



// -------------------------------- OPERATORS ------------------------------- //
// Incomplete setter operator
/// \brief          Incomplete setter operator.
/// \details        Starts recursive extraction of the missing values from the 
///                 current array.
/// \tparam         Types (Variadic types.)
/// \param[in]      values Set of constant arithmetic values.
/// \return         Copy.
template <typename Type, unsigned int Size>  
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::operator()(const Types&... values) const
{
    return operator()(values..., _data[sizeof...(values)]); 
}

// Setter operator
/// \brief          Setter operator.
/// \details        Sets the values.
/// \tparam         Types (Variadic types.)
/// \param[in]      values Set of constant arithmetic values.
/// \return         Copy.
template <typename Type, unsigned int Size>  
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::operator()(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Default getter operator
/// \brief          Default getter operator.
/// \details        Returns the first value.
/// \return         Value.
template <typename Type, unsigned int Size>  
constexpr Type Constant<Type, Size>::operator()() const
{
    return _data[0];
}

// Subscript operator
/// \brief          Subscript operator.
/// \details        Returns the specified value of the underlying array.
/// \param[in]      i Index of the value.
/// \return         Value.
template <typename Type, unsigned int Size>  
constexpr Type Constant<Type, Size>::operator[](const unsigned int i) const
{
    return _data[i];
}

// Value cast operator
/// \brief          Value cast operator.
/// \details        Implicitely returns the first value of the underlying array of
///                 constants. This should be used with caution in case of a set
///                 of constants as it only operates on the first element.
/// \return         Value.
template <typename Type, unsigned int Size>  
constexpr Constant<Type, Size>::operator Type() const
{
    return _data[0];
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Size
/// \brief          Size.
/// \details        Returns the size of the constant.
/// \return         Value of the size.
template <typename Type, unsigned int Size>  
constexpr unsigned int Constant<Type, Size>::size() const
{
    return Size;
}

// Resize
/// \brief          Resize.
/// \details        Returns a copy of the constant with a new size.
/// \return         Copy.
template <typename Type, unsigned int Size>  
template <unsigned int OtherSize> 
constexpr const Constant<Type, OtherSize> Constant<Type, Size>::resize() const
{
    return transmute<Type, OtherSize>(_data);
}

// Nullify
/// \brief          Nullify.
/// \details        Returns a copy of the constant filled with null values.
/// \return         Copy.
template <typename Type, unsigned int Size>  
constexpr const Constant<Type, Size> Constant<Type, Size>::nullify() const
{
    return std::array<Type, Size>({{Type()}});
}

// Incomplete value setter
/// \brief          Incomplete value setter.
/// \details        Starts recursive extraction of the missing values from the 
///                 current array.
/// \tparam         Types (Variadic types.)
/// \param[in]      values Set of constant arithmetic values.
/// \return         Copy.
template <typename Type, unsigned int Size>  
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::set(const Types&... values) const
{
    return set(values..., _data[sizeof...(values)]); 
}

// Value setter
/// \brief          Value setter.
/// \details        Sets the values.
/// \tparam         Types (Variadic types.)
/// \param[in]      values Set of constant arithmetic values.
/// \return         Copy.
template <typename Type, unsigned int Size>  
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::set(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Value getter
/// \brief          Value getter.
/// \details        Returns the specified value casted to the desired type.
/// \tparam         OtherType Other data type.
/// \param[in]      i Index of the value.
/// \return         Value.
template <typename Type, unsigned int Size>  
template <typename OtherType> 
constexpr OtherType Constant<Type, Size>::get(const unsigned int i) const
{
    return _data[i];
}

// Access a value
/// \brief          Access a value.
/// \details        Returns the specified value casted to the desired type.
/// \tparam         OtherType Other data type.
/// \param[in]      i Index of the value.
/// \return         Value.
template <typename Type, unsigned int Size>  
template <typename OtherType> 
constexpr OtherType Constant<Type, Size>::value(const unsigned int i) const
{
    return _data[i];
}

// Access data
/// \brief          Access data.
/// \details        Returns the internal constant array casted to the desired 
///                 type.
/// \tparam         OtherType Other data type.
/// \return         Copy of the array.
template <typename Type, unsigned int Size>  
template <typename OtherType> 
constexpr const std::array<OtherType, Size> Constant<Type, Size>::data() const
{
    return transmute<OtherType, Size>(_data);
}

// Cast
/// \brief          Cast.
/// \details        Casts to another constant type.
/// \tparam         OtherType Other data type.
/// \return         Copy.
template <typename Type, unsigned int Size>  
template <typename OtherType> 
constexpr const Constant<OtherType, Size> Constant<Type, Size>::cast() const
{
    return transmute<OtherType, Size>(_data);
}

// Copy
/// \brief          Copy.
/// \details        Copies the constant.
/// \return         Copy.
template <typename Type, unsigned int Size>  
constexpr const Constant<Type, Size> Constant<Type, Size>::copy() const
{
    return _data;
}
// -------------------------------------------------------------------------- //



// ---------------------------- UNARY OPERATIONS ---------------------------- //
// Start calculation of the inverse
/// \brief          Start calculation of the inverse.
/// \details        Starts recursive calculation of the inverse constant value :
///                 \f$\frac{1}{x}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::inv(const Types&... values) const
{
    return inv(values..., Type(1)/_data[sizeof...(values)]); 
}

// Stop calculation of the inverse
/// \brief          Stop calculation of the inverse.
/// \details        Stops recursive calculation of the inverse constant value :
///                 \f$\frac{1}{x}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::inv(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the opposite
/// \brief          Start calculation of the opposite.
/// \details        Starts recursive calculation of the opposite constant value :
///                 \f$-x\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::opp(const Types&... values) const
{
    return opp(values..., -_data[sizeof...(values)]); 
}

// Stop calculation of opposite
/// \brief          Stop calculation of opposite.
/// \details        Stops recursive calculation of the opposite constant value :
///                 \f$-x\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::opp(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the absolute value
/// \brief          Start calculation of the absolute value.
/// \details        Starts recursive calculation of the absolute constant value 
///                 : \f$\left|x\right|\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::abs(const Types&... values) const
{
    return abs(values..., ((_data[sizeof...(values)] == Type()) || (_data[sizeof...(values)] == -Type())) ? (Type()) : ((_data[sizeof...(values)] >= Type()) ? (_data[sizeof...(values)]) : (-_data[sizeof...(values)]))); 
}

// Stop calculation of the absolute value
/// \brief          Stop calculation of the absolute value.
/// \details        Stops recursive calculation of the absolute constant value 
///                 : \f$\left|x\right|\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::abs(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the signum function
/// \brief          Start calculation of the the signum function.
/// \details        Starts recursive calculation of the signum constant value  
///                 : \f$sgn\left(x = 0\right) = 0 \f$, \f$sgn\left(x < 0
///                 \right) = -1 \f$, \f$sgn\left(x > 0\right) = 1\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::sgn(const Types&... values) const
{
    return abs(values..., ((_data[sizeof...(values)] == Type()) || (_data[sizeof...(values)] == -Type())) ? (Type()) : ((_data[sizeof...(values)] >= Type()) ? (Type(1)) : (Type(-1)))); 
}

// Stop calculation of the signum function
/// \brief          Stop calculation of the the signum function.
/// \details        Stops recursive calculation of the signum constant value 
///                 : \f$sgn\left(x = 0\right) = 0 \f$, \f$sgn\left(x < 0
///                 \right) = -1 \f$, \f$sgn\left(x > 0\right) = 1\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::sgn(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the square
/// \brief          Start calculation of the square.
/// \details        Starts recursive calculation of the squared constant value  
///                 : \f$x^{2}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::sq(const Types&... values) const
{
    return sq(values..., metapow(_data[sizeof...(values)], 2));
}

// Stop calculation of the square
/// \brief          Stop calculation of the square.
/// \details        Stops recursive calculation of the squared constant value  
///                 : \f$x^{2}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::sq(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the cube
/// \brief          Start calculation of the cube.
/// \details        Starts recursive calculation of the cubed constant value  
///                 : \f$x^{3}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::cb(const Types&... values) const
{
    return cb(values..., metapow(_data[sizeof...(values)], 3));
}

// Stop calculation of the cube
/// \brief          Stop calculation of the cube.
/// \details        Stops recursive calculation of the cubed constant value  
///                 : \f$x^{3}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::cb(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the integer square root
/// \brief          Start calculation of the integer square root.
/// \details        Starts recursive calculation of the integer square root 
///                 constant value : \f$\left[\sqrt{\left[x\right]}\right]\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::isqrt(const Types&... values) const
{
    return isqrt(values..., metairt(_data[sizeof...(values)], 2));
}

// Stop calculation of the integer square root
/// \brief          Stop calculation of the integer square root.
/// \details        Stops recursive calculation of the integer square root 
///                 constant value : \f$\left[\sqrt{\left[x\right]}\right]\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::isqrt(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the integer cube root
/// \brief          Start calculation of the integer cube root.
/// \details        Starts recursive calculation of the integer cube root 
///                 constant value : \f$\left[\sqrt[3]{\left[x\right]}\right]
///                 \f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::icbrt(const Types&... values) const
{
    return icbrt(values..., metairt(_data[sizeof...(values)], 3));
}

// Stop calculation of the integer cube root
/// \brief          Stop calculation of the integer cube root.
/// \details        Stops recursive calculation of the integer cube root 
///                 constant value : \f$\left[\sqrt[3]{\left[x\right]}\right]
///                 \f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::icbrt(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the integer base 2 logarithm
/// \brief          Start calculation of the integer base 2 logarithm.
/// \details        Starts recursive calculation of the integer base 2 
///                 logarithm constant value : \f$\left[\log_{2}\left(\left[x
///                 \right]\right)\right]\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ilog2(const Types&... values) const
{
    return ilog2(values..., metailog(_data[sizeof...(values)], 2));
}

// Stop calculation of the integer base 2 logarithm
/// \brief          Stop calculation of the integer base 2 logarithm.
/// \details        Stops recursive calculation of the integer base 2 
///                 logarithm constant value : \f$\left[\log_{2}\left(\left[x
///                 \right]\right)\right]\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ilog2(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start calculation of the integer base 10 logarithm
/// \brief          Start calculation of the integer base 10 logarithm.
/// \details        Starts recursive calculation of the integer base 10 
///                 logarithm constant value : \f$\left[\log_{10}\left(\left[x
///                 \right]\right)\right]\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ilog10(const Types&... values) const
{
    return ilog10(values..., metailog(_data[sizeof...(values)], 10));
}

// Stop calculation of the integer base 10 logarithm
/// \brief          Stop calculation of the integer base 10 logarithm.
/// \details        Stops recursive calculation of the integer base 10 
///                 logarithm constant value : \f$\left[\log_{10}\left(\left[x
///                 \right]\right)\right]\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ilog10(const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start multiplication by a SI prefix
/// \brief          Start multiplication by a SI prefix.
/// \details        Starts recursive multiplication by a SI prefix ratio : 
///                 \f$\frac{x\times num}{den}\f$.
/// \tparam         Ratio Standard ratio type.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <class Ratio, typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::si(const Types&... values) const
{
    return si<Ratio>(values..., (_data[sizeof...(values)]*Ratio::num)/Ratio::den);
}

// Stop multiplication by a SI prefix
/// \brief          Stop multiplication by a SI prefix.
/// \details        Stops recursive multiplication by a SI prefix ratio : 
///                 \f$\frac{x\times num}{den}\f$.
/// \tparam         Ratio Standard ratio type.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <class Ratio, typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::si(const Types... values) const
{
    return Constant<Type, Size>(values...);
}
// -------------------------------------------------------------------------- //



// -------------------------- ARITHMETIC OPERATIONS ------------------------- //
// Start addition
/// \brief          Start addition.
/// \details        Starts recursive addition of a value : \f$x+y\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::add(const Type rhs, const Types&... values) const
{
    return add(rhs, values..., _data[sizeof...(values)]+rhs); 
}

// Stop addition
/// \brief          Stop addition.
/// \details        Stops recursive addition of a value : \f$x+y\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::add(const Type, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start substraction
/// \brief          Start substraction.
/// \details        Starts recursive substraction of a value : \f$x-y\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::sub(const Type rhs, const Types&... values) const
{
    return sub(rhs, values..., _data[sizeof...(values)]-rhs); 
}

// Stop substraction
/// \brief          Stop substraction.
/// \details        Stops recursive substraction of a value : \f$x-y\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::sub(const Type, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start multiplication
/// \brief          Start multiplication.
/// \details        Starts recursive multiplication by a value : \f$x \times y\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::mul(const Type rhs, const Types&... values) const
{
    return mul(rhs, values..., _data[sizeof...(values)]*rhs); 
}

// Stop multiplication
/// \brief          Stop multiplication.
/// \details        Stops recursive multiplication by a value : \f$x \times y\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::mul(const Type, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start division
/// \brief          Start division.
/// \details        Starts recursive division by a value : \f$\frac{x}{y}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::div(const Type rhs, const Types&... values) const
{
    return div(rhs, values..., _data[sizeof...(values)]/rhs); 
}

// Stop division
/// \brief          Stop division.
/// \details        Stops recursive division by a value : \f$\frac{x}{y}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::div(const Type, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start modulo calculation
/// \brief          Start modulo calculation.
/// \details        Starts recursive modulo calculation : \f$x\%y\f$. For 
///                 floating-point types, it computes the floating-point modulo.
/// \tparam         Types (Variadic types.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::mod(const Type rhs, const Types&... values) const
{
    return mod(rhs, values..., (std::is_floating_point<Type>::value) ? (_data[sizeof...(values)]-static_cast<long long int>(_data[sizeof...(values)]/rhs)*rhs) 
                                                                     : (static_cast<typename std::conditional<std::is_floating_point<Type>::value, long long int, Type>::type>(_data[sizeof...(values)])
                                                                       %static_cast<typename std::conditional<std::is_floating_point<Type>::value, long long int, Type>::type>(rhs))); 
}

// Stop modulo calculation
/// \brief          Stop modulo calculation.
/// \details        Stops recursive modulo calculation : \f$x\%y\f$. For 
///                 floating-point types, it computes the floating-point modulo.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::mod(const Type, const Types... values) const
{
    return Constant<Type, Size>(values...);
}
// -------------------------------------------------------------------------- //



// ------------------------- MATHEMATICAL FUNCTIONS ------------------------- //
// Start ratio multiplication
/// \brief          Start ratio multiplication.
/// \details        Starts recursive ratio multiplication : \f$\frac{x\times num
///                 }{den}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      num Numerator.
/// \param[in]      den Denominator.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ratio(const Type num, const Type den, const Types&... values) const
{
    return ratio(num, den, values..., (_data[sizeof...(values)]*num)/den);
}

// Stop ratio multiplication
/// \brief          Stop ratio multiplication.
/// \details        Stops recursive ratio multiplication : \f$\frac{x\times num
///                 }{den}\f$.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ratio(const Type, const Type, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start power calculation
/// \brief          Start power calculation.
/// \details        Starts recursive power calculation for integral exponent : 
///                 \f$x^{n}\f$. For negative exponents, the function returns 
///                 the inverse of the integral power.
/// \tparam         Types (Variadic types.)
/// \param[in]      n Exponent of the power.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::pow(const int n, const Types&... values) const
{
    return pow(n, values..., metapow(_data[sizeof...(values)], n));
}

// Stop power calculation
/// \brief          Stop power calculation.
/// \details        Stops recursive power calculation for integral exponent : 
///                 \f$x^{n}\f$. For negative exponents, the function returns 
///                 the inverse of the integral power.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::pow(const int, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start integer root calculation
/// \brief          Start integer part root calculation.
/// \details        Starts recursive integer root calculation for integral 
///                 arguments : \f$\left[\sqrt[n]{\left[x\right]}\right]\f$. 
///                 For non-compatible arguments, the function returns 0.
/// \tparam         Types (Variadic types.)
/// \param[in]      n Root order.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class>
constexpr const Constant<Type, Size> Constant<Type, Size>::irt(const int n, const Types&... values) const
{
    return irt(n, values..., metairt(_data[sizeof...(values)], n));
}

// Stop integer root calculation
/// \brief          Stop integer root calculation.
/// \details        Stops recursive integer root calculation for integral 
///                 arguments : \f$\left[\sqrt[n]{\left[x\right]}\right]\f$. 
///                 For non-compatible arguments, the function returns 0.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::irt(const int, const Types... values) const
{
    return Constant<Type, Size>(values...);
}

// Start integer logarithm calculation
/// \brief          Start integer logarithm calculation.
/// \details        Starts recursive integer logarithm calculation in the given 
///                 base for integral arguments : \f$\left[\log_{b}\left(\left[
///                 x\right]\right)\right]\f$. For non-compatible arguments, the 
///                 function returns 0.
/// \tparam         Types (Variadic types.)
/// \param[in]      b Logarithm integral base.
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ilog(const int b, const Types&... values) const
{
    return ilog(b, values..., metailog(_data[sizeof...(values)], b));
}

// Stop integer logarithm calculation
/// \brief          Stop integer logarithm calculation.
/// \details        Stops recursive integer logarithm calculation in the given 
///                 base for integral arguments : \f$\left[\log_{b}\left(\left[
///                 x\right]\right)\right]\f$. For non-compatible arguments, the 
///                 function returns 0.
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Copy.
template <typename Type, unsigned int Size> 
template <typename... Types, class> 
constexpr const Constant<Type, Size> Constant<Type, Size>::ilog(const int, const Types... values) const
{
    return Constant<Type, Size>(values...);
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Adds each element to the stream using the <tt>fill()</tt>
///                 character as a separator.
/// \tparam         SelfType (Numerical type of the constant.)
/// \tparam         SelfSize (Size of the set of constants.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side constant.
/// \return         Output stream.
template <typename SelfType, unsigned int SelfSize>
std::ostream& operator<<(std::ostream& lhs, const Constant<SelfType, SelfSize>& rhs)
{
    const char separator = lhs.fill();
    for (int i = 0; i < std::integral_constant<int, static_cast<int>(SelfSize)-1>::value; ++i) {
        lhs<<rhs._data[i]<<separator; 
    }
    return lhs<<rhs._data[std::integral_constant<int, static_cast<int>(SelfSize)-1>::value];
}
// -------------------------------------------------------------------------- //



// --------------------------------- HELPERS -------------------------------- //
// Start conversion of array data types
/// \brief          Start conversion of array data types.
/// \details        Recursively extracts and casts data from source to construct
///                 an array of another type.
/// \tparam         OtherType Other data type.
/// \tparam         OtherSize Other size.
/// \tparam         Array (Array type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      source Array to convert.
/// \param[in]      values (Recursively extracted values.)
/// \return         Casted array.
template <typename Type, unsigned int Size> 
template <typename OtherType, unsigned int OtherSize, class Array, typename... Types, class>
constexpr const std::array<OtherType, OtherSize> Constant<Type, Size>::transmute(const Array source, const Types&... values)
{
    return transmute<OtherType, OtherSize>(source, values..., OtherType(source[sizeof...(values)]));
}

// Stop conversion of array data types
/// \brief          Stop conversion of array data types.
/// \details        Stops the recursive extraction and constructs the final 
///                 array from data.
/// \tparam         OtherType Other data type.
/// \tparam         OtherSize Other size.
/// \tparam         Array (Array type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      values (Recursively extracted values.)
/// \return         Casted array.
template <typename Type, unsigned int Size> 
template <typename OtherType, unsigned int OtherSize, class Array, typename... Types, class>
constexpr const std::array<OtherType, OtherSize> Constant<Type, Size>::transmute(const Array, const Types... values)
{
    return std::array<OtherType, OtherSize>({{values...}});
}

// Compile-time exponentiation
/// \brief          Compile-time exponentiation.
/// \details        Returns the value of the integral power of the argument : 
///                 \f$x^{n}\f$. For negative exponents, the function returns 
///                 the inverse of the integral power.
/// \tparam         OtherType (Other data type.)
/// \param[in]      x Value.
/// \param[in]      n Exponent of the power.
/// \return         Value of the result.
template <typename Type, unsigned int Size> 
template <typename OtherType, class>
constexpr OtherType Constant<Type, Size>::metapow(const OtherType x, const int n)
{
    return (n != 0) ? ((n > 0) ? (x*metapow(x, n-1)) : (OtherType(1)/metapow(x, -n))) : (OtherType(1));
}

// Compile-time integer root
/// \brief          Compile-time integer root.
/// \details        Returns the value of the integer root for integral arguments 
///                 at compile-time : \f$\left[\sqrt[n]{\left[x\right]}\right]
///                 \f$. For non-compatible arguments, the function returns 0.
/// \tparam         OtherType (Other data type.)
/// \tparam         IntegralType (Integer type.)
/// \param[in]      x Value.
/// \param[in]      n Root order.
/// \param[in]      low (Lower bound for binary recursive search.)
/// \param[in]      high (Higher bound for binary recursive search.)
/// \param[in]      mid (Mid value for binary recursive search.)
/// \return         Value of the result.
template <typename Type, unsigned int Size> 
template <typename OtherType, typename IntegralType, class> 
constexpr OtherType Constant<Type, Size>::metairt(const OtherType x, const int n, const IntegralType low, const IntegralType high, const IntegralType mid)
{
    return (std::is_same<IntegralType, std::true_type>::value) ? ((static_cast<long long int>(x) && n && ((static_cast<long long int>(x) > 0) || (2*(n/2) != n)))
                                                               ? ((n > 0)
                                                               ? ((static_cast<long long int>(x) > 0)
                                                               ? (OtherType(metairt<long long int, long long int>(x, n, 1, x, 1+x/2)))
                                                               : (-metairt(-x, n))) 
                                                               : (OtherType(1)/metairt(x, -n)))
                                                               : (OtherType()))
                                                               : ((low != high) 
                                                               ? ((x < metapow(static_cast<double>(mid), n)) 
                                                               ? (metairt<long long int, long long int>(x, n, low, mid-1, (low+mid)/2)) 
                                                               : (metairt<long long int, long long int>(x, n, mid, high, (mid+high+1)/2))) 
                                                               : (mid));
}

// Compile-time integer logarithm
/// \brief          Compile-time integer logarithm.
/// \details        Returns the value of the integer logarithm in the given base
///                 for integral arguments at compile-time : \f$\left[\log_{b}
///                 \left(\left[x\right]\right)\right]\f$. For non-compatible 
///                 arguments, the function returns 0.
/// \tparam         OtherType (Other data type.)
/// \param[in]      x Value.
/// \param[in]      b Logarithm integral base.
/// \return         Value of the result.
template <typename Type, unsigned int Size> 
template <typename OtherType, class>
constexpr OtherType Constant<Type, Size>::metailog(const OtherType x, const int b)
{
    return ((static_cast<long long int>(x) >= b) && (b > 1)) ? (OtherType(1+metailog(static_cast<long long int>(x)/b, b))) : (OtherType());
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Constant.
/// \return         0 if no error.
template <typename Type, unsigned int Size> 
int Constant<Type, Size>::example()
{
    // Initialize
    std::cout<<"BEGIN = Constant::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    const std::array<double, 2> piarray({{3.1415926, 6.2831853}});
    
    // Construction
    Constant<int> icst(42);
    Constant<double> pi(3.1415926);
    Constant<double, 2> pi2(3.1415926, 6.2831853);
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42) : "                         <<Constant<int>(42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(3.1415926) : "               <<Constant<double>(3.1415926)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double, 2>(3.14159, 6.2831) : "      <<Constant<double, 2>(3.14159, 6.2831)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42.) : "                        <<Constant<int>(42.)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42) : "                      <<Constant<double>(42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(pi) : "                      <<Constant<double>(pi)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double, 2>(pi2) : "                  <<Constant<double, 2>(pi2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(icst) : "                    <<Constant<double>(icst)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double, 2>(piarray) : "              <<Constant<double, 2>(piarray)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int, 2>(piarray) : "                 <<Constant<int, 2>(piarray)<<std::endl;
    
    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"pi2(1.)[0] : "                                <<pi2(1.)[0]<<std::endl;
    std::cout<<std::setw(width)<<"pi2(1.)[1] : "                                <<pi2(1.)[1]<<std::endl;
    std::cout<<std::setw(width)<<"pi2(1., 2.)[0] : "                            <<pi2(1., 2.)[0]<<std::endl;
    std::cout<<std::setw(width)<<"pi2(1., 2.)[1] : "                            <<pi2(1., 2.)[1]<<std::endl;
    std::cout<<std::setw(width)<<"pi() : "                                      <<pi()<<std::endl;
    std::cout<<std::setw(width)<<"pi2[0] : "                                    <<pi2[0]<<std::endl;
    std::cout<<std::setw(width)<<"pi2[1] : "                                    <<pi2[1]<<std::endl;
    std::cout<<std::setw(width)<<"double(pi2) : "                               <<double(pi2)<<std::endl;
    std::cout<<std::setw(width)<<"int(pi2) : "                                  <<int(pi2)<<std::endl;
    std::cout<<std::setw(width)<<"static_cast<long double>(pi2) : "             <<static_cast<long double>(pi2)<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"pi2.size() : "                                <<pi2.size()<<std::endl;
    std::cout<<std::setw(width)<<"pi2.resize<1>() : "                           <<pi2.resize<1>()<<std::endl;
    std::cout<<std::setw(width)<<"pi2.nullify() : "                             <<pi2.nullify()<<std::endl;
    std::cout<<std::setw(width)<<"pi2.set(1.)[0] : "                            <<pi2.set(1.)[0]<<std::endl;
    std::cout<<std::setw(width)<<"pi2.set(1.)[1] : "                            <<pi2.set(1.)[1]<<std::endl;
    std::cout<<std::setw(width)<<"pi2.set(1., 2.)[0] : "                        <<pi2.set(1., 2.)[0]<<std::endl;
    std::cout<<std::setw(width)<<"pi2.set(1., 2.)[1] : "                        <<pi2.set(1., 2.)[1]<<std::endl;
    std::cout<<std::setw(width)<<"pi2.get() : "                                 <<pi2.get()<<std::endl;
    std::cout<<std::setw(width)<<"pi2.get(1) : "                                <<pi2.get(1)<<std::endl;
    std::cout<<std::setw(width)<<"pi2.get<int>(1) : "                           <<pi2.get<int>(1)<<std::endl;
    std::cout<<std::setw(width)<<"pi2.value() : "                               <<pi2.value()<<std::endl;
    std::cout<<std::setw(width)<<"pi2.value(1) : "                              <<pi2.value(1)<<std::endl;
    std::cout<<std::setw(width)<<"pi2.value<int>(1) : "                         <<pi2.value<int>(1)<<std::endl;
    std::cout<<std::setw(width)<<"pi2.data()[1] : "                             <<pi2.data()[1]<<std::endl;
    std::cout<<std::setw(width)<<"pi2.data<int>()[1] : "                        <<pi2.data<int>()[1]<<std::endl;
    std::cout<<std::setw(width)<<"pi2.cast<int>() : "                           <<pi2.cast<int>()<<std::endl;
    std::cout<<std::setw(width)<<"pi2.copy() : "                                <<pi2.copy()<<std::endl;

    // Unary operations
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Unary operations : "                          <<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).inv() : "             <<Constant<double>(42.42).inv()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).inv() : "                   <<Constant<int>(42.42).inv()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).opp() : "             <<Constant<double>(42.42).opp()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42.42).opp() : "            <<Constant<double>(-42.42).opp()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).abs() : "             <<Constant<double>(42.42).abs()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(0.).abs() : "                <<Constant<double>(0.).abs()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-0.).abs() : "               <<Constant<double>(-0.).abs()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42.42).abs() : "            <<Constant<double>(-42.42).abs()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).sgn() : "             <<Constant<double>(42.42).sgn()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(0.).sgn() : "                <<Constant<double>(0.).sgn()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-0.).sgn() : "               <<Constant<double>(-0.).sgn()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42.42).sgn() : "            <<Constant<double>(-42.42).sgn()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).sq() : "              <<Constant<double>(42.42).sq()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).cb() : "              <<Constant<double>(42.42).cb()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).isqrt() : "           <<Constant<double>(42.42).isqrt()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).icbrt() : "           <<Constant<double>(42.42).icbrt()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).ilog2() : "           <<Constant<double>(42.42).ilog2()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).ilog10() : "          <<Constant<double>(42.42).ilog10()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).si<std::mega>() : "      <<Constant<double>(42).si<std::mega>()<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).si<std::micro>() : "     <<Constant<double>(42).si<std::micro>()<<std::endl;

    // Arithmetic operations
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Arithmetic operations : "                     <<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).add(0.42) : "         <<Constant<double>(42.42).add(0.42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).sub(0.42) : "         <<Constant<double>(42.42).sub(0.42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).mul(0.42) : "         <<Constant<double>(42.42).mul(0.42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).div(0.42) : "         <<Constant<double>(42.42).div(0.42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).mod(0.42) : "         <<Constant<double>(42.42).mod(0.42)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42.42).mod(21) : "           <<Constant<double>(42.42).mod(21)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).mod(20) : "                 <<Constant<int>(42).mod(20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).mod(-20) : "                <<Constant<int>(42).mod(-20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(-42).mod(20) : "                <<Constant<int>(-42).mod(20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(-42).mod(-20) : "               <<Constant<int>(-42).mod(-20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).mod(20) : "              <<Constant<double>(42).mod(20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).mod(-20) : "             <<Constant<double>(42).mod(-20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42).mod(20) : "             <<Constant<double>(-42).mod(20)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42).mod(-20) : "            <<Constant<double>(-42).mod(-20)<<std::endl;

    // Mathematical functions
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Mathematical functions : "                    <<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).ratio(3, 2) : "             <<Constant<int>(42).ratio(3, 2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).pow(2) : "                  <<Constant<int>(42).pow(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).pow(0) : "                  <<Constant<int>(42).pow(0)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<int>(42).pow(-2) : "                 <<Constant<int>(42).pow(-2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).pow(2) : "               <<Constant<double>(42).pow(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).pow(0) : "               <<Constant<double>(42).pow(0)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).pow(-2) : "              <<Constant<double>(42).pow(-2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).irt(2) : "               <<Constant<double>(42).irt(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(42).irt(0) : "               <<Constant<double>(42).irt(0)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42).irt(0) : "              <<Constant<double>(-42).irt(0)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42).irt(1) : "              <<Constant<double>(-42).irt(1)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42).irt(2) : "              <<Constant<double>(-42).irt(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-42).irt(3) : "              <<Constant<double>(-42).irt(3)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(-1).ilog(2) : "              <<Constant<double>(-1).ilog(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(0).ilog(2) : "               <<Constant<double>(0).ilog(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(1).ilog(2) : "               <<Constant<double>(1).ilog(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(2).ilog(2) : "               <<Constant<double>(2).ilog(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(3).ilog(2) : "               <<Constant<double>(3).ilog(2)<<std::endl;
    std::cout<<std::setw(width)<<"Constant<double>(4).ilog(2) : "               <<Constant<double>(4).ilog(2)<<std::endl;

    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, pi2) : "                <<pi2<<std::endl;

    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Helpers : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"icst.transmute<int, 2>(piarray)[0] : "        <<icst.transmute<int, 2>(piarray)[0]<<std::endl;
    std::cout<<std::setw(width)<<"icst.metapow(42.42, 2) : "                    <<icst.metapow(42.42, 2)<<std::endl;
    std::cout<<std::setw(width)<<"icst.metairt(42.42, 2) : "                    <<icst.metairt(42.42, 2)<<std::endl;
    std::cout<<std::setw(width)<<"icst.metailog(42.42, 2) : "                   <<icst.metailog(42.42, 2)<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Constant::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // CONSTANT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
