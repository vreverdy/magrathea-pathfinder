/* ***************************** ABSTRACTNARRAY ***************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          AbstractNArray
// DESCRIPTION :    Abstract base class of n-dimensional mathematical arrays
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           abstractnarray.h
/// \brief          Abstract base class of n-dimensional mathematical arrays
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABSTRACTNARRAY_H_INCLUDED
#define ABSTRACTNARRAY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <stdexcept>
#include <initializer_list>
#include <array>
#include <iterator>
#include <random>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <functional>
#include <utility>
// Include libs
// Include project
#include "staticvectorizer.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Abstract base class of n-dimensional mathematical arrays
/// \brief          Abstract base class of n-dimensional mathematical arrays.
/// \details        Provides common base for n-dimensional order 1 mathematical 
///                 containers thanks to the curiously recurring template 
///                 pattern (CRTP) trick. The class derives from 
///                 StaticVectorizer, provides the storage and adds functions 
///                 specific to arithmetic types. To use it, one has to derive 
///                 from this class and pass the derived class itself as the 
///                 CRTP parameter. The derived classes have to satisfy the 
///                 conditions required by the Vectorizer base class.
/// \tparam         Kind Kind of arguments.
/// \tparam         Size Number of elements.
/// \tparam         Crtp Derived CRTP class.
/// \tparam         Type Arithmetic data type.
/// \tparam         Parameters List of parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
class AbstractNArray
: public StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>
{
    // Setup
    static_assert(std::is_arithmetic<Type>::value, "ERROR = AbstractNArray : non arithmetic type");
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator=;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::min;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::max;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::null;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::eq;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::ne;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::modify;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::apply;
    public: using StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::reduce;

    // Protected lifecycle
    /// \name           Protected lifecycle
    //@{
    protected: 
        inline ~AbstractNArray();
    //@}
    
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        inline AbstractNArray();
        template <typename FundamentalType = Type, class = typename std::enable_if<std::is_fundamental<FundamentalType>::value>::type> inline AbstractNArray(const AbstractNArray<Kind, Size, Crtp, FundamentalType, Parameters...>& source);
        template <typename OtherType = Type, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline AbstractNArray(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc, class = typename std::enable_if<sizeof...(Misc) != 0>::type> explicit inline AbstractNArray(const Misc&... misc);
    //@}
    
    // Vectorization
    /// \name           Vectorization
    //@{
    public:
        inline Type& operator[](const unsigned int i);
        inline const Type& operator[](const unsigned int i) const;
    //@}

    // Access
    /// \name           Access
    //@{
    public:
        inline Type* data();
        inline const Type* data() const;
    //@}
    
    // Iterators
    /// \name           Iterators
    //@{
    public:
        inline Type* begin();
        inline const Type* begin() const;
        inline const Type* cbegin() const;
        inline Type* end();
        inline const Type* end() const;
        inline const Type* cend() const;
        inline std::reverse_iterator<Type*> rbegin();
        inline std::reverse_iterator<const Type*> rbegin() const;
        inline std::reverse_iterator<const Type*> crbegin() const;
        inline std::reverse_iterator<Type*> rend();
        inline std::reverse_iterator<const Type*> rend() const;
        inline std::reverse_iterator<const Type*> crend() const;
        template <typename IteratorType> inline unsigned int index(const IteratorType& it, typename std::iterator_traits<IteratorType>::iterator_category* = nullptr) const;
    //@}

    // Comparison
    /// \name           Comparison
    //@{
    public:
        inline bool null(const Type& tolerance) const;
        template <class GenericType> inline bool eq(const GenericType& rhs, const Type& tolerance) const; 
        template <class GenericType> inline bool ne(const GenericType& rhs, const Type& tolerance) const; 
    //@}

    // Statistics
    /// \name           Statistics
    //@{
    public:
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline const Type& amin(const Mask& bitmask = Mask()) const;
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline const Type& amax(const Mask& bitmask = Mask()) const;
        template <typename Return = typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type, class Coefficient = Type, class Mask = std::true_type, class = typename std::enable_if<(std::is_arithmetic<Return>::value) && ((std::is_base_of<Vectorizer, Coefficient>::value) || (std::is_convertible<Coefficient, Type>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Return mean(const Coefficient& coefficient = Coefficient(1), const Mask& bitmask = Mask()) const;
        template <typename Return = typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type, typename Correction = Return, class Coefficient = Type, class Mask = std::true_type, class = typename std::enable_if<(std::is_arithmetic<Return>::value) && (std::is_arithmetic<Correction>::value) && ((std::is_base_of<Vectorizer, Coefficient>::value) || (std::is_convertible<Coefficient, Type>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Return sigma(const Correction& correction = Correction(0), const Coefficient& coefficient = Coefficient(1), const Mask& bitmask = Mask()) const;
    //@}
    
    // Application
    /// \name           Application
    //@{
    public:
        template <typename Return = Type, typename First = Type, class... Types, class... Args, class = typename std::enable_if<(std::is_convertible<typename std::decay<First>::type, typename std::decay<Type>::type>::value) && (sizeof...(Types) == sizeof...(Args))>::type> inline Crtp<Type, Parameters...>& modify(Return (*f)(First, Types...), Args&&... args);
        template <typename Return = Type, typename First = Type, class... Types, class... Args, class Mask, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) && (std::is_convertible<typename std::decay<First>::type, typename std::decay<Type>::type>::value) && (sizeof...(Types) == sizeof...(Args))>::type> inline Crtp<Type, Parameters...>& modify(const Mask& bitmask, Return (*f)(First, Types...), Args&&... args);
        template <typename Return = Type, typename First = Type, class... Types, class... Args, class = typename std::enable_if<(std::is_convertible<typename std::decay<First>::type, typename std::decay<Type>::type>::value) && (sizeof...(Types) == sizeof...(Args))>::type> inline Crtp<Return, Parameters...> apply(Return (*f)(First, Types...), Args&&... args) const;
        template <typename Return = Type, typename First = Type, class... Types, class... Args, class Mask, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) && (std::is_convertible<typename std::decay<First>::type, typename std::decay<Type>::type>::value) && (sizeof...(Types) == sizeof...(Args))>::type> inline Crtp<Return, Parameters...> apply(const Mask& bitmask, Return (*f)(First, Types...), Args&&... args) const;
    //@}

    // Count
    /// \name           Count
    //@{
    public:
        template <class Function = std::equal_to<Type>, class Mask = std::true_type, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool unicity(Function&& f = Function(), const Mask& bitmask = Mask()) const;
        template <class Function = std::equal_to<Type>, class Mask = std::true_type, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline unsigned int distinct(Function&& f = Function(), const Mask& bitmask = Mask()) const;
    //@}
    
    // Sort
    /// \name           Sort
    //@{
    public:
        template <class Function = std::less<Type>, class Indexes = Type, class Mask = std::true_type, class Pair = std::pair<Type, typename std::remove_const<typename std::remove_reference<decltype(Vectorizer::get(std::declval<Indexes>(), 0))>::type>::type>, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Crtp<Type, Parameters...>& arrange(Function&& f = Function(), Indexes&& indexes = Indexes(), const Mask& bitmask = Mask());        
        template <class Function = std::less<Type>, class Indexes = Type, class Mask = std::true_type, class Pair = std::pair<Type, typename std::remove_const<typename std::remove_reference<decltype(Vectorizer::get(std::declval<Indexes>(), 0))>::type>::type>, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Crtp<Type, Parameters...> sort(Function&& f = Function(), Indexes&& indexes = Indexes(), const Mask& bitmask = Mask()) const;    
        template <class Function = std::less<Type>, class Indexes = Type, class Mask = std::true_type, class Pair = std::pair<Type, typename std::remove_const<typename std::remove_reference<decltype(Vectorizer::get(std::declval<Indexes>(), 0))>::type>::type>, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Crtp<Type, Parameters...>& aarrange(Function&& f = Function(), Indexes&& indexes = Indexes(), const Mask& bitmask = Mask());   
        template <class Function = std::less<Type>, class Indexes = Type, class Mask = std::true_type, class Pair = std::pair<Type, typename std::remove_const<typename std::remove_reference<decltype(Vectorizer::get(std::declval<Indexes>(), 0))>::type>::type>, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Crtp<Type, Parameters...> asort(Function&& f = Function(), Indexes&& indexes = Indexes(), const Mask& bitmask = Mask()) const;   
        template <class Indexes = Type, class Pair = std::pair<typename std::remove_const<typename std::remove_reference<decltype(Vectorizer::get(std::declval<Indexes>(), 0))>::type>::type, unsigned int> > inline Crtp<Type, Parameters...>& rearrange(const Indexes& indexes = Indexes());
        template <class Indexes = Type, class Pair = std::pair<typename std::remove_const<typename std::remove_reference<decltype(Vectorizer::get(std::declval<Indexes>(), 0))>::type>::type, unsigned int> > inline Crtp<Type, Parameters...> resort(const Indexes& indexes = Indexes()) const;
    //@}
    
    // Mathematical functions
    /// \name           Mathematical functions
    //@{
    public:
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline Crtp<Type, Parameters...> abs(const Mask& bitmask = Mask()) const;
        template <int Exponent = 1> inline Crtp<Type, Parameters...> pow() const;
        template <class GenericType, class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline Crtp<Type, Parameters...> pow(const GenericType& rhs, const Mask& bitmask = Mask()) const;
        template <int Degree = 1> inline Crtp<Type, Parameters...> rt() const;
        template <class GenericType, class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline Crtp<Type, Parameters...> rt(const GenericType& rhs, const Mask& bitmask = Mask()) const;
        template <unsigned int Base = 0> inline Crtp<Type, Parameters...> log() const;
        template <class GenericType, class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline Crtp<Type, Parameters...> log(const GenericType& rhs, const Mask& bitmask = Mask()) const;
    //@}

    // Norm
    /// \name           Norm
    //@{
    public:
        template <unsigned int Degree = 2, typename NormType = typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type, class Mask = std::true_type, class = typename std::enable_if<(std::is_arithmetic<NormType>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline NormType norm(const Mask& bitmask = Mask()) const;
        template <unsigned int Degree = 2, typename NormType = typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type, class Mask = std::true_type, class = typename std::enable_if<(std::is_arithmetic<NormType>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Crtp<Type, Parameters...>& renormalize(const Mask& bitmask = Mask());
        template <unsigned int Degree = 2, typename NormType = typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type, class Mask = std::true_type, class = typename std::enable_if<(std::is_arithmetic<NormType>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline Crtp<Type, Parameters...> normalize(const Mask& bitmask = Mask()) const;
    //@}
    
    // Predefined
    /// \name           Predefined
    //@{
    public:
        static inline Crtp<Type, Parameters...> zero();
        static inline Crtp<Type, Parameters...> one();
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> static inline Crtp<Type, Parameters...> value(const Type& source = Type(), const Mask& bitmask = Mask());
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> static inline Crtp<Type, Parameters...> random(const Type& lowest = 0, const Type& highest = 1, const Mask& bitmask = Mask());
        template <class Engine, class Distribution, class Mask = std::true_type, class = typename std::enable_if<(std::decay<Engine>::type::min() != std::decay<Engine>::type::max()) && (!std::is_void<typename std::decay<Distribution>::type::result_type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> static inline Crtp<Type, Parameters...> random(Engine&& engine, Distribution&& distribution, const Mask& bitmask = Mask());
        template <class Function = unsigned int&& (*)(unsigned int&), class Mask = std::true_type, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> static inline Crtp<Type, Parameters...> indexed(Function&& f = std::forward<unsigned int>, const Mask& bitmask = Mask());
        template <class Function = std::plus<Type>, class Mask = std::true_type, class = typename std::enable_if<(!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> static inline Crtp<Type, Parameters...> progressive(const Type& first = Type(0), const Type& step = Type(1), Function&& f = Function(), const Mask& bitmask = Mask());
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



// -------------------------- PROTECTED LIFECYCLE --------------------------- //
// Protected destructor
/// \brief          Protected destructor.
/// \details        Does nothing.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::~AbstractNArray()
= default;
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit empty constructor
/// \brief          Implicit empty constructor.
/// \details        Does nothing.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::AbstractNArray()
: StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>()
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
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename FundamentalType, class> 
inline AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::AbstractNArray(const AbstractNArray<Kind, Size, Crtp, FundamentalType, Parameters...>& source)
: StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>()
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
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class... Misc, class>
inline AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::AbstractNArray(const std::initializer_list<OtherType>& source, const Misc&... misc)
: StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>()
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
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Misc, class> 
inline AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::AbstractNArray(const Misc&... misc)
: StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>()
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
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::operator[](const unsigned int i)
{
    return _data[i];
}

// Immutable direct access to the element
/// \brief          Immutable direct access to the element.
/// \details        Provides a constant direct access to the specified element.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::operator[](const unsigned int i) const
{
    return _data[i];
}
// -------------------------------------------------------------------------- //



// --------------------------------- ACCESS --------------------------------- //
// Direct access to the underlying storage 
/// \brief          Direct access to the underlying storage.
/// \details        Returns a pointer to the first element of the underlying 
///                 array.
/// \return         Pointer to the first element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::data()
{
    return &Vectorizer::get(*this, 0);
}

// Immutable direct access to the underlying storage 
/// \brief          Immutable direct access to the underlying storage.
/// \details        Returns a pointer to the first element of the underlying 
///                 array.
/// \return         Immutable pointer to the first element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::data() const
{
    return &Vectorizer::get(*this, 0);
}
// -------------------------------------------------------------------------- //



// ------------------------------- ITERATORS -------------------------------- //
// Iterator to the beginning
/// \brief          Iterator to the beginning.
/// \details        Returns a pointer to the first element.
/// \return         Pointer to the beginning.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::begin()
{
    return &Vectorizer::get(*this, 0);
}

// Immutable iterator to the beginning
/// \brief          Immutable iterator to the beginning.
/// \details        Returns a pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::begin() const
{
    return &Vectorizer::get(*this, 0);
}

// Forced immutable iterator to the beginning
/// \brief          Forced immutable iterator to the beginning.
/// \details        Returns a constant pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::cbegin() const
{
    return &Vectorizer::get(*this, 0);
}

// Iterator to the end
/// \brief          Iterator to the end.
/// \details        Returns a pointer to the position after the last element.
/// \return         Pointer to the end.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::end()
{
    return &Vectorizer::get(*this, 0)+Size;
}

// Immutable iterator to the end
/// \brief          Immutable iterator to the end.
/// \details        Returns a pointer to the position after the last element.
/// \return         Immutable pointer to the end.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::end() const
{
    return &Vectorizer::get(*this, 0)+Size;
}

// Forced immutable iterator to the end
/// \brief          Forced immutable iterator to the end.
/// \details        Returns a constant pointer to the position after the last 
///                 element.
/// \return         Immutable pointer to the end.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type* AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::cend() const
{
    return &Vectorizer::get(*this, 0)+Size;
}

// Reverse iterator to the beginning
/// \brief          Reverse iterator to the beginning.
/// \details        Returns a reversed pointer to the position after the last 
///                 element.
/// \return         Pointer to the end.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline std::reverse_iterator<Type*> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rbegin()
{
    return std::reverse_iterator<Type*>(&Vectorizer::get(*this, 0)+Size);
}

// Immutable reverse iterator to the beginning
/// \brief          Immutable reverse iterator to the beginning.
/// \details        Returns a reversed pointer to the position after the last 
///                 element.
/// \return         Immutable pointer to the end.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline std::reverse_iterator<const Type*> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rbegin() const
{
    return std::reverse_iterator<const Type*>(&Vectorizer::get(*this, 0)+Size);
}

// Forced immutable reverse iterator to the beginning
/// \brief          Forced immutable reverse iterator to the beginning.
/// \details        Returns a constant reversed pointer to the position after 
///                 the last element.
/// \return         Immutable pointer to the end.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline std::reverse_iterator<const Type*> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::crbegin() const
{
    return std::reverse_iterator<const Type*>(&Vectorizer::get(*this, 0)+Size);
}

// Reverse iterator to the end
/// \brief          Reverse iterator to the end.
/// \details        Returns a reversed pointer to the first element.
/// \return         Pointer to the beginning.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline std::reverse_iterator<Type*> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rend()
{
    return std::reverse_iterator<Type*>(&Vectorizer::get(*this, 0));
}

// Immutable reverse iterator to the end
/// \brief          Immutable reverse iterator to the end.
/// \details        Returns a reversed pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline std::reverse_iterator<const Type*> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rend() const
{
    return std::reverse_iterator<const Type*>(&Vectorizer::get(*this, 0));
}

// Forced immutable reverse iterator to the end
/// \brief          Forced immutable reverse iterator to the end.
/// \details        Returns a constant reversed pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline std::reverse_iterator<const Type*> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::crend() const
{
    return std::reverse_iterator<const Type*>(&Vectorizer::get(*this, 0));
}

// Index of an iterator in the array
/// \brief          Index of an iterator in the array.
/// \details        Returns the index of the element pointed by an iterator or 
///                 a pointer.
/// \tparam         IteratorType (Pointer or iterator type.)
/// \param[in]      it Iterator to the element.
/// \return         Copy of the value of the index of the element.
/// \exception      std::out_of_range Out of range.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename IteratorType> 
inline unsigned int AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::index(const IteratorType& it, typename std::iterator_traits<IteratorType>::iterator_category*) const
{
    if ((&*it < &Vectorizer::get(*this, 0)) || (&*it >= &Vectorizer::get(*this, 0)+Size)) {
        throw std::out_of_range("ERROR = AbstractNArray::index() : out of range");
    }
    return &*it-&Vectorizer::get(*this, 0);
}
// -------------------------------------------------------------------------- //



// ------------------------------- COMPARISON ------------------------------- //
// Check whether all elements are approximately null
/// \brief          Check whether all elements are approximately null.
/// \details        Returns true if all elements are approximately set to their 
///                 default value, returns false otherwise. The condition is 
///                 that their absolute value is less or equal the absolute 
///                 value of the tolerance.
/// \param[in]      tolerance Tolerance of approximation.
/// \return         Copy of the result of the test.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline bool AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::null(const Type& tolerance) const
{
    const Type abstol = std::abs(tolerance);
    bool result = true;
    for (unsigned int i = 0; i < Size; ++i) {
        result &= static_cast<bool>(std::abs(Vectorizer::get(*this, i)) <= abstol);
    }
    return result;
}

// Compare for approximate equality
/// \brief          Compare for approximate equality.
/// \details        Returns true if all elements of the containers are 
///                 approximately equal, returns false otherwise. The condition
///                 is that the difference of all elements of the containers has
///                 to be less or equal the absolute value of the tolerance.
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \param[in]      tolerance Tolerance of approximation.
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType> 
inline bool AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::eq(const GenericType& rhs, const Type& tolerance) const
{
    const Type abstol = std::abs(tolerance);
    bool result = true;
    Vectorizer::check(*this, rhs);
    if (!std::is_unsigned<Type>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result &= static_cast<bool>(std::abs(Vectorizer::get(*this, i)-Vectorizer::get(rhs, i)) <= abstol);
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            result &= static_cast<bool>(Vectorizer::get(*this, i) < Vectorizer::get(rhs, i)) ? static_cast<bool>(Vectorizer::get(rhs, i)-Vectorizer::get(*this, i) <= abstol)
                                                                                             : static_cast<bool>(Vectorizer::get(*this, i)-Vectorizer::get(rhs, i) <= abstol);
        }
    }
    return result;
}

// Compare for noticeable difference
/// \brief          Compare for noticeable difference.
/// \details        Returns true if at least one element is noticeably different
///                 in the two containers, returns false otherwise. The 
///                 condition is that at least one difference of two elements of 
///                 the containers has to be strictly greater than the absolute 
///                 value of the tolerance.
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \param[in]      tolerance Tolerance of approximation.
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType> 
inline bool AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::ne(const GenericType& rhs, const Type& tolerance) const
{
    const Type abstol = std::abs(tolerance);
    bool result = false;
    Vectorizer::check(*this, rhs);
    if (!std::is_unsigned<Type>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result |= static_cast<bool>(std::abs(Vectorizer::get(*this, i)-Vectorizer::get(rhs, i)) > abstol);
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            result |= static_cast<bool>(Vectorizer::get(*this, i) < Vectorizer::get(rhs, i)) ? static_cast<bool>(Vectorizer::get(rhs, i)-Vectorizer::get(*this, i) > abstol)
                                                                                             : static_cast<bool>(Vectorizer::get(*this, i)-Vectorizer::get(rhs, i) > abstol);
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- STATISTICS ------------------------------- //
// Absolute minimum element
/// \brief          Absolute minimum element.
/// \details        Returns a reference to the absolute minimum element of the 
///                 container or masked container.
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \return         Immutable reference to the element.
/// \exception      std::runtime_error Empty search.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class>
inline const Type& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::amin(const Mask& bitmask) const
{
    const Type* result = (std::is_same<std::true_type, Mask>::value) ? (static_cast<const Type*>(&Vectorizer::get(*this, 0))) : (static_cast<const Type*>(nullptr));
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 1; i < Size; ++i) {
            if (static_cast<bool>(std::abs(Vectorizer::get(*this, i)) < std::abs(*result))) {
                result = &Vectorizer::get(*this, i);
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && ((result == nullptr) || static_cast<bool>(std::abs(Vectorizer::get(*this, i)) < std::abs(*result)))) {
                result = &Vectorizer::get(*this, i);
            }
        }
        if (result == nullptr) {
            throw std::runtime_error("ERROR = AbstractNArray::amin() : empty search");
        } 
    }
    return *result;
}

// Maximum element
/// \brief          Maximum element.
/// \details        Returns a reference to the absolute maximum element of the 
///                 container or masked container.
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \return         Immutable reference to the element.
/// \exception      std::runtime_error Empty search.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class>
inline const Type& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::amax(const Mask& bitmask) const
{
    const Type* result = (std::is_same<std::true_type, Mask>::value) ? (static_cast<const Type*>(&Vectorizer::get(*this, 0))) : (static_cast<const Type*>(nullptr));
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 1; i < Size; ++i) {
            if (static_cast<bool>(std::abs(Vectorizer::get(*this, i)) > std::abs(*result))) {
                result = &Vectorizer::get(*this, i);
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && ((result == nullptr) || static_cast<bool>(std::abs(Vectorizer::get(*this, i)) > std::abs(*result)))) {
                result = &Vectorizer::get(*this, i);
            }
        }
        if (result == nullptr) {
            throw std::runtime_error("ERROR = AbstractNArray::amax() : empty search");
        } 
    }
    return *result;
}

// Mean value
/// \brief          Mean value.
/// \details        Computes the mean value of the container or masked 
///                 container. If a coefficient container is provided, the
///                 weighted mean is computed : \f$\bar{x}=\frac{\sum_{i=1}^{n}
///                 w_{i}x_{i}}{\sum_{i=1}^{n}w_{i}}\f$.
/// \tparam         Return Return type.
/// \tparam         Coefficient (Coefficient type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      coefficient Coefficients of weighted mean.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the mean value.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Coefficient, class Mask, class> 
inline Return AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::mean(const Coefficient& coefficient, const Mask& bitmask) const
{
    static const Return nothing = Return();
    Return result = nothing;
    Return total = nothing;
    Vectorizer::check(*this, coefficient);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result += Vectorizer::get(coefficient, i)*Vectorizer::get(*this, i);
            total += Vectorizer::get(coefficient, i); 
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                result += Vectorizer::get(coefficient, i)*Vectorizer::get(*this, i);
                total += Vectorizer::get(coefficient, i); 
            }
        }
    }
    return (total != nothing) ? (result/total) : (nothing);
}

// Standard deviation
/// \brief          Standard deviation.
/// \details        Computes the standard deviation of the container or masked 
///                 container. The correction, oftenly equals to \f$0\f$ or 
///                 \f$-1\f$, allows to correct the standard deviation from the
///                 bias. If a coefficient container is provided, the standard 
///                 deviation of the weighted mean is computed : \f$\sigma=\sqrt
///                 {\frac{\sum_{i=1}^{n}w_{i}\left(x_{i}-\bar{x}\right)^{2}}
///                 {\sum_{i=1}^{n}w_{i}+c}}\f$.
/// \tparam         Return Return type.
/// \tparam         Correction (Correction arithmetic type.)
/// \tparam         Coefficient (Coefficient type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      correction Additive term at denominator for bias control.
/// \param[in]      coefficient Coefficients of weighted mean.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the mean value.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, typename Correction, class Coefficient, class Mask, class> 
inline Return AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::sigma(const Correction& correction, const Coefficient& coefficient, const Mask& bitmask) const
{
    static const Return nothing = Return();
    Return result = nothing;
    Return average = nothing;
    Return total = nothing;
    Vectorizer::check(*this, coefficient);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            average += Vectorizer::get(coefficient, i)*Vectorizer::get(*this, i);
            total += Vectorizer::get(coefficient, i); 
        }
        if (total != nothing) {
            average /= total;
            for (unsigned int i = 0; i < Size; ++i) {
                result += Vectorizer::get(coefficient, i)*(Vectorizer::get(*this, i)-average)*(Vectorizer::get(*this, i)-average);
            }
            result = sqrt(result/(total+correction));
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                average += Vectorizer::get(coefficient, i)*Vectorizer::get(*this, i);
                total += Vectorizer::get(coefficient, i); 
            }
        }
        if (total != nothing) {
            average /= total;
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                    result += Vectorizer::get(coefficient, i)*(Vectorizer::get(*this, i)-average)*(Vectorizer::get(*this, i)-average);
                }
            }
            result = sqrt(result/(total+correction));
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- APPLICATION ------------------------------ //
// Modification by a function pointer
/// \brief          Modification by a function pointer.
/// \details        Modifies the container by applying a function pointer to each 
///                 element. For an element <tt>x</tt>, a function <tt>f</tt> 
///                 and for extra arguments <tt>args...</tt>, an equivalent 
///                 expression is : <tt>x = f(x, args...)</tt>. The return type 
///                 is used as an internal cast before affectation.
/// \tparam         Return Return type.
/// \tparam         First Type of the first argument of the function.
/// \tparam         Types Extra types of the function.
/// \tparam         Args (Extra types.)
/// \param[in]      f Function pointer <tt>Return(First, Types...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Self reference.
/// \warning        In case of implicit cast or extra arguments, types related
///                 to the function have to be specified manually as template
///                 parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, typename First, class... Types, class... Args, class> 
inline Crtp<Type, Parameters...>& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::modify(Return (*f)(First, Types...), Args&&... args)
{
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) = static_cast<Return>(f(Vectorizer::get(*this, i), Vectorizer::get(args, i)...));
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Masked modification by a function pointer
/// \brief          Masked modification by a function pointer.
/// \details        Modifies the container by applying a function pointer to 
///                 each element where the mask is true. For an element 
///                 <tt>x</tt>, a function <tt>f</tt> and for extra arguments 
///                 <tt>args...</tt>, an equivalent expression is : 
///                 <tt>x = f(x, args...)</tt>. The return type is used as an 
///                 internal cast before affectation.
/// \tparam         Return Return type.
/// \tparam         First Type of the first argument of the function.
/// \tparam         Types Extra types of the function.
/// \tparam         Args (Extra types.)
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \param[in]      f Function pointer <tt>Return(First, Types...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Self reference.
/// \warning        In case of implicit cast or extra arguments, types related
///                 to the function have to be specified manually as template
///                 parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, typename First, class... Types, class... Args, class Mask, class> 
inline Crtp<Type, Parameters...>& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::modify(const Mask& bitmask, Return (*f)(First, Types...), Args&&... args)
{
    Vectorizer::check(*this, bitmask);
    for (unsigned int i = 0; i < Size; ++i) {
        if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
            Vectorizer::get(*this, i) = static_cast<Return>(f(Vectorizer::get(*this, i), Vectorizer::get(args, i)...));
        } else {
            Vectorizer::get(*this, i) = static_cast<Return>(Vectorizer::get(*this, i));
        }
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Application of a function pointer
/// \brief          Application of a function pointer.
/// \details        Applies a function pointer to each element of the container 
///                 and returns a copy of the result. For a result <tt>y</tt>, 
///                 an element <tt>x</tt>, a function <tt>f</tt> and for extra 
///                 arguments <tt>args...</tt>, an equivalent expression is : 
///                 <tt>y = f(x, args...)</tt>.
/// \tparam         Return Return type.
/// \tparam         First Type of the first argument of the function.
/// \tparam         Types Extra types of the function.
/// \tparam         Args (Extra types.)
/// \param[in]      f Function pointer <tt>Return(First, Types...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Copy.
/// \warning        In case of implicit cast or extra arguments, types related
///                 to the function have to be specified manually as template
///                 parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, typename First, class... Types, class... Args, class> 
inline Crtp<Return, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::apply(Return (*f)(First, Types...), Args&&... args) const
{
    Crtp<Return, Parameters...> result; 
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = static_cast<Return>(f(Vectorizer::get(*this, i), Vectorizer::get(args, i)...));
    }
    return result;
}

// Masked application of a function pointer
/// \brief          Masked application of a function pointer.
/// \details        Applies a function pointer to each element of the container 
///                 where the mask is true and returns a copy of the result. For
///                 a result <tt>y</tt>, an element <tt>x</tt>, a function 
///                 <tt>f</tt> and for extra arguments <tt>args...</tt>, an 
///                 equivalent expression is : <tt>y = f(x, args...)</tt>.
/// \tparam         Return Return type.
/// \tparam         First Type of the first argument of the function.
/// \tparam         Types Extra types of the function.
/// \tparam         Args (Extra types.)
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \param[in]      f Function pointer <tt>Return(First, Types...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Copy.
/// \warning        In case of implicit cast or extra arguments, types related
///                 to the function have to be specified manually as template
///                 parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, typename First, class... Types, class... Args, class Mask, class> 
inline Crtp<Return, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::apply(const Mask& bitmask, Return (*f)(First, Types...), Args&&... args) const
{
    Crtp<Return, Parameters...> result; 
    Vectorizer::check(*this, bitmask);
    for (unsigned int i = 0; i < Size; ++i) {
        if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
            Vectorizer::get(result, i) = static_cast<Return>(f(Vectorizer::get(*this, i), Vectorizer::get(args, i)...));
        } else {
            Vectorizer::get(result, i) = static_cast<Return>(Vectorizer::get(*this, i));
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- COUNT --------------------------------- //
// Unicity of elements
/// \brief          Unicity of elements.
/// \details        Checks for unicity of elements using the provided function
///                 with the <tt>std::unique()</tt> algorithm. If two elements
///                 are equal, the function returns false, and true otherwise.
/// \tparam         Function (Function type : <tt>bool(Type, Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Function object <tt>bool(Type, Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class> 
inline bool AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::unicity(Function&& f, const Mask& bitmask) const
{
    Type sorted[Size];
    long long int counter = 0;
    bool result = true;
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            sorted[i] = Vectorizer::get(*this, i); 
        }
        std::sort(&sorted[0], &sorted[0]+Size);
        result = (std::is_same<std::equal_to<Type>, Function>::value) ? (std::distance(&sorted[0], std::unique(&sorted[0], &sorted[0]+Size)) == Size) 
                                                                      : (std::distance(&sorted[0], std::unique(&sorted[0], &sorted[0]+Size, f)) == Size);
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                sorted[counter++] = Vectorizer::get(*this, i); 
            }
        }
        if (counter > 0) {
            std::sort(&sorted[0], &sorted[0]+counter);
            result = (std::distance(&sorted[0], std::unique(&sorted[0], &sorted[0]+counter, f)) == counter);
        }
    }
    return result; 
}

// Number of distincts elements
/// \brief          Number of distincts elements.
/// \details        Counts the number of distincts elements using the provided 
///                 function with the <tt>std::unique()</tt> algorithm.
/// \tparam         Function (Function type : <tt>bool(Type, Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Function object <tt>bool(Type, Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class> 
inline unsigned int AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::distinct(Function&& f, const Mask& bitmask) const
{
    Type sorted[Size];
    long long int counter = 0;
    unsigned int result = 0;
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            sorted[i] = Vectorizer::get(*this, i); 
        }
        std::sort(&sorted[0], &sorted[0]+Size);
        result = (std::is_same<std::equal_to<Type>, Function>::value) ? (std::distance(&sorted[0], std::unique(&sorted[0], &sorted[0]+Size))) 
                                                                      : (std::distance(&sorted[0], std::unique(&sorted[0], &sorted[0]+Size, f)));
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                sorted[counter++] = Vectorizer::get(*this, i); 
            }
        }
        if (counter > 0) {
            std::sort(&sorted[0], &sorted[0]+counter);
            result = std::distance(&sorted[0], std::unique(&sorted[0], &sorted[0]+counter, f));
        }
    }
    return result; 
}
// -------------------------------------------------------------------------- //



// ---------------------------------- SORT ---------------------------------- //
// Arrange the container
/// \brief          Arrange the container.
/// \details        Sorts the contents and returns a sorted reference using 
///                 <tt>std::sort</tt>. A comparison function, a list of 
///                 indexes and a boolean mask can be passed. The list of 
///                 indexes is a vectorized container which will be sorted in
///                 the same order as the array. If a scalar is passed, nothing
///                 is done on it. This is equivalent to the <tt>sort()</tt> 
///                 function except that it works on the contents instead of a 
///                 copy.   
/// \tparam         Function (Function type : <tt>bool(Type, Type)</tt>.)
/// \tparam         Indexes (Vectorized container type or dummy scalar.)
/// \tparam         Mask (Mask type.)
/// \tparam         Pair (Inner pair type for indexes.)
/// \param[in]      f Function object <tt>bool(Type, Type)</tt>.
/// \param[in,out]  indexes Vectorized container or dummy scalar to be sorted.
/// \param[in]      bitmask Boolean mask.
/// \return         Self reference.  
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Indexes, class Mask, class Pair, class> 
inline Crtp<Type, Parameters...>& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::arrange(Function&& f, Indexes&& indexes, const Mask& bitmask)
{
    Pair* paired = nullptr;
    unsigned int counter = 0;
    Vectorizer::check(*this, indexes);
    Vectorizer::check(*this, bitmask);
    if ((std::is_same<std::less<Type>, Function>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Indexes>::type>::value) && (!std::is_base_of<Vectorizer, Mask>::value)) {
        std::sort(&Vectorizer::get(*this, 0), &Vectorizer::get(*this, 0)+Size);
    } else if ((!std::is_base_of<Vectorizer, typename std::decay<Indexes>::type>::value) && (!std::is_base_of<Vectorizer, Mask>::value)) {
        std::sort(&Vectorizer::get(*this, 0), &Vectorizer::get(*this, 0)+Size, f); 
    } else {
        paired = new Pair[Size];
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                paired[counter++] = Pair(Vectorizer::get(*this, i), Vectorizer::get(indexes, i));
            }
        }
        std::sort(&paired[0], &paired[0]+counter, [=, &f](const Pair& lhs, const Pair& rhs){return f(lhs.first, rhs.first);});
        counter = 0;
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(indexes, i) = paired[counter].second;
                Vectorizer::get(*this, i) = paired[counter++].first;
            }
        }
        if (paired != nullptr) {
            delete[] paired;
            paired = nullptr;
        }
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Sort
/// \brief          Sort.
/// \details        Sorts the contents and returns a sorted copy using 
///                 <tt>std::sort</tt>. A comparison function, a list of 
///                 indexes and a boolean mask can be passed. The list of 
///                 indexes is a vectorized container which will be sorted in
///                 the same order as the array. If a scalar is passed, nothing
///                 is done on it. This is equivalent to the <tt>arrange()</tt> 
///                 function except that it works on a copy instead of the 
///                 contents. 
/// \tparam         Function (Function type : <tt>bool(Type, Type)</tt>.)
/// \tparam         Indexes (Vectorized container type or dummy scalar.)
/// \tparam         Mask (Mask type.)
/// \tparam         Pair (Inner pair type for indexes.)
/// \param[in]      f Function object <tt>bool(Type, Type)</tt>.
/// \param[in,out]  indexes Vectorized container or dummy scalar to be sorted.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.  
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Indexes, class Mask, class Pair, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::sort(Function&& f, Indexes&& indexes, const Mask& bitmask) const
{
    Crtp<Type, Parameters...> result(*this); 
    Pair* paired = nullptr;
    unsigned int counter = 0;
    Vectorizer::check(*this, indexes);
    Vectorizer::check(*this, bitmask);
    if ((std::is_same<std::less<Type>, Function>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Indexes>::type>::value) && (!std::is_base_of<Vectorizer, Mask>::value)) {
        std::sort(&Vectorizer::get(result, 0), &Vectorizer::get(result, 0)+Size);
    } else if ((!std::is_base_of<Vectorizer, typename std::decay<Indexes>::type>::value) && (!std::is_base_of<Vectorizer, Mask>::value)) {
        std::sort(&Vectorizer::get(result, 0), &Vectorizer::get(result, 0)+Size, f); 
    } else {
        paired = new Pair[Size];
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                paired[counter++] = Pair(Vectorizer::get(result, i), Vectorizer::get(indexes, i));
            }
        }
        std::sort(&paired[0], &paired[0]+counter, [=, &f](const Pair& lhs, const Pair& rhs){return f(lhs.first, rhs.first);});
        counter = 0;
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(indexes, i) = paired[counter].second;
                Vectorizer::get(result, i) = paired[counter++].first;
            }
        }
        if (paired != nullptr) {
            delete[] paired;
            paired = nullptr;
        }
    }
    return result;
}

// Absolute arrange the container
/// \brief          Absolute arrange the container.
/// \details        Sorts the contents by its absolute value and returns a 
///                 sorted reference using <tt>std::sort</tt>. A comparison 
///                 function, a list of indexes and a boolean mask can be 
///                 passed. The list of indexes is a vectorized container which 
///                 will be sorted in the same order as the array. If a scalar 
///                 is passed, nothing is done on it. This is equivalent to the
///                 <tt>asort()</tt> function except that it works on the 
///                 contents instead of a copy.   
/// \tparam         Function (Function type : <tt>bool(Type, Type)</tt>.)
/// \tparam         Indexes (Vectorized container type or dummy scalar.)
/// \tparam         Mask (Mask type.)
/// \tparam         Pair (Inner pair type for indexes.)
/// \param[in]      f Function object <tt>bool(Type, Type)</tt>.
/// \param[in,out]  indexes Vectorized container or dummy scalar to be sorted.
/// \param[in]      bitmask Boolean mask.
/// \return         Self reference.  
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Indexes, class Mask, class Pair, class> 
inline Crtp<Type, Parameters...>& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::aarrange(Function&& f, Indexes&& indexes, const Mask& bitmask)
{
    Pair* paired = nullptr;
    unsigned int counter = 0;
    Vectorizer::check(*this, indexes);
    Vectorizer::check(*this, bitmask);
    if ((!std::is_base_of<Vectorizer, typename std::decay<Indexes>::type>::value) && (!std::is_base_of<Vectorizer, Mask>::value)) {
        std::sort(&Vectorizer::get(*this, 0), &Vectorizer::get(*this, 0)+Size, [=, &f](const Type& lhs, const Type& rhs){return f(std::abs(lhs), std::abs(rhs));}); 
    } else {
        paired = new Pair[Size];
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                paired[counter++] = Pair(Vectorizer::get(*this, i), Vectorizer::get(indexes, i));
            }
        }
        std::sort(&paired[0], &paired[0]+counter, [=, &f](const Pair& lhs, const Pair& rhs){return f(std::abs(lhs.first), std::abs(rhs.first));});
        counter = 0;
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(indexes, i) = paired[counter].second;
                Vectorizer::get(*this, i) = paired[counter++].first;
            }
        }
        if (paired != nullptr) {
            delete[] paired;
            paired = nullptr;
        }
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Absolute sort
/// \brief          Absolute sort.
/// \details        Sorts the contents by its absolute value and returns a 
///                 sorted copy using <tt>std::sort</tt>. A comparison function,
///                 a list of indexes and a boolean mask can be passed. The 
///                 list of indexes is a vectorized container which will be 
///                 sorted in the same order as the array. If a scalar is 
///                 passed, nothing is done on it. This is equivalent to the
///                 <tt>aarrange()</tt> function except that it works on a copy 
///                 instead of the contents.  
/// \tparam         Function (Function type : <tt>bool(Type, Type)</tt>.)
/// \tparam         Indexes (Vectorized container type or dummy scalar.)
/// \tparam         Mask (Mask type.)
/// \tparam         Pair (Inner pair type for indexes.)
/// \param[in]      f Function object <tt>bool(Type, Type)</tt>.
/// \param[in,out]  indexes Vectorized container or dummy scalar to be sorted.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.  
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Indexes, class Mask, class Pair, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::asort(Function&& f, Indexes&& indexes, const Mask& bitmask) const
{
    Crtp<Type, Parameters...> result(*this); 
    Pair* paired = nullptr;
    unsigned int counter = 0;
    Vectorizer::check(*this, indexes);
    Vectorizer::check(*this, bitmask);
    if ((!std::is_base_of<Vectorizer, typename std::decay<Indexes>::type>::value) && (!std::is_base_of<Vectorizer, Mask>::value)) {
        std::sort(&Vectorizer::get(result, 0), &Vectorizer::get(result, 0)+Size, [=, &f](const Type& lhs, const Type& rhs){return f(std::abs(lhs), std::abs(rhs));}); 
    } else {
        paired = new Pair[Size];
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                paired[counter++] = Pair(Vectorizer::get(result, i), Vectorizer::get(indexes, i));
            }
        }
        std::sort(&paired[0], &paired[0]+counter, [=, &f](const Pair& lhs, const Pair& rhs){return f(std::abs(lhs.first), std::abs(rhs.first));});
        counter = 0;
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(indexes, i) = paired[counter].second;
                Vectorizer::get(result, i) = paired[counter++].first;
            }
        }
        if (paired != nullptr) {
            delete[] paired;
            paired = nullptr;
        }
    }
    return result;
}

// Re-arrange as a function of indexes
/// \brief          Re-arrange as a function of indexes.
/// \details        Sorts the contents by the provided indexes previously 
///                 generated by an <tt>arrange()</tt> or a <tt>sort()</tt> 
///                 function. It can be used to sort several arrays in the
///                 same order of a single one. If no indexes or a scalar is
///                 provided, the container is randomly shuffled. This is 
///                 equivalent to the <tt>resort()</tt> function except that it 
///                 works on the contents instead of a copy.   
/// \tparam         Indexes (Vectorized container type or dummy scalar.)
/// \tparam         Pair (Inner pair type for indexes.)
/// \param[in]      indexes Vectorized container of indexes.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Indexes, class Pair> 
inline Crtp<Type, Parameters...>& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rearrange(const Indexes& indexes)
{
    Type* original = nullptr;
    Pair* paired = nullptr;
    Vectorizer::check(*this, indexes);
    if (!std::is_base_of<Vectorizer, Indexes>::value) {
        std::random_shuffle(&Vectorizer::get(*this, 0), &Vectorizer::get(*this, 0)+Size);
    } else {
        original = new Type[Size];
        paired = new Pair[Size];
        for (unsigned int i = 0; i < Size; ++i) {
            original[i] = Vectorizer::get(*this, i);
            paired[i] = Pair(Vectorizer::get(indexes, i), i);
        }
        std::sort(&paired[0], &paired[0]+Size, [](const Pair& lhs, const Pair& rhs){return lhs.first < rhs.first;});
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(*this, paired[i].second) = original[i];
        }
        if (original != nullptr) {
            delete[] original;
            original = nullptr;
        }
        if (paired != nullptr) {
            delete[] paired;
            paired = nullptr;
        }
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Re-sort as a function of indexes
/// \brief          Re-sort as a function of indexes.
/// \details        Sorts a copy by the provided indexes previously generated 
///                 by an <tt>arrange()</tt> or a <tt>sort()</tt> function. It 
///                 can be used to sort several arrays in the same order of a 
///                 single one. If no indexes or a scalar is provided, the 
///                 container is randomly shuffled. This is equivalent to the 
///                 <tt>rearrange()</tt> function except that it works on a 
///                 copy instead of the contents.   
/// \tparam         Indexes (Vectorized container type or dummy scalar.)
/// \tparam         Pair (Inner pair type for indexes.)
/// \param[in]      indexes Vectorized container of indexes.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Indexes, class Pair> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::resort(const Indexes& indexes) const
{
    Crtp<Type, Parameters...> result(*this); 
    Pair* paired = nullptr;
    Vectorizer::check(*this, indexes);
    if (!std::is_base_of<Vectorizer, Indexes>::value) {
        std::random_shuffle(&Vectorizer::get(result, 0), &Vectorizer::get(result, 0)+Size);
    } else {
        paired = new Pair[Size];
        for (unsigned int i = 0; i < Size; ++i) {
            paired[i] = Pair(Vectorizer::get(indexes, i), i);
        }
        std::sort(&paired[0], &paired[0]+Size, [](const Pair& lhs, const Pair& rhs){return lhs.first < rhs.first;});
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, paired[i].second) = Vectorizer::get(*this, i);
        }
        if (paired != nullptr) {
            delete[] paired;
            paired = nullptr;
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------- MATHEMATICAL FUNCTIONS ------------------------- //
// Absolute value
/// \brief          Absolute value.
/// \details        Applies the <tt>std::abs()</tt> function to each element. 
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class>
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::abs(const Mask& bitmask) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::abs(Vectorizer::get(*this, i));
        } 
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(result, i) = std::abs(Vectorizer::get(*this, i));
            } else {
                Vectorizer::get(result, i) = Vectorizer::get(*this, i);
            }
        }
    }
    return result;
}

// Power metafunction
/// \brief          Power metafunction.
/// \details        Recursively calls the multiplication operator at 
///                 compile-time computing the integer exponentiation of the 
///                 array. Negative exponent inverses the result.
/// \tparam         Exponent Exponent of the power function.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <int Exponent> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::pow() const
{
    return (Exponent > 1) ? (static_cast<const Crtp<Type, Parameters...>&>(*this)*pow<(Exponent-1)*(Exponent > 1)>())
                          : ((Exponent < 0) ? (Crtp<Type, Parameters...>(1)/pow<(-Exponent)*(Exponent < 0)>())
                                            : ((Exponent == 1) ? (static_cast<const Crtp<Type, Parameters...>&>(*this))
                                                               : (Crtp<Type, Parameters...>(1))));
}

// Power
/// \brief          Power.
/// \details        Computes the power for the specified exponent of each 
///                 element using the <tt>std::pow()</tt> function.
/// \tparam         GenericType (Value or vectorized type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType, class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::pow(const GenericType& rhs, const Mask& bitmask) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::pow(Vectorizer::get(*this, i), Vectorizer::get(rhs, i));
        } 
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(result, i) = std::pow(Vectorizer::get(*this, i), Vectorizer::get(rhs, i));
            } else {
                Vectorizer::get(result, i) = Vectorizer::get(*this, i);
            }
        }
    }
    return result;
}

// Root metafunction
/// \brief          Root metafunction.
/// \details        Calls the inverse of the power function for each element
///                 at compile-time. For 2 and 3 exponents, the 
///                 <tt>std::sqrt()</tt> and <tt>std::cbrt()</tt> functions are
///                 called, otherwise <tt>std::pow()</tt> with the inverse 
///                 exponent is called. At the end, the result is casted to the
///                 array type.
/// \tparam         Degree Root degree.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <int Degree> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rt() const
{
    static_assert(Degree != 0, "ERROR = AbstractNArray::root() : null root");
    static const typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type num = 1;
    Crtp<Type, Parameters...> result; 
    if (Degree == 1) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = Vectorizer::get(*this, i);
        } 
    } else if (Degree == 2) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::sqrt(Vectorizer::get(*this, i));
        } 
    } else if (Degree == 3) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::cbrt(Vectorizer::get(*this, i));
        } 
    } else if (Degree == -2) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = 1/std::sqrt(Vectorizer::get(*this, i));
        } 
    } else if (Degree == -3) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = 1/std::cbrt(Vectorizer::get(*this, i));
        } 
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::pow(Vectorizer::get(*this, i), num/Degree);
        }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    }
    return result;
}

// Root
/// \brief          Root.
/// \details        Computes the root for the specified degree of each element 
///                 using the <tt>std::pow()</tt> function.
/// \tparam         GenericType (Value or vectorized type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType, class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::rt(const GenericType& rhs, const Mask& bitmask) const
{
    Crtp<Type, Parameters...> result; 
    static const typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type num = 1;
    Vectorizer::check(*this, rhs);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::pow(Vectorizer::get(*this, i), num/Vectorizer::get(rhs, i));
        } 
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(result, i) = std::pow(Vectorizer::get(*this, i), num/Vectorizer::get(rhs, i));
            } else {
                Vectorizer::get(result, i) = Vectorizer::get(*this, i);
            }
        }
    }
    return result;
}

// Log metafunction
/// \brief          Log metafunction.
/// \details        Computes the logarithm with compile-time options. If the 
///                 base is 0, then the natural base logarithm is called. If the
///                 data type is a floating point, and if the base is 2 or 10, 
///                 then the <tt>std::log2()</tt> or <tt>std::log10</tt> 
///                 function is called. For other bases, the division implying 
///                 the logarithm of the base is used. If the data type is not a 
///                 floating point, the exact integral part of the logarithm in
///                 the specified base is computed using loops.
/// \tparam         Base Logarithm base.
/// \return         Copy.
/// \exception      std::domain_error Logarithm of a negative integer undefined.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <unsigned int Base> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::log() const
{
    static_assert(Base != 1, "ERROR = AbstractNArray::log() : unknown base");
    static const Type logbase = std::log(static_cast<Type>(Base+2*(Base <= 2)))+(Base <= 2);
    Crtp<Type, Parameters...> result; 
    int base = 0;
    int counter = 0;
    if (Base == 0) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::log(Vectorizer::get(*this, i));
        } 
    } else if (std::is_floating_point<Type>::value) {
        if (Base == 2) {
            for (unsigned int i = 0; i < Size; ++i) {
                Vectorizer::get(result, i) = std::log2(Vectorizer::get(*this, i));
            } 
        } else if (Base == 10) {
            for (unsigned int i = 0; i < Size; ++i) {
                Vectorizer::get(result, i) = std::log10(Vectorizer::get(*this, i));
            } 
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                Vectorizer::get(result, i) = std::log(Vectorizer::get(*this, i))/logbase;
            } 
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (Vectorizer::get(*this, i) <= 0) {
                throw std::domain_error("ERROR = AbstractNArray::log() : logarithm of a negative integer undefined");
            }
            base = Base;
            counter = 0;
            while (static_cast<Type>(base) <= Vectorizer::get(*this, i)) {
                base *= Base;
                ++counter;
            }
            Vectorizer::get(result, i) = counter;
        }
    }
    return result;
}

// Logarithm
/// \brief          Logarithm.
/// \details        Computes the logarithm in the specified base of each element
///                 using the <tt>std::log()</tt> function.
/// \tparam         GenericType (Value or vectorized type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType, class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::log(const GenericType& rhs, const Mask& bitmask) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) = std::log(Vectorizer::get(*this, i))/std::log(Vectorizer::get(rhs, i));
        } 
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(result, i) = std::log(Vectorizer::get(*this, i))/std::log(Vectorizer::get(rhs, i));
            } else {
                Vectorizer::get(result, i) = Vectorizer::get(*this, i);
            }
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- NORM ---------------------------------- //
// Norm
/// \brief          Norm.
/// \details        Computes the p-norm of the array : \f$\|x\|_{p} = \left(
///                 \sum_{i=1}^n\left|x_{i}\right|^p\right)^{1/p}\f$ where 
///                 \f$p\f$ is specified by the <tt>NormType</tt> template 
///                 parameter. If it is equal to zero, then the infinite norm is
///                 taken.
/// \tparam         Degree Norm degree.
/// \tparam         NormType Norm type.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the norm.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <unsigned int Degree, typename NormType, class Mask, class> 
inline NormType AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::norm(const Mask& bitmask) const
{
    static const typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type inv = static_cast<typename std::conditional<std::is_floating_point<NormType>::value, NormType, double>::type>(1)/(Degree+(Degree <= 0));
    NormType result = NormType();
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        if (Degree > 0) {
            for (unsigned int i = 0; i < Size; ++i) {
                result += std::pow(std::abs(Vectorizer::get(*this, i)), Degree);
            }
            result = std::pow(result, inv);
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(std::abs(Vectorizer::get(*this, i)) > result)) {
                    result = std::abs(Vectorizer::get(*this, i));
                }
            }
        }
    } else {
        if (Degree > 0) {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                    result += std::pow(std::abs(Vectorizer::get(*this, i)), Degree);
                }
            }
            result = std::pow(result, inv);
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(std::abs(Vectorizer::get(*this, i)) > result)) {
                    result = std::abs(Vectorizer::get(*this, i));
                }
            }
        }
    }
    return result; 
}

// Renormalize
/// \brief          Renormalize.
/// \details        Renormalizes the contents using the the p-norm of the array 
///                 : \f$\|x\|_{p} = \left(\sum_{i=1}^n\left|x_{i}\right|^p
///                 \right)^{1/p}\f$ where \f$p\f$ is specified by the 
///                 <tt>NormType</tt> template parameter. If it is equal to 
///                 zero, then the infinite norm is taken. This is equivalent to
///                 the <tt>normalize()</tt> function except that it works on 
///                 the contents instead of a copy.
/// \tparam         Degree Norm degree.
/// \tparam         NormType Norm type.
/// \param[in]      bitmask Boolean mask.
/// \return         Self reference.
/// \exception      std::domain_error The norm is not normal.
/// \exception      std::domain_error The norm is not normal.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <unsigned int Degree, typename NormType, class Mask, class> 
inline Crtp<Type, Parameters...>& AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::renormalize(const Mask& bitmask)
{
    static const typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type inv = static_cast<typename std::conditional<std::is_floating_point<NormType>::value, NormType, double>::type>(1)/(Degree+(Degree <= 0));
    NormType n = NormType();
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        if (Degree > 0) {
            for (unsigned int i = 0; i < Size; ++i) {
                n += std::pow(std::abs(Vectorizer::get(*this, i)), Degree);
            }
            n = std::pow(n, inv);
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(std::abs(Vectorizer::get(*this, i)) > n)) {
                    n = std::abs(Vectorizer::get(*this, i));
                }
            }
        }
        if (!std::isnormal(n)) {
            throw std::domain_error("ERROR = AbstractNArray::renormalize() : the norm is not normal");
        }
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(*this, i) /= n;
        }
    } else {
        if (Degree > 0) {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                    n += std::pow(std::abs(Vectorizer::get(*this, i)), Degree);
                }
            }
            n = std::pow(n, inv);
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(std::abs(Vectorizer::get(*this, i)) > n)) {
                    n = std::abs(Vectorizer::get(*this, i));
                }
            }
        }
        if (!std::isnormal(n)) {
            throw std::domain_error("ERROR = AbstractNArray::renormalize() : the norm is not normal");
        }
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(*this, i) /= n;
            }
        }
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Normalize
/// \brief          Normalize.
/// \details        Normalizes the contents using the the p-norm of the array 
///                 : \f$\|x\|_{p} = \left(\sum_{i=1}^n\left|x_{i}\right|^p
///                 \right)^{1/p}\f$ where \f$p\f$ is specified by the 
///                 <tt>NormType</tt> template parameter. If it is equal to 
///                 zero, then the infinite norm is taken. This is equivalent to
///                 the <tt>renormalize()</tt> function except that it works on 
///                 a copy instead of the contents.
/// \tparam         Degree Norm degree.
/// \tparam         NormType Norm type.
/// \param[in]      bitmask Boolean mask.
/// \return         Self reference.
/// \exception      std::domain_error The norm is not normal.
/// \exception      std::domain_error The norm is not normal.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <unsigned int Degree, typename NormType, class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::normalize(const Mask& bitmask) const      
{
    static const typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type inv = static_cast<typename std::conditional<std::is_floating_point<NormType>::value, NormType, double>::type>(1)/(Degree+(Degree <= 0));
    Crtp<Type, Parameters...> result(*this); 
    NormType n = NormType();
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        if (Degree > 0) {
            for (unsigned int i = 0; i < Size; ++i) {
                n += std::pow(std::abs(Vectorizer::get(result, i)), Degree);
            }
            n = std::pow(n, inv);
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(std::abs(Vectorizer::get(result, i)) > n)) {
                    n = std::abs(Vectorizer::get(result, i));
                }
            }
        }
        if (!std::isnormal(n)) {
            throw std::domain_error("ERROR = AbstractNArray::normalize() : the norm is not normal");
        }
        for (unsigned int i = 0; i < Size; ++i) {
            Vectorizer::get(result, i) /= n;
        }
    } else {
        if (Degree > 0) {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                    n += std::pow(std::abs(Vectorizer::get(result, i)), Degree);
                }
            }
            n = std::pow(n, inv);
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(std::abs(Vectorizer::get(result, i)) > n)) {
                    n = std::abs(Vectorizer::get(result, i));
                }
            }
        }
        if (!std::isnormal(n)) {
            throw std::domain_error("ERROR = AbstractNArray::normalize() : the norm is not normal");
        }
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                Vectorizer::get(result, i) /= n;
            }
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Zero array
/// \brief          Zero array.
/// \details        Creates an array filled with zeroes.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::zero()
{
    return Crtp<Type, Parameters...>(static_cast<Type>(0));
}

// One array
/// \brief          One array.
/// \details        Creates an array filled with ones.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::one()
{
    return Crtp<Type, Parameters...>(static_cast<Type>(1));
}

// Array from value
/// \brief          Array from value.
/// \details        Creates an array filled with values.
/// \tparam         Mask (Mask type.)
/// \param[in]      source Source value.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::value(const Type& source, const Mask& bitmask)
{
    Crtp<Type, Parameters...> result;
    Vectorizer::check(result, bitmask);
    return (!std::is_base_of<Vectorizer, Mask>::value) ? (Vectorizer::set(result, source)) : (Vectorizer::set(result, source, bitmask)); 
}

// Basic random array creation
/// \brief          Basic random array creation.
/// \details        Creates a random array filled with random values generated
///                 with a uniform distribution over the provided interval.
/// \tparam         Mask (Mask type.)
/// \param[in]      minimum Minimum value of the distribution.
/// \param[in]      maximum Maximum value of the distribution.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
/// \warning        As the internal engine is a static one, do not use this 
///                 function in parallel.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::random(const Type& minimum, const Type& maximum, const Mask& bitmask)
{
    static std::random_device device;
    static std::mt19937 engine(device());
    std::uniform_real_distribution<typename std::conditional<std::is_floating_point<Type>::value, Type, double>::type> rdistribution(minimum, maximum);
    std::uniform_int_distribution<typename std::conditional<!std::is_floating_point<Type>::value, Type, int>::type> idistribution(minimum, maximum);
    Crtp<Type, Parameters...> result;
    Vectorizer::check(result, bitmask);
    if (std::is_floating_point<Type>::value) {
        if (!std::is_base_of<Vectorizer, Mask>::value) {
            for (unsigned int i = 0; i < Size; ++i) {
                result[i] = rdistribution(engine);
            }
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                    result[i] = rdistribution(engine);
                }
            }
        }
    } else {
        if (!std::is_base_of<Vectorizer, Mask>::value) {
            for (unsigned int i = 0; i < Size; ++i) {
                result[i] = idistribution(engine);
            }
        } else {
            for (unsigned int i = 0; i < Size; ++i) {
                if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                    result[i] = idistribution(engine);
                }
            }
        }
    }
    return result;
}

// Generic random array creation
/// \brief          Generic random array creation.
/// \details        Creates a random array filled with random values generated
///                 from the provided distribution and engine.
/// \tparam         Engine (Random engine type.)
/// \tparam         Distribution (Random distribution type.)
/// \tparam         Mask (Mask type.)
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Engine, class Distribution, class Mask, class> 
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::random(Engine&& engine, Distribution&& distribution, const Mask& bitmask)
{
    Crtp<Type, Parameters...> result;
    Vectorizer::check(result, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result[i] = distribution(engine);
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                result[i] = distribution(engine);
            }
        }
    }
    return result;
}

// Indexed array
/// \brief          Indexed array.
/// \details        Creates an array based on indexes. The passed function 
///                 should take an index value and return the associated
///                 element value.
/// \tparam         Function (Function type : <tt>Type(unsigned int)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Function object <tt>Type(unsigned int)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class>
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::indexed(Function&& f, const Mask& bitmask)
{
    Crtp<Type, Parameters...> result;
    Vectorizer::check(result, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result[i] = f(i);
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                result[i] = f(i);
            }
        }
    }
    return result;
}

// Progressive array
/// \brief          Progressive array.
/// \details        Creates an array filled with progressive values. The next 
///                 non-masked element is equal to the result of the function
///                 applied to the previous non-masked element and the step.
/// \tparam         Function (Function type : <tt>Type(Type, Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      init Initial value.
/// \param[in]      step Step of progression.
/// \param[in]      f Function object <tt>Type(Type, Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class>
inline Crtp<Type, Parameters...> AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::progressive(const Type& init, const Type& step, Function&& f, const Mask& bitmask)
{
    Crtp<Type, Parameters...> result;
    Type* previous = nullptr;
    Vectorizer::check(result, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        result[0] = init;
        for (unsigned int i = 1; i < Size; ++i) {
            result[i] = f(result[i-1], step);
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
                result[i] = (previous != nullptr) ? (static_cast<Type>(f(*previous, step))) : (static_cast<Type>(init));
                previous = &result[i];
            }
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of AbstractNArray.
/// \return         0 if no error.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
int AbstractNArray<Kind, Size, Crtp, Type, Parameters...>::example()
{
    std::cout<<"BEGIN = AbstractNArray::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = AbstractNArray::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = AbstractNArray::example()"<<std::endl;
    return 1;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABSTRACTNARRAY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

