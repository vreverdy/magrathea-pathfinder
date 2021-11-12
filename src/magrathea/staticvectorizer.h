/* **************************** STATICVECTORIZER **************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          StaticVectorizer
// DESCRIPTION :    Helper class for generic constant size vectorization
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           staticvectorizer.h
/// \brief          Helper class for generic constant size vectorization
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef STATICVECTORIZER_H_INCLUDED
#define STATICVECTORIZER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <iterator>
#include <initializer_list>
#include <vector>
#include <bitset>
// Include libs
// Include project
#include "vectorizer.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Helper class for generic constant size vectorization
/// \brief          Helper class for generic constant size vectorization.
/// \details        Provides vectorization for constant size containers thanks 
///                 to the curiously recurring template pattern (CRTP) trick. 
///                 To use it, one has to derive from this class and pass the 
///                 derived class itself as the CRTP parameter. The derived 
///                 classes have to satisfy the conditions required by the 
///                 Vectorizer base class  and have to implement the following 
///                 functions required by CRTP : <ul>
///                 <li><tt>operator[](const unsigned int)</tt></li></ul>
///                 One can also modify members like <tt>operator()</tt> to 
///                 change the behaviour of the function. 
/// \tparam         Kind Kind of arguments.
/// \tparam         Size Number of elements.
/// \tparam         Crtp Derived CRTP class.
/// \tparam         Type Data type.
/// \tparam         Parameters List of parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
class StaticVectorizer
: public Vectorizer
{
    // Setup
    static_assert(Size > 0, "ERROR = StaticVectorizer : null size");
    
    // Protected lifecycle
    /// \name           Protected lifecycle
    //@{
    protected: 
        inline ~StaticVectorizer();
    //@}
    
    // Vectorization
    /// \name           Vectorization
    //@{
    public:
        inline Type& operator[](const unsigned int i);
        inline const Type& operator[](const unsigned int i) const;
        inline Crtp<Type, Parameters...>& resize(const unsigned int n);
    //@}

    // Operators : assignment
    /// \name           Operators : assignment
    //@{
    public:
        inline Crtp<Type, Parameters...>& operator=(const StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& rhs);
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<Type, Parameters...>& operator=(const std::initializer_list<OtherType>& rhs);
        template <class Misc> inline Crtp<Type, Parameters...>& operator=(const Misc& rhs);
    //@}
    
    // Operators : compound assignement
    /// \name           Operators : compound assignement
    //@{
    public:
        template <class GenericType> inline Crtp<Type, Parameters...>& operator+=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator-=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator*=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator/=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator%=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator&=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator|=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator^=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator<<=(const GenericType& rhs);
        template <class GenericType> inline Crtp<Type, Parameters...>& operator>>=(const GenericType& rhs);
    //@}

    // Operators : main
    /// \name           Operators : main
    //@{
    public:
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator+(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator-(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator*(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator/(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator%(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator&(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator|(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator^(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<Type, Parameters...> operator<<(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<Type, Parameters...> operator>>(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator&&(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator||(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator==(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator!=(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator>(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator<(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator>=(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
        template <typename OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<bool, Parameters...> operator<=(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const;
    //@}

    // Operators : with rhs value
    /// \name           Operators : with rhs value
    //@{
    public:
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator+(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator-(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator*(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator/(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator%(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator&(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator|(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> operator^(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<Type, Parameters...> operator<<(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<Type, Parameters...> operator>>(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator&&(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator||(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator==(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator!=(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator>(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator<(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator>=(const OtherType& rhs) const;
        template <typename OtherType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, Type>::value)>::type> inline Crtp<bool, Parameters...> operator<=(const OtherType& rhs) const;
    //@}

    // Operators : with lhs value
    /// \name           Operators : with lhs value
    //@{
    public:
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator+(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator-(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator*(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator/(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator%(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator&(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator|(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator^(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<OtherType, SelfParameters...> operator<<(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<OtherType, SelfParameters...> operator>>(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator&&(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);                
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator||(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs); 
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator==(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator!=(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator>(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator<(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator>=(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class> friend inline SelfCrtp<bool, SelfParameters...> operator<=(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
    //@}
    
    // Operators : unary
    /// \name           Operators : unary
    //@{
    public:
        inline Crtp<bool, Parameters...> operator!() const;
        inline Crtp<Type, Parameters...> operator~() const;
        inline Crtp<Type, Parameters...> operator+() const;
        inline Crtp<Type, Parameters...> operator-() const;
        inline Crtp<Type, Parameters...>& operator++();
        inline Crtp<Type, Parameters...>& operator--();
        inline Crtp<Type, Parameters...> operator++(int);
        inline Crtp<Type, Parameters...> operator--(int);
    //@}

    // Access
    /// \name           Access
    //@{
    public:
        inline StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& operator()();
        inline const StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& operator()() const;
        inline Type& operator()(const unsigned int i);
        inline const Type& operator()(const unsigned int i) const;  
        inline Type& at(const unsigned int i);
        inline const Type& at(const unsigned int i) const;  
        inline Type& front(const unsigned int i = 0);
        inline const Type& front(const unsigned int i = 0) const;  
        inline Type& back(const unsigned int i = 0);
        inline const Type& back(const unsigned int i = 0) const;  
        inline Type& cycle(const int i);
        inline const Type& cycle(const int i) const; 
    //@}

    // Assignment
    /// \name           Assignment
    //@{
    public:
        template <typename OtherType, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<Type, Parameters...>& assign(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc> inline Crtp<Type, Parameters...>& assign(const Misc&... misc);
        template <typename OtherType, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<Type, Parameters...>& fill(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc> inline Crtp<Type, Parameters...>& fill(const Misc&... misc);
        template <typename OtherType, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<Type, Parameters...> replace(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc> inline Crtp<Type, Parameters...> replace(const Misc&... misc);
        template <class GenericType> inline Crtp<Type, Parameters...>& put(const GenericType& source, const unsigned int pos, const unsigned int num = 1);
        template <class GenericType> inline Crtp<Type, Parameters...> change(const GenericType& source, const unsigned int pos, const unsigned int num = 1);
    //@}
    
    // Management
    /// \name           Management
    //@{
    public:
        inline Crtp<Type, Parameters...>& reserve(const unsigned int n);
        inline Crtp<Type, Parameters...>& clear();
        template <class... Location> inline Crtp<Type, Parameters...>& nullify(const Location&... location);
        template <class... Location> inline Crtp<Type, Parameters...>& swap(StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& rhs, const Location&... location);
        inline Crtp<Type, Parameters...> copy() const;
        template <typename OtherType = Type, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Crtp<OtherType, Parameters...> cast() const;
    //@}
    
    // Comparison
    /// \name           Comparison
    //@{
    public:
        inline bool null() const;
        template <class GenericType> inline bool eq(const GenericType& rhs) const; 
        template <class GenericType> inline bool ne(const GenericType& rhs) const; 
    //@}

    // Statistics
    /// \name           Statistics
    //@{
    public:
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline const Type& min(const Mask& bitmask = Mask()) const;
        template <class Mask = std::true_type, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value)>::type> inline const Type& max(const Mask& bitmask = Mask()) const;
    //@}
    
    // Application
    /// \name           Application
    //@{
    public:
        template <typename Return = Type, class Function, class... Args, class = typename std::enable_if<!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value>::type> inline Crtp<Type, Parameters...>& modify(Function&& f, Args&&... args);
        template <typename Return = Type, class Mask, class Function, class... Args, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)>::type> inline Crtp<Type, Parameters...>& modify(const Mask& bitmask, Function&& f, Args&&... args);
        template <typename Return = Type, class Function, class... Args, class = typename std::enable_if<!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value>::type> inline Crtp<Return, Parameters...> apply(Function&& f, Args&&... args) const;
        template <typename Return = Type, class Mask, class Function, class... Args, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)>::type> inline Crtp<Return, Parameters...> apply(const Mask& bitmask, Function&& f, Args&&... args) const;
        template <typename Return = Type, class Function = std::plus<Type>, class = typename std::enable_if<!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value>::type> inline Return reduce(Function&& f = Function(), const Return& init = Return()) const;
        template <typename Return = Type, class Mask, class Function = std::plus<Type>, class = typename std::enable_if<(std::is_base_of<Vectorizer, Mask>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)>::type> inline Return reduce(const Mask& bitmask, Function&& f = Function(), const Return& init = Return()) const;
        template <typename Result = void, class Function, class Arg, class... Args, typename Return = typename std::conditional<std::is_void<Result>::value, typename std::common_type<Type, Arg, Args...>::type, Result>::type> inline Crtp<Return, Parameters...> combine(Function&& f, const StaticVectorizer<Kind, Size, Crtp, Arg, Parameters...>& arg, const StaticVectorizer<Kind, Size, Crtp, Args, Parameters...>&... args) const;
        template <typename Return = Type, class Function> inline Crtp<Return, Parameters...> combine(Function&&) const;
    //@}
    
    // Count
    /// \name           Count
    //@{
    public:
        template <class Reference = bool, class Mask = std::true_type, class = typename std::enable_if<((std::is_convertible<Type, Reference>::value) || (std::is_base_of<Vectorizer, Reference>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline unsigned int count(const Reference& r = true, const Mask& bitmask = Mask()) const;
        template <class Function, class Mask = std::true_type, class = typename std::enable_if<((!std::is_convertible<Type, typename std::decay<Function>::type>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline unsigned int count(Function&& f, const Mask& bitmask = Mask()) const;
        template <class Reference = bool, class Mask = std::true_type, class = typename std::enable_if<((std::is_convertible<Type, Reference>::value) || (std::is_base_of<Vectorizer, Reference>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool all(const Reference& r = true, const Mask& bitmask = Mask()) const;
        template <class Function, class Mask = std::true_type, class = typename std::enable_if<((!std::is_convertible<Type, typename std::decay<Function>::type>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool all(Function&& f, const Mask& bitmask = Mask()) const;
        template <class Reference = bool, class Mask = std::true_type, class = typename std::enable_if<((std::is_convertible<Type, Reference>::value) || (std::is_base_of<Vectorizer, Reference>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool any(const Reference& r = true, const Mask& bitmask = Mask()) const;
        template <class Function, class Mask = std::true_type, class = typename std::enable_if<((!std::is_convertible<Type, typename std::decay<Function>::type>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool any(Function&& f, const Mask& bitmask = Mask()) const;
        template <class Reference = bool, class Mask = std::true_type, class = typename std::enable_if<((std::is_convertible<Type, Reference>::value) || (std::is_base_of<Vectorizer, Reference>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool none(const Reference& r = true, const Mask& bitmask = Mask()) const;
        template <class Function, class Mask = std::true_type, class = typename std::enable_if<((!std::is_convertible<Type, typename std::decay<Function>::type>::value) && (!std::is_base_of<Vectorizer, typename std::decay<Function>::type>::value)) && ((std::is_base_of<Vectorizer, Mask>::value) || (std::is_same<std::true_type, Mask>::value))>::type> inline bool none(Function&& f, const Mask& bitmask = Mask()) const;
    //@}

    // Stream
    /// \name           Stream
    //@{
    public:
        template <typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters> friend std::ostream& operator<<(std::ostream& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
        template <typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters> friend std::istream& operator>>(std::istream& lhs, StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs);
    //@}

    // Static vectorization
    /// \name           Static vectorization
    //@{
    public:
        static constexpr unsigned int size();
        static constexpr bool constant();
        static constexpr bool boolean();
        static constexpr std::array<Kind, sizeof...(Parameters)> parameters();
        static inline Type type();
    //@}
    
    // Size
    /// \name           Size
    //@{
    public:
        static constexpr bool empty();
        static constexpr unsigned int capacity();
        static constexpr unsigned int tbytes();
        static constexpr unsigned long long int bytes();
        static inline unsigned long long int space();
    //@}
    
    // Predefined
    /// \name           Predefined
    //@{
    public:
        static inline Crtp<bool, Parameters...> mask(const bool value = true);
        template <class Container, class = typename std::enable_if<(!std::is_base_of<Vectorizer, Container>::value) && (std::is_convertible<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type, bool>::value) && ((std::is_void<decltype(std::declval<Container>().flip())>::value) || (std::is_reference<decltype(std::declval<Container>().flip())>::value))>::type> static inline Crtp<bool, Parameters...> mask(const Container& container);
        template <typename OtherType, class... Misc, class = typename std::enable_if<std::is_convertible<OtherType, bool>::value>::type> static inline Crtp<bool, Parameters...> mask(const std::initializer_list<OtherType>& source, const Misc&... misc);
        template <class... Misc, class = typename std::enable_if<sizeof...(Misc) != 0>::type> static inline Crtp<bool, Parameters...> mask(const Misc&... misc);
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
/// \details        Does nothing.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::~StaticVectorizer()
= default;
// -------------------------------------------------------------------------- //



// ----------------------------- VECTORIZATION ------------------------------ //
// Direct access to the element
/// \brief          Direct access to the element.
/// \details        Provides a direct access to the specified element. This 
///                 function is required by the vectorization mechanism. 
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator[](const unsigned int i)
{
    return static_cast<Crtp<Type, Parameters...>&>(*this)[i];
}

// Immutable direct access to the element
/// \brief          Immutable direct access to the element.
/// \details        Provides a constant direct access to the specified element.
///                 This function is required by the vectorization mechanism.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator[](const unsigned int i) const
{
    return static_cast<const Crtp<Type, Parameters...>&>(*this)[i];
}

// Resize the container
/// \brief          Resize the container.
/// \details        Resizes the container and returns a reference to it. This 
///                 function is required by the vectorization mechanism. 
/// \param[in]      n New size.
/// \return         Self reference.
/// \exception      std::length_error The container cannot be resized.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::resize(const unsigned int n)
{
    if (n != Size) {
        throw std::length_error("ERROR = StaticVectorizer::resize() : the container cannot be resized");
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}
// -------------------------------------------------------------------------- //



// ------------------------- OPERATORS : ASSIGNMENT ------------------------- //
// Copy assignment operator
/// \brief          Copy assignment operator.
/// \details        Copies the contents of the container. 
/// \param[in]      rhs Right-hand side.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator=(const StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& rhs)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, rhs));
}

// Initializer list assignment operator
/// \brief          Initializer list assignment operator.
/// \details        Provides an initializer list assignment. The assignment is 
///                 delegated to the <tt>set()</tt> helper function thanks to 
///                 the following call : <tt>set(*this, rhs)</tt>.
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator=(const std::initializer_list<OtherType>& rhs)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, rhs));
}

// Generic assignment operator
/// \brief          Copy assignment operator.
/// \details        Provides a generic assignment. The assignment is delegated 
///                 to the <tt>set()</tt> helper function. Conversion assignment, 
///                 value assignment and assignment are provided through this
///                 function thanks to the following call : <tt>set(*this,
///                  misc...)</tt>.
/// \tparam         Misc (Miscellaneous type.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Misc> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator=(const Misc& rhs)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, rhs));
}
// -------------------------------------------------------------------------- //



// --------------------- OPERATORS : COMPOUND ASSIGNMENT -------------------- //
// Addition assignment 
/// \brief          Addition assignment.
/// \details        Applies the addition assignment operator to each element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator+=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) += Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Substraction assignment 
/// \brief          Substraction assignment. 
/// \details        Applies the substraction assignment operator to each 
///                 element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator-=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) -= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Multiplication assignment 
/// \brief          Multiplication assignment. 
/// \details        Applies the multiplication assignment operator to each 
///                 element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator*=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) *= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Division assignment 
/// \brief          Division assignment. 
/// \details        Applies the division assignment operator to each element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator/=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) /= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Modulo assignment 
/// \brief          Modulo assignment. 
/// \details        Applies the modulo assignment operator to each element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator%=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) %= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Bitwise AND assignment 
/// \brief          Bitwise AND assignment. 
/// \details        Applies the bitwise AND assignment operator to each element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator&=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) &= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Bitwise OR assignment 
/// \brief          Bitwise OR assignment.
/// \details        Applies the bitwise OR assignment operator to each element. 
/// \tparam         GenericType (Value or vectorized type.) 
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator|=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) |= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Bitwise XOR assignment 
/// \brief          Bitwise XOR assignment. 
/// \details        Applies the bitwise XOR assignment operator to each element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator^=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) ^= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Bitwise left shift assignment 
/// \brief          Bitwise left shift assignment. 
/// \details        Applies the bitwise left shift assignment operator to each 
///                 element. 
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<<=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) <<= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Bitwise right shift assignment 
/// \brief          Bitwise right shift assignment. 
/// \details        Applies the bitwise right shift assignment operator to each 
///                 element.
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Self reference
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>>=(const GenericType& rhs)
{
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) >>= Vectorizer::get(rhs, i);
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}
// -------------------------------------------------------------------------- //



// ---------------------------- OPERATORS : MAIN ---------------------------- //
// Addition
/// \brief          Addition.
/// \details        Applies the addition operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator+(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) + Vectorizer::get(rhs, i));
    } 
    return result;
}

// Substraction
/// \brief          Substraction.
/// \details        Applies the substraction operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Common type copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator-(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) - Vectorizer::get(rhs, i));
    } 
    return result;
}

// Multiplication
/// \brief          Multiplication.
/// \details        Applies the multiplication operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator*(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) * Vectorizer::get(rhs, i));
    } 
    return result;
}

// Division
/// \brief          Division.
/// \details        Applies the division operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator/(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) / Vectorizer::get(rhs, i));
    } 
    return result;
}

// Modulo
/// \brief          Modulo.
/// \details        Applies the modulo operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator%(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) % Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise AND
/// \brief          Bitwise AND.
/// \details        Applies the bitwise AND operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator&(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) & Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise OR
/// \brief          Bitwise OR.
/// \details        Applies the bitwise OR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator|(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) | Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise XOR
/// \brief          Bitwise XOR.
/// \details        Applies the bitwise XOR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator^(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) ^ Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise left shift
/// \brief          Bitwise left shift.
/// \details        Applies the bitwise left shift operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<<(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) << Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise right shift
/// \brief          Bitwise right shift.
/// \details        Applies the bitwise right shift operator to each element.
/// \tparam         OtherType (Other data type.) 
/// \param[in]      rhs Right-hand side. 
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>>(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) >> Vectorizer::get(rhs, i));
    } 
    return result;
}

// Logical AND
/// \brief          Logical AND.
/// \details        Applies the logical AND operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator&&(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) && Vectorizer::get(rhs, i));
    } 
    return result;
}

// Logical OR
/// \brief          Logical OR.
/// \details        Applies the logical OR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator||(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) || Vectorizer::get(rhs, i));
    } 
    return result;
}

// Equal to
/// \brief          Equal to.
/// \details        Applies the equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator==(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) == Vectorizer::get(rhs, i));
    } 
    return result;
}

// Not equal to
/// \brief          Not equal to.
/// \details        Applies the not equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator!=(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) != Vectorizer::get(rhs, i));
    } 
    return result;
}

// Greater than
/// \brief          Greater than.
/// \details        Applies the greater than operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) > Vectorizer::get(rhs, i));
    } 
    return result;
}

// Less than
/// \brief          Less than.
/// \details        Applies the less than operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) < Vectorizer::get(rhs, i));
    } 
    return result;
}

// Greater than or equal to
/// \brief          Greater than or equal to.
/// \details        Applies the greater than or equal to operator to each 
///                 element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>=(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) >= Vectorizer::get(rhs, i));
    } 
    return result;
}

// Less than or equal to
/// \brief          Less than or equal to.
/// \details        Applies the less than or equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<=(const StaticVectorizer<Kind, Size, Crtp, OtherType, Parameters...>& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) <= Vectorizer::get(rhs, i));
    } 
    return result;
}
// -------------------------------------------------------------------------- //



// ----------------------- OPERATORS : WITH RHS VALUE ----------------------- //
// Addition with rhs value
/// \brief          Addition with rhs value.
/// \details        Applies the addition operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator+(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) + Vectorizer::get(rhs, i));
    } 
    return result;
}

// Substraction with rhs value
/// \brief          Substraction with rhs value.
/// \details        Applies the substraction operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Common type copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator-(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) - Vectorizer::get(rhs, i));
    } 
    return result;
}

// Multiplication with rhs value
/// \brief          Multiplication with rhs value.
/// \details        Applies the multiplication operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator*(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) * Vectorizer::get(rhs, i));
    } 
    return result;
}

// Division with rhs value
/// \brief          Division with rhs value.
/// \details        Applies the division operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator/(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) / Vectorizer::get(rhs, i));
    } 
    return result;
}

// Modulo with rhs value
/// \brief          Modulo with rhs value.
/// \details        Applies the modulo operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator%(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) % Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise AND with rhs value
/// \brief          Bitwise AND with rhs value.
/// \details        Applies the bitwise AND operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator&(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) & Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise OR with rhs value
/// \brief          Bitwise OR with rhs value.
/// \details        Applies the bitwise OR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator|(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) | Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise XOR with rhs value
/// \brief          Bitwise XOR with rhs value.
/// \details        Applies the bitwise XOR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Common type copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator^(const OtherType& rhs) const
{
    Crtp<typename std::common_type<Type, OtherType>::type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) ^ Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise left shift with rhs value
/// \brief          Bitwise left shift with rhs value.
/// \details        Applies the bitwise left shift operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<<(const OtherType& rhs) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) << Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise right shift with rhs value
/// \brief          Bitwise right shift with rhs value.
/// \details        Applies the bitwise right shift operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>>(const OtherType& rhs) const
{
    Crtp<Type, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) >> Vectorizer::get(rhs, i));
    } 
    return result;
}

// Logical AND with rhs value
/// \brief          Logical AND with rhs value.
/// \details        Applies the logical AND operator to each element.
/// \tparam         OtherType (Other data type.) 
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator&&(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) && Vectorizer::get(rhs, i));
    } 
    return result;
}

// Logical OR with rhs value
/// \brief          Logical OR with rhs value.
/// \details        Applies the logical OR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator||(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) || Vectorizer::get(rhs, i));
    } 
    return result;
}

// Equal to with rhs value
/// \brief          Equal to with rhs value.
/// \details        Applies the equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator==(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) == Vectorizer::get(rhs, i));
    } 
    return result;
}

// Not equal to with rhs value
/// \brief          Not equal to with rhs value.
/// \details        Applies the not equal to operator to each element.
/// \tparam         OtherType (Other data type.) 
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator!=(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) != Vectorizer::get(rhs, i));
    } 
    return result;
}

// Greater than with rhs value
/// \brief          Greater than with rhs value.
/// \details        Applies the greater than operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) > Vectorizer::get(rhs, i));
    } 
    return result;
}

// Less than with rhs value
/// \brief          Less than with rhs value.
/// \details        Applies the less than operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) < Vectorizer::get(rhs, i));
    } 
    return result;
}

// Greater than or equal to with rhs value
/// \brief          Greater than or equal to with rhs value.
/// \details        Applies the greater than or equal to operator to each 
///                 element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator>=(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) >= Vectorizer::get(rhs, i));
    } 
    return result;
}

// Less than or equal to with rhs value
/// \brief          Less than or equal to with rhs value.
/// \details        Applies the less than or equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator<=(const OtherType& rhs) const
{
    Crtp<bool, Parameters...> result; 
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(*this, i) <= Vectorizer::get(rhs, i));
    } 
    return result;
}
// -------------------------------------------------------------------------- //



// ----------------------- OPERATORS : WITH LHS VALUE ----------------------- //
// Addition with lhs value
/// \brief          Addition with lhs value.
/// \details        Applies the addition operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type> 
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator+(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) + Vectorizer::get(rhs, i));
    } 
    return result;
}

// Substraction with lhs value
/// \brief          Substraction with lhs value.
/// \details        Applies the substraction operator to each element.
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side. 
/// \param[in]      rhs Right-hand side. 
/// \return         Common type copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type> 
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator-(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) - Vectorizer::get(rhs, i));
    } 
    return result;
}

// Multiplication with lhs value
/// \brief          Multiplication with lhs value.
/// \details        Applies the multiplication operator to each element.
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator*(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) * Vectorizer::get(rhs, i));
    } 
    return result;
}

// Division with lhs value
/// \brief          Division with lhs value.
/// \details        Applies the division operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator/(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) / Vectorizer::get(rhs, i));
    } 
    return result;
}

// Modulo with lhs value
/// \brief          Modulo with lhs value.
/// \details        Applies the modulo operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator%(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) % Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise AND with lhs value
/// \brief          Bitwise AND with lhs value.
/// \details        Applies the bitwise AND operator to each element.
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side. 
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator&(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) & Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise OR with lhs value
/// \brief          Bitwise OR with lhs value.
/// \details        Applies the bitwise OR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator|(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) | Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise XOR with lhs value
/// \brief          Bitwise XOR with lhs value.
/// \details        Applies the bitwise XOR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side. 
/// \return         Common type copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> operator^(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs) 
{
    SelfCrtp<typename std::common_type<SelfType, OtherType>::type, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) ^ Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise left shift with lhs value
/// \brief          Bitwise left shift with lhs value.
/// \details        Applies the bitwise left shift operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side. 
/// \return         Copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<OtherType, SelfParameters...> operator<<(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<OtherType, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) << Vectorizer::get(rhs, i));
    } 
    return result;
}

// Bitwise right shift with lhs value
/// \brief          Bitwise right shift with lhs value.
/// \details        Applies the bitwise right shift operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side. 
/// \return         Copy.
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type> 
inline SelfCrtp<OtherType, SelfParameters...> operator>>(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<OtherType, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) >> Vectorizer::get(rhs, i));
    } 
    return result;
}

// Logical AND with lhs value
/// \brief          Logical AND with lhs value.
/// \details        Applies the logical AND operator to each element.
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side. 
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator&&(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) && Vectorizer::get(rhs, i));
    } 
    return result;
}

// Logical OR with lhs value
/// \brief          Logical OR with lhs value.
/// \details        Applies the logical OR operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator||(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) || Vectorizer::get(rhs, i));
    } 
    return result;
}

// Equal to with lhs value
/// \brief          Equal to with lhs value.
/// \details        Applies the equal to operator to each element.
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side. 
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator==(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) == Vectorizer::get(rhs, i));
    } 
    return result;
}

// Not equal to with lhs value
/// \brief          Not equal to with lhs value.
/// \details        Applies the not equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type> 
inline SelfCrtp<bool, SelfParameters...> operator!=(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) != Vectorizer::get(rhs, i));
    } 
    return result;
}

// Greater than with lhs value
/// \brief          Greater than with lhs value.
/// \details        Applies the greater than operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator>(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) > Vectorizer::get(rhs, i));
    } 
    return result;
}

// Less than with lhs value
/// \brief          Less than with lhs value.
/// \details        Applies the less than operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator<(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) < Vectorizer::get(rhs, i));
    } 
    return result;
}

// Greater than or equal to with lhs value
/// \brief          Greater than or equal to with lhs value.
/// \details        Applies the greater than or equal to operator to each 
///                 element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side. 
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator>=(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) >= Vectorizer::get(rhs, i));
    } 
    return result;
}

// Less than or equal to with lhs value
/// \brief          Less than or equal to with lhs value.
/// \details        Applies the less than or equal to operator to each element. 
/// \tparam         OtherType (Other data type.)
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in]      lhs Left-hand side.
/// \param[in]      rhs Right-hand side.
/// \return         Boolean copy. 
template <typename OtherType, typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters, class = typename std::enable_if<(!std::is_base_of<Vectorizer, OtherType>::value) && (std::is_convertible<OtherType, SelfType>::value)>::type>
inline SelfCrtp<bool, SelfParameters...> operator<=(const OtherType& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    SelfCrtp<bool, SelfParameters...> result; 
    Vectorizer::check(lhs, rhs);
    for (unsigned int i = 0; i < SelfSize; ++i) {
        Vectorizer::get(result, i) = (Vectorizer::get(lhs, i) <= Vectorizer::get(rhs, i));
    } 
    return result;
}
// -------------------------------------------------------------------------- //



// --------------------------- OPERATORS : UNARY ---------------------------- //
// Logical NOT 
/// \brief          Logical NOT.
/// \details        Applies the logical NOT operator to each element.
/// \return         Boolean copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator!() const
{
    Crtp<bool, Parameters...> result; 
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = !Vectorizer::get(*this, i);
    } 
    return result;
}

// Bitwise NOT 
/// \brief          Bitwise NOT.
/// \details        Applies the bitwise NOT operator to each element.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator~() const
{
    Crtp<Type, Parameters...> result; 
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = ~Vectorizer::get(*this, i);
    } 
    return result;
}

// Integer promotion
/// \brief          Integer promotion.
/// \details        Applies the integer promotion operator to each element.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator+() const
{
    Crtp<Type, Parameters...> result; 
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = +Vectorizer::get(*this, i);
    } 
    return result;
}

// Additive inverse
/// \brief          Additive inverse.
/// \details        Applies the additive inverse operator to each element.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator-() const
{
    Crtp<Type, Parameters...> result; 
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = -Vectorizer::get(*this, i);
    } 
    return result;
}

// Increment prefix
/// \brief          Increment prefix.
/// \details        Applies the increment prefix operator to each element.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator++()
{
    for (unsigned int i = 0; i < Size; ++i) {
        ++Vectorizer::get(*this, i);
    } 
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Decrement prefix
/// \brief          Decrement prefix.
/// \details        Applies the decrement prefix operator to each element.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator--()
{
    for (unsigned int i = 0; i < Size; ++i) {
        --Vectorizer::get(*this, i);
    } 
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Increment suffix
/// \brief          Increment suffix.
/// \details        Applies the increment suffix operator to each element.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator++(int)
{
    Crtp<Type, Parameters...> result(*this);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i)++;
    }
    return result;
}

// Decrement suffix
/// \brief          Decrement suffix.
/// \details        Applies the decrement suffix operator to each element.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator--(int)
{
    Crtp<Type, Parameters...> result(*this);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i)--;
    }
    return result;
}
// -------------------------------------------------------------------------- //



// --------------------------------- ACCESS --------------------------------- //
// Abstract class access
/// \brief          Abstract class access.
/// \details        Casts to the abstract class.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator()()
{
    return *this;
}

// Immutable abstract class access
/// \brief          Immutable abstract class access.
/// \details        Casts to the abstract class.
/// \return         Immutable self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator()() const
{
    return *this;
}

// Monodimensional access operator
/// \brief          Monodimensional access operator.
/// \details        Provides a monodimensional access to the element. For a 
///                 monodimensional array, it is equivalent to the <tt>[]</tt>
///                 operator.
/// \param[in]      i Element index.
/// \return         Reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator()(const unsigned int i)
{
    return Vectorizer::get(*this, i);
}

// Immutable monodimensional access operator
/// \brief          Immutable monodimensional access operator.
/// \details        Provides a monodimensional access to the element. For a 
///                 monodimensional array, it is equivalent to the <tt>[]</tt>
///                 operator.
/// \param[in]      i Element index.
/// \return         Immutable reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::operator()(const unsigned int i) const
{
    return Vectorizer::get(*this, i);
}

// Monodimensional access with range-check
/// \brief          Monodimensional access with range-check.
/// \details        Provides a monodimensional access to the element with a 
///                 range-check. Due to the check this function is slower than
///                 the <tt>[]</tt> or the <tt>()</tt> operator.
/// \param[in]      i Element index.
/// \return         Reference to the element.
/// \exception      std::out_of_range Out of range.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::at(const unsigned int i)
{
    if (i >= Size) {
        throw std::out_of_range("ERROR = StaticVectorizer::at() : out of range");
    }
    return Vectorizer::get(*this, i);
}

// Immutable monodimensional access with range-check
/// \brief          Immutable monodimensional access with range-check.
/// \details        Provides a monodimensional access to the element with a 
///                 range-check. Due to the check this function is slower than
///                 the <tt>[]</tt> or the <tt>()</tt> operator.
/// \param[in]      i Element index.
/// \return         Immutable reference to the element.
/// \exception      std::out_of_range Out of range.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::at(const unsigned int i) const
{
    if (i >= Size) {
        throw std::out_of_range("ERROR = StaticVectorizer::at() : out of range");
    }
    return Vectorizer::get(*this, i);
}

// Mnodimensional access to the i-th element from the beginning
/// \brief          Monodimensional access to the i-th element from the 
///                 beginning.
/// \details        Returns a reference to the i-th first element in the 
///                 container without doing any range check.
/// \param[in]      i Element index.
/// \return         Reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::front(const unsigned int i)
{
    return Vectorizer::get(*this, i);
}

// Immutable monodimensional access to the i-th element from the beginning
/// \brief          Immutable monodimensional access to the i-th element from 
///                 the beginning.
/// \details        Returns a reference to the i-th first element in the 
///                 container without doing any range check.
/// \param[in]      i Element index.
/// \return         Immutable reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::front(const unsigned int i) const
{
    return Vectorizer::get(*this, i);
}

// Monodimensional access to the i-th element from the end
/// \brief          Monodimensional access to the i-th element from the end.
/// \details        Returns a reference to the i-th last element in the 
///                 container without doing any range check.
/// \param[in]      i Element index.
/// \return         Reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::back(const unsigned int i)
{
    return Vectorizer::get(*this, static_cast<int>(Size)-1-static_cast<int>(i));
}

// Immutable monodimensional access to the i-th element from the end
/// \brief          Immutable monodimensional access to the i-th element from 
///                 the end.
/// \details        Returns a reference to the i-th last element in the 
///                 container without doing any range check.
/// \param[in]      i Element index.
/// \return         Immutable reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::back(const unsigned int i) const
{
    return Vectorizer::get(*this, static_cast<int>(Size)-1-static_cast<int>(i));
}

// Cyclic monodimensional access to the contents
/// \brief          Cyclic monodimensional access to the contents.
/// \details        Provides a cyclic access to the contents, using the index
///                 modulo. Negative indexes are supported. It allows to iterate
///                 several times over the contents just by incrementing the
///                 provided index.
/// \param[in]      i Element index.
/// \return         Reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::cycle(const int i)
{
    static const int isize = Size;
    return Vectorizer::get(*this, (i >= 0) ? (i%isize) : ((isize-((-i)%isize))%isize));
}

// Immutable cyclic monodimensional access to the contents
/// \brief          Immutable cyclic monodimensional access to the contents.
/// \details        Provides a cyclic access to the contents, using the index
///                 modulo. Negative indexes are supported. It allows to iterate
///                 several times over the contents just by incrementing the
///                 provided index.
/// \param[in]      i Element index.
/// \return         Immutable reference to the element.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::cycle(const int i) const
{
    static const int isize = Size;
    return Vectorizer::get(*this, (i >= 0) ? (i%isize) : ((isize-((-i)%isize))%isize));
}
// -------------------------------------------------------------------------- //



// ------------------------------- ASSIGNMENT ------------------------------- //
// Initializer list assignment
/// \brief          Initializer list assignment.
/// \details        Provides an assignment from an initializer list equivalent 
///                 to a call to a constructor. Before any operation, the 
///                 contents is reinitialized to its default value.
/// \tparam         OtherType (Other data type.)
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      source Source of the copy. 
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class... Misc, class> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::assign(const std::initializer_list<OtherType>& source, const Misc&... misc)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(Vectorizer::set(*this, Type()), source, misc...));
}

// Generic assignment
/// \brief          Generic assignment.
/// \details        Provides a generic assignment equivalent to a call to a 
///                 constructor. Before any operation, the contents is 
///                 reinitialized to its default value.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Misc> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::assign(const Misc&... misc)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(Vectorizer::set(*this, Type()), misc...));
}

// Initializer list fill
/// \brief          Initializer list fill.
/// \details        Fills the contents using the <tt>set()</tt> function and
///                 without reinitializing the contents before any operation.
///                 This is equivalent to the <tt>replace()</tt> function but
///                 it operates on the container itself and not on a copy. 
/// \tparam         OtherType (Other data type.)
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      source Source of the copy. 
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class... Misc, class> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::fill(const std::initializer_list<OtherType>& source, const Misc&... misc)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, source, misc...));
}

// Generic fill
/// \brief          Generic fill.
/// \details        Fills the contents using the <tt>set()</tt> function and
///                 without reinitializing the contents before any operation.
///                 This is equivalent to the <tt>replace()</tt> function but
///                 it operates on the container itself and not on a copy. 
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Misc> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::fill(const Misc&... misc)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, misc...));
}

// Initializer list replace
/// \brief          Initializer list replace.
/// \details        Replaces the contents using the <tt>set()</tt> function and
///                 without reinitializing the contents before any operation.
///                 This is equivalent to the <tt>fill()</tt> function but
///                 it operates on a copy and not on the container itself. 
/// \tparam         OtherType (Other data type.)
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      source Source of the copy. 
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Self reference.

template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class... Misc, class> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::replace(const std::initializer_list<OtherType>& source, const Misc&... misc)
{
    Crtp<Type, Parameters...> result(*this);
    return Vectorizer::set(result, source, misc...);
}

// Generic replace
/// \brief          Generic replace.
/// \details        Replaces the contents using the <tt>set()</tt> function and
///                 without reinitializing the contents before any operation.
///                 This is equivalent to the <tt>fill()</tt> function but
///                 it operates on a copy and not on the container itself. 
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Misc> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::replace(const Misc&... misc)
{
    Crtp<Type, Parameters...> result(*this);
    return Vectorizer::set(result, misc...);
}

// Put an element in the container
/// \brief          Put an element in the container.
/// \details        Provides a single element assignment to modify the contents.
///                 This function is well suited for chaining assignment. This 
///                 is equivalent to the <tt>change()</tt> function but it 
///                 operates on the container itself and not on a copy.  
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      source Source of the copy.
/// \param[in]      pos Starting position of the copy.
/// \param[in]      num Number of elements to copy.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::put(const GenericType& source, const unsigned int pos, const unsigned int num)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, source, pos, num));
}

// Change an element of the container
/// \brief          Change an element of the container.
/// \details        Provides a single element assignment to modify the contents.
///                 This function is well suited for chaining assignment. This 
///                 is equivalent to the <tt>put()</tt> function but it operates
///                 on a copy and not on the container itself.  
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      source Source of the copy.
/// \param[in]      pos Starting position of the copy.
/// \param[in]      num Number of elements to copy.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType> 
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::change(const GenericType& source, const unsigned int pos, const unsigned int num)
{
    Crtp<Type, Parameters...> result(*this);
    return Vectorizer::set(result, source, pos, num);
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Reserve new space for the container
/// \brief          Reserve new space for the container.
/// \details        Reserves new space for the container in order to optimize
///                 future resize calls.
/// \param[in]      n New size for reservation.
/// \return         Self reference.
/// \exception      std::length_error The container cannot be resized.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::reserve(const unsigned int n)
{
    if (n != Size) {
        throw std::length_error("ERROR = StaticVectorizer::reserve() : the container cannot be resized");
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Clear contents
/// \brief          Clear contents.
/// \details        Clear the whole contents and sets the size to zero.
/// \return         Self reference.
/// \exception      std::length_error The container cannot be resized.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::clear()
{
    if (Size != 0) {
        throw std::length_error("ERROR = StaticVectorizer::clear() : the container cannot be resized");
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Set elements to their default values.
/// \brief          Set elements to their default values.
/// \details        Sets the specified elements to their default value using
///                 the <tt>set()</tt> function for the provided location.
/// \tparam         Location (Mask or position specifiers.)
/// \param[in]      location Boolean mask or starting position and number of
///                 elements.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Location>
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::nullify(const Location&... location)
{
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, Type(), location...));
}

// Swap elements by copy
/// \brief          Swap elements by copy.
/// \details        Swaps the elements of the two containers by copy at the 
///                 provided location and returns a reference to the left-hand
///                 side. The resulting operation is not optimal because it 
///                 implies a temporary copy.
/// \tparam         Location (Mask or position specifiers.)
/// \param[in]      rhs Right-hand side.
/// \param[in]      location Boolean mask or starting position and number of
///                 elements.
/// \return         Self reference.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Location> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::swap(StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>& rhs, const Location&... location)
{
    const Crtp<Type, Parameters...> tmp(rhs);
    Vectorizer::set(rhs, *this, location...);
    return static_cast<Crtp<Type, Parameters...>&>(Vectorizer::set(*this, tmp, location...));
}

// Copy
/// \brief          Copy.
/// \details        Returns a copy of the container.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<Type, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::copy() const
{
    return Crtp<Type, Parameters...>(*this);
}

// Cast to a different data type
/// \brief          Cast to a different data type.
/// \details        Returns a copy of the container converted to another data
///                 type.
/// \tparam         OtherType Other data type.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class> 
inline Crtp<OtherType, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::cast() const
{
    return Crtp<OtherType, Parameters...>(*this);
}
// -------------------------------------------------------------------------- //



// ------------------------------- COMPARISON ------------------------------- //
// Check whether all elements are null
/// \brief          Check whether all elements are null.
/// \details        Returns true if all elements are set to their default value,
///                 returns false otherwise.
/// \return         Copy of the result of the test.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::null() const
{
    static const Type nothing = Type();
    bool result = true;
    for (unsigned int i = 0; i < Size; ++i) {
        result &= static_cast<bool>(Vectorizer::get(*this, i) == nothing);
    }
    return result;
}

// Compare for equality
/// \brief          Compare for equality.
/// \details        Returns true if all elements of the containers are equal, 
///                 returns false otherwise.
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::eq(const GenericType& rhs) const
{
    bool result = true;
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        result &= static_cast<bool>(Vectorizer::get(*this, i) == Vectorizer::get(rhs, i));
    }
    return result;
}

// Compare for difference
/// \brief          Compare for difference.
/// \details        Returns true if at least one element is different in the two
///                 containers, returns false otherwise.
/// \tparam         GenericType (Value or vectorized type.)
/// \param[in]      rhs Right-hand side. 
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class GenericType> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::ne(const GenericType& rhs) const
{
    bool result = false;
    Vectorizer::check(*this, rhs);
    for (unsigned int i = 0; i < Size; ++i) {
        result |= static_cast<bool>(Vectorizer::get(*this, i) != Vectorizer::get(rhs, i));
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- STATISTICS ------------------------------- //
// Minimum element
/// \brief          Minimum element.
/// \details        Returns a reference to the minimum element of the container 
///                 or masked container.
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \return         Immutable reference to the element.
/// \exception      std::runtime_error Empty search.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::min(const Mask& bitmask) const
{
    const Type* result = (std::is_same<std::true_type, Mask>::value) ? (static_cast<const Type*>(&Vectorizer::get(*this, 0))) : (static_cast<const Type*>(nullptr));
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 1; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(*this, i) < *result)) {
                result = &Vectorizer::get(*this, i);
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && ((result == nullptr) || static_cast<bool>(Vectorizer::get(*this, i) < *result))) {
                result = &Vectorizer::get(*this, i);
            }
        }
        if (result == nullptr) {
            throw std::runtime_error("ERROR = StaticVectorizer::min() : empty search");
        } 
    }
    return *result;
}

// Maximum element
/// \brief          Maximum element.
/// \details        Returns a reference to the maximum element of the container 
///                 or masked container.
/// \tparam         Mask (Mask type.)
/// \param[in]      bitmask Boolean mask.
/// \return         Immutable reference to the element.
/// \exception      std::runtime_error Empty search.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Mask, class>
inline const Type& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::max(const Mask& bitmask) const
{
    const Type* result = (std::is_same<std::true_type, Mask>::value) ? (static_cast<const Type*>(&Vectorizer::get(*this, 0))) : (static_cast<const Type*>(nullptr));
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 1; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(*this, i) > *result)) {
                result = &Vectorizer::get(*this, i);
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && ((result == nullptr) || static_cast<bool>(Vectorizer::get(*this, i) > *result))) {
                result = &Vectorizer::get(*this, i);
            }
        }
        if (result == nullptr) {
            throw std::runtime_error("ERROR = StaticVectorizer::max() : empty search");
        } 
    }
    return *result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- APPLICATION ------------------------------ //
// Modification by a function object
/// \brief          Modification by a function object.
/// \details        Modifies the container by applying a function object to each 
///                 element. For an element <tt>x</tt>, a function <tt>f</tt> 
///                 and for extra arguments <tt>args...</tt>, an equivalent 
///                 expression is : <tt>x = f(x, args...)</tt>. The return type
///                 is used as an internal cast before affectation.
/// \tparam         Return Return type.
/// \tparam         Function (Function type : <tt>Type(Type, Args...)</tt>.)
/// \tparam         Args (Extra types.)
/// \param[in]      f Function object <tt>Type(Type, Args...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Self reference.
/// \warning        The size of the extra arguments is not checked.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Function, class... Args, class> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::modify(Function&& f, Args&&... args)
{
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(*this, i) = static_cast<Return>(f(Vectorizer::get(*this, i), Vectorizer::get(args, i)...));
    }
    return static_cast<Crtp<Type, Parameters...>&>(*this);
}

// Masked modification by a function object
/// \brief          Masked modification by a function object.
/// \details        Modifies the container by applying a function object to each 
///                 element where the mask is true. For an element <tt>x</tt>,
///                 a function <tt>f</tt> and for extra arguments <tt>args...
///                 </tt>, an equivalent expression is : <tt>x = f(x, args...)
///                 </tt>. The return type is used as an internal cast before 
///                 affectation.
/// \tparam         Return Return type.
/// \tparam         Mask (Mask type.)
/// \tparam         Function (Function type : <tt>Type(Type, Args...)</tt>.)
/// \tparam         Args (Extra types.)
/// \param[in]      bitmask Boolean mask.
/// \param[in]      f Function object <tt>Type(Type, Args...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Self reference.
/// \warning        The size of the extra arguments is not checked.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Mask, class Function, class... Args, class> 
inline Crtp<Type, Parameters...>& StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::modify(const Mask& bitmask, Function&& f, Args&&... args)
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

// Application of a function object
/// \brief          Application of a function object.
/// \details        Applies a function object to each element of the container 
///                 and returns a copy of the result. For a result <tt>y</tt>, 
///                 an element <tt>x</tt>, a function <tt>f</tt> and for extra 
///                 arguments <tt>args...</tt>, an equivalent expression is : 
///                 <tt>y = f(x, args...)</tt>.
/// \tparam         Return Return type.
/// \tparam         Function (Function type : <tt>Return(Type, Args...)</tt>.)
/// \tparam         Args (Extra types.)
/// \param[in]      f Function object <tt>Return(Type, Args...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Copy.
/// \warning        The size of the extra arguments is not checked.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Function, class... Args, class> 
inline Crtp<Return, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::apply(Function&& f, Args&&... args) const
{
    Crtp<Return, Parameters...> result; 
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = static_cast<Return>(f(Vectorizer::get(*this, i), Vectorizer::get(args, i)...));
    }
    return result;
}

// Masked application of a function object
/// \brief          Masked application of a function object.
/// \details        Applies a function object to each element of the container 
///                 where the mask is true and returns a copy of the result. For
///                 a result <tt>y</tt>, an element <tt>x</tt>, a function 
///                 <tt>f</tt> and for extra arguments <tt>args...</tt>, an 
///                 equivalent expression is : <tt>y = f(x, args...)</tt>.
/// \tparam         Return Return type.
/// \tparam         Mask (Mask type.)
/// \tparam         Function (Function type : <tt>Return(Type, Args...)</tt>.)
/// \tparam         Args (Extra types.)
/// \param[in]      bitmask Boolean mask.
/// \param[in]      f Function object <tt>Return(Type, Args...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Copy.
/// \warning        The size of the extra arguments is not checked.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Mask, class Function, class... Args, class> 
inline Crtp<Return, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::apply(const Mask& bitmask, Function&& f, Args&&... args) const
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

// Reduction by a function object
/// \brief          Reduction by a function object.
/// \details        Reduces the contents using a binary function object 
///                 initialized to the <tt>init</tt> value. For each reduced
///                 element <tt>x</tt>, the equivalent expression is <tt>result
///                 = f(result, x)</tt>.
/// \tparam         Return Return type.
/// \tparam         Function (Function type : <tt>Return(Return, Type)</tt>.)
/// \param[in]      f Function object <tt>Return(Return, Type)</tt>.
/// \param[in]      init Initial value for the reduction.
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Function, class> 
inline Return StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::reduce(Function&& f, const Return& init) const
{
    Return result = init;
    for (unsigned int i = 0; i < Size; ++i) {
        result = f(result, Vectorizer::get(*this, i));
    }
    return result;
}

// Masked reduction by a function object
/// \brief          Reduction by a function object.
/// \details        Reduces the contents using a binary function object 
///                 initialized to the <tt>init</tt> value where the mask is 
///                 true. For each reduced element <tt>x</tt>, the equivalent 
///                 expression is <tt>result = f(result, x)</tt>.
/// \tparam         Return Return type.
/// \tparam         Mask (Mask type.)
/// \tparam         Function (Function type : <tt>Return(Return, Type)</tt>.)
/// \param[in]      bitmask Boolean mask.
/// \param[in]      f Function object <tt>Return(Return, Type)</tt>.
/// \param[in]      init Initial value for the reduction.
/// \return         Copy of the result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Mask, class Function, class> 
inline Return StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::reduce(const Mask& bitmask, Function&& f, const Return& init) const
{
    Return result = init;
    for (unsigned int i = 0; i < Size; ++i) {
        if (static_cast<bool>(Vectorizer::get(bitmask, i))) {
            result = f(result, Vectorizer::get(*this, i));
        }
    }
    return result;
}

// Combination by a function object
/// \brief          Combination by a function object.
/// \details        Combines several vectorized containers using a function 
///                 object. It is equivalent to a transversal reduction 
///                 operation.
/// \tparam         Result Specified return type.
/// \tparam         Function (Function type : <tt>Return(Return, Type)</tt>.)
/// \tparam         Arg (First argument type.)
/// \tparam         Args (Other argument types.)
/// \tparam         Return Return type.
/// \param[in]      f Function object <tt>Return(Return, Type)</tt>.
/// \param[in]      arg First argument.
/// \param[in]      args Other arguments.
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Result, class Function, class Arg, class... Args, typename Return> 
inline Crtp<Return, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::combine(Function&& f, const StaticVectorizer<Kind, Size, Crtp, Arg, Parameters...>& arg, const StaticVectorizer<Kind, Size, Crtp, Args, Parameters...>&... args) const
{
    Crtp<Return, Parameters...> result(*this);
    for (unsigned int i = 0; i < Size; ++i) {
        Vectorizer::get(result, i) = f(Vectorizer::get(result, i), Vectorizer::get(arg, i));
    }
    return result.template combine<Return>(f, args...);
}

// Unique combination by a function object
/// \brief          Unique combination by a function object.
/// \details        Does nothing.
/// \tparam         Return Return type.
/// \tparam         Function (Function type : <tt>Return(Return, Type)</tt>.)
/// \return         Copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename Return, class Function> 
inline Crtp<Return, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::combine(Function&&) const
{
    return Crtp<Return, Parameters...>(*this);
}
// -------------------------------------------------------------------------- //



// ---------------------------------- COUNT --------------------------------- //
// Count values
/// \brief          Count values.
/// \details        Count the number of valid comparisons with the provided
///                 reference in the specified region. Note that before any 
///                 comparison, the values in the container are casted to the
///                 reference data type. With no argument, this function counts 
///                 the number of non-null values (true) in the container.
/// \tparam         Reference (Reference type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      r Reference for comparison : value or vectorized container.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the valid counts.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Reference, class Mask, class> 
inline unsigned int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::count(const Reference& r, const Mask& bitmask) const
{
    unsigned int result = 0;
    Vectorizer::check(*this, r);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result += static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i));
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            result += (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i)));
        }
    }
    return result;
}

// Count with predicate
/// \brief          Count with predicate.
/// \details        Count the numbers of times the predicate is true in the
///                 specified region. 
/// \tparam         Function (Function type : <tt>bool(Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Predicate <tt>bool(Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the valid counts.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class> 
inline unsigned int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::count(Function&& f, const Mask& bitmask) const
{
    unsigned int result = 0;
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            result += static_cast<bool>(f(Vectorizer::get(*this, i)));
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            result += (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(f(Vectorizer::get(*this, i))));
        }
    }
    return result;
}

// All values equal
/// \brief          All values equal.
/// \details        Checks if the comparison with the provided reference returns
///                 true for all elements in the specified region. Note that 
///                 before any comparison, the values in the container are 
///                 casted to the reference data type. With no argument, this
///                 function returns true if the whole contents is non-null 
///                 (true). It returns true if the container is empty.
/// \tparam         Reference (Reference type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      r Reference for comparison : value or vectorized container.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the boolean result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Reference, class Mask, class> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::all(const Reference& r, const Mask& bitmask) const
{
    bool result = true;
    Vectorizer::check(*this, r);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            if (!static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i))) {
                result = false;
                break;
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && !static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i))) {
                result = false;
                break;
            }
        }
    }
    return result;
}

// All values satisfying the predicate
/// \brief          All values satisfying the predicate.
/// \details        Checks if the unary predicate returns true for all elements
///                 in the specified region. It returns true if the container is 
///                 empty.
/// \tparam         Function (Function type : <tt>bool(Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Predicate <tt>bool(Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the boolean result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::all(Function&& f, const Mask& bitmask) const
{
    bool result = true;
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            if (!static_cast<bool>(f(Vectorizer::get(*this, i)))) {
                result = false;
                break;
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && !static_cast<bool>(f(Vectorizer::get(*this, i)))) {
                result = false;
                break;
            }
        }
    }
    return result;
}

// Any value equal
/// \brief          Any value equal.
/// \details        Checks if the comparison with the provided reference returns 
///                 true for any element in the specified region. Note that 
///                 before any comparison, the values in the container are 
///                 casted to the reference data type. With no argument, this
///                 function returns true if the whole contents is non-null 
///                 (true). It returns false if the container is empty.
/// \tparam         Reference (Reference type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      r Reference for comparison : value or vectorized container.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the boolean result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Reference, class Mask, class> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::any(const Reference& r, const Mask& bitmask) const
{
    bool result = false;
    Vectorizer::check(*this, r);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i))) {
                result = true;
                break;
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i))) {
                result = true;
                break;
            }
        }
    }
    return result;
}

// Any value satisfying the predicate
/// \brief          Any value satisfying the predicate.
/// \details        Checks if the unary predicate returns true for any element
///                 in the specified region. It returns false if the container 
///                 is empty.
/// \tparam         Function (Function type : <tt>bool(Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Predicate <tt>bool(Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the boolean result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::any(Function&& f, const Mask& bitmask) const
{
    bool result = false;
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(f(Vectorizer::get(*this, i)))) {
                result = true;
                break;
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(f(Vectorizer::get(*this, i)))) {
                result = true;
                break;
            }
        }
    }
    return result;
}

// No value equal
/// \brief          No value equal.
/// \details        Checks if the comparison with the provided reference returns 
///                 true for no element in the specified region. Note that 
///                 before any comparison, the values in the container are 
///                 casted to the reference data type. With no argument, this
///                 function returns true if the whole contents is null (false).
///                 It returns true if the container is empty.
/// \tparam         Reference (Reference type.)
/// \tparam         Mask (Mask type.)
/// \param[in]      r Reference for comparison : value or vectorized container.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the boolean result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Reference, class Mask, class> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::none(const Reference& r, const Mask& bitmask) const
{
    bool result = true;
    Vectorizer::check(*this, r);
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i))) {
                result = false;
                break;
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(static_cast<typename std::remove_reference<decltype(Vectorizer::get(r, i))>::type>(Vectorizer::get(*this, i)) == Vectorizer::get(r, i))) {
                result = false;
                break;
            }
        }
    }
    return result;
}

// No value satisfying the predicate
/// \brief          No value satisfying the predicate.
/// \details        Checks if the unary predicate returns true for no element
///                 in the specified region. It returns true if the container 
///                 is empty.
/// \tparam         Function (Function type : <tt>bool(Type)</tt>.)
/// \tparam         Mask (Mask type.)
/// \param[in]      f Predicate <tt>bool(Type)</tt>.
/// \param[in]      bitmask Boolean mask.
/// \return         Copy of the boolean result.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Function, class Mask, class> 
inline bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::none(Function&& f, const Mask& bitmask) const
{
    bool result = true;
    Vectorizer::check(*this, bitmask);
    if (!std::is_base_of<Vectorizer, Mask>::value) {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(f(Vectorizer::get(*this, i)))) {
                result = false;
                break;
            }
        }
    } else {
        for (unsigned int i = 0; i < Size; ++i) {
            if (static_cast<bool>(Vectorizer::get(bitmask, i)) && static_cast<bool>(f(Vectorizer::get(*this, i)))) {
                result = false;
                break;
            }
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Adds each element to the stream using the <tt>fill()</tt>
///                 character as a separator.
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side container.
/// \return         Output stream.
template <typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters>
std::ostream& operator<<(std::ostream& lhs, const StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    const char separator = lhs.fill();
    for (int i = 0; i < std::integral_constant<int, static_cast<int>(SelfSize)-1>::value; ++i) {
        lhs<<Vectorizer::get(rhs, i)<<separator; 
    }
    return lhs<<Vectorizer::get(rhs, std::integral_constant<int, static_cast<int>(SelfSize)-1>::value);
}

// Input stream operator
/// \brief          Input stream operator.
/// \details        Fills each element from the stream.
/// \tparam         SelfKind (Kind of arguments.)
/// \tparam         SelfSize (Number of elements.)
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfType (Data type.)
/// \tparam         SelfParameters (List of parameters.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side container.
/// \return         Input stream.
template <typename SelfKind, unsigned int SelfSize, template <typename, SelfKind...> class SelfCrtp, typename SelfType, SelfKind... SelfParameters> 
std::istream& operator>>(std::istream& lhs, StaticVectorizer<SelfKind, SelfSize, SelfCrtp, SelfType, SelfParameters...>& rhs)
{
    for (int i = 0; i < std::integral_constant<int, static_cast<int>(SelfSize)-1>::value; ++i) {
        lhs>>Vectorizer::get(rhs, i);
    }
    return lhs>>Vectorizer::get(rhs, std::integral_constant<int, static_cast<int>(SelfSize)-1>::value);
}
// -------------------------------------------------------------------------- //



// -------------------------- STATIC VECTORIZATION -------------------------- //
// Get the size of the container
/// \brief          Get the size of the container.
/// \details        Returns the current number of elements. This function is 
///                 required by the vectorization mechanism. 
/// \return         Copy of the size.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr unsigned int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::size()
{
    return Size;
}

// Get whether the container has a constant size
/// \brief          Get whether the container has a constant size.
/// \details        Returns true if the container has a constant size, false
///                 otherwise. This function is required by the vectorization 
///                 mechanism. 
/// \return         Copy of true.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::constant()
{
    return true;
}

// Get whether the container has a boolean type
/// \brief          Get whether the container has a boolean type.
/// \details        Returns true if the container has a boolean type, false
///                 otherwise. This function is required by the vectorization 
///                 mechanism. 
/// \return         Copy of true if the container has a boolean type.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::boolean()
{
    return std::is_same<bool, Type>::value;
}

// Get the template parameters
/// \brief          Get the template parameters.
/// \details        Returns an array containing the template parameters.
/// \return         Copy of an array of parameters.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr std::array<Kind, sizeof...(Parameters)> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::parameters()
{
    return std::array<Kind, sizeof...(Parameters)>{{Parameters...}};
}

// Get the data type
/// \brief          Get the data type.
/// \details        Returns a copy of the default value of the data type.
/// \return         Copy of the default value of the data type.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Type StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::type()
{
    return Type();
}
// -------------------------------------------------------------------------- //



// ---------------------------------- SIZE ---------------------------------- //  
// Get whether the container is empty
/// \brief          Get whether the container is empty.
/// \details        Returns the result of the comparison between the size and 
///                 zero.
/// \return         Copy of the result of the test.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr bool StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::empty()
{
    return (Size == 0);
}

// Get the capacity of the container
/// \brief          Get the capacity of the container.
/// \details        Returns the capacity of the container, which is equal to its
///                 size.
/// \return         Copy of the capacity.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr unsigned int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::capacity()
{
    return Size;
}

// Get the size of the data type
/// \brief          Get the size of the data type.
/// \details        Returns the size of the element type in bytes.
/// \return         Copy of the size of the data type.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr unsigned int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::tbytes()
{
    return sizeof(Type);
}

// Get the size in bytes
/// \brief          Get the size in bytes.
/// \details        Returns the total size in bytes which is the number of 
///                 elements multiplied by the size in bytes of an element. 
///                 This does not take into account alignment bytes : to get the
///                 real memory imprint use the <tt>sizeof()</tt> function.
/// \return         Copy of the size in bytes.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
constexpr unsigned long long int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::bytes()
{
    return Size*sizeof(Type);
}

// Get the maximum available space
/// \brief          Get the maximum available space.
/// \details        Returns the <tt>max_size()</tt> of a <tt>std::vector</tt>
///                 of the same type.
/// \return         Copy of the maximum size.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline unsigned long long int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::space()
{
    return std::vector<Type>().max_size();
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Default mask creation
/// \brief          Default mask creation.
/// \details        Creates a mask from a boolean value, which is true per 
///                 default.
/// \param[in]      value Boolean value used to create the mask.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::mask(const bool value)
{
    return Crtp<bool, Parameters...>(value);
}

// Standard boolean container mask creation
/// \brief          Standard boolean container mask creation.
/// \details        Creates a mask from a standard boolean container, which can
///                 be a <tt>std::bitset</tt> or a <tt>std::vector<bool></tt>.
/// \tparam         Container (Boolean container type.)
/// \param[in]      container Boolean container.
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class Container, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::mask(const Container& container)
{
    const unsigned int lim = (container.size() > Size) ? (Size) : (container.size());
    Crtp<bool, Parameters...> result;
    for (unsigned int i = 0; i < lim; ++i) {
        Vectorizer::get(result, i) = container[i];
    }
    return result;
}

// Initializer list mask creation
/// \brief          Initializer list mask creation.
/// \details        Creates a mask from an initializer list by calling the
///                 associated constructor.
/// \tparam         OtherType (Other data type.)
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      source Source of the copy. 
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Boolean copy.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <typename OtherType, class... Misc, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::mask(const std::initializer_list<OtherType>& source, const Misc&... misc)
{
    return Crtp<bool, Parameters...>(source, misc...);
}

// Generic mask creation
/// \brief          Generic mask creation.
/// \details        Creates a mask from generic arguments by calling the 
///                 associated constructor.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.  
/// \return         Boolean copy. 
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
template <class... Misc, class> 
inline Crtp<bool, Parameters...> StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::mask(const Misc&... misc)
{
    Crtp<bool, Parameters...> result;
    return Crtp<bool, Parameters...>(misc...);
}
// -------------------------------------------------------------------------- //


        
// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of StaticVectorizer.
/// \return         0 if no error.
template <typename Kind, unsigned int Size, template <typename, Kind...> class Crtp, typename Type, Kind... Parameters>
int StaticVectorizer<Kind, Size, Crtp, Type, Parameters...>::example()
{
    std::cout<<"BEGIN = StaticVectorizer::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = StaticVectorizer::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = StaticVectorizer::example()"<<std::endl;
    return 1;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // STATICVECTORIZER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
