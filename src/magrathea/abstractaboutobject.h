/* *************************** ABSTRACTABOUTOBJECT ************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA-PATHFINDER 
// TITLE :          AbstractAboutObject 
// DESCRIPTION :    Tuple abstraction of generic about object 
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           abstractaboutobject.h 
/// \brief          Tuple abstraction of generic about object 
/// \author         Vincent Reverdy (vince.rev@gmail.com) 
/// \date           2012-2013 
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABSTRACTABOUTOBJECT_H_INCLUDED
#define ABSTRACTABOUTOBJECT_H_INCLUDED
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
// Misc
namespace magrathea {
//--------------------------------------------------------------------------- //



//---------------------------------- CLASS ---------------------------------- //
// Tuple abstraction of generic about object 
/// \brief          Tuple abstraction of generic about object. 
/// \details        This class is an abstraction of an about object with 
///                 common functions, like assignment, copy, getters and 
///                 setters. An about object is basically a standard structure 
///                 with additional features in order to be able to group 
///                 information about a software, its authors, its license and 
///                 other information. 
/// \tparam         Crtp Derived CRTP class.
/// \tparam         Types Variadic list of components types.
template <class Crtp, class... Types>
class AbstractAboutObject
{
    // Protected lifecycle 
    /// \name           Protected lifecycle 
    //@{
    protected: 
        inline ~AbstractAboutObject(); 
    //@}

    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        inline AbstractAboutObject(); 
        template <class OtherCrtp, class... OtherTypes> explicit inline AbstractAboutObject(const AbstractAboutObject<OtherCrtp, OtherTypes...>& source); 
        template <class... OtherTypes> explicit inline AbstractAboutObject(const std::tuple<OtherTypes...>& source); 
        template <class... OtherTypes, class = typename std::enable_if<(sizeof...(OtherTypes) != 0) && (std::is_constructible<typename std::tuple_element<0, typename std::conditional<sizeof...(Types) != 0, std::tuple<Types...>, std::tuple<std::true_type> >::type>::type, typename std::tuple_element<0, typename std::conditional<sizeof...(OtherTypes) != 0, std::tuple<OtherTypes...>, std::tuple<std::true_type> >::type>::type>::value)>::type> explicit inline AbstractAboutObject(const OtherTypes&... source); 
    //@}

    // Operators 
    /// \name           Operators 
    //@{
    public: 
        inline Crtp& operator=(const AbstractAboutObject<Crtp, Types...>& rhs); 
        template <class OtherCrtp, class... OtherTypes> inline Crtp& operator=(const AbstractAboutObject<OtherCrtp, OtherTypes...>& rhs); 
        template <class... OtherTypes> inline Crtp& operator=(const std::tuple<OtherTypes...>& rhs); 
        template <class OtherCrtp, class... OtherTypes> inline bool operator==(const AbstractAboutObject<OtherCrtp, OtherTypes...>& rhs) const; 
        template <class OtherCrtp, class... OtherTypes> inline bool operator!=(const AbstractAboutObject<OtherCrtp, OtherTypes...>& rhs) const; 
    //@}

    // Assignment 
    /// \name           Assignment 
    //@{
    public: 
        inline Crtp& assign(); 
        inline Crtp& assign(const AbstractAboutObject<Crtp, Types...>& source); 
        template <class OtherCrtp, class... OtherTypes> inline Crtp& assign(const AbstractAboutObject<OtherCrtp, OtherTypes...>& source); 
        template <class... OtherTypes> inline Crtp& assign(const std::tuple<OtherTypes...>& source); 
        template <class... OtherTypes, class = typename std::enable_if<(sizeof...(OtherTypes) != 0) && (std::is_constructible<typename std::tuple_element<0, typename std::conditional<sizeof...(Types) != 0, std::tuple<Types...>, std::tuple<std::true_type> >::type>::type, typename std::tuple_element<0, typename std::conditional<sizeof...(OtherTypes) != 0, std::tuple<OtherTypes...>, std::tuple<std::true_type> >::type>::type>::value)>::type> inline Crtp& assign(const OtherTypes&... source); 
    //@}

    // Management 
    /// \name           Management 
    //@{
    public: 
        inline Crtp& nullify(); 
        inline Crtp copy() const; 
        template <class OtherCrtp = Crtp, class = typename std::enable_if<std::is_constructible<OtherCrtp, Crtp>::value>::type> inline OtherCrtp cast() const; 
    //@}

    // Data 
    /// \name           Data 
    //@{
    public: 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type> inline Type& data(Dummy...); 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type> inline const Type& data(Dummy...) const; 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::value>::type> inline Crtp& data(const Type& value); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::value>::type> inline Type& data(Dummy...); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::value>::type> inline const Type& data(Dummy...) const; 
        template <unsigned int Index, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::type>::value>::type> inline Crtp& data(const Type& value); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::type>::value>::type> inline Type& data(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>::type>::value>::type> inline const Type& data(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::type>::type>::value>::type> inline Crtp& data(const Type& value); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& data(const Subscript subscript, Dummy...); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Types...>, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& data(const Subscript subscript, Dummy...) const; 
        template <unsigned int Index, typename Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Types...>, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& data(const Subscript subscript, const Type& value); 
    //@}

    // Getters 
    /// \name           Getters 
    //@{
    public: 
        template <class Type = std::tuple<Types...>, class = typename std::enable_if<std::is_convertible<Type, std::tuple<Types...> >::value>::type> inline const Type& get() const; 
        template <unsigned int Index, class Type = typename std::tuple_element<Index, std::tuple<Types...> >::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Index, std::tuple<Types...> >::type>::value)>::type> inline const Type& get(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Types...> >::type>::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Types...> >::type>::type>::value)>::type> inline const Type& get(const Dummy&...) const; 
        template <unsigned int Index, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Types...> >::type>()[0])>::type>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Types...> >::type>()[0])>::type>::type>::value>::type> inline const Type& get(const unsigned int subscript) const; 
    //@}

    // Setters 
    /// \name           Setters 
    //@{
    public: 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, std::tuple<Types...> >::value>::type> inline Crtp& set(const Type& value); 
        template <unsigned int Index, class Type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Index, std::tuple<Types...> >::type>::value)>::type> inline Crtp& set(const Type& value, Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Types...> >::type>::type>::value)>::type> inline Crtp& set(const Type& value, const Dummy&...); 
        template <unsigned int Index, class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Types...> >::type>()[0])>::type>::type>::value>::type> inline Crtp& set(const unsigned int subscript, const Type& value); 
    //@}

    // Stream 
    /// \name           Stream 
    //@{
    public: 
        template <class SelfCrtp, class... SelfTypes> friend std::ostream& operator<<(std::ostream& lhs, const AbstractAboutObject<SelfCrtp, SelfTypes...>& rhs); 
    //@}

    // Types 
    /// \name           Types 
    //@{
    public: 
        template <class Type = std::tuple<Types...>, class = typename std::enable_if<std::is_convertible<Type, std::tuple<Types...> >::value>::type> static constexpr Type type(); 
        template <unsigned int Index, class Type = typename std::tuple_element<Index, std::tuple<Types...> >::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Index, std::tuple<Types...> >::type>::value)>::type> static constexpr Type type(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Types...> >::type>::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Types...> >::type>::type>::value)>::type> static constexpr Type type(const Dummy&...); 
        template <unsigned int Index, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Types...> >::type>()[0])>::type>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Types...> >::type>()[0])>::type>::type>::value>::type> static constexpr Type type(const unsigned int subscript); 
    //@}

    // Properties 
    /// \name           Properties 
    //@{
    public: 
        static constexpr unsigned int types(); 
    //@}

    // Helpers 
    /// \name           Helpers 
    //@{
    public: 
        template <class Output> static constexpr Output transmute(); 
        template <class Output, class Input, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<Input>::type>::type, typename std::remove_cv<typename std::remove_reference<Output>::type>::type>::value>::type> static constexpr Input transmute(Input&& input); 
        template <class Output, class... Input, class = typename std::enable_if<(sizeof...(Input) <= std::tuple_size<Output>::value) && (sizeof...(Input) != std::tuple_size<Output>::value)>::type> static constexpr Output transmute(const Input&... input); 
        template <class Output, class... Input, class = typename std::enable_if<(sizeof...(Input) == std::tuple_size<Output>::value) && (sizeof...(Input) != 0)>::type> static constexpr Output transmute(Input&&... input); 
        template <class Input> static constexpr bool printable(typename std::enable_if<!std::is_void<decltype(std::declval<std::ostream&>()<<std::declval<Input>())>::value>::type* = nullptr); 
        template <class Input, class... Dummy, class = typename std::enable_if<sizeof...(Dummy) == 0>::type> static constexpr bool printable(Dummy...); 
        template <class Input, class = typename std::enable_if<printable<typename std::remove_cv<typename std::remove_reference<Input>::type>::type>()>::type> static inline bool print(std::ostream& stream, Input&& input); 
        template <class Input, class = typename std::enable_if<(!printable<Input>()) && (!std::is_void<decltype(std::declval<Input>()[0])>::value) && (!std::is_void<decltype(std::declval<Input>().size())>::value)>::type> static inline bool print(std::ostream& stream, const Input& input); 
        template <unsigned int Current = 0, class Input, class... Dummy, class = typename std::enable_if<(!printable<Input>()) && (Current <= std::tuple_size<Input>::value) && (Current != std::tuple_size<Input>::value) && (sizeof...(Dummy) == 0)>::type> static inline bool print(std::ostream& stream, const Input& input, Dummy...); 
        template <unsigned int Current = 0, class... Dummy, class = typename std::enable_if<sizeof...(Dummy) <= 1>::type> static inline bool print(std::ostream& stream, const Dummy&...); 
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
        std::tuple<Types...> _data;                                             ///< Internal tuple container. 
    //@}
};
//--------------------------------------------------------------------------- //



//--------------------------- PROTECTED LIFECYCLE --------------------------- //
// Protected destructor 
/// \brief          Protected destructor. 
/// \details        Avoids direct instantiation of the class, and only allows 
///                 it through its derived children. 
template <class Crtp, class... Types>
inline AbstractAboutObject<Crtp, Types...>::~AbstractAboutObject()
= default;
//--------------------------------------------------------------------------- //



//-------------------------------- LIFECYCLE -------------------------------- //
// Implicit empty constructor 
/// \brief          Implicit empty constructor. 
/// \details        Provides an implicit construction of an object initialized 
///                 to its default value. 
template <class Crtp, class... Types>
inline AbstractAboutObject<Crtp, Types...>::AbstractAboutObject()
: _data(std::tuple<Types...>())
{
    ;
}

// Explicit conversion constructor 
/// \brief          Explicit conversion constructor. 
/// \details        Provides an explicit construction from another type of 
///                 object. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherTypes (Other variadic list of components types.)
/// \param[in]      source Source of the copy.
template <class Crtp, class... Types>
template <class OtherCrtp, class... OtherTypes>
inline AbstractAboutObject<Crtp, Types...>::AbstractAboutObject(const AbstractAboutObject<OtherCrtp, OtherTypes...>& source)
: _data(transmute<std::tuple<Types...> >(source.data()))
{
    ;
}

// Explicit data constructor 
/// \brief          Explicit data constructor. 
/// \details        Provides an explicit construction from data. 
/// \tparam         OtherTypes (Other variadic list of object property types.)
/// \param[in]      source Source of the copy.
template <class Crtp, class... Types>
template <class... OtherTypes>
inline AbstractAboutObject<Crtp, Types...>::AbstractAboutObject(const std::tuple<OtherTypes...>& source)
: _data(transmute<std::tuple<Types...> >(source))
{
    ;
}

// Explicit components constructor 
/// \brief          Explicit components constructor. 
/// \details        Provides an explicit construction from components. 
/// \tparam         OtherTypes (Other variadic list of object property types.)
/// \param[in]      source Source of the copy.
template <class Crtp, class... Types>
template <class... OtherTypes, class>
inline AbstractAboutObject<Crtp, Types...>::AbstractAboutObject(const OtherTypes&... source)
: _data(transmute<std::tuple<Types...> >(source...))
{
    ;
}
//--------------------------------------------------------------------------- //



//-------------------------------- OPERATORS -------------------------------- //
// Copy assignment operator 
/// \brief          Copy assignment operator. 
/// \details        Assigns contents from the same type of object. 
/// \param[in]      rhs Right-hand side.
/// \return         Self reference. 
template <class Crtp, class... Types>
inline Crtp& AbstractAboutObject<Crtp, Types...>::operator=(const AbstractAboutObject<Crtp, Types...>& rhs)
{
    _data = rhs._data;
    return static_cast<Crtp&>(*this);
}

// Conversion assignment operator 
/// \brief          Conversion assignment operator. 
/// \details        Assigns contents from another type of object. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherTypes (Other variadic list of components types.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class OtherCrtp, class... OtherTypes>
inline Crtp& AbstractAboutObject<Crtp, Types...>::operator=(const AbstractAboutObject<OtherCrtp, OtherTypes...>& rhs)
{
    _data = transmute<std::tuple<Types...> >(rhs.data());
    return static_cast<Crtp&>(*this);
}

// Data assignment operator 
/// \brief          Data assignment operator. 
/// \details        Assigns contents from data. 
/// \tparam         OtherTypes (Other variadic list of object property types.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class... OtherTypes>
inline Crtp& AbstractAboutObject<Crtp, Types...>::operator=(const std::tuple<OtherTypes...>& rhs)
{
    _data = transmute<std::tuple<Types...> >(rhs);
    return static_cast<Crtp&>(*this);
}

// Equal to 
/// \brief          Equal to. 
/// \details        Compares for equality and returns true if the contents is 
///                 equal. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherTypes (Other variadic list of components types.)
/// \param[in]      rhs Right-hand side.
/// \return         True if equal, false if not equal. 
template <class Crtp, class... Types>
template <class OtherCrtp, class... OtherTypes>
inline bool AbstractAboutObject<Crtp, Types...>::operator==(const AbstractAboutObject<OtherCrtp, OtherTypes...>& rhs) const
{
    return (_data == rhs.data());
}

// Not equal to 
/// \brief          Not equal to. 
/// \details        Compares for difference and returns true if the contents 
///                 is different. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherTypes (Other variadic list of components types.)
/// \param[in]      rhs Right-hand side.
/// \return         True if not equal, false if equal. 
template <class Crtp, class... Types>
template <class OtherCrtp, class... OtherTypes>
inline bool AbstractAboutObject<Crtp, Types...>::operator!=(const AbstractAboutObject<OtherCrtp, OtherTypes...>& rhs) const
{
    return (_data != rhs.data());
}
//--------------------------------------------------------------------------- //



//-------------------------------- ASSIGNMENT ------------------------------- //
// Empty assignment 
/// \brief          Empty assignment. 
/// \details        Assigns contents from an object initialized to its default 
///                 value. 
/// \return         Self reference. 
template <class Crtp, class... Types>
inline Crtp& AbstractAboutObject<Crtp, Types...>::assign()
{
    _data = std::tuple<Types...>();
    return static_cast<Crtp&>(*this);;
}

// Copy assignment 
/// \brief          Copy assignment. 
/// \details        Assigns contents from the same type of object. 
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <class Crtp, class... Types>
inline Crtp& AbstractAboutObject<Crtp, Types...>::assign(const AbstractAboutObject<Crtp, Types...>& source)
{
    _data = source._data;
    return static_cast<Crtp&>(*this);;
}

// Conversion assignment 
/// \brief          Conversion assignment. 
/// \details        Assigns contents from another type of object. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherTypes (Other variadic list of components types.)
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class OtherCrtp, class... OtherTypes>
inline Crtp& AbstractAboutObject<Crtp, Types...>::assign(const AbstractAboutObject<OtherCrtp, OtherTypes...>& source)
{
    _data = transmute<std::tuple<Types...> >(source.data());
    return static_cast<Crtp&>(*this);
}

// Data assignment 
/// \brief          Data assignment. 
/// \details        Assigns contents from data. 
/// \tparam         OtherTypes (Other variadic list of object property types.)
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class... OtherTypes>
inline Crtp& AbstractAboutObject<Crtp, Types...>::assign(const std::tuple<OtherTypes...>& source)
{
    _data = transmute<std::tuple<Types...> >(source);
    return static_cast<Crtp&>(*this);
}

// Components assignment 
/// \brief          Components assignment. 
/// \details        Assigns contents from components. 
/// \tparam         OtherTypes (Other variadic list of object property types.)
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class... OtherTypes, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::assign(const OtherTypes&... source)
{
    _data = transmute<std::tuple<Types...> >(source...);
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//-------------------------------- MANAGEMENT ------------------------------- //
// Nullify 
/// \brief          Nullify. 
/// \details        Resets all data members to their default values. 
/// \return         Self reference. 
template <class Crtp, class... Types>
inline Crtp& AbstractAboutObject<Crtp, Types...>::nullify()
{
    _data = std::tuple<Types...>();
    return static_cast<Crtp&>(*this);
}

// Copy 
/// \brief          Copy. 
/// \details        Generates a copy of the object. 
/// \return         Copy. 
template <class Crtp, class... Types>
inline Crtp AbstractAboutObject<Crtp, Types...>::copy() const
{
    return static_cast<const Crtp&>(*this);
}

// Cast 
/// \brief          Cast. 
/// \details        Casts contents to another object type. 
/// \tparam         OtherCrtp Other derived CRTP class.
/// \return         Casted copy. 
template <class Crtp, class... Types>
template <class OtherCrtp, class>
inline OtherCrtp AbstractAboutObject<Crtp, Types...>::cast() const
{
    return OtherCrtp(*this);
}
//--------------------------------------------------------------------------- //



//----------------------------------- DATA ---------------------------------- //
// Unified data access 
/// \brief          Unified data access. 
/// \details        Provides a direct access to the data. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Data std::tuple<Types...> type.)
/// \return         Reference to the data. 
template <class Crtp, class... Types>
template <class... Dummy, class Type, class, class>
inline Type& AbstractAboutObject<Crtp, Types...>::data(Dummy...)
{
    return _data;
}

// Unified data getter 
/// \brief          Unified data getter. 
/// \details        Gets the data. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Data std::tuple<Types...> type.)
/// \return         Immutable reference to the data. 
template <class Crtp, class... Types>
template <class... Dummy, class Type, class, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::data(Dummy...) const
{
    return _data;
}

// Unified data setter 
/// \brief          Unified data setter. 
/// \details        Sets the data. 
/// \param[in]      value Data value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class Type, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::data(const Type& value)
{
    _data = value;
    return static_cast<Crtp&>(*this);
}

// Unified data component access 
/// \brief          Unified data component access. 
/// \details        Provides a direct access to the specified component of the 
///                 data. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Reference to the component of the data. 
template <class Crtp, class... Types>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractAboutObject<Crtp, Types...>::data(Dummy...)
{
    return std::get<Index>(_data);
}

// Unified data component getter 
/// \brief          Unified data component getter. 
/// \details        Gets the specified component of the data. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Immutable reference to the component of the data. 
template <class Crtp, class... Types>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::data(Dummy...) const
{
    return std::get<Index>(_data);
}

// Unified data component setter 
/// \brief          Unified data component setter. 
/// \details        Sets the specified component of the data. 
/// \tparam         Index Index of the component.
/// \param[in]      value Component value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::data(const Type& value)
{
    std::get<Index>(_data) = value;
    return static_cast<Crtp&>(*this);
}

// Unified data inner component access 
/// \brief          Unified data inner component access. 
/// \details        Provides a direct access to the specified inner component 
///                 of the specified component of the data. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Reference to the inner component of the data. 
template <class Crtp, class... Types>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractAboutObject<Crtp, Types...>::data(Dummy...)
{
    return std::get<Subscript>(std::get<Index>(_data));
}

// Unified data inner component getter 
/// \brief          Unified data inner component getter. 
/// \details        Gets the specified inner component of the specified 
///                 component of the data. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Immutable reference to the inner component of the data. 
template <class Crtp, class... Types>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::data(Dummy...) const
{
    return std::get<Subscript>(std::get<Index>(_data));
}

// Unified data inner component setter 
/// \brief          Unified data inner component setter. 
/// \details        Sets the specified inner component of the specified 
///                 component of the data. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \param[in]      value Inner component value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <unsigned int Index, unsigned int Subscript, class Type, class, class, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::data(const Type& value)
{
    std::get<Subscript>(std::get<Index>(_data)) = value;
    return static_cast<Crtp&>(*this);
}

// Unified data component element access 
/// \brief          Unified data component element access. 
/// \details        Provides a direct access to the specified element of the 
///                 specified component of the data. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Reference to the component element of the data. 
template <class Crtp, class... Types>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractAboutObject<Crtp, Types...>::data(const Subscript subscript, Dummy...)
{
    return std::get<Index>(_data)[subscript];
}

// Unified data component element getter 
/// \brief          Unified data component element getter. 
/// \details        Gets the specified element of the specified component of 
///                 the data. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Immutable reference to the component element of the data. 
template <class Crtp, class... Types>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::data(const Subscript subscript, Dummy...) const
{
    return std::get<Index>(_data)[subscript];
}

// Unified data component element setter 
/// \brief          Unified data component element setter. 
/// \details        Sets the specified element of the specified component of 
///                 the data. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \param[in]      subscript Subscript of the component element.
/// \param[in]      value Component element value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <unsigned int Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::data(const Subscript subscript, const Type& value)
{
    std::get<Index>(_data)[subscript] = value;
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//--------------------------------- GETTERS --------------------------------- //
// Data getter 
/// \brief          Data getter. 
/// \details        Gets the underlying tuple. 
/// \tparam         Type (Underlying tuple type.)
/// \return         Immutable reference to the tuple. 
template <class Crtp, class... Types>
template <class Type, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::get() const
{
    return _data;
}

// Component getter 
/// \brief          Component getter. 
/// \details        Gets the specified component of the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Type (Component type.)
/// \tparam         Dummy (Dummy types.)
/// \return         Immutable reference to the component. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class... Dummy, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::get(Dummy...) const
{
    return std::get<Index>(_data);
}

// Inner component getter 
/// \brief          Inner component getter. 
/// \details        Gets the specified inner component of the specified 
///                 component of the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Index of the inner component.
/// \tparam         Type (Inner component type.)
/// \tparam         Dummy (Dummy types.)
/// \return         Immutable reference to the inner component. 
template <class Crtp, class... Types>
template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::get(const Dummy&...) const
{
    return std::get<Subscript>(std::get<Index>(_data));
}

// Component element getter 
/// \brief          Component element getter. 
/// \details        Gets the specified element of the specified component of 
///                 the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Type (Element type.)
/// \param[in]      subscript Index of the element.
/// \return         Immutable reference to the component element. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class>
inline const Type& AbstractAboutObject<Crtp, Types...>::get(const unsigned int subscript) const
{
    return std::get<Index>(_data)[subscript];
}
//--------------------------------------------------------------------------- //



//--------------------------------- SETTERS --------------------------------- //
// Data setter 
/// \brief          Data setter. 
/// \details        Sets the underlying tuple. 
/// \tparam         Type (Tuple type.)
/// \param[in]      value Tuple value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <class Type, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::set(const Type& value)
{
    _data = value;
    return static_cast<Crtp&>(*this);
}

// Component setter 
/// \brief          Component setter. 
/// \details        Sets the specified component of the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Type (Component type.)
/// \tparam         Dummy (Dummy types.)
/// \param[in]      value Component value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class... Dummy, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::set(const Type& value, Dummy...)
{
    std::get<Index>(_data) = value;
    return static_cast<Crtp&>(*this);
}

// Inner component setter 
/// \brief          Inner component setter. 
/// \details        Sets the specified inner component of the specified 
///                 component of the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Index of the inner component.
/// \tparam         Type (Inner component type.)
/// \tparam         Dummy (Dummy types.)
/// \param[in]      value Inner component value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::set(const Type& value, const Dummy&...)
{
    std::get<Subscript>(std::get<Index>(_data)) = value;
    return static_cast<Crtp&>(*this);
}

// Component element setter 
/// \brief          Component element setter. 
/// \details        Sets the specified element of the specified component of 
///                 the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Type (Element type.)
/// \param[in]      subscript Index of the element.
/// \param[in]      value Element value.
/// \return         Self reference. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class>
inline Crtp& AbstractAboutObject<Crtp, Types...>::set(const unsigned int subscript, const Type& value)
{
    std::get<Index>(_data)[subscript] = value;
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//---------------------------------- STREAM --------------------------------- //
// Output stream operator 
/// \brief          Output stream operator. 
/// \details        Adds each element to the stream. 
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfTypes (Variadic list of components types.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side object.
/// \return         Output stream. 
template <class SelfCrtp, class... SelfTypes>
std::ostream& operator<<(std::ostream& lhs, const AbstractAboutObject<SelfCrtp, SelfTypes...>& rhs)
{
    rhs.print(lhs, rhs._data);
    return lhs;
}
//--------------------------------------------------------------------------- //



//---------------------------------- TYPES ---------------------------------- //
// Data type 
/// \brief          Data type. 
/// \details        Returns the default value of the type of the underlying 
///                 tuple. 
/// \tparam         Type (Underlying tuple type.)
/// \return         Default tuple value. 
template <class Crtp, class... Types>
template <class Type, class>
constexpr Type AbstractAboutObject<Crtp, Types...>::type()
{
    return Type();
}

// Component type 
/// \brief          Component type. 
/// \details        Returns the default value of the type of the specified 
///                 component of the underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Type (Component type.)
/// \tparam         Dummy (Dummy types.)
/// \return         Default component value. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class... Dummy, class>
constexpr Type AbstractAboutObject<Crtp, Types...>::type(Dummy...)
{
    return Type();
}

// Inner component type 
/// \brief          Inner component type. 
/// \details        Returns the default value of the type of the specified 
///                 inner component of the specified component of the 
///                 underlying tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Index of the inner component.
/// \tparam         Type (Inner component type.)
/// \tparam         Dummy (Dummy types.)
/// \return         Default inner component value. 
template <class Crtp, class... Types>
template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class>
constexpr Type AbstractAboutObject<Crtp, Types...>::type(const Dummy&...)
{
    return Type();
}

// Component element type 
/// \brief          Component element type. 
/// \details        Returns the default value of the type of the specified 
///                 element of the specified component of the underlying 
///                 tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Type (Element type.)
/// \param[in]      subscript Index of the element.
/// \return         Default component element value. 
template <class Crtp, class... Types>
template <unsigned int Index, class Type, class>
constexpr Type AbstractAboutObject<Crtp, Types...>::type(const unsigned int subscript)
{
    return (subscript != 0) ? (Type()) : (Type());
}
//--------------------------------------------------------------------------- //



//-------------------------------- PROPERTIES ------------------------------- //
// Total number of types 
/// \brief          Total number of types. 
/// \details        Counts the number of internal components. 
/// \return         Number of components. 
template <class Crtp, class... Types>
constexpr unsigned int AbstractAboutObject<Crtp, Types...>::types()
{
    return std::tuple_size<std::tuple<Types...> >::value;
}
//--------------------------------------------------------------------------- //



//--------------------------------- HELPERS --------------------------------- //
// Default transmutation 
/// \brief          Default transmutation. 
/// \details        Transmutes the input into the specified output type by 
///                 creating a new default value. 
/// \tparam         Output Output type.
/// \return         Default output value. 
template <class Crtp, class... Types>
template <class Output>
constexpr Output AbstractAboutObject<Crtp, Types...>::transmute()
{
    return Output();
}

// Standard transmutation 
/// \brief          Standard transmutation. 
/// \details        Transmutes the input into the specified output type by 
///                 returning the provided convertible value. 
/// \tparam         Output Output type.
/// \tparam         Input (Input type.)
/// \param[in]      input Input value.
/// \return         Forwarded input value. 
template <class Crtp, class... Types>
template <class Output, class Input, class>
constexpr Input AbstractAboutObject<Crtp, Types...>::transmute(Input&& input)
{
    return std::forward<Input>(input);
}

// Recursive transmutation 
/// \brief          Recursive transmutation. 
/// \details        Transmutes the input into the specified output type by 
///                 recursively creating new values to fill a tuple. 
/// \tparam         Output Output type.
/// \tparam         Input (Input types.)
/// \param[in]      input Input values.
/// \return         Forwarded tuple. 
template <class Crtp, class... Types>
template <class Output, class... Input, class>
constexpr Output AbstractAboutObject<Crtp, Types...>::transmute(const Input&... input)
{
    return transmute<Output>(input..., typename std::tuple_element<sizeof...(Input), Output>::type());
}

// Forwarding transmutation 
/// \brief          Forwarding transmutation. 
/// \details        Transmutes the input into the specified output type by 
///                 forwarding the provided values as a tuple. 
/// \tparam         Output Output type.
/// \tparam         Input (Input types.)
/// \param[in]      input Input values.
/// \return         Forwarded tuple. 
template <class Crtp, class... Types>
template <class Output, class... Input, class>
constexpr Output AbstractAboutObject<Crtp, Types...>::transmute(Input&&... input)
{
    return transmute<Output>(std::forward_as_tuple(std::forward<Input>(input)...));
}

// Is printable 
/// \brief          Is printable. 
/// \details        True overload to check at compile time whether a type can 
///                 be inserted in an ostream. 
/// \tparam         Input Input type.
/// \return         True. 
template <class Crtp, class... Types>
template <class Input>
constexpr bool AbstractAboutObject<Crtp, Types...>::printable(typename std::enable_if<!std::is_void<decltype(std::declval<std::ostream&>()<<std::declval<Input>())>::value>::type*)
{
    return true;
}

// Is not printable 
/// \brief          Is not printable. 
/// \details        False overload to check at compile time whether a type 
///                 cannot be inserted in an ostream. 
/// \tparam         Input Input type.
/// \tparam         Dummy (Dummy types.)
/// \return         False. 
template <class Crtp, class... Types>
template <class Input, class... Dummy, class>
constexpr bool AbstractAboutObject<Crtp, Types...>::printable(Dummy...)
{
    return false;
}

// Standard print 
/// \brief          Standard print. 
/// \details        Prints a printable object to the specified stream. 
/// \tparam         Input (Printable type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      input Printable value.
/// \return         Stream status. 
template <class Crtp, class... Types>
template <class Input, class>
inline bool AbstractAboutObject<Crtp, Types...>::print(std::ostream& stream, Input&& input)
{
    return (stream<<std::forward<Input>(input)).good();
}

// Iterative print 
/// \brief          Iterative print. 
/// \details        Prints each element of an array to the specified stream. 
/// \tparam         Input (Array type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      input Array value.
/// \return         Stream status. 
template <class Crtp, class... Types>
template <class Input, class>
inline bool AbstractAboutObject<Crtp, Types...>::print(std::ostream& stream, const Input& input)
{
    const unsigned int length = input.size();
    const char separator = stream.fill();
    if (length > 0) {
        print(stream, input[0]);
        for (unsigned int i = 1; i < length; ++i) {
            print(stream<<separator, input[i]);
        }
    }
    return stream.good();
}

// Recursive print 
/// \brief          Recursive print. 
/// \details        Prints each element of a tuple to the specified stream. 
/// \tparam         Current (Current level of recursion.)
/// \tparam         Input (Tuple type.)
/// \tparam         Dummy (Dummy types.)
/// \param[in,out]  stream Output stream.
/// \param[in]      input Tuple value.
/// \return         Stream status. 
template <class Crtp, class... Types>
template <unsigned int Current, class Input, class... Dummy, class>
inline bool AbstractAboutObject<Crtp, Types...>::print(std::ostream& stream, const Input& input, Dummy...)
{
    print((Current > 0) ? (stream<<stream.fill()) : (stream), std::get<Current>(input));
    return print<Current+1>(stream, input);
}

// Ineffective print 
/// \brief          Ineffective print. 
/// \details        Prints nothing to the specified stream. 
/// \tparam         Current (Current level of recursion.)
/// \tparam         Dummy (Dummy types.)
/// \param[in,out]  stream Output stream.
/// \return         Stream status. 
template <class Crtp, class... Types>
template <unsigned int Current, class... Dummy, class>
inline bool AbstractAboutObject<Crtp, Types...>::print(std::ostream& stream, const Dummy&...)
{
    return stream.good();
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of AbstractAboutObject. 
/// \return         0 if no error. 
template <class Crtp, class... Types>
int AbstractAboutObject<Crtp, Types...>::example()
{
    std::cout<<"BEGIN = AbstractAboutObject::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = AbstractAboutObject::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = AbstractAboutObject::example()"<<std::endl;
    return 1;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABSTRACTABOUTOBJECT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

