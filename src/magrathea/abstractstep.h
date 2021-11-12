/* ****************************** ABSTRACTSTEP ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA-PATHFINDER 
// TITLE :          AbstractStep 
// DESCRIPTION :    Abstraction of an evolution step 
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           abstractstep.h 
/// \brief          Abstraction of an evolution step 
/// \author         Vincent Reverdy (vince.rev@gmail.com) 
/// \date           2012-2013 
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABSTRACTSTEP_H_INCLUDED
#define ABSTRACTSTEP_H_INCLUDED
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
// Abstraction of an evolution step 
/// \brief          Abstraction of an evolution step. 
/// \details        This class is an abstraction of an evolution step with 
///                 common functions, like assignment, copy, getters and 
///                 setters. A step object is basically a standard structure 
///                 with additional features. The internal behaviour is based 
///                 on three groups of quantities : <ul><li>id : a scalar id 
///                 of the step</li> <li>core : an array of basic 
///                 quantities</li> <li>extra : a tuple of derived 
///                 quantities</li></ul>. 
/// \tparam         Crtp Derived CRTP class.
/// \tparam         Scalar Scalar type of id.
/// \tparam         Array Array type of core quantities.
/// \tparam         Tuple Tuple type of extra quantities.
template <class Crtp, class Scalar, class Array, class Tuple>
class AbstractStep
{
    // Protected lifecycle 
    /// \name           Protected lifecycle 
    //@{
    protected: 
        inline ~AbstractStep(); 
    //@}

    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        inline AbstractStep(); 
        template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple> explicit inline AbstractStep(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& source); 
        template <class OtherScalar, class OtherArray = Array, class OtherTuple = Tuple, class = typename std::enable_if<(std::is_constructible<Scalar, OtherScalar>::value) && (std::is_constructible<Array, OtherArray>::value) && (std::is_constructible<Tuple, OtherTuple>::value)>::type> explicit inline AbstractStep(const OtherScalar& scalar, const OtherArray& array = OtherArray(), const OtherTuple& tuple = OtherTuple()); 
    //@}

    // Operators 
    /// \name           Operators 
    //@{
    public: 
        inline Crtp& operator=(const AbstractStep<Crtp, Scalar, Array, Tuple>& rhs); 
        template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple> inline Crtp& operator=(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& rhs); 
        template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple> inline bool operator==(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& rhs) const; 
        template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple> inline bool operator!=(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& rhs) const; 
    //@}

    // Assignment 
    /// \name           Assignment 
    //@{
    public: 
        inline Crtp& assign(); 
        inline Crtp& assign(const AbstractStep<Crtp, Scalar, Array, Tuple>& source); 
        template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple> inline Crtp& assign(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& source); 
        template <class OtherScalar, class OtherArray = Array, class OtherTuple = Tuple, class = typename std::enable_if<(std::is_constructible<Scalar, OtherScalar>::value) && (std::is_constructible<Array, OtherArray>::value) && (std::is_constructible<Tuple, OtherTuple>::value)>::type> inline Crtp& assign(const OtherScalar& scalar, const OtherArray& array = OtherArray(), const OtherTuple& tuple = OtherTuple()); 
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
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type> inline Type& data(Dummy...); 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type> inline const Type& data(Dummy...) const; 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type> inline Crtp& data(const Type& value); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::value>::type> inline Type& data(Dummy...); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::value>::type> inline const Type& data(Dummy...) const; 
        template <unsigned int Index, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::value>::type> inline Crtp& data(const Type& value); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::type>::value>::type> inline Type& data(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::type>::value>::type> inline const Type& data(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::type>::type>::value>::type> inline Crtp& data(const Type& value); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& data(const Subscript subscript, Dummy...); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, std::tuple<Scalar, Array, Tuple>, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& data(const Subscript subscript, Dummy...) const; 
        template <unsigned int Index, typename Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, std::tuple<Scalar, Array, Tuple>, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& data(const Subscript subscript, const Type& value); 
    //@}

    // Getters 
    /// \name           Getters 
    //@{
    public: 
        template <class Type = std::tuple<Scalar, Array, Tuple>, class = typename std::enable_if<std::is_convertible<Type, std::tuple<Scalar, Array, Tuple> >::value>::type> inline const Type& get() const; 
        template <unsigned int Index, class Type = typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::value)>::type> inline const Type& get(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::type>::value)>::type> inline const Type& get(const Dummy&...) const; 
        template <unsigned int Index, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>()[0])>::type>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>()[0])>::type>::type>::value>::type> inline const Type& get(const unsigned int subscript) const; 
    //@}

    // Setters 
    /// \name           Setters 
    //@{
    public: 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, std::tuple<Scalar, Array, Tuple> >::value>::type> inline Crtp& set(const Type& value); 
        template <unsigned int Index, class Type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::value)>::type> inline Crtp& set(const Type& value, Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::type>::value)>::type> inline Crtp& set(const Type& value, const Dummy&...); 
        template <unsigned int Index, class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>()[0])>::type>::type>::value>::type> inline Crtp& set(const unsigned int subscript, const Type& value); 
    //@}

    // Id 
    /// \name           Id 
    //@{
    public: 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type> inline Type& id(Dummy...); 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type> inline const Type& id(Dummy...) const; 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::value>::type> inline Crtp& id(const Type& value); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::value>::type> inline Type& id(Dummy...); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::value>::type> inline const Type& id(Dummy...) const; 
        template <unsigned int Index, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::type>::value>::type> inline Crtp& id(const Type& value); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::type>::value>::type> inline Type& id(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>::type>::value>::type> inline const Type& id(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::type>::type>::value>::type> inline Crtp& id(const Type& value); 
        template <typename Index, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type>::value>::type> inline Type& id(const Index index, Dummy...); 
        template <typename Index, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type>::value>::type> inline const Type& id(const Index index, Dummy...) const; 
        template <typename Index, class Type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>()[Index()])>::type>::type>::value>::type> inline Crtp& id(const Index index, const Type& value); 
        template <typename Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& id(const Index index, const Subscript subscript, Dummy...); 
        template <typename Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& id(const Index index, const Subscript subscript, Dummy...) const; 
        template <typename Index, typename Subscript, class Type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& id(const Index index, const Subscript subscript, const Type& value); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& id(const Subscript subscript, Dummy...); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Scalar, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& id(const Subscript subscript, Dummy...) const; 
        template <unsigned int Index, typename Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Scalar, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& id(const Subscript subscript, const Type& value); 
    //@}

    // Core 
    /// \name           Core 
    //@{
    public: 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type> inline Type& core(Dummy...); 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type> inline const Type& core(Dummy...) const; 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::value>::type> inline Crtp& core(const Type& value); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::value>::type> inline Type& core(Dummy...); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::value>::type> inline const Type& core(Dummy...) const; 
        template <unsigned int Index, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::type>::value>::type> inline Crtp& core(const Type& value); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::type>::value>::type> inline Type& core(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>::type>::value>::type> inline const Type& core(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::type>::type>::value>::type> inline Crtp& core(const Type& value); 
        template <typename Index, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type>::value>::type> inline Type& core(const Index index, Dummy...); 
        template <typename Index, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type>::value>::type> inline const Type& core(const Index index, Dummy...) const; 
        template <typename Index, class Type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<!std::is_void<Type>::value, Array, void>::type>()[Index()])>::type>::type>::value>::type> inline Crtp& core(const Index index, const Type& value); 
        template <typename Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& core(const Index index, const Subscript subscript, Dummy...); 
        template <typename Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& core(const Index index, const Subscript subscript, Dummy...) const; 
        template <typename Index, typename Subscript, class Type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<!std::is_void<Type>::value, Array, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& core(const Index index, const Subscript subscript, const Type& value); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& core(const Subscript subscript, Dummy...); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Array, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& core(const Subscript subscript, Dummy...) const; 
        template <unsigned int Index, typename Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Array, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& core(const Subscript subscript, const Type& value); 
    //@}

    // Extra 
    /// \name           Extra 
    //@{
    public: 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type> inline Type& extra(Dummy...); 
        template <class... Dummy, class Type = typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type> inline const Type& extra(Dummy...) const; 
        template <class Type, class = typename std::enable_if<std::is_convertible<Type, typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::value>::type> inline Crtp& extra(const Type& value); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::value>::type> inline Type& extra(Dummy...); 
        template <unsigned int Index, class... Dummy, class Type = typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::value>::type> inline const Type& extra(Dummy...) const; 
        template <unsigned int Index, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::type>::value>::type> inline Crtp& extra(const Type& value); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::type>::value>::type> inline Type& extra(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>::type>::value>::type> inline const Type& extra(Dummy...) const; 
        template <unsigned int Index, unsigned int Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::value>::type, class = typename std::enable_if<Subscript+1 <= std::tuple_size<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::type>::type>::value>::type> inline Crtp& extra(const Type& value); 
        template <typename Index, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type>::value>::type> inline Type& extra(const Index index, Dummy...); 
        template <typename Index, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type>::value>::type> inline const Type& extra(const Index index, Dummy...) const; 
        template <typename Index, class Type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>()[Index()])>::type>::type>::value>::type> inline Crtp& extra(const Index index, const Type& value); 
        template <typename Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& extra(const Index index, const Subscript subscript, Dummy...); 
        template <typename Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& extra(const Index index, const Subscript subscript, Dummy...) const; 
        template <typename Index, typename Subscript, class Type, class = typename std::enable_if<std::is_convertible<Index, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>()[Index()])>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& extra(const Index index, const Subscript subscript, const Type& value); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Type& extra(const Subscript subscript, Dummy...); 
        template <unsigned int Index, typename Subscript, class... Dummy, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>()[Subscript()])>::type>::type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<sizeof...(Dummy) == 0>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<sizeof...(Dummy) == 0, Tuple, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline const Type& extra(const Subscript subscript, Dummy...) const; 
        template <unsigned int Index, typename Subscript, class Type, class = typename std::enable_if<Index+1 <= std::tuple_size<typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::value>::type, class = typename std::enable_if<std::is_convertible<Subscript, unsigned int>::value>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, typename std::conditional<!std::is_void<Type>::value, Tuple, void>::type>::type>()[Subscript()])>::type>::type>::value>::type> inline Crtp& extra(const Subscript subscript, const Type& value); 
    //@}

    // Stream 
    /// \name           Stream 
    //@{
    public: 
        template <class SelfCrtp, class SelfScalar, class SelfArray, class SelfTuple> friend std::ostream& operator<<(std::ostream& lhs, const AbstractStep<SelfCrtp, SelfScalar, SelfArray, SelfTuple>& rhs); 
    //@}

    // Types 
    /// \name           Types 
    //@{
    public: 
        template <class Type = std::tuple<Scalar, Array, Tuple>, class = typename std::enable_if<std::is_convertible<Type, std::tuple<Scalar, Array, Tuple> >::value>::type> static constexpr Type type(); 
        template <unsigned int Index, class Type = typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::value)>::type> static constexpr Type type(Dummy...); 
        template <unsigned int Index, unsigned int Subscript, class Type = typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Type, typename std::tuple_element<Subscript, typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>::type>::value)>::type> static constexpr Type type(const Dummy&...); 
        template <unsigned int Index, class Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>()[0])>::type>::type, class = typename std::enable_if<std::is_convertible<Type, typename std::remove_cv<typename std::remove_reference<decltype(std::declval<typename std::tuple_element<Index, std::tuple<Scalar, Array, Tuple> >::type>()[0])>::type>::type>::value>::type> static constexpr Type type(const unsigned int subscript); 
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
        template <int Exponent = 1, class Coefficient = std::ratio<1>, class Type, class = typename std::enable_if<Coefficient::num+Coefficient::den != Coefficient::num>::type> static constexpr Type monomial(Type&& value); 
        template <class Coefficient = std::ratio<1>, int Exponent = 1, class Type, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (Coefficient::num+Coefficient::den != Coefficient::num)>::type> static constexpr Type monomial(Type&& value, Dummy...); 
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
        std::tuple<Scalar, Array, Tuple> _data;                                 ///< Internal tuple container. 
    //@}
};
//--------------------------------------------------------------------------- //



//--------------------------- PROTECTED LIFECYCLE --------------------------- //
// Protected destructor 
/// \brief          Protected destructor. 
/// \details        Avoids direct instantiation of the class, and only allows 
///                 it through its derived children. 
template <class Crtp, class Scalar, class Array, class Tuple>
inline AbstractStep<Crtp, Scalar, Array, Tuple>::~AbstractStep()
= default;
//--------------------------------------------------------------------------- //



//-------------------------------- LIFECYCLE -------------------------------- //
// Implicit empty constructor 
/// \brief          Implicit empty constructor. 
/// \details        Provides an implicit construction of an object initialized 
///                 to its default value. 
template <class Crtp, class Scalar, class Array, class Tuple>
inline AbstractStep<Crtp, Scalar, Array, Tuple>::AbstractStep()
: _data(std::tuple<Scalar, Array, Tuple>())
{
    ;
}

// Explicit conversion constructor 
/// \brief          Explicit conversion constructor. 
/// \details        Provides an explicit construction from another type of 
///                 object. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      source Source of the copy.
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple>
inline AbstractStep<Crtp, Scalar, Array, Tuple>::AbstractStep(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& source)
: _data(transmute<std::tuple<Scalar, Array, Tuple> >(source.data()))
{
    ;
}

// Explicit step constructor 
/// \brief          Explicit step constructor. 
/// \details        Provides an explicit construction from step components. 
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      scalar Source of the scalar.
/// \param[in]      array Source of the array.
/// \param[in]      tuple Source of the tuple.
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherScalar, class OtherArray, class OtherTuple, class>
inline AbstractStep<Crtp, Scalar, Array, Tuple>::AbstractStep(const OtherScalar& scalar, const OtherArray& array, const OtherTuple& tuple)
: _data(transmute<std::tuple<Scalar, Array, Tuple> >(scalar, array, tuple))
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
template <class Crtp, class Scalar, class Array, class Tuple>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::operator=(const AbstractStep<Crtp, Scalar, Array, Tuple>& rhs)
{
    _data = rhs._data;
    return static_cast<Crtp&>(*this);
}

// Conversion assignment operator 
/// \brief          Conversion assignment operator. 
/// \details        Assigns contents from another type of object. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::operator=(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& rhs)
{
    _data = transmute<std::tuple<Scalar, Array, Tuple> >(rhs.data());
    return static_cast<Crtp&>(*this);
}

// Equal to 
/// \brief          Equal to. 
/// \details        Compares for equality and returns true if the contents is 
///                 equal. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      rhs Right-hand side.
/// \return         True if equal, false if not equal. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple>
inline bool AbstractStep<Crtp, Scalar, Array, Tuple>::operator==(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& rhs) const
{
    return (_data == rhs.data());
}

// Not equal to 
/// \brief          Not equal to. 
/// \details        Compares for difference and returns true if the contents 
///                 is different. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      rhs Right-hand side.
/// \return         True if not equal, false if equal. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple>
inline bool AbstractStep<Crtp, Scalar, Array, Tuple>::operator!=(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& rhs) const
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
template <class Crtp, class Scalar, class Array, class Tuple>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::assign()
{
    _data = std::tuple<Scalar, Array, Tuple>();
    return static_cast<Crtp&>(*this);;
}

// Copy assignment 
/// \brief          Copy assignment. 
/// \details        Assigns contents from the same type of object. 
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::assign(const AbstractStep<Crtp, Scalar, Array, Tuple>& source)
{
    _data = source._data;
    return static_cast<Crtp&>(*this);;
}

// Conversion assignment 
/// \brief          Conversion assignment. 
/// \details        Assigns contents from another type of object. 
/// \tparam         OtherCrtp (Other derived CRTP class.)
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherCrtp, class OtherScalar, class OtherArray, class OtherTuple>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::assign(const AbstractStep<OtherCrtp, OtherScalar, OtherArray, OtherTuple>& source)
{
    _data = transmute<std::tuple<Scalar, Array, Tuple> >(source.data());
    return static_cast<Crtp&>(*this);
}

// Step assignment 
/// \brief          Step assignment. 
/// \details        Assigns contents from step components. 
/// \tparam         OtherScalar (Other scalar type of id.)
/// \tparam         OtherArray (Other array type of core quantities.)
/// \tparam         OtherTuple (Other tuple type of extra quantities.)
/// \param[in]      scalar Source of the scalar.
/// \param[in]      array Source of the array.
/// \param[in]      tuple Source of the tuple.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherScalar, class OtherArray, class OtherTuple, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::assign(const OtherScalar& scalar, const OtherArray& array, const OtherTuple& tuple)
{
    _data = transmute<std::tuple<Scalar, Array, Tuple> >(scalar, array, tuple);
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//-------------------------------- MANAGEMENT ------------------------------- //
// Nullify 
/// \brief          Nullify. 
/// \details        Resets all data members to their default values. 
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::nullify()
{
    _data = std::tuple<Scalar, Array, Tuple>();
    return static_cast<Crtp&>(*this);
}

// Copy 
/// \brief          Copy. 
/// \details        Generates a copy of the object. 
/// \return         Copy. 
template <class Crtp, class Scalar, class Array, class Tuple>
inline Crtp AbstractStep<Crtp, Scalar, Array, Tuple>::copy() const
{
    return static_cast<const Crtp&>(*this);
}

// Cast 
/// \brief          Cast. 
/// \details        Casts contents to another object type. 
/// \tparam         OtherCrtp Other derived CRTP class.
/// \return         Casted copy. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class OtherCrtp, class>
inline OtherCrtp AbstractStep<Crtp, Scalar, Array, Tuple>::cast() const
{
    return OtherCrtp(*this);
}
//--------------------------------------------------------------------------- //



//----------------------------------- DATA ---------------------------------- //
// Unified data access 
/// \brief          Unified data access. 
/// \details        Provides a direct access to the data. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Data std::tuple<Scalar, Array, Tuple> type.)
/// \return         Reference to the data. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(Dummy...)
{
    return _data;
}

// Unified data getter 
/// \brief          Unified data getter. 
/// \details        Gets the data. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Data std::tuple<Scalar, Array, Tuple> type.)
/// \return         Immutable reference to the data. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(Dummy...) const
{
    return _data;
}

// Unified data setter 
/// \brief          Unified data setter. 
/// \details        Sets the data. 
/// \param[in]      value Data value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::data(const Type& value)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(Dummy...) const
{
    return std::get<Index>(_data);
}

// Unified data component setter 
/// \brief          Unified data component setter. 
/// \details        Sets the specified component of the data. 
/// \tparam         Index Index of the component.
/// \param[in]      value Component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::data(const Type& value)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(Dummy...) const
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::data(const Type& value)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(const Subscript subscript, Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::data(const Subscript subscript, Dummy...) const
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::data(const Subscript subscript, const Type& value)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::get() const
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class... Dummy, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::get(Dummy...) const
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::get(const Dummy&...) const
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::get(const unsigned int subscript) const
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::set(const Type& value)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class... Dummy, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::set(const Type& value, Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::set(const Type& value, const Dummy&...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::set(const unsigned int subscript, const Type& value)
{
    std::get<Index>(_data)[subscript] = value;
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//------------------------------------ ID ----------------------------------- //
// Id access 
/// \brief          Id access. 
/// \details        Provides a direct access to the id scalar. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Id scalar type.)
/// \return         Reference to the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(Dummy...)
{
    return std::get<0>(_data);
}

// Id getter 
/// \brief          Id getter. 
/// \details        Gets the id scalar. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Id scalar type.)
/// \return         Immutable reference to the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(Dummy...) const
{
    return std::get<0>(_data);
}

// Id setter 
/// \brief          Id setter. 
/// \details        Sets the id scalar. 
/// \param[in]      value Id scalar value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Type& value)
{
    std::get<0>(_data) = value;
    return static_cast<Crtp&>(*this);
}

// Id component access 
/// \brief          Id component access. 
/// \details        Provides a direct access to the specified component of the 
///                 id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Reference to the component of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(Dummy...)
{
    return std::get<Index>(std::get<0>(_data));
}

// Id component getter 
/// \brief          Id component getter. 
/// \details        Gets the specified component of the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Immutable reference to the component of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(Dummy...) const
{
    return std::get<Index>(std::get<0>(_data));
}

// Id component setter 
/// \brief          Id component setter. 
/// \details        Sets the specified component of the id scalar. 
/// \tparam         Index Index of the component.
/// \param[in]      value Component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Type& value)
{
    std::get<Index>(std::get<0>(_data)) = value;
    return static_cast<Crtp&>(*this);
}

// Id inner component access 
/// \brief          Id inner component access. 
/// \details        Provides a direct access to the specified inner component 
///                 of the specified component of the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Reference to the inner component of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(Dummy...)
{
    return std::get<Subscript>(std::get<Index>(std::get<0>(_data)));
}

// Id inner component getter 
/// \brief          Id inner component getter. 
/// \details        Gets the specified inner component of the specified 
///                 component of the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Immutable reference to the inner component of the id 
///                 scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(Dummy...) const
{
    return std::get<Subscript>(std::get<Index>(std::get<0>(_data)));
}

// Id inner component setter 
/// \brief          Id inner component setter. 
/// \details        Sets the specified inner component of the specified 
///                 component of the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \param[in]      value Inner component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Type& value)
{
    std::get<Subscript>(std::get<Index>(std::get<0>(_data))) = value;
    return static_cast<Crtp&>(*this);
}

// Id element access 
/// \brief          Id element access. 
/// \details        Provides a direct access to the specified element of the 
///                 id scalar. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Element type.)
/// \param[in]      index Index of the component.
/// \return         Reference to the element of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Index index, Dummy...)
{
    return std::get<0>(_data)[index];
}

// Id element getter 
/// \brief          Id element getter. 
/// \details        Gets the specified element of the id scalar. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Element type.)
/// \param[in]      index Index of the component.
/// \return         Immutable reference to the element of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Index index, Dummy...) const
{
    return std::get<0>(_data)[index];
}

// Id element setter 
/// \brief          Id element setter. 
/// \details        Sets the specified element of the id scalar. 
/// \tparam         Index (Index type of the component.)
/// \param[in]      index Index of the component.
/// \param[in]      value Element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Index index, const Type& value)
{
    std::get<0>(_data)[index] = value;
    return static_cast<Crtp&>(*this);
}

// Id inner element access 
/// \brief          Id inner element access. 
/// \details        Provides a direct access to the specified inner element of 
///                 the specified element of the id scalar. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner element type.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \return         Reference to the inner element of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Index index, const Subscript subscript, Dummy...)
{
    return std::get<0>(_data)[index][subscript];
}

// Id inner element getter 
/// \brief          Id inner element getter. 
/// \details        Gets the specified inner element of the specified element 
///                 of the id scalar. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner element type.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \return         Immutable reference to the inner element of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Index index, const Subscript subscript, Dummy...) const
{
    return std::get<0>(_data)[index][subscript];
}

// Id inner element setter 
/// \brief          Id inner element setter. 
/// \details        Sets the specified inner element of the specified element 
///                 of the id scalar. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \param[in]      value Inner element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Index index, const Subscript subscript, const Type& value)
{
    std::get<0>(_data)[index][subscript] = value;
    return static_cast<Crtp&>(*this);
}

// Id component element access 
/// \brief          Id component element access. 
/// \details        Provides a direct access to the specified element of the 
///                 specified component of the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Reference to the component element of the id scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Subscript subscript, Dummy...)
{
    return std::get<Index>(std::get<0>(_data))[subscript];
}

// Id component element getter 
/// \brief          Id component element getter. 
/// \details        Gets the specified element of the specified component of 
///                 the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Immutable reference to the component element of the id 
///                 scalar. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Subscript subscript, Dummy...) const
{
    return std::get<Index>(std::get<0>(_data))[subscript];
}

// Id component element setter 
/// \brief          Id component element setter. 
/// \details        Sets the specified element of the specified component of 
///                 the id scalar. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \param[in]      subscript Subscript of the component element.
/// \param[in]      value Component element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::id(const Subscript subscript, const Type& value)
{
    std::get<Index>(std::get<0>(_data))[subscript] = value;
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//----------------------------------- CORE ---------------------------------- //
// Core access 
/// \brief          Core access. 
/// \details        Provides a direct access to the core array. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Core array type.)
/// \return         Reference to the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(Dummy...)
{
    return std::get<1>(_data);
}

// Core getter 
/// \brief          Core getter. 
/// \details        Gets the core array. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Core array type.)
/// \return         Immutable reference to the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(Dummy...) const
{
    return std::get<1>(_data);
}

// Core setter 
/// \brief          Core setter. 
/// \details        Sets the core array. 
/// \param[in]      value Core array value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Type& value)
{
    std::get<1>(_data) = value;
    return static_cast<Crtp&>(*this);
}

// Core component access 
/// \brief          Core component access. 
/// \details        Provides a direct access to the specified component of the 
///                 core array. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Reference to the component of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(Dummy...)
{
    return std::get<Index>(std::get<1>(_data));
}

// Core component getter 
/// \brief          Core component getter. 
/// \details        Gets the specified component of the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Immutable reference to the component of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(Dummy...) const
{
    return std::get<Index>(std::get<1>(_data));
}

// Core component setter 
/// \brief          Core component setter. 
/// \details        Sets the specified component of the core array. 
/// \tparam         Index Index of the component.
/// \param[in]      value Component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Type& value)
{
    std::get<Index>(std::get<1>(_data)) = value;
    return static_cast<Crtp&>(*this);
}

// Core inner component access 
/// \brief          Core inner component access. 
/// \details        Provides a direct access to the specified inner component 
///                 of the specified component of the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Reference to the inner component of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(Dummy...)
{
    return std::get<Subscript>(std::get<Index>(std::get<1>(_data)));
}

// Core inner component getter 
/// \brief          Core inner component getter. 
/// \details        Gets the specified inner component of the specified 
///                 component of the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Immutable reference to the inner component of the core 
///                 array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(Dummy...) const
{
    return std::get<Subscript>(std::get<Index>(std::get<1>(_data)));
}

// Core inner component setter 
/// \brief          Core inner component setter. 
/// \details        Sets the specified inner component of the specified 
///                 component of the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \param[in]      value Inner component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Type& value)
{
    std::get<Subscript>(std::get<Index>(std::get<1>(_data))) = value;
    return static_cast<Crtp&>(*this);
}

// Core element access 
/// \brief          Core element access. 
/// \details        Provides a direct access to the specified element of the 
///                 core array. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Element type.)
/// \param[in]      index Index of the component.
/// \return         Reference to the element of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Index index, Dummy...)
{
    return std::get<1>(_data)[index];
}

// Core element getter 
/// \brief          Core element getter. 
/// \details        Gets the specified element of the core array. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Element type.)
/// \param[in]      index Index of the component.
/// \return         Immutable reference to the element of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Index index, Dummy...) const
{
    return std::get<1>(_data)[index];
}

// Core element setter 
/// \brief          Core element setter. 
/// \details        Sets the specified element of the core array. 
/// \tparam         Index (Index type of the component.)
/// \param[in]      index Index of the component.
/// \param[in]      value Element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Index index, const Type& value)
{
    std::get<1>(_data)[index] = value;
    return static_cast<Crtp&>(*this);
}

// Core inner element access 
/// \brief          Core inner element access. 
/// \details        Provides a direct access to the specified inner element of 
///                 the specified element of the core array. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner element type.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \return         Reference to the inner element of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Index index, const Subscript subscript, Dummy...)
{
    return std::get<1>(_data)[index][subscript];
}

// Core inner element getter 
/// \brief          Core inner element getter. 
/// \details        Gets the specified inner element of the specified element 
///                 of the core array. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner element type.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \return         Immutable reference to the inner element of the core 
///                 array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Index index, const Subscript subscript, Dummy...) const
{
    return std::get<1>(_data)[index][subscript];
}

// Core inner element setter 
/// \brief          Core inner element setter. 
/// \details        Sets the specified inner element of the specified element 
///                 of the core array. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \param[in]      value Inner element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Index index, const Subscript subscript, const Type& value)
{
    std::get<1>(_data)[index][subscript] = value;
    return static_cast<Crtp&>(*this);
}

// Core component element access 
/// \brief          Core component element access. 
/// \details        Provides a direct access to the specified element of the 
///                 specified component of the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Reference to the component element of the core array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Subscript subscript, Dummy...)
{
    return std::get<Index>(std::get<1>(_data))[subscript];
}

// Core component element getter 
/// \brief          Core component element getter. 
/// \details        Gets the specified element of the specified component of 
///                 the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Immutable reference to the component element of the core 
///                 array. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Subscript subscript, Dummy...) const
{
    return std::get<Index>(std::get<1>(_data))[subscript];
}

// Core component element setter 
/// \brief          Core component element setter. 
/// \details        Sets the specified element of the specified component of 
///                 the core array. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \param[in]      subscript Subscript of the component element.
/// \param[in]      value Component element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::core(const Subscript subscript, const Type& value)
{
    std::get<Index>(std::get<1>(_data))[subscript] = value;
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//---------------------------------- EXTRA ---------------------------------- //
// Extra access 
/// \brief          Extra access. 
/// \details        Provides a direct access to the extra tuple. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Extra tuple type.)
/// \return         Reference to the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(Dummy...)
{
    return std::get<2>(_data);
}

// Extra getter 
/// \brief          Extra getter. 
/// \details        Gets the extra tuple. 
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Extra tuple type.)
/// \return         Immutable reference to the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class... Dummy, class Type, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(Dummy...) const
{
    return std::get<2>(_data);
}

// Extra setter 
/// \brief          Extra setter. 
/// \details        Sets the extra tuple. 
/// \param[in]      value Extra tuple value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Type& value)
{
    std::get<2>(_data) = value;
    return static_cast<Crtp&>(*this);
}

// Extra component access 
/// \brief          Extra component access. 
/// \details        Provides a direct access to the specified component of the 
///                 extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Reference to the component of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(Dummy...)
{
    return std::get<Index>(std::get<2>(_data));
}

// Extra component getter 
/// \brief          Extra component getter. 
/// \details        Gets the specified component of the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component type.)
/// \return         Immutable reference to the component of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(Dummy...) const
{
    return std::get<Index>(std::get<2>(_data));
}

// Extra component setter 
/// \brief          Extra component setter. 
/// \details        Sets the specified component of the extra tuple. 
/// \tparam         Index Index of the component.
/// \param[in]      value Component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Type& value)
{
    std::get<Index>(std::get<2>(_data)) = value;
    return static_cast<Crtp&>(*this);
}

// Extra inner component access 
/// \brief          Extra inner component access. 
/// \details        Provides a direct access to the specified inner component 
///                 of the specified component of the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Reference to the inner component of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(Dummy...)
{
    return std::get<Subscript>(std::get<Index>(std::get<2>(_data)));
}

// Extra inner component getter 
/// \brief          Extra inner component getter. 
/// \details        Gets the specified inner component of the specified 
///                 component of the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner component type.)
/// \return         Immutable reference to the inner component of the extra 
///                 tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(Dummy...) const
{
    return std::get<Subscript>(std::get<Index>(std::get<2>(_data)));
}

// Extra inner component setter 
/// \brief          Extra inner component setter. 
/// \details        Sets the specified inner component of the specified 
///                 component of the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript Subscript of the inner component.
/// \param[in]      value Inner component value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Type& value)
{
    std::get<Subscript>(std::get<Index>(std::get<2>(_data))) = value;
    return static_cast<Crtp&>(*this);
}

// Extra element access 
/// \brief          Extra element access. 
/// \details        Provides a direct access to the specified element of the 
///                 extra tuple. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Element type.)
/// \param[in]      index Index of the component.
/// \return         Reference to the element of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class... Dummy, class Type, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Index index, Dummy...)
{
    return std::get<2>(_data)[index];
}

// Extra element getter 
/// \brief          Extra element getter. 
/// \details        Gets the specified element of the extra tuple. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Element type.)
/// \param[in]      index Index of the component.
/// \return         Immutable reference to the element of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class... Dummy, class Type, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Index index, Dummy...) const
{
    return std::get<2>(_data)[index];
}

// Extra element setter 
/// \brief          Extra element setter. 
/// \details        Sets the specified element of the extra tuple. 
/// \tparam         Index (Index type of the component.)
/// \param[in]      index Index of the component.
/// \param[in]      value Element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, class Type, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Index index, const Type& value)
{
    std::get<2>(_data)[index] = value;
    return static_cast<Crtp&>(*this);
}

// Extra inner element access 
/// \brief          Extra inner element access. 
/// \details        Provides a direct access to the specified inner element of 
///                 the specified element of the extra tuple. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner element type.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \return         Reference to the inner element of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Index index, const Subscript subscript, Dummy...)
{
    return std::get<2>(_data)[index][subscript];
}

// Extra inner element getter 
/// \brief          Extra inner element getter. 
/// \details        Gets the specified inner element of the specified element 
///                 of the extra tuple. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Inner element type.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \return         Immutable reference to the inner element of the extra 
///                 tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Index index, const Subscript subscript, Dummy...) const
{
    return std::get<2>(_data)[index][subscript];
}

// Extra inner element setter 
/// \brief          Extra inner element setter. 
/// \details        Sets the specified inner element of the specified element 
///                 of the extra tuple. 
/// \tparam         Index (Index type of the component.)
/// \tparam         Subscript (Subscript type of the inner element.)
/// \param[in]      index Index of the component.
/// \param[in]      subscript Subscript of the inner element.
/// \param[in]      value Inner element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <typename Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Index index, const Subscript subscript, const Type& value)
{
    std::get<2>(_data)[index][subscript] = value;
    return static_cast<Crtp&>(*this);
}

// Extra component element access 
/// \brief          Extra component element access. 
/// \details        Provides a direct access to the specified element of the 
///                 specified component of the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Reference to the component element of the extra tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Subscript subscript, Dummy...)
{
    return std::get<Index>(std::get<2>(_data))[subscript];
}

// Extra component element getter 
/// \brief          Extra component element getter. 
/// \details        Gets the specified element of the specified component of 
///                 the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \tparam         Dummy (Dummy types.)
/// \tparam         Type (Component element type.)
/// \param[in]      subscript Subscript of the component element.
/// \return         Immutable reference to the component element of the extra 
///                 tuple. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class... Dummy, class Type, class, class, class, class>
inline const Type& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Subscript subscript, Dummy...) const
{
    return std::get<Index>(std::get<2>(_data))[subscript];
}

// Extra component element setter 
/// \brief          Extra component element setter. 
/// \details        Sets the specified element of the specified component of 
///                 the extra tuple. 
/// \tparam         Index Index of the component.
/// \tparam         Subscript (Subscript type of the component element.)
/// \param[in]      subscript Subscript of the component element.
/// \param[in]      value Component element value.
/// \return         Self reference. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, typename Subscript, class Type, class, class, class>
inline Crtp& AbstractStep<Crtp, Scalar, Array, Tuple>::extra(const Subscript subscript, const Type& value)
{
    std::get<Index>(std::get<2>(_data))[subscript] = value;
    return static_cast<Crtp&>(*this);
}
//--------------------------------------------------------------------------- //



//---------------------------------- STREAM --------------------------------- //
// Output stream operator 
/// \brief          Output stream operator. 
/// \details        Adds each element to the stream. 
/// \tparam         SelfCrtp (Derived CRTP class.)
/// \tparam         SelfScalar (Scalar type of id.)
/// \tparam         SelfArray (Array type of core quantities.)
/// \tparam         SelfTuple (Tuple type of extra quantities.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side object.
/// \return         Output stream. 
template <class SelfCrtp, class SelfScalar, class SelfArray, class SelfTuple>
std::ostream& operator<<(std::ostream& lhs, const AbstractStep<SelfCrtp, SelfScalar, SelfArray, SelfTuple>& rhs)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Type, class>
constexpr Type AbstractStep<Crtp, Scalar, Array, Tuple>::type()
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class... Dummy, class>
constexpr Type AbstractStep<Crtp, Scalar, Array, Tuple>::type(Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, unsigned int Subscript, class Type, class... Dummy, class>
constexpr Type AbstractStep<Crtp, Scalar, Array, Tuple>::type(const Dummy&...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Index, class Type, class>
constexpr Type AbstractStep<Crtp, Scalar, Array, Tuple>::type(const unsigned int subscript)
{
    return (subscript != 0) ? (Type()) : (Type());
}
//--------------------------------------------------------------------------- //



//-------------------------------- PROPERTIES ------------------------------- //
// Total number of types 
/// \brief          Total number of types. 
/// \details        Counts the number of internal components. 
/// \return         Number of components. 
template <class Crtp, class Scalar, class Array, class Tuple>
constexpr unsigned int AbstractStep<Crtp, Scalar, Array, Tuple>::types()
{
    return std::tuple_size<std::tuple<Scalar, Array, Tuple> >::value;
}
//--------------------------------------------------------------------------- //



//--------------------------------- HELPERS --------------------------------- //
// Exponent coefficient monomial 
/// \brief          Exponent coefficient monomial. 
/// \details        Computes a monomial of a value with exponent priority. 
/// \tparam         Exponent Value of the exponent \f$n\f$.
/// \tparam         Coefficient Value of the coefficient \f$c\f$.
/// \tparam         Type (Type of the value \f$x\f$.)
/// \param[in]      value Value \f$x\f$.
/// \return         Value of \f$ \left(x^{n}\right) \times c\f$. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <int Exponent, class Coefficient, class Type, class>
constexpr Type AbstractStep<Crtp, Scalar, Array, Tuple>::monomial(Type&& value)
{
    return (Exponent > 0) ? (std::forward<Type>(value)*monomial<(Exponent-1)*(Exponent > 1), Coefficient>(std::forward<Type>(value))) : ((Exponent < 0) ? (Type(1)/monomial<(-Exponent)*(Exponent < 0), Coefficient>(std::forward<Type>(value))) : (Type(Coefficient::num)/Type(Coefficient::den)));
}

// Coefficient exponent monomial 
/// \brief          Coefficient exponent monomial. 
/// \details        Computes a monomial of a value with coefficient priority. 
/// \tparam         Coefficient Value of the coefficient \f$c\f$.
/// \tparam         Exponent Value of the exponent \f$n\f$.
/// \tparam         Type (Type of the value \f$x\f$.)
/// \tparam         Dummy (Dummy types.)
/// \param[in]      value Value \f$x\f$.
/// \return         Value of \f$\left(c \times x\right)^{n}\f$. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Coefficient, int Exponent, class Type, class... Dummy, class>
constexpr Type AbstractStep<Crtp, Scalar, Array, Tuple>::monomial(Type&& value, Dummy...)
{
    return (Exponent > 1) ? (std::forward<Type>(value)*(Type(Coefficient::num)/Type(Coefficient::den))*monomial<(Exponent-1)*(Exponent > 1), Coefficient>(std::forward<Type>(value))) : ((Exponent < 0) ? (Type(1)/monomial<(-Exponent)*(Exponent < 0), Coefficient>(std::forward<Type>(value))) : ((Exponent == 1) ? (std::forward<Type>(value)*(Type(Coefficient::num)/Type(Coefficient::den))) : (Type(1))));
}

// Default transmutation 
/// \brief          Default transmutation. 
/// \details        Transmutes the input into the specified output type by 
///                 creating a new default value. 
/// \tparam         Output Output type.
/// \return         Default output value. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Output>
constexpr Output AbstractStep<Crtp, Scalar, Array, Tuple>::transmute()
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Output, class Input, class>
constexpr Input AbstractStep<Crtp, Scalar, Array, Tuple>::transmute(Input&& input)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Output, class... Input, class>
constexpr Output AbstractStep<Crtp, Scalar, Array, Tuple>::transmute(const Input&... input)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Output, class... Input, class>
constexpr Output AbstractStep<Crtp, Scalar, Array, Tuple>::transmute(Input&&... input)
{
    return transmute<Output>(std::forward_as_tuple(std::forward<Input>(input)...));
}

// Is printable 
/// \brief          Is printable. 
/// \details        True overload to check at compile time whether a type can 
///                 be inserted in an ostream. 
/// \tparam         Input Input type.
/// \return         True. 
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Input>
constexpr bool AbstractStep<Crtp, Scalar, Array, Tuple>::printable(typename std::enable_if<!std::is_void<decltype(std::declval<std::ostream&>()<<std::declval<Input>())>::value>::type*)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Input, class... Dummy, class>
constexpr bool AbstractStep<Crtp, Scalar, Array, Tuple>::printable(Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Input, class>
inline bool AbstractStep<Crtp, Scalar, Array, Tuple>::print(std::ostream& stream, Input&& input)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <class Input, class>
inline bool AbstractStep<Crtp, Scalar, Array, Tuple>::print(std::ostream& stream, const Input& input)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Current, class Input, class... Dummy, class>
inline bool AbstractStep<Crtp, Scalar, Array, Tuple>::print(std::ostream& stream, const Input& input, Dummy...)
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
template <class Crtp, class Scalar, class Array, class Tuple>
template <unsigned int Current, class... Dummy, class>
inline bool AbstractStep<Crtp, Scalar, Array, Tuple>::print(std::ostream& stream, const Dummy&...)
{
    return stream.good();
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of AbstractStep. 
/// \return         0 if no error. 
template <class Crtp, class Scalar, class Array, class Tuple>
int AbstractStep<Crtp, Scalar, Array, Tuple>::example()
{
    std::cout<<"BEGIN = AbstractStep::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = AbstractStep::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = AbstractStep::example()"<<std::endl;
    return 1;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABSTRACTSTEP_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

