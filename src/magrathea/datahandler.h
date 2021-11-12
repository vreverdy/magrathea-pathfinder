/* ******************************* DATAHANDLER ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA-PATHFINDER
// TITLE :          DataHandler
// DESCRIPTION :    Set of basic operations on binary data related to IO
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           datahandler.h
/// \brief          Set of basic operations on binary data related to IO
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef DATAHANDLER_H_INCLUDED
#define DATAHANDLER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <type_traits>
#include <algorithm>
#include <limits>
#include <iterator>
#include <array>
#include <tuple>
#include <vector>
#include <bitset>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Set of basic operations on binary data related to IO
/// \brief          Set of basic operations on binary data related to IO.
/// \details        Provides a wide range of overloaded utilities to handle and 
///                 format data for reading and writing tasks : byteswap, size, 
///                 extraction from tuples...
class DataHandler final
{
    // Utilities
    /// \name           Utilities
    //@{
    public:
        template <typename Type> static constexpr bool array(const Type&); 
        template <typename Type, std::size_t Size> static constexpr bool array(const std::array<Type, Size>&); 
        template <typename Type> static constexpr bool tuple(const Type&); 
        template <typename... Types> static constexpr bool tuple(const std::tuple<Types...>&); 
    //@}

    // Size
    /// \name           Size
    //@{
    public:
        template <typename... Types> static constexpr unsigned int size();
        static constexpr unsigned int size(const std::tuple<>&);
        template <typename Type> static constexpr unsigned int size(const Type& variable);
        template <typename Type, typename... Types> static constexpr unsigned int size(const Type& variable, const Types&... variables);
        template <typename Type, std::size_t Size> static constexpr unsigned int size(const std::array<Type, Size>& container);
        template <unsigned int Current = 0, typename... Types> static constexpr unsigned int size(const std::tuple<Types...>& container);
        static inline unsigned long long int rsize(const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <typename Type> static inline unsigned long long int rsize(const Type* const first, const Type* const last);
        template <typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline unsigned long long int rsize(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Bytesize
    /// \name           Bytesize
    //@{
    public:
        template <typename Type = void> static constexpr unsigned int bytesize(const std::nullptr_t = nullptr);
        template <typename Type, typename... Types, class = typename std::enable_if<sizeof...(Types) != 0>::type> static constexpr unsigned int bytesize();
        static constexpr unsigned int bytesize(const std::tuple<>&);
        template <typename Type> static constexpr unsigned int bytesize(const Type& variable);
        template <typename Type, typename... Types> static constexpr unsigned int bytesize(const Type& variable, const Types&... variables);
        template <typename Type, std::size_t Size> static constexpr unsigned int bytesize(const std::array<Type, Size>& container);
        template <unsigned int Current = 0, typename... Types> static constexpr unsigned int bytesize(const std::tuple<Types...>& container);
        static inline unsigned long long int rbytesize(const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <typename Type> static inline unsigned long long int rbytesize(const Type* const first, const Type* const last);
        template <typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline unsigned long long int rbytesize(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}
    
    // Byteswap
    /// \name           Byteswap
    //@{
    public:
        template <bool Byteswap, typename... Types, class = typename std::enable_if<!Byteswap>::type> static inline unsigned int byteswap(const Types&...);
        template <bool Byteswap, typename... Types, class = typename std::enable_if<!Byteswap>::type> static inline unsigned int rbyteswap(const Types&...);
        template <bool Byteswap = true, class = typename std::enable_if<Byteswap>::type> static inline unsigned int byteswap();
        template <bool Byteswap = true, class = typename std::enable_if<Byteswap>::type> static inline unsigned int byteswap(std::tuple<>&);
        template <bool Byteswap = true, typename Type, class = typename std::enable_if<Byteswap>::type> static inline unsigned int byteswap(Type& variable);
        template <bool Byteswap = true, typename Type, typename... Types, class = typename std::enable_if<Byteswap>::type> static inline unsigned int byteswap(Type& variable, Types&... variables);
        template <bool Byteswap = true, typename Type, std::size_t Size, class = typename std::enable_if<Byteswap>::type> static inline unsigned int byteswap(std::array<Type, Size>& container);
        template <bool Byteswap = true, unsigned int Current = 0, typename... Types, class = typename std::enable_if<Byteswap>::type> static inline unsigned int byteswap(std::tuple<Types...>& container);
        template <bool Byteswap = true, class = typename std::enable_if<Byteswap>::type> static inline unsigned long long int rbyteswap(const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <bool Byteswap = true, typename Type, class = typename std::enable_if<Byteswap>::type> static inline unsigned long long int rbyteswap(Type* const first, Type* const last);
        template <bool Byteswap = true, typename Type, class = typename std::enable_if<(Byteswap) && (!std::is_pointer<Type>::value)>::type> static inline unsigned long long int rbyteswap(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}
    
    // Write stream
    /// \name           Write stream
    //@{
    public:
        template <bool Byteswap = false> static inline bool write(std::ostream& stream);
        template <bool Byteswap = false> static inline bool write(std::ostream& stream, const std::tuple<>&);
        template <bool Byteswap = false, typename Type> static inline bool write(std::ostream& stream, const Type& variable);
        template <bool Byteswap = false, typename Type, typename... Types> static inline bool write(std::ostream& stream, const Type& variable, const Types&... variables);
        template <bool Byteswap = false, typename Type, std::size_t Size> static inline bool write(std::ostream& stream, const std::array<Type, Size>& container);
        template <bool Byteswap = false, unsigned int Current = 0, typename... Types> static inline bool write(std::ostream& stream, const std::tuple<Types...>& container);
        template <bool Byteswap = false> static inline bool rwrite(std::ostream& stream, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <bool Byteswap = false, typename Type> static inline bool rwrite(std::ostream& stream, const Type* const first, const Type* const last);
        template <bool Byteswap = false, typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline bool rwrite(std::ostream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Read stream
    /// \name           Read stream
    //@{
    public:
        template <bool Byteswap = false> static inline bool read(std::istream& stream);
        template <bool Byteswap = false> static inline bool read(std::istream& stream, std::tuple<>&);
        template <bool Byteswap = false, typename Type> static inline bool read(std::istream& stream, Type& variable);
        template <bool Byteswap = false, typename Type, typename... Types> static inline bool read(std::istream& stream, Type& variable, Types&... variables);
        template <bool Byteswap = false, typename Type, std::size_t Size> static inline bool read(std::istream& stream, std::array<Type, Size>& container);
        template <bool Byteswap = false, unsigned int Current = 0, typename... Types> static inline bool read(std::istream& stream, std::tuple<Types...>& container);
        template <bool Byteswap = false> static inline bool rread(std::istream& stream, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <bool Byteswap = false, typename Type> static inline bool rread(std::istream& stream, Type* const first, Type* const last);
        template <bool Byteswap = false, typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline bool rread(std::istream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Write buffer
    /// \name           Write buffer
    //@{
    public:
        template <bool Byteswap = false> static inline char*& write(char*& buffer);
        template <bool Byteswap = false> static inline char*& write(char*& buffer, const std::tuple<>&);
        template <bool Byteswap = false, typename Type> static inline char*& write(char*& buffer, const Type& variable);
        template <bool Byteswap = false, typename Type, typename... Types> static inline char*& write(char*& buffer, const Type& variable, const Types&... variables);
        template <bool Byteswap = false, typename Type, std::size_t Size> static inline char*& write(char*& buffer, const std::array<Type, Size>& container);
        template <bool Byteswap = false, unsigned int Current = 0, typename... Types> static inline char*& write(char*& buffer, const std::tuple<Types...>& container);
        template <bool Byteswap = false> static inline char*& rwrite(char*& buffer, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <bool Byteswap = false, typename Type> static inline char*& rwrite(char*& buffer, const Type* const first, const Type* const last);
        template <bool Byteswap = false, typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline char*& rwrite(char*& buffer, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Read buffer
    /// \name           Read buffer
    //@{
    public:
        template <bool Byteswap = false> static inline char*& read(char*& buffer);
        template <bool Byteswap = false> static inline char*& read(char*& buffer, std::tuple<>&);
        template <bool Byteswap = false, typename Type> static inline char*& read(char*& buffer, Type& variable);
        template <bool Byteswap = false, typename Type, typename... Types> static inline char*& read(char*& buffer, Type& variable, Types&... variables);
        template <bool Byteswap = false, typename Type, std::size_t Size> static inline char*& read(char*& buffer, std::array<Type, Size>& container);
        template <bool Byteswap = false, unsigned int Current = 0, typename... Types> static inline char*& read(char*& buffer, std::tuple<Types...>& container);
        template <bool Byteswap = false> static inline char*& rread(char*& buffer, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <bool Byteswap = false, typename Type> static inline char*& rread(char*& buffer, Type* const first, Type* const last);
        template <bool Byteswap = false, typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline char*& rread(char*& buffer, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}
    
    // Nullification
    /// \name           Nullification
    //@{
    public:
        static inline unsigned int nullify();
        static inline unsigned int nullify(std::tuple<>&);
        template <typename Type> static inline unsigned int nullify(Type& variable);
        template <typename Type, typename... Types> static inline unsigned int nullify(Type& variable, Types&... variables);
        template <typename Type, std::size_t Size> static inline unsigned int nullify(std::array<Type, Size>& container);
        template <unsigned int Current = 0, typename... Types> static inline unsigned int nullify(std::tuple<Types...>& container);
        static inline unsigned long long int rnullify(const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <typename Type> static inline unsigned long long int rnullify(Type* const first, Type* const last);
        template <typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline unsigned long long int rnullify(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Equalization
    /// \name           Equalization
    //@{
    public:
        template <typename Reference> static inline unsigned int equalize(const Reference& reference);
        template <typename Reference> static inline unsigned int equalize(const Reference& reference, std::tuple<>&);
        template <typename Reference, typename Type> static inline unsigned int equalize(const Reference& reference, Type& variable);
        template <typename Reference, typename Type, typename... Types> static inline unsigned int equalize(const Reference& reference, Type& variable, Types&... variables);
        template <typename Reference, typename Type, std::size_t Size> static inline unsigned int equalize(const Reference& reference, std::array<Type, Size>& container);
        template <unsigned int Current = 0, typename Reference, typename... Types> static inline unsigned int equalize(const Reference& reference, std::tuple<Types...>& container);
        template <typename Reference> static inline unsigned long long int requalize(const Reference& reference, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <typename Reference, typename Type> static inline unsigned long long int requalize(const Reference& reference, Type* const first, Type* const last);
        template <typename Reference, typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static inline unsigned long long int requalize(const Reference& reference, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Hexification
    /// \name           Hexification
    //@{
    public:
        template <bool Byteswap = false, bool Uppercase = false> static std::string hexify();
        template <bool Byteswap = false, bool Uppercase = false> static std::string hexify(const std::tuple<>&);
        template <bool Byteswap = false, bool Uppercase = false, typename Type> static std::string hexify(const Type& variable);
        template <bool Byteswap = false, bool Uppercase = false, typename Type, typename... Types> static std::string hexify(const Type& variable, const Types&... variables);
        template <bool Byteswap = false, bool Uppercase = false, typename Type, std::size_t Size> static std::string hexify(const std::array<Type, Size>& container);
        template <bool Byteswap = false, bool Uppercase = false, unsigned int Current = 0, typename... Types> static std::string hexify(const std::tuple<Types...>& container);
        template <bool Byteswap = false, bool Uppercase = false> static std::string rhexify(const std::nullptr_t = nullptr, const std::nullptr_t = nullptr, const std::string& separator = " ");
        template <bool Byteswap = false, bool Uppercase = false, typename Type> static std::string rhexify(const Type* const first, const Type* const last, const std::string& separator = " ");
        template <bool Byteswap = false, bool Uppercase = false, typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static std::string rhexify(const Type& first, const Type& last, const std::string& separator = " ", typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Stringification
    /// \name           Stringification
    //@{
    public:
        template <unsigned int Base = 10, char Leading = char()> static std::string stringify();
        template <unsigned int Base = 10, char Leading = char()> static std::string stringify(const std::tuple<>&);
        template <unsigned int Base = 10, char Leading = char(), typename Type> static std::string stringify(const Type& variable);
        template <unsigned int Base = 10, char Leading = char(), typename Type, typename... Types> static std::string stringify(const Type& variable, const Types&... variables);
        template <unsigned int Base = 10, char Leading = char(), typename Type, std::size_t Size> static std::string stringify(const std::array<Type, Size>& container);
        template <unsigned int Base = 10, char Leading = char(), unsigned int Current = 0, typename... Types> static std::string stringify(const std::tuple<Types...>& container);
        template <unsigned int Base = 10, char Leading = char()> static std::string rstringify(const std::nullptr_t = nullptr, const std::nullptr_t = nullptr, const std::string& separator = " ");
        template <unsigned int Base = 10, char Leading = char(), typename Type> static std::string rstringify(const Type* const first, const Type* const last, const std::string& separator = " ");
        template <unsigned int Base = 10, char Leading = char(), typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static std::string rstringify(const Type& first, const Type& last, const std::string& separator = " ", typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Print
    /// \name           Print
    //@{
    public:
        static bool print(std::ostream& stream);
        static bool print(std::ostream& stream, const std::tuple<>&);
        template <typename Type> static bool print(std::ostream& stream, const Type& variable);
        template <typename Type, typename... Types> static bool print(std::ostream& stream, const Type& variable, const Types&... variables);
        template <typename Type, std::size_t Size> static bool print(std::ostream& stream, const std::array<Type, Size>& container);
        template <unsigned int Current = 0, typename... Types> static bool print(std::ostream& stream, const std::tuple<Types...>& container);
        static bool rprint(std::ostream& stream, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <typename Type> static bool rprint(std::ostream& stream, const Type* const first, const Type* const last);
        template <typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static bool rprint(std::ostream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}

    // Scan
    /// \name           Scan
    //@{
    public:
        static bool scan(std::istream& stream);
        static bool scan(std::istream& stream, std::tuple<>&);
        template <typename Type> static bool scan(std::istream& stream, Type& variable);
        template <typename Type, typename... Types> static bool scan(std::istream& stream, Type& variable, Types&... variables);
        template <typename Type, std::size_t Size> static bool scan(std::istream& stream, std::array<Type, Size>& container);
        template <unsigned int Current = 0, typename... Types> static bool scan(std::istream& stream, std::tuple<Types...>& container);
        static bool rscan(std::istream& stream, const std::nullptr_t = nullptr, const std::nullptr_t = nullptr);
        template <typename Type> static bool rscan(std::istream& stream, Type* const first, Type* const last);
        template <typename Type, class = typename std::enable_if<!std::is_pointer<Type>::value>::type> static bool rscan(std::istream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category* = nullptr);
    //@}
    
    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// -------------------------------- UTILITIES ------------------------------- //
// Is not an array
/// \brief          Is not an array.
/// \details        Returns false as the passed argument is not a standard 
///                 array.
/// \return         Whether the argument is a standard array.
template <typename Type> 
constexpr bool DataHandler::array(const Type&)
{
    return false;
}

// Is an array
/// \brief          Is an array.
/// \details        Returns true as the passed argument is a standard array.
/// \return         Whether the argument is a standard array.
template <typename Type, std::size_t Size>
constexpr bool DataHandler::array(const std::array<Type, Size>&)
{
    return true;
}

// Is not a tuple
/// \brief          Is not a tuple.
/// \details        Returns false as the passed argument is not a standard 
///                 tuple.
/// \return         Whether the argument is a standard tuple.
template <typename Type> 
constexpr bool DataHandler::tuple(const Type&)
{
    return false;
}

// Is a tuple
/// \brief          Is a tuple.
/// \details        Returns true as the passed argument is a standard tuple.
/// \return         Whether the argument is a standard tuple.
template <typename... Types>
constexpr bool DataHandler::tuple(const std::tuple<Types...>&)
{
    return true;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- SIZE ---------------------------------- //
// Number of elements of several types
/// \brief          Number of elements of several types.
/// \details        Computes the number of elements of the passed types. Note 
///                 that if the passed type is a templated type, like a tuple, 
///                 the inner arguments are not extracted, and the result is 
///                 simply the same as <tt>sizeof...()</tt> .
/// \return         Number of elements.
template<typename... Types> 
constexpr unsigned int DataHandler::size()
{
    return sizeof...(Types);
}

// Number of elements of an empty tuple
/// \brief          Number of elements of an empty tuple.
/// \details        Function overload that does nothing.
/// \return         Number of elements.
constexpr unsigned int DataHandler::size(const std::tuple<>&)
{
    return size();
}

// Number of elements of a single variable
/// \brief          Number of elements of a single variable.
/// \details        Computes the number of elements of the passed variable.
/// \tparam         Type (Variable type.)
/// \param[in]      variable Variable.
/// \return         Number of elements.
template <typename Type>
constexpr unsigned int DataHandler::size(const Type& variable)
{
    return size<Type>()+sizeof(variable)*0;
}

// Number of elements of several variables
/// \brief          Number of elements of several variables.
/// \details        Computes the number of elements of the passed variables.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      variable Variable.
/// \param[in]      variables Variables.
/// \return         Number of elements.
template <typename Type, typename... Types> 
constexpr unsigned int DataHandler::size(const Type& variable, const Types&... variables)
{
    return size(variable)+size(variables...);
}

// Number of elements of an array
/// \brief          Number of elements of an array.
/// \details        Computes the number of elements of the passed array.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in]      container Array.
/// \return         Number of elements.
template <typename Type, std::size_t Size> 
constexpr unsigned int DataHandler::size(const std::array<Type, Size>& container)
{
    return (Size > 0) ? (Size*size(container[0])) : (size());
}

// Number of elements of a tuple
/// \brief          Number of elements of a tuple.
/// \details        Computes the number of elements of the passed tuple.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in]      container Tuple.
/// \return         Number of elements.
template <unsigned int Current, typename... Types> 
constexpr unsigned int DataHandler::size(const std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (size(std::get<Current>(container))+size<((Current+1) < sizeof...(Types))*(Current+1)>(container)) : (size(std::get<Current>(container)));
}

// Number of elements of an empty range
/// \brief          Number of elements of an empty range.
/// \details        Function overload that does nothing.
/// \return         Number of elements.
inline unsigned long long int DataHandler::rsize(const std::nullptr_t, const std::nullptr_t)
{
    return size();
}

// Number of elements of a range between pointers
/// \brief          Number of elements of a range between pointers.
/// \details        Computes the number of elements of the range between the 
///                 passed pointers.
/// \tparam         Type (Pointer type.)
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Number of elements.
template <typename Type>
inline unsigned long long int DataHandler::rsize(const Type* const first, const Type* const last)
{
    return (first <= last) ? (static_cast<unsigned long long int>((last-first)*size(*first))) : (static_cast<unsigned long long int>(size()));
}

// Number of elements of a range between iterators
/// \brief          Number of elements of a range between iterators.
/// \details        Computes the number of elements of the range between the 
///                 passed iterators.
/// \tparam         Type (Iterator type.)
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Number of elements.
template <typename Type, class>
inline unsigned long long int DataHandler::rsize(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    return (first <= last) ? (static_cast<unsigned long long int>((last-first)*size(*first))) : (static_cast<unsigned long long int>(size()));
}
// -------------------------------------------------------------------------- //



// -------------------------------- BYTESIZE -------------------------------- //
// Size in bytes of a single type
/// \brief          Size in bytes of a single type.
/// \details        Computes the size in bytes of the passed type. The returned
///                 size is either the result of <tt>sizeof()</tt> or 0 if the
///                 type is void. Note that if the passed type is a templated
///                 type, like a tuple, the inner arguments are not extracted, 
///                 and the result is simply the same as <tt>sizeof()</tt> .
/// \return         Size in bytes.
template <typename Type> 
constexpr unsigned int DataHandler::bytesize(const std::nullptr_t)
{
    return (!std::is_void<Type>::value) ? (sizeof(typename std::conditional<!std::is_void<Type>::value, Type, char>::type)) : (0);
}

// Size in bytes of several types
/// \brief          Size in bytes of a single type.
/// \details        Computes the size in bytes of the passed types. The returned
///                 size is the sum of the sizes of each type.
/// \return         Size in bytes.
template <typename Type, typename... Types, class> 
constexpr unsigned int DataHandler::bytesize()
{
    return bytesize<Type>()+bytesize<Types...>();
}

// Size in bytes of an empty tuple
/// \brief          Size in bytes of an empty tuple.
/// \details        Function overload that does nothing.
/// \return         Size in bytes.
constexpr unsigned int DataHandler::bytesize(const std::tuple<>&)
{
    return bytesize();
}

// Size in bytes of a single variable
/// \brief          Size in bytes of a single variable.
/// \details        Computes the size in bytes of the passed variable.
/// \tparam         Type (Variable type.)
/// \param[in]      variable Variable.
/// \return         Size in bytes.
template <typename Type>
constexpr unsigned int DataHandler::bytesize(const Type& variable)
{
    return bytesize<Type>()+sizeof(variable)*0;
}

// Size in bytes of several variables
/// \brief          Size in bytes of several variables.
/// \details        Computes the sum of the size in bytes of the passed 
///                 variables.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      variable Variable.
/// \param[in]      variables Variables.
/// \return         Size in bytes.
template <typename Type, typename... Types> 
constexpr unsigned int DataHandler::bytesize(const Type& variable, const Types&... variables)
{
    return bytesize(variable)+bytesize(variables...);
}

// Size in bytes of an array
/// \brief          Size in bytes of an array.
/// \details        Computes the sum of the size in bytes of elements of the 
///                 passed array.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in]      container Array.
/// \return         Size in bytes.
template <typename Type, std::size_t Size> 
constexpr unsigned int DataHandler::bytesize(const std::array<Type, Size>& container)
{
    return (Size > 0) ? (Size*bytesize(container[0])) : (bytesize());
}

// Size in bytes of a tuple
/// \brief          Size in bytes of a tuple.
/// \details        Computes the sum of the size in bytes of elements of the 
///                 passed tuple.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in]      container Tuple.
/// \return         Size in bytes.
template <unsigned int Current, typename... Types> 
constexpr unsigned int DataHandler::bytesize(const std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (bytesize(std::get<Current>(container))+bytesize<((Current+1) < sizeof...(Types))*(Current+1)>(container)) : (bytesize(std::get<Current>(container)));
}

// Size in bytes of an empty range
/// \brief          Size in bytes of an empty range.
/// \details        Function overload that does nothing.
/// \return         Size in bytes.
inline unsigned long long int DataHandler::rbytesize(const std::nullptr_t, const std::nullptr_t)
{
    return bytesize();
}

// Size in bytes of a range between pointers
/// \brief          Size in bytes of a range between pointers.
/// \details        Computes the sum of the size in bytes of the elements of 
///                 the range between the passed pointers.
/// \tparam         Type (Pointer type.)
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Size in bytes.
template <typename Type>
inline unsigned long long int DataHandler::rbytesize(const Type* const first, const Type* const last)
{
    return (first <= last) ? (static_cast<unsigned long long int>((last-first)*bytesize(*first))) : (static_cast<unsigned long long int>(bytesize()));
}

// Size in bytes of a range between iterators
/// \brief          Size in bytes of a range between iterators.
/// \details        Computes the sum of the size in bytes of the range between
///                 the passed iterators.
/// \tparam         Type (Iterator type.)
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Size in bytes.
template <typename Type, class>
inline unsigned long long int DataHandler::rbytesize(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    return (first <= last) ? (static_cast<unsigned long long int>((last-first)*bytesize(*first))) : (static_cast<unsigned long long int>(bytesize()));
}
// -------------------------------------------------------------------------- //



// -------------------------------- BYTESWAP -------------------------------- //
// Do not swap bytes
/// \brief          Do not swap bytes.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Do not swap endianness if false.
/// \tparam         Types (Variadic types.)
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename... Types, class> 
inline unsigned int DataHandler::byteswap(const Types&...)
{
    return 0;
}

// Do not swap bytes of a range
/// \brief          Do not swap bytes of a range.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Do not swap endianness if false.
/// \tparam         Types (Variadic types.)
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename... Types, class>
inline unsigned int DataHandler::rbyteswap(const Types&...)
{
    return 0;
}

// Swap bytes of nothing
/// \brief          Swap bytes of nothing.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness if true.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, class>
inline unsigned int DataHandler::byteswap()
{
    return 0;
}

// Swap bytes of an empty tuple
/// \brief          Swap bytes of an empty tuple.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness if true.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, class>
inline unsigned int DataHandler::byteswap(std::tuple<>&)
{
    return byteswap<Byteswap>();
}

// Swap bytes of a single variable
/// \brief          Swap bytes of a single variable.
/// \details        Inverts the order of bytes of the passed variable to change 
///                 the endianness.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Type (Variable type.)
/// \param[in,out]  variable Variable.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename Type, class>
inline unsigned int DataHandler::byteswap(Type& variable)
{
    static const int half = sizeof(variable)/2;
    static const int end = sizeof(variable)-1;
    char* buffer = reinterpret_cast<char*>(&variable);
    char temporary;
    for (int i = 0; i < half; ++i) {
        temporary = buffer[i];
        buffer[i] = buffer[end-i];
        buffer[end-i] = temporary;
    }
    return 1;
}

// Swap bytes of several variables
/// \brief          Swap bytes of several variables.
/// \details        Inverts the order of bytes of the passed variables to change 
///                 the endianness.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  variable Variable.
/// \param[in,out]  variables Variables.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename Type, typename... Types, class> 
inline unsigned int DataHandler::byteswap(Type& variable, Types&... variables)
{
    return byteswap<Byteswap>(variable)+byteswap<Byteswap>(variables...);
}

// Swap bytes of an array
/// \brief          Swap bytes of an array.
/// \details        Inverts the order of bytes of each element of the passed 
///                 array to change the endianness.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  container Array.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename Type, std::size_t Size, class> 
inline unsigned int DataHandler::byteswap(std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        byteswap<Byteswap>(container[i]);
    }
    return Size;
}

// Swap bytes of a tuple
/// \brief          Swap bytes of a tuple.
/// \details        Inverts the order of bytes of each element of the passed 
///                 tuple to change the endianness.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in, out] container Tuple.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, unsigned int Current, typename... Types, class> 
inline unsigned int DataHandler::byteswap(std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (byteswap<Byteswap>(std::get<Current>(container))+byteswap<Byteswap, ((Current+1) < sizeof...(Types))*(Current+1)>(container)) : (byteswap<Byteswap>(std::get<Current>(container)));
}

// Swap bytes of an empty range
/// \brief          Swap bytes of an empty range.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness if true.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, class>
inline unsigned long long int DataHandler::rbyteswap(const std::nullptr_t, const std::nullptr_t)
{
    return byteswap<Byteswap>();
}

// Swap bytes of a range between pointers
/// \brief          Swap bytes of a range between pointers.
/// \details        Inverts the order of bytes of each element of the range
///                 between the passed pointers.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Type (Pointer type.)
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename Type, class>
inline unsigned long long int DataHandler::rbyteswap(Type* const first, Type* const last)
{
    unsigned long long int n = 0;
    for (Type* ptr = first; ptr < last; ++ptr) {
        n += byteswap<Byteswap>(*ptr);
    }
    return n;
}

// Swap bytes of a range between iterators
/// \brief          Swap bytes of a range between iterators.
/// \details        Inverts the order of bytes of each element of the range
///                 between the passed iterators.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Type (Iterator type.)
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Number of elements successfully byteswapped.
template <bool Byteswap, typename Type, class>
inline unsigned long long int DataHandler::rbyteswap(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    unsigned long long int n = 0;
    for (Type it = first; it != last; ++it) {
        n += byteswap<Byteswap>(*it);
    }
    return n;
}
// -------------------------------------------------------------------------- //



// ------------------------------ WRITE STREAM ------------------------------ //
// Write nothing to stream
/// \brief          Write nothing to stream.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \param[in,out]  stream Output stream.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap>
inline bool DataHandler::write(std::ostream& stream)
{
    return stream.good();
}

// Write an empty tuple to stream
/// \brief          Write an empty tuple to stream.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \param[in,out]  stream Output stream.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap>
inline bool DataHandler::write(std::ostream& stream, const std::tuple<>&)
{
    return write<Byteswap>(stream);
}

// Write a single variable to stream
/// \brief          Write a single variable to stream.
/// \details        Writes the next variable to the stream and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Variable type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      variable Variable.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type>
inline bool DataHandler::write(std::ostream& stream, const Type& variable)
{
    Type temporary = (Byteswap) ? (variable) : (Type());
    if (Byteswap) {
        byteswap<Byteswap>(temporary);
        stream.write(reinterpret_cast<char*>(&temporary), sizeof(temporary));
    } else {
        stream.write(reinterpret_cast<const char*>(&variable), sizeof(variable));
    }
    return stream.good();
}

// Write several variables to stream
/// \brief          Write several variables to stream.
/// \details        Writes the next variables to the stream and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Output stream.
/// \param[in]      variable Variable.
/// \param[in]      variables Variables.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type, typename... Types> 
inline bool DataHandler::write(std::ostream& stream, const Type& variable, const Types&... variables)
{
    write<Byteswap>(stream, variable);
    return write<Byteswap>(stream, variables...);
}

// Write an array to stream
/// \brief          Write an array to stream.
/// \details        Writes each element of the passed array to the stream and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  stream Output stream.
/// \param[in]      container Array.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type, std::size_t Size> 
inline bool DataHandler::write(std::ostream& stream, const std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        write<Byteswap>(stream, container[i]);
    }
    return stream.good();
}

// Write a tuple to stream
/// \brief          Write a tuple to stream.
/// \details        Writes each element of the passed tuple to the stream and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Output stream.
/// \param[in]      container Tuple.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, unsigned int Current, typename... Types> 
inline bool DataHandler::write(std::ostream& stream, const std::tuple<Types...>& container)
{
    write<Byteswap>(stream, std::get<Current>(container));
    return ((Current+1) < sizeof...(Types)) ? (write<Byteswap, ((Current+1) < sizeof...(Types))*(Current+1)>(stream, container)) : (write(stream));
}

// Write an empty range to stream
/// \brief          Write an empty range to stream.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \param[in,out]  stream Output stream.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap>
inline bool DataHandler::rwrite(std::ostream& stream, const std::nullptr_t, const std::nullptr_t)
{
    return write<Byteswap>(stream);
}

// Write a range between pointers to stream
/// \brief          Write a range between pointers to stream.
/// \details        Writes each element of the range between the passed pointers 
///                 to the stream and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Pointer type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type>
inline bool DataHandler::rwrite(std::ostream& stream, const Type* const first, const Type* const last)
{
    if (std::is_fundamental<Type>::value) {
        if (Byteswap) {
            for (const Type* ptr = first; ptr < last; ++ptr) {
                write<Byteswap>(stream, *ptr);
            }
        } else {
            stream.write(reinterpret_cast<const char*>(first), (last-first)*sizeof(*first));
        }
    } else {
        for (const Type* ptr = first; ptr < last; ++ptr) {
            write<Byteswap>(stream, *ptr);
        }
    }
    return stream.good();
}

// Write a range between iterators to stream
/// \brief          Write a range between iterators to stream.
/// \details        Writes each element of the range between the passed 
///                 iterators to the stream and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Iterator type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type, class>
inline bool DataHandler::rwrite(std::ostream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    for (Type it = first; it != last; ++it) {
        write<Byteswap>(stream, *it);
    }
    return stream.good();
}
// -------------------------------------------------------------------------- //



// ------------------------------- READ STREAM ------------------------------ //
// Read nothing from stream
/// \brief          Read nothing from stream.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \param[in,out]  stream Input stream.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap>
inline bool DataHandler::read(std::istream& stream)
{
    return stream.good();
}

// Read an empty tuple from stream
/// \brief          Read an empty tuple from stream.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \param[in,out]  stream Input stream.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap>
inline bool DataHandler::read(std::istream& stream, std::tuple<>&)
{
    return read<Byteswap>(stream);
}

// Read a single variable from stream
/// \brief          Read a single variable from stream.
/// \details        Reads the next variable from the stream and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Variable type.)
/// \param[in,out]  stream Input stream.
/// \param[out]     variable Variable.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type>
inline bool DataHandler::read(std::istream& stream, Type& variable)
{
// HERE
    //stream.read(reinterpret_cast<char*>(&variable), sizeof(variable));
    stream.read(reinterpret_cast<char*>(&variable), SIZEOFOCTREE);
    byteswap<Byteswap>(variable);
    return stream.good();
}

// Read several variables from stream
/// \brief          Read several variables from stream.
/// \details        Reads the next variables from the stream and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Input stream.
/// \param[out]     variable Variable.
/// \param[out]     variables Variables.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type, typename... Types> 
inline bool DataHandler::read(std::istream& stream, Type& variable, Types&... variables)
{
    read<Byteswap>(stream, variable);
    return read<Byteswap>(stream, variables...);
}

// Read an array from stream
/// \brief          Read an array from stream.
/// \details        Reads each element of the passed array from the stream and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  stream Input stream.
/// \param[out]     container Array.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type, std::size_t Size> 
inline bool DataHandler::read(std::istream& stream, std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        read<Byteswap>(stream, container[i]);
    }
    return stream.good();
}

// Read a tuple from stream
/// \brief          Read a tuple from stream.
/// \details        Reads each element of the passed tuple from the stream and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Input stream.
/// \param[out]     container Tuple.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, unsigned int Current, typename... Types> 
inline bool DataHandler::read(std::istream& stream, std::tuple<Types...>& container)
{
    read<Byteswap>(stream, std::get<Current>(container));
    return ((Current+1) < sizeof...(Types)) ? (read<Byteswap, ((Current+1) < sizeof...(Types))*(Current+1)>(stream, container)) : (read(stream));
}

// Read an empty range from stream
/// \brief          Read an empty range from stream.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \param[in,out]  stream Input stream.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap>
inline bool DataHandler::rread(std::istream& stream, const std::nullptr_t, const std::nullptr_t)
{
    return read<Byteswap>(stream);
}

// Read a range between pointers from stream
/// \brief          Read a range between pointers from stream.
/// \details        Reads each element of the range between the passed pointers 
///                 from the stream and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Pointer type.)
/// \param[in,out]  stream Input stream.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type>
inline bool DataHandler::rread(std::istream& stream, Type* const first, Type* const last)
{ 
    if (std::is_fundamental<Type>::value) {
        stream.read(reinterpret_cast<char*>(first), (last-first)*sizeof(*first));
        rbyteswap<Byteswap>(first, last);
    } else {
	int g = 0;
        for (Type* ptr = first; ptr < last; ++ptr) {
            read<Byteswap>(stream, *ptr);
        }
	g++;
    } 
    return stream.good();
}

// Read a range between iterators from stream
/// \brief          Read a range between iterators from stream.
/// \details        Reads each element of the range between the passed iterators 
///                 from the stream and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Iterator type.)
/// \param[in,out]  stream Input stream.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         True if the stream has no failure, false otherwise.
template <bool Byteswap, typename Type, class>
inline bool DataHandler::rread(std::istream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    for (Type it = first; it != last; ++it) {
        read<Byteswap>(stream, *it);
    }
    return stream.good();
}
// -------------------------------------------------------------------------- //



// ------------------------------ WRITE BUFFER ------------------------------ //
// Write nothing to buffer
/// \brief          Write nothing to buffer.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \param[in,out]  buffer Output buffer.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap>
inline char*& DataHandler::write(char*& buffer)
{
    return buffer;
}

// Write an empty tuple to buffer
/// \brief          Write an empty tuple to buffer.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \param[in,out]  buffer Output buffer.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap>
inline char*& DataHandler::write(char*& buffer, const std::tuple<>&)
{
    return write<Byteswap>(buffer);
}

// Write a single variable to buffer
/// \brief          Write a single variable to buffer.
/// \details        Writes the next variable to the buffer and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Variable type.)
/// \param[in,out]  buffer Output buffer.
/// \param[in]      variable Variable.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type>
inline char*& DataHandler::write(char*& buffer, const Type& variable)
{
    *reinterpret_cast<Type*>(buffer) = variable;
    byteswap<Byteswap>(*reinterpret_cast<Type*>(buffer));
    return buffer += sizeof(variable);
}

// Write several variables to buffer
/// \brief          Write several variables to buffer.
/// \details        Writes the next variables to the buffer and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  buffer Output buffer.
/// \param[in]      variable Variable.
/// \param[in]      variables Variables.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type, typename... Types> 
inline char*& DataHandler::write(char*& buffer, const Type& variable, const Types&... variables)
{
    write<Byteswap>(buffer, variable);
    return write<Byteswap>(buffer, variables...);
}

// Write an array to buffer
/// \brief          Write an array to buffer.
/// \details        Writes each element of the passed array to the buffer and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  buffer Output buffer.
/// \param[in]      container Array.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type, std::size_t Size> 
inline char*& DataHandler::write(char*& buffer, const std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        write<Byteswap>(buffer, container[i]);
    }
    return buffer;
}

// Write a tuple to buffer
/// \brief          Write a tuple to buffer.
/// \details        Writes each element of the passed tuple to the buffer and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  buffer Output buffer.
/// \param[in]      container Tuple.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, unsigned int Current, typename... Types> 
inline char*& DataHandler::write(char*& buffer, const std::tuple<Types...>& container)
{
    write<Byteswap>(buffer, std::get<Current>(container));
    return ((Current+1) < sizeof...(Types)) ? (write<Byteswap, ((Current+1) < sizeof...(Types))*(Current+1)>(buffer, container)) : (write(buffer));
}

// Write an empty range to buffer
/// \brief          Write an empty range to buffer.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \param[in,out]  buffer Output buffer.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap>
inline char*& DataHandler::rwrite(char*& buffer, const std::nullptr_t, const std::nullptr_t)
{
    return write<Byteswap>(buffer);
}

// Write a range between pointers to buffer
/// \brief          Write a range between pointers to buffer.
/// \details        Writes each element of the range between the passed pointers 
///                 to the buffer and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Pointer type.)
/// \param[in,out]  buffer Output buffer.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type>
inline char*& DataHandler::rwrite(char*& buffer, const Type* const first, const Type* const last)
{
    if (std::is_fundamental<Type>::value) {
        std::memcpy(buffer, first, (last-first)*sizeof(*first));
        rbyteswap<Byteswap>(reinterpret_cast<Type*>(buffer), reinterpret_cast<Type*>(buffer)+(last-first)*sizeof(*first));
        buffer += (last-first)*sizeof(*first);
    } else {
        for (const Type* ptr = first; ptr < last; ++ptr) {
            write<Byteswap>(buffer, *ptr);
        }
    }
    return buffer;
}

// Write a range between iterators to buffer
/// \brief          Write a range between iterators to buffer.
/// \details        Writes each element of the range between the passed 
///                 iterators to the buffer and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness for writing if true.
/// \tparam         Type (Iterator type.)
/// \param[in,out]  buffer Output buffer.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type, class>
inline char*& DataHandler::rwrite(char*& buffer, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    for (Type it = first; it != last; ++it) {
        write<Byteswap>(buffer, *it);
    }
    return buffer;
}
// -------------------------------------------------------------------------- //



// ------------------------------- READ BUFFER ------------------------------ //
// Read nothing from buffer
/// \brief          Read nothing from buffer.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \param[in,out]  buffer Input buffer.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap>
inline char*& DataHandler::read(char*& buffer)
{
    return buffer;
}

// Read an empty tuple from buffer
/// \brief          Read an empty tuple from buffer.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \param[in,out]  buffer Input buffer.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap>
inline char*& DataHandler::read(char*& buffer, std::tuple<>&)
{
    return read<Byteswap>(buffer);
}

// Read a single variable from buffer
/// \brief          Read a single variable from buffer.
/// \details        Reads the next variable from the buffer and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Variable type.)
/// \param[in,out]  buffer Input buffer.
/// \param[out]     variable Variable.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type>
inline char*& DataHandler::read(char*& buffer, Type& variable)
{
    variable = *reinterpret_cast<Type*>(buffer);
    byteswap<Byteswap>(variable);
    return buffer += sizeof(variable);
}

// Read several variables from buffer
/// \brief          Read several variables from buffer.
/// \details        Reads the next variables from the buffer and swap bytes if 
///                 necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  buffer Input buffer.
/// \param[out]     variable Variable.
/// \param[out]     variables Variables.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type, typename... Types> 
inline char*& DataHandler::read(char*& buffer, Type& variable, Types&... variables)
{
    read<Byteswap>(buffer, variable);
    return read<Byteswap>(buffer, variables...);
}

// Read an array from buffer
/// \brief          Read an array from buffer.
/// \details        Reads each element of the passed array from the buffer and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  buffer Input buffer.
/// \param[out]     container Array.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type, std::size_t Size> 
inline char*& DataHandler::read(char*& buffer, std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        read<Byteswap>(buffer, container[i]);
    }
    return buffer;
}

// Read a tuple from buffer
/// \brief          Read a tuple from buffer.
/// \details        Reads each element of the passed tuple from the buffer and 
///                 swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  buffer Input buffer.
/// \param[out]     container Tuple.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, unsigned int Current, typename... Types> 
inline char*& DataHandler::read(char*& buffer, std::tuple<Types...>& container)
{
    read<Byteswap>(buffer, std::get<Current>(container));
    return ((Current+1) < sizeof...(Types)) ? (read<Byteswap, ((Current+1) < sizeof...(Types))*(Current+1)>(buffer, container)) : (read(buffer));
}

// Read an empty range from buffer
/// \brief          Read an empty range from buffer.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \param[in,out]  buffer Input buffer.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap>
inline char*& DataHandler::rread(char*& buffer, const std::nullptr_t, const std::nullptr_t)
{
    return read<Byteswap>(buffer);
}

// Read a range between pointers from buffer
/// \brief          Read a range between pointers from buffer.
/// \details        Reads each element of the range between the passed pointers 
///                 from the buffer and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Pointer type.)
/// \param[in,out]  buffer Input buffer.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type>
inline char*& DataHandler::rread(char*& buffer, Type* const first, Type* const last)
{
    if (std::is_fundamental<Type>::value) {
        std::memcpy(first, buffer, (last-first)*sizeof(*first));
        rbyteswap<Byteswap>(first, last);
        buffer += (last-first)*sizeof(*first);
    } else {
        for (Type* ptr = first; ptr < last; ++ptr) {
            read<Byteswap>(buffer, *ptr);
        }
    }
    return buffer;
}

// Read a range between iterators from buffer
/// \brief          Read a range between iterators from buffer.
/// \details        Reads each element of the range between the passed iterators 
///                 from the buffer and swap bytes if necessary.
/// \tparam         Byteswap Swap endianness after reading if true.
/// \tparam         Type (Iterator type.)
/// \param[in,out]  buffer Input buffer.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Pointer to the new position in the buffer.
template <bool Byteswap, typename Type, class>
inline char*& DataHandler::rread(char*& buffer, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    for (Type it = first; it != last; ++it) {
        read<Byteswap>(buffer, *it);
    }
    return buffer;
}
// -------------------------------------------------------------------------- //



// ------------------------------ NULLIFICATION ----------------------------- //
// Nullify nothing
/// \brief          Nullify nothing.
/// \details        Function overload that does nothing.
/// \return         Number of elements successfully set to zero.
inline unsigned int DataHandler::nullify()
{
    return 0;
}

// Nullify an empty tuple
/// \brief          Nullify an empty tuple.
/// \details        Function overload that does nothing.
/// \return         Number of elements successfully set to zero.
inline unsigned int DataHandler::nullify(std::tuple<>&)
{
    return nullify();
}

// Nullify a single variable
/// \brief          Nullify a single variable.
/// \details        Calls the constructor of the passed variable to nullify it.
/// \tparam         Type (Variable type.)
/// \param[out]     variable Variable.
/// \return         Number of elements successfully set to zero.
template <typename Type>
inline unsigned int DataHandler::nullify(Type& variable)
{
    variable = Type();
    return 1;
}

// Nullify several variables
/// \brief          Nullify several variables.
/// \details        Calls the constructor of the passed variables to nullify 
///                 them.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[out]     variable Variable.
/// \param[out]     variables Variables.
/// \return         Number of elements successfully set to zero.
template <typename Type, typename... Types> 
inline unsigned int DataHandler::nullify(Type& variable, Types&... variables)
{
    return nullify(variable)+nullify(variables...);
}

// Nullify an array
/// \brief          Nullify an array.
/// \details        Calls the constructor of each element of the passed array to
///                 nullify it.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[out]     container Array.
/// \return         Number of elements successfully set to zero.
template <typename Type, std::size_t Size> 
inline unsigned int DataHandler::nullify(std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        nullify(container[i]);
    }
    return Size;
}

// Nullify a tuple
/// \brief          Nullify a tuple.
/// \details        Calls the constructor of each element of the passed tuple to
///                 nullify it.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[out]     container Tuple.
/// \return         Number of elements successfully set to zero.
template <unsigned int Current, typename... Types> 
inline unsigned int DataHandler::nullify(std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (nullify(std::get<Current>(container))+nullify<((Current+1) < sizeof...(Types))*(Current+1)>(container)) : (nullify(std::get<Current>(container)));
}

// Nullify an empty range
/// \brief          Nullify an empty range.
/// \details        Function overload that does nothing.
/// \return         Number of elements successfully set to zero.
inline unsigned long long int DataHandler::rnullify(const std::nullptr_t, const std::nullptr_t)
{
    return nullify();
}

// Nullify a range between pointers
/// \brief          Nullify a range between pointers.
/// \details        Calls the constructor of each element of the range between
///                 the passed pointers.
/// \tparam         Type (Pointer type.)
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Number of elements successfully set to zero.
template <typename Type>
inline unsigned long long int DataHandler::rnullify(Type* const first, Type* const last)
{
    unsigned long long int n = 0;
    for (Type* ptr = first; ptr < last; ++ptr) {
        n += nullify(*ptr);
    }
    return n;
}

// Nullify a range between iterators
/// \brief          Nullify a range between iterators.
/// \details        Calls the constructor of each element of the range between
///                 the passed iterators.
/// \tparam         Type (Iterator type.)
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Number of elements successfully set to zero.
template <typename Type, class>
inline unsigned long long int DataHandler::rnullify(const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    unsigned long long int n = 0;
    for (Type it = first; it != last; ++it) {
        n += nullify(*it);
    }
    return n;
}
// -------------------------------------------------------------------------- //



// ------------------------------ EQUALIZATION ------------------------------ //
// Equalize nothing
/// \brief          Equalize nothing.
/// \details        Function overload that does nothing.
/// \tparam         Reference (Reference type.)
/// \param[in]      reference Reference value.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference>
inline unsigned int DataHandler::equalize(const Reference& reference)
{
    (void)(reference);
    return 0;
}

// Equalize an empty tuple
/// \brief          Equalize an empty tuple.
/// \details        Function overload that does nothing.
/// \tparam         Reference (Reference type.)
/// \param[in]      reference Reference value.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference>
inline unsigned int DataHandler::equalize(const Reference& reference, std::tuple<>&)
{
    return equalize(reference);
}

// Equalize a single variable
/// \brief          Equalize a single variable.
/// \details        Set the passed variable to the reference value.
/// \tparam         Reference (Reference type.)
/// \tparam         Type (Variable type.)
/// \param[in]      reference Reference value.
/// \param[out]     variable Variable.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference, typename Type>
inline unsigned int DataHandler::equalize(const Reference& reference, Type& variable)
{
    variable = reference;
    return 1;
}

// Equalize several variables
/// \brief          Equalize several variables.
/// \details        Set the passed variables to the reference value.
/// \tparam         Reference (Reference type.)
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      reference Reference value.
/// \param[out]     variable Variable.
/// \param[out]     variables Variables.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference, typename Type, typename... Types> 
inline unsigned int DataHandler::equalize(const Reference& reference, Type& variable, Types&... variables)
{
    return equalize(reference, variable)+equalize(reference, variables...);
}

// Equalize an array
/// \brief          Equalize an array.
/// \details        Set each element of the passed array to the reference value.
/// \tparam         Reference (Reference type.)
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in]      reference Reference value.
/// \param[out]     container Array.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference, typename Type, std::size_t Size> 
inline unsigned int DataHandler::equalize(const Reference& reference, std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        equalize(reference, container[i]);
    }
    return Size;
}

// Equalize a tuple
/// \brief          Equalize a tuple.
/// \details        Set each element of the passed tuple to the reference value.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Reference (Reference type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      reference Reference value.
/// \param[out]     container Tuple.
/// \return         Number of elements successfully set to the reference value.
template <unsigned int Current, typename Reference, typename... Types> 
inline unsigned int DataHandler::equalize(const Reference& reference, std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (equalize(reference, std::get<Current>(container))+equalize<((Current+1) < sizeof...(Types))*(Current+1)>(reference, container)) : (equalize(reference, std::get<Current>(container)));
}

// Equalize an empty range
/// \brief          Equalize an empty range.
/// \details        Function overload that does nothing.
/// \tparam         Reference (Reference type.)
/// \param[in]      reference Reference value.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference>
inline unsigned long long int DataHandler::requalize(const Reference& reference, const std::nullptr_t, const std::nullptr_t)
{
    return equalize(reference);
}

// Equalize a range between pointers
/// \brief          Equalize a range between pointers.
/// \details        Set each element of the range between the passed pointers 
///                 to the reference value.
/// \tparam         Reference (Reference type.)
/// \tparam         Type (Pointer type.)
/// \param[in]      reference Reference value.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference, typename Type>
inline unsigned long long int DataHandler::requalize(const Reference& reference, Type* const first, Type* const last)
{
    unsigned long long int n = 0;
    for (Type* ptr = first; ptr < last; ++ptr) {
        n += equalize(reference, *ptr);
    }
    return n;
}

// Equalize a range between iterators
/// \brief          Equalize a range between iterators.
/// \details        Set each element of the range between the passed iterators 
///                 to the reference value.
/// \tparam         Reference (Reference type.)
/// \tparam         Type (Iterator type.)
/// \param[in]      reference Reference value.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         Number of elements successfully set to the reference value.
template <typename Reference, typename Type, class>
inline unsigned long long int DataHandler::requalize(const Reference& reference, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    unsigned long long int n = 0;
    for (Type it = first; it != last; ++it) {
        n += equalize(reference, *it);
    }
    return n;
}
// -------------------------------------------------------------------------- //



// ------------------------------ HEXIFICATION ------------------------------ //
// Convert nothing to a hexadecimal string
/// \brief          Convert nothing to a hexadecimal string.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase>
std::string DataHandler::hexify()
{
    return std::string();
}

// Convert an empty tuple to a hexadecimal string
/// \brief          Convert an empty tuple to a hexadecimal string.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase>
std::string DataHandler::hexify(const std::tuple<>&)
{
    return hexify<Byteswap, Uppercase>();
}

// Convert a single variable to a hexadecimal string
/// \brief          Convert a single variable to a hexadecimal string.
/// \details        Converts bytes of the passed variable to their hexadecimal
///                 representation and returns the associated string.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \tparam         Type (Variable type.)
/// \param[in]      variable Variable.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase, typename Type>
std::string DataHandler::hexify(const Type& variable)
{
    static const int half = sizeof(variable);
    static const int end = 2*sizeof(variable)-1;
    static const std::array<unsigned char, 16> alphabet = (Uppercase) ? (std::array<unsigned char, 16>({{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'}})) : (std::array<unsigned char, 16>({{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'}}));
    const unsigned char* const ptr = reinterpret_cast<unsigned char const*>(&variable);
    std::string result(end+1, 0);
    if (Byteswap) {
        for (int i = 0; i < half; ++i)
        {
            result[end-2*i-1] = alphabet[ptr[i]/16];
            result[end-2*i] = alphabet[ptr[i]%16];
        }
    } else {
        for (int i = 0; i < half; ++i)
        {
            result[2*i] = alphabet[ptr[i]/16];
            result[2*i+1] = alphabet[ptr[i]%16];
        }
    }
    return result;
}

// Convert several variables to a hexadecimal string
/// \brief          Convert several variables to a hexadecimal string.
/// \details        Converts bytes of the passed variables to their hexadecimal
///                 representation and returns the associated string.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in]      variable Variable.
/// \param[in]      variables Variables.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase, typename Type, typename... Types> 
std::string DataHandler::hexify(const Type& variable, const Types&... variables)
{
    return hexify<Byteswap, Uppercase>(variable)+std::string((bytesize(variable) > 0) ? (" ") : (""))+hexify<Byteswap, Uppercase>(variables...);
}

// Convert an array to a hexadecimal string
/// \brief          Convert an array to a hexadecimal string.
/// \details        Converts bytes of each element of the passed array to their 
///                 hexadecimal representation and returns the associated 
///                 string.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in]      container Array.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase, typename Type, std::size_t Size> 
std::string DataHandler::hexify(const std::array<Type, Size>& container)
{
    static const std::string separator(" ");
    std::string result;
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        result += hexify<Byteswap, Uppercase>(container[i])+separator;
    }
    if ((!result.empty()) && (result.size() >= separator.size())) {
        result.erase(result.size()-separator.size());
    }
    return result;
}

// Convert a tuple to a hexadecimal string
/// \brief          Convert a tuple to a hexadecimal string.
/// \details        Converts bytes of each element of the passed tuple to their 
///                 hexadecimal representation and returns the associated 
///                 string.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in]      container Tuple.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase, unsigned int Current, typename... Types> 
std::string DataHandler::hexify(const std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (hexify<Byteswap, Uppercase>(std::get<Current>(container))+std::string((bytesize(std::get<Current>(container)) > 0) ? (" ") : (""))+hexify<Byteswap, Uppercase, ((Current+1) < sizeof...(Types))*(Current+1)>(container)) : (hexify<Byteswap, Uppercase>(std::get<Current>(container)));
}

// Convert an empty range to a hexadecimal string
/// \brief          Convert an empty range to a hexadecimal string.
/// \details        Function overload that does nothing.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \param[in]      separator String used as a separator between elements.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase>
std::string DataHandler::rhexify(const std::nullptr_t, const std::nullptr_t, const std::string& separator)
{
    (void)(separator);
    return hexify<Byteswap, Uppercase>();
}

// Convert a range between pointers to a hexadecimal string
/// \brief          Convert a range between pointers to a hexadecimal string.
/// \details        Converts bytes of each element of the range between the
///                 passed pointers to their hexadecimal representation and 
///                 returns the associated string.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \tparam         Type (Pointer type.)
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \param[in]      separator String used as a separator between elements.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase, typename Type>
std::string DataHandler::rhexify(const Type* const first, const Type* const last, const std::string& separator)
{
    std::string result;
    for (const Type* ptr = first; ptr < last; ++ptr) {
        result += hexify<Byteswap, Uppercase>(*ptr)+separator;
    }
    if ((!result.empty()) && (result.size() >= separator.size())) {
        result.erase(result.size()-separator.size());
    }
    return result;
}

// Convert a range between iterators to a hexadecimal string
/// \brief          Convert a range between iterators to a hexadecimal string.
/// \details        Converts bytes of each element of the range between the
///                 passed iterators to their hexadecimal representation and 
///                 returns the associated string.
/// \tparam         Byteswap Swap bytes of hexadecimal conversion if true.
/// \tparam         Uppercase Force uppercase letters for hexadecimal 
///                 representation if true.
/// \tparam         Type (Iterator type.)
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \param[in]      separator String used as a separator between elements.
/// \return         Hexadecimal representation of the input.
template <bool Byteswap, bool Uppercase, typename Type, class>
std::string DataHandler::rhexify(const Type& first, const Type& last, const std::string& separator, typename std::iterator_traits<Type>::iterator_category*)
{
    std::string result;
    for (Type it = first; it != last; ++it) {
        result += hexify<Byteswap, Uppercase>(*it)+separator;
    }
    if ((!result.empty()) && (result.size() >= separator.size())) {
        result.erase(result.size()-separator.size());
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ----------------------------- STRINGIFICATION ---------------------------- //
// Convert nothing to a string
/// \brief          Convert nothing to a string.
/// \details        Function overload that does nothing.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading>
std::string DataHandler::stringify()
{
    static_assert((Base == 2) || (Base == 8) || (Base == 10) || (Base == 16), "ERROR = DataHandler::stringify() : incompatible base");
    static_assert((Leading == static_cast<char>(0)) || ((Leading >= static_cast<char>(9)) && (Leading <= static_cast<char>(13))) || ((Leading >= static_cast<char>(32)) && (Leading <= static_cast<char>(126))), "ERROR = DataHandler::stringify() : incompatible leading character");
    return std::string();
}

// Convert an empty tuple to a string
/// \brief          Convert an empty tuple to a string.
/// \details        Function overload that does nothing.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading>
std::string DataHandler::stringify(const std::tuple<>&)
{
    return stringify<Base, Leading>();
}

// Convert a single variable to a string
/// \brief          Convert a single variable to a string.
/// \details        Converts the passed variable to a string. The base of 
///                 conversion can be specified for integral numbers, and it
///                 is possible to fill-in the result with leading characters.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \tparam         Type (Variable type.)
/// \param[in]      variable Variable.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading, typename Type>
std::string DataHandler::stringify(const Type& variable)
{
    static const unsigned int max = (Base == 10) ? (std::max(static_cast<int>(std::to_string(std::numeric_limits<typename std::conditional<std::is_integral<Type>::value, Type, int>::type>::min()).size()), static_cast<int>(std::to_string(std::numeric_limits<typename std::conditional<std::is_integral<Type>::value, Type, int>::type>::max()).size())+1)) : (static_cast<int>((16*sizeof(variable))/Base));
    static const char zero = '0';
    static const char plus = '+';
    static const std::string signs("-+");
    std::string result(stringify<Base, Leading>());
    std::ostringstream oss;
    unsigned int mark = 0;
    if (std::is_integral<Type>::value) {
        oss<<((Base == 8) ? (std::oct) : ((Base == 16) ? (std::hex) : (std::dec)))<<(typename std::conditional<(std::is_integral<Type>::value) && (Base == 2), std::bitset<8*sizeof(variable)>, typename std::conditional<(sizeof(variable) > 1), Type, int>::type>::type(variable));
        result = oss.str();
        if (result.find_first_not_of(zero) != std::string::npos) {
            result.erase(0, result.find_first_not_of(zero));
        }
        if (result.empty()) {
            result.push_back(zero);
        }
        if ((Leading != char()) && (result.size() < max)) {
            if (Base == 10) {
                if (result.find_first_of(signs) == std::string::npos) {
                    result.insert(0, 1, plus);
                }
                if ((Leading >= static_cast<char>(48)) && (Leading <= static_cast<char>(57)) && (result.find_first_of(signs) != std::string::npos) && (result.find_first_of(signs)+1 < result.size())) {
                    mark = result.find_first_of(signs)+1;
                }
            }
            result.insert(mark, max-result.size(), Leading);
        }
    } else {
        oss<<variable;
        result = oss.str();
    }
    return result;
}

// Convert several variables to a string
/// \brief          Convert several variables to a string.
/// \details        Converts the passed variables to a string.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[out]     variable Variable.
/// \param[out]     variables Variables.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading, typename Type, typename... Types> 
std::string DataHandler::stringify(const Type& variable, const Types&... variables)
{
    return stringify<Base, Leading>(variable)+std::string((bytesize(variable) > 0) ? (" ") : (""))+stringify<Base, Leading>(variables...);
}

// Convert an array to a string
/// \brief          Convert an array to a string.
/// \details        Converts each element of the passed array to a string.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[out]     container Array.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading, typename Type, std::size_t Size> 
std::string DataHandler::stringify(const std::array<Type, Size>& container)
{
    static const std::string separator(" ");
    std::string result;
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        result += stringify<Base, Leading>(container[i])+separator;
    }
    if ((!result.empty()) && (result.size() >= separator.size())) {
        result.erase(result.size()-separator.size());
    }
    return result;
}

// Convert a tuple to a string
/// \brief          Convert a tuple to a string.
/// \details        Converts each element of the passed tuple to a string.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[out]     container Tuple.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading, unsigned int Current, typename... Types> 
std::string DataHandler::stringify(const std::tuple<Types...>& container)
{
    return ((Current+1) < sizeof...(Types)) ? (stringify<Base, Leading>(std::get<Current>(container))+std::string((bytesize(std::get<Current>(container)) > 0) ? (" ") : (""))+stringify<Base, Leading, ((Current+1) < sizeof...(Types))*(Current+1)>(container)) : (stringify<Base, Leading>(std::get<Current>(container)));
}

// Convert an empty range to a string
/// \brief          Convert an empty range to a string.
/// \details        Function overload that does nothing.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \param[in]      separator String used as a separator between elements.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading>
std::string DataHandler::rstringify(const std::nullptr_t, const std::nullptr_t, const std::string& separator)
{
    (void)(separator);
    return stringify<Base, Leading>();
}

// Convert a range between pointers to a string
/// \brief          Convert a range between pointers to a string.
/// \details        Converts each element of the range between the passed 
///                 pointers to a string.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \tparam         Type (Pointer type.)
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \param[in]      separator String used as a separator between elements.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading, typename Type>
std::string DataHandler::rstringify(const Type* const first, const Type* const last, const std::string& separator)
{
    std::string result;
    for (const Type* ptr = first; ptr < last; ++ptr) {
        result += stringify<Base, Leading>(*ptr)+separator;
    }
    if ((!result.empty()) && (result.size() >= separator.size())) {
        result.erase(result.size()-separator.size());
    }
    return result;
}

// Convert a range between iterators to a string
/// \brief          Convert a range between iterators to a string.
/// \details        Converts each element of the range between the passed 
///                 iterators to a string.
/// \tparam         Base Base representation for integers : 2, 8, 10 or 16.
/// \tparam         Leading Leading character for integers : 0, 9-13 or 32-126.
/// \tparam         Type (Iterator type.)
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \param[in]      separator String used as a separator between elements.
/// \return         String conversion of arguments.
template <unsigned int Base, char Leading, typename Type, class>
std::string DataHandler::rstringify(const Type& first, const Type& last, const std::string& separator, typename std::iterator_traits<Type>::iterator_category*)
{
    std::string result;
    for (Type it = first; it != last; ++it) {
        result += stringify<Base, Leading>(*it)+separator;
    }
    if ((!result.empty()) && (result.size() >= separator.size())) {
        result.erase(result.size()-separator.size());
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- PRINT --------------------------------- //
// Print nothing to stream
/// \brief          Print nothing to stream.
/// \details        Function overload that does nothing.
/// \param[in,out]  stream Output stream.
/// \return         True if the stream has no failure, false otherwise.
inline bool DataHandler::print(std::ostream& stream)
{
    return stream.good();
}

// Print an empty tuple to stream
/// \brief          Print an empty tuple to stream.
/// \details        Function overload that does nothing.
/// \param[in,out]  stream Output stream.
/// \return         True if the stream has no failure, false otherwise.
inline bool DataHandler::print(std::ostream& stream, const std::tuple<>&)
{
    return print(stream);
}

// Print a single variable to stream
/// \brief          Print a single variable to stream.
/// \details        Prints the next variable to the stream.
/// \tparam         Type (Variable type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      variable Variable.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type>
inline bool DataHandler::print(std::ostream& stream, const Type& variable)
{
    return (stream<<variable).good();
}

// Print several variables to stream
/// \brief          Print several variables to stream.
/// \details        Prints the next variables to the stream.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Output stream.
/// \param[in]      variable Variable.
/// \param[in]      variables Variables.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type, typename... Types> 
inline bool DataHandler::print(std::ostream& stream, const Type& variable, const Types&... variables)
{
    print(stream, variable);
    return print((bytesize(variable) > 0) ? (stream<<stream.fill()) : (stream), variables...);
}

// Print an array to stream
/// \brief          Print an array to stream.
/// \details        Prints each element of the passed array to the stream.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  stream Output stream.
/// \param[in]      container Array.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type, std::size_t Size> 
inline bool DataHandler::print(std::ostream& stream, const std::array<Type, Size>& container)
{
    if (bytesize(container) > 0) {
        print(stream, container[0]);
        for (unsigned int i = 1+(Size == 0); i < Size+(Size == 0); ++i) {
            print(stream<<stream.fill(), container[i]);
        }
    }
    return stream.good();
}

// Print a tuple to stream
/// \brief          Print a tuple to stream.
/// \details        Prints each element of the passed tuple to the stream.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Output stream.
/// \param[in]      container Tuple.
/// \return         True if the stream has no failure, false otherwise.
template <unsigned int Current, typename... Types> 
inline bool DataHandler::print(std::ostream& stream, const std::tuple<Types...>& container)
{
    print(stream, std::get<Current>(container));
    return ((Current+1) < sizeof...(Types)) ? (print<((Current+1) < sizeof...(Types))*(Current+1)>((bytesize(std::get<Current>(container)) > 0) ? (stream<<stream.fill()) : (stream), container)) : (print(stream));
}

// Print an empty range to stream
/// \brief          Print an empty range to stream.
/// \details        Function overload that does nothing.
/// \param[in,out]  stream Output stream.
/// \return         True if the stream has no failure, false otherwise.
inline bool DataHandler::rprint(std::ostream& stream, const std::nullptr_t, const std::nullptr_t)
{
    return print(stream);
}

// Print a range between pointers to stream
/// \brief          Print a range between pointers to stream.
/// \details        Prints each element of the range between the passed pointers 
///                 to the stream.
/// \tparam         Type (Pointer type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type>
inline bool DataHandler::rprint(std::ostream& stream, const Type* const first, const Type* const last)
{
    if (rbytesize(first, last) > 0) {
        print(stream, *first);
        for (const Type* ptr = first+1; ptr < last; ++ptr) {
            print(stream<<stream.fill(), *ptr);
        }
    }
    return stream.good();
}

// Print a range between iterators to stream
/// \brief          Print a range between iterators to stream.
/// \details        Prints each element of the range between the passed 
///                 iterators to the stream.
/// \tparam         Type (Iterator type.)
/// \param[in,out]  stream Output stream.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type, class>
inline bool DataHandler::rprint(std::ostream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    if (rbytesize(first, last) > 0) {
        print(stream, *first);
        for (Type it = first+1; it != last; ++it) {
            print(stream<<stream.fill(), *it);
        }
    }
    return stream.good();
}
// -------------------------------------------------------------------------- //



// ---------------------------------- SCAN ---------------------------------- //
// Scan nothing from stream
/// \brief          Scan nothing from stream.
/// \details        Function overload that does nothing.
/// \param[in,out]  stream Input stream.
/// \return         True if the stream has no failure, false otherwise.
inline bool DataHandler::scan(std::istream& stream)
{
    return stream.good();
}

// Scan an empty tuple from stream
/// \brief          Scan an empty tuple from stream.
/// \details        Function overload that does nothing.
/// \param[in,out]  stream Input stream.
/// \return         True if the stream has no failure, false otherwise.
inline bool DataHandler::scan(std::istream& stream, std::tuple<>&)
{
    return scan(stream);
}

// Scan a single variable from stream
/// \brief          Scan a single variable from stream.
/// \details        Scans the next variable from the stream.
/// \tparam         Type (Variable type.)
/// \param[in,out]  stream Input stream.
/// \param[out]     variable Variable.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type>
inline bool DataHandler::scan(std::istream& stream, Type& variable)
{
    return (stream>>variable).good();
}

// Scan several variables from stream
/// \brief          Scan several variables from stream.
/// \details        Scans the next variables from the stream.
/// \tparam         Type (Variable type.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Input stream.
/// \param[out]     variable Variable.
/// \param[out]     variables Variables.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type, typename... Types> 
inline bool DataHandler::scan(std::istream& stream, Type& variable, Types&... variables)
{
    scan(stream, variable);
    return scan(stream, variables...);
}

// Scan an array from stream
/// \brief          Scan an array from stream.
/// \details        Scans each element of the passed array from the stream.
/// \tparam         Type (Array type.)
/// \tparam         Size (Array size.)
/// \param[in,out]  stream Input stream.
/// \param[out]     container Array.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type, std::size_t Size> 
inline bool DataHandler::scan(std::istream& stream, std::array<Type, Size>& container)
{
    for (unsigned int i = (Size == 0); i < Size+(Size == 0); ++i) {
        scan(stream, container[i]);
    }
    return stream.good();
}

// Scan a tuple from stream
/// \brief          Scan a tuple from stream.
/// \details        Scans each element of the passed tuple from the stream.
/// \tparam         Current (Current level of recursion.)
/// \tparam         Types (Variadic types.)
/// \param[in,out]  stream Input stream.
/// \param[out]     container Tuple.
/// \return         True if the stream has no failure, false otherwise.
template <unsigned int Current, typename... Types> 
inline bool DataHandler::scan(std::istream& stream, std::tuple<Types...>& container)
{
    scan(stream, std::get<Current>(container));
    return ((Current+1) < sizeof...(Types)) ? (scan<((Current+1) < sizeof...(Types))*(Current+1)>(stream, container)) : (scan(stream));
}

// Scan an empty range from stream
/// \brief          Scan an empty range from stream.
/// \details        Function overload that does nothing.
/// \param[in,out]  stream Input stream.
/// \return         True if the stream has no failure, false otherwise.
inline bool DataHandler::rscan(std::istream& stream, const std::nullptr_t, const std::nullptr_t)
{
    return scan(stream);
}

// Scan a range between pointers from stream
/// \brief          Scan a range between pointers from stream.
/// \details        Scans each element of the range between the passed pointers 
///                 from the stream.
/// \tparam         Type (Pointer type.)
/// \param[in,out]  stream Input stream.
/// \param[in]      first Pointer to the beginning.
/// \param[in]      last Pointer to the end.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type>
inline bool DataHandler::rscan(std::istream& stream, Type* const first, Type* const last)
{
    for (Type* ptr = first; ptr < last; ++ptr) {
        scan(stream, *ptr);
    }
    return stream.good();
}

// Scan a range between iterators from stream
/// \brief          Scan a range between iterators from stream.
/// \details        Scans each element of the range between the passed iterators 
///                 from the stream.
/// \tparam         Type (Iterator type.)
/// \param[in,out]  stream Input stream.
/// \param[in]      first Iterator to the beginning.
/// \param[in]      last Iterator to the end.
/// \return         True if the stream has no failure, false otherwise.
template <typename Type, class>
inline bool DataHandler::rscan(std::istream& stream, const Type& first, const Type& last, typename std::iterator_traits<Type>::iterator_category*)
{
    for (Type it = first; it != last; ++it) {
        scan(stream, *it);
    }
    return stream.good();
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Timer.
/// \return         0 if no error.
int DataHandler::example()
{
    // Initialize
    std::cout<<"BEGIN = DataHandler::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    const std::string name(std::tmpnam(nullptr));
    std::fstream stream;
    char buffer[std::numeric_limits<short unsigned int>::max()] = {static_cast<char>(0)};
    char* ptrbuf = &buffer[0];
    char c = 'v';
    int i = 42;
    double d = 42.42;
    std::array<char, 0> ac0;
    std::array<int, 4> ai4 = {{0, 1, 2, 3}};
    std::tuple<> t = std::tuple<>();
    std::tuple<char, int, char> tcic('4', 42, '2'); 
    std::tuple<int, char, char> ticc(42, '4', '2');
    std::vector<int> vi(std::begin(ai4), std::end(ai4));
    int* ptrvib = &*vi.begin();
    int* ptrvie = &*vi.end();
    std::vector<int>::iterator itvib = vi.begin();
    std::vector<int>::iterator itvie = vi.end();
    
    // Construction
    DataHandler handler;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"DataHandler() : "                             ; DataHandler(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"handler = DataHandler() : "                   ; handler = DataHandler(); std::cout<<std::endl;

    // Utilities
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Utilities : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"handler.array(c) : "                          <<handler.array(c)<<std::endl;
    std::cout<<std::setw(width)<<"handler.array(ac0) : "                        <<handler.array(ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.array(ai4) : "                        <<handler.array(ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.array(t) : "                          <<handler.array(t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.array(tcic) : "                       <<handler.array(tcic)<<std::endl;
    std::cout<<std::setw(width)<<"handler.tuple(c) : "                          <<handler.tuple(c)<<std::endl;
    std::cout<<std::setw(width)<<"handler.tuple(ac0) : "                        <<handler.tuple(ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.tuple(ai4) : "                        <<handler.tuple(ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.tuple(t) : "                          <<handler.tuple(t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.tuple(tcic) : "                       <<handler.tuple(tcic)<<std::endl;

    // Size
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Size : "                                                                            <<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size() : "                                                                  <<handler.size()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size<void>() : "                                                            <<handler.size<void>()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size<double>() : "                                                          <<handler.size<double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size<std::tuple<char, double> >() : "                                       <<handler.size<std::tuple<char, double> >()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size<char, double>() : "                                                    <<handler.size<char, double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(t) : "                                                                 <<handler.size(t)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(i) : "                                                                 <<handler.size(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(c, i) : "                                                              <<handler.size(c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(ac0) : "                                                               <<handler.size(ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(ai4) : "                                                               <<handler.size(ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(tcic) : "                                                              <<handler.size(tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.size(tcic, ticc) : "                                                        <<handler.size(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rsize() : "                                                                 <<handler.rsize()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rsize(ptrvib, ptrvie) : "                                                   <<handler.rsize(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rsize(itvib, itvie) : "                                                     <<handler.rsize(itvib, itvie)<<std::endl;
    
    // Bytesize
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Bytesize : "                                                                        <<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize() : "                                                              <<handler.bytesize()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize<void>() : "                                                        <<handler.bytesize<void>()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize<double>() : "                                                      <<handler.bytesize<double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize<std::tuple<char, double> >() : "                                   <<handler.bytesize<std::tuple<char, double> >()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize<char, double>() : "                                                <<handler.bytesize<char, double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(t) : "                                                             <<handler.bytesize(t)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(i) : "                                                             <<handler.bytesize(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(c, i) : "                                                          <<handler.bytesize(c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(ac0) : "                                                           <<handler.bytesize(ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(ai4) : "                                                           <<handler.bytesize(ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(tcic) : "                                                          <<handler.bytesize(tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.bytesize(tcic, ticc) : "                                                    <<handler.bytesize(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rbytesize() : "                                                             <<handler.rbytesize()<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rbytesize(ptrvib, ptrvie) : "                                               <<handler.rbytesize(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rbytesize(itvib, itvie) : "                                                 <<handler.rbytesize(itvib, itvie)<<std::endl;

    // Byteswap
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Byteswap : "                                                                        <<std::endl;
    std::cout<<std::setw(width*2)<<"i / handler.byteswap<false>(i) / i : "                                              <<handler.hexify(i); std::cout<<" / "<<handler.byteswap<false>(i)<<" / "; std::cout<<handler.hexify(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"ptrvib ptrvie / handler.rbyteswap<false>(ptrvib, ptrvie) / ptrvib ptrvie : "        <<handler.rhexify(ptrvib, ptrvie); std::cout<<" / "<<handler.rbyteswap<false>(ptrvib, ptrvie)<<" / "; std::cout<<handler.rhexify(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<" / handler.byteswap() / : "                                                         <<handler.hexify(); std::cout<<" / "<<handler.byteswap()<<" / "; std::cout<<handler.hexify()<<std::endl;
    std::cout<<std::setw(width*2)<<"t / handler.byteswap(t) / t : "                                                     <<handler.hexify(t); std::cout<<" / "<<handler.byteswap(t)<<" / "; std::cout<<handler.hexify(t)<<std::endl;
    std::cout<<std::setw(width*2)<<"i / handler.byteswap(i) / i : "                                                     <<handler.hexify(i); std::cout<<" / "<<handler.byteswap(i)<<" / "; std::cout<<handler.hexify(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"c i / handler.byteswap(c, i) / c i : "                                              <<handler.hexify(c, i); std::cout<<" / "<<handler.byteswap(c, i)<<" / "; std::cout<<handler.hexify(c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"ac0 / handler.byteswap(ac0) / ac0 : "                                               <<handler.hexify(ac0); std::cout<<" / "<<handler.byteswap(ac0)<<" / "; std::cout<<handler.hexify(ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"ai4 / handler.byteswap(ai4) / ai4 : "                                               <<handler.hexify(ai4); std::cout<<" / "<<handler.byteswap(ai4)<<" / "; std::cout<<handler.hexify(ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"tcic / handler.byteswap(tcic) / tcic : "                                            <<handler.hexify(tcic); std::cout<<" / "<<handler.byteswap(tcic)<<" / "; std::cout<<handler.hexify(tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"tcic ticc / handler.byteswap(tcic, ticc) / tcic ticc : "                            <<handler.hexify(tcic, ticc); std::cout<<" / "<<handler.byteswap(tcic, ticc)<<" / "; std::cout<<handler.hexify(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<" / handler.rbyteswap() / : "                                                        <<handler.rhexify(); std::cout<<" / "<<handler.rbyteswap()<<" / "; std::cout<<handler.rhexify()<<std::endl;
    std::cout<<std::setw(width*2)<<"ptrvib ptrvie / handler.rbyteswap(ptrvib, ptrvie) / ptrvib ptrvie : "               <<handler.rhexify(ptrvib, ptrvie); std::cout<<" / "<<handler.rbyteswap(ptrvib, ptrvie)<<" / "; std::cout<<handler.rhexify(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"itvib itvie / handler.rbyteswap(itvib, itvie) / itvib itvie : "                     <<handler.rhexify(itvib, itvie); std::cout<<" / "<<handler.rbyteswap(itvib, itvie)<<" / "; std::cout<<handler.rhexify(itvib, itvie)<<std::endl;
 
    // Write stream
    stream.open(name, std::ios::out|std::ios::trunc|std::ios::binary);
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Write stream : "                                                                    <<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream) : "                                                           <<handler.write(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, t) : "                                                        <<handler.write(stream, t)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, i) : "                                                        <<handler.write(stream, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, c, i) : "                                                     <<handler.write(stream, c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, ac0) : "                                                      <<handler.write(stream, ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, ai4) : "                                                      <<handler.write(stream, ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, tcic) : "                                                     <<handler.write(stream, tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(stream, tcic, ticc) : "                                               <<handler.write(stream, tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite(stream) : "                                                          <<handler.rwrite(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite(stream, ptrvib, ptrvie) : "                                          <<handler.rwrite(stream, ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite(stream, itvib, itvie) : "                                            <<handler.rwrite(stream, itvib, itvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream) : "                                                     <<handler.write<true>(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, t) : "                                                  <<handler.write<true>(stream, t)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, i) : "                                                  <<handler.write<true>(stream, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, c, i) : "                                               <<handler.write<true>(stream, c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, ac0) : "                                                <<handler.write<true>(stream, ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, ai4) : "                                                <<handler.write<true>(stream, ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, tcic) : "                                               <<handler.write<true>(stream, tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(stream, tcic, ticc) : "                                         <<handler.write<true>(stream, tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite<true>(stream) : "                                                    <<handler.rwrite<true>(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite<true>(stream, ptrvib, ptrvie) : "                                    <<handler.rwrite<true>(stream, ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite<true>(stream, itvib, itvie) : "                                      <<handler.rwrite<true>(stream, itvib, itvie)<<std::endl;
    stream.close();
    
    // Read stream
    stream.open(name, std::ios::in|std::ios::binary);
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Read stream : "                                                                     <<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream) : "                                                            <<handler.read(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, t) : "                                                         <<handler.read(stream, t)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, i) : "                                                         <<handler.read(stream, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, c, i) : "                                                      <<handler.read(stream, c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, ac0) : "                                                       <<handler.read(stream, ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, ai4) : "                                                       <<handler.read(stream, ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, tcic) : "                                                      <<handler.read(stream, tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(stream, tcic, ticc) : "                                                <<handler.read(stream, tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread(stream) : "                                                           <<handler.rread(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread(stream, ptrvib, ptrvie) : "                                           <<handler.rread(stream, ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread(stream, itvib, itvie) : "                                             <<handler.rread(stream, itvib, itvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream) : "                                                      <<handler.read<true>(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, t) : "                                                   <<handler.read<true>(stream, t)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, i) : "                                                   <<handler.read<true>(stream, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, c, i) : "                                                <<handler.read<true>(stream, c, i)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, ac0) : "                                                 <<handler.read<true>(stream, ac0)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, ai4) : "                                                 <<handler.read<true>(stream, ai4)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, tcic) : "                                                <<handler.read<true>(stream, tcic)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(stream, tcic, ticc) : "                                          <<handler.read<true>(stream, tcic, ticc)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread<true>(stream) : "                                                     <<handler.rread<true>(stream)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread<true>(stream, ptrvib, ptrvie) : "                                     <<handler.rread<true>(stream, ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread<true>(stream, itvib, itvie) : "                                       <<handler.rread<true>(stream, itvib, itvie)<<std::endl;
    stream.close();
    
    // Write buffer
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Write buffer : "                                                                    <<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf) : "                                                           ; handler.write(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, t) : "                                                        ; handler.write(ptrbuf, t); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, i) : "                                                        ; handler.write(ptrbuf, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, c, i) : "                                                     ; handler.write(ptrbuf, c, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, ac0) : "                                                      ; handler.write(ptrbuf, ac0); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, ai4) : "                                                      ; handler.write(ptrbuf, ai4); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, tcic) : "                                                     ; handler.write(ptrbuf, tcic); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write(ptrbuf, tcic, ticc) : "                                               ; handler.write(ptrbuf, tcic, ticc); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite(ptrbuf) : "                                                          ; handler.rwrite(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite(ptrbuf, ptrvib, ptrvie) : "                                          ; handler.rwrite(ptrbuf, ptrvib, ptrvie); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite(ptrbuf, itvib, itvie) : "                                            ; handler.rwrite(ptrbuf, itvib, itvie); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf) : "                                                     ; handler.write<true>(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, t) : "                                                  ; handler.write<true>(ptrbuf, t); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, i) : "                                                  ; handler.write<true>(ptrbuf, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, c, i) : "                                               ; handler.write<true>(ptrbuf, c, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, ac0) : "                                                ; handler.write<true>(ptrbuf, ac0); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, ai4) : "                                                ; handler.write<true>(ptrbuf, ai4); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, tcic) : "                                               ; handler.write<true>(ptrbuf, tcic); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.write<true>(ptrbuf, tcic, ticc) : "                                         ; handler.write<true>(ptrbuf, tcic, ticc); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite<true>(ptrbuf) : "                                                    ; handler.rwrite<true>(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite<true>(ptrbuf, ptrvib, ptrvie) : "                                    ; handler.rwrite<true>(ptrbuf, ptrvib, ptrvie); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rwrite<true>(ptrbuf, itvib, itvie) : "                                      ; handler.rwrite<true>(ptrbuf, itvib, itvie); std::cout<<std::endl;
    
    // Read buffer
    ptrbuf = &buffer[0];
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Read buffer : "                                                                     <<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf) : "                                                            ; handler.read(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, t) : "                                                         ; handler.read(ptrbuf, t); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, i) : "                                                         ; handler.read(ptrbuf, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, c, i) : "                                                      ; handler.read(ptrbuf, c, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, ac0) : "                                                       ; handler.read(ptrbuf, ac0); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, ai4) : "                                                       ; handler.read(ptrbuf, ai4); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, tcic) : "                                                      ; handler.read(ptrbuf, tcic); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read(ptrbuf, tcic, ticc) : "                                                ; handler.read(ptrbuf, tcic, ticc); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread(ptrbuf) : "                                                           ; handler.rread(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread(ptrbuf, ptrvib, ptrvie) : "                                           ; handler.rread(ptrbuf, ptrvib, ptrvie); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread(ptrbuf, itvib, itvie) : "                                             ; handler.rread(ptrbuf, itvib, itvie); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf) : "                                                      ; handler.read<true>(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, t) : "                                                   ; handler.read<true>(ptrbuf, t); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, i) : "                                                   ; handler.read<true>(ptrbuf, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, c, i) : "                                                ; handler.read<true>(ptrbuf, c, i); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, ac0) : "                                                 ; handler.read<true>(ptrbuf, ac0); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, ai4) : "                                                 ; handler.read<true>(ptrbuf, ai4); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, tcic) : "                                                ; handler.read<true>(ptrbuf, tcic); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.read<true>(ptrbuf, tcic, ticc) : "                                          ; handler.read<true>(ptrbuf, tcic, ticc); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread<true>(ptrbuf) : "                                                     ; handler.rread<true>(ptrbuf); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread<true>(ptrbuf, ptrvib, ptrvie) : "                                     ; handler.rread<true>(ptrbuf, ptrvib, ptrvie); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"handler.rread<true>(ptrbuf, itvib, itvie) : "                                       ; handler.rread<true>(ptrbuf, itvib, itvie); std::cout<<std::endl;
    
    // Nullification
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Nullification : "                             <<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify() : "                         <<handler.nullify()<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(t) : "                        <<handler.nullify(t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(i) : "                        <<handler.nullify(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(c, i) : "                     <<handler.nullify(c, i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(ac0) : "                      <<handler.nullify(ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(ai4) : "                      <<handler.nullify(ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(tcic) : "                     <<handler.nullify(tcic)<<std::endl;
    std::cout<<std::setw(width)<<"handler.nullify(tcic, ticc) : "               <<handler.nullify(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rnullify() : "                        <<handler.rnullify()<<std::endl;
    std::cout<<std::setw(width)<<"handler.rnullify(ptrvib, ptrvie) : "          <<handler.rnullify(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rnullify(itvib, itvie) : "            <<handler.rnullify(itvib, itvie)<<std::endl;

    // Equalization
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Equalization : "                              <<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42) : "                      <<handler.equalize(42)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, t) : "                   <<handler.equalize(42, t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, i) : "                   <<handler.equalize(42, i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, c, i) : "                <<handler.equalize(42, c, i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, ac0) : "                 <<handler.equalize(42, ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, ai4) : "                 <<handler.equalize(42, ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, tcic) : "                <<handler.equalize(42, tcic)<<std::endl;
    std::cout<<std::setw(width)<<"handler.equalize(42, tcic, ticc) : "          <<handler.equalize(42, tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.requalize(42) : "                     <<handler.requalize(42)<<std::endl;
    std::cout<<std::setw(width)<<"handler.requalize(42, ptrvib, ptrvie) : "     <<handler.requalize(42, ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.requalize(42, itvib, itvie) : "       <<handler.requalize(42, itvib, itvie)<<std::endl;

    // Hexification
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Hexification : "                              <<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify() : "                          <<handler.hexify()<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(t) : "                         <<handler.hexify(t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(i) : "                         <<handler.hexify(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(c, i, d) : "                   <<handler.hexify(c, i, d)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(ac0) : "                       <<handler.hexify(ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(ai4) : "                       <<handler.hexify(ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(tcic) : "                      <<handler.hexify(tcic)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify(tcic, ticc) : "                <<handler.hexify(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rhexify() : "                         <<handler.rhexify()<<std::endl;
    std::cout<<std::setw(width)<<"handler.rhexify(ptrvib, ptrvie) : "           <<handler.rhexify(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rhexify(itvib, itvie) : "             <<handler.rhexify(itvib, itvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify<true>(tcic, ticc) : "          <<handler.hexify<true>(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.hexify<true, true>(tcic, ticc) :"     <<handler.hexify<true, true>(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rhexify(itvib, itvie, \", \") : "     <<handler.rhexify(itvib, itvie, ", ")<<std::endl;

    // Stringification
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stringification : "                           <<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify() : "                       <<handler.stringify()<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(t) : "                      <<handler.stringify(t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(i) : "                      <<handler.stringify(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(c, i, d) : "                <<handler.stringify(c, i, d)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(ac0) : "                    <<handler.stringify(ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(ai4) : "                    <<handler.stringify(ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(tcic) : "                   <<handler.stringify(tcic)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify(tcic, ticc) : "             <<handler.stringify(tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rstringify() : "                      <<handler.rstringify()<<std::endl;
    std::cout<<std::setw(width)<<"handler.rstringify(ptrvib, ptrvie) : "        <<handler.rstringify(ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rstringify(itvib, itvie) : "          <<handler.rstringify(itvib, itvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<16>(c) : "                  <<handler.stringify<16>(c)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<10, '0'>(c) : "             <<handler.stringify<10, '0'>(c)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<16, '0'>(c) : "             <<handler.stringify<16, '0'>(c)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<10, '0'>(i) : "             <<handler.stringify<10, '0'>(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<16, '0'>(i) : "             <<handler.stringify<16, '0'>(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<2>(i) : "                   <<handler.stringify<2>(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<10, '0'>(i) : "             <<handler.stringify<10, '0'>(i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.stringify<10, '0'>(-i) : "            <<handler.stringify<10, '0'>(-i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rstringify(itvib, itvie, \", \") :"   <<handler.rstringify(itvib, itvie, ", ")<<std::endl;

    // Print
    stream.open(name, std::ios::out|std::ios::trunc);
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Print : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream) : "                     <<handler.print(stream)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, t) : "                  <<handler.print(stream, t)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, i) : "                  <<handler.print(stream, i)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, c, i, d) : "            <<handler.print(stream, c, i, d)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, ac0) : "                <<handler.print(stream, ac0)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, ai4) : "                <<handler.print(stream, ai4)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, tcic) : "               <<handler.print(stream, tcic)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.print(stream, tcic, ticc) : "         <<handler.print(stream, tcic, ticc)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.rprint(stream) : "                    <<handler.rprint(stream)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.rprint(stream, ptrvib, ptrvie) :"     <<handler.rprint(stream, ptrvib, ptrvie)<<std::endl; stream<<std::endl;
    std::cout<<std::setw(width)<<"handler.rprint(stream, itvib, itvie) : "      <<handler.rprint(stream, itvib, itvie)<<std::endl; stream<<std::endl;
    stream.close();
    
    // Scan
    stream.open(name, std::ios::in);
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Scan : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream) : "                      <<handler.scan(stream)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, t) : "                   <<handler.scan(stream, t)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, i) : "                   <<handler.scan(stream, i)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, c, i, d) : "             <<handler.scan(stream, c, i, d)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, ac0) : "                 <<handler.scan(stream, ac0)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, ai4) : "                 <<handler.scan(stream, ai4)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, tcic) : "                <<handler.scan(stream, tcic)<<std::endl;
    std::cout<<std::setw(width)<<"handler.scan(stream, tcic, ticc) : "          <<handler.scan(stream, tcic, ticc)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rscan(stream) : "                     <<handler.rscan(stream)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rscan(stream, ptrvib, ptrvie) : "     <<handler.rscan(stream, ptrvib, ptrvie)<<std::endl;
    std::cout<<std::setw(width)<<"handler.rscan(stream, itvib, itvie) : "       <<handler.rscan(stream, itvib, itvie)<<std::endl;
    stream.close();
    
    // Finalize
    std::remove(name.data());
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = DataHandler::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // DATAHANDLER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
