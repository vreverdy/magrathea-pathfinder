/* ******************************* VECTORIZER ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Vectorizer
// DESCRIPTION :    Helper base class for generic vectorization
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           vectorizer.h
/// \brief          Helper base class for generic vectorization
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef VECTORIZER_H_INCLUDED
#define VECTORIZER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <stdexcept>
#include <functional>
#include <initializer_list>
#include <vector>
#include <array>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Helper base class for generic vectorization
/// \brief          Helper base class for generic vectorization.
/// \details        Provides global functions for vectorization implementation.
///                 Vectorizer helpers (like StaticVectorizer or 
///                 DynamicVectorizer) are derived from this class and have to
///                 implement the following functions required by CRTP : <ul>
///                 <li><tt>operator[]</tt>
///                 <li><tt>resize()</tt>
///                 <li><tt>size()</tt>
///                 <li><tt>constant()</tt>
///                 <li><tt>boolean()</tt>
///                 <li><tt>parameters()</tt>
///                 <li><tt>type()</tt>
///                 </ul>
class Vectorizer
{
    // Protected lifecycle
    /// \name           Protected lifecycle
    //@{
    protected: 
        inline ~Vectorizer();
    //@}
    
    // Vectorization
    /// \name           Vectorization
    //@{
    public:
        inline Vectorizer& operator[](const unsigned int i);
        inline const Vectorizer& operator[](const unsigned int i) const;
        inline Vectorizer& resize(const unsigned int n);
        inline unsigned int size() const;
        inline bool constant() const;
        inline bool boolean() const;
        inline std::array<unsigned int, 0> parameters() const;
        inline void type() const;
    //@}
    
    // Check
    /// \name           Check
    //@{
    public:
        template <bool Exception = true, class FirstType, class SecondType, class = typename std::enable_if<(!std::is_base_of<Vectorizer, FirstType>::value) || (!std::is_base_of<Vectorizer, SecondType>::value)>::type> static constexpr bool check(const FirstType&, const SecondType&);
        template <bool Exception = true, class FirstType, class SecondType, class = typename std::enable_if<(std::is_base_of<Vectorizer, FirstType>::value) && (std::is_base_of<Vectorizer, SecondType>::value)>::type, class = typename std::enable_if<(FirstType::constant()) && (SecondType::constant())>::type> static constexpr bool check(const FirstType&, const SecondType&);
        template <bool Exception = true, class FirstType, class SecondType, class = typename std::enable_if<(std::is_base_of<Vectorizer, FirstType>::value) && (std::is_base_of<Vectorizer, SecondType>::value)>::type, class = typename std::enable_if<(!FirstType::constant()) || (!SecondType::constant())>::type, class = void> static inline bool check(const FirstType& first, const SecondType& second);
    //@}
    
    // Getters
    /// \name           Getters
    //@{
    public:
        template <typename IntegralType, IntegralType Value, typename DummyType = unsigned int, class = typename std::enable_if<(std::is_integral<IntegralType>::value) && (std::is_convertible<DummyType, unsigned int>::value)>::type> static constexpr IntegralType get(const std::integral_constant<IntegralType, Value>, const DummyType = DummyType());
        template <typename DataType, typename DummyType = unsigned int, class = typename std::enable_if<(!std::is_base_of<Vectorizer, DataType>::value) && (std::is_convertible<DummyType, unsigned int>::value)>::type> static inline DataType& get(DataType& source, const DummyType = DummyType());
        template <typename DataType, typename DummyType = unsigned int, class = typename std::enable_if<(!std::is_base_of<Vectorizer, DataType>::value) && (std::is_convertible<DummyType, unsigned int>::value)>::type> static inline const DataType& get(const DataType& source, const DummyType = DummyType());
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, typename DataType = typename std::remove_reference<decltype(VectorizedType::type())>::type> static inline DataType& get(VectorizedType& source, const unsigned int i);
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, typename DataType = typename std::remove_reference<decltype(VectorizedType::type())>::type> static inline const DataType& get(const VectorizedType& source, const unsigned int i);
    //@}
    
    // Setters
    /// \name           Setters
    //@{
    public:
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type> static inline VectorizedType& set(VectorizedType& destination);
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, class GenericType> static inline VectorizedType& set(VectorizedType& destination, const GenericType& source);
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, class GenericType, class First, class Second, class... Others, class = typename std::enable_if<(sizeof...(Others) != 0) && (sizeof...(Others) != 1)>::type> static inline VectorizedType& set(VectorizedType& destination, const GenericType& source, const First& first, const Second& second, const Others&... others);
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, class GenericType, typename SizeType = std::true_type, class = typename std::enable_if<(std::is_same<SizeType, std::true_type>::value) || (std::is_convertible<SizeType, unsigned int>::value)>::type> static inline VectorizedType& set(VectorizedType& destination, const GenericType& source, const unsigned int pos, const SizeType num = SizeType());
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, typename DataType = typename std::remove_reference<decltype(VectorizedType::type())>::type, typename SizeType = std::true_type, class = typename std::enable_if<(std::is_convertible<DataType, typename std::remove_reference<decltype(std::declval<VectorizedType>()[0])>::type>::value) && ((std::is_same<SizeType, std::true_type>::value) || (std::is_convertible<SizeType, unsigned int>::value))>::type> static inline VectorizedType& set(VectorizedType& destination, const std::initializer_list<DataType>& source, const unsigned int pos = 0, const SizeType num = SizeType());
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, typename IteratorType, typename SizeType = std::true_type, class = typename std::enable_if<(std::is_same<SizeType, std::true_type>::value) || (std::is_convertible<SizeType, unsigned int>::value)>::type> static inline VectorizedType& set(VectorizedType& destination, const IteratorType& first, const IteratorType& last, const unsigned int pos = 0, const SizeType num = SizeType(), typename std::iterator_traits<IteratorType>::iterator_category* = nullptr);
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, class GenericType, class MaskType, typename DummyType = unsigned int, class = typename std::enable_if<(std::is_base_of<Vectorizer, MaskType>::value) && (std::is_convertible<DummyType, unsigned int>::value)>::type> static inline VectorizedType& set(VectorizedType& destination, const GenericType& source, const MaskType& mask, const DummyType = DummyType());
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, typename DataType = typename std::remove_reference<decltype(VectorizedType::type())>::type, class MaskType, typename DummyType = unsigned int, class = typename std::enable_if<(std::is_base_of<Vectorizer, MaskType>::value) && (std::is_convertible<DummyType, unsigned int>::value)>::type> static inline VectorizedType& set(VectorizedType& destination, const std::initializer_list<DataType>& source, const MaskType& mask, const DummyType = DummyType());
        template <class VectorizedType, class = typename std::enable_if<std::is_base_of<Vectorizer, VectorizedType>::value>::type, typename IteratorType, class MaskType, typename DummyType = unsigned int, class = typename std::enable_if<(std::is_base_of<Vectorizer, MaskType>::value) && (std::is_convertible<DummyType, unsigned int>::value)>::type> static inline VectorizedType& set(VectorizedType& destination, const IteratorType& first, const IteratorType& last, const MaskType& mask, const DummyType = DummyType(), typename std::iterator_traits<IteratorType>::iterator_category* = nullptr);
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
inline Vectorizer::~Vectorizer()
= default;
// -------------------------------------------------------------------------- //



// ----------------------------- VECTORIZATION ------------------------------ //
// Direct access to the element
/// \brief          Direct access to the element.
/// \details        Provides a direct access to the specified element. This 
///                 function is required by the vectorization mechanism. 
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline Vectorizer& Vectorizer::operator[](const unsigned int i)
{
    throw std::logic_error("ERROR = Vectorizer::operator[]() : this function should be overloaded by derived classes");
    (void)(i);
    return *this;
}

// Immutable direct access to the element
/// \brief          Immutable direct access to the element.
/// \details        Provides a constant direct access to the specified element.
///                 This function is required by the vectorization mechanism.
/// \param[in]      i Index of the element.
/// \return         Const reference to the element.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline const Vectorizer& Vectorizer::operator[](const unsigned int i) const
{
    throw std::logic_error("ERROR = Vectorizer::operator[]() : this function should be overloaded by derived classes");
    (void)(i);
    return *this;
}

// Resize the container
/// \brief          Resize the container.
/// \details        Resizes the container and returns a reference to it. This 
///                 function is required by the vectorization mechanism. 
/// \param[in]      n New size.
/// \return         Self reference.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline Vectorizer& Vectorizer::resize(const unsigned int n)
{
    throw std::logic_error("ERROR = Vectorizer::resize() : this function should be overloaded by derived classes");
    (void)(n);
    return *this;
}

// Get the size of the container
/// \brief          Get the size of the container.
/// \details        Returns the current number of elements. This function is 
///                 required by the vectorization mechanism. 
/// \return         Copy of the size.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline unsigned int Vectorizer::size() const
{
    throw std::logic_error("ERROR = Vectorizer::size() : this function should be overloaded by derived classes");
    return 0;
}

// Get whether the container has a constant size
/// \brief          Get whether the container has a constant size.
/// \details        Returns true if the container has a constant size, false
///                 otherwise. This function is required by the vectorization 
///                 mechanism. 
/// \return         Copy of true if the container has a constant size.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline bool Vectorizer::constant() const
{
    throw std::logic_error("ERROR = Vectorizer::constant() : this function should be overloaded by derived classes");
    return true;
}

// Get whether the container has a boolean type
/// \brief          Get whether the container has a boolean type.
/// \details        Returns true if the container has a boolean type size, false
///                 otherwise. This function is required by the vectorization 
///                 mechanism. 
/// \return         Copy of true if the container has a boolean type.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline bool Vectorizer::boolean() const
{
    throw std::logic_error("ERROR = Vectorizer::boolean() : this function should be overloaded by derived classes");
    return false;
}

// Get the template parameters
/// \brief          Get the template parameters.
/// \details        Returns an array containing the template parameters. This 
///                 function is required by the vectorization mechanism.
/// \return         Copy of an array of parameters.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline std::array<unsigned int, 0> Vectorizer::parameters() const
{
    throw std::logic_error("ERROR = Vectorizer::parameters() : this function should be overloaded by derived classes");
    return std::array<unsigned int, 0>();
}

// Get the data type
/// \brief          Get the data type.
/// \details        Returns a copy of the default value of the data type.
/// \return         Copy of the default value of the data type.
/// \exception      std::logic_error This function should be overloaded by 
///                 derived classes.
inline void Vectorizer::type() const
{
    throw std::logic_error("ERROR = Vectorizer::type() : this function should be overloaded by derived classes");
}
// -------------------------------------------------------------------------- //



// ---------------------------------- CHECK --------------------------------- //
// Check compatibility with at least one non-vectorized type
/// \brief          Check compatibility with at least one non-vectorized type.
/// \details        Checks whether the two provided vectorizer have compatible
///                 properties.
/// \tparam         Exception Throw exception or assertion on error.
/// \tparam         FirstType (First type.)
/// \tparam         SecondType (Second type.)
/// \return         True whether the two types are compatible, false otherwise.
template <bool Exception, class FirstType, class SecondType, class> 
constexpr bool Vectorizer::check(const FirstType&, const SecondType&)
{
    return true;
}

// Check compatibility between two static vectorized types
/// \brief          Check compatibility between two static vectorized types.
/// \details        Checks whether the two provided vectorizer have compatible
///                 properties.
/// \tparam         Exception Throw exception or assertion on error.
/// \tparam         FirstType (First type.)
/// \tparam         SecondType (Second type.)
/// \return         True whether the two types are compatible, false otherwise.
template <bool Exception, class FirstType, class SecondType, class, class> 
constexpr bool Vectorizer::check(const FirstType&, const SecondType&)
{
    static_assert((!Exception) || (FirstType::size() == SecondType::size()), "ERROR = Vectorizer::check() : vectorizer sizes are not equal");
    return (FirstType::size() == SecondType::size());
}

// Check compatibility between two vectorized types
/// \brief          Check compatibility between two vectorized types
/// \details        Checks whether the two provided vectorizer have compatible
///                 properties.
/// \tparam         Exception Throw exception or assertion on error.
/// \tparam         FirstType (First type.)
/// \tparam         SecondType (Second type.)
/// \param[in]      first First argument.
/// \param[in]      second Second argument.
/// \return         True whether the two types are compatible, false otherwise.
/// \exception      std::length_error Vectorizer sizes are not equal.
template <bool Exception, class FirstType, class SecondType, class, class, class>
inline bool Vectorizer::check(const FirstType& first, const SecondType& second)
{
    if ((Exception) && (first.size() != second.size())) {
        throw std::length_error("ERROR = Vectorizer::check() : vectorizer sizes are not equal");
    }
    return (first.size() == second.size());
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// Integral constant getter
/// \brief          Integral constant getter.
/// \details        Returns the value of the provided integral constant.
/// \tparam         IntegralType (Integral type.)
/// \tparam         DummyType (Dummy parameter type.)
/// \return         Value of the integral constant.
template <typename IntegralType, IntegralType Value, typename DummyType, class> 
constexpr IntegralType Vectorizer::get(const std::integral_constant<IntegralType, Value>, const DummyType)
{
    return std::integral_constant<IntegralType, Value>::value;
}

// Non-vectorized getter
/// \brief          Non-vectorized getter.
/// \details        Returns a reference to the provided non-vectorized source.
/// \tparam         DataType (Data type.)
/// \tparam         DummyType (Dummy parameter type.)
/// \param[in]      source Argument to get.
/// \return         Reference to the argument.
template <typename DataType, typename DummyType, class> 
inline DataType& Vectorizer::get(DataType& source, const DummyType)
{
    return source;
}

// Immutable non-vectorized getter
/// \brief          Immutable non-vectorized getter.
/// \details        Returns a constant reference to the provided non-vectorized 
///                 source.
/// \tparam         DataType (Data type.)
/// \tparam         DummyType (Dummy parameter type.)
/// \param[in]      source Argument to get.
/// \return         Const reference to the argument.
template <typename DataType, typename DummyType, class> 
inline const DataType& Vectorizer::get(const DataType& source, const DummyType)
{
    return source;
}

// Vector element getter
/// \brief          Vector element getter.
/// \details        Returns a reference to the i-th element of the provided 
///                 vectorized source.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         DataType (Data type.)
/// \param[in]      source Source container.
/// \param[in]      i Index of the element.
/// \return         Reference to the specified element.
template <class VectorizedType, class, typename DataType> 
inline DataType& Vectorizer::get(VectorizedType& source, const unsigned int i)
{
    return source[i];
}

// Immutable vector element getter
/// \brief          Immutable vector element getter.
/// \details        Returns a constant reference to the i-th element of the 
///                 provided vectorized source.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         DataType (Data type.)
/// \param[in]      source Source container.
/// \param[in]      i Index of the element.
/// \return         Const reference to the specified element.
template <class VectorizedType, class, typename DataType> 
inline const DataType& Vectorizer::get(const VectorizedType& source, const unsigned int i)
{
    return source[i];
}
// -------------------------------------------------------------------------- //



// --------------------------------- SETTERS -------------------------------- //
// Empty setter
/// \brief          Empty setter.
/// \details        Does nothing.
/// \tparam         VectorizedType (Vectorized type.)
/// \param[in,out]  destination Destination of the copy.
/// \return         Reference to the destination. 
template <class VectorizedType, class> 
inline VectorizedType& Vectorizer::set(VectorizedType& destination)
{
    return destination;
}

// Copy setter
/// \brief          Copy setter.
/// \details        Copies the whole contents of the source to the destination.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         GenericType (Generic type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      source Source of the copy.
/// \return         Reference to the destination.
template <class VectorizedType, class, class GenericType> 
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const GenericType& source)
{
    const unsigned int lim = check(destination, source)*destination.size();
    for (unsigned int i = 0; i < lim; ++i) {
        get(destination, i) = get(source, i);
    }
    return destination;
}

// Variadic setter
/// \brief          Variadic setter.
/// \details        Calls recursively the setters for a long list of arguments.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         GenericType (Generic type.)
/// \tparam         First (First type.)
/// \tparam         Second (Second type.)
/// \tparam         Others (Other types.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      source Source of the copy.
/// \param[in]      first First extra argument.
/// \param[in]      second Second extra argument.
/// \param[in]      others Other arguments.
/// \return         Reference to the destination.
template <class VectorizedType, class, class GenericType, class First, class Second, class... Others, class> 
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const GenericType& source, const First& first, const Second& second, const Others&... others)
{
    return set(set(destination, source, first, second), source, others...);
}

// Partial setter
/// \brief          Partial setter.
/// \details        Copies the contents of the source to a part of the 
///                 destination.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         GenericType (Generic type.)
/// \tparam         SizeType (Size type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      source Source of the copy.
/// \param[in]      pos Starting position of the copy.
/// \param[in]      num Number of elements to copy.
/// \return         Reference to the destination.
template <class VectorizedType, class, class GenericType, typename SizeType, class>
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const GenericType& source, const unsigned int pos, const SizeType num)
{
    const unsigned int lim = check(destination, source)*(((std::is_same<SizeType, std::true_type>::value) || (pos+get(num) > destination.size())) ? (destination.size()) : (pos+get(num)));
    for (unsigned int i = pos; i < lim; ++i) {
        get(destination, i) = get(source, i);
    }
    return destination;
}

// Partial list setter
/// \brief          Partial list setter.
/// \details        Copies the contents of the source to a part of the 
///                 destination. The first element of the list is copied at the
///                 provided position, and the next elements are copied after 
///                 it. If the list is too small, empty values are added to its
///                 end.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         DataType (Data type.)
/// \tparam         SizeType (Size type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      source Source of the copy.
/// \param[in]      pos Starting position of the copy.
/// \param[in]      num Number of elements to copy.
/// \return         Reference to the destination.
template <class VectorizedType, class, typename DataType, typename SizeType, class>
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const std::initializer_list<DataType>& source, const unsigned int pos, const SizeType num)
{
    const unsigned int lim = ((std::is_same<SizeType, std::true_type>::value) || (pos+get(num) > destination.size())) ? (destination.size()) : (pos+get(num));
    const DataType* last = source.end();
    const DataType* it = source.begin();
    for (unsigned int i = pos; i < lim; ++i) {
        if (it < last) {
            get(destination, i) = *(it++);
        } else {
            get(destination, i) = DataType();
        }
    }
    return destination;
}

// Partial range setter
/// \brief          Partial range setter.
/// \details        Copies the values from the range to a part of the 
///                 destination. The first element of the range is copied at the
///                 provided position, and the next elements are copied after 
///                 it. The copy stops as soon as the end of the range is 
///                 encountered or if the number of elements to copy is reached.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         IteratorType (Iterator or pointer type.)
/// \tparam         SizeType (Size type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      first Iterator to the beginning of the range.
/// \param[in]      last Iterator to the end of the range.
/// \param[in]      pos Starting position of the copy.
/// \param[in]      num Number of elements to copy.
/// \return         Reference to the destination.
template <class VectorizedType, class, typename IteratorType, typename SizeType, class>
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const IteratorType& first, const IteratorType& last, const unsigned int pos, const SizeType num, typename std::iterator_traits<IteratorType>::iterator_category*)
{
    const unsigned int mark = ((std::is_same<SizeType, std::true_type>::value) || (pos+get(num) > destination.size())) ? (destination.size()) : (pos+get(num));
    const unsigned int lim = pos+(static_cast<bool>(last > first)*(last-first)) < mark ? pos+(static_cast<bool>(last > first)*(last-first)) : mark;
    IteratorType it = first;
    for (unsigned int i = pos; i < lim; ++i) {
        get(destination, i) = *(it++);
    }
    return destination;
}

// Masked setter
/// \brief          Masked setter.
/// \details        Copies elements of the source to the destination using a 
///                 mask of boolean values : the values are copied only where 
///                 the mask is true.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         GenericType (Generic type.)
/// \tparam         MaskType (Mask Type.)
/// \tparam         DummyType (Dummy parameter type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      source Source of the copy.
/// \param[in]      mask Boolean mask.
/// \return         Reference to the destination.
template <class VectorizedType, class, class GenericType, class MaskType, typename DummyType, class> 
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const GenericType& source, const MaskType& mask, const DummyType)
{
    const unsigned int lim = check(destination, mask)*destination.size();
    for (unsigned int i = 0; i < lim; ++i) {
        if (static_cast<bool>(get(mask, i))) {
            get(destination, i) = get(source, i);
        }
    }
    return destination;
}

// Masked list setter
/// \brief          Masked list setter.
/// \details        Copies the contents of the source to the destination using a
///                 mask of boolean values : the values are copied only where 
///                 the mask is true. The iteration over values in the 
///                 destination and in the source list are independant : the 
///                 n-th element of the list is copied to the n-th true element
///                 of the destination. If the list is too small, empty values 
///                 are added to its end.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         DataType (Data type.)
/// \tparam         MaskType (Mask Type.)
/// \tparam         DummyType (Dummy parameter type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      source Source of the copy.
/// \param[in]      mask Boolean mask.
/// \return         Reference to the destination.
template <class VectorizedType, class, typename DataType, class MaskType, typename DummyType, class> 
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const std::initializer_list<DataType>& source, const MaskType& mask, const DummyType)
{
    const unsigned int lim = check(destination, mask)*destination.size();
    const DataType* last = source.end();
    const DataType* it = source.begin();
    for (unsigned int i = 0; i < lim; ++i) {
        if (static_cast<bool>(get(mask, i))) {
            if (it < last) {
                get(destination, i) = *(it++);
            } else {
                get(destination, i) = DataType();
            }
        }
    }
    return destination;
}

// Masked range setter
/// \brief          Masked range setter.
/// \details        Copies the values from the range to the destination using a
///                 mask of boolean values : the values are copied only where 
///                 the mask is true. The iteration over values in the 
///                 destination and in the range list are independant : the 
///                 n-th element of the range is copied to the n-th true element
///                 of the destination. The copy stops as soon as the end of the 
///                 range is encountered or if the number of elements to copy is
///                 reached.
/// \tparam         VectorizedType (Vectorized type.)
/// \tparam         Iterator (Iterator or pointer type.)
/// \tparam         MaskType (Mask Type.)
/// \tparam         DummyType (Dummy parameter type.)
/// \param[in,out]  destination Destination of the copy.
/// \param[in]      first Iterator to the beginning of the range.
/// \param[in]      last Iterator to the end of the range.
/// \param[in]      mask Boolean mask.
/// \return         Reference to the destination.
template <class VectorizedType, class, typename IteratorType, class MaskType, typename DummyType, class> 
inline VectorizedType& Vectorizer::set(VectorizedType& destination, const IteratorType& first, const IteratorType& last, const MaskType& mask, const DummyType, typename std::iterator_traits<IteratorType>::iterator_category*)
{
    const unsigned int lim = check(destination, mask)*destination.size();
    IteratorType it = first;
    for (unsigned int i = 0; (i < lim) && (it < last); ++i) {
        if (static_cast<bool>(get(mask, i))) {
            get(destination, i) = *(it++);
        }
    }
    return destination;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Vectorizer.
/// \return         0 if no error.
int Vectorizer::example()
{
    std::cout<<"BEGIN = Vectorizer::example()"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"ERROR = Vectorizer::example() : no example is provided for an abstract class"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"END = Vectorizer::example()"<<std::endl;
    return 1;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // VECTORIZER_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
