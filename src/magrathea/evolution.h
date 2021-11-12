/* ******************************** EVOLUTION ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          Evolution 
// DESCRIPTION :    Resizable container of steps dedicated to integration
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           evolution.h 
/// \brief          Resizable container of steps dedicated to integration
/// \date           2012-2013 
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef EVOLUTION_H_INCLUDED
#define EVOLUTION_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



//------------------------------- PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <utility>
#include <type_traits>
#include <initializer_list>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <vector>
#include <string>
// Include libs
// Include project
#include "step.h"
// Misc
namespace magrathea {
//--------------------------------------------------------------------------- //



//---------------------------------- CLASS ---------------------------------- //
// Resizable container of steps dedicated to integration
/// \brief          Resizable container of steps dedicated to integration. 
/// \details        An evolution container is an accumulator of integration 
///                 steps. It has the standard functions of containers plus
///                 additional functions to interpolate and sort the results.
/// \tparam         Type Step type.
/// \tparam         Container Container type.
template <class Type = Step<>, class Container = std::vector<Type> >
class Evolution final
{
    // Setup
    static_assert(std::is_same<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, Type>::value, "ERROR = Evolution : non matching element and container types");
    
    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        template <class OtherType> inline Evolution(const std::initializer_list<OtherType>& source);
        template <class OtherType, class OtherContainer> explicit inline Evolution(const Evolution<OtherType, OtherContainer>& source);
        template <class... Misc, class = typename std::enable_if<std::is_constructible<Container, Misc...>::value>::type> explicit inline Evolution(Misc&&... misc); 
    //@}
    
    // Operators
    /// \name           Operators
    //@{
    public:
        template <class OtherType> inline Evolution<Type, Container>& operator=(const std::initializer_list<OtherType>& source);
        template <class OtherType, class OtherContainer> inline Evolution<Type, Container>& operator=(const Evolution<OtherType, OtherContainer>& source);
        template <class Misc, class = typename std::enable_if<std::is_convertible<Container, typename std::remove_cv<typename std::remove_reference<Misc>::type>::type>::value>::type> inline Evolution<Type, Container>& operator=(Misc&& misc);
        inline bool operator==(const Evolution<Type, Container>& rhs) const;
        inline bool operator!=(const Evolution<Type, Container>& rhs) const;
        inline bool operator>(const Evolution<Type, Container>& rhs) const;
        inline bool operator<(const Evolution<Type, Container>& rhs) const;
        inline bool operator>=(const Evolution<Type, Container>& rhs) const;
        inline bool operator<=(const Evolution<Type, Container>& rhs) const;
        inline Type& operator[](const unsigned int i);
        inline const Type& operator[](const unsigned int i) const;
    //@}
    
    // Assignment 
    /// \name           Assignment 
    //@{
    public: 
        inline Evolution<Type, Container>& assign(const Evolution<Type, Container>& source);
        template <class OtherType> inline Evolution<Type, Container>& assign(const std::initializer_list<OtherType>& source);
        template <class OtherType, class OtherContainer> inline Evolution<Type, Container>& assign(const Evolution<OtherType, OtherContainer>& source);
        template <class... Misc, class = typename std::enable_if<std::is_constructible<Container, Misc...>::value>::type> inline Evolution<Type, Container>& assign(Misc&&... misc); 
    //@}

    // Management
    /// \name           Management
    //@{
    public:
        inline Evolution<Type, Container>& nullify();
        inline Evolution<Type, Container> copy() const;
        template <class OtherType = Type, class OtherContainer = std::vector<OtherType> > inline Evolution<OtherType, OtherContainer> cast() const;
    //@}

    // Access 
    /// \name           Access 
    //@{
    public: 
        inline Type& at(const unsigned int i);
        inline const Type& at(const unsigned int i) const;  
        inline Type& front(const unsigned int i = 0);
        inline const Type& front(const unsigned int i = 0) const;  
        inline Type& back(const unsigned int i = 0);
        inline const Type& back(const unsigned int i = 0) const;  
        inline Type& cycle(const int i);
        inline const Type& cycle(const int i) const; 
        inline Container& container();
        inline const Container& container() const;
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
        template <typename Iterator> inline unsigned int index(const Iterator& it) const;
    //@}

    // Capacity 
    /// \name           Capacity 
    //@{
    public: 
        inline bool empty() const;
        inline unsigned int size() const;
        inline unsigned int capacity() const;
        inline Evolution<Type, Container>& reserve(const unsigned int n);
        inline Evolution<Type, Container>& shrink();
        inline unsigned long long int space() const;
    //@}

    // Modifiers  
    /// \name           Modifiers  
    //@{
    public: 
        inline Evolution<Type, Container>& fullclear();
        inline Evolution<Type, Container>& clear();
        template <class... Misc> inline Evolution<Type, Container>& resize(const Misc&... misc);
        inline Evolution<Type, Container>& pop();
        template <class OtherType, class = typename std::enable_if<std::is_convertible<OtherType, Type>::value>::type> inline Evolution<Type, Container>& append(const OtherType& value);
    //@}

    // Stream 
    /// \name           Stream 
    //@{
    public: 
        template <class SelfType, class SelfContainer> friend std::ostream& operator<<(std::ostream& lhs, const Evolution<SelfType, SelfContainer>& rhs); 
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
        Container _container;                                                   ///< Internal container. 
    //@}
};
//--------------------------------------------------------------------------- //



//-------------------------------- LIFECYCLE -------------------------------- //
// Implicit initializer list constructor
/// \brief          Implicit initializer list constructor.
/// \details        Provides an implicit conversion from an initializer list.
/// \tparam         OtherType (Other step type.)
/// \param[in]      source Source of the copy.  
template <class Type, class Container>
template <class OtherType>
inline Evolution<Type, Container>::Evolution(const std::initializer_list<OtherType>& source)
: _container(source)
{
    ;
}

// Explicit conversion constructor
/// \brief          Explicit conversion constructor.
/// \details        Provides an explicit construction from another type of
///                 evolution container.
/// \tparam         OtherType (Other step type.)
/// \tparam         OtherContainer (Other container type.)
/// \param[in]      source Source of the copy. 
template <class Type, class Container> 
template <class OtherType, class OtherContainer>
inline Evolution<Type, Container>::Evolution(const Evolution<OtherType, OtherContainer>& source)
: _container()
{
    const unsigned int n = source.container().size();
    _container.resize(n);
    for (unsigned int i = 0; i < n; ++i) {
        _container[i] = source.container()[i];
    }
}

// Explicit generic constructor 
/// \brief          Explicit generic constructor. 
/// \details        Provides a generic interface to all constructors of the 
///                 container.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
template <class Type, class Container>
template <class... Misc, class> 
inline Evolution<Type, Container>::Evolution(Misc&&... misc)
: _container(std::forward<Misc>(misc)...)
{
    ;
}
//--------------------------------------------------------------------------- //



//-------------------------------- OPERATORS -------------------------------- //
// Initializer list assignment operator
/// \brief          Initializer list assignment operator.
/// \details        Provides an initializer list assignment.
/// \tparam         OtherType (Other step type.)
/// \param[in]      source Source of the copy.  
/// \return         Self reference.
template <class Type, class Container>
template <class OtherType> 
inline Evolution<Type, Container>& Evolution<Type, Container>::operator=(const std::initializer_list<OtherType>& source)
{
    _container = source;
    return *this;
}

// Conversion assignment operator  
/// \brief          Conversion assignment operator. 
/// \details        Provides an assignment from another type of evolution 
///                 container.
/// \tparam         OtherType (Other step type.)
/// \tparam         OtherContainer (Other container type.)
/// \param[in]      source Source of the copy. 
/// \return         Self reference. 
template <class Type, class Container>
template <class OtherType, class OtherContainer>
inline Evolution<Type, Container>& Evolution<Type, Container>::operator=(const Evolution<OtherType, OtherContainer>& source)
{
    const unsigned int n = source.container().size();
    _container.resize(n);
    for (unsigned int i = 0; i < n; ++i) {
        _container[i] = source.container()[i];
    }
    return *this;
}

// Generic assignment operator  
/// \brief          Generic assignment operator. 
/// \details        Provides a generic interface to all assignment operators of 
///                 the container
/// \tparam         Misc (Miscellaneous type.)
/// \param[in]      misc Miscellaneous argument.
/// \return         Self reference. 
template <class Type, class Container>
template <class Misc, class> 
inline Evolution<Type, Container>& Evolution<Type, Container>::operator=(Misc&& misc)
{
    _container = std::forward<Misc>(misc);
    return *this;
}

// Equal to
/// \brief          Equal to.
/// \details        Lexicographically compares the values in the container for
///                 equality. 
/// \param[in]      rhs Right-hand side.
/// \return         True if equal, false if not equal. 
template <class Type, class Container>
inline bool Evolution<Type, Container>::operator==(const Evolution<Type, Container>& rhs) const
{
    return (_container == rhs._container);
}

// Not equal to
/// \brief          Not equal to.
/// \details        Lexicographically compares the values in the container for
///                 difference. 
/// \param[in]      rhs Right-hand side.
/// \return         True if not equal, false if equal. 
template <class Type, class Container>
inline bool Evolution<Type, Container>::operator!=(const Evolution<Type, Container>& rhs) const
{
    return (_container != rhs._container);
}

// Greater than
/// \brief          Greater than.
/// \details        Lexicographically compares the values in the container for
///                 the greater than comparison operator. 
/// \param[in]      rhs Right-hand side.
/// \return         True if greater than, false otherwise.  
template <class Type, class Container>
inline bool Evolution<Type, Container>::operator>(const Evolution<Type, Container>& rhs) const
{
    return (_container > rhs._container);
}

// Less than
/// \brief          Less than.
/// \details        Lexicographically compares the values in the container for
///                 the less than comparison operator. 
/// \param[in]      rhs Right-hand side.
/// \return         True if less than, false otherwise. 
template <class Type, class Container>
inline bool Evolution<Type, Container>::operator<(const Evolution<Type, Container>& rhs) const
{
    return (_container < rhs._container);
}

// Greater than or equal to
/// \brief          Greater than or equal to.
/// \details        Lexicographically compares the values in the container for
///                 the greater than or equal to comparison operator. 
/// \param[in]      rhs Right-hand side.
/// \return         True if greater than or equal to, false otherwise. 
template <class Type, class Container>
inline bool Evolution<Type, Container>::operator>=(const Evolution<Type, Container>& rhs) const
{
    return (_container >= rhs._container);
}

// Less than or equal to
/// \brief          Less than or equal to.
/// \details        Lexicographically compares the values in the container for
///                 the less than or equal to comparison operator. 
/// \param[in]      rhs Right-hand side.
/// \return         True if less than or equal to, false otherwise. 
template <class Type, class Container>
inline bool Evolution<Type, Container>::operator<=(const Evolution<Type, Container>& rhs) const
{
    return (_container <= rhs._container);
}

// Direct access to the element
/// \brief          Direct access to the element.
/// \details        Provides a direct access to the specified element.
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <class Type, class Container>
inline Type& Evolution<Type, Container>::operator[](const unsigned int i)
{
    return _container[i];
}

// Immutable direct access to the element
/// \brief          Immutable direct access to the element.
/// \details        Provides a direct access to the specified element.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <class Type, class Container>
inline const Type& Evolution<Type, Container>::operator[](const unsigned int i) const
{
    return _container[i];
}
//--------------------------------------------------------------------------- //



//-------------------------------- ASSIGNMENT ------------------------------- //
// Copy assignment
/// \brief          Copy assignment.
/// \details        Provides a copy assignment from another evolution container.
/// \param[in]      source Source of the copy.  
/// \return         Self reference.
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::assign(const Evolution<Type, Container>& source)
{
    _container = source._container;
    return *this;
}

// Initializer list assignment
/// \brief          Initializer list assignment.
/// \details        Provides an assignment from an initializer list equivalent 
///                 to a call to the constructor.
/// \tparam         OtherType (Other step type.)
/// \param[in]      source Source of the copy.
/// \return         Self reference.
template <class Type, class Container>
template <class OtherType> 
inline Evolution<Type, Container>& Evolution<Type, Container>::assign(const std::initializer_list<OtherType>& source)
{
    _container = source;
    return *this;
}

// Conversion assignment  
/// \brief          Conversion assignment. 
/// \details        Provides an assignment from another type of evolution 
///                 container.
/// \tparam         OtherType (Other step type.)
/// \tparam         OtherContainer (Other container type.)
/// \param[in]      source Source of the copy. 
/// \return         Self reference. 
template <class Type, class Container>
template <class OtherType, class OtherContainer>
inline Evolution<Type, Container>& Evolution<Type, Container>::assign(const Evolution<OtherType, OtherContainer>& source)
{
    const unsigned int n = source.container().size();
    _container.resize(n);
    for (unsigned int i = 0; i < n; ++i) {
        _container[i] = source.container()[i];
    }
    return *this;
}

// Generic assignment 
/// \brief          Generic assignment. 
/// \details        Provides a generic interface to all assignments of the 
///                 container.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Self reference. 
template <class Type, class Container>
template <class... Misc, class> 
inline Evolution<Type, Container>& Evolution<Type, Container>::assign(Misc&&... misc)
{
    _container = Container(std::forward<Misc>(misc)...);
    return *this;
}
//--------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Nullify 
/// \brief          Nullify. 
/// \details        Resets all data members to their default values. 
/// \return         Self reference. 
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::nullify()
{
    std::fill(_container.begin(), _container.end(), Type());
    return *this;
}

// Copy 
/// \brief          Copy. 
/// \details        Generates a copy of the object. 
/// \return         Copy. 
template <class Type, class Container>
inline Evolution<Type, Container> Evolution<Type, Container>::copy() const
{
    return *this;
}

// Cast 
/// \brief          Cast. 
/// \details        Casts contents to another object type. 
/// \tparam         OtherType Other step type.
/// \tparam         OtherContainer Other container type.
/// \return         Casted copy. 
template <class Type, class Container>
template <class OtherType, class OtherContainer> 
inline Evolution<OtherType, OtherContainer> Evolution<Type, Container>::cast() const
{
    return Evolution<OtherType, OtherContainer>(*this);
}
//--------------------------------------------------------------------------- //



//---------------------------------- ACCESS --------------------------------- //
// Access with range-check
/// \brief          Access with range-check.
/// \details        Provides an access to the element with a range-check.
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <class Type, class Container>
inline Type& Evolution<Type, Container>::at(const unsigned int i)
{
    return _container.at(i);
}

// Immutable access with range-check
/// \brief          Immutable access with range-check.
/// \details        Provides an access to the element with a range-check.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <class Type, class Container>
inline const Type& Evolution<Type, Container>::at(const unsigned int i) const
{
    return _container.at(i);
}

// Access to the i-th element from the beginning
/// \brief          Access to the i-th element from the beginning.
/// \details        Returns a reference to the i-th first element in the 
///                 container without doing any range check.
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <class Type, class Container>
inline Type& Evolution<Type, Container>::front(const unsigned int i)
{
    return _container[i];
}

// Immutable access to the i-th element from the beginning
/// \brief          Immutable access to the i-th element from the beginning.
/// \details        Returns a reference to the i-th first element in the 
///                 container without doing any range check.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <class Type, class Container>
inline const Type& Evolution<Type, Container>::front(const unsigned int i) const 
{
    return _container[i];
}

// Access to the i-th element from the end
/// \brief          Access to the i-th element from the end.
/// \details        Returns a reference to the i-th last element in the 
///                 container without doing any range check.
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <class Type, class Container>
inline Type& Evolution<Type, Container>::back(const unsigned int i)
{
    return _container[static_cast<int>(_container.size())-1-static_cast<int>(i)];
}

// Immutable access to the i-th element from the end
/// \brief          Immutable access to the i-th element from the end.
/// \details        Returns a reference to the i-th last element in the 
///                 container without doing any range check.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <class Type, class Container>
inline const Type& Evolution<Type, Container>::back(const unsigned int i) const
{
    return _container[static_cast<int>(_container.size())-1-static_cast<int>(i)];
}

// Cyclic access to elements
/// \brief          Cyclic access to elements.
/// \details        Provides a cyclic access to the elements, using the index
///                 modulo. Negative indexes are supported. It allows to iterate
///                 several times over the contents just by incrementing the
///                 provided index.
/// \param[in]      i Index of the element.
/// \return         Reference to the element.
template <class Type, class Container>
inline Type& Evolution<Type, Container>::cycle(const int i)
{
    const int isize = _container.size();
    return _container[(i >= 0) ? (i%isize) : ((isize-((-i)%isize))%isize)];
}

// Immutable cyclic access to elements
/// \brief          Immutable cyclic access to elements.
/// \details        Provides a cyclic access to the elements, using the index
///                 modulo. Negative indexes are supported. It allows to iterate
///                 several times over the contents just by incrementing the
///                 provided index.
/// \param[in]      i Index of the element.
/// \return         Immutable reference to the element.
template <class Type, class Container>
inline const Type& Evolution<Type, Container>::cycle(const int i) const
{
    const int isize = _container.size();
    return _container[(i >= 0) ? (i%isize) : ((isize-((-i)%isize))%isize)];
}

// Direct access to the underlying container 
/// \brief          Direct access to the underlying container.
/// \details        Provides a direct access to the underlying container by 
///                 returning a reference to it.
/// \return         Reference to the underlying container.
template <class Type, class Container>
inline Container& Evolution<Type, Container>::container()
{
    return _container;
}

// Immutable direct access to the underlying container 
/// \brief          Direct access to the underlying container.
/// \details        Provides a direct access to the underlying container by 
///                 returning a reference to it.
/// \return         Immutable reference to the underlying container.
template <class Type, class Container>
inline const Container& Evolution<Type, Container>::container() const
{
    return _container;
}

// Direct access to the underlying array 
/// \brief          Direct access to the underlying array.
/// \details        Provides a direct access to the underlying array by 
///                 returning a pointer to the first element of storage.
/// \return         Pointer to the underlying element storage.
template <class Type, class Container>
inline Type* Evolution<Type, Container>::data()
{
    return _container.data();
}

// Immutable direct access to the underlying array 
/// \brief          Immutable direct access to the underlying array.
/// \details        Provides a direct access to the underlying array by 
///                 returning a pointer to the first element of storage.
/// \return         Immutable pointer to the underlying element storage.
template <class Type, class Container>
inline const Type* Evolution<Type, Container>::data() const
{
    return _container.data();
}
//--------------------------------------------------------------------------- //



// ------------------------------- ITERATORS -------------------------------- //
// Iterator to the beginning
/// \brief          Iterator to the beginning.
/// \details        Returns a pointer to the first element.
/// \return         Pointer to the beginning.
template <class Type, class Container>
inline Type* Evolution<Type, Container>::begin()
{
    return _container.data();
}

// Immutable iterator to the beginning
/// \brief          Immutable iterator to the beginning.
/// \details        Returns a pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <class Type, class Container>
inline const Type* Evolution<Type, Container>::begin() const
{
    return _container.data();
}

// Forced immutable iterator to the beginning
/// \brief          Forced immutable iterator to the beginning.
/// \details        Returns a constant pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <class Type, class Container>
inline const Type* Evolution<Type, Container>::cbegin() const
{
    return _container.data();
}

// Iterator to the end
/// \brief          Iterator to the end.
/// \details        Returns a pointer to the position after the last element.
/// \return         Pointer to the end.
template <class Type, class Container>
inline Type* Evolution<Type, Container>::end()
{
    return _container.data()+_container.size();
}

// Immutable iterator to the end
/// \brief          Immutable iterator to the end.
/// \details        Returns a pointer to the position after the last element.
/// \return         Immutable pointer to the end.
template <class Type, class Container>
inline const Type* Evolution<Type, Container>::end() const
{
    return _container.data()+_container.size();
}

// Forced immutable iterator to the end
/// \brief          Forced immutable iterator to the end.
/// \details        Returns a constant pointer to the position after the last 
///                 element.
/// \return         Immutable pointer to the end.
template <class Type, class Container>
inline const Type* Evolution<Type, Container>::cend() const
{
    return _container.data()+_container.size();
}

// Reverse iterator to the beginning
/// \brief          Reverse iterator to the beginning.
/// \details        Returns a reversed pointer to the position after the last 
///                 element.
/// \return         Pointer to the end.
template <class Type, class Container>
inline std::reverse_iterator<Type*> Evolution<Type, Container>::rbegin()
{
    return std::reverse_iterator<Type*>(_container.data()+_container.size());
}

// Immutable reverse iterator to the beginning
/// \brief          Immutable reverse iterator to the beginning.
/// \details        Returns a reversed pointer to the position after the last 
///                 element.
/// \return         Immutable pointer to the end.
template <class Type, class Container>
inline std::reverse_iterator<const Type*> Evolution<Type, Container>::rbegin() const
{
    return std::reverse_iterator<const Type*>(_container.data()+_container.size());
}

// Forced immutable reverse iterator to the beginning
/// \brief          Forced immutable reverse iterator to the beginning.
/// \details        Returns a constant reversed pointer to the position after 
///                 the last element.
/// \return         Immutable pointer to the end.
template <class Type, class Container>
inline std::reverse_iterator<const Type*> Evolution<Type, Container>::crbegin() const
{
    return std::reverse_iterator<const Type*>(_container.data()+_container.size());
}

// Reverse iterator to the end
/// \brief          Reverse iterator to the end.
/// \details        Returns a reversed pointer to the first element.
/// \return         Pointer to the beginning.
template <class Type, class Container>
inline std::reverse_iterator<Type*> Evolution<Type, Container>::rend()
{
    return std::reverse_iterator<Type*>(_container.data());
}

// Immutable reverse iterator to the end
/// \brief          Immutable reverse iterator to the end.
/// \details        Returns a reversed pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <class Type, class Container>
inline std::reverse_iterator<const Type*> Evolution<Type, Container>::rend() const
{
    return std::reverse_iterator<const Type*>(_container.data());
}

// Forced immutable reverse iterator to the end
/// \brief          Forced immutable reverse iterator to the end.
/// \details        Returns a constant reversed pointer to the first element.
/// \return         Immutable pointer to the beginning.
template <class Type, class Container>
inline std::reverse_iterator<const Type*> Evolution<Type, Container>::crend() const
{
    return std::reverse_iterator<const Type*>(_container.data());
}

// Index of an iterator in the container
/// \brief          Index of an iterator in the container.
/// \details        Returns the index of the element pointed by an iterator or 
///                 a pointer.
/// \tparam         Iterator (Pointer or iterator type.)
/// \param[in]      it Iterator to the element.
/// \return         Index of the element.
/// \exception      std::out_of_range Out of range.
template <class Type, class Container>
template <typename Iterator> 
inline unsigned int Evolution<Type, Container>::index(const Iterator& it) const
{
    if ((&*it < _container.data()) || (&*it >= _container.data()+_container.size())) {
        throw std::out_of_range("ERROR = Evolution::index() : out of range");
    }
    return &*it-_container.data();
}
// -------------------------------------------------------------------------- //



//--------------------------------- CAPACITY -------------------------------- //
// Emptiness checking
/// \brief          Emptiness checking.
/// \details        Checks if the container has no elements.
/// \return         True if empty, false otherwise.
template <class Type, class Container>
inline bool Evolution<Type, Container>::empty() const
{
    return _container.empty();
}

// Number of elements 
/// \brief          Number of elements.
/// \details        Returns the distance between the first and the last element.
/// \return         The number of elements in the container.
template <class Type, class Container>
inline unsigned int Evolution<Type, Container>::size() const
{
    return _container.size();
}

// Capacity of the underlying storage
/// \brief          Capacity of the underlying storage.
/// \details        Returns the number of elements that the container has 
///                 currently allocated space for. 
/// \return         Capacity of the currently allocated storage. 
template <class Type, class Container>
inline unsigned int Evolution<Type, Container>::capacity() const
{
    return _container.capacity();
}

// Storage reservation
/// \brief          Storage reservation.
/// \details        Increases the capacity of the underlying storage. Existing
///                 elements are protected so it could not invalidate the 
///                 actual contents.
/// \param[in]      n New size.
/// \return         Self reference.
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::reserve(const unsigned int n)
{
    _container.reserve((n < _container.size()) ? (_container.size()) : (n));
    return *this;
}

// Storage shrinking
/// \brief          Storage shrinking.
/// \details        Reduces memory usage by freeing unused memory.
/// \return         Self reference.
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::shrink()
{
    _container.shrink_to_fit();
    return *this;
}

// Available space
/// \brief          Available space.
/// \details        Returns the maximum possible number of elements.
/// \return         Maximum number of elements. 
template <class Type, class Container>
inline unsigned long long int Evolution<Type, Container>::space() const
{
    return _container.max_size();
}
//--------------------------------------------------------------------------- //        



//--------------------------------- MODIFIERS ------------------------------- //
// Clear the contents
/// \brief          Full clear the container.
/// \details        Erase the container.
/// \return         Self reference.
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::fullclear()
{
    Container().swap(_container);
    return *this;
}


// Clear the contents
/// \brief          Clear the contents.
/// \details        Removes all elements from the container.
/// \return         Self reference.
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::clear()
{
    _container.clear();
    return *this;
}

// Resize
/// \brief          Resize.
/// \details        Resizes the container to contain a new number of elements. 
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Self reference.
template <class Type, class Container>
template <class... Misc> 
inline Evolution<Type, Container>& Evolution<Type, Container>::resize(const Misc&... misc)
{
    _container.resize(misc...);
    return *this;
}

// Pop back
/// \brief          Pop back.
/// \details        Removes the last element of the container.
/// \return         Self reference.
template <class Type, class Container>
inline Evolution<Type, Container>& Evolution<Type, Container>::pop()
{
    _container.pop_back();
    return *this;
}

// Append an element
/// \brief          Append an element.
/// \details        Appends an element to the end of the container.
/// \tparam         OtherType (Other step type.)
/// \param[in]      value Value of the element to be pushed back.
/// \return         Self reference.
template <class Type, class Container>
template <class OtherType, class> 
inline Evolution<Type, Container>& Evolution<Type, Container>::append(const OtherType& value)
{
    _container.push_back(value);
    return *this;
}
//--------------------------------------------------------------------------- //



//---------------------------------- STREAM --------------------------------- //
// Output stream operator 
/// \brief          Output stream operator. 
/// \details        Adds each element to the stream using the <tt>fill()</tt>
///                 character as a separator.
/// \tparam         SelfType (Step type.)
/// \tparam         SelfContainer (Container type.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side object.
/// \return         Output stream. 
template <class SelfType, class SelfContainer>
std::ostream& operator<<(std::ostream& lhs, const Evolution<SelfType, SelfContainer>& rhs)
{
    const char separator = lhs.fill();
    const int isize = static_cast<int>(rhs._container.size())-1;
    for (int i = 0; i < isize; ++i) {
        lhs<<rhs._container[i]<<separator; 
    }
    return (isize < 0) ? (lhs) : (lhs<<rhs[isize]);
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of Evolution. 
/// \return         0 if no error. 
template <class Type, class Container>
int Evolution<Type, Container>::example()
{
    // Initialize
    std::cout<<"BEGIN = Evolution::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    unsigned int ui = 42;
    std::array<std::array<double, 3>, 2> arr({{std::array<double, 3>({{4., 8., 15.}}), std::array<double, 3>({{16., 23., 42.}})}});
    std::tuple<std::string> tup("string");
    Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > step(ui, arr, tup);
    Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > step2(ui*2, arr, tup);
    Step<int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > istep(ui, arr, tup);
    std::vector<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > > vec({step, step2});
    
    // Construction
    Evolution<double> basic({static_cast<double>(ui)});
    Evolution<double> basic2({static_cast<double>(ui*2)});
    Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > > evolution({step});
    Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > > evolution2({step, step});
    Evolution<Step<int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > > ievolution({istep});
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Lifecycle : "                                                                                                               <<std::endl;
    std::cout<<std::setw(width*3)<<"Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >() : "                        <<Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >()<<std::endl;
    std::cout<<std::setw(width*3)<<"Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >(evolution) : "               <<Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >(evolution)<<std::endl;    
    std::cout<<std::setw(width*3)<<"Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >(ievolution) : "              <<Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >(ievolution)<<std::endl;    
    std::cout<<std::setw(width*3)<<"Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >({step}) : "                  <<Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >({step})<<std::endl;
    std::cout<<std::setw(width*3)<<"Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >(2) : "                       <<Evolution<Step<unsigned int, std::array<std::array<double, 3>, 2>, std::tuple<std::string> > >(2)<<std::endl;  

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"evolution = evolution2"                       <<(evolution = evolution2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution = {step}"                           <<(evolution = {step2})<<std::endl;
    std::cout<<std::setw(width)<<"evolution = vec"                              <<(evolution = vec)<<std::endl;
    std::cout<<std::setw(width)<<"evolution == evolution2"                      <<(evolution == evolution2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution != evolution2"                      <<(evolution != evolution2)<<std::endl;
    std::cout<<std::setw(width)<<"basic > basic2"                               <<(basic > basic2)<<std::endl;
    std::cout<<std::setw(width)<<"basic < basic2"                               <<(basic < basic2)<<std::endl;
    std::cout<<std::setw(width)<<"basic >= basic2"                              <<(basic >= basic2)<<std::endl;
    std::cout<<std::setw(width)<<"basic <= basic2"                              <<(basic <= basic2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution[0] = step"                          <<(evolution[0] = step)<<std::endl;
    std::cout<<std::setw(width)<<"evolution[0]"                                 <<evolution[0]<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"evolution.assign(evolution2)"                 <<evolution.assign(evolution2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.assign({step})"                     <<evolution.assign({step})<<std::endl;
    std::cout<<std::setw(width)<<"evolution.assign(ievolution)"                 <<evolution.assign(ievolution)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.assign(4)"                          <<evolution.assign(4)<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"evolution.nullify()"                          <<evolution.nullify()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.copy()"                             <<evolution.copy()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.cast()"                             <<evolution.cast()<<std::endl;

    // Access
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Access : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"evolution.at(0) = step"                       <<(evolution.at(0) = step)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.at(0)"                              <<evolution.at(0)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.front(1) = step"                    <<(evolution.front(1) = step)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.front(1)"                           <<evolution.front(1)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.back(1) = step2"                    <<(evolution.back(1) = step2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.back(1)"                            <<evolution.back(1)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.cycle(-1) = step"                   <<(evolution.cycle(-1) = step)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.cycle(-1)"                          <<evolution.cycle(-1)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.container()[0] = step2"             <<(evolution.container()[0] = step2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.container()[0]"                     <<evolution.container()[0]<<std::endl;
    std::cout<<std::setw(width)<<"evolution.data()[0] = step"                   <<(evolution.data()[0] = step)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.data()[0]"                          <<evolution.data()[0]<<std::endl;

    // Iterators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Iterators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.begin() : "                       <<&*evolution.begin()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.cbegin() : "                      <<&*evolution.cbegin()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.end() : "                         <<&*evolution.end()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.cend() : "                        <<&*evolution.cend()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.rbegin() : "                      <<&*evolution.rbegin()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.crbegin() : "                     <<&*evolution.crbegin()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.rend() : "                        <<&*evolution.rend()<<std::endl;
    std::cout<<std::setw(width)<<"&*evolution.crend() : "                       <<&*evolution.crend()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.index(&evolution[3]) : "            <<evolution.index(&evolution[3])<<std::endl;
    
    // Capacity
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Capacity : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"evolution.empty() : "                         <<evolution.empty()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.size() : "                          <<evolution.size()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.capacity() : "                      <<evolution.capacity()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.reserve(42) : "                     <<evolution.reserve(42)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.shrink() : "                        <<evolution.shrink()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.space() : "                         <<evolution.space()<<std::endl;

    // Modifiers
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Modifiers : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"evolution.clear() : "                         <<evolution.clear()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.resize(2) : "                       <<evolution.resize(2)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.resize(4, step) : "                 <<evolution.resize(4, step)<<std::endl;
    std::cout<<std::setw(width)<<"evolution.pop() : "                           <<evolution.pop()<<std::endl;
    std::cout<<std::setw(width)<<"evolution.append() : "                        <<evolution.append(step)<<std::endl;
    
    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, evolution) : "          <<evolution<<std::endl;
        
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Evolution::example()"<<std::endl;
    return 0;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // EVOLUTION_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

