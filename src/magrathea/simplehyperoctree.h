/* *************************** SIMPLEHYPEROCTREE **************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          SimpleHyperOctree 
// DESCRIPTION :    A simple hyperoctree based on bit manipulations
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2021), Michel-Andrès Breton (2015-2016)]
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           simplehyperoctree.h 
/// \brief          A simple hyperoctree based on bit manipulations
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2021
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef SIMPLEHYPEROCTREE_H_INCLUDED
#define SIMPLEHYPEROCTREE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

 

//------------------------------- PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <sstream>
#include <type_traits>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <algorithm>
#include <vector>
#include <bitset>
#include <array>
#include <tuple>
// Include libs
#include "simplehyperoctreeindex.h"
#include "../utility.h"
// Include project
// Misc
namespace magrathea {
//--------------------------------------------------------------------------- //



//---------------------------------- CLASS ---------------------------------- //
// A simple hyperoctree based on bit manipulations
/// \brief          A simple hyperoctree based on bit manipulations.
/// \details        Implementation of a simple and easy-to-use and hyperoctree 
///                 structure in arbitrary dimension. It provides basic find
///                 and search algorithms based on indices relying on bit 
///                 manipulations.
/// \tparam         Type Scalar position type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         Element Underlying element type.
/// \tparam         Container Underlying container type.
template <typename Type = double, class Index = SimpleHyperOctreeIndex<unsigned long long int, 3>, class Data = double, unsigned int Dimension = Index::dimension(), class Position = std::ratio<0>, class Extent = std::ratio<1>, class Element = std::pair<Index, Data>, class Container = std::vector<Element> >
class SimpleHyperOctree final
{
    // Setup
    static_assert(std::is_floating_point<Type>::value, "ERROR = SimpleHyperOctree : non compliant type");
    static_assert(Index::dimension() == Dimension, "ERROR = SimpleHyperOctree : non compliant index");
    static_assert(Dimension > 0, "ERROR = SimpleHyperOctree : non compliant dimension");
    static_assert(Position::den != 0, "ERROR = SimpleHyperOctree : non compliant position");
    static_assert(Extent::num >= 0, "ERROR = SimpleHyperOctree : non compliant extent");
    static_assert(std::is_same<typename std::remove_cv<typename std::remove_reference<decltype(*std::begin(std::declval<Container>()))>::type>::type, Element>::value, "ERROR = SimpleHyperOctree : non compliant container");

    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        inline SimpleHyperOctree();
        explicit inline SimpleHyperOctree(const unsigned int ilvl, const unsigned int nref = 0);
    //@}

    // Operators
    /// \name           Operators
    //@{
    public:
        inline bool operator==(const SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& rhs) const;
        inline bool operator!=(const SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& rhs) const;
        inline Element& operator[](const unsigned long long int ielem);
        inline const Element& operator[](const unsigned long long int ielem) const;
        inline Container& operator()();
        inline const Container& operator()() const;
        inline Element& operator()(const Index& idx);
        inline const Element& operator()(const Index& idx) const;
        template <typename... Types, class = typename std::enable_if<(sizeof...(Types) != 0)>::type, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Element& operator()(Types&&... iposs);
        template <typename... Types, class = typename std::enable_if<(sizeof...(Types) != 0)>::type, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline const Element& operator()(Types&&... iposs) const;
    //@}

    // Assignment 
    /// \name           Assignment 
    //@{
    public: 
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& assign();
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& assign(const unsigned int ilvl, const unsigned int nref = 0);
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& assign(const SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& source);
    //@}

    // Management 
    /// \name           Management 
    //@{
    public: 
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& nullify();
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container> copy() const;
        template <class Template = SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>, class = typename std::enable_if<std::is_constructible<Template, SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container> >::value>::type> inline Template cast() const;
    //@}

    // Access 
    /// \name           Access 
    //@{
    public: 
        inline Element& at(const unsigned long long int ielem);
        inline const Element& at(const unsigned long long int ielem) const;  
        inline Element& front(const unsigned long long int ielem = 0);
        inline const Element& front(const unsigned long long int ielem = 0) const;  
        inline Element& back(const unsigned long long int ielem = 0);
        inline const Element& back(const unsigned long long int ielem = 0) const;  
        inline Element& cycle(const long long int ielem);
        inline const Element& cycle(const long long int ielem) const; 
        inline Container& container();
        inline const Container& container() const;
        inline Element* data();
        inline const Element* data() const;

    //@}
    
    // Iterators
    /// \name           Iterators
    //@{
    public:
        template <typename Iterator = decltype(std::declval<Container>().begin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator begin();
        template <typename Iterator = decltype(std::declval<const Container>().begin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator begin() const;
        template <typename Iterator = decltype(std::declval<Container>().cbegin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator cbegin() const;
        template <typename Iterator = decltype(std::declval<Container>().end()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator end();
        template <typename Iterator = decltype(std::declval<const Container>().end()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator end() const;
        template <typename Iterator = decltype(std::declval<Container>().cend()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator cend() const;
        template <typename Iterator = decltype(std::declval<Container>().rbegin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator rbegin();
        template <typename Iterator = decltype(std::declval<const Container>().rbegin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator rbegin() const;
        template <typename Iterator = decltype(std::declval<Container>().crbegin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator crbegin() const;
        template <typename Iterator = decltype(std::declval<Container>().rend()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator rend();
        template <typename Iterator = decltype(std::declval<const Container>().rend()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator rend() const;
        template <typename Iterator = decltype(std::declval<Container>().crend()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator crend() const;
        template <typename Iterator, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline unsigned long long int index(const Iterator& it) const;
    //@}
    
    // Search
    /// \name           Search
    //@{
    public:
        template <typename Iterator = decltype(std::declval<Container>().begin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator find(const Index& idx);
        template <typename Iterator = decltype(std::declval<const Container>().begin()), class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline Iterator find(const Index& idx) const;
        template <typename Iterator = decltype(std::declval<Container>().begin()), typename... Types, class = typename std::enable_if<(std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value) && ((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Iterator locate(Types&&... iposs);
        template <typename Iterator = decltype(std::declval<const Container>().begin()), typename... Types, class = typename std::enable_if<(std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value) && ((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Iterator locate(Types&&... iposs) const;
    //@}

    // Capacity 
    /// \name           Capacity 
    //@{
    public: 
        inline bool empty() const;
        inline unsigned long long int size() const;
        inline unsigned long long int capacity() const;
        inline unsigned long long int space() const;
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& shrink();
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& reserve(const unsigned long long int nelem);
    //@}

    // Modifiers  
    /// \name           Modifiers  
    //@{
    public: 
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& fullclear();
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& clear();
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& pop();
        template <class... Misc> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& append(Misc&&... misc);
        template <class... Misc> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& resize(Misc&&... misc);

	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& rho(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& rho(Misc&&... misc) const;
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& phi(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& phi(Misc&&... misc) const;
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidxyz(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidxyz(Misc&&... misc) const;
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidx(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidx(Misc&&... misc) const;
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidy(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidy(Misc&&... misc) const;
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidz(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& dphidz(Misc&&... misc) const;
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& a(Misc&&... misc);
	template <unsigned int... Values,class... Misc,class Template, class> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& a(Misc&&... misc) const;

    //@}
    
    // Refinement
    /// \name           Refinement
    //@{
    public:
        inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& update();
        template <typename Iterator, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline bool root(const Iterator& it);
        template <typename Iterator, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline bool leaf(const Iterator& it);
        template <typename Iterator, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& refine(const Iterator& it);
        template <typename Iterator, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<decltype(*std::declval<Iterator>())>::type>::type, Element>::value>::type> inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& coarsen(const Iterator& it);
    //@}

    // Interpolation
    /// \name           Interpolation
    //@{
    public:
        template <typename... Types, class = typename std::enable_if<(sizeof...(Types) != 0)>::type, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Data ngp(Types&&... iposs) const; 
        template <typename... Types, class = typename std::enable_if<(sizeof...(Types) != 0)>::type, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Data cic(Types&&... iposs) const; 
        template <typename... Types, class = typename std::enable_if<(sizeof...(Types) != 0)>::type, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Data tsc(Types&&... iposs) const; 
        template <typename... Types, class = typename std::enable_if<(sizeof...(Types) != 0)>::type, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1)) && (!std::is_same<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Index>::value)>::type> inline Data tsc(std::vector<Element>& elemsTsc, Types&&... iposs) const; 
    //@}

    // Stream 
    /// \name           Stream 
    //@{
    public: 
        template <typename SelfType, class SelfIndex, class SelfData, unsigned int SelfDimension, class SelfPosition, class SelfExtent, class SelfElement, class SelfContainer> friend std::ostream& operator<<(std::ostream& lhs, const SimpleHyperOctree<SelfType, SelfIndex, SelfData, SelfDimension, SelfPosition, SelfExtent, SelfElement, SelfContainer>& rhs); 
    //@}

    // Properties 
    /// \name           Properties 
    //@{
    public: 
        static constexpr Type type(); 
        static constexpr Position position(); 
        static constexpr Extent extent(); 
        static constexpr Element element(); 
        static constexpr unsigned int dimension();
    //@}
    
    // Helpers
    /// \name           Helpers
    //@{
    public:
        template <class Template, class Tuple, typename... Kinds, class = typename std::enable_if<(sizeof...(Kinds)+1 <= std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value)>::type> static constexpr Template tupleify(Tuple&& tuple, Kinds&&... values);
        template <class Template, class Tuple, typename... Kinds, class = typename std::enable_if<(sizeof...(Kinds) == std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value)>::type> static constexpr Template tupleify(const Tuple& tuple, Kinds&&... values);
        template <class Template, typename... Kinds, class = typename std::enable_if<std::is_convertible<typename std::tuple_element<0, typename std::conditional<static_cast<bool>(sizeof...(Kinds)), std::tuple<typename std::remove_cv<typename std::remove_reference<Kinds>::type>::type...>, std::tuple<Template> >::type>::type, typename std::tuple_element<0, typename std::conditional<static_cast<bool>(sizeof...(Kinds)), Template, std::tuple<Template> >::type>::type>::value>::type> static constexpr Template tupleify(Kinds&&... values);
        template <class Template, class Function, class Tuple, typename... Kinds, class = typename std::enable_if<(sizeof...(Kinds)+1 <= std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value)>::type> static constexpr Template variadify(Function&& f, Tuple&& tuple, Kinds&&... values);
        template <class Template, class Function, class Tuple, typename... Kinds, class = typename std::enable_if<(sizeof...(Kinds) == std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value)>::type> static constexpr Template variadify(Function&& f, const Tuple& tuple, Kinds&&... values);
        template <class Template, class Function, typename... Kinds, class = typename std::enable_if<(!std::is_void<typename std::result_of<Function(Kinds...)>::type>::value) || (std::is_void<typename std::result_of<Function(Kinds...)>::type>::value)>::type> static constexpr Template variadify(Function&& f, Kinds&&... values);
        template <unsigned int Step = 0, class Object, typename Kind, class = typename std::enable_if<(Step <= std::tuple_size<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Object>().data())>::type>::type>::value)>::type> static constexpr Object& mac(Object& accumulator, const Object& factor, const Kind coefficient);
        template <unsigned int Step = 0, class Tuple, typename Kind, class Dummy = std::nullptr_t, class = typename std::enable_if<(Step <= std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value)>::type> static constexpr Tuple& mac(Tuple& accumulator, const Tuple& factor, const Kind coefficient, Dummy&& = Dummy());
        template <class Regular, typename Kind, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (std::is_convertible<Kind, Regular>::value)>::type> static constexpr Regular& mac(Regular& accumulator, const Regular& factor, const Kind coefficient, Dummy...);
        template <class Irregular, typename Kind, class... Dummy, class = typename std::enable_if<(sizeof...(Dummy) == 0) && (!std::is_convertible<Kind, Irregular>::value)>::type> static constexpr Irregular& mac(Irregular& accumulator, const Irregular& factor, const Kind coefficient, const Dummy...);
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
// Implicit empty constructor 
/// \brief          Implicit empty constructor. 
/// \details        Provides an implicit construction of an empty hyperoctree.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::SimpleHyperOctree()
: _container(Container())
{
    ;
}

// Explicit level constructor 
/// \brief          Explicit level constructor. 
/// \details        Provides an explicit construction of a fixed mesh starting
///                 from the specified level and refining the provided amount of
///                 times.
/// \param[in]      ilvl Index of the first level.
/// \param[in]      nref Number of refinements of the first level.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::SimpleHyperOctree(const unsigned int ilvl, const unsigned int nref)
: _container(Container())
{
    static const Element zero = Element(Index(), Data());
    Index idx = Index();
    idx = idx.previous(ilvl, nref);
    _container.resize(Index::indices(ilvl, nref), zero);
    std::for_each(std::begin(_container), std::end(_container), [=, &ilvl, &nref, &idx](Element& elem){idx = idx.next(ilvl, nref); std::get<0>(elem) = idx;});
} 
//--------------------------------------------------------------------------- //



//-------------------------------- OPERATORS -------------------------------- //
// Equal to
/// \brief          Equal to.
/// \details        Compares the hyperoctrees for equality.
/// \param[in]      rhs Right-hand side.
/// \return         True if equal, false if not equal. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline bool SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator==(const SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& rhs) const
{
    return (_container == rhs._container);
}

// Not equal to
/// \brief          Not equal to.
/// \details        Compares the hyperoctrees for difference.
/// \param[in]      rhs Right-hand side.
/// \return         True if not equal, false if equal. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline bool SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator!=(const SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& rhs) const
{
    return (_container != rhs._container);
}

// Element access operator
/// \brief          Element access operator.
/// \details        Provides a direct access to the specified element.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator[](const unsigned long long int ielem)
{
    return _container[ielem];
}

// Immutable element access operator
/// \brief          Immutable element access operator.
/// \details        Provides a direct access to the specified element.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator[](const unsigned long long int ielem) const
{
    return _container[ielem];
}

// Container access operator
/// \brief          Container access operator.
/// \details        Provides a direct access to the underlying container.
/// \return         Reference to the container.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Container& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator()()
{
    return _container;
}

// Immutable container access operator
/// \brief          Immutable container access operator.
/// \details        Provides a direct access to the underlying container.
/// \return         Immutable reference to the container.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Container& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator()() const
{
    return _container;
}

// Element access operator from hyperoctree index
/// \brief          Element access operator from hyperoctree index.
/// \details        Provides a direct access to the element corresponding to
///                 the specified index.
/// \param[in]      idx Index value.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator()(const Index& idx)
{
    return *(std::upper_bound(std::begin(_container), std::end(_container), Element(idx, Data()), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);})-1);
}

// Immutable element access operator from hyperoctree index
/// \brief          Immutable element access operator from hyperoctree index.
/// \details        Provides a direct access to the element corresponding to
///                 the specified index.
/// \param[in]      idx Index value.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator()(const Index& idx) const
{
    return *(std::upper_bound(std::begin(_container), std::end(_container), Element(idx, Data()), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);})-1);
}

// Element access operator from position
/// \brief          Element access operator from position.
/// \details        Provides a direct access to the element corresponding to
///                 the specified position.
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename... Types, class, class> 
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator()(Types&&... iposs)
{
    return *(std::upper_bound(std::begin(_container), std::end(_container), Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data()), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);})-1);
}

// Immutable element access operator from position
/// \brief          Immutable element access operator from position.
/// \details        Provides a direct access to the element corresponding to
///                 the specified position.
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename... Types, class, class> 
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::operator()(Types&&... iposs) const
{
    return *(std::upper_bound(std::begin(_container), std::end(_container), Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data()), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);})-1);
}
//--------------------------------------------------------------------------- //



//-------------------------------- ASSIGNMENT ------------------------------- //
// Empty assignment 
/// \brief          Empty assignment. 
/// \details        Provides an empty assignment of an empty hyperoctree.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::assign()
{
    _container.clear();
    return *this;
}

// Level assignment 
/// \brief          Level assignment . 
/// \details        Provides an assignment of a fixed mesh starting from the
///                 specified level and refining the provided amount of times.
/// \param[in]      ilvl Index of the first level.
/// \param[in]      nref Number of refinements of the first level.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::assign(const unsigned int ilvl, const unsigned int nref)
{
    static const Element zero = Element(Index(), Data());
    Index idx = Index();
    idx = idx.previous(ilvl, nref);
    _container.resize(Index::indices(ilvl, nref), zero);
    std::for_each(std::begin(_container), std::end(_container), [=, &ilvl, &nref, &idx](Element& elem){idx = idx.next(ilvl, nref); std::get<0>(elem) = idx;});
    return *this;
}

// Copy assignment 
/// \brief          Copy assignment . 
/// \details        Provides a copy assignment from another hyperoctree.
/// \param[in]      source Source of the copy.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::assign(const SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& source)
{
    _container = source._container;
    return *this;
}
//--------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Nullify 
/// \brief          Nullify. 
/// \details        Resets all data to their default values but keeping the
///                 same tree structure. 
/// \return         Self reference. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::nullify()
{
    const Data zero = Data();
    std::for_each(std::begin(_container), std::end(_container), [=, &zero](Element& elem){std::get<1>(elem) = zero;});
    return *this;
}

// Copy 
/// \brief          Copy. 
/// \details        Generates a copy of the hyperoctree. 
/// \return         Copy. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container> SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::copy() const
{
    return *this;
}

// Cast 
/// \brief          Cast. 
/// \details        Casts contents to another object type. 
/// \tparam         Template Output type.
/// \return         Casted copy. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, class> 
inline Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cast() const
{
    return Template(*this);
}
//--------------------------------------------------------------------------- //



//---------------------------------- ACCESS --------------------------------- //
// Access to an element with range-check
/// \brief          Access to an element with range-check.
/// \details        Provides an access to the specified element with a 
///                 range-check.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::at(const unsigned long long int ielem)
{
    return _container.at(ielem);
}

// Immutable access to an element with range-check
/// \brief          Immutable access to an element with range-check.
/// \details        Provides an access to the specified element with a 
///                 range-check.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::at(const unsigned long long int ielem) const
{
    return _container.at(ielem);
}

// Access to the i-th element from the beginning
/// \brief          Access to the i-th element from the beginning.
/// \details        Returns a reference to the i-th first element in the 
///                 underlying container without doing any range check.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::front(const unsigned long long int ielem)
{
    return *(std::begin(_container)+ielem);
}

// Immutable access to the i-th element from the beginning
/// \brief          Immutable access to the i-th element from the beginning.
/// \details        Returns a reference to the i-th first element in the 
///                 underlying container without doing any range check.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::front(const unsigned long long int ielem) const 
{
    return *(std::begin(_container)+ielem);
}

// Access to the i-th element from the end
/// \brief          Access to the i-th element from the end.
/// \details        Returns a reference to the i-th last element in the 
///                 underlying container without doing any range check.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::back(const unsigned long long int ielem)
{
    return *(std::begin(_container)+std::distance(std::begin(_container), std::end(_container))-1-static_cast<long long int>(ielem));
}

// Immutable access to the i-th element from the end
/// \brief          Immutable access to the i-th element from the end.
/// \details        Returns a reference to the i-th last element in the 
///                 underlying container without doing any range check.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::back(const unsigned long long int ielem) const
{
    return *(std::begin(_container)+std::distance(std::begin(_container), std::end(_container))-1-static_cast<long long int>(ielem));
}

// Cyclic access to an element
/// \brief          Cyclic access to an element.
/// \details        Provides a cyclic access to the elements, using the index
///                 modulo. Negative indexes are supported. It allows to iterate
///                 several times over the contents just by incrementing the
///                 provided index.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cycle(const long long int ielem)
{
    const long long int length = std::distance(std::begin(_container), std::end(_container));
    return _container[(ielem >= 0) ? (ielem%length) : ((length-((-ielem)%length))%length)];
}

// Immutable cyclic access to an element
/// \brief          Immutable cyclic access to an element.
/// \details        Provides a cyclic access to the elements, using the index
///                 modulo. Negative indexes are supported. It allows to iterate
///                 several times over the contents just by incrementing the
///                 provided index.
/// \param[in]      ielem Index of the element in the underlying container.
/// \return         Immutable reference to the element.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cycle(const long long int ielem) const
{
    const long long int length = std::distance(std::begin(_container), std::end(_container));
    return _container[(ielem >= 0) ? (ielem%length) : ((length-((-ielem)%length))%length)];
}

// Direct access to the underlying container 
/// \brief          Direct access to the underlying container.
/// \details        Provides a direct access to the underlying container by 
///                 returning a reference to it.
/// \return         Reference to the underlying container.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Container& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::container()
{
    return _container;
}

// Immutable direct access to the underlying container 
/// \brief          Direct access to the underlying container.
/// \details        Provides a direct access to the underlying container by 
///                 returning a reference to it.
/// \return         Immutable reference to the underlying container.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Container& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::container() const
{
    return _container;
}

// Direct access to the underlying array 
/// \brief          Direct access to the underlying array.
/// \details        Provides a direct access to the underlying array by 
///                 returning a pointer to the first element of storage.
/// \return         Pointer to the underlying element storage.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline Element* SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::data()
{
    return &(*std::begin(_container));
}

// Immutable direct access to the underlying array 
/// \brief          Immutable direct access to the underlying array.
/// \details        Provides a direct access to the underlying array by 
///                 returning a pointer to the first element of storage.
/// \return         Immutable pointer to the underlying element storage.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline const Element* SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::data() const
{
    return &(*std::begin(_container));
}


//--------------------------------------------------------------------------- //



// ------------------------------- ITERATORS -------------------------------- //
// Iterator to the beginning
/// \brief          Iterator to the beginning.
/// \details        Returns a pointer to the first element.
/// \tparam         Iterator (Iterator type.)
/// \return         Pointer to the beginning.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::begin()
{
    return _container.begin();
}

// Immutable iterator to the beginning
/// \brief          Immutable iterator to the beginning.
/// \details        Returns a pointer to the first element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the beginning.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::begin() const
{
    return _container.begin();
}

// Forced immutable iterator to the beginning
/// \brief          Forced immutable iterator to the beginning.
/// \details        Returns a constant pointer to the first element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the beginning.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cbegin() const
{
    return _container.cbegin();
}

// Iterator to the end
/// \brief          Iterator to the end.
/// \details        Returns a pointer to the position after the last element.
/// \tparam         Iterator (Iterator type.)
/// \return         Pointer to the end.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::end()
{
    return _container.end();
}

// Immutable iterator to the end
/// \brief          Immutable iterator to the end.
/// \details        Returns a pointer to the position after the last element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the end.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::end() const
{
    return _container.end();
}

// Forced immutable iterator to the end
/// \brief          Forced immutable iterator to the end.
/// \details        Returns a constant pointer to the position after the last 
///                 element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the end.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cend() const
{
    return _container.cend();
}

// Reverse iterator to the beginning
/// \brief          Reverse iterator to the beginning.
/// \details        Returns a reversed pointer to the position after the last 
///                 element.
/// \tparam         Iterator (Iterator type.)
/// \return         Pointer to the end.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::rbegin()
{
    return _container.rbegin();
}

// Immutable reverse iterator to the beginning
/// \brief          Immutable reverse iterator to the beginning.
/// \details        Returns a reversed pointer to the position after the last 
///                 element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the end.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::rbegin() const
{
    return _container.rbegin();
}

// Forced immutable reverse iterator to the beginning
/// \brief          Forced immutable reverse iterator to the beginning.
/// \details        Returns a constant reversed pointer to the position after 
///                 the last element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the end.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::crbegin() const
{
    return _container.crbegin();
}

// Reverse iterator to the end
/// \brief          Reverse iterator to the end.
/// \details        Returns a reversed pointer to the first element.
/// \tparam         Iterator (Iterator type.)
/// \return         Pointer to the beginning.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::rend()
{
    return _container.rend();
}

// Immutable reverse iterator to the end
/// \brief          Immutable reverse iterator to the end.
/// \details        Returns a reversed pointer to the first element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the beginning.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::rend() const
{
    return _container.rend();
}

// Forced immutable reverse iterator to the end
/// \brief          Forced immutable reverse iterator to the end.
/// \details        Returns a constant reversed pointer to the first element.
/// \tparam         Iterator (Iterator type.)
/// \return         Immutable pointer to the beginning.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::crend() const
{
    return _container.crend();
}

// Index of an iterator in the underlying container
/// \brief          Index of an iterator in the underlying container.
/// \details        Returns the index of the element pointed by an iterator.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      it Iterator to an element.
/// \return         Index of the element in the underlying container.
/// \exception      std::out_of_range Out of range.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline unsigned long long int SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::index(const Iterator& it) const
{
    const decltype(std::begin(_container)) iter(&*it);
    if ((iter < std::begin(_container)) || (iter >= std::end(_container))) {
        throw std::out_of_range("ERROR = SimpleHyperOctree::index() : out of range");
    }
    return std::distance(std::begin(_container), iter);
}
// -------------------------------------------------------------------------- //



// --------------------------------- SEARCH --------------------------------- //
// Find element from index
/// \name           Find element from index.
/// \details        Finds the most refined cell that contains the provided 
///                 index and returns an iterator to it.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      idx Index.
/// \return         Iterator to the element containing the specified index.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class> 
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::find(const Index& idx)
{
    const Element elem = Element(idx, Data());
    const unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    return (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::begin(_container)+marker-(marker > 0)) : (std::end(_container));
}

// Find immutable element from index
/// \name           Find immutable element from index.
/// \details        Finds the most refined cell that contains the provided 
///                 index and returns an immutable iterator to it.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      idx Index.
/// \return         Immutable iterator to the element containing the specified 
///                 index.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class> 
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::find(const Index& idx) const
{
    const Element elem = Element(idx, Data());
    const unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    return (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::begin(_container)+marker-(marker > 0)) : (std::end(_container));
}

// Locate element from position
/// \name           Locate element from position.
/// \details        Locates the most refined cell at the provided position and
///                 returns an iterator to it.
/// \tparam         Iterator (Iterator type.)
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Iterator to the element found at the specified position.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, typename... Types, class> 
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::locate(Types&&... iposs)
{
    const Element elem = Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data());
    const unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    return (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::begin(_container)+marker-(marker > 0)) : (std::end(_container));
}

// Locate immutable element from position
/// \name           Locate immutable element from position.
/// \details        Locates the most refined cell at the provided position and
///                 returns an immutable iterator to it.
/// \tparam         Iterator (Iterator type.)
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Immutable iterator to the element found at the specified 
///                 position.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, typename... Types, class> 
inline Iterator SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::locate(Types&&... iposs) const
{
    const Element elem = Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data());
    const unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    return (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::begin(_container)+marker-(marker > 0)) : (std::end(_container));
}
// -------------------------------------------------------------------------- //



//--------------------------------- CAPACITY -------------------------------- //
// Emptiness checking
/// \brief          Emptiness checking.
/// \details        Checks if the hyperoctree has no elements.
/// \return         True if empty, false otherwise.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline bool SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::empty() const
{
    return (std::begin(_container) == std::end(_container));
}

// Number of elements 
/// \brief          Number of elements.
/// \details        Returns the distance between the first and the last element.
/// \return         The number of elements in the hyperoctree.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline unsigned long long int SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::size() const
{
    return std::distance(std::begin(_container), std::end(_container));
}

// Capacity of the underlying storage
/// \brief          Capacity of the underlying storage.
/// \details        Returns the number of elements that the hyperoctree has 
///                 currently allocated space for. 
/// \return         Capacity of the currently allocated underlying storage. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline unsigned long long int SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::capacity() const
{
    return _container.capacity();
}

// Available space
/// \brief          Available space.
/// \details        Returns the maximum possible number of elements.
/// \return         Maximum number of elements. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline unsigned long long int SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::space() const
{
    return _container.max_size();
}

// Shrink reserved storage
/// \brief          Shrink reserved storage
/// \details        Reduces memory usage by freeing unused memory.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::shrink()
{
    _container.shrink_to_fit();
    return *this;
}

// Increase reserved storage
/// \brief          Increase reserved storage.
/// \details        Increases the capacity of the underlying storage. Existing
///                 elements are protected so it could not invalidate the 
///                 actual contents.
/// \param[in]      nelem New reserved number of elements.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::reserve(const unsigned long long int nelem)
{
    _container.reserve(nelem);
    return *this;
}
//--------------------------------------------------------------------------- //



//--------------------------------- MODIFIERS ------------------------------- //
// Clear the contents
/// \brief          Full clear the content.
/// \details        Erase the hyperoctree.
/// \return         Self reference. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::fullclear()
{
    Container().swap(_container);
    return *this;
}

// Clear the contents
/// \brief          Clear the contents.
/// \details        Removes all elements from the hyperoctree.
/// \return         Self reference. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::clear()
{
    _container.clear();
    return *this;
}

// Pop back
/// \brief          Pop back.
/// \details        Removes the last element of the hyperoctree.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::pop()
{
    _container.pop_back();
    return *this;
}

// Append an element
/// \brief          Append an element.
/// \details        Appends an element by constructing it from the arguments 
///                 and emplacing it at the end of the hyperoctree.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class... Misc> 
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::append(Misc&&... misc)
{
    _container.emplace_back(std::forward<Misc>(misc)...);
    return *this;
}

// Resize
/// \brief          Resize.
/// \details        Resizes the hyperoctree to contain a new number of elements.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class... Misc> 
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::resize(Misc&&... misc)
{
    _container.resize(std::forward<Misc>(misc)...);
    return *this;
}

// Access to the rho data
/// \brief          Access to the rho data. 
/// \details        Provides an access to the rho data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::rho(Misc&&... misc)
{
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the rho data
/// \brief          Access to the rho data. 
/// \details        Provides an access to the rho data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::rho(Misc&&... misc) const
{
    return this->template data<0, Values...>(std::forward<Misc>(misc)...);
}

// Access to the phi data
/// \brief          Access to the phi data. 
/// \details        Provides an access to the phi data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::phi(Misc&&... misc)
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the phi data
/// \brief          Immutable access to the phi data. 
/// \details        Provides an access to the phi data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::phi(Misc&&... misc) const
{
    return this->template data<1, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dphidxyz data
/// \brief          Access to the dphidxyz data. 
/// \details        Provides an access to the dphidxyz data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidxyz(Misc&&... misc)
{
    return this->template data<2, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dphidxyz data
/// \brief          Immutable access to the dphidxyz data. 
/// \details        Provides an access to the dphidxyz data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidxyz(Misc&&... misc) const
{
    return this->template data<2, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dphidx data
/// \brief          Access to the dphidx data. 
/// \details        Provides an access to the dphidx data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidx(Misc&&... misc)
{
    return this->template data<2, Values...>(0 ,std::forward<Misc>(misc)...);
}

// Immutable access to the dphidx data
/// \brief          Immutable access to the dphidx data. 
/// \details        Provides an access to the dphidx data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidx(Misc&&... misc) const
{
    return this->template data<2, Values...>(0 ,std::forward<Misc>(misc)...);
}

// Access to the dphidy data
/// \brief          Access to the dphidy data. 
/// \details        Provides an access to the dphidy data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidy(Misc&&... misc)
{
    return this->template data<2, Values...>(1 ,std::forward<Misc>(misc)...);
}

// Immutable access to the dphidy data
/// \brief          Immutable access to the dphidy data. 
/// \details        Provides an access to the dphidy data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidy(Misc&&... misc) const
{
    return this->template data<2, Values...>(1 ,std::forward<Misc>(misc)...);
}

// Access to the dphidz data
/// \brief          Access to the dphidz data. 
/// \details        Provides an access to the dphidz data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidz(Misc&&... misc)
{
    return this->template data<2, Values...>(2 ,std::forward<Misc>(misc)...);
}

// Immutable access to the dphidz data
/// \brief          Immutable access to the dphidz data. 
/// \details        Provides an access to the dphidz data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dphidz(Misc&&... misc) const
{
    return this->template data<2, Values...>(2 ,std::forward<Misc>(misc)...);
}

// Access to the a data
/// \brief          Access to the a data. 
/// \details        Provides an access to the a data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::a(Misc&&... misc)
{
    return this->template data<3, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the a data
/// \brief          Immutable access to the a data. 
/// \details        Provides an access to the a data by forwarding 
///                 parameters to the unified base accessor member. 
/// \tparam         Values List of template values. 
/// \tparam         Misc (Miscellaneous types.) 
/// \tparam         Template (Deduced template type.) 
/// \param[in]      misc Miscellaneous arguments. 
/// \return         Forwarded result.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int... Values,class... Misc,class Template, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::a(Misc&&... misc) const
{
    return this->template data<3, Values...>(std::forward<Misc>(misc)...);
}
//--------------------------------------------------------------------------- //



//-------------------------------- REFINEMENT ------------------------------- 
// Update refinement
/// \brief          Update refinement.
/// \details        Updates the octree refinement by removing coarsened cells
///                 , by sorting refined ones and by removing cells with the
///                 same index.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::update()
{ 
    _container.resize(std::distance(std::begin(_container), std::remove_if(std::begin(_container), std::end(_container), [](const Element& elem){return std::get<0>(elem).invalidated();})));
#ifdef MEMORY_SAVING
    std::sort(std::begin(_container), std::end(_container), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);});
#else
    unsigned long long int marker = std::distance(std::begin(_container), std::is_sorted_until(std::begin(_container), std::end(_container), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    if (marker < static_cast<unsigned long long int>(std::distance(std::begin(_container), std::end(_container)))) {
        std::sort(std::begin(_container)+marker, std::end(_container), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);});
        std::inplace_merge(std::begin(_container), std::begin(_container)+marker, std::end(_container), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);});
    }
#endif     

    unsigned int size = _container.size();
    std::vector<unsigned int> count(size);
    Utility::parallelize(size, [=, &count](const unsigned int i){
	auto j = std::adjacent_find(std::begin(_container)+i, std::end(_container), [](const Element& first, const Element& second){return std::get<0>(first) != std::get<0>(second);});
	count[i] = std::distance(std::begin(_container)+i, j);
    });

    Utility::parallelize(size, [=](const unsigned int i){
	if(count[i]){
	    if(!((i > 0) && count[i-1])){ // Previous index must not be in a serie of identical indices.
	        std::vector<Data> tmpvector;
		for(unsigned int k = 0; k < count[i]+1; ++k){
		    tmpvector.push_back(std::get<1>(_container.data()[i+k]));
		}
	 	std::sort(std::begin(tmpvector), std::end(tmpvector), [](const Data& first, const Data& second){return first.a() < second.a();});
		std::get<1>(_container.data()[i]) = tmpvector[0];
		if(std::abs(tmpvector.front().a() - tmpvector.back().a()) > 1e-6){ // To avoid numerical noise
		    std::get<1>(_container.data()[i]).dphidt() = (tmpvector.front().phi()-tmpvector.back().phi())/(tmpvector.front().a()-tmpvector.back().a());
		}
	    }
	}
    });

    _container.resize(std::distance(std::begin(_container), std::unique(std::begin(_container), std::end(_container), [](const Element& first, const Element& second){	return std::get<0>(first) == std::get<0>(second);})));
#ifndef MEMORY_SAVING
    _container.shrink_to_fit();
#endif
    return *this;
}

// Root level
/// \brief          Root level.
/// \details        Checks whether the given element corresponds to the root 
///                 level.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      it Iterator to an element.
/// \return         True if the element is at the root level, false otherwise.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class> 
inline bool SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::root(const Iterator& it)
{
    const Index idx = std::get<0>(*it);
    return std::none_of(std::begin(_container), decltype(std::begin(_container))(it), [=, &idx](const Element& elem){return std::get<0>(elem).containing(idx);});
}

// Leaf level
/// \brief          Leaf level.
/// \details        Checks whether the given element corresponds to a leaf 
///                 level. If the element does not have any child, it is 
///                 considered at the leaf level.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      it Iterator to an element.
/// \return         True if the element is not refined, false otherwise.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class> 
inline bool SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::leaf(const Iterator& it)
{
    return (decltype(std::begin(_container))(it)+1 < std::end(_container)) ? (!(std::get<0>(*it).containing(std::get<0>(*(decltype(std::begin(_container))(it)+1))))) : (true);
}

// Refine an element
/// \brief          Refine an element.
/// \details        Refines a given element by inserting the related children 
///                 at the end of the container. To become valid and accessible,
///                 the container should be updated.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      it Iterator to an element.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::refine(const Iterator& it)
{
    const Index idx = std::get<0>(*it);
    const Data dat = std::get<1>(*it);
    for (unsigned int isite = 0; isite < Index::sites(); ++isite) {
        _container.emplace_back(idx.child(isite), dat);
    }
    return *this;
}

// Coarsen an element
/// \brief          Coarsen an element.
/// \details        Coarsens a given element by invalidating all the underlying 
///                 children. To become valid again, the container should be 
///                 updated.
/// \tparam         Iterator (Iterator type.)
/// \param[in]      it Iterator to an element.
/// \return         Self reference.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename Iterator, class>
inline SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::coarsen(const Iterator& it)
{
    const decltype(std::begin(_container)) iter(&*it);
    const Index idx = std::get<0>(*iter);
    std::for_each(iter+1, std::find_if_not(iter+1, std::end(_container), [=, &idx](const Element& elem){return idx.containing(std::get<0>(elem));}), [](Element& elem){std::get<0>(elem).invalidate();}); 
    return *this;
}
//--------------------------------------------------------------------------- //



// ------------------------------ INTERPOLATION ----------------------------- //
// Nearest grid point interpolation
/// \brief          Nearest grid point interpolation.
/// \details        Computes the value of the data at the provided position 
///                 using a nearest grid point interpolation scheme.
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Value of the data at the provided position.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename... Types, class, class> 
inline Data SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::ngp(Types&&... iposs) const
{
    const Element elem = Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data());
    const unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    return (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::get<1>(*(std::begin(_container)+marker-1))) : (Data());
}

// Cloud-in-cell interpolation
/// \brief          Cloud in cell interpolation.
/// \details        Computes the value of the data at the provided position 
///                 using a cloud in cell interpolation scheme.
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Value of the data at the provided position.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename... Types, class, class> 
inline Data SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cic(Types&&... iposs) const
{
    static const Type one = Type(1);
    static const Type two = Type(2);
    const std::array<Type, sizeof...(iposs)> point = tupleify<std::array<Type, sizeof...(iposs)> >(std::forward<Types>(iposs)...);
    Element elem = Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data());
    unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    Index idx = (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::get<0>(*(std::begin(_container)+marker-1))) : (Index::invalid());
    unsigned int ilvl = idx.level();
    unsigned int lvl = ilvl;
    Type half = Type(Extent::num)/Type(Extent::den)*idx.extent()/two;  
    std::array<Type, sizeof...(iposs)> vertex = point; 
    std::array<Type, sizeof...(iposs)> dist = point; 
    std::bitset<sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits> site;
    std::bitset<sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits> bits;
    std::bitset<sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits> diff;
    Type weight = one;
    Data result = std::get<1>(elem);
    do {
        ilvl = lvl;
        if (!idx.invalidated()) {
            for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                dist[idim] = std::abs((idx.template center<Type, Position, Extent>(idim)-point[idim])/(half+half));
            }
            for (unsigned int isite = 0; isite < Index::sites(); ++isite) {
                if (idx.child(isite).containing(std::get<0>(elem))) {
                    site = isite;
                    isite = Index::sites();
                }
            }
            for (unsigned int ivertex = 0; ivertex < Index::sites(); ++ivertex) {
                weight = one;
                bits = ivertex;
                diff = site ^ bits; 
                for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                    vertex[idim] = point[idim]+((bits[idim])*half)-((!bits[idim])*half);
                    weight *= ((diff[idim])*(one-dist[idim]))+((!diff[idim])*dist[idim]);
                }
                std::get<0>(elem) = variadify<Index>([=, &ilvl](const typename std::remove_cv<typename std::remove_reference<Types>::type>::type&... input){return Index::template compute<Type, Position, Extent>(ilvl, input...);}, vertex); 
                marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
                if (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) {
                    lvl = std::min(lvl, std::get<0>(*(std::begin(_container)+marker-1)).level());
                    mac(result, std::get<1>(*(std::begin(_container)+marker-1)), weight);
                } else {
                    result = Data();
                    ivertex = Index::sites();
                }
            }
        } 
        if (lvl < ilvl) {
	    elem = Element(Index::template compute<Type, Position, Extent>(lvl+1, std::forward<Types>(iposs)...), Data());
	    idx = std::get<0>(elem).parent();
            result = std::get<1>(elem);
	    half = Type(Extent::num)/Type(Extent::den)*idx.extent()/two;
        }
    } while (lvl < ilvl);
    return result;
}


// Triangular Shaped Cloud interpolation
/// \brief          Triangular Shaped Cloud interpolation.
/// \details        Computes the value of the data at the provided position 
///                 using a triangular shaped cloud interpolation scheme.
/// \tparam         Types (Scalar position types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Value of the data at the provided position.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename... Types, class, class> 
inline Data SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::tsc(Types&&... iposs) const
{
    static const Type one = Type(1);
    static const Type two = Type(2);
    const std::array<Type, sizeof...(iposs)> point = tupleify<std::array<Type, sizeof...(iposs)> >(std::forward<Types>(iposs)...);
    Element elem = Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data());
    unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    Index idx = (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::get<0>(*(std::begin(_container)+marker-1))) : (Index::invalid());
    unsigned int ilvl = idx.level();
    unsigned int lvl = ilvl;
    Type half = Type(Extent::num)/Type(Extent::den)*idx.extent()/two;  
    std::array<Type, sizeof...(iposs)> dist = point; 
    Data result = std::get<1>(elem);
    Type weight = 0;
    if(point == std::array<Type, sizeof...(iposs)>()){ // if position is the center of the box, perform CIC
	return SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cic(0.,0.,0.);
    } else{
        do {
            ilvl = lvl;
            if (!idx.invalidated()) {
	        const double twohalves = two*half;
	        for (int ix = -one; ix <= one; ++ix) {
		    int aix = abs(ix);
		    for (int iy = -one; iy <= one; ++iy) {
		        int aiy = abs(iy);
		        for (int iz = -one; iz <= one; ++iz) {
		   	    int aiz = abs(iz);
	 	            std::get<0>(elem) = Index::template compute<Type, Position, Extent>(ilvl, point[0]+ix*twohalves, point[1]+iy*twohalves, point[2]+iz*twohalves);
			    marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
                    	    if (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) {
			        for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                    	            dist[idim] = std::abs((std::get<0>(elem).template center<Type, Position, Extent>(idim)-point[idim])/(twohalves));
			        }
				weight = (aix*0.5*(1.5-dist[0])*(1.5-dist[0]) + (1-aix)*(0.75-dist[0]*dist[0]))
					*(aiy*0.5*(1.5-dist[1])*(1.5-dist[1]) + (1-aiy)*(0.75-dist[1]*dist[1]))
					*(aiz*0.5*(1.5-dist[2])*(1.5-dist[2]) + (1-aiz)*(0.75-dist[2]*dist[2]));
                        	lvl = std::min(lvl, std::get<0>(*(std::begin(_container)+marker-1)).level());
				if(lvl < ilvl){
				    result = Data();
				    ix = iy = iz = one;
				} else{
                        	    mac(result, std::get<1>(*(std::begin(_container)+marker-1)), weight);
				}
                    	    } else {
                        	result = Data();
                        	ix = iy = iz = one;
                    	    }
		        }
		    }
	        }	
            } 
            if (lvl < ilvl) {
	        elem = Element(Index::template compute<Type, Position, Extent>(lvl+1, std::forward<Types>(iposs)...), Data());
	        idx = std::get<0>(elem).parent();
                result = std::get<1>(elem);
	        half = Type(Extent::num)/Type(Extent::den)*idx.extent()/two;
            }
        } while (lvl < ilvl);
    	return result;
    }
}

// Triangular Shaped Cloud interpolation
/// \brief          Triangular Shaped Cloud interpolation.
/// \details        Computes the value of the data at the provided position 
///                 using a triangular shaped cloud interpolation scheme. 
///		    Also use previsouly computed neighbouring cells
/// \tparam         Types (Scalar position types.)
/// \param[in]      elemsTsc Neighbouring cells previsouly computed
/// \param[in]      iposs Real positions along each dimension.
/// \return         Value of the data at the provided position.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <typename... Types, class, class> 
inline Data SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::tsc(std::vector<Element>& elemsTsc, Types&&... iposs) const
{
    static const Type one = Type(1);
    static const Type two = Type(2);
    const std::array<Type, sizeof...(iposs)> point = tupleify<std::array<Type, sizeof...(iposs)> >(std::forward<Types>(iposs)...);
    Element elem = Element(Index::template cipher<Type, Position, Extent>(std::forward<Types>(iposs)...), Data());
    unsigned long long int marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
    Index idx = (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) ? (std::get<0>(*(std::begin(_container)+marker-1))) : (Index::invalid());
    unsigned int ilvl = idx.level();
    unsigned int lvl = ilvl;
    Type half = Type(Extent::num)/Type(Extent::den)*idx.extent()/two;  
    std::array<Type, sizeof...(iposs)> dist = point; 
    Data result = std::get<1>(elem);
    Type weight = 0;
    std::vector<Element> elemsTsctmp(27);
    unsigned int ic = 0;
    if(point == std::array<Type, sizeof...(iposs)>()){ // if position is the center of the box, perform CIC
	return SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::cic(0.,0.,0.);
    } else{
        do {
            ilvl = lvl;
	    const double twohalves = two*half;
	    if(std::get<0>(elemsTsc[13]) == idx){
		ic = 0;
	        for (int ix = -one; ix <= one; ++ix) {
		    int aix = abs(ix);
		    for (int iy = -one; iy <= one; ++iy) {
		        int aiy = abs(iy);
		        for (int iz = -one; iz <= one; ++iz) {
		            int aiz = abs(iz);
		            for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                        	dist[idim] = std::abs((std::get<0>(elemsTsc[ic]).template center<Type, Position, Extent>(idim)-point[idim])/(twohalves));
		    	    }
		    	    weight =  (aix*0.5*(1.5-dist[0])*(1.5-dist[0]) + (1-aix)*(0.75-dist[0]*dist[0]))
				     *(aiy*0.5*(1.5-dist[1])*(1.5-dist[1]) + (1-aiy)*(0.75-dist[1]*dist[1]))
				     *(aiz*0.5*(1.5-dist[2])*(1.5-dist[2]) + (1-aiz)*(0.75-dist[2]*dist[2]));  
                    	    mac(result, std::get<1>(elemsTsc[ic]), weight);
			    ic++;
			}
		    }
		}
		return result;
	    } else{
                if (!idx.invalidated()) {
		    ic = 0;
	            for (int ix = -one; ix <= one; ++ix) {
		        int aix = abs(ix);
		        for (int iy = -one; iy <= one; ++iy) {
		            int aiy = abs(iy);
		            for (int iz = -one; iz <= one; ++iz) {
		   	        int aiz = abs(iz);
	 	                std::get<0>(elem) = Index::template compute<Type, Position, Extent>(ilvl, point[0]+ix*twohalves, point[1]+iy*twohalves, point[2]+iz*twohalves);
			        marker = std::distance(std::begin(_container), std::upper_bound(std::begin(_container), std::end(_container), elem, [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
                    	        if (std::get<0>(*(std::begin(_container)+marker-(marker > 0))).containing(std::get<0>(elem))) {
				    elemsTsctmp[ic] = *(std::begin(_container)+marker-(marker > 0));
			            for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                    	                dist[idim] = std::abs((std::get<0>(elem).template center<Type, Position, Extent>(idim)-point[idim])/(twohalves));
			            }
		    		    weight =  (aix*0.5*(1.5-dist[0])*(1.5-dist[0]) + (1-aix)*(0.75-dist[0]*dist[0]))
					     *(aiy*0.5*(1.5-dist[1])*(1.5-dist[1]) + (1-aiy)*(0.75-dist[1]*dist[1]))
					     *(aiz*0.5*(1.5-dist[2])*(1.5-dist[2]) + (1-aiz)*(0.75-dist[2]*dist[2])); 
                        	    lvl = std::min(lvl, std::get<0>(elemsTsctmp[ic]).level());
				    if(lvl < ilvl){
				        result = Data();
				        ix = iy = iz = one;
				    } else{
                        	        mac(result, std::get<1>(elemsTsctmp[ic]), weight);
					ic++;
				    }
                    	        } else {
                        	    result = Data();
                        	    ix = iy = iz = one;
                    	        }
		            }
		        }
	            }	
		    if(27 == ic){
			elemsTsc = elemsTsctmp;
		    }
                }
	    } 
            if (lvl < ilvl) {
	        elem = Element(Index::template compute<Type, Position, Extent>(lvl+1, std::forward<Types>(iposs)...), Data());
	        idx = std::get<0>(elem).parent();
                result = std::get<1>(elem);
	        half = Type(Extent::num)/Type(Extent::den)*idx.extent()/two;
            }		
        } while (lvl < ilvl);
    	return result;
    }
}

//--------------------------------------------------------------------------- //



//---------------------------------- STREAM --------------------------------- //
// Output stream operator 
/// \brief          Output stream operator. 
/// \details        Displays the whole structure of the octree.
/// \tparam         SelfType (Scalar position type.)
/// \tparam         SelfIndex (Index type.)
/// \tparam         SelfData (Data type.)
/// \tparam         SelfDimension (Number of dimensions.)
/// \tparam         SelfPosition (Position of the hyperoctree center.)
/// \tparam         SelfExtent (Extent of the hyperoctree.)
/// \tparam         SelfElement (Underlying element type.)
/// \tparam         SelfContainer (Underlying container type.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side object.
/// \return         Output stream. 
template <typename SelfType, class SelfIndex, class SelfData, unsigned int SelfDimension, class SelfPosition, class SelfExtent, class SelfElement, class SelfContainer>
std::ostream& operator<<(std::ostream& lhs, const SimpleHyperOctree<SelfType, SelfIndex, SelfData, SelfDimension, SelfPosition, SelfExtent, SelfElement, SelfContainer>& rhs)
{
    static const unsigned int iwidth = std::ceil((sizeof(SelfIndex::type())*std::numeric_limits<unsigned char>::digits)*std::log10(2));
    static const unsigned int lwidth = std::to_string(SelfIndex::instructions()).size();
    static const unsigned int fwidth = std::numeric_limits<SelfType>::max_digits10;
    static const unsigned int twidth = std::numeric_limits<SelfType>::max_digits10+std::numeric_limits<unsigned char>::digits;
    static const unsigned int ten = 10;
    static const char czero = '0';
    const unsigned long long int length = std::distance(std::begin(rhs._container), std::end(rhs._container));
    const char separator = lhs.fill();
    const std::string separators(std::numeric_limits<unsigned char>::digits, separator);
    const std::ios::fmtflags flags = lhs.flags();
    SelfIndex idx = SelfIndex();
    for (unsigned long long int ielem = 0; ielem < length; ++ielem) {
        idx = std::get<0>(*(std::begin(rhs._container)+ielem));
        lhs<<std::setw(iwidth)<<ielem<<separator;
        lhs<<std::setfill(czero)<<std::setw(iwidth)<<std::dec<<idx.stringify(ten)<<separator;
        lhs<<std::setfill(separator)<<std::setw(lwidth)<<std::dec<<idx.level()<<separators;
        lhs.flags(flags);
        lhs<<std::setfill(separator)<<std::scientific<<std::setprecision(fwidth)<<std::right;
        for (unsigned int idim = 0; idim < SelfIndex::dimension(); ++idim) {
            lhs<<std::setw(twidth)<<idx.template center<SelfType, SelfPosition, SelfExtent>(idim)<<separator;
        }
        lhs<<std::setw(twidth)<<idx.template length<SelfType, SelfPosition, SelfExtent>()<<separators;
        lhs.flags(flags);
        lhs<<std::get<1>(*(std::begin(rhs._container)+ielem))<<std::endl;
    }
    lhs<<std::setfill(separator);
    lhs.flags(flags);
    return lhs;
}
//--------------------------------------------------------------------------- //



//-------------------------------- PROPERTIES ------------------------------- //
// Type
/// \brief          Type. 
/// \details        Returns a copy of the default type value.
/// \return         Copy of the default type value.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
constexpr Type SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::type()
{
    return Type();
}

// Position
/// \brief          Position. 
/// \details        Returns a copy of the default position value.
/// \return         Copy of the default position value.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
constexpr Position SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::position()
{
    return Position();
}

// Extent
/// \brief          Extent. 
/// \details        Returns a copy of the default extent value.
/// \return         Copy of the default extent value.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
constexpr Extent SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::extent()
{
    return Extent();
}

// Element
/// \brief          Element. 
/// \details        Returns a copy of the default element value.
/// \return         Copy of the default element value.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
constexpr Element SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::element()
{
    return Element();
}

// Number of dimensions
/// \brief          Number of dimensions. 
/// \details        Returns the number of spatial dimensions.
/// \return         Copy of the number of dimensions.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
constexpr unsigned int SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::dimension()
{
    return Dimension;
}
//--------------------------------------------------------------------------- //



//--------------------------------- HELPERS --------------------------------- //
// Convert a tuple to another tuple by extraction
/// \brief          Convert a tuple to another tuple by extraction.
/// \details        Converts the provided tuple to another one by recursively
///                 extracting the underlying values.
/// \tparam         Template Output tuple type.
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in]      tuple Input tuple.
/// \param[in]      values (Extracted values.)
/// \return         Tuple constructed from input.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, class Tuple, typename... Kinds, class> 
constexpr Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::tupleify(Tuple&& tuple, Kinds&&... values)
{
    return tupleify<Template>(std::forward<Tuple>(tuple), std::forward<Kinds>(values)..., std::get<sizeof...(Kinds)>(std::forward<Tuple>(tuple)));
}

// Convert a tuple to another tuple
/// \brief          Convert a tuple to another tuple.
/// \details        Converts the provided tuple to another one.
/// \tparam         Template Output tuple type.
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in]      tuple Input tuple.
/// \param[in]      values (Extracted values.)
/// \return         Tuple constructed from input.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, class Tuple, typename... Kinds, class> 
constexpr Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::tupleify(const Tuple& tuple, Kinds&&... values)
{
    return tupleify<typename std::conditional<static_cast<bool>(sizeof(tuple)), Template, Template>::type>(std::forward<Kinds>(values)...);
}

// Convert values to a tuple
/// \brief          Convert values to a tuple.
/// \details        Converts the provided values to a tuple.
/// \tparam         Template Output tuple type.
/// \tparam         Kinds (Arithmetic types.)
/// \param[in]      values Input values.
/// \return         Tuple constructed from input.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, typename... Kinds, class>
constexpr Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::tupleify(Kinds&&... values)
{
    return Template{std::forward<Kinds>(values)...};
}

// Forward a tuple to a variadic function by extraction
/// \brief          Forward a tuple to a variadic function by extraction.
/// \details        Forwards the underlying values of the provided tuple to a 
///                 variadic function by recursively extracting them.
/// \tparam         Template Output type.
/// \tparam         Function (Function type.)
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in,out]  f Variadic function object.
/// \param[in,out]  tuple Input tuple.
/// \param[in,out]  values (Extracted values.)
/// \return         Function result converted to the output type.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, class Function, class Tuple, typename... Kinds, class> 
constexpr Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::variadify(Function&& f, Tuple&& tuple, Kinds&&... values)
{
    return variadify<Template>(std::forward<Function>(f), std::forward<Tuple>(tuple), std::forward<Kinds>(values)..., std::get<sizeof...(Kinds)>(std::forward<Tuple>(tuple)));
}

// Forward a tuple to a variadic function
/// \brief          Forward a tuple to a variadic function.
/// \details        Forwards the underlying values of the provided tuple to a 
///                 variadic function.
/// \tparam         Template Output type.
/// \tparam         Function (Function type.)
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in,out]  f Variadic function object.
/// \param[in,out]  tuple Input tuple.
/// \param[in,out]  values (Extracted values.)
/// \return         Function result converted to the output type.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, class Function, class Tuple, typename... Kinds, class> 
constexpr Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::variadify(Function&& f, const Tuple& tuple, Kinds&&... values)
{
    return variadify<typename std::conditional<static_cast<bool>(sizeof(tuple)), Template, Template>::type>(std::forward<Function>(f), std::forward<Kinds>(values)...);
}

// Forward values to a variadic function
/// \brief          Forward values to a variadic function.
/// \details        Forwards the values to a variadic function.
/// \tparam         Template Output type.
/// \tparam         Function (Function type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in,out]  f Variadic function object.
/// \param[in,out]  values Input values.
/// \return         Function result converted to the output type.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Template, class Function, typename... Kinds, class> 
constexpr Template SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::variadify(Function&& f, Kinds&&... values)
{
    return Template(std::forward<Function>(f)(std::forward<Kinds>(values)...));
}

// Multiply-accumulate objects
/// \brief          Multiply-accumulate objects.
/// \details        Executes a multiply-accumulate operation on each element of
///                 the provided object and returns the accumulator.
/// \tparam         Step (Recursion step.)
/// \tparam         Object (Object type with a tuple data member.)
/// \tparam         Kind (Arithmetic type.)
/// \param[in,out]  accumulator Accumulator object reference.
/// \param[in]      factor Factor object value.
/// \param[in]      coefficient Multiplier value.
/// \return         Reference to the accumulator object.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int Step, class Object, typename Kind, class> 
constexpr Object& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::mac(Object& accumulator, const Object& factor, const Kind coefficient)
{
    return (mac(accumulator.data(), factor.data(), coefficient*(sizeof(Step)/sizeof(Step))), accumulator);
}

// Multiply-accumulate tuples
/// \brief          Multiply-accumulate tuples.
/// \details        Executes a multiply-accumulate operation on each element of
///                 the provided tuples and returns the accumulator.
/// \tparam         Step (Recursion step.)
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kind (Arithmetic type.)
/// \tparam         Dummy (Dummy type.)
/// \param[in,out]  accumulator Accumulator tuple reference.
/// \param[in]      factor Factor tuple value.
/// \param[in]      coefficient Multiplier value.
/// \return         Reference to the accumulator tuple.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <unsigned int Step, class Tuple, typename Kind, class Dummy, class> 
constexpr Tuple& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::mac(Tuple& accumulator, const Tuple& factor, const Kind coefficient, Dummy&&)
{
    return (Step < std::tuple_size<Tuple>::value) ? (mac<Step+(Step < std::tuple_size<Tuple>::value)>(accumulator, factor, coefficient, mac(std::get<Step*(Step < std::tuple_size<Tuple>::value)>(accumulator), std::get<Step*(Step < std::tuple_size<Tuple>::value)>(factor), coefficient))) : (accumulator);
}

// Multiply-accumulate arithmetic types
/// \brief          Multiply-accumulate arithmetic types.
/// \details        Adds the factor times the coefficient to the accumulator
///                 and returns it.
/// \tparam         Regular (Arithmetic type.)
/// \tparam         Kind (Arithmetic type.)
/// \tparam         Dummy (Dummy types.)
/// \param[in,out]  accumulator Accumulator reference.
/// \param[in]      factor Factor value.
/// \param[in]      coefficient Multiplier value.
/// \return         Reference to the modified accumulator.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Regular, typename Kind, class... Dummy, class> 
constexpr Regular& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::mac(Regular& accumulator, const Regular& factor, const Kind coefficient, Dummy...)
{
    return (accumulator += factor*coefficient);
}

// Multiply-accumulate non arithmetic types
/// \brief          Multiply-accumulate non arithmetic types.
/// \details        Does nothing and returns the accumulator.
/// \tparam         Irregular (Non arithmetic type.)
/// \tparam         Kind (Arithmetic type.)
/// \tparam         Dummy (Dummy types.)
/// \param[in,out]  accumulator Accumulator reference.
/// \param[in]      factor Factor value.
/// \param[in]      coefficient Multiplier value.
/// \return         Reference to the unmodified accumulator.
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
template <class Irregular, typename Kind, class... Dummy, class> 
constexpr Irregular& SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::mac(Irregular& accumulator, const Irregular& factor, const Kind coefficient, const Dummy...)
{
    return (std::integral_constant<bool, static_cast<bool>(sizeof(factor)+sizeof(coefficient))>::value) ? (accumulator) : (accumulator);
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of SimpleHyperOctree. 
/// \return         0 if no error. 
template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
int SimpleHyperOctree<Type, Index, Data, Dimension, Position, Extent, Element, Container>::example()
{
    // Initialize
    std::cout<<"BEGIN = SimpleHyperOctree::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::ostringstream stream;
    SimpleHyperOctreeIndex<unsigned long long int, 3> i(42);
    std::tuple<int, int> t(4, 8);
    std::string s = "forty-two";

    // Construction
    SimpleHyperOctree<double, SimpleHyperOctreeIndex<unsigned long long int, 3>, double> octree(0, 4);
    std::pair<SimpleHyperOctreeIndex<unsigned long long int, 3>, double> element{};

    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Lifecycle : "                                                                                                               <<std::endl;
    std::cout<<std::setw(width*3)<<"SimpleHyperOctree<double, SimpleHyperOctreeIndex<unsigned long long int, 3>, double>() : "                                  ; SimpleHyperOctree<double, SimpleHyperOctreeIndex<unsigned long long int, 3>, double>(); std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"SimpleHyperOctree<double, SimpleHyperOctreeIndex<unsigned long long int, 3>, double>(0, 4) : "                              ; SimpleHyperOctree<double, SimpleHyperOctreeIndex<unsigned long long int, 3>, double>(0, 4); std::cout<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Operators : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"octree == octree : "                                                                <<(octree == octree)<<std::endl;
    std::cout<<std::setw(width*2)<<"octree != octree : "                                                                <<(octree != octree)<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<1>(octree[4]) = 8 : "                                                      <<(std::get<1>(octree[4]) = 8)<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<1>(octree[15]) : "                                                         <<std::get<1>(octree[15])<<std::endl;
    std::cout<<std::setw(width*2)<<"octree().reserve(16) : "                                                            ; octree().reserve(16); std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"octree().size() : "                                                                 <<octree().size()<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<1>(octree(SimpleHyperOctreeIndex<unsigned long long int, 3>())) = 23 : "   <<(std::get<1>(octree(SimpleHyperOctreeIndex<unsigned long long int, 3>())) = 23)<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<1>(octree(SimpleHyperOctreeIndex<unsigned long long int, 3>())) : "        <<std::get<1>(octree(SimpleHyperOctreeIndex<unsigned long long int, 3>()))<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<1>(octree(0.25, 0.25, 0.25)) = 42. : "                                     <<(std::get<1>(octree(0.25, 0.25, 0.25)) = 42.)<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<1>(octree(0.25, 0.25, 0.25)) : "                                           <<std::get<1>(octree(0.25, 0.25, 0.25))<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"octree.assign() : "                           ; octree.assign(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"octree.assign(0, 4) : "                       ; octree.assign(0, 4); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"octree.assign(octree) : "                     ; octree.assign(octree); std::cout<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"octree.nullify() : "                          ; octree.nullify(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"octree.copy() : "                             ; octree.copy(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"octree.cast() : "                             ; octree.cast(); std::cout<<std::endl;

    // Access
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Access : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.at(4)) = 8 : "             <<(std::get<1>(octree.at(4)) = 8)<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.at(15)) : "                <<std::get<1>(octree.at(15))<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.front(16)) = 23 : "        <<(std::get<1>(octree.front(16)) = 23)<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.front(42)) : "             <<std::get<1>(octree.front(42))<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.back(16)) = 23 : "         <<(std::get<1>(octree.back(16)) = 23)<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.back(15)) : "              <<std::get<1>(octree.back(15))<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.cycle(8)) = 4 : "          <<(std::get<1>(octree.cycle(8)) = 4)<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.cycle(4)) : "              <<std::get<1>(octree.cycle(4))<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.container()[8]) = 15 : "   <<(std::get<1>(octree.container()[8]) = 15)<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.container()[23]) : "       <<std::get<1>(octree.container()[23])<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.data()[42]) = 23 : "       <<(std::get<1>(octree.data()[42]) = 23)<<std::endl;
    std::cout<<std::setw(width)<<"std::get<1>(octree.data()[16]) : "            <<std::get<1>(octree.data()[16])<<std::endl;

    // Iterators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Iterators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"&*octree.begin() : "                          <<&*octree.begin()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.cbegin() : "                         <<&*octree.cbegin()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.end() : "                            <<&*octree.end()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.cend() : "                           <<&*octree.cend()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.rbegin() : "                         <<&*octree.rbegin()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.crbegin() : "                        <<&*octree.crbegin()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.rend() : "                           <<&*octree.rend()<<std::endl;
    std::cout<<std::setw(width)<<"&*octree.crend() : "                          <<&*octree.crend()<<std::endl;
    std::cout<<std::setw(width)<<"octree.index(&octree[42]) : "                 <<octree.index(&octree[42])<<std::endl;

    // Search
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Search : "                                    <<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<0>(*octree.find(std::get<0>(*std::begin(octree)+20))) : "                  <<std::get<0>(*octree.find(std::get<0>(*(std::begin(octree)+20))))<<std::endl;
    std::cout<<std::setw(width*2)<<"std::get<0>(*octree.locate(0.4, 0.8, 0.15)) : "                                     <<std::get<0>(*octree.locate(0.4, 0.8, 0.15))<<std::endl;
    
    // Capacity
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Capacity : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"octree.empty() : "                            <<octree.empty()<<std::endl;
    std::cout<<std::setw(width)<<"octree.size() : "                             <<octree.size()<<std::endl;
    std::cout<<std::setw(width)<<"octree.capacity() : "                         <<octree.capacity()<<std::endl;
    std::cout<<std::setw(width)<<"octree.size() : "                             <<octree.size()<<std::endl;
    std::cout<<std::setw(width)<<"octree.shrink().capacity() : "                <<octree.shrink().capacity()<<std::endl;
    std::cout<<std::setw(width)<<"octree.reserve(42).capacity() : "             <<octree.reserve(42).capacity()<<std::endl;

    // Modifiers
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Modifiers : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"octree.copy().clear().size() : "                                                    <<octree.copy().clear().size()<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.copy().pop().size() : "                                                      <<octree.copy().pop().size()<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.copy().append(octree[42]).size() : "                                         <<octree.copy().append(octree[42]).size()<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.copy().resize(42).size() : "                                                 <<octree.copy().resize(42).size()<<std::endl;

    // Refinement
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Refinement : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"octree.update().size() : "                                                          <<octree.update().size()<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.root(octree.begin()) : "                                                     <<octree.root(octree.begin())<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.leaf(octree.begin()) : "                                                     <<octree.leaf(octree.begin())<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.refine(octree.begin()).update().size() : "                                   <<octree.refine(octree.begin()).update().size()<<std::endl;
    std::cout<<std::setw(width*2)<<"octree.coarsen(octree.begin()+1).update().size() : "                                <<octree.coarsen(octree.begin()+1).update().size()<<std::endl;

    // Interpolation
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Interpolation : "                             <<std::endl;
    std::cout<<std::setw(width)<<"octree.ngp(0.25, 0.25, 0.25) : "              <<octree.ngp(0.25, 0.25, 0.25)<<std::endl;
    std::cout<<std::setw(width)<<"octree.cic(0.25, 0.25, 0.25) : "              <<octree.cic(0.25, 0.25, 0.25)<<std::endl;

    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(stream, octree) : "                ; stream<<octree; std::cout<<std::endl;

    // Properties
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Properties : "                                <<std::endl;
    std::cout<<std::setw(width)<<"octree.type() : "                             <<octree.type()<<std::endl;
    std::cout<<std::setw(width)<<"octree.position().num : "                     <<octree.position().num<<std::endl;
    std::cout<<std::setw(width)<<"octree.extent().num : "                       <<octree.extent().num<<std::endl;
    std::cout<<std::setw(width)<<"std::get<0>(octree.element()) : "             <<std::get<0>(octree.element())<<std::endl;
    std::cout<<std::setw(width)<<"octree.dimension() : "                        <<octree.dimension()<<std::endl;

    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Helpers : "                                                                                                                 <<std::endl;
    std::cout<<std::setw(width*3)<<"std::get<0>(octree.tupleify<std::tuple<int, int> >(std::make_tuple(4, 8))) : "                                              <<std::get<0>(octree.tupleify<std::tuple<int, int> >(std::make_tuple(4, 8)))<<std::endl;
    std::cout<<std::setw(width*3)<<"std::get<0>(octree.tupleify<std::tuple<int, int> >(std::array<int, 2>()) : "                                                <<std::get<0>(octree.tupleify<std::tuple<int, int> >(std::array<int, 2>()))<<std::endl;
    std::cout<<std::setw(width*3)<<"std::get<0>(octree.tupleify<std::tuple<int, int> >(15, 16) : "                                                              <<std::get<0>(octree.tupleify<std::tuple<int, int> >(15, 16))<<std::endl;
    std::cout<<std::setw(width*3)<<"octree.variadify<int>([](int i, int j){return i+j;}, std::make_tuple(23, 42)) : "                                           <<octree.variadify<int>([](int i, int j){return i+j;}, std::make_tuple(23, 42))<<std::endl;
    std::cout<<std::setw(width*3)<<"octree.variadify<int>([](int i, int j){return i+j;}, std::array<int, 2>()) : "                                              <<octree.variadify<int>([](int i, int j){return i+j;}, std::array<int, 2>())<<std::endl;
    std::cout<<std::setw(width*3)<<"octree.variadify<int>([](int i, int j){return i+j;}, 4, 8) : "                                                              <<octree.variadify<int>([](int i, int j){return i+j;}, 4, 8)<<std::endl;
    std::cout<<std::setw(width*3)<<"octree.mac(i, i, 15) : "                                                                                                    <<octree.mac(i, i, 15)<<std::endl;
    std::cout<<std::setw(width*3)<<"std::get<1>(octree.mac(t, std::make_tuple(16, 23), 4)) : "                                                                  <<std::get<1>(octree.mac(t, std::make_tuple(16, 23), 4))<<std::endl;
    std::cout<<std::setw(width*3)<<"octree.mac(std::get<0>(t), 42, 4) : "                                                                                       <<octree.mac(std::get<0>(t), 42, 4)<<std::endl;
    std::cout<<std::setw(width*3)<<"octree.mac(s, std::string(\"twenty three\"), 4) : "                                                                         <<octree.mac(s, std::string("twenty three"), 4)<<std::endl;
        
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = SimpleHyperOctree::example()"<<std::endl;
    return 0;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // SIMPLEHYPEROCTREE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
