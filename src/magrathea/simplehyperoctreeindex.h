/* ************************* SIMPLEHYPEROCTREEINDEX ************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          SimpleHyperOctreeIndex 
// DESCRIPTION :    A simple hyperoctree index based on an integer
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           simplehyperoctreeindex.h 
/// \brief          A simple hyperoctree index based on an integer
/// \author         Vincent Reverdy (vince.rev@gmail.com) 
/// \date           2012-2013 
/// \copyright      CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef SIMPLEHYPEROCTREEINDEX_H_INCLUDED
#define SIMPLEHYPEROCTREEINDEX_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



//------------------------------- PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <typeinfo>
#include <limits>
#include <cstdint>
#include <string>
#include <utility>
#include <algorithm>
#include <bitset>
#include <array>
#include <tuple>
#include <ratio>
// Include libs
// Include project
// Misc
namespace magrathea {
//--------------------------------------------------------------------------- //



//---------------------------------- CLASS ---------------------------------- //
// A simple hyperoctree index based on an integer
/// \brief          A simple hyperoctree index based on an integer.
/// \details        Implements a simple hyperoctree index with no dependency 
///                 using a single integer and providing all standard 
///                 operations to easily recover level and position.
/// \tparam         Type Unsigned integer type.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Bits Size of the type in bits.
template <typename Type = unsigned long long int, unsigned int Dimension = 3, unsigned int Bits = sizeof(Type)*std::numeric_limits<unsigned char>::digits> 
class SimpleHyperOctreeIndex final
{
    // Setup
    static_assert(((std::is_integral<Type>::value) ? (std::is_unsigned<Type>::value) : (std::is_convertible<Type, int>::value)) && (!std::is_floating_point<Type>::value), "ERROR = SimpleHyperOctreeIndex : non compliant type");
    static_assert((Dimension > 0) && (Dimension < Bits), "ERROR = SimpleHyperOctreeIndex : non compliant dimension");
    static_assert(Bits == sizeof(Type)*std::numeric_limits<unsigned char>::digits, "ERROR = SimpleHyperOctreeIndex : non compliant number of bits");
    
    // Lifecycle 
    /// \name           Lifecycle 
    //@{
    public: 
        inline SimpleHyperOctreeIndex();
        explicit inline SimpleHyperOctreeIndex(const Type source);
        template <class String, class = typename std::enable_if<(std::is_convertible<String, std::string>::value) && (std::is_convertible<decltype(std::declval<String>()[0]), char>::value)>::type> explicit inline SimpleHyperOctreeIndex(const String& source);
    //@}
    
    // Operators
    /// \name           Operators
    //@{
    public:
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& operator=(const Type rhs);
        template <class String, class = typename std::enable_if<(std::is_convertible<String, std::string>::value) && (std::is_convertible<decltype(std::declval<String>()[0]), char>::value)>::type> inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& operator=(const String& rhs);
        inline operator Type() const;
        inline Type& operator()();
        inline const Type& operator()() const;
        inline bool operator[](const unsigned int ibit) const;
    //@}

    // Assignment 
    /// \name           Assignment 
    //@{
    public: 
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& assign(); 
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& assign(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& source); 
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& assign(const Type source); 
        template <class String, class = typename std::enable_if<(std::is_convertible<String, std::string>::value) && (std::is_convertible<decltype(std::declval<String>()[0]), char>::value)>::type> inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& assign(const String& source); 
    //@}
    
    // Management
    /// \name           Management
    //@{
    public: 
        inline Type& data();
        inline const Type& data() const;
        inline Type get() const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& set(const Type value);
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& nullify();
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> copy() const;
        template <class OtherType = SimpleHyperOctreeIndex<Type, Dimension, Bits>, class = typename std::enable_if<std::is_constructible<OtherType, Type>::value>::type> inline OtherType cast() const;
        template <typename Base = std::true_type, class = typename std::enable_if<((std::is_same<Base, std::true_type>::value) || (std::is_convertible<Base, int>::value)) && (!std::is_floating_point<Base>::value)>::type> std::string stringify(const Base base = Base()) const;
    //@}
    
    // Core
    /// \name           Core
    //@{
    public:
        inline unsigned int level() const;
        inline bool coarsest() const;
        inline bool finest() const;
        inline bool limited() const;
        inline bool check() const;
        inline bool invalidated() const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& fix();
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& invalidate();
    //@}
    
    // Tree
    /// \name           Tree
    //@{
    public:
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> parent() const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> child(const unsigned int isite) const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> brother(const unsigned int isite) const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> preceding() const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> following() const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> previous(const unsigned int ilvl = 0, const unsigned int nref = Bits/(Dimension+1)) const;
        inline SimpleHyperOctreeIndex<Type, Dimension, Bits> next(const unsigned int ilvl = 0, const unsigned int nref = Bits/(Dimension+1)) const;
    //@}
    
    // Curve
    /// \name           Curve
    //@{
    public:
        template <unsigned int From = 0, unsigned int To = Bits/(Dimension+1), unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class Level = std::true_type, class = typename std::enable_if<(From <= Mid) && (Mid <= To) && (To <= Bits/(Dimension+1)) && (((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value))>::type> inline Type linear(const Level ilvl = Level()) const;
        inline Type zcurve() const;
        template <class Level = std::true_type, class = typename std::enable_if<((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value)>::type> inline Type site(const Level ilvl = Level()) const;
        template <unsigned int From = 0, unsigned int To = Dimension-1, unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class = typename std::enable_if<(From <= Mid) && (Mid <= To) && (To <= Dimension)>::type> inline Type coordinate(const unsigned int idim) const;
    //@}

    // Position
    /// \name           Position
    //@{
    public:
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, unsigned int From = 0, unsigned int To = Dimension-1, unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class Reference = typename std::ratio_subtract<Position, typename std::ratio_divide<Extent, std::ratio<2> > >, class Level = std::true_type, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0) && (From <= Mid) && (Mid <= To) && (To <= Dimension) && (((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value))>::type> inline Kind position(const unsigned int idim, const Level ilvl = Level()) const;
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, unsigned int From = 0, unsigned int To = Dimension-1, unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class Reference = typename std::ratio_subtract<Position, typename std::ratio_divide<Extent, std::ratio<2> > >, class Level = std::true_type, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0) && (From <= Mid) && (Mid <= To) && (To <= Dimension) && (((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value))>::type> inline Kind center(const unsigned int idim, const Level ilvl = Level()) const;
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, unsigned int From = 0, unsigned int To = Dimension-1, unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class Reference = typename std::ratio_subtract<Position, typename std::ratio_divide<Extent, std::ratio<2> > >, class Level = std::true_type, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0) && (From <= Mid) && (Mid <= To) && (To <= Dimension) && (((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value))>::type> inline Kind minimum(const unsigned int idim, const Level ilvl = Level()) const;
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, unsigned int From = 0, unsigned int To = Dimension-1, unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class Reference = typename std::ratio_subtract<Position, typename std::ratio_divide<Extent, std::ratio<2> > >, class Level = std::true_type, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0) && (From <= Mid) && (Mid <= To) && (To <= Dimension) && (((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value))>::type> inline Kind maximum(const unsigned int idim, const Level ilvl = Level()) const;
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0)>::type> inline Kind extent() const;
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0)>::type> inline Kind length() const;
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0)>::type> inline Kind volume() const;
    //@}
    
    // Space
    /// \name           Space
    //@{
    public:
        inline bool containing(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const;
        inline bool contained(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const;
        inline bool intersecting(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const;
        inline bool adjoining(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const;
    //@}

    // Stream 
    /// \name           Stream 
    //@{
    public: 
        template <typename SelfType, unsigned int SelfDimension, unsigned int SelfBits> friend std::ostream& operator<<(std::ostream& lhs, const SimpleHyperOctreeIndex<SelfType, SelfDimension, SelfBits>& rhs); 
    //@}
    
    // Computation
    /// \name           Computation
    //@{
    public:
        template <unsigned int From = 0, unsigned int To = Bits/(Dimension+1), unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class = typename std::enable_if<(From <= Mid) && (Mid <= To) && (To <= Bits/(Dimension+1))>::type> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> delinearize(const unsigned int ilvl, const Type ilinear);
        template <unsigned int From = 0, unsigned int To = Bits/(Dimension+1), unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class = typename std::enable_if<(From <= Mid) && (Mid <= To) && (To <= Bits/(Dimension+1))>::type> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> dezcurvify(const unsigned int ilvl, const Type izcurve);
        template <typename... Types, class = typename std::enable_if<((std::is_convertible<typename std::tuple_element<0, typename std::conditional<static_cast<bool>(sizeof...(Types)), std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...>, std::tuple<Type> >::type>::type, Type>::value) ? (true) : (sizeof...(Types) == 1))>::type> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> navigate(Types&&... isites);
        template <unsigned int From = 0, unsigned int To = Bits/(Dimension+1), unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), typename... Types, class = typename std::enable_if<(From <= Mid) && (Mid <= To) && (To <= Bits/(Dimension+1)) && ((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Types>::type>::type...> >::type, Type>::value) ? (sizeof...(Types) == Dimension) : (sizeof...(Types) == 1))>::type> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> make(const unsigned int ilvl, Types&&... icoords);
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, unsigned int From = 0, unsigned int To = Bits/(Dimension+1), unsigned int Mid = From+(To-From)/2, bool Condition = (From != To), class Reference = typename std::ratio_subtract<Position, typename std::ratio_divide<Extent, std::ratio<2> > >, typename... Kinds, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0) && (From <= Mid) && (Mid <= To) && (To <= Bits/(Dimension+1)) && ((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Kinds>::type>::type...> >::type, Type>::value) ? (sizeof...(Kinds) == Dimension) : (sizeof...(Kinds) == 1))>::type> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> compute(const unsigned int ilvl, Kinds&&... iposs);
        template <typename Kind = double, class Position = std::ratio<0>, class Extent = std::ratio<1>, class Reference = typename std::ratio_subtract<Position, typename std::ratio_divide<Extent, std::ratio<2> > >, typename... Kinds, class = typename std::enable_if<(std::is_floating_point<Kind>::value) && (Position::den != 0) && (Extent::num >= 0) && ((std::is_convertible<typename std::tuple_element<0, std::tuple<typename std::remove_cv<typename std::remove_reference<Kinds>::type>::type...> >::type, Type>::value) ? (sizeof...(Kinds) == Dimension) : (sizeof...(Kinds) == 1))>::type> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> cipher(Kinds&&... iposs);
        template <class... Indices> static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> merge(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx, Indices&&... idcs);
    //@}

    // Properties 
    /// \name           Properties 
    //@{
    public: 
        static constexpr Type type(); 
        static constexpr unsigned int types(); 
        static constexpr unsigned int dimension(); 
        static constexpr unsigned int bits(); 
        static constexpr unsigned int size();
        static constexpr unsigned int instructions();
        static constexpr unsigned int opcodes();
        static constexpr unsigned int operands();
        static constexpr unsigned int remainder();
        static constexpr unsigned int refinements();
        static constexpr unsigned int sites();
        template <class Level = std::true_type, class = typename std::enable_if<((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value)>::type> static constexpr Type indices(const Level ilvl = Level(), const unsigned int nref = 0);
        template <class Level = std::true_type, class = typename std::enable_if<((std::is_same<Level, std::true_type>::value) || (std::is_convertible<Level, int>::value)) && (!std::is_floating_point<Level>::value)>::type> static constexpr Type subdivisions(const Level ilvl = Level());
        static constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> invalid();
    //@}
    
    // Helpers 
    /// \name           Helpers 
    //@{
    public:
        template <class Template, class Function, class Tuple, typename... Kinds, class = typename std::enable_if<(sizeof...(Kinds)+1 <= std::tuple_size<Template>::value)>::type> static constexpr Template apply(Function&& f, Tuple&& tuple, Kinds&&... values);
        template <class Template, class Function, class Tuple, typename... Kinds, class = typename std::enable_if<(sizeof...(Kinds) == std::tuple_size<Template>::value)>::type> static constexpr Template apply(Function&& f, const Tuple& tuple, Kinds&&... values);
        template <class Template, class Function, typename Kind, class = typename std::enable_if<std::is_convertible<Template, int>::value>::type> static constexpr Template apply(Function&& f, Kind&& value);
        template <typename Kind, class Operation = std::plus<Kind>, unsigned int Index = 0, class Tuple, bool Condition = (Index+2 <= std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value), class = typename std::enable_if<(std::is_convertible<Kind, int>::value) && (std::is_convertible<typename std::result_of<Operation(Kind, Kind)>::type, Kind>::value) && (std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value >= 1)>::type> static constexpr Kind accumulate(Tuple&& tuple);
        template <typename Kind, class Operation = std::plus<Kind>, typename... Kinds, class = typename std::enable_if<(std::is_convertible<Kind, int>::value) && (std::is_convertible<typename std::result_of<Operation(Kind, Kind)>::type, Kind>::value)>::type> static constexpr Kind accumulate(const Kind value, const Kinds... values);
        template <typename Kind, class Operation = std::plus<Kind>, class = typename std::enable_if<(std::is_convertible<Kind, int>::value) && (std::is_convertible<typename std::result_of<Operation(Kind, Kind)>::type, Kind>::value)>::type> static constexpr Kind accumulate(const Kind value);
        template <typename Kind, int Exponent = 1, int One = 1, bool Greater = (Exponent > 1), bool Less = (Exponent < 0), bool Equal = (Exponent == 1), class = typename std::enable_if<std::is_convertible<Kind, int>::value>::type> static constexpr Kind pow(const Kind value);
        template <typename Integer, Integer Zero = Integer(), Integer One = Integer(1), Integer Ones = ~Zero, Integer Size = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer block(const Integer location = Zero, const Integer length = Size);
        template <typename Integer, Integer Index = Integer(), Integer Zero = Integer(), Integer One = Integer(1), Integer Condition = (Index+One <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer periodic(const Integer period = One, const Integer offset = Zero);
        template <typename Integer, Integer Index = Integer(), Integer Zero = Integer(), Integer One = Integer(1), Integer Size = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, Integer Condition = (Index+One <= Size), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer comb(const Integer period = One, const Integer offset = Zero);
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = Integer(1), Integer One = Integer(1), Integer Ones = ~Integer(), Integer Condition = (Step+One <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer nhp(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = Integer(1), Integer One = Integer(1), Integer Ones = ~Integer(), Integer Condition = (Step+One <= nhp<Integer>(sizeof(Integer)*std::numeric_limits<unsigned char>::digits)), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer bhsmask(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = nhp<Integer>(sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Zero = Integer(), Integer One = Integer(1), Integer Ones = ~Zero, Integer Size = nhp<Integer>(sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Temporary = block<Integer>(Step, Step), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer lzcnt(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = nhp<Integer>(sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Zero = Integer(), Integer One = Integer(1), Integer Ones = ~Zero, Integer Size = nhp<Integer>(sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Temporary = periodic<Integer>(Step*(One+One)+(Step == Zero)), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer tzcnt(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = Integer(), Integer One = Integer(1), Integer Condition = (Step+One <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Temporary = ((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer popcnt(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = Integer(), Integer Shift = Integer(), Integer One = Integer(1), Integer Condition = (Step+One <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Temporary = ((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer pext(const Integer value); 
        template <typename Integer, Integer Mask = ~Integer(), Integer Step = Integer(), Integer Shift = Integer(), Integer One = Integer(1), Integer Condition = (Step+One <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits), Integer Temporary = ((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer pdep(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Period = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, Integer Step = Integer(), Integer Count = Integer(), Integer Zero = Integer(), Integer One = Integer(1), Integer Size = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, Integer Condition = (Step+One <= Size), Integer Population = popcnt<Integer>(Mask)+((Period-(popcnt<Integer>(Mask)%Period))*(popcnt<Integer>(Mask)%Period != Zero)), Integer Destination = (((((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)) ? (Count) : (Zero))%((Population+(Period*(Population+One <= Period)))/Period))*Period)+((((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)) ? (Count) : (Zero))/((Population+(Period*(Population+One <= Period)))/Period)), Integer Temporary = ((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)) && (Destination < Size), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer itlc(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Period = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, Integer Step = Integer(), Integer Count = Integer(), Integer Zero = Integer(), Integer One = Integer(1), Integer Size = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, Integer Condition = (Step+One <= Size), Integer Population = popcnt<Integer>(Mask)+((Period-(popcnt<Integer>(Mask)%Period))*(popcnt<Integer>(Mask)%Period != Zero)), Integer Destination = (((((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)) ? (Count) : (Zero))%Period)*(Population/Period))+((((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)) ? (Count) : (Zero))/Period), Integer Temporary = ((Condition) ? ((Mask >> Step) & (Condition)) : (Condition)) && (Destination < Size), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer dtlc(const Integer value);
        template <typename Integer, Integer Mask = ~Integer(), Integer Length = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, bool Msb = false, std::size_t Step = Integer(), Integer Zero = Integer(), Integer Direction = ((!Msb) || ((Length*(Step+1)) <= (sizeof(Integer)*std::numeric_limits<unsigned char>::digits))), Integer Left = ((Msb) ? (sizeof(Integer)*std::numeric_limits<unsigned char>::digits-(Length*(Step+1))) : (Length*Step))*(Direction), Integer Right = ((Msb) ? ((Length*(Step+1))-sizeof(Integer)*std::numeric_limits<unsigned char>::digits) : (Length*Step))*(!Direction), Integer Condition = ((Left+1 <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits) && (Right+1 <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits) && (Right+1 <= Length)), class Tuple, Integer Count = ((std::tuple_size<typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::value)-(Step+1)), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value) && (std::is_convertible<typename std::tuple_element<Step, typename std::remove_cv<typename std::remove_reference<Tuple>::type>::type>::type, Integer>::value)>::type> static constexpr Integer glue(Tuple&& tuple);
        template <typename Integer, Integer Mask = ~Integer(), Integer Length = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, bool Msb = false, std::size_t Step = Integer(), Integer Zero = Integer(), Integer Direction = ((!Msb) || ((Length*(Step+1)) <= (sizeof(Integer)*std::numeric_limits<unsigned char>::digits))), Integer Left = ((Msb) ? (sizeof(Integer)*std::numeric_limits<unsigned char>::digits-(Length*(Step+1))) : (Length*Step))*(Direction), Integer Right = ((Msb) ? ((Length*(Step+1))-sizeof(Integer)*std::numeric_limits<unsigned char>::digits) : (Length*Step))*(!Direction), Integer Condition = ((Left+1 <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits) && (Right+1 <= sizeof(Integer)*std::numeric_limits<unsigned char>::digits) && (Right+1 <= Length)), typename... Integers, Integer Count = sizeof...(Integers), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer glue(const Integer value, const Integers... values);
        template <typename Integer, Integer Mask = ~Integer(), Integer Length = sizeof(Integer)*std::numeric_limits<unsigned char>::digits, bool Msb = false, std::size_t Step = Integer(), Integer Zero = Integer(), class = typename std::enable_if<((std::is_integral<Integer>::value) ? (std::is_unsigned<Integer>::value) : (std::is_convertible<Integer, int>::value)) && (!std::is_floating_point<Integer>::value)>::type> static constexpr Integer glue();
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
        Type _data;                                                             ///< Internal index value. 
    //@}
};
//--------------------------------------------------------------------------- //



//-------------------------------- LIFECYCLE -------------------------------- //
// Implicit empty constructor 
/// \brief          Implicit empty constructor. 
/// \details        Provides an implicit construction of the index initialized 
///                 to its default value. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>::SimpleHyperOctreeIndex()
: _data(Type())
{
    ;
}

// Explicit value constructor 
/// \brief          Explicit value constructor. 
/// \details        Provides an explicit construction of the index initialized 
///                 to a particular value. 
/// \param[in]      source Source of the copy.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>::SimpleHyperOctreeIndex(const Type source)
: _data(source)
{
    ;
} 

// Explicit string constructor
/// \brief          Explicit string constructor.
/// \details        Provide an explicit construction from a string of zeros and
///                 ones. The index is filled from the most significant bit.
/// \tparam         String String type.
/// \param[in]      source Source of the copy.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class String, class>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>::SimpleHyperOctreeIndex(const String& source)
: _data()
{
    static const char czero = '0';
    static const char cone = '1';
    static const Type zero = Type();
    static const Type one = Type(1);
    const unsigned int length = std::end(source)-std::begin(source);
    Type ibit = zero;
    for (unsigned int ichar = 0; (ichar < length) && (ibit < Bits); ++ichar) {
        (source[ichar] == cone) ? (_data |= (one << (Bits-(++ibit)))) : (_data |= zero*(ibit += (source[ichar] == czero)));
    }
} 
//--------------------------------------------------------------------------- //



//-------------------------------- OPERATORS -------------------------------- //
// Value assignment operator 
/// \brief          Value assignment operator. 
/// \details        Assigns data from a value. 
/// \param[in]      rhs Right-hand side.
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::operator=(const Type rhs)
{
    _data = rhs;
    return *this;
}

// String assignment operator 
/// \brief          String assignment operator. 
/// \details        Assigns data from a string of zeros and ones. The index is 
///                 filled from the most significant bit.
/// \tparam         String String type.
/// \param[in]      rhs Right-hand side.
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class String, class>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::operator=(const String& rhs)
{
    static const char czero = '0';
    static const char cone = '1';
    static const Type zero = Type();
    static const Type one = Type(1);
    const unsigned int length = std::end(rhs)-std::begin(rhs);
    Type ibit = zero;
    _data = zero;
    for (unsigned int ichar = 0; (ichar < length) && (ibit < Bits); ++ichar) {
        (rhs[ichar] == cone) ? (_data |= (one << (Bits-(++ibit)))) : (_data |= zero*(ibit += (rhs[ichar] == czero)));
    }
    return *this;
}

// Implicit cast operator
/// \brief          Implicit cast operator.
/// \details        Implicitely converts the index to an integer.
/// \return         Copy of the underlying data.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>::operator Type() const
{
    return _data;
}

// Data access operator
/// \brief          Data access operator. 
/// \details        Provides direct access to internal data. 
/// \return         Reference to underlying data. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline Type& SimpleHyperOctreeIndex<Type, Dimension, Bits>::operator()()
{
    return _data;
}

// Immutable data access operator
/// \brief          Immutable data access operator. 
/// \details        Provides an immutable direct access to internal data. 
/// \return         Immutable reference to underlying data. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline const Type& SimpleHyperOctreeIndex<Type, Dimension, Bits>::operator()() const
{
    return _data;
}

// Immutable bit access
/// \brief          Immutable bit access. 
/// \details        Provides an immutable access to the i-th bit. 
/// \param[in]      ibit Bit index.
/// \return         Copy of the i-th bit. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::operator[](const unsigned int ibit) const
{
    return (_data & (std::integral_constant<Type, 1>::value << Type(ibit))) && (ibit < Bits);
}
//--------------------------------------------------------------------------- //



//-------------------------------- ASSIGNMENT ------------------------------- //
// Empty assignment 
/// \brief          Empty assignment. 
/// \details        Assigns contents from an index initialized to its default 
///                 value. 
/// \return         Self reference.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::assign()
{
    _data = Type();
    return *this;
}

// Copy assignment 
/// \brief          Copy assignment. 
/// \details        Assigns contents from the same type of index. 
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::assign(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& source)
{
    _data = source._data;
    return *this;
}

// Value assignment 
/// \brief          Value assignment. 
/// \details        Assigns contents from an index value. 
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::assign(const Type source)
{
    _data = source;
    return *this;
}

// String assignment 
/// \brief          String assignment. 
/// \details        Assigns contents from a string of zeros and ones. The index 
///                 is filled from the most significant bit.
/// \tparam         String String type.
/// \param[in]      source Source of the copy.
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class String, class> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::assign(const String& source) 
{
    static const char czero = '0';
    static const char cone = '1';
    static const Type zero = Type();
    static const Type one = Type(1);
    const unsigned int length = std::end(source)-std::begin(source);
    Type ibit = zero;
    _data = zero;
    for (unsigned int ichar = 0; (ichar < length) && (ibit < Bits); ++ichar) {
        (source[ichar] == cone) ? (_data |= (one << (Bits-(++ibit)))) : (_data |= zero*(ibit += (source[ichar] == czero)));
    }
    return *this;
}
//--------------------------------------------------------------------------- //



//-------------------------------- MANAGEMENT ------------------------------- //
// Data access
/// \brief          Data access. 
/// \details        Provides direct access to internal data. 
/// \return         Reference to underlying data. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline Type& SimpleHyperOctreeIndex<Type, Dimension, Bits>::data()
{
    return _data;
}

// Immutable data access
/// \brief          Immutable data access. 
/// \details        Provides an immutable direct access to internal data. 
/// \return         Immutable reference to underlying data. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline const Type& SimpleHyperOctreeIndex<Type, Dimension, Bits>::data() const
{
    return _data;
}

// Data getter
/// \brief          Data getter. 
/// \details        Returns a copy of the internal data. 
/// \return         Copy of the underlying data. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::get() const
{
    return _data;
}

// Data setter
/// \brief          Data setter. 
/// \details        Sets the internal data to the provided value.
/// \param[in]      value Input value.
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::set(const Type value)
{
    _data = value;
    return *this;
}

// Nullify 
/// \brief          Nullify. 
/// \details        Resets the internal data to its default value. 
/// \return         Self reference. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::nullify()
{
    _data = std::integral_constant<Type, 0>::value;
    return *this;
}

// Copy 
/// \brief          Copy. 
/// \details        Generates a copy of the index. 
/// \return         Copy. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::copy() const
{
    return *this;
}

// Cast 
/// \brief          Cast. 
/// \details        Casts the index to another object type. 
/// \tparam         OtherType Other data type.
/// \return         Casted copy. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class OtherType, class> 
inline OtherType SimpleHyperOctreeIndex<Type, Dimension, Bits>::cast() const
{
    return OtherType(_data);
}

// Stringify
/// \brief          Stringify.
/// \details        Converts the index to a string. If no base is specified, the
///                 standard display is used.
/// \tparam         Base Integral base type.
/// \param[in]      base Integral base value.
/// \return         String corresponding to the index.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Base, class> 
std::string SimpleHyperOctreeIndex<Type, Dimension, Bits>::stringify(const Base base) const
{
    static const unsigned int n = std::integral_constant<unsigned int, Bits-1>::value;
    static const unsigned int d = std::integral_constant<unsigned int, Dimension+1>::value;
    static const char czero = '0';
    static const char cone = '1';
    static const char calpha = 'A';
    static const char separator = ' ';
    static const char ten = static_cast<char>(10);
    static const Type one = Type(1);
    std::string result;
    Type number = _data;
    char remainder = 0;
    if (std::is_same<Base, std::true_type>::value) {
        for (unsigned int i = 0; i <= n; ++i) {
            (((i > 0) && (i%d == 0)) ? (result += separator) : (result));
            result += ((_data & (one<<(n-i))) ? (cone) : (czero));
            (((i < n) && (i%d == 0)) ? (result += separator) : (result));
        }
    } else if (base > 0) {
        do {
            remainder = number%base;
            result += (remainder < ten) ? (czero+remainder) : (calpha+remainder-ten);
            number /= base;
        } while (number > 0);
        std::reverse(std::begin(result), std::end(result));
    }
    return result;
}
//--------------------------------------------------------------------------- //



//----------------------------------- CORE ---------------------------------- //
// Level
/// \brief          Level.
/// \details        Computes the level of refinement according to the number of
///                 trailing zeros.
/// \return         Refinement level.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::level() const
{
    return ((Bits-tzcnt<Type, (comb<Type>(Dimension+1, Dimension) << (Bits%(Dimension+1)))>(_data))/std::integral_constant<Type, Dimension+1>::value)+(_data != std::integral_constant<Type, 0>::value);
}

// Coarsest level
/// \brief          Coarsest level.
/// \details        Checks whether the index corresponds to the coarsest, 
///                 unrefined level.
/// \return         True if coarsest level, false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::coarsest() const
{
    return (_data == std::integral_constant<Type, 0>::value);
}

// Finest level
/// \brief          Finest level.
/// \details        Checks whether the index corresponds to the finest, 
///                 most refined level.
/// \return         True if finest level, false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::finest() const
{
    return (_data & std::integral_constant<Type, (std::integral_constant<Type, 1>::value << (Dimension+(Bits%(Dimension+1))))>::value);
}

// Limited 
/// \brief          Limited.
/// \details        Checks whether the index corresponds to the coarsest
///                 or finest level of refinement.
/// \return         True if coarsest or finest level, false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::limited() const
{
    return ((_data & std::integral_constant<Type, (std::integral_constant<Type, 1>::value << (Dimension+(Bits%(Dimension+1))))>::value) || (_data == std::integral_constant<Type, 0>::value)); 
}

// Check
/// \brief          Check.
/// \details        Checks that the underlying integer represents an index with
///                 no error.
/// \return         True if no error, false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::check() const
{
    return !((~block<Type>(Bits-(popcnt<Type, (comb<Type>(Dimension+1, Dimension) << (Bits%(Dimension+1)))>(_data)*std::integral_constant<Type, Dimension+1>::value))) & _data);
}

// Invalidated
/// \brief          Invalidated.
/// \details        Checks whether the index is invalidated.
/// \return         True if the index is invalidated, false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::invalidated() const
{
    return (_data == std::integral_constant<Type, block<Type>(0, Bits-1)>::value);
}

// Fix
/// \brief          Fix.
/// \details        Fixes the underlying integer if it does not represents a 
///                 correct index. To do so, every bit set after the last 
///                 correct level is cleared.
/// \return         Self reference.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::fix()
{
    static const Type zero = Type();
    static const Type one = Type(1);
    static const Type mask = std::integral_constant<Type, comb<Type>(Dimension+one, Dimension)>::value;
    static const Type first = std::integral_constant<Type, block<Type>(Bits-one)>::value;
    static const Type dimensionp = std::integral_constant<Type, Dimension+one>::value;
    const Type marker = Bits-(popcnt<Type, (mask << (Bits%(dimensionp)))>(_data)*dimensionp);
    const bool error = (~block<Type>(marker)) & _data;
    Type temporary = zero;
    if (error) {
        _data &= block<Type>(marker);
        temporary = (mask << marker)*(marker < Bits);
        temporary = (_data & temporary) ^ temporary;
        if (temporary) {
            _data &= (~(nhp<Type>(temporary)-one))*((_data & first) != zero);
        }
    }
    return *this;
}

// Invalidate
/// \brief          Invalidate.
/// \details        Invalidates the index.
/// \return         Self reference.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits>& SimpleHyperOctreeIndex<Type, Dimension, Bits>::invalidate()
{
    _data = std::integral_constant<Type, block<Type>(0, Bits-1)>::value;
    return *this;
}
//--------------------------------------------------------------------------- //



//----------------------------------- TREE ---------------------------------- //
// Parent index
/// \brief          Parent index.
/// \details        Computes the parent index in the tree.
/// \return         Index of the parent.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::parent() const
{
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((_data) & (block<Type>(std::integral_constant<Type, Bits+Dimension+1>::value-(level()*std::integral_constant<Type, Dimension+1>::value))));
}

// Child index
/// \brief          Child index.
/// \details        Computes the child index in the tree.
/// \param[in]      isite Site of the child.
/// \return         Index of the specified child.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::child(const unsigned int isite) const
{
    const Type ilvl = level();
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((_data | ((((Type(isite) & (std::integral_constant<Type, block<Type>(0, Dimension)>::value)) | std::integral_constant<Type, (std::integral_constant<Type, 1>::value << Dimension)>::value) << (Bits-((ilvl+std::integral_constant<Type, 1>::value)*std::integral_constant<Type, Dimension+1>::value)))*(ilvl < std::integral_constant<Type, Bits/(Dimension+1)>::value))));
}

// Brother index
/// \brief          Brother index.
/// \details        Computes the brother index in the tree, which is an index
///                 with the same final parent.
/// \param[in]      isite Site of the brother.
/// \return         Index of the specified brother.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::brother(const unsigned int isite) const
{
    const Type ilvl = Type(level())-std::integral_constant<Type, 1>::value;
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((ilvl <= std::integral_constant<Type, Bits/(Dimension+1)>::value) ? (((_data) & (block<Type>(std::integral_constant<Type, Bits>::value-(ilvl*std::integral_constant<Type, Dimension+1>::value)))) | ((((Type(isite) & (std::integral_constant<Type, block<Type>(0, Dimension)>::value)) | std::integral_constant<Type, (std::integral_constant<Type, 1>::value << Dimension)>::value) << (Bits-((ilvl+std::integral_constant<Type, 1>::value)*std::integral_constant<Type, Dimension+1>::value)))*(ilvl < std::integral_constant<Type, Bits/(Dimension+1)>::value))) : (_data));
}

// Preceding brother index
/// \brief          Preceding brother index.
/// \details        Computes the index of the preceding brother in the tree, 
///                 which is an index with the same final parent. A cyclic 
///                 operation is performed if the first brother of this final
///                 parent has already been reached.
/// \return         Index of the preceding brother.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::preceding() const
{
    const Type shift = std::integral_constant<Type, Bits>::value-(level()*std::integral_constant<Type, Dimension+1>::value);
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((_data) ^ (((_data) ^ ((_data)-((std::integral_constant<Type, 1>::value << shift)*(shift < std::integral_constant<Type, Bits>::value)))) & (block<Type>(shift, std::integral_constant<Type, Dimension>::value))));
}

// Following brother index
/// \brief          Following brother index.
/// \details        Computes the index of the following brother in the tree, 
///                 which is an index with the same final parent. A cyclic 
///                 operation is performed if the last brother of this final
///                 parent has already been reached.
/// \return         Index of the following brother.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::following() const
{
    const Type shift = std::integral_constant<Type, Bits>::value-(level()*std::integral_constant<Type, Dimension+1>::value);
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((_data) ^ (((_data) ^ ((_data)+((std::integral_constant<Type, 1>::value << shift)*(shift < std::integral_constant<Type, Bits>::value)))) & (block<Type>(shift, std::integral_constant<Type, Dimension>::value))));
}

// Previous index
/// \brief          Previous index.
/// \details        Computes the previous index in the tree from the specified 
///                 level to the specified number of refinements. A cyclic 
///                 operation is performed if the first index has been reached.
/// \param[in]      ilvl Index of the first level.
/// \param[in]      nref Number of refinements of the first level.
/// \return         Index of the previous element.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::previous(const unsigned int ilvl, const unsigned int nref) const
{
    const Type flvl = (ilvl <= std::integral_constant<Type, Bits/(Dimension+1)>::value) ? (Type(ilvl)) : (std::integral_constant<Type, Bits/(Dimension+1)>::value);
    const Type mlvl = (flvl+nref <= std::integral_constant<Type, Bits/(Dimension+1)>::value) ? (Type(flvl+nref)) : (std::integral_constant<Type, Bits/(Dimension+1)>::value);
    const Type clvl = level();
    const Type tlvl = (clvl <= mlvl) ? (clvl) : (mlvl);
    const Type shift = ((_data & block<Type>(std::integral_constant<Type, Bits>::value-(clvl*std::integral_constant<Type, Dimension+1>::value), std::integral_constant<Type, Dimension>::value))+(!_data)) ? (std::integral_constant<Type, Bits>::value-(mlvl*std::integral_constant<Type, Dimension+1>::value)) : (std::integral_constant<Type, Bits+Dimension>::value-(clvl*std::integral_constant<Type, Dimension+1>::value));
    const SimpleHyperOctreeIndex<Type, Dimension, Bits> result((_data & block<Type>(shift))-((std::integral_constant<Type, 1>::value << shift)*(clvl == tlvl)*(mlvl > std::integral_constant<Type, 0>::value)));
    return (result.level() >= flvl) ? (result) : (result.previous(flvl, mlvl-flvl));
}

// Next index
/// \brief          Next index.
/// \details        Computes the next index in the tree from the specified 
///                 level to the specified number of refinements. A cyclic 
///                 operation is performed if the last index has been reached.
/// \param[in]      ilvl Index of the first level.
/// \param[in]      nref Number of refinements of the first level.
/// \return         Index of the next element.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::next(const unsigned int ilvl, const unsigned int nref) const
{
    const Type flvl = (ilvl <= std::integral_constant<Type, Bits/(Dimension+1)>::value) ? (Type(ilvl)) : (std::integral_constant<Type, Bits/(Dimension+1)>::value);
    const Type mlvl = (flvl+nref <= std::integral_constant<Type, Bits/(Dimension+1)>::value) ? (Type(flvl+nref)) : (std::integral_constant<Type, Bits/(Dimension+1)>::value);
    const Type clvl = level();
    const Type tlvl = (clvl <= mlvl) ? (clvl) : (mlvl);
    const Type shift = std::integral_constant<Type, Bits>::value-(tlvl*std::integral_constant<Type, Dimension+1>::value)-(tlvl < mlvl);
    const SimpleHyperOctreeIndex<Type, Dimension, Bits> result((_data & block<Type>(shift))+((std::integral_constant<Type, 1>::value << shift)*(clvl == tlvl)));
    return (result.level() >= flvl) ? (result) : (result.next(flvl, mlvl-flvl));
}
//--------------------------------------------------------------------------- //



//----------------------------------- CURVE --------------------------------- //
// Linear index
/// \brief          Linear index.
/// \details        Computes the linear index of the current level from the 
///                 tree.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Level (Level index type.)
/// \param[in]      ilvl (Recursion level index.)
/// \return         Linear index at the current level. 
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class Level, class> 
inline Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::linear(const Level ilvl) const
{
    return (!std::is_same<Level, std::true_type>::value) ? ((Condition) ? ((ilvl <= Mid) ? (linear<From, Mid>(ilvl)) : (linear<Mid+(Mid != To), To>(ilvl))) : (dtlc<Type, (Mid > 0) ? ((~comb<Type>(Dimension+1, Dimension)) << (Bits-(Mid*(Dimension+1)))) : (Mid), Dimension>(_data))) : (linear<From, To>(level()));
}
        
// Z-curve index
/// \brief          Z-curve index.
/// \details        Computes the Z-curve index of the current level from the 
///                 tree index by removing the interleaving opcodes.
/// \return         Z-curve index at the current level. 
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::zcurve() const
{
    return ((pext<Type, ((~comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>(_data)) >> ((std::integral_constant<Type, Bits/(Dimension+1)>::value-level())*(std::integral_constant<Type, Dimension>::value)));
}

// Site index
/// \brief          Site index.
/// \details        Computes the index of the site at the specified level, 
///                 which is the child index at this level. If no level is 
///                 specified, it returns the site index at the current level.
/// \tparam         Level (Level index type.)
/// \param[in]      ilvl Level index.
/// \return         Site index at the specified level. 
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <class Level, class> 
inline Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::site(const Level ilvl) const
{
    const Type shift = (std::is_same<Level, std::true_type>::value) ? (std::integral_constant<Type, Bits>::value-(level()*std::integral_constant<Type, Dimension+1>::value)) : (std::integral_constant<Type, Bits>::value-((ilvl*std::integral_constant<Type, Dimension+1>::value)*(Type(ilvl) <= std::integral_constant<Type, Bits/(Dimension+1)>::value)));
    return ((_data >> shift) & (std::integral_constant<Type, block<Type>(0, Dimension)>::value))*(shift < std::integral_constant<Type, Bits>::value);
}

// Coordinate index
/// \brief          Coordinate index.
/// \details        Computes the coordinate along the specified dimension at the
///                 current level.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      idim Dimension index.
/// \return         Monodimensional index at the current level.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class> 
inline Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::coordinate(const unsigned int idim) const
{
    return (Condition) ? ((idim <= Mid) ? (coordinate<From, Mid>(idim)) : (coordinate<Mid+(Mid != To), To>(idim))) : ((pext<Type, (comb<Type>(Dimension+1, Mid) << (Bits%(Dimension+1)))>(_data)) >> (std::integral_constant<Type, Bits/(Dimension+1)>::value-level()));
}
//--------------------------------------------------------------------------- //



//--------------------------------- POSITION -------------------------------- //
// Position
/// \brief          Position.
/// \details        Computes the position of the element from the tree index 
///                 along the specified dimension and using the provided 
///                 topology.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Reference (Reference position.)
/// \tparam         Level (Level index type.)
/// \param[in]      idim Dimension index.
/// \param[in]      ilvl (Recursion level index.)
/// \return         Position of the element.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class Reference, class Level, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::position(const unsigned int idim, const Level ilvl) const
{
    return (!std::is_same<Level, std::true_type>::value) ? ((Condition) ? ((idim <= Mid) ? (position<Kind, Position, Extent, From, Mid>(idim, ilvl)) : (position<Kind, Position, Extent, Mid+(Mid != To), To>(idim, ilvl))) : ((Kind(Reference::num)/Kind(Reference::den))+((Kind(Extent::num)/Kind(Extent::den))*(Kind((std::integral_constant<Type, 2>::value*((pext<Type, (comb<Type>(Dimension+1, Mid) << (Bits%(Dimension+1)))>(_data)) >> (std::integral_constant<Type, Bits/(Dimension+1)>::value-ilvl)))+std::integral_constant<Type, 1>::value)/Kind((std::integral_constant<Type, 1>::value) << (ilvl+std::integral_constant<Type, 1>::value)))))) : (position<Kind, Position, Extent, From, To>(idim, level()));
}

// Center
/// \brief          Center.
/// \details        Computes the position of the center of the element from 
///                 the tree index along the specified dimension and using the 
///                 provided topology.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Reference (Reference position.)
/// \tparam         Level (Level index type.)
/// \param[in]      idim Dimension index.
/// \param[in]      ilvl (Recursion level index.)
/// \return         Position of the center of the element.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class Reference, class Level, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::center(const unsigned int idim, const Level ilvl) const
{
    return (!std::is_same<Level, std::true_type>::value) ? ((Condition) ? ((idim <= Mid) ? (center<Kind, Position, Extent, From, Mid>(idim, ilvl)) : (center<Kind, Position, Extent, Mid+(Mid != To), To>(idim, ilvl))) : ((Kind(Reference::num)/Kind(Reference::den))+((Kind(Extent::num)/Kind(Extent::den))*(Kind((std::integral_constant<Type, 2>::value*((pext<Type, (comb<Type>(Dimension+1, Mid) << (Bits%(Dimension+1)))>(_data)) >> (std::integral_constant<Type, Bits/(Dimension+1)>::value-ilvl)))+std::integral_constant<Type, 1>::value)/Kind((std::integral_constant<Type, 1>::value) << (ilvl+std::integral_constant<Type, 1>::value)))))) : (center<Kind, Position, Extent, From, To>(idim, level()));
}

// Minimum
/// \brief          Minimum.
/// \details        Computes the position of the minimum of the element from 
///                 the tree index along the specified dimension and using the 
///                 provided topology.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Reference (Reference position.)
/// \tparam         Level (Level index type.)
/// \param[in]      idim Dimension index.
/// \param[in]      ilvl (Recursion level index.)
/// \return         Position of the minimum of the element.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class Reference, class Level, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::minimum(const unsigned int idim, const Level ilvl) const
{
    return (!std::is_same<Level, std::true_type>::value) ? ((Condition) ? ((idim <= Mid) ? (minimum<Kind, Position, Extent, From, Mid>(idim, ilvl)) : (minimum<Kind, Position, Extent, Mid+(Mid != To), To>(idim, ilvl))) : ((Kind(Reference::num)/Kind(Reference::den))+((Kind(Extent::num)/Kind(Extent::den))*(Kind((pext<Type, (comb<Type>(Dimension+1, Mid) << (Bits%(Dimension+1)))>(_data)) >> (std::integral_constant<Type, Bits/(Dimension+1)>::value-ilvl))/Kind(std::integral_constant<Type, 1>::value << ilvl))))) : (minimum<Kind, Position, Extent, From, To>(idim, level()));
}

// Maximum
/// \brief          Maximum.
/// \details        Computes the position of the maximum of the element from 
///                 the tree index along the specified dimension and using the 
///                 provided topology.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Reference (Reference position.)
/// \tparam         Level (Level index type.)
/// \param[in]      idim Dimension index.
/// \param[in]      ilvl (Recursion level index.)
/// \return         Position of the maximum of the element.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class Reference, class Level, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::maximum(const unsigned int idim, const Level ilvl) const
{
    return (!std::is_same<Level, std::true_type>::value) ? ((Condition) ? ((idim <= Mid) ? (maximum<Kind, Position, Extent, From, Mid>(idim, ilvl)) : (maximum<Kind, Position, Extent, Mid+(Mid != To), To>(idim, ilvl))) : ((Kind(Reference::num)/Kind(Reference::den))+((Kind(Extent::num)/Kind(Extent::den))*(Kind(std::integral_constant<Type, 1>::value+((pext<Type, (comb<Type>(Dimension+1, Mid) << (Bits%(Dimension+1)))>(_data)) >> (std::integral_constant<Type, Bits/(Dimension+1)>::value-ilvl)))/Kind(std::integral_constant<Type, 1>::value << ilvl))))) : (maximum<Kind, Position, Extent, From, To>(idim, level()));
}

// Extent
/// \brief          Extent.
/// \details        Computes the extent of the element from the tree index.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \return         Extent.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::extent() const
{
    return Kind(Extent::num)/(Kind(Extent::den)*Kind(std::integral_constant<Type, 1>::value << level())); // WARNING : Extent not taken into account here
}

// Length
/// \brief          Length.
/// \details        Computes the edge length of the element from the tree index.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \return         Element length.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::length() const
{
    return Kind(Extent::num)/(Kind(Extent::den)*Kind(std::integral_constant<Type, 1>::value << level()));
}

// Volume
/// \brief          Volume.
/// \details        Computes the volume of the element from the tree index.
/// \tparam         Kind Output type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \return         Element volume.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename Kind, class Position, class Extent, class> 
inline Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::volume() const
{
    return pow<Kind, Dimension>(Kind(Extent::num)/(Kind(Extent::den)*Kind(std::integral_constant<Type, 1>::value << level())));
}
//--------------------------------------------------------------------------- //



//---------------------------------- SPACE ---------------------------------- //
// Containing an index
/// \brief          Containing an index.
/// \details        Checks whether the index contains another index. It returns
///                 true if the index is a parent of the provided index at any
///                 less or equal level of refinement. 
/// \param[in]      idx Other index.
/// \return         True if the index contains the provided index, false 
///                 otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::containing(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const
{
    return (((std::integral_constant<Type, ((~comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data)) <= ((std::integral_constant<Type, ((comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data))) && (_data <= idx._data);
}

// Contained by an index
/// \brief          Contained by an index.
/// \details        Checks whether the index is contained by another index. It 
///                 returns true if the index is a child of the provided index
///                 at any greater or equal level of refinement.
/// \param[in]      idx Other index.
/// \return         True if the index is contained by the provided index, false
///                 otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::contained(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const
{
    return (((std::integral_constant<Type, ((~comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data)) <= ((std::integral_constant<Type, ((comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data))) && (idx._data <= _data);
}

// Intersecting an index
/// \brief          Intersecting an index.
/// \details        Checks whether the index is intersecting another index. It
///                 returns true if the index is a child or a parent of the
///                 provided index.
/// \param[in]      idx Other index.
/// \return         True if the index intersection is not null, false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::intersecting(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const
{
    return (((std::integral_constant<Type, ((~comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data)) <= ((std::integral_constant<Type, ((comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data)));
}

// Adjoining an index
/// \brief          Adjoining an index.
/// \details        Checks whether the index is adjoining another index. It 
///                 returns true if the index has the same parent as the 
///                 provided index.
/// \param[in]      idx Other index.
/// \return         True if the index and the provided index are adjoining, 
///                 false otherwise.
template <typename Type, unsigned int Dimension, unsigned int Bits>
inline bool SimpleHyperOctreeIndex<Type, Dimension, Bits>::adjoining(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx) const
{
    const Type msk = block<Type>(std::integral_constant<Type, Bits+Dimension+1>::value-(level()*std::integral_constant<Type, Dimension+1>::value));
    return ((!((std::integral_constant<Type, ((comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1)))>::value) & (_data ^ idx._data))) && ((_data & msk) == (idx._data & msk)));
}
//--------------------------------------------------------------------------- //



//---------------------------------- STREAM --------------------------------- //
// Output stream operator 
/// \brief          Output stream operator. 
/// \details        Displays the bits of the index.
/// \tparam         SelfType (Unsigned integer type.)
/// \tparam         SelfDimension (Number of dimensions.)
/// \tparam         SelfBits (Size of the type in bits.)
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side object.
/// \return         Output stream. 
template <typename SelfType, unsigned int SelfDimension, unsigned int SelfBits> 
std::ostream& operator<<(std::ostream& lhs, const SimpleHyperOctreeIndex<SelfType, SelfDimension, SelfBits>& rhs)
{
    static const unsigned int n = std::integral_constant<unsigned int, SelfBits-1>::value;
    static const unsigned int d = std::integral_constant<unsigned int, SelfDimension+1>::value;
    static const char czero = '0';
    static const char cone = '1';
    static const SelfType one = SelfType(1);
    const char separator = lhs.fill();
    for (unsigned int i = 0; i <= n; ++i) {
        (((i > 0) && (i%d == 0)) ? (lhs<<separator) : (lhs));
        lhs<<((rhs._data & (one<<(n-i))) ? (cone) : (czero));
        (((i < n) && (i%d == 0)) ? (lhs<<separator) : (lhs));
    }
    return lhs;
}
//--------------------------------------------------------------------------- //



//-------------------------------- COMPUTATION ------------------------------ //
// Index from linear index
/// \brief          Index from linear index.
/// \details        Creates an index from a linear index provided at a given 
///                 level of refinement.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      ilvl Level index.
/// \param[in]      ilinear Linear index.
/// \return         Index value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class> 
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::delinearize(const unsigned int ilvl, const Type ilinear)
{
    return (Condition) ? ((ilvl <= Mid) ? (delinearize<From, Mid>(ilvl, ilinear)) : (delinearize<Mid+(Mid != To), To>(ilvl, ilinear))) : (SimpleHyperOctreeIndex<Type, Dimension, Bits>(((itlc<Type, block<Type>(0, (Dimension+1)*Mid), Dimension+1>((ilinear) | (std::integral_constant<Type, block<Type>(Mid*Dimension, Mid)>::value))) << (std::integral_constant<Type, Bits-(Mid*(Dimension+1))>::value))*std::integral_constant<Type, (Mid > 0)>::value));
}

// Index from Z-curve
/// \brief          Index from Z-curve.
/// \details        Creates an index from a Z-curve provided at a given level of 
///                 refinement.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      ilvl Level index.
/// \param[in]      izcurve Z-curve index.
/// \return         Index value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class> 
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::dezcurvify(const unsigned int ilvl, const Type izcurve)
{
    return (Condition) ? ((ilvl <= Mid) ? (dezcurvify<From, Mid>(ilvl, izcurve)) : (dezcurvify<Mid+(Mid != To), To>(ilvl, izcurve))) : (SimpleHyperOctreeIndex<Type, Dimension, Bits>((pdep<Type, ((~comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1))) & block<Type>(Bits-(Mid*(Dimension+1)))>(izcurve)) | (std::integral_constant<Type, ((comb<Type>(Dimension+1, Dimension)) << (Bits%(Dimension+1))) & block<Type>(Bits-(Mid*(Dimension+1)))>::value))); 
}

// Index from navigation
/// \brief          Index from navigation.
/// \details        Produces an index from a serie of operands representing
///                 navigation through consecutive sites.
/// \tparam         Types (Unsigned integer types.)
/// \param[in]      isites Local site indices.
/// \return         Index value.
template <typename Type, unsigned int Dimension, unsigned int Bits>
template <typename... Types, class> 
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::navigate(Types&&... isites)
{
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((glue<Type, block<Type>(0, Dimension+1), Dimension+1, true>((std::forward<Types>(isites) | std::integral_constant<Type, (std::integral_constant<Type, 1>::value << Dimension)>::value)...)) & (std::integral_constant<Type, block<Type>(Bits%(Dimension+1))>::value));
}

// Index from coordinates
/// \brief          Index from coordinates.
/// \details        Makes an index from the provided coordinates along each 
///                 dimension at the specified level of refinement.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Types (Unsigned integer types.)
/// \param[in]      ilvl Level index.
/// \param[in]      icoords Integral coordinates along each dimension.
/// \return         Index value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <unsigned int From, unsigned int To, unsigned int Mid, bool Condition, typename... Types, class>
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::make(const unsigned int ilvl, Types&&... icoords)
{
    return (Condition) ? ((ilvl <= Mid) ? (make<From, Mid>(ilvl, std::forward<Types>(icoords)...)) : (make<Mid+(Mid != To), To>(ilvl, std::forward<Types>(icoords)...))) : (SimpleHyperOctreeIndex<Type, Dimension, Bits>(((itlc<Type, block<Type>(0, (Dimension+1)*Mid), Dimension+1>((glue<Type, block<Type>(0, Mid), Mid, false>(std::forward<Types>(icoords)...)) | (std::integral_constant<Type, block<Type>(Mid*Dimension, Mid)>::value))) << (std::integral_constant<Type, Bits-(Mid*(Dimension+1))>::value))*std::integral_constant<Type, (Mid > 0)>::value));
}

// Index from position
/// \brief          Index from position.
/// \details        Computes an index from the provided real positions along 
///                 each dimension at the specified level of refinement and 
///                 using the provided topology.
/// \tparam         Kind Position type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         From (Lower search bound.)
/// \tparam         To (Upper search bound.)
/// \tparam         Mid (Middle value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Reference (Reference position.)
/// \tparam         Kinds (Input types.)
/// \param[in]      ilvl Level index.
/// \param[in]      iposs Real positions along each dimension.
/// \return         Index value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Kind, class Position, class Extent, unsigned int From, unsigned int To, unsigned int Mid, bool Condition, class Reference, typename... Kinds, class>
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::compute(const unsigned int ilvl, Kinds&&... iposs)
{
    return (Condition) ? ((ilvl <= Mid) ? (compute<Kind, Position, Extent, From, Mid>(ilvl, std::forward<Kinds>(iposs)...)) : (compute<Kind, Position, Extent, Mid+(Mid != To), To>(ilvl, std::forward<Kinds>(iposs)...))) : (SimpleHyperOctreeIndex<Type, Dimension, Bits>(((itlc<Type, block<Type>(0, (Dimension+1)*Mid), Dimension+1>((glue<Type, block<Type>(0, Mid), Mid, false>(Type((Kind(Extent::den)/Kind(Extent::num))*(Kind(std::integral_constant<Type, (std::integral_constant<Type, 1>::value << std::integral_constant<Type, Mid>::value)>::value))*(Kind(std::forward<Kinds>(iposs))-(Kind(Reference::num)/Kind(Reference::den))))...)) | (std::integral_constant<Type, block<Type>(Mid*Dimension, Mid)>::value))) << (std::integral_constant<Type, Bits-(Mid*(Dimension+1))>::value))*std::integral_constant<Type, (Mid > 0)>::value));
}

// Index from position ciphering
/// \brief          Index from position ciphering.
/// \details        Ciphers an index from the provided real positions along 
///                 each dimension at the maximum level of refinement and 
///                 using the provided topology.
/// \tparam         Kind Position type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         Reference (Reference position.)
/// \tparam         Kinds (Input types.)
/// \param[in]      iposs Real positions along each dimension.
/// \return         Index value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Kind, class Position, class Extent, class Reference, typename... Kinds, class>
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::cipher(Kinds&&... iposs)
{
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((itlc<Type, block<Type>(0, (Dimension+1)*(Bits/(Dimension+1))), Dimension+1>((glue<Type, block<Type>(0, Bits/(Dimension+1)), Bits/(Dimension+1), false>(Type((Kind(Extent::den)/Kind(Extent::num))*(Kind(std::integral_constant<Type, (std::integral_constant<Type, 1>::value << std::integral_constant<Type, Bits/(Dimension+1)>::value)>::value))*(Kind(std::forward<Kinds>(iposs))-(Kind(Reference::num)/Kind(Reference::den))))...)) | (std::integral_constant<Type, block<Type>((Bits/(Dimension+1))*Dimension, Bits/(Dimension+1))>::value))) << (std::integral_constant<Type, Bits-((Bits/(Dimension+1))*(Dimension+1))>::value));
}

// Index from merging to lowest common ancestor
/// \brief          Index from merging to lowest common ancestor.
/// \details        Computes the index of the lowest common ancestor of a serie
///                 of indices.
/// \tparam         Indices (Indices types.)
/// \param[in]      idx Index value.
/// \param[in]      idcs Indices values.
/// \return         Lowest common ancestor index.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class... Indices> 
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::merge(const SimpleHyperOctreeIndex<Type, Dimension, Bits>& idx, Indices&&... idcs)
{
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>((idx._data) & (block<Type>((std::integral_constant<Type, Bits>::value)-(((std::integral_constant<Type, Bits>::value-tzcnt<Type>(~bhsmask<Type>(~((accumulate<Type, std::bit_and<Type> >(idx._data, std::forward<Indices>(idcs)._data...)) | (accumulate<Type, std::bit_and<Type> >(~idx._data, ~std::forward<Indices>(idcs)._data...))))))/std::integral_constant<Type, Dimension+1>::value)*(std::integral_constant<Type, Dimension+1>::value)))));
}
//--------------------------------------------------------------------------- //



//-------------------------------- PROPERTIES ------------------------------- //
// Type
/// \brief          Type. 
/// \details        Returns a copy of the default type value.
/// \return         Copy of the default type value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::type()
{
    return Type();
}

// Total number of types 
/// \brief          Total number of types. 
/// \details        Counts the number of internal integers. 
/// \return         One. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::types()
{
    return std::integral_constant<unsigned int, 1>::value;
}

// Number of dimensions
/// \brief          Number of dimensions. 
/// \details        Returns the number of spatial dimensions.
/// \return         Copy of the number of dimensions.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::dimension()
{
    return Dimension;
}

// Number of bits
/// \brief          Number of bits. 
/// \details        Returns the number of bits used for the index.
/// \return         Copy of the number of bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::bits()
{
    return Bits;
}

// Size
/// \brief          Size. 
/// \details        Returns the size of the type in bits.
/// \return         Copy of the size of the type in bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::size()
{
    return std::integral_constant<unsigned int, sizeof(Type)*std::numeric_limits<unsigned char>::digits>::value;
}

// Total number of instructions
/// \brief          Total number of instructions. 
/// \details        Returns the maximum number of instructions in the index.
/// \return         Maximum number of instructions.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::instructions()
{
    return std::integral_constant<unsigned int, Bits/(Dimension+1)>::value;
}

// Number of opcodes
/// \brief          Number of opcodes. 
/// \details        Returns the number of opcodes in the mini-language.
/// \return         Number of opcodes.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::opcodes()
{
    return std::integral_constant<unsigned int, (static_cast<unsigned int>(1) << static_cast<unsigned int>(1))>::value;
}

// Number of operands
/// \brief          Number of operands. 
/// \details        Returns the number of operands in the mini-language.
/// \return         Number of operands.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::operands()
{
    return std::integral_constant<unsigned int, (static_cast<unsigned int>(1) << Dimension)>::value;
}

// Remainder
/// \brief          Remainder. 
/// \details        Returns the number of unused bits at the end.
/// \return         Number of unused bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::remainder()
{
    return std::integral_constant<unsigned int, Bits%(Dimension+1)>::value;
}

// Total number of refinements
/// \brief          Total number of refinements. 
/// \details        Returns the maximum number of refinements in the index.
/// \return         Maximum number of refinements.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::refinements()
{
    return std::integral_constant<unsigned int, Bits/(Dimension+1)>::value;
}

// Number of sites.
/// \brief          Number of sites.
/// \details        Returns the total number of different sites which 
///                 corresponds to the number of children.
/// \return         Number of sites.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr unsigned int SimpleHyperOctreeIndex<Type, Dimension, Bits>::sites()
{
    return std::integral_constant<unsigned int, (static_cast<unsigned int>(1) << Dimension)>::value;
}

// Total number of indices
/// \brief          Total number of indices.
/// \details        Returns the total number of possible indices across a given
///                 number of levels. If no level is provided, it returns the 
///                 total number of possible indices accross all levels.
/// \tparam         Level (Level index type.)
/// \param[in]      ilvl Index of the first level.
/// \param[in]      nref Number of refinements of the first level.
/// \return         Number of indices.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class Level, class> 
constexpr Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::indices(const Level ilvl, const unsigned int nref)
{
    return (std::is_same<Level, std::true_type>::value) ? (std::integral_constant<Type, ((comb<Type>(Dimension)) & (block<Type>(0, ((Bits/(Dimension+1))*Dimension)+1)))>::value) : ((std::integral_constant<Type, ((comb<Type>(Dimension)) & (block<Type>(0, ((Bits/(Dimension+1))*Dimension)+1)))>::value) & (block<Type>(ilvl*std::integral_constant<Type, Dimension>::value, (nref+std::integral_constant<Type, 1>::value)*std::integral_constant<Type, Dimension>::value)));
}

// Monodimensional subdivisions
/// \brief          Monodimensional subdivisions.
/// \details        Returns the number of monodimensional subdivisions at a
///                 given level. If no level is provided, it returns the 
///                 number of monodimensional subdivisions at the greatest
///                 possible level.
/// \tparam         Level (Level index type.)
/// \param[in]      ilvl Level index.
/// \return         Number of subdivisions.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class Level, class> 
constexpr Type SimpleHyperOctreeIndex<Type, Dimension, Bits>::subdivisions(const Level ilvl)
{
    return (std::is_same<Level, std::true_type>::value) ? (std::integral_constant<Type, ((std::integral_constant<Type, 1>::value) << (Bits/(Dimension+1)))>::value) : (((std::integral_constant<Type, 1>::value) << (ilvl))*(static_cast<unsigned int>(ilvl) <= std::integral_constant<unsigned int, Bits/(Dimension+1)>::value));
}

// Invalid index
/// \brief          Invalid index.
/// \details        Produces an invalidated index.
/// \return         Invalid index.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
constexpr SimpleHyperOctreeIndex<Type, Dimension, Bits> SimpleHyperOctreeIndex<Type, Dimension, Bits>::invalid()
{
    return SimpleHyperOctreeIndex<Type, Dimension, Bits>(std::integral_constant<Type, block<Type>(0, Bits-1)>::value);
}
//--------------------------------------------------------------------------- //



//--------------------------------- HELPERS --------------------------------- //
// Apply function to tuple extraction
/// \brief          Apply function to tuple extraction.
/// \details        Applies the provided function to each element of the 
///                 specified tuple by extracting them recursively.
/// \tparam         Template Output type.
/// \tparam         Function (Function type.)
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in,out]  f Function object taking a single argument.
/// \param[in,out]  tuple Input tuple.
/// \param[in,out]  values (Extracted values.)
/// \return         Pack of function results.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class Template, class Function, class Tuple, typename... Kinds, class> 
constexpr Template SimpleHyperOctreeIndex<Type, Dimension, Bits>::apply(Function&& f, Tuple&& tuple, Kinds&&... values)
{
    return apply<Template>(std::forward<Function>(f), std::forward<Tuple>(tuple), std::forward<Kinds>(values)..., std::get<sizeof...(Kinds)>(std::forward<Tuple>(tuple)));
}

// Apply function to tuple
/// \brief          Apply function to tuple.
/// \details        Applies the provided function to each element of the 
///                 specified tuple.
/// \tparam         Template Output type.
/// \tparam         Function (Function type.)
/// \tparam         Tuple (Tuple type.)
/// \tparam         Kinds (Arithmetic types.)
/// \param[in,out]  f Function object taking a single argument.
/// \param[in,out]  tuple Input tuple.
/// \param[in,out]  values (Extracted values.)
/// \return         Pack of function results.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class Template, class Function, class Tuple, typename... Kinds, class> 
constexpr Template SimpleHyperOctreeIndex<Type, Dimension, Bits>::apply(Function&& f, const Tuple& tuple, Kinds&&... values)
{
    return Template{apply<typename std::conditional<static_cast<bool>(sizeof(tuple)), typename std::remove_cv<typename std::remove_reference<Kinds>::type>::type, typename std::remove_cv<typename std::remove_reference<Kinds>::type>::type>::type>(std::forward<Function>(f), std::forward<Kinds>(values))...};
}

// Apply function to value
/// \brief          Apply function to value.
/// \details        Applies the provided function to the specified value.
/// \tparam         Template Output type.
/// \tparam         Function (Function type.)
/// \tparam         Kind (Arithmetic type.)
/// \param[in,out]  f Function object taking a single argument.
/// \param[in,out]  value Input value.
/// \return         Function result.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <class Template, class Function, typename Kind, class> 
constexpr Template SimpleHyperOctreeIndex<Type, Dimension, Bits>::apply(Function&& f, Kind&& value)
{
    return Template(std::forward<Function>(f)(std::forward<Kind>(value)));
}
        
// Accumulate values from tuple
/// \brief          Accumulate values from tuple.
/// \details        Recursively applies the function object to consecutive 
///                 pairs of values, eventually accumulating them.
/// \tparam         Kind Arithmetic type.
/// \tparam         Operation Commutative binary function object.
/// \tparam         Index (Recursion index.)
/// \tparam         Tuple (Tuple of arithmetic types.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      tuple Input tuple.
/// \return         Accumulation of values from tuple.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Kind, class Operation, unsigned int Index, class Tuple, bool Condition, class>
constexpr Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::accumulate(Tuple&& tuple)
{
    return (Condition) ? (Operation()(std::get<Index>(std::forward<Tuple>(tuple)), accumulate<Kind, Operation, Index+Condition>(std::forward<Tuple>(tuple)))) : (std::get<Index>(std::forward<Tuple>(tuple)));
}

// Accumulate several values
/// \brief          Accumulate several values.
/// \details        Recursively applies the function object to consecutive 
///                 pairs of values, eventually accumulating them.
/// \tparam         Kind Arithmetic type.
/// \tparam         Operation Commutative binary function object.
/// \tparam         Kinds (Arithmetic types.)
/// \param[in]      value Input value.
/// \param[in]      values Input values.
/// \return         Accumulation of several values.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Kind, class Operation, typename... Kinds, class> 
constexpr Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::accumulate(const Kind value, const Kinds... values)
{
    return Operation()(value, accumulate<Kind, Operation>(values...));
}

// Accumulate a value
/// \brief          Accumulate a value.
/// \details        Accumulates a single value to end tail recursion, returning 
///                 the input value itself.
/// \tparam         Kind Arithmetic type.
/// \tparam         Operation Commutative binary function object.
/// \param[in]      value Input value.
/// \return         Accumulation of a value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Kind, class Operation, class> 
constexpr Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::accumulate(const Kind value)
{
    return value;
}
        
// Integral exponentiation
/// \brief          Integral exponentiation.
/// \details        Computes the integral exponentiation of the specified value.
/// \tparam         Kind Arithmetic type.
/// \tparam         Exponent Exponent value.
/// \tparam         One (Value of one.)
/// \tparam         Greater (Greater than one recursion branching.)
/// \tparam         Less (Less than zero recursion branching.)
/// \tparam         Equal (Equal to one recursion branching.)
/// \param[in]      value Input value.
/// \return         Exponentiated value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Kind, int Exponent, int One, bool Greater, bool Less, bool Equal, class>
constexpr Kind SimpleHyperOctreeIndex<Type, Dimension, Bits>::pow(const Kind value)
{
    return (Greater) ? (value*pow<Kind, (Exponent-One)*Greater>(value)) : ((Less) ? (Kind(One)/pow<Kind, (-Exponent)*Less>(value)) : ((Equal) ? (value) : (Kind(One))));
}

// Block mask
/// \brief          Block mask.
/// \details        Produces a block mask with a given number of bits set after
///                 a given position.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Index (Recursion index.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Ones (Complementary of zero.)
/// \tparam         Size (Total number of bits.)
/// \param[in]      location Position of the first bit set.
/// \param[in]      length Number of bit set.
/// \return         Value of the block mask.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Zero, Integer One, Integer Ones, Integer Size, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::block(const Integer location, const Integer length)
{
    return (location < Size) ? (((Ones*(Size <= length)) | (((length < Size) ? (One << length) : (Zero))-One)) << location) : (Zero);
}

// Periodic mask
/// \brief          Periodic mask.
/// \details        Produces a periodic mask with each period composed of ones 
///                 and zeros. This kind of mask is used a lot for bit twiddling
///                 hacks with a power of two period.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Index (Recursion index.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      period Value of the period.
/// \param[in]      offset Optional value of the offset.
/// \return         Value of the periodic mask.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Index, Integer Zero, Integer One, Integer Condition, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::periodic(const Integer period, const Integer offset)
{
    return (Condition) ? (((One << Index)*((((Index+period)-(offset%period))%period) < (period+(period%(One+One)))/(One+One))) | (periodic<Integer, Index+Condition>(period, offset))) : (Condition);
}

// Comb mask
/// \brief          Comb mask.
/// \details        Produces a comb mask with one bit set at the beginning of 
///                 each period.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Index (Recursion index.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Size (Total number of bits.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      period Value of the period.
/// \param[in]      offset Optional value of the offset.
/// \return         Value of the comb mask.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Index, Integer Zero, Integer One, Integer Size, Integer Condition, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::comb(const Integer period, const Integer offset)
{
    return (offset < Size)*(((Condition) ? (((One << Index)*(Index%period == Zero)) | (comb<Integer, Index+Condition>(period, offset))) : (Condition)) << ((Index == Zero)*(offset%period)));
}

// Next higher power of two
/// \brief          Next higher power of two.
/// \details        Computes recursively the power of two that is greater or
///                 equal to the specified value masked by the provided mask. 
///                 If the specified value is equal to zero or if the result 
///                 cannot be stored in the specified integer type, the 
///                 function returns zero.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Step (Recursion step.)
/// \tparam         One (Value of one.)
/// \tparam         Ones (Complementary of zero.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      value Input value.
/// \return         Power of two greater or equal to the specified value.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer One, Integer Ones, Integer Condition, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::nhp(const Integer value)
{
    return (Condition) ? (nhp<Integer, Ones, (Step << One)*Condition>(((Step != One) ? (value) : ((Mask == Ones) ? (value-One) : ((value & Mask)-One))) | (((Step != One) ? (value) : ((Mask == Ones) ? (value-One) : ((value & Mask)-One))) >> Step))) : (value+One);
}

// Mask up to highest set bit
/// \brief          Mask up to highest set bit.
/// \details        Computes a mask of bits set up to the highest set bit in 
///                 the specified value masked by the provided mask. If the 
///                 specified value is equal to zero, the function returns
///                 zero.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Step (Recursion step.)
/// \tparam         One (Value of one.)
/// \tparam         Ones (Complementary of zero.)
/// \tparam         Condition (Recursion branching.)
/// \param[in]      value Input value.
/// \return         Number with all bits set up to the highest set bit.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer One, Integer Ones, Integer Condition, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::bhsmask(const Integer value)
{
    return (Condition) ? (bhsmask<Integer, Ones, (Step << One)*Condition>(((Mask == Ones) ? (value) : (value & Mask)) | (((Mask == Ones) ? (value) : (value & Mask)) >> Step))) : (value); 
}

// Leading zero count
/// \brief          Leading zero count.
/// \details        Computes the number of leading zero bits in the specified 
///                 value masked by the provided mask. If the specified value 
///                 is equal to zero, the result is set to the number of bits 
///                 of the integer.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Step (Recursion step.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Ones (Complementary of zero.)
/// \tparam         Size (Total number of bits.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Number of leading zero bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer Zero, Integer One, Integer Ones, Integer Size, Integer Temporary, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::lzcnt(const Integer value)
{
    return (Step) ? ((Step != Size) ? ((Step*((value & Temporary) != Zero)) | (lzcnt<Integer, Ones, (Step >> One)>(value >> (Step*((value & Temporary) != Zero))))) : ((Mask == Ones) ? (Size-(value != Zero)-lzcnt<Integer, Ones, (Step >> One)>(value)) : (lzcnt<Integer, Ones, Step>(Mask & value)))) : (Step);
}

// Trailing zero count
/// \brief          Trailing zero count.
/// \details        Computes the number of trailing zero bits in the specified 
///                 value masked by the provided mask. If the specified value 
///                 is equal to zero, the result is set to the number of bits 
///                 of the integer.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Step (Recursion step.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Ones (Complementary of zero.)
/// \tparam         Size (Total number of bits.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Number of trailing zero bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer Zero, Integer One, Integer Ones, Integer Size, Integer Temporary, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::tzcnt(const Integer value)
{
    return (Step) ? ((Step != Size) ? (((!(Temporary & value)) ? (Zero) : (Step))+tzcnt<Integer, Ones, (Step >> One)>(value)) : ((Mask == Ones) ? (((value) ? (Size-One) : (Size))-tzcnt<Integer, Ones, (Step >> One)>(value & -value)) : (tzcnt<Integer, Ones, Step>(Mask & value)))) : (Step); 
}

// Population count
/// \brief          Population count.
/// \details        Computes the number of bit set in the specified value 
///                 masked by the provided mask. 
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Step (Recursion step.)
/// \tparam         One (Value of one.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Number of bit set.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer One, Integer Condition, Integer Temporary, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::popcnt(const Integer value)
{
    return (Condition) ? (((Temporary) && ((value >> Step) & Condition))+popcnt<Integer, Mask, Step+Condition>(value)) : (Condition);
}

// Parallel bits extract
/// \brief          Parallel bits extract.
/// \details        Extracts the bits of the specified value corresponding to 
///                 the bits set in the provided mask.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask for extraction.
/// \tparam         Step (Recursion step.)
/// \tparam         Shift (Recursion shift value.)
/// \tparam         One (Value of one.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Extracted bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer Shift, Integer One, Integer Condition, Integer Temporary, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::pext(const Integer value)
{
    return (Condition) ? ((!Temporary) ? (pext<Integer, Mask, Step+Condition, Shift+Temporary>(value)) : ((((value >> Step) & Condition) << Shift) | (pext<Integer, Mask, Step+Condition, Shift+Temporary>(value)))) : (Condition);
}

// Parallel bits deposit
/// \brief          Parallel bits deposit.
/// \details        Deposits the bits of the specified value corresponding to 
///                 the bits set in the provided mask.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask for deposit.
/// \tparam         Step (Recursion step.)
/// \tparam         Shift (Recursion shift value.)
/// \tparam         One (Value of one.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Deposited bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Step, Integer Shift, Integer One, Integer Condition, Integer Temporary, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::pdep(const Integer value)
{
    return (Condition) ? ((!Temporary) ? (pdep<Integer, Mask, Step+Condition, Shift+Temporary>(value)) : ((((value >> Shift) & Condition) << Step) | (pdep<Integer, Mask, Step+Condition, Shift+Temporary>(value)))) : (Condition);
}

// Interlace bits
/// \brief          Interlace bits.
/// \details        Interlaces bits of the specified value corresponding to 
///                 the bits set in the provided mask and using the provided
///                 period. For example if a period of N is specified, the bits
///                 of the original value masked by one are redistributed in 
///                 sequences of N bits.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask for interlacing.
/// \tparam         Period Interlacing period.
/// \tparam         Step (Recursion step.)
/// \tparam         Count (Recursion set bit count of the mask.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Size (Total number of bits.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Population (Maximized population count of the mask.)
/// \tparam         Destination (Extracted bit destination.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Interlaced bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Period, Integer Step, Integer Count, Integer Zero, Integer One, Integer Size, Integer Condition, Integer Population, Integer Destination, Integer Temporary, class>
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::itlc(const Integer value)
{
    return (Condition) ? ((Temporary) ? ((((value >> Step) & Condition) << Destination) | itlc<Integer, Mask, Period, Step+Condition, (Count+(((Step+Condition < Size) ? (Mask >> (Step+Condition)) : (Zero)) & Condition))*(Step >= tzcnt<Integer>(Mask))>(value)) : (itlc<Integer, Mask, Period, Step+Condition, (Count+(((Step+Condition < Size) ? (Mask >> (Step+Condition)) : (Zero)) & Condition))*(Step >= tzcnt<Integer>(Mask))>(value))) : (Condition);
}

// Deinterlace bits
/// \brief          Deinterlace bits.
/// \details        Deinterlaces bits of the specified value corresponding to 
///                 the bits set in the provided mask and using the provided
///                 period. For example if a period of N is specified, the bits
///                 of the original value masked by one are redistributed in N
///                 separate regions.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask for deinterlacing.
/// \tparam         Period Deinterlacing period.
/// \tparam         Step (Recursion step.)
/// \tparam         Count (Recursion set bit count of the mask.)
/// \tparam         Zero (Value of zero.)
/// \tparam         One (Value of one.)
/// \tparam         Size (Total number of bits.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Population (Maximized population count of the mask.)
/// \tparam         Destination (Extracted bit destination.)
/// \tparam         Temporary (Temporary mask.)
/// \param[in]      value Input value.
/// \return         Deinterlaced bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Period, Integer Step, Integer Count, Integer Zero, Integer One, Integer Size, Integer Condition, Integer Population, Integer Destination, Integer Temporary, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::dtlc(const Integer value)
{
    return (Condition) ? ((Temporary) ? ((((value >> Step) & Condition) << Destination) | dtlc<Integer, Mask, Period, Step+Condition, (Count+(((Step+Condition < Size) ? (Mask >> (Step+Condition)) : (Zero)) & Condition))*(Step >= tzcnt<Integer>(Mask))>(value)) : (dtlc<Integer, Mask, Period, Step+Condition, (Count+(((Step+Condition < Size) ? (Mask >> (Step+Condition)) : (Zero)) & Condition))*(Step >= tzcnt<Integer>(Mask))>(value))) : (Condition);
}

// Glue bit fields from tuple
/// \brief          Glue bit fields from tuple.
/// \details        Glues bit fields of the specified tuple one after another
///                 starting either from the least significant bit or the most
///                 significant bit. Each value is first masked using the 
///                 provided mask and the first remaining bits of the specified
///                 length are put one after another.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Length Bit field length.
/// \tparam         Msb Start from the most significant bit or not.
/// \tparam         Step (Recursion step.)
/// \tparam         Zero (Value of zero.)
/// \tparam         Direction (Shift direction.)
/// \tparam         Left (Recursion left shift value.)
/// \tparam         Right (Recursion right shift value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Tuple (Tuple of unsigned integer types.)
/// \tparam         Count (Count of remaining integers.)
/// \param[in]      tuple Input tuple.
/// \return         Glued bits from tuple.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Length, bool Msb, std::size_t Step, Integer Zero, Integer Direction, Integer Left, Integer Right, Integer Condition, class Tuple, Integer Count, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::glue(Tuple&& tuple)
{
    return (Condition) ? ((Count) ? (((Direction) ? ((Mask & Integer(std::get<Step>(std::forward<Tuple>(tuple)))) << Left) : ((Mask & Integer(std::get<Step>(std::forward<Tuple>(tuple)))) >> Right)) | (glue<Integer, Mask, Length, Msb, (Step+Condition)*(Count > Zero)>(std::forward<Tuple>(tuple)))) : ((Direction) ? ((Mask & Integer(std::get<Step>(std::forward<Tuple>(tuple)))) << Left) : ((Mask & Integer(std::get<Step>(std::forward<Tuple>(tuple)))) >> Right))) : (Zero);
}

// Glue bit fields from values
/// \brief          Glue bit fields from values.
/// \details        Glues bit fields of the specified values one after another
///                 starting either from the least significant bit or the most
///                 significant bit. Each value is first masked using the 
///                 provided mask and the first remaining bits of the specified
///                 length are put one after another.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Length Bit field length.
/// \tparam         Msb Start from the most significant bit or not.
/// \tparam         Step (Recursion step.)
/// \tparam         Zero (Value of zero.)
/// \tparam         Direction (Shift direction.)
/// \tparam         Left (Recursion left shift value.)
/// \tparam         Right (Recursion right shift value.)
/// \tparam         Condition (Recursion branching.)
/// \tparam         Integers (Unsigned integer types.)
/// \tparam         Count (Count of remaining integers.)
/// \param[in]      value Input value.
/// \param[in]      values Input values.
/// \return         Glued bits.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Length, bool Msb, std::size_t Step, Integer Zero, Integer Direction, Integer Left, Integer Right, Integer Condition, typename... Integers, Integer Count, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::glue(const Integer value, const Integers... values)
{
    return (Condition) ? ((Count) ? (((Direction) ? ((Mask & value) << Left) : ((Mask & value) >> Right)) | (glue<Integer, Mask, Length, Msb, Step+Condition>(values...))) : ((Direction) ? ((Mask & value) << Left) : ((Mask & value) >> Right))) : (Zero);
}

// Glue bit fields from nothing
/// \brief          Glue bit fields from nothing.
/// \details        Glues nothing to end tail recursion and returns zero.
/// \tparam         Integer Unsigned integer type.
/// \tparam         Mask Bit mask applied with bitwise and operation.
/// \tparam         Length Bit field length.
/// \tparam         Msb Start from the most significant bit or not.
/// \tparam         Step (Recursion step.)
/// \tparam         Zero (Value of zero.)
/// \return         Value of zero.
template <typename Type, unsigned int Dimension, unsigned int Bits> 
template <typename Integer, Integer Mask, Integer Length, bool Msb, std::size_t Step, Integer Zero, class> 
constexpr Integer SimpleHyperOctreeIndex<Type, Dimension, Bits>::glue()
{
    return Zero;
}
//--------------------------------------------------------------------------- //



//----------------------------------- TEST ---------------------------------- //
// Example function 
/// \brief          Example function. 
/// \details        Tests and demonstrates the use of SimpleHyperOctreeIndex. 
/// \return         0 if no error. 
template <typename Type, unsigned int Dimension, unsigned int Bits> 
int SimpleHyperOctreeIndex<Type, Dimension, Bits>::example()
{
    // Initialize
    std::cout<<"BEGIN = SimpleHyperOctreeIndex::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;

    // Construction
    SimpleHyperOctreeIndex<unsigned long long int, 3> i(4);
    SimpleHyperOctreeIndex<unsigned long long int, 3> j(8);
    SimpleHyperOctreeIndex<unsigned long long int, 2> k("1 01 1 01 1 11 1 10 1 11");
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"SimpleHyperOctreeIndex<unsigned long long int, 3>() : "                             <<SimpleHyperOctreeIndex<unsigned long long int, 3>()<<std::endl;
    std::cout<<std::setw(width*2)<<"SimpleHyperOctreeIndex<unsigned long long int, 3>(16) : "                           <<SimpleHyperOctreeIndex<unsigned long long int, 3>(16)<<std::endl;
    std::cout<<std::setw(width*2)<<"SimpleHyperOctreeIndex<unsigned char, 3>(23) : "                                    <<SimpleHyperOctreeIndex<unsigned char, 3>(23)<<std::endl;
    std::cout<<std::setw(width*2)<<"SimpleHyperOctreeIndex<unsigned long long int, 3>(i) : "                            <<SimpleHyperOctreeIndex<unsigned long long int, 3>(i)<<std::endl;
    std::cout<<std::setw(width*2)<<"SimpleHyperOctreeIndex<unsigned long long int, 3>(\"1 001 1 100 0\") : "            <<SimpleHyperOctreeIndex<unsigned long long int, 3>("1 001 1 100 0")<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"i = j : "                                     <<(i = j)<<std::endl;
    std::cout<<std::setw(width)<<"i = 42 : "                                    <<(i = 42)<<std::endl;
    std::cout<<std::setw(width)<<"i = \"1 001 1 100 0\" : "                     <<(i = "1 001 1 100 0")<<std::endl;
    std::cout<<std::setw(width)<<"i = 42 : "                                    <<(i = 42)<<std::endl;
    std::cout<<std::setw(width)<<"i == j : "                                    <<(i == j)<<std::endl;
    std::cout<<std::setw(width)<<"i() : "                                       <<i()<<std::endl;
    std::cout<<std::setw(width)<<"i[0] : "                                      <<i[0]<<std::endl;
    std::cout<<std::setw(width)<<"i[1] : "                                      <<i[1]<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Assignment : "                                <<std::endl;
    std::cout<<std::setw(width)<<"j.assign()"                                   <<j.assign()<<std::endl;
    std::cout<<std::setw(width)<<"j.assign(i)"                                  <<j.assign(i)<<std::endl;
    std::cout<<std::setw(width)<<"j.assign(42)"                                 <<j.assign(42)<<std::endl;
    std::cout<<std::setw(width)<<"j.assign(\"1 001 1 100 0\")"                  <<j.assign("1 001 1 100 0")<<std::endl;
    
    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"j.data() = 42"                                <<(j.data() = 42)<<std::endl;
    std::cout<<std::setw(width)<<"j.data()"                                     <<j.data()<<std::endl;
    std::cout<<std::setw(width)<<"j.get()"                                      <<j.get()<<std::endl;
    std::cout<<std::setw(width)<<"j.set()"                                      <<j.set(42)<<std::endl;
    std::cout<<std::setw(width)<<"j.nullify()"                                  <<j.nullify()<<std::endl;
    std::cout<<std::setw(width)<<"j.copy()"                                     <<j.copy()<<std::endl;
    std::cout<<std::setw(width)<<"j.cast()"                                     <<j.cast()<<std::endl;
    std::cout<<std::setw(width)<<"i.stringify(2)"                               <<i.stringify(2)<<std::endl;

    // Core
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Core : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"k.level()"                                    <<k.level()<<std::endl;
    std::cout<<std::setw(width)<<"k.coarsest()"                                 <<k.coarsest()<<std::endl;
    std::cout<<std::setw(width)<<"k.finest()"                                   <<k.finest()<<std::endl;
    std::cout<<std::setw(width)<<"k.limited()"                                  <<k.limited()<<std::endl;
    std::cout<<std::setw(width)<<"k.check()"                                    <<k.check()<<std::endl;
    std::cout<<std::setw(width)<<"i.invalidated()"                              <<i.invalidated()<<std::endl;
    std::cout<<std::setw(width)<<"k.fix()"                                      <<k.fix()<<std::endl;
    std::cout<<std::setw(width)<<"i.invalidate()"                               <<i.invalidate()<<std::endl;

    // Tree
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Tree : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"k.parent()"                                   <<k.parent()<<std::endl;
    std::cout<<std::setw(width)<<"k.child(4)"                                   <<k.child(4)<<std::endl;
    std::cout<<std::setw(width)<<"k.brother(4)"                                 <<k.brother(4)<<std::endl;
    std::cout<<std::setw(width)<<"k.preceding()"                                <<k.preceding()<<std::endl;
    std::cout<<std::setw(width)<<"k.following()"                                <<k.following()<<std::endl;
    std::cout<<std::setw(width)<<"k.previous(0, 8)"                             <<k.previous(0, 8)<<std::endl;
    std::cout<<std::setw(width)<<"k.next(0, 8)"                                 <<k.next(0, 8)<<std::endl;

    // Curve
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Curve : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"k.linear()"                                   <<k.linear()<<std::endl;
    std::cout<<std::setw(width)<<"k.zcurve()"                                   <<k.zcurve()<<std::endl;
    std::cout<<std::setw(width)<<"k.site(4)"                                    <<k.site(4)<<std::endl;
    std::cout<<std::setw(width)<<"k.coordinate(0)"                              <<k.coordinate(0)<<std::endl;

    // Position
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Position : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"k.position(0)"                                <<k.position(0)<<std::endl;
    std::cout<<std::setw(width)<<"k.center(0)"                                  <<k.center(0)<<std::endl;
    std::cout<<std::setw(width)<<"k.minimum(0)"                                 <<k.minimum(0)<<std::endl;
    std::cout<<std::setw(width)<<"k.maximum(0)"                                 <<k.maximum(0)<<std::endl;
    std::cout<<std::setw(width)<<"k.extent()"                                   <<k.extent()<<std::endl;
    std::cout<<std::setw(width)<<"k.length()"                                   <<k.length()<<std::endl;
    std::cout<<std::setw(width)<<"k.volume()"                                   <<k.volume()<<std::endl;

    // Space
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Space : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"i.containing(j) : "                           <<i.containing(j)<<std::endl;
    std::cout<<std::setw(width)<<"i.contained(j) : "                            <<i.contained(j)<<std::endl;
    std::cout<<std::setw(width)<<"i.intersecting(j) : "                         <<i.intersecting(j)<<std::endl;
    std::cout<<std::setw(width)<<"i.adjoining(j) : "                            <<i.adjoining(j)<<std::endl;
        
    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, i) : "                  <<i<<std::endl;

    // Computation
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Computation : "                               <<std::endl;
    std::cout<<std::setw(width)<<"i.delinearize(3, 424) : "                     <<i.delinearize(3, 424)<<std::endl;
    std::cout<<std::setw(width)<<"i.dezcurvify(3, 424) : "                      <<i.dezcurvify(3, 424)<<std::endl;
    std::cout<<std::setw(width)<<"i.navigate(0, 1, 2, 3, 4, 5, 6, 7) : "        <<i.navigate(0, 1, 2, 3, 4, 5, 6, 7)<<std::endl;
    std::cout<<std::setw(width)<<"i.make(4, 4, 8, 15) : "                       <<i.make(4, 4, 8, 15)<<std::endl;
    std::cout<<std::setw(width)<<"i.compute(4, 0.8, 0.15, 0.16) : "             <<i.compute(4, 0.8, 0.15, 0.16)<<std::endl;
    std::cout<<std::setw(width)<<"i.cipher(0.8, 0.15, 0.16) : "                 <<i.cipher(0.8, 0.15, 0.16)<<std::endl;
    std::cout<<std::setw(width)<<"i.merge(i, j) : "                             <<i.merge(i, j)<<std::endl;
    
    // Properties
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Properties : "                                <<std::endl;
    std::cout<<std::setw(width)<<"i.type() : "                                  <<i.type()<<std::endl;
    std::cout<<std::setw(width)<<"i.types() : "                                 <<i.types()<<std::endl;
    std::cout<<std::setw(width)<<"i.dimension() : "                             <<i.dimension()<<std::endl;
    std::cout<<std::setw(width)<<"i.bits() : "                                  <<i.bits()<<std::endl;
    std::cout<<std::setw(width)<<"i.size() : "                                  <<i.size()<<std::endl;
    std::cout<<std::setw(width)<<"i.instructions() : "                          <<i.instructions()<<std::endl;
    std::cout<<std::setw(width)<<"i.opcodes() : "                               <<i.opcodes()<<std::endl;
    std::cout<<std::setw(width)<<"i.operands() : "                              <<i.operands()<<std::endl;
    std::cout<<std::setw(width)<<"i.remainder() : "                             <<i.remainder()<<std::endl;
    std::cout<<std::setw(width)<<"i.refinements() : "                           <<i.refinements()<<std::endl;
    std::cout<<std::setw(width)<<"i.sites() : "                                 <<i.sites()<<std::endl;
    std::cout<<std::setw(width)<<"i.indices(4, 8) : "                           <<i.indices(4, 8)<<std::endl;
    std::cout<<std::setw(width)<<"i.subdivisions(15) : "                        <<i.subdivisions(15)<<std::endl;
    std::cout<<std::setw(width)<<"i.invalid() : "                               <<i.invalid()<<std::endl;
    
    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Helpers : "                                                                                                                 <<std::endl;
    std::cout<<std::setw(width*3)<<"std::get<0>(i.apply<std::tuple<int> >([](const int& x){return x*2;}, std::make_tuple(4))) : "                               <<std::get<0>(i.apply<std::tuple<int> >([](const int& x){return x*2;}, std::make_tuple(4)))<<std::endl;
    std::cout<<std::setw(width*3)<<"std::get<1>(i.apply<std::tuple<int, int> >([](const int& x){return x*2;}, std::make_tuple(4, 8))) : "                       <<std::get<1>(i.apply<std::tuple<int, int> >([](const int& x){return x*2;}, std::make_tuple(4, 8)))<<std::endl;
    std::cout<<std::setw(width*3)<<"i.apply<int>([](const int& x){return x*2;}, 4) : "                                                                          <<i.apply<int>([](const int& x){return x*2;}, 4)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.accumulate<unsigned int, std::plus<unsigned int> >(std::make_tuple(4, 8)) : "                                             <<i.accumulate<unsigned int, std::plus<unsigned int> >(std::make_tuple(4, 8))<<std::endl;
    std::cout<<std::setw(width*3)<<"i.accumulate<unsigned int, std::plus<unsigned int> >(4, 8, 15, 16, 23, 42) : "                                              <<i.accumulate<unsigned int, std::plus<unsigned int> >(4, 8, 15, 16, 23, 42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.accumulate<unsigned int, std::plus<unsigned int> >(42) : "                                                                <<i.accumulate<unsigned int, std::plus<unsigned int> >(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.pow<double, 2>(42) : "                                                                                                    <<i.pow<double, 2>(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.block<unsigned int>(3, 4) : "                                                                                             <<i.block<unsigned int>(3, 4)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.periodic<unsigned int>(4) : "                                                                                             <<i.periodic<unsigned int>(4)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.comb<unsigned int>(4, 3) : "                                                                                              <<i.comb<unsigned int>(4, 3)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.nhp<unsigned int>(42) : "                                                                                                 <<i.nhp<unsigned int>(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.bhsmask<unsigned int>(42) : "                                                                                             <<i.bhsmask<unsigned int>(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.lzcnt<unsigned int>(42) : "                                                                                               <<i.lzcnt<unsigned int>(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.tzcnt<unsigned int>(42) : "                                                                                               <<i.tzcnt<unsigned int>(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.popcnt<unsigned int>(42) : "                                                                                              <<i.popcnt<unsigned int>(42)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.pext<unsigned int, 42>(4242) : "                                                                                          <<i.pext<unsigned int, 42>(42424242)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.pdep<unsigned int, 42>(4242) : "                                                                                          <<i.pdep<unsigned int, 42>(42424242)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.itlc<unsigned int, 42, 2>(4242) : "                                                                                       <<i.itlc<unsigned int, 42, 2>(42424242)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.dtlc<unsigned int, 42, 2>(4242) : "                                                                                       <<i.dtlc<unsigned int, 42, 2>(42424242)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.glue<unsigned int, 63, 6, false>(std::make_tuple(4, 8, 15, 16, 23)) : "                                                   <<i.glue<unsigned int, 63, 6, false>(std::make_tuple(4, 8, 15, 16, 23))<<std::endl;
    std::cout<<std::setw(width*3)<<"i.glue<unsigned int, 63, 6, false>(4, 8, 15, 16, 23) : "                                                                    <<i.glue<unsigned int, 63, 6, false>(4, 8, 15, 16, 23)<<std::endl;
    std::cout<<std::setw(width*3)<<"i.glue<unsigned int, 63, 6, true>() : "                                                                                     <<i.glue<unsigned int, 63, 6, true>()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = SimpleHyperOctreeIndex::example()"<<std::endl;
    return 0;
}
//--------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // SIMPLEHYPEROCTREEINDEX_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
