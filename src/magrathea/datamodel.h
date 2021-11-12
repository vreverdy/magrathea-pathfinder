/* ******************************* DATAMODEL ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          DataModel
// DESCRIPTION :    Management of fundamental types representation
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           datamodel.h
/// \brief          Management of fundamental types representation
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef DATAMODEL_H_INCLUDED
#define DATAMODEL_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <stdexcept>
#include <limits>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Management of fundamental types representation
/// \brief          Management of fundamental types representation.
/// \details        Class to hold the data representation of fundamental types 
///                 on a system. The information is encoded in an <tt>unsigned 
///                 long long int</tt> in the following way where <tt>[BXbY-Z]
///                 </tt> means the information starts at the bit <tt>Y</tt> of
///                 the byte <tt>X</tt> and its size is <tt>Z</tt> bits : <ul>
///                 <li><tt>[B0b0-1]</tt> endianness</li>
///                 <li><tt>[B0b1-1]</tt> complement</li>
///                 <li><tt>[B1b0-1]</tt> <tt>float</tt> compatibility to the 
///                 IEEE-754 <tt>binary32</tt> representation</li>
///                 <li><tt>[B1b1-1]</tt> <tt>double</tt> compatibility to the 
///                 IEEE-754 <tt>binary64</tt> representation</li>
///                 <li><tt>[B1b2-1]</tt> <tt>long double</tt> compatibility to 
///                 the IEEE-754 <tt>binary128</tt> representation</li>
///                 <li><tt>[B2b0-4]</tt> <tt>void*</tt> size</li>
///                 <li><tt>[B2b4-4]</tt> <tt>bool</tt> size</li>
///                 <li><tt>[B3b0-4]</tt> <tt>char</tt> size</li>
///                 <li><tt>[B3b4-4]</tt> <tt>short int</tt> size</li>
///                 <li><tt>[B4b0-4]</tt> <tt>int</tt> size</li>
///                 <li><tt>[B4b4-4]</tt> <tt>long int</tt> size</li>
///                 <li><tt>[B5b0-4]</tt> <tt>long long int</tt> size</li>
///                 <li><tt>[B6b0-4]</tt> <tt>float</tt> size</li>
///                 <li><tt>[B6b4-4]</tt> <tt>double</tt> size</li>
///                 <li><tt>[B7b0-4]</tt> <tt>long double</tt> size</li></ul>
///                 Specified sizes cannot be equal to zero and are defaulted to
///                 one. Futhermore, the IEEE-754 compatibility corresponds to 
///                 correct byte size, IEC-559 compatibility, denormalization, 
///                 correct radix and correct number of mantissa digits.
class DataModel final
{
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        explicit inline DataModel(const unsigned long long int source = 0);
        explicit inline DataModel(const bool big, const bool twos, const bool fieee754, const bool dieee754, const bool ldieee754, const unsigned int psize, const unsigned int bsize, const unsigned int csize, const unsigned int sisize, const unsigned int isize, const unsigned int lisize, const unsigned int llisize, const unsigned int fsize, const unsigned int dsize, const unsigned int ldsize);
    //@}

    // Operators
    /// \name           Operators
    //@{
    public:
        inline unsigned long long int operator()();
        inline bool operator==(const DataModel& rhs);
        inline bool operator!=(const DataModel& rhs);
    //@}
    
    // Assignment
    /// \name           Assignment
    //@{
    public:
        inline DataModel& assign(const DataModel& source);
        inline DataModel& assign(const unsigned long long int source = 0);
        inline DataModel& assign(const bool big, const bool twos, const bool fieee754, const bool dieee754, const bool ldieee754, const unsigned int psize, const unsigned int bsize, const unsigned int csize, const unsigned int sisize, const unsigned int isize, const unsigned int lisize, const unsigned int llisize, const unsigned int fsize, const unsigned int dsize, const unsigned int ldsize);
    //@}
    
    // Management
    /// \name           Management
    //@{
    public:
        inline unsigned int size() const;
        inline const unsigned long long int& data() const;
        inline DataModel& clear();
        inline DataModel copy() const;
        template <typename Type = DataModel> inline Type cast() const;
        inline bool check();
    //@}
    
    // Getters
    /// \name           Getters
    //@{
    public:
        inline unsigned long long int get() const;
        inline bool endianness() const;
        inline bool complement() const;
        template <typename Type, class = typename std::enable_if<(!std::is_pointer<Type>::value) && (std::is_floating_point<typename std::decay<Type>::type>::value)>::type> inline bool ieee754() const;
        template <typename Type, class = typename std::enable_if<(std::is_pointer<Type>::value) || (std::is_integral<typename std::decay<Type>::type>::value) || (std::is_floating_point<typename std::decay<Type>::type>::value)>::type> inline unsigned int size() const; 
    //@}
    
    // Setters
    /// \name           Setters
    //@{
    public:
        inline DataModel& set(const bool big, const bool twos, const bool fieee754, const bool dieee754, const bool ldieee754, const unsigned int psize, const unsigned int bsize, const unsigned int csize, const unsigned int sisize, const unsigned int isize, const unsigned int lisize, const unsigned int llisize, const unsigned int fsize, const unsigned int dsize, const unsigned int ldsize);
        inline DataModel& endianness(const bool value);
        inline DataModel& complement(const bool value);
        template <typename Type, class = typename std::enable_if<(!std::is_pointer<Type>::value) && (std::is_floating_point<typename std::decay<Type>::type>::value)>::type> inline DataModel& ieee754(const bool value);
        template <typename Type, class = typename std::enable_if<(std::is_pointer<Type>::value) || (std::is_integral<typename std::decay<Type>::type>::value) || (std::is_floating_point<typename std::decay<Type>::type>::value)>::type> inline DataModel& size(const unsigned int value); 
    //@}
    
    // Stream
    /// \name           Stream
    //@{
    public:
        friend std::ostream& operator<<(std::ostream& lhs, const DataModel& rhs);
    //@}
    
    // Helpers
    /// \name           Helpers
    //@{
    public:
        template <typename Type, typename Decayed = typename std::decay<Type>::type, typename Signed = typename std::make_signed<typename std::conditional<(std::is_integral<Decayed>::value) && (!std::is_same<Decayed, bool>::value), Decayed, int>::type>::type, typename Fundamental = typename std::conditional<std::is_pointer<Type>::value, void*, typename std::conditional<(std::is_integral<Decayed>::value) && (!std::is_same<Decayed, bool>::value), typename std::conditional<std::is_same<Signed, signed char>::value, char, Signed>::type, Decayed>::type>::type, class = typename std::enable_if<(std::is_pointer<Fundamental>::value) || (std::is_integral<Fundamental>::value) || (std::is_floating_point<Fundamental>::value)>::type> static constexpr Fundamental fundamental();
        template <typename Type, class = typename std::enable_if<(!std::is_pointer<Type>::value) && (std::is_floating_point<typename std::decay<Type>::type>::value)>::type> static constexpr bool control754();
        static constexpr bool control();
    //@}
    
    // Predefined
    /// \name           Predefined
    //@{
    public:
        static inline const DataModel& system();
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
        unsigned long long int _code;                                           ///< Internal encoding of data types.
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Explicit value constructor
/// \brief          Explicit value constructor.
/// \details        Explicitely constructs the data model from an <tt>unsigned 
///                 long long int</tt> code.
/// \param[in]      source Code to be used for construction.
inline DataModel::DataModel(const unsigned long long int source)
: _code(source)
{
    ;
}

// Explicit detailed constructor
/// \brief          Explicit detailed constructor.
/// \details        Explicitely constructs the data model using all the needed 
///                 values.
/// \param[in]      big False for little-endian, true for big-endian.
/// \param[in]      twos True for two's complement, false otherwise.
/// \param[in]      fieee754 IEEE-754 <tt>binary32</tt> compatibility of 
///                 <tt>float</tt>.
/// \param[in]      dieee754 IEEE-754 <tt>binary64</tt> compatibility of 
///                 <tt>double</tt>.
/// \param[in]      ldieee754 IEEE-754 <tt>binary128</tt> compatibility of 
///                 <tt>long double</tt>.
/// \param[in]      psize Byte size of pointers.
/// \param[in]      bsize Byte size of <tt>bool</tt>.
/// \param[in]      csize Byte size of <tt>char</tt>.
/// \param[in]      sisize Byte size of <tt>short int</tt>.
/// \param[in]      isize Byte size of <tt>int</tt>.
/// \param[in]      lisize Byte size of <tt>long int</tt>.
/// \param[in]      llisize Byte size of <tt>long long int</tt>.
/// \param[in]      fsize Byte size of <tt>float</tt>.
/// \param[in]      dsize Byte size of <tt>double</tt>.
/// \param[in]      ldsize Byte size of <tt>long double</tt>.
inline DataModel::DataModel(const bool big, const bool twos, const bool fieee754, const bool dieee754, const bool ldieee754, const unsigned int psize, const unsigned int bsize, const unsigned int csize, const unsigned int sisize, const unsigned int isize, const unsigned int lisize, const unsigned int llisize, const unsigned int fsize, const unsigned int dsize, const unsigned int ldsize)
: _code(0)
{
    set(big, twos, fieee754, dieee754, ldieee754, psize, bsize, csize, sisize, isize, lisize, llisize, fsize, dsize, ldsize);
}
// -------------------------------------------------------------------------- //



// -------------------------------- OPERATORS ------------------------------- //
// Access operator
/// \brief          Access operator.
/// \details        Returns a copy of the underlying complete code of the data 
///                 model.
/// \return         Copy of the code.
inline unsigned long long int DataModel::operator()()
{
    return _code;
}

// Equal to comparison
/// \brief          Equal to comparison.
/// \details        Compares two data models for equality.
/// \return         True if the two data models are equal, false otherwise.
inline bool DataModel::operator==(const DataModel& rhs)
{
    return _code == rhs._code;
}

// Not equal to comparison
/// \brief          Not equal to comparison.
/// \details        Compares two data models for difference.
/// \return         True if the two data models are not equal, false otherwise.
inline bool DataModel::operator!=(const DataModel& rhs)
{
    return _code != rhs._code;
}
// -------------------------------------------------------------------------- //



// ------------------------------- ASSIGNMENT ------------------------------- //
// Copy assignment
/// \brief          Copy assignment.
/// \details        Assign the code from another data model.
/// \param[in]      source Source of the copy.
/// \return         Self reference.
inline DataModel& DataModel::assign(const DataModel& source)
{
    _code = source._code;
    return *this;
}

// Value assignment
/// \brief          Value assignment.
/// \details        Assigns a code to the data model.
/// \param[in]      source Source of the copy.
/// \return         Self reference.
inline DataModel& DataModel::assign(const unsigned long long int source)
{
    _code = source;
    return *this;
}

// Detailed assignment
/// \brief          Detailed assignment.
/// \details        Assigns the contents of the data model using all the needed 
///                 values.
/// \param[in]      big False for little-endian, true for big-endian.
/// \param[in]      twos True for two's complement, false otherwise.
/// \param[in]      fieee754 IEEE-754 <tt>binary32</tt> compatibility of 
///                 <tt>float</tt>.
/// \param[in]      dieee754 IEEE-754 <tt>binary64</tt> compatibility of 
///                 <tt>double</tt>.
/// \param[in]      ldieee754 IEEE-754 <tt>binary128</tt> compatibility of 
///                 <tt>long double</tt>.
/// \param[in]      psize Byte size of pointers.
/// \param[in]      bsize Byte size of <tt>bool</tt>.
/// \param[in]      csize Byte size of <tt>char</tt>.
/// \param[in]      sisize Byte size of <tt>short int</tt>.
/// \param[in]      isize Byte size of <tt>int</tt>.
/// \param[in]      lisize Byte size of <tt>long int</tt>.
/// \param[in]      llisize Byte size of <tt>long long int</tt>.
/// \param[in]      fsize Byte size of <tt>float</tt>.
/// \param[in]      dsize Byte size of <tt>double</tt>.
/// \param[in]      ldsize Byte size of <tt>long double</tt>.
/// \return         Self reference.
inline DataModel& DataModel::assign(const bool big, const bool twos, const bool fieee754, const bool dieee754, const bool ldieee754, const unsigned int psize, const unsigned int bsize, const unsigned int csize, const unsigned int sisize, const unsigned int isize, const unsigned int lisize, const unsigned int llisize, const unsigned int fsize, const unsigned int dsize, const unsigned int ldsize)
{
    return set(big, twos, fieee754, dieee754, ldieee754, psize, bsize, csize, sisize, isize, lisize, llisize, fsize, dsize, ldsize);
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Get the size of the code
/// \brief          Get the size of the code.
/// \details        Returns the result of the <tt>sizeof</tt> operator on the
///                 underlying code.
/// \return         Size of the code in bytes.
inline unsigned int DataModel::size() const
{
    return sizeof(_code);
}

// Access data
/// \brief          Access data.
/// \details        Returns a constant reference to the internal underlying data
///                 which is the data model code.
/// \return         Immutable reference to the code.
inline const unsigned long long int& DataModel::data() const
{
    return _code;
}

// Clear code
/// \brief          Clear code.
/// \details        Clears the whole contents and sets the flags to zero, and 
///                 the sizes to one.
/// \return         Self reference.
inline DataModel& DataModel::clear()
{
    _code = 0;
    return *this;
}

// Copy
/// \brief          Copy.
/// \details        Returns a copy of the data model.
/// \return         Copy.
inline DataModel DataModel::copy() const
{
    return *this;
}

// Cast
/// \brief          Cast.
/// \details        Returns a copy of the data model casted to the provided 
///                 type.
/// \tparam         Type Data type.
/// \return         Casted copy.
template <typename Type> 
inline Type DataModel::cast() const
{
    return Type(_code);
}

// Check if standard
/// \brief          Check if standard.
/// \details        Checks whether the data model is a standard one. It means
///                 that it has all the following properties : <ul>
///                 <li>it uses two's complement</li>
///                 <li><tt>float</tt> and <tt>double</tt> are IEEE-754
///                 compliant</li>
///                 <li>pointer size is 4 or 8</li>
///                 <li><tt>bool</tt> size is 1</li>
///                 <li><tt>char</tt> size is 1</li>
///                 <li><tt>short int</tt> size is 2</li>
///                 <li><tt>int</tt> size is 4</li>
///                 <li><tt>long int</tt> size is 4 or 8</li>
///                 <li><tt>long long int</tt> size is 8</li>
///                 <li><tt>float</tt> size is 4</li>
///                 <li><tt>double</tt> size is 8</li>
///                 <li><tt>long double</tt> size is 8, 10, 12 or 16</li></ul>
/// \return         True if the data model is standard, false if not.
inline bool DataModel::check()
{
    return ((complement()) && (ieee754<float>()) && (ieee754<double>()) && ((size<void*>() == 4) || (size<void*>() == 8)) && (size<bool>() == 1) && (size<char>() == 1) && (size<short int>() == 2) && (size<int>() == 4) && ((size<long int>() == 4) || (size<long int>() == 8)) && (size<long long int>() == 8) && (size<float>() == 4) && (size<double>() == 8) && ((size<long double>() == 8) || (size<long double>() == 10) || (size<long double>() == 12) || (size<long double>() == 16)));
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// Global getter
/// \brief          Global getter.
/// \details        Returns a copy of the underlying complete code of the data 
///                 model.
/// \return         Copy of the code.
inline unsigned long long int DataModel::get() const
{
    return _code;
}

// Get endianness
/// \brief          Get endianness.
/// \details        Returns the endianness from data model.
/// \return         False for little-endian, true for big-endian.
inline bool DataModel::endianness() const
{
    return (((*(reinterpret_cast<const unsigned char*>(&_code)+0))>>0)&1);
}

// Get complement
/// \brief          Get complement.
/// \details        Returns whether the system use two's complement encoding or
///                 not according to the data model.
/// \return         True for two's complement, false otherwise.
inline bool DataModel::complement() const
{
    return (((*(reinterpret_cast<const unsigned char*>(&_code)+0))>>1)&1);
}

// Get IEEE-754 compatibility
/// \brief          Get IEEE-754 compatibility.
/// \details        Returns whether the specified floating-point type is 
///                 compliant to the IEEE-754 standard according to the
///                 data model.
/// \tparam         Type Floating-point type.
/// \return         True if compliant, false otherwise. 
template <typename Type, class> 
inline bool DataModel::ieee754() const
{
    return ((std::is_same<decltype(fundamental<Type>()), float>::value) ? (((*(reinterpret_cast<const unsigned char*>(&_code)+1))>>0)&1)
         : ((std::is_same<decltype(fundamental<Type>()), double>::value) ? (((*(reinterpret_cast<const unsigned char*>(&_code)+1))>>1)&1)
         : ((std::is_same<decltype(fundamental<Type>()), long double>::value) ? (((*(reinterpret_cast<const unsigned char*>(&_code)+1))>>2)&1)
         : (0))));
}

// Get type size
/// \brief          Get type size.
/// \details        Returns the size of the provided fundamental type according
///                 to the data model. 
/// \tparam         Type Pointer, integral of floating-point type.
/// \return         Type size in bytes. 
template <typename Type, class> 
inline unsigned int DataModel::size() const
{
    return ((std::is_same<decltype(fundamental<Type>()), void*>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+2))>>0)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), bool>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+2))>>4)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), char>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+3))>>0)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), short int>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+3))>>4)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), int>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+4))>>0)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), long int>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+4))>>4)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), long long int>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+5))>>0)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), float>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+6))>>0)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), double>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+6))>>4)&((1<<4)-1))+1)
         : ((std::is_same<decltype(fundamental<Type>()), long double>::value) ? ((((*(reinterpret_cast<const unsigned char*>(&_code)+7))>>0)&((1<<4)-1))+1)
         : (0)))))))))));
}
// -------------------------------------------------------------------------- //



// --------------------------------- SETTERS -------------------------------- //
// Global setter
/// \brief          Global setter.
/// \details        Sets the content using all the needed values.
/// \param[in]      big False for little-endian, true for big-endian.
/// \param[in]      twos True for two's complement, false otherwise.
/// \param[in]      fieee754 IEEE-754 <tt>binary32</tt> compatibility of 
///                 <tt>float</tt>.
/// \param[in]      dieee754 IEEE-754 <tt>binary64</tt> compatibility of 
///                 <tt>double</tt>.
/// \param[in]      ldieee754 IEEE-754 <tt>binary128</tt> compatibility of 
///                 <tt>long double</tt>.
/// \param[in]      psize Byte size of pointers.
/// \param[in]      bsize Byte size of <tt>bool</tt>.
/// \param[in]      csize Byte size of <tt>char</tt>.
/// \param[in]      sisize Byte size of <tt>short int</tt>.
/// \param[in]      isize Byte size of <tt>int</tt>.
/// \param[in]      lisize Byte size of <tt>long int</tt>.
/// \param[in]      llisize Byte size of <tt>long long int</tt>.
/// \param[in]      fsize Byte size of <tt>float</tt>.
/// \param[in]      dsize Byte size of <tt>double</tt>.
/// \param[in]      ldsize Byte size of <tt>long double</tt>.
/// \return         Self reference.
inline DataModel& DataModel::set(const bool big, const bool twos, const bool fieee754, const bool dieee754, const bool ldieee754, const unsigned int psize, const unsigned int bsize, const unsigned int csize, const unsigned int sisize, const unsigned int isize, const unsigned int lisize, const unsigned int llisize, const unsigned int fsize, const unsigned int dsize, const unsigned int ldsize)
{
    return endianness(big).complement(twos).ieee754<float>(fieee754).ieee754<double>(dieee754).ieee754<long double>(ldieee754).size<void*>(psize).size<bool>(bsize).size<char>(csize).size<short int>(sisize).size<int>(isize).size<long int>(lisize).size<long long int>(llisize).size<float>(fsize).size<double>(dsize).size<long double>(ldsize);
}

// Set endianness
/// \brief          Set endianness.
/// \details        Sets the endianness of the data model.
/// \param[in]      value False for little-endian, true for big-endian.
/// \return         Self reference.
inline DataModel& DataModel::endianness(const bool value)
{
    ((*(reinterpret_cast<unsigned char*>(&_code)+0) &= ~(1<<0)) |= static_cast<unsigned char>(value)<<0);
    return *this;
}

// Set complement
/// \brief          Set complement.
/// \details        Sets the complement of the data model.
/// \param[in]      value True for two's complement, false otherwise.
/// \return         Self reference.
inline DataModel& DataModel::complement(const bool value)
{
    ((*(reinterpret_cast<unsigned char*>(&_code)+0) &= ~(1<<1)) |= static_cast<unsigned char>(value)<<1);
    return *this;
}

// Set IEEE-754 compatibility
/// \brief          Set IEEE-754 compatibility.
/// \details        Sets whether the specified floating-point type is 
///                 compliant to the IEEE-754 standard.
/// \tparam         Type Floating-point type.
/// \param[in]      value True if compliant, false otherwise.
/// \return         Self reference.
template <typename Type, class> 
inline DataModel& DataModel::ieee754(const bool value)
{
   (void)((std::is_same<decltype(fundamental<Type>()), float>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+1) &= ~(1<<0)) |= static_cast<unsigned char>(value)<<0)
       : ((std::is_same<decltype(fundamental<Type>()), double>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+1) &= ~(1<<1)) |= static_cast<unsigned char>(value)<<1)
       : ((std::is_same<decltype(fundamental<Type>()), long double>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+1) &= ~(1<<2)) |= static_cast<unsigned char>(value)<<2)
       : (0))));
    return *this;
}

// Set type size
/// \brief          Set type size.
/// \details        Sets the size of the provided fundamental type of the data
///                 model.
/// \tparam         Type Pointer, integral of floating-point type.
/// \param[in]      value Type size in bytes.
/// \return         Self reference. 
template <typename Type, class> 
inline DataModel& DataModel::size(const unsigned int value)
{
    (void)((std::is_same<decltype(fundamental<Type>()), void*>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+2) &= ~(((1<<4)-1)<<0)) |= static_cast<unsigned char>(value-(value > 0))<<0)
        : ((std::is_same<decltype(fundamental<Type>()), bool>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+2) &= ~(((1<<4)-1)<<4)) |= static_cast<unsigned char>(value-(value > 0))<<4)
        : ((std::is_same<decltype(fundamental<Type>()), char>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+3) &= ~(((1<<4)-1)<<0)) |= static_cast<unsigned char>(value-(value > 0))<<0)
        : ((std::is_same<decltype(fundamental<Type>()), short int>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+3) &= ~(((1<<4)-1)<<4)) |= static_cast<unsigned char>(value-(value > 0))<<4)
        : ((std::is_same<decltype(fundamental<Type>()), int>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+4) &= ~(((1<<4)-1)<<0)) |= static_cast<unsigned char>(value-(value > 0))<<0)
        : ((std::is_same<decltype(fundamental<Type>()), long int>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+4) &= ~(((1<<4)-1)<<4)) |= static_cast<unsigned char>(value-(value > 0))<<4)
        : ((std::is_same<decltype(fundamental<Type>()), long long int>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+5) &= ~(((1<<4)-1)<<0)) |= static_cast<unsigned char>(value-(value > 0))<<0)
        : ((std::is_same<decltype(fundamental<Type>()), float>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+6) &= ~(((1<<4)-1)<<0)) |= static_cast<unsigned char>(value-(value > 0))<<0)
        : ((std::is_same<decltype(fundamental<Type>()), double>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+6) &= ~(((1<<4)-1)<<4)) |= static_cast<unsigned char>(value-(value > 0))<<4)
        : ((std::is_same<decltype(fundamental<Type>()), long double>::value) ? ((*(reinterpret_cast<unsigned char*>(&_code)+7) &= ~(((1<<4)-1)<<0)) |= static_cast<unsigned char>(value-(value > 0))<<0)
        : (0)))))))))));
    return *this;
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Prints out the data model.
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side data model.
/// \return         Output stream.
std::ostream& operator<<(std::ostream& lhs, const DataModel& rhs)
{
    const char separator = lhs.fill();
    return lhs<<rhs.endianness()<<separator<<rhs.complement()<<separator<<rhs.ieee754<float>()<<separator<<rhs.ieee754<double>()<<separator<<rhs.ieee754<long double>()<<separator<<rhs.size<void*>()<<separator<<rhs.size<bool>()<<separator<<rhs.size<char>()<<separator<<rhs.size<short int>()<<separator<<rhs.size<int>()<<separator<<rhs.size<long int>()<<separator<<rhs.size<long long int>()<<separator<<rhs.size<float>()<<separator<<rhs.size<double>()<<separator<<rhs.size<long double>();
}
// -------------------------------------------------------------------------- //



// --------------------------------- HELPERS -------------------------------- //
// Get the associated fundamental type
/// \brief          Get the associated fundamental type.
/// \details        Converts the type passed as first template argument to an
///                 associated fundamental type : <ul>
///                 <li><tt>void*</tt> if the type is a pointer</li>
///                 <li><tt>bool</tt> if the type is a cv-qualified boolean</li>
///                 <li><tt>char</tt> if the type is a cv-qualified plain, 
///                 signed or unsigned character</li>
///                 <li>signed decayed type if the type is an integral type</li>
///                 <li>decayed type if the type is a floating-point type</li>
///                 </ul>.
/// \tparam         Type Type to convert.
/// \tparam         Decayed (Decayed version of the type.)
/// \tparam         Signed (Signed version of the type.)
/// \tparam         Fundamental (Fundamental version of the type.)
/// \return         Default value of the fundamental type.
template <typename Type, typename Decayed, typename Signed, typename Fundamental, class> 
constexpr Fundamental DataModel::fundamental()
{
    return Fundamental();
}

// Control IEEE-754 system compatibility
/// \brief          Control IEEE-754 system compatibility.
/// \details        Control whether the floating-point type is compatible with
///                 the IEEE-754 standard on the current architecture : correct 
///                 byte size, IEC-559 compatibility, denormalization, correct 
///                 radix and correct number of mantissa digits.
/// \tparam         Type Floating-point type.
/// \return         True if the type is compliant to IEEE-754, false otherwise.
template <typename Type, class> 
constexpr bool DataModel::control754()
{
    return (std::numeric_limits<Type>::is_iec559) && (std::numeric_limits<Type>::has_denorm) && (std::numeric_limits<Type>::radix == 2) && (((std::is_same<Type, float>::value) && (sizeof(Type) == 4) && (std::numeric_limits<Type>::digits == 24)) || ((std::is_same<Type, double>::value) && (sizeof(Type) == 8) && (std::numeric_limits<Type>::digits == 53)) || ((std::is_same<Type, long double>::value) && (sizeof(Type) == 16) && (std::numeric_limits<Type>::digits == 113)));
}

// Control if standard
/// \brief          Control if standard.
/// \details        Controls whether the system data model is a standard one. 
///                 It means that it has all the following properties : <ul>
///                 <li>it uses two's complement</li>
///                 <li><tt>float</tt> and <tt>double</tt> are IEEE-754
///                 compliant</li>
///                 <li><tt>long double</tt> is IEC-559 compliant, can be
///                 denormalized and has a binary radix</li>
///                 <li>pointer size is 4 or 8</li>
///                 <li><tt>bool</tt> size is 1</li>
///                 <li><tt>char</tt> size is 1</li>
///                 <li><tt>short int</tt> size is 2</li>
///                 <li><tt>int</tt> size is 4</li>
///                 <li><tt>long int</tt> size is 4 or 8</li>
///                 <li><tt>long long int</tt> size is 8</li>
///                 <li><tt>float</tt> size is 4</li>
///                 <li><tt>double</tt> size is 8</li>
///                 <li><tt>long double</tt> size is 8, 10, 12 or 16</li></ul>
/// \return         True if the data model is standard, false if not.
constexpr bool DataModel::control()
{
    return (-1 == ~0) && ((control754<float>()) && (control754<double>()) && ((std::numeric_limits<long double>::is_iec559) && (std::numeric_limits<long double>::has_denorm) && (std::numeric_limits<long double>::radix == 2))) && ((sizeof(void*) == 4) || (sizeof(void*) == 8)) && (sizeof(bool) == 1) && (sizeof(char) == 1) && (sizeof(short int) == 2) && (sizeof(int) == 4) && ((sizeof(long int) == 4) || (sizeof(long int) == 8)) && (sizeof(long long int) == 8) && (sizeof(float) == 4) && (sizeof(double) == 8) && ((sizeof(long double) == 8) || (sizeof(long double) == 10) || (sizeof(long double) == 12) || (sizeof(long double) == 16));
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// System data model
/// \brief          System data model.
/// \details        Returns an immutable reference to singleton representing the 
///                 data model of the current architecture.
/// \return         Immutable reference to system data model singleton.
inline const DataModel& DataModel::system()
{
    static const union {uint32_t i; char c[4];} x = {0x0000FEFF};
    static const DataModel singleton(!(x.c[0]), -1 == ~0, control754<float>(), control754<double>(), control754<long double>(), sizeof(void*), sizeof(bool), sizeof(char), sizeof(short int), sizeof(int), sizeof(long int), sizeof(long long int), sizeof(float), sizeof(double), sizeof(long double));
    return singleton;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of DataModel.
/// \return         0 if no error.
int DataModel::example()
{
    // Initialize
    std::cout<<"BEGIN = DataModel::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    DataModel model;
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"DataModel() : "                                                                     <<DataModel()<<std::endl;
    std::cout<<std::setw(width*2)<<"DataModel(0) : "                                                                    <<DataModel(0)<<std::endl;
    std::cout<<std::setw(width*2)<<"DataModel(0, 1, 1, 1, 0, 8, 1, 1, 2, 4, 4, 8, 4, 8, 16) : "                         <<DataModel(0, 1, 1, 1, 0, 8, 1, 1, 2, 4, 4, 8, 4, 8, 16)<<std::endl;
    std::cout<<std::setw(width*2)<<"DataModel(model) : "                                                                <<DataModel(model)<<std::endl;
    
    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"model() : "                                   <<model()<<std::endl;
    std::cout<<std::setw(width)<<"model == d.system() : "                       <<(model == model.system())<<std::endl;
    std::cout<<std::setw(width)<<"model != d.system() : "                       <<(model != model.system())<<std::endl;
    std::cout<<std::setw(width)<<"(model = model.system()) : "                  <<(model = model.system())<<std::endl;

    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Assignment : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"model.assign(0) : "                                                                 <<model.assign(0)<<std::endl;
    std::cout<<std::setw(width*2)<<"model.assign(model.system()) : "                                                    <<model.assign(model.system())<<std::endl;
    std::cout<<std::setw(width*2)<<"model.assign(0, 1, 1, 1, 0, 8, 1, 1, 2, 4, 4, 8, 4, 8, 16) : "                      <<model.assign(0, 1, 1, 1, 0, 8, 1, 1, 2, 4, 4, 8, 4, 8, 16)<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"model.size() : "                              <<model.size()<<std::endl;
    std::cout<<std::setw(width)<<"model.data() : "                              <<model.data()<<std::endl;
    std::cout<<std::setw(width)<<"model.clear() : "                             <<model.clear()<<std::endl;
    std::cout<<std::setw(width)<<"model.copy() : "                              <<model.copy()<<std::endl;
    std::cout<<std::setw(width)<<"model.cast() : "                              <<model.cast()<<std::endl;
    std::cout<<std::setw(width)<<"model.check() : "                             <<model.check()<<std::endl;

    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"model.get() : "                               <<model.get()<<std::endl;
    std::cout<<std::setw(width)<<"model.endianness() : "                        <<model.endianness()<<std::endl;
    std::cout<<std::setw(width)<<"model.complement() : "                        <<model.complement()<<std::endl;
    std::cout<<std::setw(width)<<"model.ieee754<float>() : "                    <<model.ieee754<float>()<<std::endl;
    std::cout<<std::setw(width)<<"model.size<long double>() : "                 <<model.size<long double>()<<std::endl;

    // Setters
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Setters : "                                                                         <<std::endl;
    std::cout<<std::setw(width*2)<<"model.set(0, 1, 1, 1, 0, 8, 1, 1, 2, 4, 4, 8, 4, 8, 16) : "                         <<model.set(0, 1, 1, 1, 0, 8, 1, 1, 2, 4, 4, 8, 4, 8, 16)<<std::endl;
    std::cout<<std::setw(width*2)<<"model.endianness(true) : "                                                          <<model.endianness(true)<<std::endl;
    std::cout<<std::setw(width*2)<<"model.complement(false) : "                                                         <<model.complement(false)<<std::endl;
    std::cout<<std::setw(width*2)<<"model.ieee754<float>(false) : "                                                     <<model.ieee754<float>(false)<<std::endl;
    std::cout<<std::setw(width*2)<<"model.size<long double>(10) : "                                                     <<model.size<long double>(10)<<std::endl;

    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, model) : "              <<model<<std::endl;

    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Helpers : "                                                                         <<std::endl;
    std::cout<<std::setw(width*2)<<"model.fundamental<volatile signed int>() :"                                         <<model.fundamental<volatile signed int>()<<std::endl;
    std::cout<<std::setw(width*2)<<"model.control754<float>() : "                                                       <<model.control754<float>()<<std::endl;
    std::cout<<std::setw(width*2)<<"model.control754<double>() : "                                                      <<model.control754<double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"model.control754<long double>() : "                                                 <<model.control754<long double>()<<std::endl;
    std::cout<<std::setw(width*2)<<"model.control() : "                                                                 <<model.control()<<std::endl;
    
    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Predefined : "                                <<std::endl;
    std::cout<<std::setw(width)<<"model.system() : "                            <<model.system()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = DataModel::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // DATAMODEL_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
