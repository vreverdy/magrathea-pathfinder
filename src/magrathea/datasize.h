/* ******************************** DATASIZE ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          DataSize
// DESCRIPTION :    Wrapper of binary data size and manager of unit conversion
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           datasize.h
/// \brief          Wrapper of binary data size and manager of unit conversion
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef DATASIZE_H_INCLUDED
#define DATASIZE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <stdexcept>
#include <ratio>
#include <limits>
#include <string>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Wrapper of binary data size and manager of unit conversion
/// \brief          Wrapper of binary data size and manager of unit conversion.
/// \details        Class to hold a byte count to represent the size of a file 
///                 or a chunk of memory. The byte count is stored internally 
///                 into a long long integer whose negative values are 
///                 associated with a void data size like a non-existing file. 
class DataSize final
{
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        constexpr DataSize();
        explicit constexpr DataSize(const long long int amount);
    //@}
    
    // Operators
    /// \name           Operators
    //@{
    public:
        template <typename Type> inline DataSize& operator=(const Type& rhs);
        inline DataSize& operator=(const DataSize& rhs);
        inline unsigned long long int operator()();
    //@}

    // Managemement
    /// \name           Management
    //@{
    public:
        constexpr const long long int& data() const;
        constexpr DataSize copy() const;
        template <typename Type = DataSize> constexpr Type cast() const;
    //@}
    
    // Getters
    /// \name           Getters
    //@{
    public:
        bool constexpr valid() const;
        bool constexpr empty() const;
        template <typename Type = unsigned long long int> constexpr Type bytes() const;
        template <typename Type = unsigned long long int> constexpr Type kibibytes() const;
        template <typename Type = unsigned long long int> constexpr Type mebibytes() const;
        template <typename Type = unsigned long long int> constexpr Type gibibytes() const;
        template <typename Type = unsigned long long int> constexpr Type tebibytes() const;
        template <typename Type = unsigned long long int> constexpr Type pebibytes() const;
        template <typename Type = unsigned long long int> constexpr Type exbibytes() const;
        template <typename Type = unsigned long long int> constexpr Type kilobytes() const;
        template <typename Type = unsigned long long int> constexpr Type megabytes() const;
        template <typename Type = unsigned long long int> constexpr Type gigabytes() const;
        template <typename Type = unsigned long long int> constexpr Type terabytes() const;
        template <typename Type = unsigned long long int> constexpr Type petabytes() const;
        template <typename Type = unsigned long long int> constexpr Type exabytes() const;
        template <typename Type = unsigned long long int> constexpr Type bits() const;
    //@}

    // Setters
    /// \name           Setters
    //@{
    public:
        inline DataSize& valid(const bool ok);
        inline DataSize& empty(const bool ok);
        template <typename Type> inline DataSize& bytes(const Type& amount);
        template <typename Type> inline DataSize& kibibytes(const Type& amount);
        template <typename Type> inline DataSize& mebibytes(const Type& amount);
        template <typename Type> inline DataSize& gibibytes(const Type& amount);
        template <typename Type> inline DataSize& tebibytes(const Type& amount);
        template <typename Type> inline DataSize& pebibytes(const Type& amount);
        template <typename Type> inline DataSize& exbibytes(const Type& amount);
        template <typename Type> inline DataSize& kilobytes(const Type& amount);
        template <typename Type> inline DataSize& megabytes(const Type& amount);
        template <typename Type> inline DataSize& gigabytes(const Type& amount);
        template <typename Type> inline DataSize& terabytes(const Type& amount);
        template <typename Type> inline DataSize& petabytes(const Type& amount);
        template <typename Type> inline DataSize& exabytes(const Type& amount);
        template <typename Type> inline DataSize& bits(const Type& amount);
    //@}
        
    // Predefined
    /// \name           Predefined
    //@{
    public:
        template <typename Type = unsigned long long int> static constexpr DataSize byte(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize kibi(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize mebi(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize gibi(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize tebi(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize pebi(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize exbi(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize kilo(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize mega(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize giga(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize tera(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize peta(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize exa(const Type& amount = Type(1));
        template <typename Type = unsigned long long int> static constexpr DataSize bit(const Type& amount = Type(8));
    //@}

    // Stream
    /// \name           Stream
    //@{
    public:
        friend std::ostream& operator<<(std::ostream& lhs, const DataSize& rhs);
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
        long long int _size;                                                    ///< Data size in bytes.
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit constexpr empty constructor
/// \brief          Implicit constexpr empty constructor.
/// \details        Implicitely constructs the data size and set it to an 
///                 invalid size.
constexpr DataSize::DataSize()
: _size(-1LL)
{
    ;
}

// Explicit constexpr value constructor
/// \brief          Explicit constexpr value constructor.
/// \details        Explicitely constructs the data size from a long long
///                 integer.
/// \param[in]      amount Value to be used for construction.
constexpr DataSize::DataSize(const long long int amount)
: _size((amount >= 0) ? (static_cast<long long int>(amount)) : (-1LL))
{
    ;
}
// -------------------------------------------------------------------------- //



// -------------------------------- OPERATORS ------------------------------- //
// Conversion assignment operator
/// \brief          Conversion assignment operator.
/// \details        Convert the provided amount of bytes to a data size.
/// \tparam         Type (Data type.)
/// \param[in]      rhs Right-hand side.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::operator=(const Type& rhs)
{
    _size = rhs;
    return *this;
}

// Copy assignment operator
/// \brief          Copy assignment operator.
/// \details        Copies the contents of another data size.
/// \param[in]      rhs Right-hand side.
/// \return         Self reference.
inline DataSize& DataSize::operator=(const DataSize& rhs)
{
    _size = rhs._size;
    return *this;
}

// Conversion operator
/// \brief          Conversion operator.
/// \details        Returns the number of bytes or throw an error if the size
///                 is not defined.
/// \return         The size in bytes.
/// \exception      std::underflow_error Undefined size.
inline unsigned long long int DataSize::operator()()
{
    if (_size < 0) {
        throw std::underflow_error("ERROR = DataSize::operator()() : undefined size");
    }
    return _size;
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Access data
/// \brief          Access data.
/// \details        Returns a constant reference to the internal underlying data 
///                 which is the size in bytes as a long long integer, negative 
///                 if undefined.
/// \return         Const reference to the data.
constexpr const long long int& DataSize::data() const
{
    return _size;
}

// Copy
/// \brief          Copy.
/// \details        Returns a copy of the data size.
/// \return         Copy.
constexpr DataSize DataSize::copy() const
{
    return DataSize(_size);
}

// Cast
/// \brief          Cast.
/// \details        Returns a copy of the data size casted to the provided type.
///                 If the size is undefined and the type has a signaling NaN, 
///                 then this signaling NaN is returned. 
/// \tparam         Type Data type.
/// \return         Copy.
template <typename Type> 
constexpr Type DataSize::cast() const
{
    return ((_size >= 0) || (!std::numeric_limits<Type>::has_signaling_NaN)) ? (Type(_size)) : (Type((std::numeric_limits<Type>::has_signaling_NaN) ? (std::numeric_limits<typename std::conditional<std::numeric_limits<Type>::has_signaling_NaN, Type, double>::type>::signaling_NaN()) : (typename std::conditional<std::numeric_limits<Type>::has_signaling_NaN, Type, double>::type(_size))));
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// Get the validity of the data size
/// \brief          Get the validity of the data size.
/// \details        Returns whether the data exists or not. The existence of
///                 data is defined as a size greater or equal to zero.
/// \return         Data validity.
constexpr bool DataSize::valid() const
{
    return (_size >= 0);
}

// Get whether the size is not null
/// \brief          Get whether the size is not null.
/// \details        Returns true if the size is equal to zero or if data is not
///                 valid.
/// \return         Whether the size is not strictly greater than zero.
constexpr bool DataSize::empty() const
{
    return (_size <= 0);
}

// Get data size in bytes
/// \brief          Get data size in bytes.
/// \details        Returns the data size in bytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of bytes.
template <typename Type> 
constexpr Type DataSize::bytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)) : (static_cast<Type>(0));
}

// Get data size in kibibytes
/// \brief          Get data size in kibibytes.
/// \details        Returns the data size in kibibytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of kibibytes.
template <typename Type> 
constexpr Type DataSize::kibibytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(1LL<<10LL)) : (static_cast<Type>(0));
}

// Get data size in mebibytes
/// \brief          Get data size in mebibytes.
/// \details        Returns the data size in mebibytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of mebibytes.
template <typename Type> 
constexpr Type DataSize::mebibytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(1LL<<20LL)) : (static_cast<Type>(0));
}

// Get data size in gibibytes
/// \brief          Get data size in gibibytes.
/// \details        Returns the data size in gibibytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of gibibytes.
template <typename Type> 
constexpr Type DataSize::gibibytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(1LL<<30LL)) : (static_cast<Type>(0));
}

// Get data size in tebibytes
/// \brief          Get data size in tebibytes.
/// \details        Returns the data size in tebibytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of tebibytes.
template <typename Type> 
constexpr Type DataSize::tebibytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(1LL<<40ULL)) : (static_cast<Type>(0));
}

// Get data size in pebibytes
/// \brief          Get data size in pebibytes.
/// \details        Returns the data size in pebibytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of pebibytes.
template <typename Type> 
constexpr Type DataSize::pebibytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(1LL<<50LL)) : (static_cast<Type>(0));
}

// Get data size in exbibytes
/// \brief          Get data size in exbibytes.
/// \details        Returns the data size in exbibytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of exbibytes.
template <typename Type> 
constexpr Type DataSize::exbibytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(1LL<<60LL)) : (static_cast<Type>(0));
}

// Get data size in kilobytes
/// \brief          Get data size in kilobytes.
/// \details        Returns the data size in kilobytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of kilobytes.
template <typename Type> 
constexpr Type DataSize::kilobytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(std::kilo::num)) : (static_cast<Type>(0));
}

// Get data size in megabytes
/// \brief          Get data size in megabytes.
/// \details        Returns the data size in megabytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of megabytes.
template <typename Type> 
constexpr Type DataSize::megabytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(std::mega::num)) : (static_cast<Type>(0));
}

// Get data size in gigabytes
/// \brief          Get data size in gigabytes.
/// \details        Returns the data size in gigabytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of gigabytes.
template <typename Type> 
constexpr Type DataSize::gigabytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(std::giga::num)) : (static_cast<Type>(0));
}

// Get data size in terabytes
/// \brief          Get data size in terabytes.
/// \details        Returns the data size in terabytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of terabytes.
template <typename Type> 
constexpr Type DataSize::terabytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(std::tera::num)) : (static_cast<Type>(0));
}

// Get data size in petabytes
/// \brief          Get data size in petabytes.
/// \details        Returns the data size in petabytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of petabytes.
template <typename Type> 
constexpr Type DataSize::petabytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(std::peta::num)) : (static_cast<Type>(0));
}

// Get data size in exabytes
/// \brief          Get data size in exabytes.
/// \details        Returns the data size in exabytes in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of exabytes.
template <typename Type> 
constexpr Type DataSize::exabytes() const
{
    return (_size > 0) ? (static_cast<Type>(_size)/static_cast<Type>(std::exa::num)) : (static_cast<Type>(0));
}

// Get data size in bits
/// \brief          Get data size in bits.
/// \details        Returns the data size in bits in the provided type.
/// \tparam         Type Data type.
/// \return         The amount of bits.
template <typename Type> 
constexpr Type DataSize::bits() const
{
    return (_size > 0) ? (static_cast<Type>(_size)*static_cast<Type>(8)) : (static_cast<Type>(0));
}
// -------------------------------------------------------------------------- //



// --------------------------------- SETTERS -------------------------------- //
// Set valid status
/// \brief          Set valid status.
/// \details        Sets the size to invalid if the argument is false, change it
///                 to the maximum between the current size and zero if true.
/// \param[in]      ok Status.
/// \return         Self reference.  
inline DataSize& DataSize::valid(const bool ok)
{
    _size = (ok) ? ((_size >= 0) ? (_size) : (0LL)) : (-1LL);
    return *this;
}

// Set empty status
/// \brief          Set empty status.
/// \details        Sets the size to zero if the argument is true, change it to
///                 the maximum between the current size and one byte if false.
/// \param[in]      ok Status.
/// \return         Self reference.  
inline DataSize& DataSize::empty(const bool ok)
{
    _size = (ok) ? (0LL) : ((_size > 0) ? (_size) : (1LL));
    return *this;
}

// Set data size in bytes
/// \brief          Set data size in bytes.
/// \details        Sets the current data size in bytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of bytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::bytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount)) : (-1LL);
    return *this;
}

// Set data size in kibibytes
/// \brief          Set data size in kibibytes.
/// \details        Sets the current data size in kibibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of kibibytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::kibibytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*(1LL<<10LL))) : (-1LL);
    return *this;
}

// Set data size in mebibytes
/// \brief          Set data size in mebibytes.
/// \details        Sets the current data size in mebibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of mebibytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::mebibytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*(1LL<<20LL))) : (-1LL);
    return *this;
}

// Set data size in gibibytes
/// \brief          Set data size in gibibytes.
/// \details        Sets the current data size in gibibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of gibibytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::gibibytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*(1LL<<30LL))) : (-1LL);
    return *this;
}

// Set data size in tebibytes
/// \brief          Set data size in tebibytes.
/// \details        Sets the current data size in tebibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of tebibytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::tebibytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*(1LL<<40LL))) : (-1LL);
    return *this;
}

// Set data size in pebibytes
/// \brief          Set data size in pebibytes.
/// \details        Sets the current data size in pebibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of pebibytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::pebibytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*(1LL<<50LL))) : (-1LL);
    return *this;
}

// Set data size in exbibytes
/// \brief          Set data size in exbibytes.
/// \details        Sets the current data size in exbibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of exbibytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::exbibytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*(1LL<<60LL))) : (-1LL);
    return *this;
}

// Set data size in kilobytes
/// \brief          Set data size in kilobytes.
/// \details        Sets the current data size in kilobytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of kilobytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::kilobytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*std::kilo::num)) : (-1LL);
    return *this;
}

// Set data size in megabytes
/// \brief          Set data size in megabytes.
/// \details        Sets the current data size in megabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of megabytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::megabytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*std::mega::num)) : (-1LL);
    return *this;
}

// Set data size in gigabytes
/// \brief          Set data size in gigabytes.
/// \details        Sets the current data size in gigabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of gigabytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::gigabytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*std::giga::num)) : (-1LL);
    return *this;
}

// Set data size in terabytes
/// \brief          Set data size in terabytes.
/// \details        Sets the current data size in terabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of terabytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::terabytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*std::tera::num)) : (-1LL);
    return *this;
}

// Set data size in petabytes
/// \brief          Set data size in petabytes.
/// \details        Sets the current data size in petabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of petabytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::petabytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*std::peta::num)) : (-1LL);
    return *this;
}

// Set data size in exabytes
/// \brief          Set data size in exabytes.
/// \details        Sets the current data size in exabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of exabytes.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::exabytes(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount*std::exa::num)) : (-1LL);
    return *this;
}

// Set data size in bits
/// \brief          Set data size in bits.
/// \details        Sets the current data size in bits. If the number of 
///                 provided bits is not divisible by 8, an extra byte is 
///                 added.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of bits.
/// \return         Self reference.
template <typename Type> 
inline DataSize& DataSize::bits(const Type& amount)
{
    _size = (amount >= 0) ? (static_cast<long long int>(amount)/8LL+static_cast<bool>(static_cast<long long int>(amount)%8LL)) : (-1LL);
    return *this;
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Predefined data size in bytes
/// \brief          Predefined data size in bytes.
/// \details        Constructs and returns a data size based on an amount of
///                 bytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of bytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::byte(const Type& amount)
{
    return DataSize(amount);
}

// Predefined data size in kibibytes
/// \brief          Predefined data size in kibibytes.
/// \details        Constructs and returns a data size based on an amount of
///                 kibibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of kibibytes.
/// \return         Copy of the data size.
template <typename Type>
constexpr DataSize DataSize::kibi(const Type& amount)
{
    return DataSize(amount*(1LL<<10LL));
}

// Predefined data size in mebibytes
/// \brief          Predefined data size in mebibytes.
/// \details        Constructs and returns a data size based on an amount of
///                 mebibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of mebibytes.
/// \return         Copy of the data size.
template <typename Type>
constexpr DataSize DataSize::mebi(const Type& amount)
{
    return DataSize(amount*(1LL<<20ULL));
}

// Predefined data size in gibibytes
/// \brief          Predefined data size in gibibytes.
/// \details        Constructs and returns a data size based on an amount of
///                 gibibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of gibibytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::gibi(const Type& amount)
{
    return DataSize(amount*(1LL<<30LL));
}

// Predefined data size in tebibytes
/// \brief          Predefined data size in tebibytes.
/// \details        Constructs and returns a data size based on an amount of
///                 tebibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of tebibytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::tebi(const Type& amount)
{
    return DataSize(amount*(1LL<<40LL));
}

// Predefined data size in pebibytes
/// \brief          Predefined data size in pebibytes.
/// \details        Constructs and returns a data size based on an amount of
///                 pebibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of pebibytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::pebi(const Type& amount)
{
    return DataSize(amount*(1LL<<50LL));
}

// Predefined data size in exbibytes
/// \brief          Predefined data size in exbibytes.
/// \details        Constructs and returns a data size based on an amount of
///                 exbibytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of exbibytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::exbi(const Type& amount)
{
    return DataSize(amount*(1LL<<60LL));
}

// Predefined data size in kilobytes
/// \brief          Predefined data size in kilobytes.
/// \details        Constructs and returns a data size based on an amount of
///                 kilobytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of kilobytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::kilo(const Type& amount)
{
    return DataSize(amount*std::kilo::num);
}

// Predefined data size in megabytes
/// \brief          Predefined data size in megabytes.
/// \details        Constructs and returns a data size based on an amount of
///                 megabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of megabytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::mega(const Type& amount)
{
    return DataSize(amount*std::mega::num);
}

// Predefined data size in gigabytes
/// \brief          Predefined data size in gigabytes.
/// \details        Constructs and returns a data size based on an amount of
///                 gigabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of gigabytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::giga(const Type& amount)
{
    return DataSize(amount*std::giga::num);
}

// Predefined data size in terabytes
/// \brief          Predefined data size in terabytes.
/// \details        Constructs and returns a data size based on an amount of
///                 terabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of terabytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::tera(const Type& amount)
{
    return DataSize(amount*std::tera::num);
}

// Predefined data size in petabytes
/// \brief          Predefined data size in petabytes.
/// \details        Constructs and returns a data size based on an amount of
///                 petabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of petabytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::peta(const Type& amount)
{
    return DataSize(amount*std::peta::num);
}

// Predefined data size in exabytes
/// \brief          Predefined data size in exabytes.
/// \details        Constructs and returns a data size based on an amount of
///                 exabytes.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of exabytes.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::exa(const Type& amount)
{
    return DataSize(amount*std::exa::num);
}

// Predefined data size in bits
/// \brief          Predefined data size in bits.
/// \details        Constructs and returns a data size based on an amount of
///                 bits. If the number of provided bits is not divisible by
///                 8, an extra byte is added.
/// \tparam         Type (Data type.)
/// \param[in]      amount Amount of bits.
/// \return         Copy of the data size.
template <typename Type> 
constexpr DataSize DataSize::bit(const Type& amount)
{
    return DataSize(static_cast<long long int>(amount)/8LL+static_cast<bool>(static_cast<long long int>(amount)%8LL));
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Prints out the data size in the most convenient binary unit.
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side data size.
/// \return         Output stream.
std::ostream& operator<<(std::ostream& lhs, const DataSize& rhs)
{
    static const long long int threshold[7] = {1LL<<0LL, 1LL<<10LL, 1LL<<20LL, 1LL<<30LL, 1LL<<40LL, 1LL<<50LL, 1LL<<60LL};
    static const char suffix[7] = {'B', 'K', 'M', 'G', 'T', 'P', 'E'};
    const unsigned int category = (rhs._size < threshold[3]) ? ((rhs._size < threshold[1]) ? (0) : ((rhs._size < threshold[2]) ? (1) : (2))) : ((rhs._size < threshold[5]) ? ((rhs._size < threshold[4]) ? (3) : (4)) : ((rhs._size < threshold[6]) ? (5) : (6))); 
    const std::string unit = (rhs._size < 0) ? (std::string()) : (std::string(1, suffix[category])+((category > 0) ? (std::string("iB")) : (std::string())));
    return lhs<<((rhs._size < 0) ? (std::numeric_limits<double>::quiet_NaN()) : (static_cast<double>(rhs._size)/threshold[category]))<<unit;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of DataSize.
/// \return         0 if no error.
int DataSize::example()
{
    // Initialize
    std::cout<<"BEGIN = DataSize::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    DataSize d(8192);
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"DataSize() : "                                <<DataSize()<<std::endl;
    std::cout<<std::setw(width)<<"DataSize(42) : "                              <<DataSize(42)<<std::endl;
    std::cout<<std::setw(width)<<"DataSize(d) : "                               <<DataSize(d)<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Operators : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"d = 16384 : "                                 <<(d = 16384)<<std::endl;
    std::cout<<std::setw(width)<<"d = DataSize(32768) : "                       <<(d = DataSize(32768))<<std::endl;
    std::cout<<std::setw(width)<<"d() : "                                       <<d()<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"d.data() : "                                  <<d.data()<<std::endl;
    std::cout<<std::setw(width)<<"d.copy() : "                                  <<d.copy()<<std::endl;
    std::cout<<std::setw(width)<<"d.cast() : "                                  <<d.cast()<<std::endl;
    std::cout<<std::setw(width)<<"d.cast<double>() : "                          <<d.cast<double>()<<std::endl;

    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"d.valid() : "                                 <<d.valid()<<std::endl;
    std::cout<<std::setw(width)<<"d.empty() : "                                 <<d.empty()<<std::endl;
    std::cout<<std::setw(width)<<"d.bytes() : "                                 <<d.bytes()<<std::endl;
    std::cout<<std::setw(width)<<"d.kibibytes<double>() : "                     <<d.kibibytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.mebibytes<double>() : "                     <<d.mebibytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.gibibytes<double>() : "                     <<d.gibibytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.tebibytes<double>() : "                     <<d.tebibytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.pebibytes<double>() : "                     <<d.pebibytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.exbibytes<double>() : "                     <<d.exbibytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.kilobytes<double>() : "                     <<d.kilobytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.megabytes<double>() : "                     <<d.megabytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.gigabytes<double>() : "                     <<d.gigabytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.terabytes<double>() : "                     <<d.terabytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.petabytes<double>() : "                     <<d.petabytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.exabytes<double>() : "                      <<d.exabytes<double>()<<std::endl;
    std::cout<<std::setw(width)<<"d.bits() : "                                  <<d.bits()<<std::endl;

    // Setters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Setters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"d.valid(false) : "                            <<d.valid(false)<<std::endl;
    std::cout<<std::setw(width)<<"d.empty(true) : "                             <<d.empty(true)<<std::endl;
    std::cout<<std::setw(width)<<"d.bytes(42) : "                               <<d.bytes(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.kibibytes(42) : "                           <<d.kibibytes(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.mebibytes(42) : "                           <<d.mebibytes(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.gibibytes(42) : "                           <<d.gibibytes(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.tebibytes(42) : "                           <<d.tebibytes(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.pebibytes(42) : "                           <<d.pebibytes(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.exbibytes(1.42) : "                         <<d.exbibytes(1.42)<<std::endl;
    std::cout<<std::setw(width)<<"d.kilobytes(42.) : "                          <<d.kilobytes(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.megabytes(42.) : "                          <<d.megabytes(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.gigabytes(42.) : "                          <<d.gigabytes(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.terabytes(42.) : "                          <<d.terabytes(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.petabytes(42.) : "                          <<d.petabytes(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.exabytes(1.42) : "                          <<d.exabytes(1.42)<<std::endl;
    std::cout<<std::setw(width)<<"d.bits(42) : "                                <<d.bits(42)<<std::endl;
    
    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, d) : "                  <<d<<std::endl;

    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Predefined : "                                <<std::endl;
    std::cout<<std::setw(width)<<"d.byte() : "                                  <<d.byte()<<std::endl;
    std::cout<<std::setw(width)<<"d.kibi(42) : "                                <<d.kibi(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.mebi(42) : "                                <<d.mebi(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.gibi(42) : "                                <<d.gibi(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.tebi(42) : "                                <<d.tebi(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.pebi(42) : "                                <<d.pebi(42)<<std::endl;
    std::cout<<std::setw(width)<<"d.exbi(1.42) : "                              <<d.exbi(1.42)<<std::endl;
    std::cout<<std::setw(width)<<"d.kilo(42.) : "                               <<d.kilo(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.mega(42.) : "                               <<d.mega(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.giga(42.) : "                               <<d.giga(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.tera(42.) : "                               <<d.tera(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.peta(42.) : "                               <<d.peta(42.)<<std::endl;
    std::cout<<std::setw(width)<<"d.exa(1.42) : "                               <<d.exa(1.42)<<std::endl;
    std::cout<<std::setw(width)<<"d.bit(42) : "                                 <<d.bit(42)<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = DataSize::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // DATASIZE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
