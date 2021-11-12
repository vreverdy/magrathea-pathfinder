/* ******************************** FILELIST ******************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          FileList
// DESCRIPTION :    List of files based on a function or a vector
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           filelist.h
/// \brief          List of files based on a function or a vector
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef FILELIST_H_INCLUDED
#define FILELIST_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cctype>
#include <type_traits>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <initializer_list>
#include <array>
#include <vector>
#include <string>
#include <limits>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// List of files based on a function or a vector
/// \brief          List of files based on a function or a vector.
/// \details        Class to hold a list of paths to files. Internally, the list
///                 can be stored as a container of strings, a C-like format or 
///                 a function object that will generate the file names on the 
///                 fly. An offset allows to shift the real identificator of the
///                 files regarding the list index.
class FileList final
{
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        template <class Type = std::string, class = typename std::enable_if<std::is_convertible<Type, std::string>::value>::type> inline FileList(const std::initializer_list<Type>& source = std::initializer_list<Type>());
        template <class ContainerType, class PrefixType = std::string, class = typename std::enable_if<!(std::is_convertible<ContainerType, std::string>::value) && (std::is_convertible<typename std::remove_reference<decltype(*std::begin(std::declval<ContainerType>()))>::type, std::string>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> explicit inline FileList(const ContainerType& source, const PrefixType& prefix = PrefixType());
        template <typename IteratorType, class PrefixType = std::string, class = typename std::enable_if<(std::is_convertible<typename std::remove_reference<decltype(*std::declval<IteratorType>())>::type, std::string>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> explicit inline FileList(const IteratorType& begin, const IteratorType& end, const PrefixType& prefix = PrefixType());
        template <class FormatType, class PrefixType = std::string, class = typename std::enable_if<(std::is_convertible<FormatType, std::string>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> explicit inline FileList(const FormatType& cformat, const int width = 0, const unsigned int amount = 1, const unsigned int shift = 0, const PrefixType& prefix = PrefixType());
        template <class Function, class PrefixType = std::string, class = typename std::enable_if<(std::is_convertible<typename std::remove_reference<typename std::result_of<Function(const unsigned int)>::type>::type, std::string>::value) && !(std::is_same<typename std::decay<Function>::type, FileList>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> explicit inline FileList(Function&& f, const unsigned int amount = 1, const unsigned int shift = 0, const PrefixType& prefix = PrefixType());
    //@}
    
    // Operators
    /// \name           Operators
    //@{
    public:
        template <class ContainerType, class = typename std::enable_if<!(std::is_convertible<ContainerType, std::string>::value) && (std::is_convertible<typename std::remove_reference<decltype(*std::begin(std::declval<ContainerType>()))>::type, std::string>::value)>::type> inline FileList& operator=(const ContainerType& source);
        template <class FormatType, class = typename std::enable_if<std::is_convertible<FormatType, std::string>::value>::type, class = void> inline FileList& operator=(const FormatType& cformat);
        template <class Function, class = typename std::enable_if<(std::is_convertible<typename std::remove_reference<typename std::result_of<Function(const unsigned int)>::type>::type, std::string>::value) && !(std::is_same<typename std::decay<Function>::type, FileList>::value)>::type, class = void, class = void> inline FileList& operator=(Function&& f);
        inline bool operator==(const FileList& rhs) const;
        inline bool operator!=(const FileList& rhs) const;
    //@}
    
    // Assignment
    /// \name           Assignment
    //@{
    public:
        inline FileList& assign(const FileList& source);
        template <class ContainerType = std::initializer_list<std::string>, class PrefixType = std::string, class = typename std::enable_if<!(std::is_convertible<ContainerType, std::string>::value) && (std::is_convertible<typename std::remove_reference<decltype(*std::begin(std::declval<ContainerType>()))>::type, std::string>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> inline FileList& assign(const ContainerType& source = ContainerType(), const PrefixType& prefix = PrefixType());
        template <typename IteratorType, class PrefixType = std::string, class = typename std::enable_if<(std::is_convertible<typename std::remove_reference<decltype(*std::declval<IteratorType>())>::type, std::string>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> inline FileList& assign(const IteratorType& begin, const IteratorType& end, const PrefixType& prefix = PrefixType());
        template <class FormatType, class PrefixType = std::string, class = typename std::enable_if<(std::is_convertible<FormatType, std::string>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> inline FileList& assign(const FormatType& cformat, const int width = 0, const unsigned int amount = 1, const unsigned int shift = 0, const PrefixType& prefix = PrefixType());
        template <class Function, class PrefixType = std::string, class = typename std::enable_if<(std::is_convertible<typename std::remove_reference<typename std::result_of<Function(const unsigned int)>::type>::type, std::string>::value) && !(std::is_same<typename std::decay<Function>::type, FileList>::value) && (std::is_convertible<PrefixType, std::string>::value)>::type> inline FileList& assign(Function&& f, const unsigned int amount = 1, const unsigned int shift = 0, const PrefixType& prefix = PrefixType());
    //@}
    
    // Management
    /// \name           Management
    //@{
    public:
        inline bool empty() const;
        inline bool cleared() const;
        inline unsigned int capacity() const;
        inline FileList& resize(const unsigned int amount);
        inline FileList& reserve(const unsigned int amount);
        inline FileList& shrink();
        inline FileList& clear();
        inline FileList copy() const;
        template <class Type = FileList, class = typename std::enable_if<std::is_constructible<Type, FileList>::value>::type> inline Type cast() const;
    //@}
    
    // Getters
    /// \name           Getters
    //@{
    public:
        inline const std::vector<std::string>& container() const;
        inline const std::string& format() const;
        inline const std::function<std::string(const unsigned int)>& generator() const;
        inline const std::string& root() const;
        inline const int& length() const;
        inline const unsigned int& size() const;
        inline const unsigned int& offset() const;
    //@}

    // Files
    /// \name           Files
    //@{
    public:
        inline std::string operator[](const unsigned int i) const;
        inline std::string at(const unsigned int i) const;
        inline std::string front(const unsigned int i = 0) const;
        inline std::string back(const unsigned int i = 0) const;
        inline std::string operator()(const unsigned int i) const;
        inline std::string get(const unsigned int i) const;
    //@}
    
    // Manipulations
    /// \name           Manipulations
    //@{
    public:
        template <class Type = std::vector<std::string>, class = typename std::enable_if<((std::is_void<decltype(std::declval<Type>().resize(0))>::value) || (!std::is_void<decltype(std::declval<Type>().resize(0))>::value)) && (!std::is_const<decltype(std::declval<Type>()[0])>::value) && (std::is_reference<decltype(std::declval<Type>()[0])>::value) && (std::is_convertible<decltype(std::declval<Type>()[0]), std::string>::value)>::type> Type convert() const;
        std::string common() const;
        FileList formatify() const;
    //@}
    
    // Algorithms
    /// \name           Algorithms
    //@{
    public:
        template <class Function, class... Args, class = typename std::result_of<Function(std::string, Args...)>::type> FileList& apply(Function&& f, Args&&... args);
        unsigned int count() const;
        template <class Function, class = typename std::enable_if<std::is_convertible<typename std::result_of<Function(std::string)>::type, bool>::value>::type> unsigned int count(Function&& f) const;
        template <class Function, typename FirstType = std::true_type, typename AmountType = std::true_type, class = typename std::enable_if<(std::is_convertible<typename std::result_of<Function(std::string)>::type, bool>::value) && ((std::is_same<AmountType, std::true_type>::value) || (std::is_convertible<FirstType, unsigned int>::value)) && ((std::is_same<AmountType, std::true_type>::value) || (std::is_convertible<AmountType, int>::value))>::type> unsigned int find(Function&& f, const FirstType first = FirstType(), const AmountType amount = AmountType()) const;
        template <class Function = std::less<std::string>, class = typename std::enable_if<std::is_convertible<typename std::result_of<Function(std::string, std::string)>::type, bool>::value>::type> FileList& sort(Function&& f = Function());
        template <class Function = std::less<std::string>, class = typename std::enable_if<std::is_convertible<typename std::result_of<Function(std::string, std::string)>::type, bool>::value>::type> bool sorted(Function&& f = Function()) const;
        template <class Function = std::equal_to<std::string>, class = typename std::enable_if<std::is_convertible<typename std::result_of<Function(std::string, std::string)>::type, bool>::value>::type> FileList& unique(Function&& f = Function());
        template <class Function = std::equal_to<std::string>, class = typename std::enable_if<std::is_convertible<typename std::result_of<Function(std::string, std::string)>::type, bool>::value>::type> bool unicity(Function&& f = Function()) const;
    //@}

    // Stream
    /// \name           Stream
    //@{
    public:
        friend std::ostream& operator<<(std::ostream& lhs, const FileList& rhs);
    //@}
    
    // Helpers
    /// \name           Helpers
    //@{
    public:
        template <typename IndexType = std::true_type, class = typename std::enable_if<(std::is_same<IndexType, std::true_type>::value) || (std::is_convertible<IndexType, unsigned int>::value)>::type> static inline std::string generate(const IndexType i = IndexType());
        template <class ContainerType, class = typename std::enable_if<!(std::is_convertible<ContainerType, std::string>::value) && (std::is_convertible<typename std::remove_reference<decltype(*std::begin(std::declval<ContainerType>()))>::type, std::string>::value)>::type> static inline std::string generate(const unsigned int i, const ContainerType& source);
        template <unsigned int BufferSize = std::numeric_limits<unsigned char>::max()*std::numeric_limits<unsigned char>::digits, typename WidthType = std::true_type, class = typename std::enable_if<(std::is_same<WidthType, std::true_type>::value) || (std::is_convertible<WidthType, int>::value)>::type> static inline std::string generate(const unsigned int i, const std::string& cformat, const WidthType width = WidthType());
        template <class Function, class = typename std::enable_if<(std::is_convertible<typename std::remove_reference<typename std::result_of<Function(const unsigned int)>::type>::type, std::string>::value) && !(std::is_same<typename std::decay<Function>::type, FileList>::value)>::type> static inline std::string generate(const unsigned int i, Function&& f);
        template <typename Dir = unsigned int, typename Min = unsigned int, typename Num = unsigned int, typename Max = unsigned int, class = typename std::enable_if<(std::is_convertible<Dir, unsigned int>::value) && (std::is_convertible<Min, unsigned int>::value) && (std::is_convertible<Num, unsigned int>::value) && (std::is_convertible<Max, unsigned int>::value)>::type> static unsigned int numberify(const std::string& str, unsigned int pos = std::numeric_limits<unsigned int>::max(), Dir&& dir = Dir(), Min&& min = Min(), Num&& num = Num(), Max&& max = Max());
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
    public:
        std::vector<std::string> _container;                                    ///< Container of file names.
        std::string _format;                                                    ///< Format of file names.
        std::function<std::string(const unsigned int)> _generator;              ///< Generator of file names.
        std::string _root;                                                      ///< Common root prefix.
        int _length;                                                            ///< Optional width of the format field.
        unsigned int _size;                                                     ///< Total file list size.
        unsigned int _offset;                                                   ///< Offset between file identificator and list index.
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Implicit initializer list constructor
/// \brief          Impliticit initializer list constructor.
/// \details        Implicitely constructs the file list from an initializer 
///                 list of paths.
/// \tparam         Type (Data type.)
/// \param[in]      source Source of the copy.
template <class Type, class>
inline FileList::FileList(const std::initializer_list<Type>& source)
: _container(std::begin(source), std::end(source))
, _format()
, _generator()
, _root()
, _length(0)
, _size(source.size())
, _offset(0)
{
    ;
}

// Explicit container constructor
/// \brief          Explicit container constructor.
/// \details        Explicitely constructs the file list from a container of
///                 paths.
/// \tparam         ContainerType (Container type.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      source Source of the copy.
/// \param[in]      prefix Common root prefix.
template <class ContainerType, class PrefixType, class>
inline FileList::FileList(const ContainerType& source, const PrefixType& prefix)
: _container(std::begin(source), std::end(source))
, _format()
, _generator()
, _root(prefix)
, _length(0)
, _size(std::end(source)-std::begin(source))
, _offset(0)
{
    ;
}

// Explicit iterator constructor
/// \brief          Explicit iterator constructor.
/// \details        Explicitely constructs the file list from a range of 
///                 iterators.
/// \tparam         IteratorType (Iterator or pointer type.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      begin Iterator to the beginning of the range.
/// \param[in]      end Iterator to the end of the range.
/// \param[in]      prefix Common root prefix.
template <typename IteratorType, class PrefixType, class>
inline FileList::FileList(const IteratorType& begin, const IteratorType& end, const PrefixType& prefix)
: _container(begin, end)
, _format()
, _generator()
, _root(prefix)
, _length(0)
, _size(end-begin)
, _offset(0)
{
    ;
}

// Explicit format constructor
/// \brief          Explicit format constructor.
/// \details        Explicitely constructs the file list from a C-like format
///                 string and its width if required.
/// \tparam         FormatType (Type of the format string.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      cformat Format of file names.
/// \param[in]      width Optional width of the format field.
/// \param[in]      amount Total file list size.
/// \param[in]      shift Offset between file identificator and list index.
/// \param[in]      prefix Common root prefix.
template <class FormatType, class PrefixType, class> 
inline FileList::FileList(const FormatType& cformat, const int width, const unsigned int amount, const unsigned int shift, const PrefixType& prefix)
: _container()
, _format(cformat)
, _generator()
, _root(prefix)
, _length(width)
, _size(amount)
, _offset(shift)
{
    ;
}

// Explicit generator constructor
/// \brief          Explicit generator constructor.
/// \details        Explicitely constructs the file list from a generator 
///                 function.
/// \tparam         Function (Function type : <tt>std::string(const unsigned 
///                 int)</tt>.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      f Function object <tt>std::string(const unsigned int)</tt>.
/// \param[in]      amount Total file list size.
/// \param[in]      shift Offset between file identificator and list index.
/// \param[in]      prefix Common root prefix.
template <class Function, class PrefixType, class>
inline FileList::FileList(Function&& f, const unsigned int amount, const unsigned int shift, const PrefixType& prefix)
: _container()
, _format()
, _generator(f)
, _root(prefix)
, _length(0)
, _size(amount)
, _offset(shift)
{
    ;
}
// -------------------------------------------------------------------------- //



// -------------------------------- OPERATORS ------------------------------- //
// Container assignment operator
/// \brief          Container assignment operator.
/// \details        Assigns the file list contents using a container of paths
///                 and without any root prefix.
/// \tparam         ContainerType (Container type.)
/// \param[in]      source Source of the copy.
/// \return         Self reference.
template <class ContainerType, class> 
inline FileList& FileList::operator=(const ContainerType& source)
{
    _container.assign(std::begin(source), std::end(source));
    _format.clear();
    _generator = std::function<std::string(const unsigned int)>();
    _root.clear();
    _length = 0;
    _size = std::end(source)-std::begin(source);
    _offset = 0;
    return *this;
}

// Format assignment operator
/// \brief          Format assignment operator.
/// \details        Assigns the file list contents using a C-like format string
///                 without any root prefix and with a size of one element.
/// \tparam         FormatType (Type of the format string.)
/// \param[in]      cformat Format of file names.
/// \return         Self reference.
template <class FormatType, class, class> 
inline FileList& FileList::operator=(const FormatType& cformat)
{
    _container.clear();
    _format = cformat;
    _generator = std::function<std::string(const unsigned int)>();
    _root.clear();
    _length = 0;
    _size = 1;
    _offset = 0;
    return *this;
}

// Generator assignment operator
/// \brief          Generator assignment operator.
/// \details        Assigns the file list contents usign a generator function
///                 without any root prefix and with a size of one element.
/// \tparam         Function (Function type : <tt>std::string(const unsigned 
///                 int)</tt>.)
/// \param[in]      f Function object <tt>std::string(const unsigned int)</tt>.
/// \return         Self reference.
template <class Function, class, class, class> 
inline FileList& FileList::operator=(Function&& f)
{
    _container.clear();
    _format.clear();
    _generator = f;
    _root.clear();
    _length = 0;
    _size = 1;
    _offset = 0;
    return *this;
}

// Equal to
/// \brief          Equal to.
/// \details        Compares for equality and returns true if all elements 
///                 returned by the subscript operator are equal: the 
///                 internal storage is not compared and consequently a
///                 container-based file list can be equal to a format-based
///                 file list.
/// \param[in]      rhs Right-hand side. 
/// \return         True if equal, false if not equal.
inline bool FileList::operator==(const FileList& rhs) const
{
    bool result = (_size == rhs._size);
    for (unsigned int i = 0; (i < _size) && (result); ++i) {
        result = ((*this)[i] == rhs[i]);
    }
    return result;
}

// Not equal to
/// \brief          Not equal to.
/// \details        Compares for difference and returns true if two elements 
///                 returned by the subscript operator differs: the 
///                 internal storage is not compared.
/// \param[in]      rhs Right-hand side. 
/// \return         True if not equal, false if equal.
inline bool FileList::operator!=(const FileList& rhs) const
{
    bool result = (_size == rhs._size);
    for (unsigned int i = 0; (i < _size) && (result); ++i) {
        result = ((*this)[i] == rhs[i]);
    }
    return !result;
}
// -------------------------------------------------------------------------- //



// ------------------------------- ASSIGNMENT ------------------------------- //
// Copy assignment
/// \brief          Copy assignment.
/// \details        Assigns the file list contents using another file list.
/// \param[in]      source Source of the copy.
/// \return         Self reference.
inline FileList& FileList::assign(const FileList& source)
{
    _container = source._container;
    _format = source._format;
    _generator = source._generator;
    _root = source._root;
    _length = source._length;
    _size = source._size;
    _offset = source._offset;
    return *this;
}

// Container assignment
/// \brief          Container assignment.
/// \details        Assigns the file list contents using a container of paths.
/// \tparam         ContainerType (Container type.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      source Source of the copy.
/// \param[in]      prefix Common root prefix.
/// \return         Self reference.
template <class ContainerType, class PrefixType, class> 
inline FileList& FileList::assign(const ContainerType& source, const PrefixType& prefix)
{
    _container.assign(std::begin(source), std::end(source));
    _format.clear();
    _generator = std::function<std::string(const unsigned int)>();
    _root = prefix;
    _length = 0;
    _size = std::end(source)-std::begin(source);
    _offset = 0;
    return *this;
}

// Iterator assignment
/// \brief          Iterator assignment.
/// \details        Assigns the file list contents using a range of iterators.
/// \tparam         IteratorType (Iterator or pointer type.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      begin Iterator to the beginning of the range.
/// \param[in]      end Iterator to the end of the range.
/// \param[in]      prefix Common root prefix.
/// \return         Self reference.
template <typename IteratorType, class PrefixType, class> 
inline FileList& FileList::assign(const IteratorType& begin, const IteratorType& end, const PrefixType& prefix)
{
    _container.assign(begin, end);
    _format.clear();
    _generator = std::function<std::string(const unsigned int)>();
    _root = prefix;
    _length = 0;
    _size = end-begin;
    _offset = 0;
    return *this;
}

// Format assignment
/// \brief          Format assignment.
/// \details        Assigns the file list contents using a C-like format string
///                 and its width if required.
/// \tparam         FormatType (Type of the format string.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      cformat Format of file names.
/// \param[in]      width Optional width of the format field.
/// \param[in]      amount Total file list size.
/// \param[in]      shift Offset between file identificator and list index.
/// \param[in]      prefix Common root prefix.
/// \return         Self reference.
template <class FormatType, class PrefixType, class> 
inline FileList& FileList::assign(const FormatType& cformat, const int width, const unsigned int amount, const unsigned int shift, const PrefixType& prefix)
{
    _container.clear();
    _format = cformat;
    _generator = std::function<std::string(const unsigned int)>();
    _root = prefix;
    _length = width;
    _size = amount;
    _offset = shift;
    return *this;
}

// Generator assignment
/// \brief          Generator assignment.
/// \details        Assigns the file list contents usign a generator function.
/// \tparam         Function (Function type : <tt>std::string(const unsigned 
///                 int)</tt>.)
/// \tparam         PrefixType (Prefix type.)
/// \param[in]      f Function object <tt>std::string(const unsigned int)</tt>.
/// \param[in]      amount Total file list size.
/// \param[in]      shift Offset between file identificator and list index.
/// \param[in]      prefix Common root prefix.
/// \return         Self reference.
template <class Function, class PrefixType, class> 
inline FileList& FileList::assign(Function&& f, const unsigned int amount, const unsigned int shift, const PrefixType& prefix)
{
    _container.clear();
    _format.clear();
    _generator = f;
    _root = prefix;
    _length = 0;
    _size = amount;
    _offset = shift;
    return *this;
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Empty
/// \brief          Empty.
/// \details        Returns true if the file list size is equal to zero.
/// \return         True if the file list is empty, false otherwise.
inline bool FileList::empty() const
{
    return (_size == 0);
}

// Cleared
/// \brief          Cleared.
/// \details        Returns true if the file list is equal to a cleared one.
///                 It corresponds to an empty size, an empty format and an 
///                 empty function.
/// \return         True if the file list is cleared, false otherwise.
inline bool FileList::cleared() const
{
    return ((_container.empty()) && (_format.empty()) && (!bool(_generator)) && (_root.empty()) && (_length == 0) && (_size == 0) && (_offset == 0));
}
    
// Capacity
/// \brief          Capacity.
/// \details        Returns the current capacity of the underlying container if
///                 it is used to store the file names, or returns the size if
///                 a format or a generator function is used.
/// \return         Current capacity.
inline unsigned int FileList::capacity() const
{
    return ((_generator) || (!_format.empty())) ? (_size) : (_container.capacity());
}

// Resize
/// \brief          Resize.
/// \details        Resizes the current list and returns it. In case of a
///                 container based list, the new contents is initialized as
///                 empty. Note that setting a null size, empties the list but
///                 does not clear it.
/// \param[in]      amount New size.
/// \return         Self reference.
inline FileList& FileList::resize(const unsigned int amount)
{
    if ((_generator) || (!_format.empty())) {
        _size = amount;
    } else {
        _container.resize(amount);
        _size = _container.size();
    }
    return *this;
}

// Reserve space for the container
/// \brief          Reserve space for the container.
/// \details        Increases the capacity of the container if the specified 
///                 size is greater than the current one.
/// \param[in]      amount New capacity.
/// \return         Self reference.
inline FileList& FileList::reserve(const unsigned int amount)
{
    if (amount > _container.size()) {
        _container.reserve(amount);
    }
    return *this;
}

// Shrink the container
/// \brief          Shrink the container.
/// \details        Shrinks the container capacity to fit its size.
/// \return         Self reference.
inline FileList& FileList::shrink()
{
    _container.shrink_to_fit();
    return *this;
}

// Clear
/// \brief          Clear.
/// \details        Clears the whole contents of the file list.
/// \return         Self reference.
inline FileList& FileList::clear()
{
    _container.clear();
    _format.clear();
    _generator = std::function<std::string(const unsigned int)>();
    _root.clear();
    _length = 0;
    _size = 0;
    _offset = 0;
    return *this;
}

// Copy
/// \brief          Copy.
/// \details        Returns a copy of the file list.
/// \return         Copy.
inline FileList FileList::copy() const
{
    return *this;
}

// Cast
/// \brief          Cast.
/// \details        Returns a copy of the file list casted to the provided 
///                 type.
/// \tparam         Type Type.
/// \return         Casted copy.
template <class Type, class>
inline Type FileList::cast() const
{
    return Type(*this);
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// Container getter
/// \brief          Container getter.
/// \details        Gets the container of file names.
/// \return         Immutable reference to the container.
inline const std::vector<std::string>& FileList::container() const
{
    return _container;
}

// Format getter
/// \brief          Format getter.
/// \details        Gets the format of file names.
/// \return         Immutable reference to the format.
inline const std::string& FileList::format() const
{
    return _format;
}

// Generator getter
/// \brief          Generator getter.
/// \details        Gets the generator of file names.
/// \return         Immutable reference to the generator.
inline const std::function<std::string(const unsigned int)>& FileList::generator() const
{
    return _generator;
}

// Root getter
/// \brief          Root getter.
/// \details        Gets the common root prefix.
/// \return         Immutable reference to the root.
inline const std::string& FileList::root() const
{
    return _root;
}

// Length getter
/// \brief          Length getter.
/// \details        Gets the optional width of the format field.
/// \return         Immutable reference to the length.
inline const int& FileList::length() const
{
    return _length;
}

// Size getter
/// \brief          Size getter.
/// \details        Gets the total file list size.
/// \return         Immutable reference to the size.
inline const unsigned int& FileList::size() const
{
    return _size;
}

// Offset getter
/// \brief          Offset getter.
/// \details        Gets the offset between file identificator and list index.
/// \return         Immutable reference to the offset.
inline const unsigned int& FileList::offset() const
{
    return _offset;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- FILES --------------------------------- //
// Subscript operator
/// \brief          Subscript operator.
/// \details        Returns the file name at the specified position in the file 
///                 list or returns an empty string if it does not exist.
/// \param[in]      i Index.
/// \return         Copy of the specified element.
inline std::string FileList::operator[](const unsigned int i) const
{
    return (i < _size) ? ((_generator) ? (_root+_generator(_offset+i)) : ((!_format.empty()) ? ((_length != 0) ? (_root+generate(_offset+i, _format, _length)) : (_root+generate(_offset+i, _format))) : (_root+_container[i]))) : (std::string());
}

// Range check access
/// \brief          Range check access.
/// \details        Returns the file name at the specified position in the file 
///                 list or throws an error if the index cannot be reached.
/// \param[in]      i Index.
/// \return         Copy of the specified element.
/// \exception      std::out_of_range Out of range.
inline std::string FileList::at(const unsigned int i) const
{
    if (i >= _size) {
        throw std::out_of_range("ERROR = FileList::at() : out of range");
    }
    return (_generator) ? (_root+_generator(_offset+i)) : ((!_format.empty()) ? ((_length != 0) ? (_root+generate(_offset+i, _format, _length)) : (_root+generate(_offset+i, _format))) : (_root+_container[i]));
}

// Access from the beginning
/// \brief          Access from the beginning.
/// \details        Returns the file name at the specified position in the file 
///                 list starting from the beginning and without any range 
///                 check.
/// \param[in]      i Index.
/// \return         Copy of the specified element.
inline std::string FileList::front(const unsigned int i) const
{
    return (i < _size) ? ((_generator) ? (_root+_generator(_offset+i)) : ((!_format.empty()) ? ((_length != 0) ? (_root+generate(_offset+i, _format, _length)) : (_root+generate(_offset+i, _format))) : (_root+_container[i]))) : (std::string());
}

// Access from the end
/// \brief          Access from the end.
/// \details        Returns the file name at the specified position in the file 
///                 list starting from the end and without any range check.
/// \param[in]      i Index.
/// \return         Copy of the specified element.
inline std::string FileList::back(const unsigned int i) const
{
    return (i < _size) ? ((_generator) ? (_root+_generator(static_cast<int>(_offset+_size)-static_cast<int>(i+1))) : ((!_format.empty()) ? ((_length != 0) ? (_root+generate(static_cast<int>(_offset+_size)-static_cast<int>(i+1), _format, _length)) : (_root+generate(static_cast<int>(_offset+_size)-static_cast<int>(i+1), _format))) : (_root+_container[static_cast<int>(_size)-static_cast<int>(i+1)]))) : (std::string());
}


// Generator operator
/// \brief          Generator operator.
/// \details        Generates the file name using the specified index and
///                 ignoring the root, the size limit (except if the list is
///                 based on a container) and the offset.
/// \param[in]      i Index.
/// \return         Copy of the generated element.
inline std::string FileList::operator()(const unsigned int i) const
{
    return (_generator) ? (_generator(i)) : ((!_format.empty()) ? ((_length != 0) ? (generate(i, _format, _length)) : (generate(i, _format))) : ((i < _size) ? (_container[i]) : (std::string())));
}

// Generator getter
/// \brief          Generator getter.
/// \details        Generates the file name using the specified index and
///                 ignoring the root, the size limit (except if the list is
///                 based on a container) and the offset.
/// \param[in]      i Index.
/// \return         Copy of the generated element.
inline std::string FileList::get(const unsigned int i) const
{
    return (_generator) ? (_generator(i)) : ((!_format.empty()) ? ((_length != 0) ? (generate(i, _format, _length)) : (generate(i, _format))) : ((i < _size) ? (_container[i]) : (std::string())));
}
// -------------------------------------------------------------------------- //



// ------------------------------ MANIPULATIONS ----------------------------- //
// Convert
/// \brief          Convert.
/// \details        Converts data to the specified type of container.
/// \tparam         Type Type.
/// \return         Container of file names.
template <class Type, class> 
Type FileList::convert() const
{
    const unsigned int bsize = (sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length)*(!_format.empty());
    std::vector<char> buffer(bsize);
    Type result(_size);
    if (_generator) {
        for (unsigned int i = 0; i < _size; ++i) {
            result[i] = _root+_generator(_offset+i);
        }
    } else if (!_format.empty()) {
        if (_length != 0) {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _length, _offset+i);
                result[i] = _root+buffer.data();
            }
        } else {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _offset+i);
                result[i] = _root+buffer.data();
            }
        }
    } else {
        for (unsigned int i = 0; i < _size; ++i) {
            result[i] = _root+_container[i];
        }
    }
    return result;
}

// Common file part
/// \brief          Common file part.
/// \details        Detects the longest common character sequence within the 
///                 whole list, starting from the beginning.
/// \return         Longest common part.
std::string FileList::common() const
{
    const unsigned int bsize = (sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length)*(!_format.empty());
    std::vector<char> buffer(bsize);
    std::string result;
    std::string temporary;
    unsigned int n = 0;
    if (_size > 0) {
        if (_generator) {
            result = _generator(_offset);
            n = result.size();
            for (unsigned int i = 1; (i < _size) && (n > 0); ++i) {
                temporary = _generator(_offset+i);
                if (temporary.size() < n) {
                    n = temporary.size();
                }
                for (unsigned int j = 0; j < n; ++j) {
                    if (temporary[j] != result[j]) {
                        n = j;
                    }
                }
            }
        } else if (!_format.empty()) {
            (void)((_length != 0) ? (std::snprintf(buffer.data(), bsize, _format.data(), _length, 0)) : (std::snprintf(buffer.data(), bsize, _format.data(), 0)));
            result = buffer.data();
            n = (_format.size() < result.size()) ? (_format.size()) : (result.size());
            for (unsigned int i = 0; i < n; ++i) {
                if (_format[i] != result[i]) {
                    n = i;
                }
            }
        } else {
            result = _container[0];
            n = result.size();
            for (unsigned int i = 1; (i < _size) && (n > 0); ++i) {
                if (_container[i].size() < n) {
                    n = _container[i].size();
                }
                for (unsigned int j = 0; j < n; ++j) {
                    if (_container[i][j] != result[j]) {
                        n = j;
                    }
                }
            }
        }
    }
    return _root+result.erase(n);
}

// Detect common format
/// \brief          Detect common format.
/// \details        Computes whether a common format can be deduced from the 
///                 file list and returns a format based file list on success
///                 or an empty one on failure. Note that in case of a function
///                 or a container based list with a single element, the last
///                 number of the string is converted to a format.
/// \return         Formatted file list.
FileList FileList::formatify() const
{
    unsigned int bsize = 0;
    char* buffer = nullptr;
    std::string cformat;
    int width = 0;
    unsigned int amount = 0;
    unsigned int shift = 0;
    std::string prefix;
    std::string filename;
    std::string temporary;
    unsigned int n = 0;
    unsigned int ndir = 0;
    unsigned int nmin = 0;
    unsigned int nnum = 0;
    unsigned int nmax = 0;
    if (_size > 0) {
        if (_generator) {
            filename = _generator(_offset);
            n = filename.size()-((_size == 1) && (!filename.empty()));
            for (unsigned int i = 1; (i < _size) && (n > 0); ++i) {
                temporary = _generator(_offset+i);
                if (temporary.size() < n) {
                    n = temporary.size();
                }
                for (unsigned int j = 0; j < n; ++j) {
                    if (temporary[j] != filename[j]) {
                        n = j;
                    }
                }
            }
            if (n < filename.size()) {
                if ((std::isdigit(filename[n])) || (_size == 1)) {
                    shift = numberify(filename, n, ndir, nmin, nnum, nmax);
                    if (nnum < filename.size()) {
                        cformat = filename.substr(0, nmin)+((nmin < nnum) ? (std::string("%0*d")) : (std::string("%d")))+filename.substr(nmax);
                        width = nnum-nmin+(nmin < nnum);
                        bsize = sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+cformat.size()+width;
                        buffer = new char[bsize];
                        if (width != 0) {
                            for (unsigned int i = 0; (i < _size) && (n < filename.size()); ++i) {
                                std::snprintf(buffer, bsize, cformat.data(), width, shift+i);
                                n = (_generator(_offset+i) == buffer) ? (n) : (filename.size());
                               }
                        } else {
                            for (unsigned int i = 0; (i < _size) && (n < filename.size()); ++i) {
                                std::snprintf(buffer, bsize, cformat.data(), shift+i);
                                n = (_generator(_offset+i) == buffer) ? (n) : (filename.size());
                            }
                        }
                        if (n < filename.size()) {
                            prefix.append(cformat.begin(), cformat.begin()+ndir);
                            cformat.erase(cformat.begin(), cformat.begin()+ndir);
                            amount = _size;
                        }
                    }
                }
            }
        } else if (!_format.empty()) {
            bsize = sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length;
            buffer = new char[bsize];
            (void)((_length != 0) ? (std::snprintf(buffer, bsize, _format.data(), _length, _offset)) : (std::snprintf(buffer, bsize, _format.data(), _offset)));
            filename = buffer;
            n = (_format.size() < filename.size()) ? (_format.size()) : (filename.size());
            for (unsigned int i = 0; i < n; ++i) {
                if (_format[i] != filename[i]) {
                    n = i;
                }
            }
            n -= ((_size == 1) && (n == filename.size()));
            if (n < filename.size()) {
                if ((std::isdigit(filename[n]))  || (_size == 1)) {
                    shift = numberify(filename, n, ndir, nmin, nnum, nmax);
                    if (nnum < filename.size()) {
                        cformat = filename.substr(0, nmin)+((nmin < nnum) ? (std::string("%0*d")) : (std::string("%d")))+filename.substr(nmax);
                        width = nnum-nmin+(nmin < nnum);
                        (void)((width != 0) ? (std::snprintf(buffer, bsize, cformat.data(), width, shift)) : (std::snprintf(buffer, bsize, cformat.data(), shift)));
                        if (filename == buffer) {
                            (void)((_length != 0) ? (std::snprintf(buffer, bsize, _format.data(), _length, _offset+_size-(_size > 0))) : (std::snprintf(buffer, bsize, _format.data(), _offset+_size-(_size > 0))));
                            temporary = buffer;
                            (void)((width != 0) ? (std::snprintf(buffer, bsize, cformat.data(), width, shift+_size-(_size > 0))) : (std::snprintf(buffer, bsize, cformat.data(), shift+_size-(_size > 0))));
                            if (temporary == buffer) {
                                prefix.append(cformat.begin(), cformat.begin()+ndir);
                                cformat.erase(cformat.begin(), cformat.begin()+ndir);
                                amount = _size;
                            }
                        }
                    }
                }
            }
        } else {
            filename = _container[0];
            n = filename.size()-((_size == 1) && (!filename.empty()));
            for (unsigned int i = 1; (i < _size) && (n > 0); ++i) {
                if (_container[i].size() < n) {
                    n = _container[i].size();
                }
                for (unsigned int j = 0; j < n; ++j) {
                    if (_container[i][j] != filename[j]) {
                        n = j;
                    }
                }
            }
            if (n < filename.size()) {
                if ((std::isdigit(filename[n])) || (_size == 1)) {
                    shift = numberify(filename, n, ndir, nmin, nnum, nmax);
                    if (nnum < filename.size()) {
                        cformat = filename.substr(0, nmin)+((nmin < nnum) ? (std::string("%0*d")) : (std::string("%d")))+filename.substr(nmax);
                        width = nnum-nmin+(nmin < nnum);
                        bsize = sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+cformat.size()+width;
                        buffer = new char[bsize];
                        if (width != 0) {
                            for (unsigned int i = 0; (i < _size) && (n < filename.size()); ++i) {
                                std::snprintf(buffer, bsize, cformat.data(), width, shift+i);
                                n = (_container[i] == buffer) ? (n) : (filename.size());
                               }
                        } else {
                            for (unsigned int i = 0; (i < _size) && (n < filename.size()); ++i) {
                                std::snprintf(buffer, bsize, cformat.data(), shift+i);
                                n = (_container[i] == buffer) ? (n) : (filename.size());
                            }
                        }
                        if (n < filename.size()) {
                            prefix.append(cformat.begin(), cformat.begin()+ndir);
                            cformat.erase(cformat.begin(), cformat.begin()+ndir);
                            amount = _size;
                        }
                    }
                }
            }
        } 
    }
    if (buffer != nullptr) {
        delete[] buffer;
        buffer = nullptr;
    }
    return (amount > 0) ? (FileList(cformat, width, amount, shift, prefix)) : (FileList());
}
// -------------------------------------------------------------------------- //



// ------------------------------- ALGORITHMS ------------------------------- //
// Application of a function on each file
/// \brief          Application of a function on each file.
/// \details        Applies the provided algorithm on each file of the list,
///                 passing the extra arguments as extra function parameters if
///                 specified.
/// \tparam         Function (Function type : <tt>Type(std::string, Args...)
///                 </tt>.)
/// \tparam         Args (Extra types.)
/// \param[in]      f Function object <tt>Type(std::string, Args...)</tt>.
/// \param[in]      args Extra arguments of the function.
/// \return         Self reference.
template <class Function, class... Args, class> 
FileList& FileList::apply(Function&& f, Args&&... args)
{
    const unsigned int bsize = (sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length)*(!_format.empty());
    std::vector<char> buffer(bsize);
    if (_generator) {
        for (unsigned int i = 0; i < _size; ++i) {
            f(_root+_generator(_offset+i), args...);
        }
    } else if (!_format.empty()) {
        if (_length != 0) {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _length, _offset+i);
                f(_root+buffer.data(), args...);
            }
        } else {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _offset+i);
                f(_root+buffer.data(), args...);
            }
        }
    } else {
        for (unsigned int i = 0; i < _size; ++i) {
            f(_root+_container[i], args...);
        }
    }
    return *this;
}
        
// Count existing files
/// \brief          Count existing files.
/// \details        Counts the number of existing files in the list.
/// \return         Number of existing files.
unsigned int FileList::count() const
{
    const unsigned int bsize = (sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length)*(!_format.empty());
    std::vector<char> buffer(bsize);
    unsigned int result = 0;
    FILE* file = nullptr;
    if (_generator) {
        for (unsigned int i = 0; i < _size; ++i) {
            file = std::fopen((_root+_generator(_offset+i)).data(), "rb");
            result += ((file != nullptr) ? (!std::fclose(file)) : (false));
        }
    } else if (!_format.empty()) {
        if (_length != 0) {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _length, _offset+i);
                file = std::fopen((_root+buffer.data()).data(), "rb");
                result += ((file != nullptr) ? (!std::fclose(file)) : (false));
            }
        } else {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _offset+i);
                file = std::fopen((_root+buffer.data()).data(), "rb");
                result += ((file != nullptr) ? (!std::fclose(file)) : (false));
            }
        }
    } else {
        for (unsigned int i = 0; i < _size; ++i) {
            file = std::fopen((_root+_container[i]).data(), "rb");
            result += ((file != nullptr) ? (!std::fclose(file)) : (false));
        }
    }
    return result;
}

// Count files satisfying a condition
/// \brief          Count files satisfying a condition.
/// \details        Counts the number of files in the list that satisfy the 
///                 provided condition.
/// \tparam         Function (Function type : <tt>bool(std::string)</tt>.)
/// \param[in]      f Function object <tt>bool(std::string)</tt>.
/// \return         Number of files satisfying the condition.
template <class Function, class> 
unsigned int FileList::count(Function&& f) const
{
    const unsigned int bsize = (sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length)*(!_format.empty());
    std::vector<char> buffer(bsize);
    unsigned int result = 0;
    if (_generator) {
        for (unsigned int i = 0; i < _size; ++i) {
            result += f(_root+_generator(_offset+i));
        }
    } else if (!_format.empty()) {
        if (_length != 0) {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _length, _offset+i);
                result += f(_root+buffer.data());
            }
        } else {
            for (unsigned int i = 0; i < _size; ++i) {
                std::snprintf(buffer.data(), bsize, _format.data(), _offset+i);
                result += f(_root+buffer.data());
            }
        }
    } else {
        for (unsigned int i = 0; i < _size; ++i) {
            result += f(_root+_container[i]);
        }
    }
    return result;
}

// Find first file satisfying a condition
/// \brief          Find first file satisfying a condition.
/// \details        Computes the index of the first file satisfying the 
///                 specified predicate. The first and amount arguments allows
///                 to specify the first file to test and the amount of file to 
///                 test. As an example, an amount of zero file do not test any
///                 file, and an amount of -2 file tests the first file and the
///                 previous one. If no file satisfy the predicate, the function
///                 returns the size of the list. Boundary overflows and 
///                 unsigned substractions are correctly handled.
/// \tparam         Function (Function type : <tt>bool(std::string)</tt>.)
/// \tparam         FirstType (Integral type of the first index.)
/// \tparam         AmountType (Integral type of the first index.)
/// \param[in]      f Function object <tt>bool(std::string)</tt>.
/// \param[in]      first First file index to test.
/// \param[in]      amount Amount of file to test in the corresponding 
///                 direction.
/// \return         Index of the first file satisfying the predicate.
template <class Function, typename FirstType, typename AmountType, class>
unsigned int FileList::find(Function&& f, const FirstType first, const AmountType amount) const
{
    const unsigned int bsize = (sizeof(unsigned int)*std::numeric_limits<unsigned char>::digits+_format.size()+_length)*(!_format.empty());
    std::vector<char> buffer(bsize);
    const int idiff = ((std::is_same<AmountType, std::true_type>::value) || (std::is_unsigned<AmountType>::value) || (amount > 0)) ? (1) : (-1);
    const int ibegin = ((std::is_same<FirstType, std::true_type>::value) || (static_cast<int>(first) < 0)) ? (-static_cast<int>((!std::is_same<AmountType, std::true_type>::value) && (static_cast<int>(amount) < 0))) : ((static_cast<int>(first) >= static_cast<int>(_size)) ? (static_cast<int>(_size)-static_cast<int>((!std::is_same<AmountType, std::true_type>::value) && (static_cast<int>(amount) < 0))) : (static_cast<int>(first)));
    const int iend = ((std::is_same<FirstType, std::true_type>::value) || (std::is_same<AmountType, std::true_type>::value) || (static_cast<int>(first)+static_cast<int>(amount) >= static_cast<int>(_size))) ? (static_cast<int>(_size)-static_cast<int>((!std::is_same<AmountType, std::true_type>::value) && (static_cast<int>(first)+static_cast<int>(amount) < 0))) : ((static_cast<int>(first)+static_cast<int>(amount) < 0) ? (-static_cast<int>((!std::is_same<AmountType, std::true_type>::value) && (static_cast<int>(amount) < 0))) : (static_cast<int>(first)+static_cast<int>(amount)));
    unsigned int result = _size;
    if (_generator) {
        for (int i = ibegin; i != iend; i += idiff) {
            if (f(_root+_generator(_offset+i))) {
                result = i;
                break;
            }
        }
    } else if (!_format.empty()) {
        if (_length != 0) {
            for (int i = ibegin; i != iend; i += idiff) {
                std::snprintf(buffer.data(), bsize, _format.data(), _length, _offset+i);
                if (f(_root+buffer.data())) {
                    result = i;
                    break;
                }
            }
        } else {
            for (int i = ibegin; i != iend; i += idiff) {
                std::snprintf(buffer.data(), bsize, _format.data(), _offset+i);
                if (f(_root+buffer.data())) {
                    result = i;
                    break;
                }
            }
        }
    } else {
        for (int i = ibegin; i != iend; i += idiff) {
            if (f(_root+_container[i])) {
                result = i;
                break;
            }
        }
    }
    return result;
}

// Sort file names
/// \brief          Sort file names.
/// \details        Sorts the file names using the provided comparator function
///                 if the file list is based on a container. If the file list 
///                 is based on a format, it is checked whether the two first 
///                 files are sorted and an error is thrown if they are not. 
///                 Finally, if the file list is based on a generator function, 
///                 the order is checked and an error is thrown in case of non 
///                 sorted list.
/// \tparam         Function (Function type : <tt>bool(std::string, std::string)
///                 </tt>.)
/// \param[in]      f Function object <tt>bool(std::string, std::string)</tt>.
/// \return         Self reference.
/// \exception      std::runtime_error Generator based list cannot be modified.
/// \exception      std::runtime_error Constant format based list cannot be 
///                 modified.
template <class Function, class> 
FileList& FileList::sort(Function&& f)
{
    bool result = true;
    if (_size > 1) {
        if (_generator) {
            for (unsigned int i = 1; (i < _size) && (result); ++i) {
                result = !f(_generator(_offset+i), _generator(_offset+i-1));
            }
            if (!result) {
                throw std::runtime_error("ERROR = FileList::sort() : generator based list cannot be modified");
            }
        } else if (!_format.empty()) {
            if ((_length != 0) ? (f(generate(_offset+1, _format, _length), generate(_offset, _format, _length))) : (f(generate(_offset+1, _format), generate(_offset, _format)))) {
                throw std::runtime_error("ERROR = FileList::sort() : constant format based list cannot be modified");
            }
        } else {
            (std::is_same<Function, std::less<std::string> >::value) ? (std::sort(_container.begin(), _container.end())) : (std::sort(_container.begin(), _container.end(), f));
        }
    }
    return *this;
}

// Check sorting
/// \brief          Check sorting.
/// \details        Checks whether the file names are sorted using the provided 
///                 comparator function.
/// \tparam         Function (Function type : <tt>bool(std::string, std::string)
///                 </tt>.)
/// \param[in]      f Function object <tt>bool(std::string, std::string)</tt>.
/// \return         True if the contents is sorted, false otherwise.
template <class Function, class> 
bool FileList::sorted(Function&& f) const
{
    bool result = ((_generator) || (_size <= 1)) ? (true) : ((!_format.empty()) ? ((_length != 0) ? (!f(generate(_offset+1, _format, _length), generate(_offset, _format, _length))) : (!f(generate(_offset+1, _format), generate(_offset, _format)))) : ((std::is_same<Function, std::less<std::string> >::value) ? (std::is_sorted(_container.begin(), _container.end())) : (std::is_sorted(_container.begin(), _container.end(), f))));
    if (_size > 1) {
        if (_generator) {
            for (unsigned int i = 1; (i < _size) && (result); ++i) {
                result = !f(_generator(_offset+i), _generator(_offset+i-1));
            }
        }
    }
    return result;
}

// Unique file names
/// \brief          Unique file names.
/// \details        Erases consecutive equal file names using the provided 
///                 comparator function if the file list is based on a 
///                 container. If the file list is based on a format, it is 
///                 checked whether the two first files are unique and 
///                 if they are not, the file size is set equals to one. 
///                 Finally, if the file list is based on a generator function, 
///                 the unicity is checked and an error is thrown in case of non
///                 unique names.
/// \tparam         Function (Function type : <tt>bool(std::string, std::string)
///                 </tt>.)
/// \param[in]      f Function object <tt>bool(std::string, std::string)</tt>.
/// \return         Self reference.
/// \exception      std::runtime_error Generator based list cannot be modified.
template <class Function, class> 
FileList& FileList::unique(Function&& f)
{
    bool result = true;
    if (_size > 1) {
        if (_generator) {
            for (unsigned int i = 1; (i < _size) && (result); ++i) {
                result = !f(_generator(_offset+i), _generator(_offset+i-1));
            }
            if (!result) {
                throw std::runtime_error("ERROR = FileList::unique() : generator based list cannot be modified");
            }
        } else if (!_format.empty()) {
            if ((_length != 0) ? (f(generate(_offset+1, _format, _length), generate(_offset, _format, _length))) : (f(generate(_offset+1, _format), generate(_offset, _format)))) {
                _size = 1;
            }
        } else {
            _container.erase((std::is_same<Function, std::equal_to<std::string> >::value) ? (std::unique(_container.begin(), _container.end())) : (std::unique(_container.begin(), _container.end(), f)), _container.end());
            _size = _container.size();
        }
    }
    return *this;
}

// Check unicity
/// \brief          Check unicity.
/// \details        Checks whether the file names are unique using the provided 
///                 comparator function.
/// \tparam         Function (Function type : <tt>bool(std::string, std::string)
///                 </tt>.)
/// \param[in]      f Function object <tt>bool(std::string, std::string)</tt>.
/// \return         True if the contents is unique, false otherwise.
template <class Function, class> 
bool FileList::unicity(Function&& f) const
{
    bool result = ((_generator) || (_size <= 1)) ? (true) : ((!_format.empty()) ? ((_length != 0) ? (!f(generate(_offset+1, _format, _length), generate(_offset, _format, _length))) : (!f(generate(_offset+1, _format), generate(_offset, _format)))) : (true));
    if (_size > 1) {
        if (_generator) {
            for (unsigned int i = 1; (i < _size) && (result); ++i) {
                result = !f(_generator(_offset+i), _generator(_offset+i-1));
            }
        } else if (_format.empty()) {
            for (unsigned int i = 1; (i < _size) && (result); ++i) {
                result = !f(_container[i], _container[i-1]);
            }
        }
    }
    return result;
}
// -------------------------------------------------------------------------- //



// --------------------------------- STREAM --------------------------------- //
// Output stream operator
/// \brief          Output stream operator.
/// \details        Adds each file name to the stream using the filling 
///                 character as a separator.
/// \param[in,out]  lhs Left-hand side stream.
/// \param[in]      rhs Right-hand side file list.
/// \return         Output stream.
std::ostream& operator<<(std::ostream& lhs, const FileList& rhs)
{
    const char separator = lhs.fill();
    const unsigned int n = rhs._size-(rhs._size > 0);
    for (unsigned int i = 0; i < n; ++i) {
        lhs<<rhs[i]<<separator; 
    }
    return (rhs._size > 0) ? (lhs<<rhs[n]) : (lhs);
}
// -------------------------------------------------------------------------- //



// --------------------------------- HELPERS -------------------------------- //
// Temporary file name generator
/// \brief          Temporary file name generator.
/// \details        Generates a temporary file name based on the specified 
///                 index.
/// \tparam         IndexType (Integral type corresponding to the index.)
/// \param[in]      i Index.
/// \return         File name.
template <typename IndexType, class>
inline std::string FileList::generate(const IndexType i)
{
    char buffer[L_tmpnam];
    return (std::is_same<IndexType, std::true_type>::value) ? (std::string(std::tmpnam(buffer))) : (std::string(std::tmpnam(buffer))+std::to_string(i));
}
    
// Container-based file name generator
/// \brief          Container-based file name generator.
/// \details        Gets the specified element of the passed container.
/// \tparam         ContainerType (Container type.)
/// \param[in]      i Index.
/// \param[in]      source Source container.
/// \return         File name.
template <class ContainerType, class>
inline std::string FileList::generate(const unsigned int i, const ContainerType& source)
{
    return (std::begin(source)+i < std::end(source)) ? (*(std::begin(source)+i)) : (std::string());
}

// Format-based file name generator
/// \brief          Format-based file name generator.
/// \details        Generates a file name using the specified C-style format.
///                 If this format requires a width specifier, it can be passed
///                 as an optional parameter.
/// \tparam         BufferSize Value of the internal buffer size.
/// \tparam         WidthType (Integral type corresponding to the format width.)
/// \param[in]      i Index.
/// \param[in]      cformat C-like format string.
/// \param[in]      width Optional width of the format field.
/// \return         File name.
template <unsigned int BufferSize, typename WidthType, class>
inline std::string FileList::generate(const unsigned int i, const std::string& cformat, const WidthType width)
{
    char buffer[BufferSize];
    (void)((std::is_same<WidthType, std::true_type>::value) ? (std::snprintf(buffer, BufferSize, cformat.data(), i)) : (std::snprintf(buffer, BufferSize, cformat.data(), width, i)));
    return buffer;
}

// Function-based file name generator
/// \brief          Format-based file name generator.
/// \details        Generates a file name passing the index to the specified
///                 function.
/// \tparam         Function (Function type : <tt>std::string(const unsigned 
///                 int)</tt>.)
/// \param[in]      i Index.
/// \param[in]      f Function object <tt>std::string(const unsigned int)</tt>.
/// \return         File name.
template <class Function, class>
inline std::string FileList::generate(const unsigned int i, Function&& f)
{
    return f(i);
}

// Parse number in a string
/// \brief          Parse number in a string.
/// \details        Detects positions of specific markers around a digit in a 
///                 string. If the specified position is greater than the string
///                 length, it is first adjusted to the last character of the 
///                 string. If the character at the provided position is not a 
///                 digit, the function starts by searching the last digit 
///                 before the position. In case of failure, all the markers are
///                 set equal to the end of the string.  
/// \tparam         Dir (Directory type.)
/// \tparam         Min (Minimum type.)
/// \tparam         Num (Number type.)
/// \tparam         Max (Maximum type.)
/// \param[in]      str Input string containing the number to parse.
/// \param[in]      pos Position of a digit inside the string.
/// \param[out]     dir First character after the last slash before the 
///                 specified position.
/// \param[out]     min First digit of the number around the specified position,
///                 including preceding zeroes.
/// \param[out]     num First digit of the number around the specified position,
///                 excluding preceding zeroes.
/// \param[out]     max First position after the number around the specified 
///                 position.
/// \return         Value of the number around the specified position.
template <typename Dir, typename Min, typename Num, typename Max, class>
unsigned int FileList::numberify(const std::string& str, unsigned int pos, Dir&& dir, Min&& min, Num&& num, Max&& max)
{
    const unsigned int siz = str.size();
    unsigned int result = 0;
    pos = (siz > 0) ? ((pos > siz) ? (siz-1) : (pos)) : (std::numeric_limits<unsigned int>::max());
    if (pos != std::numeric_limits<unsigned int>::max()) {
        pos = siz-(std::find_if(str.rend()-pos-1, str.rend(), [](const char& c){return std::isdigit(c);})-str.rbegin());
        pos = (pos > 0) ? ((std::isdigit(str[pos-1])) ? (pos-1) : (std::numeric_limits<unsigned int>::max())) : (std::numeric_limits<unsigned int>::max());
    }
    if (pos != std::numeric_limits<unsigned int>::max()) {
        min = siz-(std::find_if_not(str.rend()-pos-1, str.rend(), [](const char& c){return std::isdigit(c);})-str.rbegin());
        max = std::find_if_not(str.begin()+pos, str.end(), [](const char& c){return std::isdigit(c);})-str.begin();
        num = std::find_if_not(str.begin()+min, str.begin()+max, [](const char& c){return c == '0';})-str.begin();
        dir = siz-(std::find_if(str.rend()-min, str.rend(), [](const char& c){return ((c == '/') || (c == '\\'));})-str.rbegin());
        result = std::stoul(str.substr(min));
        num -= (result == 0);
    } else {
        dir = siz;
        min = siz;
        num = siz;
        max = siz;
    }
    return result;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of FileList.
/// \return         0 if no error.
int FileList::example()
{
    // Initialize
    std::cout<<"BEGIN = FileList::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    std::array<std::string, 4> array = {{"the", "answer", "is", "42"}};
    
    // Construction
    FileList list([](const unsigned int i){return std::to_string(i);}, 4);
    
    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Lifecycle : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"FileList() : "                                                                      <<FileList()<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList({{\"the\", \"answer\", \"is\", \"42\"}}) : "                               <<FileList({{"the", "answer", "is", "42"}})<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(array) : "                                                                 <<FileList(array)<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(array, \"/directory/file_\") : "                                           <<FileList(array, "/directory/file_")<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(\"file_%04d\") : "                                                         <<FileList("file_%04d")<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(\"file_%04d\", 0, 4) : "                                                   <<FileList("file_%04d", 0, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(\"file_%04d\", 0, 4, 4) : "                                                <<FileList("file_%04d", 0, 4, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(\"file_%0*d\", 4, 4, 4) : "                                                <<FileList("file_%0*d", 4, 4, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList(\"file_%0*d\", 4, 4, 4, \"/directory/\") : "                               <<FileList("file_%0*d", 4, 4, 4, "/directory/")<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList([](int i){return \"file_\"+std::to_string(i);}) : "                        <<FileList([](int i){return "file_"+std::to_string(i);})<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList([](int i){return \"file_\"+std::to_string(i);}, 4) : "                     <<FileList([](int i){return "file_"+std::to_string(i);}, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList([](int i){return \"file_\"+std::to_string(i);}, 4, 4) : "                  <<FileList([](int i){return "file_"+std::to_string(i);}, 4, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"FileList([](int i){return \"file_\"+std::to_string(i);}, 4, 4, \"/dir/\") : "       <<FileList([](int i){return "file_"+std::to_string(i);}, 4, 4, "/dir/")<<std::endl;

    // Operators
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Operators : "                                                                       <<std::endl;
    std::cout<<std::setw(width*2)<<"list = FileList() : "                                                               <<(list = FileList())<<std::endl;
    std::cout<<std::setw(width*2)<<"list = {{\"the\", \"answer\", \"is\", \"42\"}} : "                                  <<(list = {{"the", "answer", "is", "42"}})<<std::endl;
    std::cout<<std::setw(width*2)<<"list = array : "                                                                    <<(list = array)<<std::endl;
    std::cout<<std::setw(width*2)<<"list = [](int i){return \"file_\"+std::to_string(i);} : "                           <<(list = [](int i){return "file_"+std::to_string(i);})<<std::endl;
    std::cout<<std::setw(width*2)<<"list == FileList(\"file_%d\") : "                                                   <<(list == FileList("file_%d"))<<std::endl;
    std::cout<<std::setw(width*2)<<"list == FileList(\"file_%05d\") : "                                                 <<(list == FileList("file_%05d"))<<std::endl;
    std::cout<<std::setw(width*2)<<"list != FileList(\"file_%d\") : "                                                   <<(list != FileList("file_%d"))<<std::endl;
    std::cout<<std::setw(width*2)<<"list != FileList(\"file_%05d\") : "                                                 <<(list != FileList("file_%05d"))<<std::endl;
            
    // Assignment
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Assignment : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(list) : "                                                               <<list.assign(list)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(FileList()) : "                                                         <<list.assign(FileList())<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign() : "                                                                   <<list.assign()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign({{\"the\", \"answer\", \"is\", \"42\"}}) : "                            <<list.assign({{"the", "answer", "is", "42"}})<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(array) : "                                                              <<list.assign(array)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(array, \"/directory/file_\") : "                                        <<list.assign(array, "/directory/file_")<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(\"file_%04d\") : "                                                      <<list.assign("file_%04d")<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(\"file_%04d\", 0, 4) : "                                                <<list.assign("file_%04d", 0, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(\"file_%04d\", 0, 4, 4) : "                                             <<list.assign("file_%04d", 0, 4, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(\"file_%0*d\", 4, 4, 4) : "                                             <<list.assign("file_%0*d", 4, 4, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign(\"file_%0*d\", 4, 4, 4, \"/directory/\") : "                            <<list.assign("file_%0*d", 4, 4, 4, "/directory/")<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign([](int i){return \"file_\"+std::to_string(i);}) : "                     <<list.assign([](int i){return "file_"+std::to_string(i);})<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign([](int i){return \"file_\"+std::to_string(i);}, 4) : "                  <<list.assign([](int i){return "file_"+std::to_string(i);}, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign([](int i){return \"file_\"+std::to_string(i);}, 4, 4) : "               <<list.assign([](int i){return "file_"+std::to_string(i);}, 4, 4)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.assign([](int i){return \"file_\"+std::to_string(i);}, 4, 4, \"/dir/\") : "    <<list.assign([](int i){return "file_"+std::to_string(i);}, 4, 4, "/dir/")<<std::endl;

    // Management
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Management : "                                <<std::endl;
    std::cout<<std::setw(width)<<"list.empty() : "                              <<list.empty()<<std::endl;
    std::cout<<std::setw(width)<<"list.cleared() : "                            <<list.cleared()<<std::endl;
    std::cout<<std::setw(width)<<"list.capacity() : "                           <<list.capacity()<<std::endl;
    std::cout<<std::setw(width)<<"list.resize(8) : "                            <<list.resize(8)<<std::endl;
    std::cout<<std::setw(width)<<"list.reserve(list.capacity()*2) : "           <<list.reserve(list.capacity()*2)<<std::endl;
    std::cout<<std::setw(width)<<"list.container().capacity() : "               <<list.container().capacity()<<std::endl;
    std::cout<<std::setw(width)<<"list.shrink() : "                             <<list.shrink()<<std::endl;
    std::cout<<std::setw(width)<<"list.clear() : "                              <<list.clear()<<std::endl;
    std::cout<<std::setw(width)<<"list.cleared() : "                            <<list.cleared()<<std::endl;
    std::cout<<std::setw(width)<<"list.copy() : "                               <<list.copy()<<std::endl;
    std::cout<<std::setw(width)<<"list.cast() : "                               <<list.cast()<<std::endl;
    
    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"list.container().size() : "                   <<list.container().size()<<std::endl;
    std::cout<<std::setw(width)<<"list.format() : "                             <<list.format()<<std::endl;
    std::cout<<std::setw(width)<<"bool(list.generator()) : "                    <<bool(list.generator())<<std::endl;
    std::cout<<std::setw(width)<<"list.format() : "                             <<list.format()<<std::endl;
    std::cout<<std::setw(width)<<"list.root() : "                               <<list.root()<<std::endl;
    std::cout<<std::setw(width)<<"list.length() : "                             <<list.length()<<std::endl;
    std::cout<<std::setw(width)<<"list.size() : "                               <<list.size()<<std::endl;
    std::cout<<std::setw(width)<<"list.offset() : "                             <<list.offset()<<std::endl;
    std::cout<<std::setw(width)<<"list.assign(array) : "                        <<list.assign(array)<<std::endl;

    // Files
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Files : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"list[0] : "                                   <<list[0]<<std::endl;
    std::cout<<std::setw(width)<<"list[4] : "                                   <<list[4]<<std::endl;
    std::cout<<std::setw(width)<<"list.at(0) : "                                <<list.at(0)<<std::endl;
    std::cout<<std::setw(width)<<"list.front() : "                              <<list.front()<<std::endl;
    std::cout<<std::setw(width)<<"list.front(1) : "                             <<list.front(1)<<std::endl;
    std::cout<<std::setw(width)<<"list.back(1) : "                              <<list.back(1)<<std::endl;
    std::cout<<std::setw(width)<<"list.back() : "                               <<list.back()<<std::endl;
    std::cout<<std::setw(width)<<"list.back() : "                               <<list.back()<<std::endl;
    std::cout<<std::setw(width)<<"list.assign(\"file_%05d\") : "                <<list.assign("file_%05d")<<std::endl;
    std::cout<<std::setw(width)<<"list(42) : "                                  <<list(42)<<std::endl;
    std::cout<<std::setw(width)<<"list.get(42) : "                              <<list.get(42)<<std::endl;
    std::cout<<std::setw(width)<<"list.assign(\"file_%05d\", 0, 10) : "         <<list.assign("file_%05d", 0, 10)<<std::endl;
    
    // Manipulations
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Manipulations : "                             <<std::endl;
    std::cout<<std::setw(width)<<"list.convert()[0] : "                         <<list.convert()[0]<<std::endl;
    std::cout<<std::setw(width)<<"list.common() : "                             <<list.common()<<std::endl;
    std::cout<<std::setw(width)<<"list.formatify() : "                          <<list.formatify()<<std::endl;
    
    // Algorithms
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Algorithms : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"list.apply([](const std::string& filename){return filename.size();}) : "            <<list.apply([](const std::string& filename){return filename.size();})<<std::endl;
    std::cout<<std::setw(width*2)<<"list.count() : "                                                                    <<list.count()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.count([](const std::string& filename){return filename.size() > 4;}) : "        <<list.count([](const std::string& filename){return filename.size() > 4;})<<std::endl;
    std::cout<<std::setw(width*2)<<"list.find([](const std::string& filename){return filename.size() > 4;}) : "         <<list.find([](const std::string& filename){return filename.size() > 4;})<<std::endl;
    std::cout<<std::setw(width*2)<<"list.find([](const std::string& filename){return filename.size()>16;},4,-9):"       <<list.find([](const std::string& filename){return filename.size() > 16;}, 4, -9)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.sort() : "                                                                     <<list.sort()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.sort(std::less<std::string>()) : "                                             <<list.sort(std::less<std::string>())<<std::endl;
    std::cout<<std::setw(width*2)<<"list.sorted() : "                                                                   <<list.sorted()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.sorted(std::less<std::string>()) : "                                           <<list.sorted(std::less<std::string>())<<std::endl;
    std::cout<<std::setw(width*2)<<"list.unique() : "                                                                   <<list.unique()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.unique(std::equal_to<std::string>()) : "                                       <<list.unique(std::equal_to<std::string>())<<std::endl;
    std::cout<<std::setw(width*2)<<"list.unicity() : "                                                                  <<list.unicity()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.unicity(std::equal_to<std::string>()) : "                                      <<list.unicity(std::equal_to<std::string>())<<std::endl;
    
    // Stream
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Stream : "                                    <<std::endl;
    std::cout<<std::setw(width)<<"operator<<(std::cout, list) : "               <<list<<std::endl;
    
    // Helpers
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Helpers : "                                                                         <<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate() : "                                                                 <<list.generate()<<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate(42) : "                                                               <<list.generate(42)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate(3, array) : "                                                         <<list.generate(3, array)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate(4, array) : "                                                         <<list.generate(4, array)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate(42, \"file_%08d\") : "                                                <<list.generate(42, "file_%08d")<<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate(42, \"file_%0*d\", 8) : "                                             <<list.generate(42, "file_%0*d", 8)<<std::endl;
    std::cout<<std::setw(width*2)<<"list.generate(42, [](const unsigned int i){return std::to_string(i);}) : "          <<list.generate(42, [](const unsigned int i){return std::to_string(i);})<<std::endl;
    std::cout<<std::setw(width*2)<<"list.numberify(\"/directory/file_0481516_02342.txt\") : "                           <<list.numberify("/directory/file_0481516_02342.txt")<<std::endl;
    std::cout<<std::setw(width*2)<<"list.numberify(\"/directory/file_0481516_02342.txt\", 20) : "                       <<list.numberify("/directory/file_0481516_02342.txt", 20)<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = FileList::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // FILELIST_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
