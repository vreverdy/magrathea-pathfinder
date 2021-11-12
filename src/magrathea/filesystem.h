/* ******************************* FILESYSTEM ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA-PATHFINDER
// TITLE :          FileSystem
// DESCRIPTION :    Global file management
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           filesystem.h
/// \brief          Global file management
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef FILESYSTEM_H_INCLUDED
#define FILESYSTEM_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <type_traits>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdint>
#include <ctime>
#include <chrono>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <limits>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Global file management
/// \brief          Global file management.
/// \details        Provides general functions to manage files, retrieve some
///                 information and perform standard operations on files. All
///                 functions are standard-compliant, but use either C++ or C
///                 depending on performances. Most functions return whether the
///                 operation is a success and do not throw any exception. In 
///                 the class, an unit means the size in bytes of a provided 
///                 type. 
class FileSystem final
{
    // Utilities
    /// \name           Utilities
    //@{
    public:
        static inline bool endianness(); 
        template <bool Byteswap, typename Type, class = typename std::enable_if<!Byteswap>::type> static inline bool byteswap(const Type& variable);
        template <bool Byteswap = true, typename Type, class = typename std::enable_if<Byteswap>::type> static inline bool byteswap(Type& variable);
        template <typename Type = unsigned int, class = typename std::enable_if<std::is_integral<Type>::value>::type> static inline Type bom();
        template <typename Type, class = typename std::enable_if<std::is_integral<Type>::value>::type> static bool bom(const Type& mark);
    //@}

    // Size
    /// \name           Size
    //@{
    public:
        template <typename Type = char, class Stream> static inline long long int weight(Stream& stream);
        static long long int size(const std::string& filename);
        static long long int size(std::istream& stream);
        static long long int size(std::ostream& stream);
    //@}
    
    // Management
    /// \name           Management
    //@{
    public:
        static inline std::string temporary(const std::string& prefix = "\b", const std::string& suffix = "");
        static inline std::string dated(const std::string& prefix = "", const std::string& suffix = "", const std::string& format = "%Y-%m-%d-%H-%M-%S");
        static inline bool remove(const std::string& filename);
        static inline bool rename(const std::string& oldname, const std::string& newname, const std::ios::openmode& mode = std::ios::out);
        static bool copy(const std::string& oldname, const std::string& newname, const std::ios::openmode& mode = std::ios::out, const long long int chunk = -1);
    //@}

    // Split and join
    /// \name           Split and join
    //@{
    public:
        template <bool Byteswap = false, typename Marker = long long int, class Container> static unsigned int split(const std::string& filename, const Container& filenames, const std::ios::openmode& mode = std::ios::out, const long long int chunk = -1, const unsigned long long int limit = 0);
        template <bool Byteswap = false, typename Marker = long long int, class Container> static unsigned int unsplit(const Container& filenames, const std::string& filename, const std::ios::openmode& mode = std::ios::out, const long long int chunk = -1);
        template <bool Byteswap = false, typename Marker = long long int, class Container> static unsigned int join(const Container& filenames, const std::string& filename, const std::ios::openmode& mode = std::ios::out, const long long int chunk = -1);
        template <bool Byteswap = false, typename Marker = long long int, class Container> static unsigned int unjoin(const std::string& filename, const Container& filenames, const std::ios::openmode& mode = std::ios::out, const long long int chunk = -1);
    //@}
    
    // Comparison
    /// \name           Comparison
    //@{
    public:
        static bool compare(const std::string& first, const std::string& second, const long long int chunk = -1);
        static bool compare(std::istream& first, std::istream& second, const long long int chunk = -1);
        static bool compare(std::ostream& first, std::ostream& second, const long long int chunk = -1);
    //@}

    // Existence and creation
    /// \name           Existence and creation
    //@{
    public:
        static inline bool exist(const std::string& filename);
        template <typename Type = char> static bool check(const std::string& filename, const long long int min = -1, const long long int max = -1);
        template <class Contents = std::true_type> static bool create(const std::string& filename, const Contents& contents = Contents());
        template <class Contents = std::true_type> static bool initialize(const std::string& filename, const Contents& contents = Contents());
        template <class Contents = std::true_type> static bool reset(const std::string& filename, const Contents& contents = Contents());
        template <class Contents = char, class = typename std::enable_if<!std::is_convertible<Contents, std::string>::value>::type> static bool generate(const std::string& filename, const std::ios::openmode& mode, const long long int amount, const long long int chunk = -1, const Contents& contents = Contents());
        template <class Contents, class = typename std::enable_if<std::is_convertible<Contents, std::string>::value>::type, class = void> static bool generate(const std::string& filename, const std::ios::openmode& mode, const long long int amount, const long long int chunk, const Contents& contents);
        template <class Engine, class Distribution, class Contents = typename std::decay<Distribution>::type::result_type, class = typename std::enable_if<(std::decay<Engine>::type::min() != std::decay<Engine>::type::max()) && (!std::is_void<typename std::decay<Distribution>::type::result_type>::value)>::type> static bool generate(const std::string& filename, const std::ios::openmode& mode, const long long int amount, const long long int chunk, Engine&& engine, Distribution&& distribution);
    //@}

    // Size control
    /// \name           Size control
    //@{
    public:
        template <typename Type = char, class Stream> static inline bool empty(Stream& stream);
        template <typename Type = char, class Stream> static inline bool exact(Stream& stream, const long long int amount);
        template <typename Type = char, class Stream> static inline bool regular(Stream& stream, const long long int min = -1, const long long int max = -1);
    //@}
    
    // File types
    /// \name           File types
    //@{
        static bool ascii(const std::string& filename, const long long int min = -1, const long long int max = -1, const long long int chunk = -1);
        static bool eascii(const std::string& filename, const long long int min = -1, const long long int max = -1, const long long int chunk = -1);
        template <typename Type = char> static bool binary(const std::string& filename, const long long int min = -1, const long long int max = -1, const long long int chunk = -1); 
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
// Get system endianness
/// \brief          Get system endianness.
/// \details        Returns the system endianness tested with an integer.
/// \return         True for big-endian, false for little-endian.
inline bool FileSystem::endianness()
{
    static const union {uint32_t i; char c[4];} x = {0x0000FEFF};
    return !(x.c[0]);
}

// Do not swap bytes
/// \brief          Do not swap bytes.
/// \details        Does not invert the order of bytes of the passed variable to
///                 keep the endianness.
/// \tparam         Byteswap Do not swap endianness if false.
/// \tparam         Type (Variable type.)
/// \param[in]      variable Variable.
/// \return         False if the variable has not been byteswapped.
template <bool Byteswap, typename Type, class>
inline bool FileSystem::byteswap(const Type& variable)
{
    return Byteswap && sizeof(variable);
}

// Swap bytes
/// \brief          Swap bytes.
/// \details        Inverts the order of bytes of the passed variable to change 
///                 the endianness.
/// \tparam         Byteswap Swap endianness if true.
/// \tparam         Type (Variable type.)
/// \param[in,out]  variable Variable.
/// \return         True if the variable has been byteswapped.
template <bool Byteswap, typename Type, class>
inline bool FileSystem::byteswap(Type& variable)
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
    return Byteswap;
}

// Byte order mark
/// \brief          Byte order mark.
/// \details        Returns the byte order mark 0xFEFF casted to the provided
///                 type. It is used to detect the endianness : for example, for
///                 a 4-bytes integer : 00-00-FE-FF indicates big-endianness and
///                 FF-FE-00-00 indicates little-endianness.
/// \tparam         Type Integral byte order mark type.
/// \return         Cast of 0xFEFF.
template <typename Type, class> 
inline Type FileSystem::bom()
{
    return static_cast<Type>(0xFEFF);
}

// Byte order mark endianness
/// \brief          Byte order mark endianness.
/// \details        Returns endianness from byte order mark passed as parameter 
///                 : true for big-endian, false for little-endian. The byte 
///                 order mark should be of the form FE-FF.
/// \tparam         Type (Integral byte order mark type.)    
/// \param[in]      mark Byte order mark.
/// \return         True for big-endian, false for little-endian.
/// \exception      std::invalid_argument Unrecognized byte order mark.
template <typename Type, class>
bool FileSystem::bom(const Type& mark)
{
    static const char fe = static_cast<char>(0xFE);
    static const char ff = static_cast<char>(0xFF);
    union {Type i; char c[sizeof(mark)];} little = {0};
    union {Type i; char c[sizeof(mark)];} big = {0};
    little.c[sizeof(mark) > 1] = fe;
    little.c[0] = ff;
    big.c[sizeof(mark)-1-(sizeof(mark) > 1)] = fe;
    big.c[sizeof(mark)-1] = ff;
    if ((sizeof(mark) <= 1) || ((mark != little.i) && (mark != big.i))) {
        throw std::invalid_argument("ERROR = FileSystem::bom() : unrecognized byte order mark");
    }
    return (mark == big.i);
}
// -------------------------------------------------------------------------- //



// ---------------------------------- SIZE ---------------------------------- //
// Get the weight of a file in terms of a specific type
/// \brief          Get the weight of a file in terms of a specific type.
/// \details        Returns the amount of data of the specified type equivalent
///                 to the file size. This is just the file size divided by the
///                 size of the specified type rounded to the lower integer. If 
///                 the file does not exist or if a stream error is detected, 
///                 -1 is returned.
/// \tparam         Type Type representing the considered unit.
/// \tparam         Stream (String, input stream or output stream.)
/// \param[in,out]  stream Stream.
/// \return         File weight in unit or -1 if error.
template <typename Type, class Stream> 
inline long long int FileSystem::weight(Stream& stream)
{
    long long int length = size(stream);
    return (length >= 0) ? (length/static_cast<long long int>(sizeof(Type))) : (-1);
}

// Get size of a file passed by name
/// \brief          Get size of a file passed by name.
/// \details        Opens the file, computes its size and closes it. If the 
///                 file does not exist or if a stream error is detected, 
///                 -1 is returned.
/// \param[in]      filename File name.
/// \return         File size in bytes or -1 if error.
long long int FileSystem::size(const std::string& filename)
{
    std::ifstream stream(filename, std::ios::binary);
    std::streampos first = 0;
    std::streampos last = 0;
    long long int length = -1;
    if (stream) {
        stream.seekg(0, std::ios::beg);
        first = stream.tellg();
        stream.seekg(0, std::ios::end);
        last = stream.tellg();
        if ((stream) && (first != -1) && (last != -1)) {
            length = last-first;
        }
    }
    return length;
}

// Get size of a file passed by input stream
/// \brief          Get size of a file passed by input stream.
/// \details        Saves the current position, computes the size of the passed
///                 stream and returns to the original position. If the file 
///                 does not exist or if a stream error is detected, -1 is 
///                 returned.
/// \param[in,out]  stream Input stream.
/// \return         File size in bytes or -1 if error.
long long int FileSystem::size(std::istream& stream)
{
    std::streampos position = stream.tellg();
    std::streampos first = 0;
    std::streampos last = 0;
    long long int length = -1;
    if (stream) {
        stream.seekg(0, std::ios::beg);
        first = stream.tellg();
        stream.seekg(0, std::ios::end);
        last = stream.tellg();
        stream.seekg(position);
        if ((stream) && (first != -1) && (last != -1)) {
            length = last-first;
        }
    }
    return length;
}

// Get size of a file passed by output stream
/// \brief          Get size of a file passed by output stream.
/// \details        Saves the current position, computes the size of the passed
///                 stream and returns to the original position. If the file 
///                 does not exist or if a stream error is detected, -1 is 
///                 returned.
/// \param[in,out]  stream Output stream.
/// \return         File size in bytes or -1 if error.
long long int FileSystem::size(std::ostream& stream)
{
    std::streampos position = stream.tellp();
    std::streampos first = 0;
    std::streampos last = 0;
    long long int length = -1;
    if (stream) {
        stream.seekp(0, std::ios::beg);
        first = stream.tellp();
        stream.seekp(0, std::ios::end);
        last = stream.tellp();
        stream.seekp(position);
        if ((stream) && (first != -1) && (last != -1)) {
            length = last-first;
        }
    }
    return length;
}
// -------------------------------------------------------------------------- //



// ------------------------------- MANAGEMENT ------------------------------- //
// Temporary file name
/// \brief          Temporary file name.
/// \details        Generates a temporary file name. If no argument is used, 
///                 then the default location is used. If an empty prefix 
///                 is specified, the path is erased, and only the file name is 
///                 kept. If a prefix or a suffix is specified, the default 
///                 path is erased, the file name is kept and prefixed and 
///                 suffixed by the arguments.
/// \param[in]      prefix File name prefix or path.
/// \param[in]      suffix File name suffix or extension.
/// \return         Generated temporary file name.
inline std::string FileSystem::temporary(const std::string& prefix, const std::string& suffix)
{
    static const std::string defaulted = "\b";
    static const std::string directory = "/\\";
    char buffer[L_tmpnam];
    std::string filename = std::tmpnam(buffer);
    return (prefix != defaulted) ? (prefix+filename.erase(0, std::distance(std::find_first_of(filename.rbegin(), filename.rend(), directory.rbegin(), directory.rend()), filename.rend()))+suffix) : (filename);
}

// Dated file name
/// \brief          Dated file name.
/// \details        Generates a file name from the current time.
/// \param[in]      prefix File name prefix or path.
/// \param[in]      suffix File name suffix or extension.
/// \param[in]      format Date format compatible with <tt>strftime()</tt>.
/// \return         Generated dated file name.
inline std::string FileSystem::dated(const std::string& prefix, const std::string& suffix, const std::string& format)
{   
    static const unsigned int length = sizeof(unsigned long long int)*sizeof(unsigned long long int);
    const decltype(std::chrono::system_clock::to_time_t(std::chrono::system_clock::time_point())) now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char buffer[length];
    std::strftime(buffer, length, format.data(), std::localtime(&now));
    return prefix+buffer+suffix;
}

// Remove a file
/// \brief          Remove a file.
/// \details        Removes an existing file from the file system and returns 
///                 true on success.
/// \param[in]      filename File name.
/// \return         True on success, false on error.
bool FileSystem::remove(const std::string& filename)
{
    return (std::remove(filename.data()) == 0);
}

// Rename a file
/// \brief          Rename a file.
/// \details        Renames an existing file and returns true on success. If the 
///                 old and new names are equal, the function returns false 
///                 without further tests. If the new file name already exists,
///                 nothing is done and the function fails except if the 
///                 truncate open mode is specified : in that case, the existing
///                 file is erased.
/// \param[in]      oldname Old name.
/// \param[in]      newname New name.
/// \param[in]      mode Open mode.
/// \return         True on success, false on error.
bool FileSystem::rename(const std::string& oldname, const std::string& newname, const std::ios::openmode& mode)
{
    const bool ok = (newname != oldname) ? ((exist(newname)) ? (((mode) & (std::ios::trunc)) ? (std::remove(newname.data()) == 0) : (false)) : (true)) : (false);
    return (ok) ? (std::rename(oldname.data(), newname.data()) == 0) : (ok);
}

// Copy
/// \brief          Copy.
/// \details        Copy the provided file to a new location. If the old and new
///                 names are equal, the function returns false without further 
///                 tests. By default, if a standard output open mode is 
///                 specified, all existing files are protected and only a new 
///                 one can be created. This protection can be removed by 
///                 explicitely specifying the truncate open mode. Finally, the 
///                 chunk parameter allows to specify whether a buffer is used 
///                 for the copy : 0 corresponds to no buffer, a positive 
///                 integer corresponds to the amount of bytes that is put into 
///                 the buffer, and a negative integer corresponds to a buffer 
///                 of the total size of the file. 
/// \param[in]      oldname Old file name.
/// \param[in]      newname New file name.
/// \param[in]      mode Open mode.
/// \param[in]      chunk Buffer size.
/// \return         True if the copy was done without errors, false otherwise.    
bool FileSystem::copy(const std::string& oldname, const std::string& newname, const std::ios::openmode& mode, const long long int chunk)
{
    std::ifstream input(oldname, std::ios::binary);
    bool ok = (((mode) & (std::ios::trunc | std::ios::app)) || (!exist(newname))) && (newname != oldname) && (input);
    const long long int amount = (ok) ? (size(input)) : (-1);
    const long long int length = ((ok) && (amount >= 0)) ? (((chunk < 0) || (amount < chunk)) ? (amount) : (chunk)) : (0);
    const long long int niterations = (length > 0) ? (amount/length) : (0);
    char contents = char();
    std::vector<char> buffer(length);
    std::ofstream output;
    if (ok) {
        output.open(newname, mode | std::ios::binary);
        if ((input) && (output)) {
            if (length > 0) {
                for (long long int iteration = 0; iteration < niterations; ++iteration) {
                    input.read(buffer.data(), length*sizeof(*buffer.data()));
                    output.write(buffer.data(), length*sizeof(*buffer.data()));
                }
                input.read(buffer.data(), (amount%length)*sizeof(*buffer.data()));
                output.write(buffer.data(), (amount%length)*sizeof(*buffer.data()));
            } else {
                for (long long int byte = 0; byte < amount; ++byte) {
                    input.read(&contents, sizeof(contents));
                    output.write(&contents, sizeof(contents));
                }
            }
        }
        ok = (input.good()) && (output.good());
    }
    return ok;
}
// -------------------------------------------------------------------------- //



// ----------------------------- SPLIT AND JOIN ----------------------------- //
// Split a file in several files
/// \brief          Split a file in several files.
/// \details        Splits the provided file into several ones of lower sizes.
///                 At the beginning and the end of each file, the current 
///                 iterator byte position regarding the splitted file is saved 
///                 as a marker. A beginning marker equals to zero corresponds 
///                 to the first and an end marker equals to zero corresponds
///                 to the last file. The limit parameter allows to limit the 
///                 size by file to an exact number of bytes. By default, if a 
///                 standard output open mode is specified, all existing files 
///                 are protected and only new ones can be created. This 
///                 protection can be removed by explicitely specifying the 
///                 truncate open mode. Finally, the chunk parameter allows to 
///                 specify whether a buffer is used for the copy : 0 
///                 corresponds to no buffer, a positive integer corresponds to 
///                 the amount of bytes that is put into the buffer, and a 
///                 negative integer corresponds to a buffer of the total size 
///                 of the file. 
/// \tparam         Byteswap Swap endianness of markers.
/// \tparam         Marker Record marker type.
/// \tparam         Container (File container type.)
/// \param[in]      filename Input file name to split.
/// \param[in]      filenames Output file names.
/// \param[in]      mode Open mode.
/// \param[in]      chunk Buffer size.
/// \param[in]      limit Size limit of each file.
/// \return         Number of files written on success, zero on error.
template <bool Byteswap, typename Marker, class Container> 
unsigned int FileSystem::split(const std::string& filename, const Container& filenames, const std::ios::openmode& mode, const long long int chunk, const unsigned long long int limit)
{
    std::ifstream input(filename, std::ios::binary);
    unsigned int nfiles = filenames.size();
    bool ok = (input.good()) && (nfiles > 0);
    long long int total = (ok) ? (size(input)) : (0);
    long long int current = 0;
    long long int amount = 0;
    long long int length = 0;
    long long int niterations = 0;
    Marker marker = Marker();
    char contents = char();
    std::vector<char> buffer((chunk >= 0) ? (chunk) : ((limit > 0) ? (limit) : (total/static_cast<long long int>(nfiles+(nfiles == 0))+1)));
    std::ofstream output;
    if (ok) {
        if (!((mode) & (std::ios::trunc | std::ios::app))) {
            for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                ok = !exist(filenames[file]);
            }
        }
        if (ok) {
            if (limit > 0) {
                nfiles = total/static_cast<long long int>(limit)+(total%static_cast<long long int>(limit) > 0)+(total == 0);
                ok = (ok) && (nfiles > 0) && (nfiles <= filenames.size());
            }
            if (ok) {
                if (chunk < 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        output.open(filenames[file], mode | std::ios::binary);
                        amount = (limit > 0) ? ((file+1 < nfiles) ? (static_cast<long long int>(limit)) : (total-current)) : ((total/static_cast<long long int>(nfiles))+(file < total%static_cast<long long int>(nfiles)));
                        marker = current;
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        input.read(buffer.data(), amount*sizeof(*buffer.data()));
                        output.write(buffer.data(), amount*sizeof(*buffer.data()));
                        current += amount;
                        marker = current*(file+1 < nfiles);
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        output.close();
                        ok = (input.good()) && (output.good());
                    }
                } else if (chunk > 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        output.open(filenames[file], mode | std::ios::binary);
                        amount = (limit > 0) ? ((file+1 < nfiles) ? (static_cast<long long int>(limit)) : (total-current)) : ((total/static_cast<long long int>(nfiles))+(file < total%static_cast<long long int>(nfiles)));
                        length = (amount < chunk) ? (amount) : (chunk);
                        niterations = (length > 0) ? (amount/length) : (0);
                        marker = current;
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        for (long long int iteration = 0; iteration < niterations; ++iteration) {
                            input.read(buffer.data(), length*sizeof(*buffer.data()));
                            output.write(buffer.data(), length*sizeof(*buffer.data()));
                        }
                        input.read(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        output.write(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        current += amount;
                        marker = current*(file+1 < nfiles);
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        output.close();
                        ok = (input.good()) && (output.good());
                    }
                } else {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        output.open(filenames[file], mode | std::ios::binary);
                        amount = (limit > 0) ? ((file+1 < nfiles) ? (static_cast<long long int>(limit)) : (total-current)) : ((total/static_cast<long long int>(nfiles))+(file < total%static_cast<long long int>(nfiles)));
                        marker = current;
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        for (long long int byte = 0; byte < amount; ++byte) {
                            input.read(&contents, sizeof(contents));
                            output.write(&contents, sizeof(contents));
                        }
                        current += amount;
                        marker = current*(file+1 < nfiles);
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        output.close();
                        ok = (input.good()) && (output.good());
                    }
                } 
            }
        }
    }
    return ok*nfiles;
}

// Unsplit several files in a file
/// \brief          Unsplit several files in a file.
/// \details        Unsplits the provided list of files into the original one. 
///                 The original marker size should correspond to the provided 
///                 one. By default, if a standard output open mode is 
///                 specified, all existing files are protected and only a new 
///                 one can be created. This protection can be removed by 
///                 explicitely specifying the truncate open mode. Finally, the 
///                 chunk parameter allows to specify whether a buffer is used 
///                 for the copy : 0 corresponds to no buffer, a positive 
///                 integer corresponds to the amount of bytes that is put into 
///                 the buffer, and a negative integer corresponds to a buffer 
///                 of the total size of the file. 
/// \tparam         Byteswap Swap endianness of markers.
/// \tparam         Marker Record marker type.
/// \tparam         Container (File container type.)
/// \param[in]      filenames Input file names to unsplit.
/// \param[in]      filename Output file name.
/// \param[in]      mode Open mode.
/// \param[in]      chunk Buffer size.
/// \return         Number of files written on success, zero on error.
template <bool Byteswap, typename Marker, class Container> 
unsigned int FileSystem::unsplit(const Container& filenames, const std::string& filename, const std::ios::openmode& mode, const long long int chunk)
{
    std::ifstream input;
    const unsigned int nfiles = filenames.size();
    bool ok = (((mode) & (std::ios::trunc | std::ios::app)) || (!exist(filename))) && (nfiles > 0);
    long long int amount = 0;
    long long int length = 0;
    long long int niterations = 0;
    Marker marker = Marker();
    Marker mark = Marker();
    char contents = char();
    std::vector<char> buffer((chunk >= 0) ? (chunk) : (1));
    std::ofstream output;
    if (ok) {
        for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
            ok = false;
            input.open(filenames[file], std::ios::binary);
            if (input) {
                amount = size(input)-static_cast<long long int>(sizeof(marker)+sizeof(mark));
                if (amount >= 0) {
                    input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                    byteswap<Byteswap>(marker);
                    if (marker == mark) {
                        input.seekg(static_cast<long long int>(sizeof(marker))+amount, std::ios::beg);
                        input.read(reinterpret_cast<char*>(&mark), sizeof(mark));
                        byteswap<Byteswap>(mark);
                        ok = (file+1 < nfiles) ? (mark == marker+amount) : (mark == Marker());
                    }
                }
            }
            input.close();
        }
        if (ok) {
            output.open(filename, mode | std::ios::binary);
            ok = output.good();
            if (ok) {
                if (chunk < 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        input.open(filenames[file], std::ios::binary);
                        amount = size(input)-static_cast<long long int>(sizeof(marker)+sizeof(mark));
                        input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                        if (amount > static_cast<long long int>(buffer.size())) {
                            buffer.resize(amount);
                        }
                        input.read(buffer.data(), amount*sizeof(*buffer.data()));
                        output.write(buffer.data(), amount*sizeof(*buffer.data()));
                        input.close();
                        ok = (input.good()) && (output.good());
                    }
                } else if (chunk > 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        input.open(filenames[file], std::ios::binary);
                        amount = size(input)-static_cast<long long int>(sizeof(marker)+sizeof(mark));
                        length = (amount < chunk) ? (amount) : (chunk);
                        niterations = (length > 0) ? (amount/length) : (0);
                        input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                        for (long long int iteration = 0; iteration < niterations; ++iteration) {
                            input.read(buffer.data(), length*sizeof(*buffer.data()));
                            output.write(buffer.data(), length*sizeof(*buffer.data()));
                        }
                        input.read(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        output.write(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        input.close();
                        ok = (input.good()) && (output.good());
                    }
                } else {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        input.open(filenames[file], std::ios::binary);
                        amount = size(input)-static_cast<long long int>(sizeof(marker)+sizeof(mark));
                        input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                        for (long long int byte = 0; byte < amount; ++byte) {
                            input.read(&contents, sizeof(contents));
                            output.write(&contents, sizeof(contents));
                        }
                        input.close();
                        ok = (input.good()) && (output.good());
                    }
                }
            }
        }
    }
    return ok;
}

// Join several files in a file
/// \brief          Join several files in a file.
/// \details        Joins the provided list of file into a single file. At the 
///                 beginning and the end of each file a record marker of the 
///                 file size in bytes is added. By default, if a standard 
///                 output open mode is specified, all existing files are 
///                 protected and only a new one can be created. This protection 
///                 can be removed by explicitely specifying the truncate open 
///                 mode. Finally, the chunk parameter allows to specify whether 
///                 a buffer is used for the copy : 0 corresponds to no buffer,
///                 a positive integer corresponds to the amount of bytes that 
///                 is put into the buffer, and a negative integer corresponds 
///                 to a buffer of the total size of the file. 
/// \tparam         Byteswap Swap endianness of markers.
/// \tparam         Marker Record marker type.
/// \tparam         Container (File container type.)
/// \param[in]      filenames Input file names to join.
/// \param[in]      filename Output file name.
/// \param[in]      mode Open mode.
/// \param[in]      chunk Buffer size.
/// \return         Number of files written on success, zero on error.
template <bool Byteswap, typename Marker, class Container>
unsigned int FileSystem::join(const Container& filenames, const std::string& filename, const std::ios::openmode& mode, const long long int chunk)
{
    std::ifstream input;
    const unsigned int nfiles = filenames.size();
    bool ok = (((mode) & (std::ios::trunc | std::ios::app)) || (!exist(filename))) && (nfiles > 0);
    long long int amount = 0;
    long long int length = 0;
    long long int niterations = 0;
    Marker marker = Marker();
    char contents = char();
    std::vector<char> buffer((chunk >= 0) ? (chunk) : (1));
    std::ofstream output;
    if (ok) {
        for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
            ok = exist(filenames[file]);
        }
        if (ok) {
            output.open(filename, mode | std::ios::binary);
            ok = output.good();
            if (ok) {
                if (chunk < 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        input.open(filenames[file], std::ios::binary);
                        amount = size(input);
                        marker = amount;
                        if (amount > static_cast<long long int>(buffer.size())) {
                            buffer.resize(amount);
                        }
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        input.read(buffer.data(), amount*sizeof(*buffer.data()));
                        output.write(buffer.data(), amount*sizeof(*buffer.data()));
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        input.close();
                        ok = (input.good()) && (output.good());
                    }
                } else if (chunk > 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        input.open(filenames[file], std::ios::binary);
                        amount = size(input);
                        length = (amount < chunk) ? (amount) : (chunk);
                        niterations = (length > 0) ? (amount/length) : (0);
                        marker = amount;
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        for (long long int iteration = 0; iteration < niterations; ++iteration) {
                            input.read(buffer.data(), length*sizeof(*buffer.data()));
                            output.write(buffer.data(), length*sizeof(*buffer.data()));
                        }
                        input.read(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        output.write(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        input.close();
                        ok = (input.good()) && (output.good());
                    }
                } else {
                    for (unsigned int file = 0;  (file < nfiles) && (ok); ++file) {
                        input.open(filenames[file], std::ios::binary);
                        amount = size(input);
                        marker = amount;
                        byteswap<Byteswap>(marker);
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        for (long long int byte = 0; byte < amount; ++byte) {
                            input.read(&contents, sizeof(contents));
                            output.write(&contents, sizeof(contents));
                        }
                        output.write(reinterpret_cast<char*>(&marker), sizeof(marker));
                        input.close();
                        ok = (input.good()) && (output.good());
                    }
                }
            }
        }
    }
    return ok;
}

// Unjoin a file in several files
/// \brief          Unjoin a file in several files.
/// \details        Unjoins the provided file into the original ones. The 
///                 original marker size should correspond to the provided one.
///                 If output file names are not unique, the previous file is
///                 truncated or not depending on the specified open mode. By 
///                 default, if a standard output open mode is specified, all 
///                 existing files are protected and only new ones can be 
///                 created. This protection can be removed by explicitely 
///                 specifying the truncate open mode. Finally, the chunk 
///                 parameter allows to specify whether a buffer is used for 
///                 the copy : 0 corresponds to no buffer, a positive integer 
///                 corresponds to the amount of bytes that is put into the 
///                 buffer, and a negative integer corresponds to a buffer of 
///                 the total size of the file. 
/// \tparam         Byteswap Swap endianness of markers.
/// \tparam         Marker Record marker type.
/// \tparam         Container (File container type.)
/// \param[in]      filename Input file name to unjoin.
/// \param[in]      filenames Output file names.
/// \param[in]      mode Open mode.
/// \param[in]      chunk Buffer size.
/// \return         Number of files written on success, zero on error.
template <bool Byteswap, typename Marker, class Container> 
unsigned int FileSystem::unjoin(const std::string& filename, const Container& filenames, const std::ios::openmode& mode, const long long int chunk)
{
    std::ifstream input(filename, std::ios::binary);
    unsigned int nfiles = filenames.size();
    bool ok = (input) && (nfiles > 0);
    long long int amount = (ok) ? (size(input)) : (0);
    long long int length = 0;
    long long int niterations = 0;
    Marker marker = Marker();
    Marker mark = Marker();
    char contents = char();
    std::vector<char> buffer((chunk >= 0) ? (chunk) : (1));
    std::ofstream output;
    if (ok) {
        if (!((mode) & (std::ios::trunc | std::ios::app))) {
            for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                ok = !exist(filenames[file]);
            }
        }
        if (ok) {
            nfiles = 0;
            while ((input.tellg() < amount) && (ok)) {
                input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                byteswap<Byteswap>(marker);
                input.seekg(marker, std::ios::cur);
                input.read(reinterpret_cast<char*>(&mark), sizeof(mark));
                byteswap<Byteswap>(mark);
                ok = (input.good()) && (marker >= Marker()) && (marker == mark);
                ++nfiles;
            }
            input.seekg(0, std::ios::beg);
            ok = (ok) && (nfiles > 0) && (nfiles <= filenames.size());
            if (ok) {
                if (chunk < 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        output.open(filenames[file], mode | std::ios::binary);
                        input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                        byteswap<Byteswap>(marker);
                        amount = marker;
                        if (amount > static_cast<long long int>(buffer.size())) {
                            buffer.resize(amount);
                        }
                        input.read(buffer.data(), amount*sizeof(*buffer.data()));
                        output.write(buffer.data(), amount*sizeof(*buffer.data()));
                        input.read(reinterpret_cast<char*>(&mark), sizeof(mark));
                        output.close();
                        ok = (input.good()) && (output.good());
                    }
                } else if (chunk > 0) {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        output.open(filenames[file], mode | std::ios::binary);
                        input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                        byteswap<Byteswap>(marker);
                        amount = marker;
                        length = (amount < chunk) ? (amount) : (chunk);
                        niterations = (length > 0) ? (amount/length) : (0);
                        for (long long int iteration = 0; iteration < niterations; ++iteration) {
                            input.read(buffer.data(), length*sizeof(*buffer.data()));
                            output.write(buffer.data(), length*sizeof(*buffer.data()));
                        }
                        input.read(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        output.write(buffer.data(), (amount%(length+(length == 0)))*sizeof(*buffer.data()));
                        input.read(reinterpret_cast<char*>(&mark), sizeof(mark));
                        output.close();
                        ok = (input.good()) && (output.good());
                    }
                } else {
                    for (unsigned int file = 0; (file < nfiles) && (ok); ++file) {
                        output.open(filenames[file], mode | std::ios::binary);
                        input.read(reinterpret_cast<char*>(&marker), sizeof(marker));
                        byteswap<Byteswap>(marker);
                        amount = marker;
                        for (long long int byte = 0; byte < amount; ++byte) {
                            input.read(&contents, sizeof(contents));
                            output.write(&contents, sizeof(contents));
                        }
                        input.read(reinterpret_cast<char*>(&mark), sizeof(mark));
                        output.close();
                        ok = (input.good()) && (output.good());
                    }
                }
            }
        }
    }
    return ok*nfiles;
}
// -------------------------------------------------------------------------- //



// ------------------------------- COMPARISON ------------------------------- //
// Compare two files passed by names
/// \brief          Compare two files passed by names.
/// \details        Compares two files from their file names and return true if 
///                 both can be tested without errors, have the same size, and 
///                 have the same contents. A file compared to itself returns
///                 true whether it can be opened without errors. Finally, the 
///                 chunk parameter allows to specify whether a buffer is used 
///                 for the comparison : 0 corresponds to no buffer, a positive 
///                 integer corresponds to the amount of bytes that is put into 
///                 the buffer, and a negative integer corresponds to a buffer 
///                 of the total size of the file. 
/// \param[in]      first First file name.
/// \param[in]      second Second file name.
/// \param[in]      chunk Buffer size.
/// \return         True if both files compares equal without errors, false 
///                 otherwise.
bool FileSystem::compare(const std::string& first, const std::string& second, const long long int chunk)
{
    std::ifstream source(first, std::ios::binary);
    std::ifstream destination;
    const bool same = (first == second);
    const long long int amount = (!same) ? (size(source)) : (0);
    bool ok = (source.good()) && (amount >= 0);
    const long long int length = ((ok) && (!same)) ? (((chunk < 0) || (amount < chunk)) ? (amount) : (chunk)) : (0);
    long long int niterations = (length > 0) ? (amount/length) : (0);
    char original = char();
    char contents = char();
    std::vector<char> reference(length);
    std::vector<char> buffer;
    if ((ok) && (!same)) {
        destination.open(second, std::ios::binary);
        if (destination) {
            ok = (size(destination) == amount);
            if (ok) {
                if (length > 0) {
                    buffer.resize(length);
                    for (long long int iteration = 0; (iteration < niterations) && (ok); ++iteration) {
                        source.read(reference.data(), length*sizeof(*buffer.data()));
                        destination.read(buffer.data(), length*sizeof(*buffer.data()));
                        ok = std::equal(buffer.begin(), buffer.end(), reference.begin()); 
                    }
                    if (ok) {
                        niterations = amount%length;
                        source.read(reference.data(), niterations*sizeof(*buffer.data()));
                        destination.read(buffer.data(), niterations*sizeof(*buffer.data()));
                        ok = std::equal(buffer.begin(), buffer.begin()+niterations, reference.begin()); 
                    }
                } else {
                    for (long long int byte = 0; (byte < amount) && (ok); ++byte) {
                        source.read(&original, sizeof(original));
                        destination.read(&contents, sizeof(contents));
                        ok = (contents == original);
                    }
                }
            }
        }
        ok = (ok) && (source.good()) && (destination.good());
    }
    return ok;
}

// Compare two files passed by input stream
/// \brief          Compare two files passed by input stream.
/// \details        Compares two opened input streams and return true if both 
///                 can be tested without errors, have the same size, and have 
///                 the same contents. A stream compared to itself returns
///                 true whether it can be accessed without errors. Finally, the 
///                 chunk parameter allows to specify whether a buffer is used 
///                 for the comparison : 0 corresponds to no buffer, a positive 
///                 integer corresponds to the amount of bytes that is put into 
///                 the buffer, and a negative integer corresponds to a buffer 
///                 of the total size of the file. 
/// \param[in,out]  first First input stream.
/// \param[in,out]  second Second input stream.
/// \param[in]      chunk Buffer size.
/// \return         True if both input streams compares equal without errors, 
///                 false otherwise.
bool FileSystem::compare(std::istream& first, std::istream& second, const long long int chunk)
{
    std::streampos position = first.tellg();
    std::streampos marker = second.tellg();
    const bool same = (&first == &second);
    const long long int amount = size(first);
    bool ok = (first.good()) && (second.good()) && (amount >= 0) && ((!same) ? (size(second) == amount) : (true));
    const long long int length = ((ok) && (!same)) ? (((chunk < 0) || (amount < chunk)) ? (amount) : (chunk)) : (0);
    long long int niterations = (length > 0) ? (amount/length) : (0);
    char original = char();
    char contents = char();
    std::vector<char> reference(length);
    std::vector<char> buffer(length);
    if ((ok) && (!same)) {
        first.seekg(0, std::ios::beg);
        second.seekg(0, std::ios::beg);
        if (length > 0) {
            for (long long int iteration = 0; (iteration < niterations) && (ok); ++iteration) {
                first.read(reference.data(), length*sizeof(*buffer.data()));
                second.read(buffer.data(), length*sizeof(*buffer.data()));
                ok = std::equal(buffer.begin(), buffer.end(), reference.begin()); 
            }
            if (ok) {
                niterations = amount%length;
                first.read(reference.data(), niterations*sizeof(*buffer.data()));
                second.read(buffer.data(), niterations*sizeof(*buffer.data()));
                ok = std::equal(buffer.begin(), buffer.begin()+niterations, reference.begin()); 
            }
        } else {
            for (long long int byte = 0; (byte < amount) && (ok); ++byte) {
                first.read(&original, sizeof(original));
                second.read(&contents, sizeof(contents));
                ok = (contents == original);
            }
        }
        first.seekg(position);
        second.seekg(marker);
        ok = (ok) && (first.good()) && (second.good());
    }
    return ok;
}

// Compare two files passed by output stream
/// \brief          Compare two files passed by output stream.
/// \details        Compares two opened output streams and return true if both 
///                 can be tested without errors and have the same size as no 
///                 further comparison can be done with output streams. Finally, 
///                 the chunk parameter allows to specify whether a buffer is 
///                 used for the comparison : 0 corresponds to no buffer, a 
///                 positive integer corresponds to the amount of bytes that is
///                 put into the buffer, and negative integer corresponds to a 
///                 buffer of the total size of the file. 
/// \param[in,out]  first First output stream.
/// \param[in,out]  second Second output stream.
/// \param[in]      chunk Buffer size.
/// \return         True if both output streams compares equal without errors, 
///                 false otherwise.
bool FileSystem::compare(std::ostream& first, std::ostream& second, const long long int chunk)
{
    const long long int amount = size(first);
    const bool ok = (amount >= 0) && ((&first != &second) ? (size(second) == amount) : (true));
    (void)(chunk);
    return (ok) && (first) && (second);
}
// -------------------------------------------------------------------------- //



// ------------------------ EXISTENCE AND CREATION -------------------------- //
// Check file existence
/// \brief          Check file existence.
/// \details        Returns whether the file can be opened. The difference
///                 with the <tt>check()</tt> function is that <tt>exists()</tt>
///                 use a C test without checking the C++ error bits. 
///                 Consequently, it may be faster than the <tt>check()</tt> 
///                 implementation and therefore is well suited to check the 
///                 existence of a large number of files.
/// \param[in]      filename File name.
/// \return         True if the file exists, false otherwise.
inline bool FileSystem::exist(const std::string& filename)
{
    FILE* file = std::fopen(filename.data(), "rb");
    return (file != nullptr) ? (!std::fclose(file)) : (false);
}

// Check file consistency
/// \brief          Check file consistency.
/// \details        Returns whether the file can be opened successfully and 
///                 whether it contains between the minimum and the maximum
///                 data of the provided type. If the total size is not 
///                 divisible by the size of the type, the function returns 
///                 false. If the file is not empty, reading a byte is also 
///                 tested.  A value of -1 for the minimum or for the 
///                 maximum means that this boundary is not tested.
/// \tparam         Type Type representing the considered unit.
/// \param[in]      filename File name.
/// \param[in]      min Check if the file contains at least this amount of data.
/// \param[in]      max Check if the file contains at most this amount of data.
/// \return         True if the file can be opened an is compliant to the 
///                 provided parameters, false otherwise.
template <typename Type> 
bool FileSystem::check(const std::string& filename, const long long int min, const long long int max)
{
    static const long long int unit = sizeof(Type);
    std::ifstream stream(filename, std::ios::binary);
    long long int length = size(stream);
    bool ok = (length >= 0) && (length%unit == 0);
    char contents = char();
    if (ok) {
        if (length >= unit) {
            stream.read(&contents, sizeof(contents));
        }
        ok = (stream) && ((min < 0) || ((min >= 0) && (length >= min*unit))) && ((max < 0) || ((max >= 0) && (length <= max*unit)));
    }
    return ok;
}

// Create a file
/// \brief          Create a file.
/// \details        Creates a file, overwriting any previous file if needed.
///                 This is equivalent to create a new file using an 
///                 <tt>std::ofstream</tt> and the <tt>std::ios::trunc</tt> open
///                 mode. The extra parameter allows to initialize the file 
///                 with contents : if this parameter is convertible to a 
///                 string, the file is opened as a text file, otherwise it is
///                 considered as a binary file.
/// \tparam         Contents Contents type.
/// \param[in]      filename File name.
/// \param[in]      contents Contents to add to the file.
/// \return         True if the file was created without errors, false 
///                 otherwise.
template <class Contents>
bool FileSystem::create(const std::string& filename, const Contents& contents)
{
    std::ofstream stream(filename, (std::is_convertible<Contents, std::string>::value) ? (std::ios::trunc) : (std::ios::trunc | std::ios::binary));
    return (((!std::is_same<Contents, std::true_type>::value) && (stream)) ? ((std::is_convertible<Contents, std::string>::value) ? (stream<<contents) : (stream.write(reinterpret_cast<const char*>(&contents), sizeof(contents)))) : (stream)).good();
}

// Initialize a new file
/// \brief          Initialize a new file.
/// \details        Creates a new file, without overwriting any previous file.
///                 This function is well suited to create a file without 
///                 risking to erase some important existing data. The extra 
///                 parameter allows to initialize the file with contents : if 
///                 this parameter is convertible to a string, the file is 
///                 opened as a text file, otherwise it is considered as a 
///                 binary file.
/// \tparam         Contents Contents type.
/// \param[in]      filename File name.
/// \param[in]      contents Contents to add to the file.
/// \return         True if the file was created without errors, false 
///                 otherwise.
template <class Contents>
bool FileSystem::initialize(const std::string& filename, const Contents& contents)
{
    return (exist(filename)) ? (false) : (create(filename, contents));
}

// Reset an existing file
/// \brief          Reset an existing file.
/// \details        Erases the contents of an existing file without creating
///                 a new one if the specified name does not exist. This 
///                 function is well suited to avoid the unexpected creation of
///                 new files. The extra parameter allows to initialize the file 
///                 with contents : if this parameter is convertible to a 
///                 string, the file is opened as a text file, otherwise it is
///                 considered as a binary file.
/// \tparam         Contents Contents type.
/// \param[in]      filename File name.
/// \param[in]      contents Contents to add to the file.
/// \return         True if the file was created without errors, false 
///                 otherwise.
template <class Contents>
bool FileSystem::reset(const std::string& filename, const Contents& contents)
{
    return (exist(filename)) ? (create(filename, contents)) : (false);
}

// Generate a binary file based on contents
/// \brief          Generate a binary file based on contents.
/// \details        Creates a new file putting a repetition of the specified 
///                 amount of contents in it. To use the binary mode, this 
///                 contents should not be convertible to a string. By default, 
///                 if a standard output open mode is specified, all existing 
///                 files are protected and only a new one can be created. This 
///                 protection can be removed by explicitely specifying the 
///                 truncate open mode. Finally, the chunk parameter allows to 
///                 specify whether a buffer is used in binary mode : 0 
///                 corresponds to no buffer, a positive integer corresponds to 
///                 the amount of contents that is put into the buffer, and a 
///                 negative integer corresponds to a buffer of the total size 
///                 of the file.
/// \tparam         Contents (Contents type.)
/// \param[in]      filename File name.
/// \param[in]      mode Open mode.
/// \param[in]      amount Amount of contents to be put in the file.
/// \param[in]      chunk Buffer size.
/// \param[in]      contents Contents to add to the file.
/// \return         True if the file was created without errors, false 
///                 otherwise.
template <class Contents, class>
bool FileSystem::generate(const std::string& filename, const std::ios::openmode& mode, const long long int amount, const long long int chunk, const Contents& contents)
{
    bool ok = (((mode) & (std::ios::trunc | std::ios::app)) || (!exist(filename))) && (amount >= 0);
    const long long int length = (ok) ? (((chunk < 0) || (amount < chunk)) ? (amount) : (chunk)) : (0);
    const long long int niterations = (length > 0) ? (amount/length) : (0);
    const std::vector<Contents> buffer(length, contents);
    std::ofstream stream;
    if (ok) {
        stream.open(filename, mode | std::ios::binary);
        if (stream) {
            if (length > 0) {
                for (long long int iteration = 0; iteration < niterations; ++iteration) {
                    stream.write(reinterpret_cast<const char*>(buffer.data()), length*sizeof(*buffer.data()));
                }
                stream.write(reinterpret_cast<const char*>(buffer.data()), (amount%length)*sizeof(*buffer.data()));
            } else {
                for (long long int byte = 0; byte < amount; ++byte) {
                    stream.write(reinterpret_cast<const char*>(&contents), sizeof(contents));
                }
            }
        }
        ok = stream.good();
    }
    return ok;
}

// Generate a text file based on contents
/// \brief          Generate a text file based on contents.
/// \details        Creates a new text file putting a repetition of the 
///                 specified amount of contents in it. To use the text mode, 
///                 this contents should be convertible to a string. By default, 
///                 if a standard output open mode is specified, all existing 
///                 files are protected and only a new one can be created. This 
///                 protection can be removed by explicitely specifying the 
///                 truncate open mode. Finally, the chunk parameter allows to 
///                 specify whether a buffer string is used in text mode : 0 
///                 corresponds to no buffer, a positive integer corresponds to
///                 the amount of contents copies that are put into the buffer, 
///                 and a negative integer corresponds to a buffer of the total
///                 size of the file.
/// \tparam         Contents (Contents type.)
/// \param[in]      filename File name.
/// \param[in]      mode Open mode.
/// \param[in]      amount Amount of contents to be put in the file.
/// \param[in]      chunk Buffer size.
/// \param[in]      contents Contents to add to the file.
/// \return         True if the file was created without errors, false 
///                 otherwise.
template <class Contents, class, class>
bool FileSystem::generate(const std::string& filename, const std::ios::openmode& mode, const long long int amount, const long long int chunk, const Contents& contents)
{
    bool ok = (((mode) & (std::ios::trunc | std::ios::app)) || (!exist(filename))) && (amount >= 0);
    const long long int length = (ok) ? (((chunk < 0) || (amount < chunk)) ? (amount) : (chunk)) : (0);
    long long int niterations = (length > 0) ? (amount/length) : (0);
    std::string original(contents);
    std::string buffer;
    std::ofstream stream;
    if (ok) {
        stream.open(filename, mode);
        if (stream) {
            if (length > 0) {
                buffer.reserve(length*original.size());
                for (long long int line = 0; line < length; ++line) {
                    buffer += original;
                }
                for (long long int iteration = 0; iteration < niterations; ++iteration) {
                    stream<<buffer;
                }
                niterations = amount%length;
                for (long long int iteration = 0; iteration < niterations; ++iteration) {
                    stream<<original;
                }
            } else {
                for (long long int line = 0; line < amount; ++line) {
                    stream<<original;
                }
            }
        }
        ok = stream.good();
    }
    return ok;
}

// Generate a random binary file
/// \brief          Generate a random binary file.
/// \details        Creates a new random file putting a repetition of the 
///                 random numbers generated thanks to the specified engine
///                 and distribution and using the distribution result type
///                 as the contents type. By default, if a standard output open 
///                 mode is specified, all existing files are protected and 
///                 only a new one can be created. This protection can be 
///                 removed by explicitely specifying the truncate open mode. 
///                 Finally, the chunk parameter allows to specify whether a 
///                 buffer is used in binary mode : 0 corresponds to no buffer, 
///                 a positive integer corresponds to the amount of contents 
///                 that is put into the buffer, and a negative integer 
///                 corresponds to a buffer of the total size of the file.
/// \tparam         Engine (Random engine type.)
/// \tparam         Distribution (Random distribution type.)
/// \tparam         Contents (Contents type.)
/// \param[in]      filename File name.
/// \param[in]      mode Open mode.
/// \param[in]      amount Amount of contents to be put in the file.
/// \param[in]      chunk Buffer size.
/// \param[in,out]  engine Random engine.
/// \param[in,out]  distribution Random distribution.
/// \return         True if the file was created without errors, false 
///                 otherwise.
template <class Engine, class Distribution, class Contents, class>
bool FileSystem::generate(const std::string& filename, const std::ios::openmode& mode, const long long int amount, const long long int chunk, Engine&& engine, Distribution&& distribution)
{
    bool ok = (((mode) & (std::ios::trunc | std::ios::app)) || (!exist(filename))) && (amount >= 0);
    const long long int length = (ok) ? (((chunk < 0) || (amount < chunk)) ? (amount) : (chunk)) : (0);
    const long long int niterations = (length > 0) ? (amount/length) : (0);
    Contents contents = Contents();
    std::vector<Contents> buffer(length);
    std::ofstream stream;
    if (ok) {
        stream.open(filename, mode | std::ios::binary);
        if (stream) {
            if (length > 0) {
                std::for_each(buffer.begin(), buffer.end(), [=, &distribution, &engine](Contents& element){element = distribution(engine);});
                for (long long int iteration = 0; iteration < niterations; ++iteration) {
                    stream.write(reinterpret_cast<const char*>(buffer.data()), length*sizeof(*buffer.data()));
                }
                stream.write(reinterpret_cast<const char*>(buffer.data()), (amount%length)*sizeof(*buffer.data()));
            } else {
                for (long long int byte = 0; byte < amount; ++byte) {
                    contents = distribution(engine);
                    stream.write(reinterpret_cast<const char*>(&contents), sizeof(contents));
                }
            }
        }
        ok = stream.good();
    }
    return ok;
}
// -------------------------------------------------------------------------- //



// ------------------------------ SIZE CONTROL ------------------------------ //
// Empty file
/// \brief          Empty file.
/// \details        Returns whether the file exists and is empty. The unit type
///                 is here for compatibility reasons with other size control
///                 functions.
/// \tparam         Type Type representing the considered unit.
/// \tparam         Stream (String, input stream or output stream.)
/// \param[in,out]  stream Stream.
/// \return         True if the file size is null, false otherwise.
template <typename Type, class Stream> 
inline bool FileSystem::empty(Stream& stream)
{
    return (size(stream) == 0*static_cast<long long int>(sizeof(Type)));
}

// File of exact specified size
/// \brief          File of exact specified size.
/// \details        Returns whether the file size measured in the specified unit
///                 is exactly equals to the value. For example if the specified
///                 type is an integer, and the value is equal to 4, it returns
///                 true if the file contains exactly 4 integers. A value equals
///                 to -1 returns true if the file size cannot be computed.
/// \tparam         Type Type representing the considered unit.
/// \tparam         Stream (String, input stream or output stream.)
/// \param[in,out]  stream Stream.
/// \param[in]      amount Amount of contents to for comparison.
/// \return         True if the file size corresponds to the provided value, 
///                 false otherwise.
template <typename Type, class Stream> 
inline bool FileSystem::exact(Stream& stream, const long long int amount)
{
    const long long int length = size(stream);
    return (length >= 0) ? (length == amount*static_cast<long long int>(sizeof(Type))) : (length == amount);
}

// Regular file
/// \brief          Regular file.
/// \details        Tests if the file is regular regarding to the provided 
///                 options. It returns true if the file can be read without
///                 problems, if its size in bytes can be divided by the size
///                 of the provided type, and if its size in the specified unit
///                 is between the minimum and maximum provided amount of data. 
///                 A value of -1 for the minimum or for the maximum means that 
///                 this boundary is not tested.
/// \tparam         Type Type representing the considered unit.
/// \tparam         Stream (String, input stream or output stream.)
/// \param[in,out]  stream Stream.
/// \param[in]      min Check if the file contains at least this amount of data.
/// \param[in]      max Check if the file contains at most this amount of data.
/// \return         True if the file size is regular regarding the parameters,
///                 false otherwise.
template <typename Type, class Stream> 
inline bool FileSystem::regular(Stream& stream, const long long int min, const long long int max)
{
    static const long long int unit = sizeof(Type);
    const long long int length = size(stream);
    return (length >= 0) ? ((length%unit == 0) && ((min < 0) || ((min >= 0) && (length >= min*unit))) && ((max < 0) || ((max >= 0) && (length <= max*unit)))) : (false);
}
// -------------------------------------------------------------------------- //



// ------------------------------- FILE TYPES ------------------------------- //
// Ascii file
/// \brief          Ascii file.
/// \details        Tests whether the file is an ascii file containing between
///                 the minimum and maximum amount of bytes and using the chunk
///                 parameter to control the buffer size. A file is considered
///                 as ascii if all its bytes are in the ranges [9, 13] or [32, 
///                 126]. An empty file is considered as an ascii file. A value 
///                 of -1 for the minimum or for the maximum means that this 
///                 boundary is not tested. Finally, the chunk parameter allows 
///                 to specify whether a buffer is used in binary mode : 0 
///                 corresponds to no buffer, a positive integer corresponds to
///                 the amount of contents that is put into the buffer, and a 
///                 negative integer corresponds to a buffer of the total size 
///                 of the file.
/// \param[in]      filename File name.
/// \param[in]      min Check if the file contains at least this amount of data.
/// \param[in]      max Check if the file contains at most this amount of data.
/// \param[in]      chunk Buffer size.
/// \return         True if the file is an ascii file satisfying the provided 
///                 conditions, false otherwise.
bool FileSystem::ascii(const std::string& filename, const long long int min, const long long int max, const long long int chunk)
{
    static const long long int unit = sizeof(char);
    std::ifstream stream(filename, std::ios::binary);
    const long long int amount = size(stream);
    bool ok = (stream.good()) && ((amount >= 0) ? ((amount%unit == 0) && ((min < 0) || ((min >= 0) && (amount >= min*unit))) && ((max < 0) || ((max >= 0) && (amount <= max*unit)))) : (false));
    const long long int length = (ok) ? (((chunk < 0) || (amount < chunk*unit)) ? (amount) : (chunk*unit)) : (0);
    const long long int niterations = (length > 0) ? (amount/length) : (0);
    char contents = char();
    std::vector<char> buffer(length);
    if (ok) {
        if (length > 0) {
            for (long long int iteration = 0; (iteration < niterations) && (ok); ++iteration) {
                stream.read(buffer.data(), length*sizeof(*buffer.data()));
                ok = std::all_of(buffer.begin(), buffer.end(), [](const char& element){return (((element >= static_cast<char>(9)) && (element <= static_cast<char>(13))) || ((element >= static_cast<char>(32)) && (element <= static_cast<char>(126))));}); 
            }
            if (ok) {
                stream.read(buffer.data(), (amount%length)*sizeof(*buffer.data()));
                ok = std::all_of(buffer.begin(), buffer.begin()+(amount%length), [](const char& element){return (((element >= static_cast<char>(9)) && (element <= static_cast<char>(13))) || ((element >= static_cast<char>(32)) && (element <= static_cast<char>(126))));}); 
            }
        } else {
            for (long long int byte = 0; (byte < amount) && (ok); ++byte) {
                stream.read(&contents, sizeof(contents));
                ok = ((contents >= static_cast<char>(9)) && (contents <= static_cast<char>(13))) || ((contents >= static_cast<char>(32)) && (contents <= static_cast<char>(126)));
            }
        }
        ok = (ok) && (stream.good());
    }
    return ok;
}

// Extended ascii file
/// \brief          Extended ascii file.
/// \details        Tests whether the file is an extended ascii file containing 
///                 between the minimum and maximum amount of bytes and using 
///                 the chunk parameter to control the buffer size. A file is 
///                 considered as extended ascii if all its bytes are in the 
///                 ranges [9, 13] or [32, 126], or [128, 255]. An empty file is
///                 considered as an extended ascii file. A value of -1 for the 
///                 minimum or for the maximum means that this boundary is not 
///                 tested. Finally, the chunk parameter allows to specify 
///                 whether a buffer is used in binary mode : 0 corresponds to 
///                 no buffer, a positive integer corresponds to the amount of 
///                 contents that is put into the buffer, and a negative integer 
///                 corresponds to a buffer of the total size of the file.
/// \param[in]      filename File name.
/// \param[in]      min Check if the file contains at least this amount of data.
/// \param[in]      max Check if the file contains at most this amount of data.
/// \param[in]      chunk Buffer size.
/// \return         True if the file is an extended ascii file satisfying the 
///                 provided conditions, false otherwise.
bool FileSystem::eascii(const std::string& filename, const long long int min, const long long int max, const long long int chunk)
{
    static const long long int unit = sizeof(char);
    std::ifstream stream(filename, std::ios::binary);
    const long long int amount = size(stream);
    bool ok = (stream.good()) && ((amount >= 0) ? ((amount%unit == 0) && ((min < 0) || ((min >= 0) && (amount >= min*unit))) && ((max < 0) || ((max >= 0) && (amount <= max*unit)))) : (false));
    const long long int length = (ok) ? (((chunk < 0) || (amount < chunk*unit)) ? (amount) : (chunk*unit)) : (0);
    const long long int niterations = (length > 0) ? (amount/length) : (0);
    char contents = char();
    std::vector<char> buffer(length);
    if (ok) {
        if (length > 0) {
            for (long long int iteration = 0; (iteration < niterations) && (ok); ++iteration) {
                stream.read(buffer.data(), length*sizeof(*buffer.data()));
                ok = std::all_of(buffer.begin(), buffer.end(), [](const char& element){return ((element < static_cast<char>(0)) || ((element >= static_cast<char>(9)) && (element <= static_cast<char>(13))) || ((element >= static_cast<char>(32)) && (element <= static_cast<char>(126))));}); 
            }
            if (ok) {
                stream.read(buffer.data(), (amount%length)*sizeof(*buffer.data()));
                ok = std::all_of(buffer.begin(), buffer.begin()+(amount%length), [](const char& element){return ((element < static_cast<char>(0)) || ((element >= static_cast<char>(9)) && (element <= static_cast<char>(13))) || ((element >= static_cast<char>(32)) && (element <= static_cast<char>(126))));}); 
            }
        } else {
            for (long long int byte = 0; (byte < amount) && (ok); ++byte) {
                stream.read(&contents, sizeof(contents));
                ok = (contents < static_cast<char>(0)) || ((contents >= static_cast<char>(9)) && (contents <= static_cast<char>(13))) || ((contents >= static_cast<char>(32)) && (contents <= static_cast<char>(126)));
            }
        }
        ok = (ok) && (stream.good());
    }
    return ok;
}

// Binary file
/// \brief          Binary file.
/// \details        Tests whether the file is a binary file containing between 
///                 the minimum and maximum amount of data and using the chunk 
///                 parameter to control the buffer size. A file is considered 
///                 as binary if it has a byte outside of the ranges [9, 13] or 
///                 [32, 126], or [128, 255]. An empty file is not considered as 
///                 a binary file. If the file size is not divisible by the 
///                 specified unit, the function returns false. A value of -1 
///                 for the minimum or for the maximum means that this boundary 
///                 is not tested. Finally, the chunk parameter allows to 
///                 specify whether a buffer is used in binary mode : 0 
///                 corresponds to no buffer, a positive integer corresponds to 
///                 the amount of contents that is put into the buffer, and a 
///                 negative integer corresponds to a buffer of the total size 
///                 of the file.
/// \param[in]      filename File name.
/// \param[in]      min Check if the file contains at least this amount of data.
/// \param[in]      max Check if the file contains at most this amount of data.
/// \param[in]      chunk Buffer size.
/// \return         True if the file is a binary file satisfying the provided 
///                 conditions, false otherwise.
template <typename Type> 
bool FileSystem::binary(const std::string& filename, const long long int min, const long long int max, const long long int chunk)
{
    static const long long int unit = sizeof(Type);
    std::ifstream stream(filename, std::ios::binary);
    const long long int amount = size(stream);
    bool ok = (stream.good()) && ((amount >= 0) ? ((amount%unit == 0) && ((min < 0) || ((min >= 0) && (amount >= min*unit))) && ((max < 0) || ((max >= 0) && (amount <= max*unit)))) : (false));
    const long long int length = (ok) ? (((chunk < 0) || (amount < chunk*unit)) ? (amount) : (chunk*unit)) : (0);
    const long long int niterations = (length > 0) ? (amount/length) : (0);
    char contents = char();
    std::vector<char> buffer(length);
    if (ok) {
        ok = false;
        if (length > 0) {
            for (long long int iteration = 0; (iteration < niterations) && (!ok); ++iteration) {
                stream.read(buffer.data(), length*sizeof(*buffer.data()));
                ok = std::any_of(buffer.begin(), buffer.end(), [](const char& element){return (((element >= static_cast<char>(0)) && (element < static_cast<char>(9))) || ((element > static_cast<char>(13)) && (element < static_cast<char>(32))) || (element == static_cast<char>(127)));}); 
            }
            if (!ok) {
                stream.read(buffer.data(), (amount%length)*sizeof(*buffer.data()));
                ok = std::any_of(buffer.begin(), buffer.begin()+(amount%length), [](const char& element){return (((element >= static_cast<char>(0)) && (element < static_cast<char>(9))) || ((element > static_cast<char>(13)) && (element < static_cast<char>(32))) || (element == static_cast<char>(127)));}); 
            }
        } else {
            for (long long int byte = 0; (byte < amount) && (!ok); ++byte) {
                stream.read(&contents, sizeof(contents));
                ok = ((contents >= static_cast<char>(0)) && (contents < static_cast<char>(9))) || ((contents > static_cast<char>(13)) && (contents < static_cast<char>(32))) || (contents == static_cast<char>(127));
            }
        }
        ok = (ok) && (stream.good());
    }
    return ok;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of FileSystem.
/// \return         0 if no error.
int FileSystem::example()
{
    // Initialize
    std::cout<<"BEGIN = FileSystem::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    int i = 42;
    unsigned int ibom = 0;
    std::uniform_int_distribution<int> d(0, 100);
    std::mt19937 e;
    std::string filename(std::tmpnam(nullptr));
    std::string ifilename(std::tmpnam(nullptr));
    std::string ofilename(std::tmpnam(nullptr));
    std::vector<std::string> filenames = {std::string(std::tmpnam(nullptr)), std::string(std::tmpnam(nullptr)), std::string(std::tmpnam(nullptr)), std::string(std::tmpnam(nullptr))};
    std::ofstream output(filename, std::ios::binary);
    std::ifstream input(filename, std::ios::binary);
    std::ofstream out;
    std::ifstream in;
    output<<"The Ultimate Question of Life, the Universe and Everything"<<std::endl;
    
    // Construction
    FileSystem filesystem;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"FileSystem() : "                              ; FileSystem(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"filesystem = FileSystem() : "                 ; filesystem = FileSystem(); std::cout<<std::endl;
    
    // Utilities
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Utilities : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"filesystem.endianness() : "                   <<filesystem.endianness()<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.byteswap<false>(i) : "             <<filesystem.byteswap<false>(i)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.byteswap<true>(i) : "              <<filesystem.byteswap<true>(i)<<std::endl;
    std::cout<<std::setw(width)<<"ibom = filesystem.bom<unsigned int>() : "     <<(ibom = filesystem.bom<unsigned int>())<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.bom(ibom) : "                      <<filesystem.bom(ibom)<<std::endl;

    // Size
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Size : "                                      <<std::endl;
    std::cout<<std::setw(width)<<"filesystem.weight(filename) : "               <<filesystem.weight(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.weight<int>(filename) : "          <<filesystem.weight<int>(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.size(filename) : "                 <<filesystem.size(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.size(input) : "                    <<filesystem.size(input)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.size(output) : "                   <<filesystem.size(output)<<std::endl;

    // Management
    input.close();
    output.close();
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Management : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.temporary() : "                                                          <<filesystem.temporary()<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.temporary(\"prefix_\", \"_suffix\") : "                                  <<filesystem.temporary("prefix_", "_suffix")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.dated() : "                                                              <<filesystem.dated()<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.dated(\"prefix_\", \"_suffix\", \"%Y\") : "                              <<filesystem.dated("prefix_", "_suffix", "%Y")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.remove(filename+\".txt\") : "                                            <<filesystem.remove(filename+".txt")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.rename(filename, filename+\".txt\") : "                                  <<filesystem.rename(filename, filename+".txt")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.copy(filename+\".txt\", filename) : "                                    <<filesystem.copy(filename+".txt", filename)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.remove(filename+\".txt\") : "                                            <<filesystem.remove(filename+".txt")<<std::endl;
    
    // Split and join
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Split and join : "                                                                  <<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.split(filename, filenames) : "                                           <<filesystem.split(filename, filenames)<<std::endl;std::cout<<std::setw(width*2)<<"filesystem.split(filename, filenames) : "                                           <<filesystem.split(filename, filenames)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.unsplit(filenames, ifilename) : "                                        <<filesystem.unsplit(filenames, ifilename)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.split<true, int>(filename, filenames, std::ios::trunc, 8192, 15) : "     <<filesystem.split<true, int>(filename, filenames, std::ios::trunc, 8192, 15)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.unsplit<true, int>(filenames, ofilename, std::ios::trunc, 8192) : "      <<filesystem.unsplit<true, int>(filenames, ofilename, std::ios::trunc, 8192)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.join(filenames, ofilename, std::ios::trunc) : "                          <<filesystem.join(filenames, ofilename, std::ios::trunc)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.unjoin(ofilename, filenames, std::ios::trunc) : "                        <<filesystem.unjoin(ofilename, filenames, std::ios::trunc)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.join<true, int>(filenames, ofilename, std::ios::trunc) : "               <<filesystem.join<true, int>(filenames, ofilename, std::ios::trunc)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.unjoin<true, int>(ofilename, filenames, std::ios::trunc) : "             <<filesystem.unjoin<true, int>(ofilename, filenames, std::ios::trunc)<<std::endl;

    // Comparison
    input.open(filename, std::ios::binary);
    in.open(ifilename, std::ios::binary);
    output.open(filename, std::ios::app | std::ios::binary);
    out.open(ifilename, std::ios::app | std::ios::binary);
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Comparison : "                                                                      <<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.compare(filename, ifilename) : "                                         <<filesystem.compare(filename, ifilename)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.compare(input, in) : "                                                   <<filesystem.compare(input, in)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.compare(output, out) : "                                                 <<filesystem.compare(output, out)<<std::endl;
    input.close();
    in.close();
    output.close();
    out.close();
    
    // Existence and creation
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Existence and creation : "                                                          <<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.exist(filename) : "                                                      <<filesystem.exist(filename)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.exist(std::tmpnam(nullptr)) : "                                          <<filesystem.exist(std::tmpnam(nullptr))<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.check(filename) : "                                                      <<filesystem.check(filename)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.check(filename, 42, 84) : "                                              <<filesystem.check(filename, 42, 84)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.create(filename, \"The Ultimate Question 1\") : "                        <<filesystem.create(filename, "The Ultimate Question 1")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.initialize(filename, \"The Ultimate Question 2\") : "                    <<filesystem.initialize(filename, "The Ultimate Question 2")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.reset(filename, \"The Ultimate Question 3\") : "                         <<filesystem.reset(filename, "The Ultimate Question 3")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.reset(ofilename) : "                                                     <<filesystem.reset(ofilename)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.reset(ifilename, \"The Ultimate Question\") : "                          <<filesystem.reset(ifilename, "The Ultimate Question")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.generate(filename, std::ios::trunc, 1024) : "                            <<filesystem.generate(filename, std::ios::trunc, 1024)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.generate(filename, std::ios::trunc, 1024, 65536, 42) : "                 <<filesystem.generate(filename, std::ios::trunc, 1024, 65536, 42)<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.generate(filename, std::ios::trunc, 1024, 65536, \"42\") : "             <<filesystem.generate(filename, std::ios::trunc, 1024, 65536, "42")<<std::endl;
    std::cout<<std::setw(width*2)<<"filesystem.generate(filename, std::ios::trunc, 1024, 65536, e, d) : "               <<filesystem.generate(filename, std::ios::trunc, 1024, 65536, e, d)<<std::endl;

    // Size control
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Size control : "                              <<std::endl;
    std::cout<<std::setw(width)<<"filesystem.empty(filename) : "                <<filesystem.empty(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.exact(filename, 42) : "            <<filesystem.exact(filename, 42)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.regular(filename, 42, 4242): "     <<filesystem.regular(filename, 42, 4242)<<std::endl;

    // File types
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"File types : "                                <<std::endl;
    std::cout<<std::setw(width)<<"filesystem.ascii(filename) : "                <<filesystem.ascii(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.ascii(ofilename) : "               <<filesystem.ascii(ofilename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.ascii(ifilename) : "               <<filesystem.ascii(ifilename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.ascii(filename, 0, 4096) : "       <<filesystem.ascii(filename, 0, 4096)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.eascii(filename) : "               <<filesystem.eascii(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.eascii(ofilename) : "              <<filesystem.eascii(ofilename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.eascii(ifilename) : "              <<filesystem.eascii(ifilename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.eascii(filename, 0, 4096) : "      <<filesystem.eascii(filename, 0, 4096)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.binary(filename) : "               <<filesystem.binary(filename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.binary(ofilename) : "              <<filesystem.binary(ofilename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.binary(ifilename) : "              <<filesystem.binary(ifilename)<<std::endl;
    std::cout<<std::setw(width)<<"filesystem.binary(filename, 0, 4096) : "      <<filesystem.binary(filename, 0, 4096)<<std::endl;
        
    // Finalize
    std::remove(filename.data());
    std::remove(ifilename.data());
    std::remove(ofilename.data());
    std::for_each(filenames.begin(), filenames.end(), [](const std::string& s){std::remove(s.data());});
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = FileSystem::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // FILESYSTEM_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
