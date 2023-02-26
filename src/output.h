/* ********************************* OUTPUT ********************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Output
// DESCRIPTION :    Output utilities for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           output.h
/// \brief          Output utilities for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>
// Include libs
// Include project
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "utility.h"
// Octree
#ifdef VELOCITYFIELD
#include "gravity2.h" // with velocity slots
#else
#include "gravity.h"
#endif
#include "photon.h"
// Misc
// -------------------------------------------------------------------------- //

// ---------------------------------- CLASS --------------------------------- //
// Ouput utilities for raytracing
/// \brief          Output utilities for raytracing.
/// \details        Provides a list of exportation routines to save data for
///                 raytracing.
class Output final {
    // Utilities
    /// \name           Utilities
    //@{
public:
    template <class Type = std::string, class = typename std::enable_if<std::is_convertible<typename std::remove_cv<typename std::remove_reference<Type>::type>::type, std::string>::value>::type>
    static inline std::string name(Type &&value = Type());
    template <class Type, class = typename std::enable_if<!std::is_convertible<typename std::remove_cv<typename std::remove_reference<Type>::type>::type, std::string>::value>::type>
    static inline std::string name(const Type &value);
    template <template <class, class> class Type, class First, class Second, class = typename std::enable_if<(std::tuple_size<Type<First, Second>>::value == std::tuple_size<std::pair<First, Second>>::value) && (std::is_convertible<First, std::string>::value)>::type>
    static inline std::string name(const Type<First, Second> &value);
    template <class Type, class... Types, class = typename std::enable_if<sizeof...(Types) != 0>::type>
    static inline std::string name(Type &&value, Types &&...values);
    //@}

    // Save
    /// \name           Save
    //@{
public:
    template <class Octree, class = typename std::enable_if<Octree::dimension() != 0>::type>
    static bool save(std::ostream &stream, const Octree &octree, const int digits = 0);
    template <class Cosmology, class = typename std::enable_if<std::tuple_size<Cosmology>::value != 0>::type>
    static bool save(std::ostream &stream, const Cosmology &cosmology, const unsigned int digits = 0);
    template <class Trajectory, class = typename std::enable_if<!std::is_void<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Trajectory>()[0].type())>::type>::type>::value>::type>
    static bool save(std::ostream &stream, const Trajectory &trajectory, const unsigned int &digits = 0);
    template <class Container, typename Type = typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type, typename Integral = std::true_type, class = typename std::enable_if<!std::is_void<typename std::remove_cv<typename std::remove_reference<decltype(std::declval<Container>()[0])>::type>::type>::value>::type>
    static bool save(std::ostream &stream, const Container &x, const Container &y, const Container &ymean, const Container &ystd, const unsigned int digits = 0, const Integral count = Integral());
    //@}

    // Test
    /// \name           Test
    //@{
public:
    static int example();
    //@}
};
// -------------------------------------------------------------------------- //

// ------------------------------- UTILITIES -------------------------------- //
// File name from a string
/// \brief          File name from a string.
/// \details        Produces a file name from a string.
/// \tparam         Type Type convertible to a string.
/// \param[in]      value Value of the string.
/// \return         Name corresponding to the string.
template <class Type, class>
inline std::string Output::name(Type &&value) {
    return std::forward<Type>(value);
}

// File name from a number
/// \brief          File name from a number.
/// \details        Produces a file name from a number.
/// \tparam         Type Number type convertible to a string.
/// \param[in]      value Value of the number.
/// \return         Name corresponding to the number.
template <class Type, class>
inline std::string Output::name(const Type &value) {
    return std::to_string(value);
}

// File name from a format
/// \brief          File name from a format.
/// \details        Produces a file name from a format.
/// \tparam         Type Pair type.
/// \tparam         First Type convertible to a string.
/// \tparam         Second Second type associated to the format.
/// \param[in]      value Value of the pair.
/// \return         Name corresponding to the format.
template <template <class, class> class Type, class First, class Second, class>
inline std::string Output::name(const Type<First, Second> &value) {
    std::string input(std::get<0>(value));
    std::vector<char> output(input.size() + std::numeric_limits<unsigned char>::max(), static_cast<char>(0));
    std::snprintf(output.data(), output.size(), input.data(), std::get<1>(value));
    return output.data();
}

// File name from a serie of components
/// \brief          File name from a format.
/// \details        Produces a file name from a serie of components by
///                 recursively calling the name function.
/// \tparam         Type First type.
/// \tparam         Types Other types.
/// \param[in]      value First value.
/// \param[in]      values Other values.
/// \return         Name corresponding to the serie of components.
template <class Type, class... Types, class>
inline std::string Output::name(Type &&value, Types &&...values) {
    return name(value) + name(values...);
}
// -------------------------------------------------------------------------- //

// ---------------------------------- SAVE ---------------------------------- //
// Save an octree
/// \brief          Save a octree.
/// \details        Writes an octree to a text file.
/// \tparam         Octree Octree type.
/// \param[in,out]  stream Output stream.
/// \param[in]      octree Octree.
/// \param[in]      digits Optional precision.
///  \return         True on success, false otherwise.
template <class Octree, class>
bool Output::save(std::ostream &stream, const Octree &octree, const int digits) {
    if (stream) {
        stream << std::setprecision((digits > 0) ? (digits) : (stream.precision()));
        stream << octree;
    }
    return stream.good();
}

// Save a cosmology
/// \brief          Save a cosmology.
/// \details        Writes each step of a cosmology to a text file.
/// \tparam         Cosmology Cosmology type.
/// \param[in,out]  stream Output stream.
/// \param[in]      cosmology Cosmology.
/// \param[in]      digits Optional precision.
///  \return         True on success, false otherwise.
template <class Cosmology, class>
bool Output::save(std::ostream &stream, const Cosmology &cosmology, const unsigned int digits) {
    const unsigned int size = std::get<0>(cosmology).size();
    const char separator = stream.fill();
    if (stream) {
        stream << std::setprecision((digits > 0) ? (digits) : (stream.precision()));
        for (unsigned int i = 0; i < size; ++i) {
            stream << std::get<0>(cosmology)[i] << separator << std::get<1>(cosmology)[i] << separator << std::get<2>(cosmology)[i] << separator << std::get<3>(cosmology)[i] << std::endl;
        }
    }
    return stream.good();
}

// Save a trajectory
/// \brief          Save a trajectory.
/// \details        Writes each step of the trajectory to a text file.
/// \tparam         Trajectory Trajectory type.
/// \param[in,out]  stream Output stream.
/// \param[in]      trajectory Trajectory.
/// \param[in]      digits Optional precision.
///  \return         True on success, false otherwise.
template <class Trajectory, class>
bool Output::save(std::ostream &stream, const Trajectory &trajectory, const unsigned int &digits) {
    const unsigned int size = trajectory.size();
    if (stream) {
        stream << std::setprecision((digits > 0) ? (digits) : (stream.precision()));
        for (unsigned int i = 0; i < size; ++i) {
            stream << trajectory[i] << std::endl;
        }
    }
    return stream.good();
}

// Save statistics
/// \brief          Save statistics.
/// \details        Save photons statistics.
/// \tparam         Container Container type.
/// \tparam         Type Data type.
/// \tparam         Integral Integral type.
/// \param[in,out]  stream Output stream.
/// \param[in]      x Abscissae.
/// \param[in]      y Homogeneous ordinates.
/// \param[in]      ymean Mean of values.
/// \param[in]      ystd Standard deviation of values.
/// \param[in]      digits Optional precision.
/// \param[in]      count Optional count.
///  \return         True on success, false otherwise.
template <class Container, typename Type, typename Integral, class>
bool Output::save(std::ostream &stream, const Container &x, const Container &y, const Container &ymean, const Container &ystd, const unsigned int digits, const Integral count) {
    const unsigned int n = count;
    const unsigned int size = x.size();
    if (stream) {
        stream << std::setprecision((digits > 0) ? (digits) : (stream.precision()));
        if (!std::is_same<Integral, std::true_type>::value) {
            stream << "# " << n << std::endl;
        }
        for (unsigned int i = 0; i < size; ++i) {
            stream << x[i] << " " << y[i] << " " << ymean[i] << " " << ystd[i] << std::endl;
        }
    }
    return stream.good();
}
// -------------------------------------------------------------------------- //

// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Output.
/// \return         0 if no error.
int Output::example() {
    // Initialize
    std::cout << "BEGIN = Input::example()" << std::endl;
    std::cout << std::boolalpha << std::left;
    const unsigned int width = 40;
    magrathea::SimpleHyperOctree<double, magrathea::SimpleHyperOctreeIndex<unsigned long long int, 3>, unsigned int> octree(0, 1);
    std::array<std::vector<double>, 4> cosmology;
    std::vector<Photon<double, 3>> trajectory;
    std::ostringstream oss;

    // Construction
    Output output;

    // Lifecycle and operators
    std::cout << std::endl;
    std::cout << std::setw(width) << "Lifecycle and operators : " << std::endl;
    std::cout << std::setw(width) << "Output() : ";
    Output();
    std::cout << std::endl;
    std::cout << std::setw(width) << "output = Output() : ";
    output = Output();
    std::cout << std::endl;

    // Utilities
    std::cout << std::endl;
    std::cout << std::setw(width * 2) << "Utilities : " << std::endl;
    std::cout << std::setw(width * 2) << "output.name(\"filename_00000.txt\") : " << output.name("filename_00000.txt") << std::endl;
    std::cout << std::setw(width * 2) << "output.name(42) : " << output.name(42) << std::endl;
    std::cout << std::setw(width * 2) << "output.name(std::make_pair(\"filename_%05d.txt\", 42)) : " << output.name(std::make_pair("filename_%05d.txt", 42)) << std::endl;
    std::cout << std::setw(width * 2) << "output.name(\"filename_\", std::make_pair(\"%05d.txt\", 42), \".txt\") : " << output.name("filename_", std::make_pair("%05d.txt", 42), ".txt") << std::endl;

    // Save
    std::cout << std::endl;
    std::cout << std::setw(width) << "Save : " << std::endl;
    std::cout << std::setw(width) << "output.save(oss, octree) : " << output.save(oss, octree) << std::endl;
    std::cout << std::setw(width) << "output.save(oss, cosmology) : " << output.save(oss, cosmology) << std::endl;
    std::cout << std::setw(width) << "output.save(oss, trajectory) : " << output.save(oss, trajectory) << std::endl;
    std::cout << std::setw(width) << "output.save(oss, cosmology[0], cosmology[1], cosmology[2], cosmology[3]) : " << output.save(oss, cosmology[0], cosmology[1], cosmology[2], cosmology[3]) << std::endl;

    // Finalize
    std::cout << std::noboolalpha << std::right << std::endl;
    std::cout << "END = Output::example()" << std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
#endif // OUTPUT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
