/* ********************************* INPUT ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Input
// DESCRIPTION :    Input utilities for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
// (michel-andres.breton@obspm.fr) CONTRIBUTIONS :  [Vincent Reverdy
// (2012-2013), Michel-Andrès Breton (2015-2021)] LICENSE :        CECILL-B
// License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           input.h
/// \brief          Input utilities for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton
/// (michel-andres.breton@obspm.fr) \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <type_traits>
#include <vector>
// Include libs
// Include project
#include "cone.h"
#include "magrathea/constants.h"
#include "magrathea/datahandler.h"
#include "magrathea/filelist.h"
#include "magrathea/filesystem.h"
#include "magrathea/hypercube.h"
#include "magrathea/hypersphere.h"
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
// Include C++ HDF5 project
// #include "ReadHDF5.h"
#include "TReadHDF5.h"
// unordered_set
#include <unordered_set>

// Misc
// -------------------------------------------------------------------------- //

// ---------------------------------- CLASS --------------------------------- //
// Input utilities for raytracing
/// \brief          Input utilities for raytracing.
/// \details        Provides a list of importation routines to load data for
///                 raytracing.
class Input final {
    // Utilities
    /// \name           Utilities
    //@{
public:
    static inline std::string trim(const std::string &text,
                                   const std::string &comment = "#");
    static inline std::pair<std::string, std::string>
    partition(const std::string &text, const std::string &separator = "=");
    template <class Octree, class Source,
              class Data = typename std::tuple_element<
                  1, decltype(Octree::element())>::type,
              class Element = decltype(Source::element()),
              class = typename std::enable_if<
                  (Octree::dimension() == Source::dimension()) &&
                  (std::is_integral<Data>::value)>::type>
    static inline unsigned int count(Octree &octree, const Source &source);
    template <class Octree, class Sphere, class Conic,
              typename Type = decltype(Octree::type()),
              class Element = decltype(Octree::element()),
              class Index = typename std::tuple_element<0, Element>::type,
              unsigned int Dimension = Octree::dimension(),
              class Position = decltype(Octree::position()),
              class Extent = decltype(Octree::extent()),
              class = typename std::enable_if<
                  (Dimension == 3) && (Dimension == Sphere::dimension())>::type>
    static inline bool collide(const Octree &octree, const Index &index,
                               const Sphere &sphere, const Conic &conic);
    template <unsigned int Selection = 0, class Octree,
              unsigned int Dimension = Octree::dimension(),
              class Element = decltype(Octree::element()),
              class Index = typename std::tuple_element<0, Element>::type,
              class Data = typename std::tuple_element<1, Element>::type,
              class Type = decltype(Data::template type<Selection>()),
              class = typename std::enable_if<(Dimension == 3)>::type>
    static inline Type mean(const Octree &octree, const Element &element,
                            int level = -1);
    template <unsigned int Selection = 0, class Octree,
              unsigned int Dimension = Octree::dimension(),
              class Element = decltype(Octree::element()),
              class Index = typename std::tuple_element<0, Element>::type,
              class Data = typename std::tuple_element<1, Element>::type,
              class Type = decltype(Data::template type<Selection>()),
              class = typename std::enable_if<(Dimension == 3)>::type>
    static inline Data meanAll(const Octree &octree, const Element &element,
                               int level = -1);
    template <typename Type, class Cosmology = std::array<std::vector<double>, 4>,
              class = typename std::enable_if<std::is_convertible<
                  Type, typename std::remove_cv<typename std::remove_reference<
                            decltype(std::declval<Cosmology>()[0][0])>::type>::
                            type>::value>::type>
    static inline Cosmology constantify(const unsigned int size, const Type tmin,
                                        const Type tmax, const Type a = Type(1),
                                        const Type dadt = Type(),
                                        const Type d2adt2 = Type());
    template <class Parameter, class Data, typename Type,
              class = typename std::enable_if<Data::types() != 0>::type>
    static inline Data sistemize(const Parameter &parameters, const Data &data,
                                 const Type a, const Type h, const Type omegam,
                                 const Type lboxmpch);
    template <class Parameter, class Octree, typename Type,
              class Element = decltype(Octree::element()),
              class Data = typename std::tuple_element<1, Element>::type,
              class = typename std::enable_if<(!std::is_void<Data>::value) &&
                                              (Octree::dimension() != 0)>::type>
    static inline unsigned int sistemize(const Parameter &parameters,
                                         Octree &octree, const Type h,
                                         const Type omegam, const Type lboxmpch);
    template <class Data,
              class = typename std::enable_if<Data::types() != 0>::type>
    static inline Data homogenize(const Data &data);
    template <class Octree, class... Dummy,
              class Element = decltype(Octree::element()),
              class Data = typename std::tuple_element<1, Element>::type,
              class = typename std::enable_if<(!std::is_void<Data>::value) &&
                                              (Octree::dimension() != 0) &&
                                              (sizeof...(Dummy) == 0)>::type>
    static inline unsigned int homogenize(Octree &octree, Dummy...);
    template <
        class Extent = std::ratio<1>, class Data, class Vector, typename Type,
        unsigned int Dimension = 3,
        class = typename std::enable_if<
            (Data::types() != 0) &&
            (std::is_convertible<
                Type, typename std::remove_cv<typename std::remove_reference<
                          decltype(std::declval<Vector>()[0])>::type>::type>::
                 value)>::type>
    static inline Data schwarzschildify(const Data &data, const Vector &center,
                                        const Vector &position, const Type mass,
                                        const Type length);
    template <class Octree, class Vector, typename Type, class Function,
              class... Dummy, class Element = decltype(Octree::element()),
              unsigned int Dimension = Octree::dimension(),
              class Extent = decltype(Octree::extent()),
              class = typename std::enable_if<(sizeof...(Dummy) == 0) &&
                                              (Dimension == 3)>::type>
    static inline unsigned int
    schwarzschildify(Octree &octree, const Vector &position, const Type mass,
                     const Type length, Function &&refiner, Dummy...);
    //@}

    // Files
    /// \name           Files
    //@{
public:
    template <class Octree, class Element = decltype(Octree::element()),
              class Index = typename std::tuple_element<0, Element>::type,
              class Data = typename std::tuple_element<1, Element>::type,
              unsigned int Dimension = Octree::dimension(),
              class = typename std::enable_if<Dimension != 0>::type>
    static unsigned int filetree(Octree &octree, const std::string &directory,
                                 const std::string &format);
    template <class List, class Octree, class Sphere, class Conic,
              unsigned int Dimension = Octree::dimension(),
              class = typename std::enable_if<
                  (Dimension == 3) && (Dimension == Sphere::dimension())>::type>
    static bool prepare(List &list, const Octree &octree, const Sphere &sphere,
                        const Conic &conic);
    //@}

    // Data
    /// \name           Data
    //@{
public:
    template <
        typename Integral = unsigned int, typename Real = float, class Parameter,
        class Octree, class Function, class Element = decltype(Octree::element()),
        class Index = typename std::tuple_element<0, Element>::type,
        class Data = typename std::tuple_element<1, Element>::type,
        unsigned int Dimension = Octree::dimension(),
        class = typename std::enable_if<std::is_convertible<
            typename std::result_of<Function(Element)>::type, bool>::value>::type>
    static bool import(const Parameter &parameters, Octree &octree,
                       const std::string &filename, Function &&filter);

    template <
        typename Type, class Parameter,
        class Cosmology = std::array<std::vector<Type>, 4>,
        class Element = typename std::remove_cv<typename std::remove_reference<
            decltype(std::declval<Cosmology>()[0])>::type>::type,
        class = typename std::enable_if<std::is_convertible<
            Element, typename std::remove_cv<typename std::remove_reference<
                         decltype(std::declval<Cosmology>()[0])>::type>::type>::
                                            value>::type>
    static Cosmology acquire(const Parameter &parameters, Type &h, Type &omegam,
                             Type &lboxmpch,
                             const std::string &outfile = std::string());

    template <
        class Container = std::map<std::string, std::string>,
        class Element = std::pair<std::string, std::string>,
        class = typename std::enable_if<
            (std::is_convertible<Container,
                                 std::map<std::string, std::string>>::value) &&
            (std::is_convertible<
                Element, std::pair<std::string, std::string>>::value)>::type>
    static Container parse(const std::string &filename,
                           const std::string &separator = "=",
                           const std::string &comment = "#");

    template <
        class Position, class Extent, typename Integral = unsigned int,
        typename Real = float, class Parameter, class Conic, class Octree,
        class Function, class Element = decltype(Octree::element()),
        class Index = typename std::tuple_element<0, Element>::type,
        class Data = typename std::tuple_element<1, Element>::type,
        unsigned int Dimension = Octree::dimension(),
        class = typename std::enable_if<std::is_convertible<
            typename std::result_of<Function(Element)>::type, bool>::value>::type>
    static bool importhdf5(const Parameter &parameters, unsigned int &rank,
                           const std::array<std::array<double, 3>, 3> &rotm1,
                           const double &thetay, const double &thetaz,
                           const Conic &conic, Octree &octree,
                           const std::string &filename, Function &&filter);

    template <
        class Position, class Extent, typename Integral = unsigned int,
        typename Real = float, class Parameter, class Octree, class Function,
        class Element = decltype(Octree::element()),
        class Index = typename std::tuple_element<0, Element>::type,
        class Data = typename std::tuple_element<1, Element>::type,
        unsigned int Dimension = Octree::dimension(),
        class = typename std::enable_if<std::is_convertible<
            typename std::result_of<Function(Element)>::type, bool>::value>::type>
    static bool importfullhdf5(const Parameter &parameters, Octree &octree,
                               const std::string &filename, Function &&filter);

    template <
        class Position, class Extent, typename Integral = unsigned int,
        typename Real = float, class Parameter, class Octree, class Function,
        class Element = decltype(Octree::element()),
        class Index = typename std::tuple_element<0, Element>::type,
        class Data = typename std::tuple_element<1, Element>::type,
        unsigned int Dimension = Octree::dimension(),
        class = typename std::enable_if<std::is_convertible<
            typename std::result_of<Function(Element)>::type, bool>::value>::type>
    static bool importascii(const Parameter &parameters,
                            const std::array<std::array<double, 3>, 3> &rotm1,
                            Octree &octree, const std::string &filename,
                            Function &&filter);
    //@}

    // Cones
    /// \name           Cones
    //@{
public:
    template <class Octree,
              class = typename std::enable_if<Octree::dimension() != 0>::type>
    static bool save(Octree &octree, const std::string &filename);
    template <class Octree,
              class = typename std::enable_if<Octree::dimension() != 0>::type>
    static bool load(Octree &octree, const std::string &filename);
    //@}

    // Correction
    /// \name           Correction
    //@{
public:
    template <class Cosmology, class Trajectory,
              class Type = typename std::remove_cv<typename std::remove_reference<
                  decltype(std::declval<typename std::tuple_element<
                               0, Cosmology>::type>()[0])>::type>::type,
              class = typename std::enable_if<std::is_convertible<
                  Type, typename std::remove_cv<typename std::remove_reference<
                            decltype(std::declval<typename std::tuple_element<
                                         0, Cosmology>::type>()[0])>::type>::
                            type>::value>::type>
    static Cosmology correct(const Cosmology &cosmology,
                             const Trajectory &trajectory);
    template <class Parameter, int Check = 0,
              unsigned int Selection = Check *(Check >= 0), class Octree,
              typename Kind = double,
              unsigned int Dimension = Octree::dimension(),
              class Element = decltype(Octree::element()),
              class Index = typename std::tuple_element<0, Element>::type,
              class Data = typename std::tuple_element<1, Element>::type,
              class Type = decltype(Data::template type<Selection>()),
              class Position = decltype(Octree::position()),
              class Extent = decltype(Octree::extent()),
              class = typename std::enable_if<(Dimension == 3)>::type>
    static Octree &correct(const Parameter &parameters, Octree &octree,
                           Kind &&amin = Kind());
    //@}

    // Test
    /// \name           Test
    //@{
    // public:
    //    static int example();
    //@}
};
// -------------------------------------------------------------------------- //

// ------------------------------- UTILITIES -------------------------------- //
// Trim string
/// \brief          Trim string.
/// \details        Trims string from leading and trailing white spaces, from
///                 comments, and from inner multi white spaces.
/// \param[in]      text Input text.
/// \param[in]      comment Comment string.
/// \return         Trimmed string.
inline std::string Input::trim(const std::string &text,
                               const std::string &comment) {
    std::string result = text;
    result.erase((result.find(comment) == std::string::npos)
                     ? (result.end())
                     : (result.begin() + result.find(comment)),
                 result.end());
    result.erase(result.begin(),
                 std::find_if(result.begin(), result.end(),
                              [](const char c) { return std::isgraph(c); }));
    result.erase(result.end() -
                     std::distance(result.rbegin(),
                                   std::find_if(result.rbegin(), result.rend(),
                                                [](const char c) {
                                                    return std::isgraph(c);
                                                })),
                 result.end());
    result.erase(std::unique(result.begin(), result.end(),
                             [](const char x, const char y) {
                                 return ((std::isspace(x)) && (std::isspace(y)));
                             }),
                 result.end());
    return result;
}

// Partition string
/// \brief          Partition string.
/// \details        Splits the string in two parts before and after the
///                 provided separator.
/// \param[in]      text Input text.
/// \param[in]      separator Separator string.
/// \return         Partitioned string.
inline std::pair<std::string, std::string>
Input::partition(const std::string &text, const std::string &separator) {
    return (text.find(separator) == std::string::npos)
               ? (std::make_pair(text, std::string()))
               : (std::make_pair(std::string(text.begin(),
                                             text.begin() + text.find(separator)),
                                 std::string(text.begin() + text.find(separator) +
                                                 separator.size(),
                                             text.end())));
}

// Count tree
/// \brief          Count tree.
/// \details        Counts the number of input cells in each output cells.
/// \tparam         Octree Output octree type.
/// \tparam         Source Input octree type.
/// \tparam         Data Output data type.
/// \tparam         Element Underlying input element type.
/// \param[in]      octree Output octree.
/// \param[in]      source Input octree.
/// \return         Octree with count of input cells per output cell.
template <class Octree, class Source, class Data, class Element, class>
inline unsigned int Input::count(Octree &octree, const Source &source) {
    std::for_each(source.begin(), source.end(),
                  [=, &octree](const Element &elem) {
                      ++std::get<1>(*(octree.find(std::get<0>(elem))));
                  });
    return octree.size();
}

// Collision between an octree index and a sphere or a cone
/// \brief          Collision between an octree index and a sphere or a cone.
/// \details        Detects collision between an index of an octree and a
///                 sphere or a cone.
/// \tparam         Octree Octree type.
/// \tparam         Sphere Sphere type.
/// \tparam         Conic Cone type.
/// \tparam         Type Scalar position type.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \param[in]      octree Input octree.
/// \param[in]      index Index of one element.
/// \param[in]      sphere Geometrical sphere.
/// \param[in]      cone Three dimensional cone.
/// \return         True if collision, false otherwise.
template <class Octree, class Sphere, class Conic, typename Type, class Element,
          class Index, unsigned int Dimension, class Position, class Extent,
          class>
inline bool Input::collide(const Octree &octree, const Index &index,
                           const Sphere &sphere, const Conic &conic) {
    Sphere spherified = Sphere();
    for (unsigned int idim = 0; idim < Dimension; ++idim) {
        spherified.position(idim) =
            index.template position<Type, Position, Extent>(idim);
    }

    spherified.extent() = index.template extent<Type, Position, Extent>() *
                          std::sqrt(Type(3)) /
                          Type(2 * (sizeof(octree) / sizeof(octree)));

    return ((Utility::collide(spherified, sphere)) ||
            (Utility::collide(spherified, conic, sphere.extent())));
}

// Mean value over cells
/// \brief          Mean value over cells.
/// \details        Computes the average of the provided data in all
///                 surrounding cells which are normal.
/// \tparam         Selection Index of the data to be computed.
/// \tparam         Octree Octree type.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Type Selected data type.
/// \param[in]      octree Input octree.
/// \param[in]      element Input element.
/// \param[in]      level Cell level for computation.
/// \return         Averaged value.
template <unsigned int Selection, class Octree, unsigned int Dimension,
          class Element, class Index, class Data, class Type, class>
inline Type Input::mean(const Octree &octree, const Element &element,
                        int level) {
    static const int zero = 0;
    static const int one = 1;
    static const int two = 2;
    const unsigned int size = octree.size();
    const Index index =
        ((std::get<0>(element).level() > zero) &&
         (static_cast<int>(std::get<0>(element).level()) > level))
            ? (std::get<0>(element).previous(
                  zero, (level < zero) ? (std::get<0>(element).level() - one)
                                       : (level)))
            : (std::get<0>(element));
    const unsigned int lvl = index.level();
    const unsigned long long int nx = index.coordinate(zero);
    const unsigned long long int ny = index.coordinate(one);
    const unsigned long long int nz = index.coordinate(two);
    unsigned int n = zero;
    Type result = zero;
    unsigned int count = zero;
    // Loop over the neighbour cells
    for (int ix = -one; ix <= one; ++ix) {
        for (int iy = -one; iy <= one; ++iy) {
            for (int iz = -one; iz <= one; ++iz) {
                if (!((ix == zero) && (iy == zero) && (iz == zero))) {
                    // Find the cell in the octree
                    n = std::distance(
                        octree.begin(),
                        octree.find(Index::make(lvl, nx + ix, ny + iy, nz + iz)));
                    if (n < size) {
                        // Check if value of neighbour is normal
                        if (std::isnormal(
                                std::get<1>(octree[n]).template data<Selection>())) {
                            result += std::get<1>(octree[n]).template data<Selection>();
                            ++count;
                        }
                    }
                }
            }
        }
    }
    return result / (count + (count == zero));
}

// Mean value over cells
/// \brief          Mean value over cells.
/// \details        Computes the average of cell tuple information all
///                 surrounding cells which are normal.
/// \tparam         Selection Index of the data to be computed.
/// \tparam         Octree Octree type.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Type Selected data type.
/// \param[in]      octree Input octree.
/// \param[in]      element Input element.
/// \param[in]      level Cell level for computation.
/// \return         Averaged value of tuple.
template <unsigned int Selection, class Octree, unsigned int Dimension,
          class Element, class Index, class Data, class Type, class>
inline Data Input::meanAll(const Octree &octree, const Element &element,
                           int level) {
    static const int zero = 0;
    static const int one = 1;
    static const int two = 2;
    const unsigned int size = octree.size();
    const Index index =
        ((std::get<0>(element).level() > zero) &&
         (static_cast<int>(std::get<0>(element).level()) > level))
            ? (std::get<0>(element).previous(
                  zero, (level < zero) ? (std::get<0>(element).level() - one)
                                       : (level)))
            : (std::get<0>(element));
    const unsigned int lvl = index.level();
    const unsigned long long int nx = index.coordinate(zero);
    const unsigned long long int ny = index.coordinate(one);
    const unsigned long long int nz = index.coordinate(two);
    unsigned int n = zero;
    Data result = Data();
    unsigned int count = zero;
    // Loop over the neighbour cells
    for (int ix = -one; ix <= one; ++ix) {
        for (int iy = -one; iy <= one; ++iy) {
            for (int iz = -one; iz <= one; ++iz) {
                if (!((ix == zero) && (iy == zero) && (iz == zero))) {
                    // Find the cell in the octree
                    n = std::distance(
                        octree.begin(),
                        octree.find(Index::make(lvl, nx + ix, ny + iy, nz + iz)));
                    if (n < size) {
                        // Check if value of neighbour is normal
                        if (std::isnormal(
                                std::get<1>(octree[n]).template data<Selection>())) {
                            // Sum tuple of values
                            octree.mac(result, std::get<1>(octree[n]), 1);
                            ++count;
                        }
                    }
                }
            }
        }
    }
    octree.mac(result, result, 1. / (count + (count == zero)) - 1);
    return result;
}

// Constant cosmology
/// \brief          Constant cosmology.
/// \details        Produces a constant cosmology of the provided size.
/// \tparam         Cosmology Cosmology container type.
/// \tparam         Type Data type.
/// \param[in]      size Size of the cosmology tables.
/// \param[in]      tmin First value of time.
/// \param[in]      tmax Last value of time.
/// \param[in]      a Scale factor.
/// \param[in]      dadt Value of the first derivative of the scale factor.
/// \param[in]      d2adt2 Value of the second derivative of the scale factor.
/// \return         Constant cosmology.
template <typename Type, class Cosmology, class>
inline Cosmology Input::constantify(const unsigned int size, const Type tmin,
                                    const Type tmax, const Type a,
                                    const Type dadt, const Type d2adt2) {
    const std::array<Type, 4> values =
        std::array<Type, 4>({{tmin, a, dadt, d2adt2}});
    Cosmology result = Cosmology();
    Utility::parallelize(result.size(),
                         [=, &result, &size, &values](const unsigned int i) {
                             result[i].assign(size, values[i]);
                         });
    Utility::parallelize(std::get<0>(result).size(),
                         [=, &result, &size, &tmin, &tmax](const unsigned int i) {
                             std::get<0>(result)[i] +=
                                 Type(i) * (tmax - tmin) / Type(size - (size > 0));
                         });
    return result;
}

// Data conversion to SI units
/// \brief          Data conversion to SI units.
/// \details        Converts a data to one expressed in SI units.
/// \tparam         Data Data type.
/// \tparam         Type Arithmetic type.
/// \param[in]      data Input data.
/// \param[in]      a Scale factor.
/// \param[in]      h Hubble parameter.
/// \param[in]      omegam Density parameter of matter.
/// \param[in]      lboxmpch Box size in megaparsecs times h.
/// \param[in]      mpc Value of one megaparsec in SI units.
/// \param[in]      rhoch2 Value of critical density times h squared in SI
///                 units.
/// \return         Data in SI units.
template <class Parameter, class Data, typename Type, class>
inline Data Input::sistemize(const Parameter &parameters, const Data &data,
                             const Type a, const Type h, const Type omegam,
                             const Type lboxmpch) {
    const Type t =
        Type(a * a * parameters.mpc) / Type(h * std::hecto::num * std::kilo::num);
    const Type l = Type(a * lboxmpch * parameters.mpc) / Type(h);
    Data result = data;
    result.phi() *= (l * l) / (t * t);
    result.dphidx() *= (l / (t * t)) * (-a);
    result.dphidy() *= (l / (t * t)) * (-a);
    result.dphidz() *= (l / (t * t)) * (-a);
    return result;
}

// Octree conversion to SI units
/// \brief          Octree conversion to SI units.
/// \details        Converts each element of the octree to one based on data
///                 expressed in SI units.
/// \tparam         Octree Octree type.
/// \tparam         Type Arithmetic type.
/// \tparam         Element Underlying element type.
/// \tparam         Data Data type.
/// \param[in,out]  octree Octree of data.
/// \param[in]      h Hubble parameter.
/// \param[in]      omegam Density parameter of matter.
/// \param[in]      lboxmpch Box size in megaparsecs times h.
/// \param[in]      mpc Value of one megaparsec in SI units.
/// \param[in]      rhoch2 Value of critical density times h squared in SI
///                 units.
/// \return         Octree in SI units.
template <class Parameter, class Octree, typename Type, class Element,
          class Data, class>
inline unsigned int Input::sistemize(const Parameter &parameters,
                                     Octree &octree, const Type h,
                                     const Type omegam, const Type lboxmpch) {
    Utility::parallelize(
        octree.begin(), octree.end(),
        [=, &h, &omegam, &lboxmpch, &parameters](Element &element) {
            std::get<1>(element) =
                sistemize(parameters, std::get<1>(element),
                          Type(std::get<1>(element).a()), h, omegam, lboxmpch);
        });
    return octree.size();
}

// Data homogenization
/// \brief          Data homogenization.
/// \details        Converts a data to one of a homogeneous empty Universe.
/// \tparam         Data Data type.
/// \param[in]      data Input data.
/// \return         Homogeneous empty data.
template <class Data, class>
inline Data Input::homogenize(const Data &data) {
    thread_local const Data zero = Data();
    Data result = data;
    result.rho() = zero.rho();
    result.phi() = zero.phi();
    result.dphidx() = zero.dphidx();
    result.dphidy() = zero.dphidy();
    result.dphidz() = zero.dphidz();
    result.dphidt() = zero.dphidt();
    result.a() =
        (result.a() != zero.a()) ? (result.a()) : (!(result.a() != zero.a()));
    return result;
}

// Octree homogenization
/// \brief          Octree homogenization.
/// \details        Converts each data of the octree to one of a homogeneous
///                 empty Universe.
/// \tparam         Octree Octree type.
/// \tparam         Dummy Dummy type.
/// \tparam         Element Underlying element type.
/// \tparam         Data Data type.
/// \param[in,out]  octree Octree of data.
/// \return         Homogeneous empty octree.
template <class Octree, class... Dummy, class Element, class Data, class>
inline unsigned int Input::homogenize(Octree &octree, Dummy...) {
    Utility::parallelize(octree.begin(), octree.end(), [](Element &element) {
        std::get<1>(element) = homogenize(std::get<1>(element));
    });
    return octree.size();
}

// Data schwarzschild setup
/// \brief          Data schwarzschild setup.
/// \details        Converts a data to one of a schwarzschild configuration.
/// \tparam         Extent Extent of the hyperoctree.
/// \tparam         Data Data type.
/// \tparam         Vector Vector type.
/// \tparam         Type Arithmetic type.
/// \tparam         Dimension Number of dimensions.
/// \param[in]      data Input data.
/// \param[in]      center Cell center.
/// \param[in]      position Mass position.
/// \param[in]      mass Mass in SI units.
/// \param[in]      length Spatial length in SI units.
/// \return         Schwarzschild data.
template <class Extent, class Data, class Vector, typename Type,
          unsigned int Dimension, class>
inline Data Input::schwarzschildify(const Data &data, const Vector &center,
                                    const Vector &position, const Type mass,
                                    const Type length) {
    static const Type extent = Type(Extent::num) / Type(Extent::den);
    static const Type null = Type();
    thread_local const Data zero = Data();
    const Vector v = Utility::join<Dimension>(position, center);
    const Type r = Utility::distance<Dimension>(position, center);
    const Type l = length / extent;
    const Type d = r * l;
    const Type gm = magrathea::Constants<Type>::g() * mass;
    Data result = data;
    result.rho() = data.rho();
    result.phi() = (d > null) ? (-gm / d) : (std::numeric_limits<Type>::max());
    result.dphidx() = (d > null) ? ((gm / (d * d)) * (v[0] / r)) : (null);
    result.dphidy() = (d > null) ? ((gm / (d * d)) * (v[1] / r)) : (null);
    result.dphidz() = (d > null) ? ((gm / (d * d)) * (v[2] / r)) : (null);
    result.a() =
        (result.a() != zero.a()) ? (result.a()) : (!(result.a() != zero.a()));
    return result;
}

// Octree schwarzschild setup
/// \brief          Octree schwarzschild setup.
/// \details        Converts each data of the octree to one of a schwarzschild
///                 configuration.
/// \tparam         Octree Octree type.
/// \tparam         Vector Vector type.
/// \tparam         Type Arithmetic type.
/// \tparam         Function Function type.
/// \tparam         Dummy Dummy type.
/// \tparam         Element Underlying element type.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Extent Extent of the hyperoctree.
/// \param[in,out]  octree Octree of data.
/// \param[in]      position Mass position in the octree.
/// \param[in]      mass Mass in SI units.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      refiner Refinement function taking a data and a level as
///                 arguments and returning true when a refinement should be
///                 triggered
/// \return         Schwarzschild octree.
template <class Octree, class Vector, typename Type, class Function,
          class... Dummy, class Element, unsigned int Dimension, class Extent,
          class>
inline unsigned int
Input::schwarzschildify(Octree &octree, const Vector &position, const Type mass,
                        const Type length, Function &&refiner, Dummy...) {
    unsigned int size = octree.size();
    std::vector<unsigned char> refine(size);
    do {
        size = octree.size();
        Utility::parallelize(
            octree.begin(), octree.end(),
            [=, &position, &length](Element &element) {
                std::get<1>(element) = schwarzschildify<Extent>(
                    std::get<1>(element),
                    std::array<Type, Dimension>({{std::get<0>(element).position(0),
                                                  std::get<0>(element).position(1),
                                                  std::get<0>(element).position(2)}}),
                    position, mass, length);
            });
        refine.resize(size);
        Utility::parallelize(size, [=, &octree, &refiner,
                                    &refine](const unsigned int i) {
            refine[i] =
                ((octree.leaf(octree.begin() + i)) &&
                 (std::get<0>(octree[i]).level() <
                  std::get<0>(octree[i]).refinements()) &&
                 (refiner(std::get<1>(octree[i]), std::get<0>(octree[i]).level())));
        });
        for (unsigned int irefine = 0; irefine < size; ++irefine) {
            if (refine[irefine]) {
                octree.refine(octree.begin() + irefine);
            }
        }
        octree.update();
    } while (size < octree.size());
    Utility::parallelize(
        octree.begin(), octree.end(), [=, &position, &length](Element &element) {
            std::get<1>(element) = schwarzschildify<Extent>(
                std::get<1>(element),
                std::array<Type, Dimension>({{std::get<0>(element).position(0),
                                              std::get<0>(element).position(1),
                                              std::get<0>(element).position(2)}}),
                position, mass, length);
        });
    return octree.size();
}
// -------------------------------------------------------------------------- //

// ---------------------------------- FILES --------------------------------- //
// File tree
/// \brief          File tree.
/// \details        Fills in an octree with the names of the existing ramses
///                 files.
/// \tparam         Octree Octree type.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Dimension Number of dimensions.
/// \param[out]     octree Octree of file names.
/// \param[in]      directory Input directory.
/// \param[in]      format Input file format.
/// \return         Number of detected files.
template <class Octree, class Element, class Index, class Data,
          unsigned int Dimension, class>
unsigned int Input::filetree(Octree &octree, const std::string &directory,
                             const std::string &format) {
    // Initialization
    const bool zeroth = magrathea::FileList(format, 0, 1, 0, directory).count();
    const bool first = magrathea::FileList(format, 0, 1, 1, directory).count();
    magrathea::FileList list(format, 0,
                             std::numeric_limits<int>::max() * (zeroth || first),
                             (!zeroth) && (first), directory);
    unsigned int level = 0;
    unsigned int current = 0;
    unsigned int size = 0;
    std::vector<unsigned int> index;

    // Find existing files and prepare the tree
    if (zeroth || first) {
        while (magrathea::FileSystem::exist(list[Index::indices(level + 1) - 1])) {
            ++level;
        }
        list.resize(Index::indices(level));
        size = octree.assign(0, level).size();
        index.resize(size, 0);
    }

    // Computation of the file index in each cell
    for (unsigned int ielem = 0; ielem < size; ++ielem) {
        current = std::get<0>(octree[ielem]).level();
        if (current == level) {
            for (unsigned int idim = 0; idim < Dimension; ++idim) {
                index[ielem] += std::get<0>(octree[ielem]).coordinate(idim)
                                << ((Dimension - (idim + 1)) * current);
            }
            ++index[ielem];
        }
    }

    // Finalization
    Utility::parallelize(size, [=, &octree, &index, &list](const unsigned int i) {
        std::get<1>(octree[i]) = (index[i]) ? (list[index[i] - 1]) : (Data());
    });
    return list.size();
}

// File list preparation
/// \brief          File list preparation.
/// \details        Adds to the list, the octree files which intersects the
///                 provided sphere and cone.
/// \tparam         List File list type.
/// \tparam         Octree Octree type.
/// \tparam         Sphere Sphere type.
/// \tparam         Conic Cone type.
/// \tparam         Dimension Number of dimensions.
/// \param[in]      list File list.
/// \param[in]      octree Input octree.
/// \param[in]      sphere Geometrical sphere.
/// \param[in]      cone Three dimensional cone.
/// \return         True if some files have been added to the list, false
///                 otherwise.
template <class List, class Octree, class Sphere, class Conic,
          unsigned int Dimension, class>
bool Input::prepare(List &list, const Octree &octree, const Sphere &sphere,
                    const Conic &conic) {
    // Initialization
    const unsigned int size = octree.size();
    const unsigned int original = list.size();
    std::vector<unsigned int> selection(size);
    unsigned int n = original;

    // Compute files to be read
    Utility::parallelize(
        size, [=, &octree, &sphere, &conic, &selection](const unsigned int i) {
            selection[i] =
                (std::get<1>(octree[i]).empty())
                    ? (0)
                    : (collide(octree, std::get<0>(octree[i]), sphere, conic));
        });
    list.resize(original +
                std::accumulate(selection.begin(), selection.end(), 0));
    std::for_each(selection.begin(), selection.end(),
                  [=, &n](unsigned int &i) { i = i ? ++n : i; });
    Utility::parallelize(size,
                         [=, &list, &octree, &selection](const unsigned int i) {
                             if (selection[i])
                                 list[selection[i] - 1] = std::get<1>(octree[i]);
                         });

    // Finalization
    return (list.size() > original);
}
// -------------------------------------------------------------------------- //

// ---------------------------------- DATA ---------------------------------- //
// Ramses importation
/// \brief          Ramses importation.
/// \details        Imports raw data from ramses gravity files. All cells
///                 selected by the provided filter are added to the octree. As
///                 there is no way to detect the coarse level, it should be
///                 specified as an argument.
/// \tparam         Integral Integral type of the file.
/// \tparam         Real Real type of the file.
/// \tparam         Parameter Parameter type.
/// \tparam         Octree Octree type.
/// \tparam         Function Function type taking an element as argument and
///                 returning a boolean.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Dimension Number of dimensions.
/// \param[in]      parameters Parameter structure.
/// \param[in,out]  octree Octree of data.
/// \param[in]      filename Input file name.
/// \param[in]      filter Filtering algorithm of cells.
/// \return         True on success, false on error.
template <typename Integral, typename Real, class Parameter, class Octree,
          class Function, class Element, class Index, class Data,
          unsigned int Dimension, class>
bool Input::import(const Parameter &parameters, Octree &octree,
                   const std::string &filename, Function &&filter) {
    // Initialization
    Integral record = Integral();
    Integral cube = Integral();
    Integral size = Integral();
    Integral n = Integral();
    std::vector<Integral> count;
    std::array<Real, Dimension + Dimension> minmax;
    std::vector<Real> center;
    std::vector<Real> force;
    std::vector<Real> a;
    std::vector<Real> phi;
    std::vector<Real> rho;
    std::vector<Integral> son;
    std::vector<Integral> selection;
    std::vector<Integral> index;
    std::ifstream stream(filename, std::ios::binary);
    std::thread thread;
    // Read file
    if (stream) {

        // Read level count
        magrathea::DataHandler::read(stream, record);
        count.resize(record / sizeof(Integral));
        magrathea::DataHandler::rread(stream, count.data(),
                                      count.data() + count.size());
        magrathea::DataHandler::read(stream, record);

        // Read cube number
        magrathea::DataHandler::read(stream, record);
        magrathea::DataHandler::read(stream, cube);
        magrathea::DataHandler::read(stream, record);

        // Read cube min and max
        magrathea::DataHandler::read(stream, record);
        magrathea::DataHandler::read(stream, minmax);
        magrathea::DataHandler::read(stream, record);

        // Loop over levels
        const unsigned int levelmax = parameters.coarseonly
                                          ? parameters.ncoarse + 1
                                          : parameters.ncoarse + count.size();
        for (unsigned int ilevel = parameters.ncoarse; ilevel < levelmax;
             ++ilevel) {

            // Initialize the loop
            size = count[ilevel - parameters.ncoarse];
            thread = std::thread([=, &size, &force, &a, &phi, &rho, &son]() {
                force.resize(Dimension * size);
                a.resize(size);
                phi.resize(size);
                rho.resize(size);
                son.resize(size);
            });
            center.resize(Dimension * size);

            // Read center
            magrathea::DataHandler::read(stream, record);
            magrathea::DataHandler::rread(stream, center.data(),
                                          center.data() + center.size());
            magrathea::DataHandler::read(stream, record);

            // Thread management
            thread.join();
            thread = std::thread([=, &size, &selection, &index]() {
                selection.resize(size);
                std::iota(selection.begin(), selection.end(), Integral());
                index.resize(size);
                std::iota(index.begin(), index.end(), Integral());
            });

            // Read force
            magrathea::DataHandler::read(stream, record);
            magrathea::DataHandler::rread(stream, force.data(),
                                          force.data() + force.size());
            magrathea::DataHandler::read(stream, record);

            // Read a
            magrathea::DataHandler::read(stream, record);
            magrathea::DataHandler::rread(stream, a.data(), a.data() + a.size());
            magrathea::DataHandler::read(stream, record);

            // Read phi
            magrathea::DataHandler::read(stream, record);
            magrathea::DataHandler::rread(stream, phi.data(),
                                          phi.data() + phi.size());
            magrathea::DataHandler::read(stream, record);

            // Read rho
            magrathea::DataHandler::read(stream, record);
            magrathea::DataHandler::rread(stream, rho.data(),
                                          rho.data() + rho.size());
            magrathea::DataHandler::read(stream, record);

            // Read son
            magrathea::DataHandler::read(stream, record);
            magrathea::DataHandler::rread(stream, son.data(),
                                          son.data() + son.size());
            magrathea::DataHandler::read(stream, record);

            // Selection
            thread.join();
            Utility::parallelize(
                index.begin(), index.end(),
                [=, &filter, &center, &force, &a, &phi, &rho, &son,
                 &selection](Integral &i) {
                    selection[i] = filter(
                        Element(Index::compute(ilevel, center[Dimension * i],
                                               center[Dimension * i + 1],
                                               center[Dimension * i + 2]),
                                Data(rho[i], phi[i],
                                     std::array<Real, 3>({{force[Dimension * i],
                                                           force[Dimension * i + 1],
                                                           force[Dimension * i + 2]}}),
                                     a[i])));
                });

            // Destination
            n = octree.size();
            size = std::accumulate(selection.begin(), selection.end(), n);
            thread = std::thread([=, &octree, &size]() { octree.resize(size); });
            std::for_each(selection.begin(), selection.end(),
                          [=, &n](Integral &i) { i = i ? ++n : i; });
            thread.join();
            Utility::parallelize(
                index.begin(), index.end(),
                [=, &octree, &center, &force, &a, &phi, &rho, &son,
                 &selection](Integral &i) {
                    if (selection[i])
                        octree[selection[i] - 1] = Element(
                            Index::compute(ilevel, center[Dimension * i],
                                           center[Dimension * i + 1],
                                           center[Dimension * i + 2]),
                            Data(rho[i], phi[i],
                                 std::array<Real, 3>(
                                     {{force[Dimension * i], force[Dimension * i + 1],
                                       force[Dimension * i + 2]}}),
                                 a[i]));
                });
        }

        // Cleaning
        stream.close();
    }

    // Finalization
    return stream.good();
}

// Ramses importation
/// \brief          Ramses importation.
/// \details        Imports HDF5 data from ramses gravity files. All cells
///                 selected by the provided filter are added to the octree.
/// \tparam         Position Position.
/// \tparam         Extent Extent type.
/// \tparam         Integral Integral type of the file.
/// \tparam         Real Real type of the file.
/// \tparam         Conic Cone type.
/// \tparam         Octree Octree type.
/// \tparam         Function Function type taking an element as argument and
///                 returning a boolean.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Dimension Number of dimensions.
/// \param[in]      parameters Parameter structure.
/// \param[in]      rank Rank.
/// \param[in]      rotm1 rotation matrix for narrow cones
/// \param[in,out]  thetay Semi-angle for solid angle in direction y
/// \param[in,out]  thetaz Semi-angle for solid angle in direction z
/// \param[in]      conic Cone characteristics
/// \param[in,out]  octree Octree of data.
/// \param[in]      filename Input file name.
/// \param[in]      filter Filtering algorithm of cells.
/// \return         True on success, false on error.
template <class Position, class Extent, typename Integral, typename Real,
          class Parameter, class Conic, class Octree, class Function,
          class Element, class Index, class Data, unsigned int Dimension, class>
bool Input::importhdf5(const Parameter &parameters, unsigned int &rank,
                       const std::array<std::array<double, 3>, 3> &rotm1,
                       const double &thetay, const double &thetaz,
                       const Conic &conic, Octree &octree,
                       const std::string &filename, Function &&filter) {

    // Initialization
    Integral size = Integral();
    Integral n = Integral();
    Real aexp(0);
    std::thread thread;
    const int levelShift(std::log2(EXTENT));
    unsigned int levelMin(0);

    TReadHDF5::getAttribute(filename, "metadata/cone_info", "levelmin", levelMin);
    if (!(levelMin = parameters.ncoarse)) {
        std::cout << "# Coarse level problem between data and 'ncoarse' parameter"
                  << std::endl;
        std::cout << "# Error at file " << __FILE__ << " line : " << __LINE__
                  << std::endl;
        std::terminate();
    }

    std::vector<std::vector<std::string>> cubeNumber;
    std::vector<std::string> cubeNb;
    std::vector<Integral> count;

    // Select cubes containing cells inside the cone
    TReadHDF5::cellsAndCubesPerLevels(parameters, filename, count, cubeNumber,
                                      thetay, thetaz, conic);
    // If we only keep coarse cells, do not count higher levels
    if (parameters.coarseonly) {
        while (count.size() != 1) {
            count.pop_back();
        }
    }
#ifdef VERBOSE
    std::cout << "# Rank : " << rank << " Total number of cells : "
              << std::accumulate(count.begin(), count.end(), 0) << " Filename "
              << filename << std::endl;
    for (unsigned int i = 0; i < count.size(); i++) {
        std::cout << "# Rank : " << rank
                  << " Total number of cells at level coarse + " << i << " is "
                  << count[i] << std::endl;
    }
#endif
    // Loop over levels
    for (unsigned int ilevel = parameters.ncoarse;
         ilevel < parameters.ncoarse + count.size(); ++ilevel) {

        std::vector<Real> center;
        std::vector<Real> force;
        std::vector<Real> a;
        std::vector<Real> phi;
        std::vector<Real> rho;
        std::vector<Integral> selection;
        std::vector<Integral> index;
        // Initialize the loop
        const unsigned int octreelevel = ilevel + levelShift;
        size = count[ilevel - parameters.ncoarse];
        cubeNb = cubeNumber[ilevel - parameters.ncoarse];
        // Thread management
        thread = std::thread([=, &a, &size, &selection, &index]() {
            a.resize(size), selection.resize(size);
            std::iota(selection.begin(), selection.end(), Integral());
            index.resize(size);
            std::iota(index.begin(), index.end(), Integral());
        });
        thread.join();

        // Read scale factor
        TReadHDF5::getAttribute(filename, "metadata/cone_info", "aexp", aexp);
        // Give the same scale factor to all cells in same shell
        if (aexp < 1) {
            std::fill(a.begin(), a.end(), aexp);
            // The shell might be slightly (one coarse step) ahead of the lightcone,
            // in this case the scale factor might be above 1
        } else {
            std::fill(a.begin(), a.end(), 1);
        }

        // Read position, potential, force and density from data
        TReadHDF5::fillVectors_grav(filename, ilevel, ilevel, cubeNb,
                                    "position_cell", center,
                                    "gravitational_field_cell", force,
                                    "potential_cell", phi, "density_cell", rho);

        if (phi.size() != size || center.size() != Dimension * size ||
            force.size() != Dimension * size) {
            std::cout << "# Problem size data " << phi.size() << " " << size << " "
                      << center.size() << " " << Dimension * size << " "
                      << force.size() << " " << Dimension * size << std::endl;
            std::cout << "# Error at file " << __FILE__ << " line : " << __LINE__
                      << std::endl;
            std::terminate();
        }

        // For narrow cones, need to rotate position and force
        if (parameters.isfullsky == 0) {
            std::vector<Real> centertmp = center;
            std::vector<Real> forcetmp = force;
            Utility::parallelize(
                size, [=, &center, &force, &centertmp, &forcetmp](Integral i) {
                    center[3 * i] = centertmp[3 * i] * rotm1[0][0] +
                                    centertmp[3 * i + 1] * rotm1[0][1] +
                                    centertmp[3 * i + 2] * rotm1[0][2];
                    center[3 * i + 1] = centertmp[3 * i] * rotm1[1][0] +
                                        centertmp[3 * i + 1] * rotm1[1][1] +
                                        centertmp[3 * i + 2] * rotm1[1][2];
                    center[3 * i + 2] = centertmp[3 * i] * rotm1[2][0] +
                                        centertmp[3 * i + 1] * rotm1[2][1] +
                                        centertmp[3 * i + 2] * rotm1[2][2];
                    force[3 * i] = forcetmp[3 * i] * rotm1[0][0] +
                                   forcetmp[3 * i + 1] * rotm1[0][1] +
                                   forcetmp[3 * i + 2] * rotm1[0][2];
                    force[3 * i + 1] = forcetmp[3 * i] * rotm1[1][0] +
                                       forcetmp[3 * i + 1] * rotm1[1][1] +
                                       forcetmp[3 * i + 2] * rotm1[1][2];
                    force[3 * i + 2] = forcetmp[3 * i] * rotm1[2][0] +
                                       forcetmp[3 * i + 1] * rotm1[2][1] +
                                       forcetmp[3 * i + 2] * rotm1[2][2];
                });
        }

        // Selection
        Utility::parallelize(
            index.begin(), index.end(),
            [=, &filter, &center, &force, &a, &phi, &rho, &selection](Integral &i) {
                selection[i] = filter(
                    Element(Index::template compute<Real, Position, Extent>(
                                octreelevel, center[Dimension * i],
                                center[Dimension * i + 1], center[Dimension * i + 2]),
                            Data(rho[i], phi[i],
                                 std::array<Real, 3>(
                                     {{force[Dimension * i], force[Dimension * i + 1],
                                       force[Dimension * i + 2]}}),
                                 a[i])));
            });

        // Destination
        n = octree.size();
        size = std::accumulate(selection.begin(), selection.end(), n);
        thread = std::thread([=, &octree, &size]() { octree.resize(size); });
        std::for_each(selection.begin(), selection.end(),
                      [=, &n](Integral &i) { i = i ? ++n : i; });
        thread.join();
        Utility::parallelize(
            index.begin(), index.end(),
            [=, &octree, &center, &force, &a, &phi, &rho, &selection](Integral &i) {
                if (selection[i])
                    octree[selection[i] - 1] = Element(
                        Index::template compute<Real, Position, Extent>(
                            octreelevel, center[Dimension * i],
                            center[Dimension * i + 1], center[Dimension * i + 2]),
                        Data(rho[i], phi[i],
                             std::array<Real, 3>(
                                 {{force[Dimension * i], force[Dimension * i + 1],
                                   force[Dimension * i + 2]}}),
                             a[i]));
            });

    } // loop over levels

    return true;
}

// Ramses importation
/// \brief          Ramses importation.
/// \details        Imports HDF5 data from ramses gravity files. Used for small
///		    spherical buffer zone in narrow cones. All cells
///                 selected by the provided filter are added to the octree.
/// \tparam         Position Position type.
/// \tparam         Extent Extent type.
/// \tparam         Integral Integral type of the file.
/// \tparam         Real Real type of the file.
/// \tparam         Parameter Parameter type.
/// \tparam         Octree Octree type.
/// \tparam         Function Function type taking an element as argument and
///                 returning a boolean.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Dimension Number of dimensions.
/// \param[in]      parameters Parameter structure.
/// \param[in,out]  octree Octree of data.
/// \param[in]      filename Input file name.
/// \param[in]      filter Filtering algorithm of cells.
/// \return         True on success, false on error.
template <class Position, class Extent, typename Integral, typename Real,
          class Parameter, class Octree, class Function, class Element,
          class Index, class Data, unsigned int Dimension, class>
bool Input::importfullhdf5(const Parameter &parameters, Octree &octree,
                           const std::string &filename, Function &&filter) {

    // Initialization
    Integral size = Integral();
    Integral n = Integral();
    std::vector<Integral> count;
    Real aexp(0);
    std::thread thread;
    const int levelShift(std::log2(EXTENT));

    // Read file
    unsigned int levelMin(0);
    TReadHDF5::getAttribute(filename, "metadata/cone_info", "levelmin", levelMin);
    if (!(levelMin = parameters.ncoarse)) {
        std::cout << "# Coarse level problem between data and 'ncoarse' parameter"
                  << std::endl;
        std::cout << "# Error at file " << __FILE__ << " line : " << __LINE__
                  << std::endl;
        std::terminate();
    }

    // Get all cells in data
    TReadHDF5::cellsPerLevels(filename, count);
    // If we only keep coarse cells, do not count higher levels
    if (parameters.coarseonly) {
        while (count.size() != 1) {
            count.pop_back();
        }
    }

    // Loop over levels
    for (unsigned int ilevel = parameters.ncoarse;
         ilevel < parameters.ncoarse + count.size(); ++ilevel) {

        std::vector<Real> center;
        std::vector<Real> force;
        std::vector<Real> a;
        std::vector<Real> phi;
        std::vector<Real> rho;
        std::vector<Integral> selection;
        std::vector<Integral> index;
        const unsigned int octreelevel = ilevel + levelShift;
        size = count[ilevel - parameters.ncoarse];

        // Thread management
        thread = std::thread([=, &a, &size, &selection, &index]() {
            a.resize(size), selection.resize(size);
            std::iota(selection.begin(), selection.end(), Integral());
            index.resize(size);
            std::iota(index.begin(), index.end(), Integral());
        });
        thread.join();

        // Read scale factor
        TReadHDF5::getAttribute(filename, "metadata/cone_info", "aexp", aexp);
        // Give the same scale factor to all cells in same shell
        if (aexp < 1) {
            std::fill(a.begin(), a.end(), aexp);
            // The shell might be slightly (one coarse step) ahead of the lightcone,
            // in this case the scale factor might be above 1
        } else {
            std::fill(a.begin(), a.end(), 1);
        }

        // Read position, potential, force and density from data
        TReadHDF5::fillVectors_grav(
            filename, ilevel, ilevel, "gravitational_field_cell", force,
            "potential_cell", phi, "density_cell", rho, "position_cell", center);

        if (phi.size() != size || center.size() != Dimension * size ||
            force.size() != Dimension * size) {
            std::cout << "# Problem size data " << phi.size() << " " << size << " "
                      << center.size() << " " << Dimension * size << " "
                      << force.size() << " " << Dimension * size << std::endl;
            std::cout << "# Error at file " << __FILE__ << " line : " << __LINE__
                      << std::endl;
            std::terminate();
        }

        // Selection
        Utility::parallelize(
            index.begin(), index.end(),
            [=, &filter, &center, &force, &a, &phi, &rho, &selection](Integral &i) {
                selection[i] = filter(
                    Element(Index::template compute<Real, Position, Extent>(
                                octreelevel, center[Dimension * i],
                                center[Dimension * i + 1], center[Dimension * i + 2]),
                            Data(rho[i], phi[i],
                                 std::array<Real, 3>(
                                     {{force[Dimension * i], force[Dimension * i + 1],
                                       force[Dimension * i + 2]}}),
                                 a[i])));
            });
        // Destination
        n = octree.size();
        size = std::accumulate(selection.begin(), selection.end(), n);
        thread = std::thread([=, &octree, &size]() { octree.resize(size); });
        std::for_each(selection.begin(), selection.end(),
                      [=, &n](Integral &i) { i = i ? ++n : i; });
        thread.join();
        Utility::parallelize(
            index.begin(), index.end(),
            [=, &octree, &center, &force, &a, &phi, &rho, &selection](Integral &i) {
                if (selection[i])
                    octree[selection[i] - 1] = Element(
                        Index::template compute<Real, Position, Extent>(
                            octreelevel, center[Dimension * i],
                            center[Dimension * i + 1], center[Dimension * i + 2]),
                        Data(rho[i], phi[i],
                             std::array<Real, 3>(
                                 {{force[Dimension * i], force[Dimension * i + 1],
                                   force[Dimension * i + 2]}}),
                             a[i]));
            });

    } // loop over levels

    return true;
}

// Ramses importation
/// \brief          Ramses importation.
/// \details        Imports ASCII data from ramses gravity files. All cells
///                 selected by the provided filter are added to the octree. As
///                 there is no way to detect the coarse level, it should be
///                 specified as an argument.
/// \tparam         Position Position type.
/// \tparam         Extent Extent type.
/// \tparam         Integral Integral type of the file.
/// \tparam         Real Real type of the file.
/// \tparam         Octree Octree type.
/// \tparam         Function Function type taking an element as argument and
///                 returning a boolean.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Dimension Number of dimensions.
/// \param[in]      parameters Parameter structure
/// \param[in]      rotm1 rotation matrix for narrow cones
/// \param[in,out]  octree Octree of data.
/// \param[in]      filename Input file name.
/// \param[in]      filter Filtering algorithm of cells.
/// \return         True on success, false on error.
template <class Position, class Extent, typename Integral, typename Real,
          class Parameter, class Octree, class Function, class Element,
          class Index, class Data, unsigned int Dimension, class>
bool Input::importascii(const Parameter &parameters,
                        const std::array<std::array<double, 3>, 3> &rotm1,
                        Octree &octree, const std::string &filename,
                        Function &&filter) {

    Integral size = Integral();
    Integral n = Integral();
    std::thread thread;
    const Integral levelShift(std::log2(EXTENT));

    // Initialization
    std::ifstream streaming(filename.c_str());
    streaming.unsetf(std::ios_base::skipws);
    size = std::count(std::istream_iterator<char>(streaming),
                      std::istream_iterator<char>(), '\n');
    streaming.close();
    std::ifstream stream(filename.c_str());

    std::vector<Real> center;
    std::vector<Real> force;
    std::vector<Real> a;
    std::vector<Real> phi;
    std::vector<Real> rho;
    std::vector<Integral> son;
    std::vector<Integral> level;
    std::vector<Integral> index;
    std::vector<Integral> selection;

    thread =
        std::thread([=, &size, &level, &center, &force, &a, &phi, &rho, &son]() {
            level.resize(size);
            center.resize(Dimension * size);
            force.resize(Dimension * size);
            a.resize(size);
            phi.resize(size);
            rho.resize(size);
            son.resize(size);
        });
    thread.join();

    // Read all cells in files
    for (unsigned int i = 0; i < size; ++i) {
        stream >> level[i] >> center[Dimension * i] >> center[Dimension * i + 1] >>
            center[Dimension * i + 2] >> rho[i] >> phi[i] >> force[Dimension * i] >>
            force[Dimension * i + 1] >> force[Dimension * i + 2] >> a[i];
    }
    // For narrow cones, need to rotate position and force
    if (parameters.isfullsky == 0) {
        std::vector<Real> centertmp = center;
        std::vector<Real> forcetmp = force;
        Utility::parallelize(
            size, [=, &center, &force, &centertmp, &forcetmp](Integral i) {
                center[3 * i] = centertmp[3 * i] * rotm1[0][0] +
                                centertmp[3 * i + 1] * rotm1[0][1] +
                                centertmp[3 * i + 2] * rotm1[0][2];
                center[3 * i + 1] = centertmp[3 * i] * rotm1[1][0] +
                                    centertmp[3 * i + 1] * rotm1[1][1] +
                                    centertmp[3 * i + 2] * rotm1[1][2];
                center[3 * i + 2] = centertmp[3 * i] * rotm1[2][0] +
                                    centertmp[3 * i + 1] * rotm1[2][1] +
                                    centertmp[3 * i + 2] * rotm1[2][2];
                force[3 * i] = forcetmp[3 * i] * rotm1[0][0] +
                               forcetmp[3 * i + 1] * rotm1[0][1] +
                               forcetmp[3 * i + 2] * rotm1[0][2];
                force[3 * i + 1] = forcetmp[3 * i] * rotm1[1][0] +
                                   forcetmp[3 * i + 1] * rotm1[1][1] +
                                   forcetmp[3 * i + 2] * rotm1[1][2];
                force[3 * i + 2] = forcetmp[3 * i] * rotm1[2][0] +
                                   forcetmp[3 * i + 1] * rotm1[2][1] +
                                   forcetmp[3 * i + 2] * rotm1[2][2];
            });
    }

    // Correct level depending on Octree extent
    std::transform(level.begin(), level.end(), level.begin(),
                   bind2nd(std::plus<Integral>(), levelShift));

    // Selection
    thread = std::thread([=, &size, &selection, &index]() {
        selection.resize(size);
        std::iota(selection.begin(), selection.end(), Integral());
        index.resize(size);
        std::iota(index.begin(), index.end(), Integral());
    });
    thread.join();
    Utility::parallelize(
        index.begin(), index.end(),
        [=, &filter, &level, &center, &force, &a, &phi, &rho, &son,
         &selection](Integral &i) {
            selection[i] = filter(
                Element(Index::template compute<Real, Position, Extent>(
                            level[i], center[Dimension * i],
                            center[Dimension * i + 1], center[Dimension * i + 2]),
                        Data(rho[i], phi[i],
                             std::array<Real, 3>(
                                 {{force[Dimension * i], force[Dimension * i + 1],
                                   force[Dimension * i + 2]}}),
                             a[i])));
        });
    // Destination
    n = octree.size();
    size = std::accumulate(selection.begin(), selection.end(), n);
    thread = std::thread([=, &octree, &size]() { octree.resize(size); });
    std::for_each(selection.begin(), selection.end(),
                  [=, &n](Integral &i) { i = i ? ++n : i; });
    thread.join();
    Utility::parallelize(
        index.begin(), index.end(),
        [=, &octree, &level, &center, &force, &a, &phi, &rho, &son,
         &selection](Integral &i) {
            if (selection[i])
                octree[selection[i] - 1] =
                    Element(Index::template compute<Real, Position, Extent>(
                                level[i], center[Dimension * i],
                                center[Dimension * i + 1], center[Dimension * i + 2]),
                            Data(rho[i], phi[i],
                                 std::array<Real, 3>(
                                     {{force[Dimension * i], force[Dimension * i + 1],
                                       force[Dimension * i + 2]}}),
                                 a[i]));
        });

    // If we only want coarse level, remove refined cells
    if (parameters.coarseonly) {
        const unsigned int coarseShift = parameters.ncoarse + levelShift;
        octree.resize(std::distance(
            std::begin(octree),
            std::remove_if(std::begin(octree), std::end(octree),
                           [=, &coarseShift](const Element &elem) {
                               return std::get<0>(elem).level() != coarseShift;
                           })));
    }

    stream.close();
    return stream.good();
}

// Cosmology acquisition
/// \brief          Cosmology acquisition.
/// \details        Acquires cosmology parameters and cosmology from
///                 input files: the Hubble parameter, the density
///                 parameter of matter and the size of the box. It
///                 returns a container with t, a(t), its first and
///                 second derivatives (where t can be either the
///                 cosmic time or the conformal time eta).
/// \tparam         Type Arithmetic type.
/// \tparam         Cosmology Cosmology container type.
/// \tparam         Element Element type of the container.
/// \param[in]      simfile Name of a simulation file or directory to extract
///                 the box size.
/// \param[in]      paramfile Name of a file of cosmological parameters.
/// \param[in]      evolfile Name of a file of cosmological evolution.
/// \param[in]      h Hubble parameter.
/// \param[in]      omegam Density parameter of matter.
/// \param[in]      lboxmpch Box size in megaparsecs times h.
/// \param[in]      mpc Value of one megaparsec in SI units.
/// \param[in]      outfile Name of an output file for debugging purposes.
/// \return         Container of cosmology evolution and derivatives. Contains
/// t, a, H, int Hdot, Hdot
template <typename Type, class Parameter, class Cosmology, class Element, class>
Cosmology Input::acquire(const Parameter &parameters, Type &h, Type &omegam,
                         Type &lboxmpch, const std::string &outfile) {
    // Initialization
    static const bool reverse = true;
    static const Cosmology empty = Cosmology();
    static const Type zero = Type();
    static const Type one = Type(1);
    static const Type two = one + one;
    static const Type year = Type(365 * 24 * 60 * 60) + Type(24 * 60 * 60 / 4);
    static const Type twob2 = two * std::pow(one / Type(std::kilo::num), two);
    static const Type nan = std::numeric_limits<Type>::quiet_NaN();
    static const unsigned int nsmooth = 500;
    static const unsigned int nderive = 1000;
    static const unsigned int precision = std::numeric_limits<Type>::max_digits10;
    static const unsigned int width =
        precision + std::numeric_limits<unsigned char>::digits;
    static const unsigned int length = 5;
    static const std::string digits = "0123456789";
    static const char minus = '-';
    static const char point = '.';
    static const char comment = '#';
    static const char space = ' ';
    static const char endl = '\n';
    static const std::string boxlen = "boxlen";
    static const std::string htext = "h";
    static const std::string omegamtext = "Omega_m";
    static const std::string separator = "=";
    unsigned int size = zero;
    std::vector<std::string> line;
    // Put evolution data in vector
    std::thread thread([=, &line]() {
        std::ifstream input(parameters.evolfile);
        if (input) {
            do {
                line.emplace_back();
            } while (std::getline(input, line.back()));
            input.close();
        } else {
            std::cout << "# WARNING : evolfile not found: " << parameters.evolfile
                      << std::endl;
            std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                      << std::endl;
            std::terminate();
        }
    });
    std::ifstream stream;
    std::ofstream outstream;
    std::string string(
        std::begin(parameters.celldir) +
            std::min(parameters.celldir.rfind(boxlen), parameters.celldir.size()),
        std::end(parameters.celldir));
    std::string temporary;
    Type hubble = zero;
    std::array<std::vector<Type>, length> input;
    std::array<std::vector<Type>, length> output;
    Cosmology result = empty;

    // Find the box length
    lboxmpch = zero;
    string.assign(string.begin() + std::min(boxlen.size(), string.size()),
                  string.end());
    string.assign(string.begin(),
                  string.begin() +
                      std::min(string.find_first_not_of(digits + minus + "."),
                               string.size()));
    std::replace(string.begin(), string.end(), minus, point);
    lboxmpch = (!string.empty()) ? (std::stold(string)) : (zero);

    // Find h and omegam
    h = zero;
    omegam = zero;
    stream.open(parameters.paramfile);
    if (stream) {
        while (std::getline(stream, string) &&
               (!((h < zero) || (h > zero)) ||
                !((omegam < zero) || (omegam > zero)))) {
            if (string.find(separator) != std::string::npos) {
                temporary.assign(string.begin(),
                                 string.begin() + string.find(separator));
                temporary.assign(
                    std::find_if(temporary.begin(), temporary.end(),
                                 [](const char c) { return std::isgraph(c); }),
                    temporary.end());
                temporary.assign(
                    temporary.begin(),
                    temporary.end() -
                        std::distance(temporary.rbegin(),
                                      std::find_if(temporary.rbegin(), temporary.rend(),
                                                   [](const char c) {
                                                       return std::isgraph(c);
                                                   })));
                string.assign(
                    string.begin() + string.find(separator) +
                        (std::distance(string.begin() + string.find(separator),
                                       string.end()) > 0),
                    string.end());
                string.assign(
                    std::find_if(string.begin(), string.end(),
                                 [](const char c) { return std::isgraph(c); }),
                    string.end());
                string.assign(
                    string.begin(),
                    string.end() -
                        std::distance(string.rbegin(),
                                      std::find_if(string.rbegin(), string.rend(),
                                                   [](const char c) {
                                                       return std::isgraph(c);
                                                   })));
                if (temporary == htext) {
                    h = std::stold(string);
                    hubble = (h * std::hecto::num * std::kilo::num) / parameters.mpc;
                } else if (temporary == omegamtext) {
                    omegam = std::stold(string);
                }
            }
        }
        stream.close();
    } else {
        std::cout << "# WARNING : paramfile not found: " << parameters.paramfile
                  << std::endl;
        std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__
                  << std::endl;
        std::terminate();
    }
    // Find evolution
    thread.join();
    line.erase(std::remove_if(line.begin(), line.end(),
                              [](const std::string &s) { return s.empty(); }),
               line.end());
    size = line.size();
    if (size > 0) {
        Utility::parallelize(input.begin(), input.end(),
                             [=, &size](std::vector<Type> &v) { v.resize(size); });
        // Put evolution from string to float in 'input'
        Utility::parallelize(size, [=, &input](const unsigned int i) {
            std::istringstream iss(line[i]);
            iss >> std::get<0>(input)[i];
            iss >> std::get<1>(input)[i];
            iss >> std::get<2>(input)[i];
            iss >> std::get<3>(input)[i];
            iss >> std::get<4>(input)[i];
        });
        Utility::parallelize(output.begin(), output.end(),
                             [=, &size](std::vector<Type> &v) { v.resize(size); });
        // Convert Ramses Units to SI and modify order of columns in 'output'
        Utility::parallelize(
            size, [=, &input, &output, &hubble](const unsigned int i) {
                std::get<0>(output)[i] = std::get<4>(input)[i] / hubble;
                std::get<2>(output)[i] = std::get<0>(input)[i];
                std::get<3>(output)[i] = std::get<1>(input)[i] * hubble;
            });
        // Reverse order of data from file (increasing from observer toward the
        // past)
        std::get<0>(output) = Utility::reverse(
            std::get<0>(output), *std::min_element(std::get<0>(output).begin(),
                                                   std::get<0>(output).end()));
        // Smooth evolution
        std::get<3>(output) = Utility::smooth(
            std::get<2>(output), std::get<3>(output),
            [](Type xi, Type xj) {
                return std::exp(-((xi - xj) * (xi - xj)) / twob2);
            },
            size / nsmooth);
        std::get<1>(output)[0] = zero;
        // Check if each line is good and monotonically increasing
        Utility::parallelize(size, [=, &output](const unsigned int i) {
            if (i) {
                if (!std::isnormal(std::abs(std::get<2>(output)[i] -
                                            std::get<2>(output)[i - 1]))) {
                    std::get<0>(output)[i] = nan;
                }
            }
        });
        // If one value is NaN, then all line is NaN
        Utility::parallelize(size, [=, &output](const unsigned int i) {
            if (std::isnan(std::get<0>(output)[i] + std::get<1>(output)[i] +
                           std::get<2>(output)[i] + std::get<3>(output)[i] +
                           std::get<4>(output)[i])) {
                std::get<0>(output)[i] = nan;
                std::get<1>(output)[i] = nan;
                std::get<2>(output)[i] = nan;
                std::get<3>(output)[i] = nan;
                std::get<4>(output)[i] = nan;
            }
        });
        // Remove bad lines
        Utility::parallelize(
            output.begin(), output.end(), [](std::vector<Type> &v) {
                v.erase(std::remove_if(v.begin(), v.end(),
                                       [](const Type x) { return std::isnan(x); }),
                        v.end());
            });
        size = std::get<0>(output).size();
        thread = std::thread([=, &result, &size]() {
            std::get<0>(result).resize(size);
            std::get<1>(result).resize(size);
            std::get<2>(result).resize(size);
            std::get<3>(result).resize(size);
        });
        for (unsigned int i = 1; i < size; ++i) {
            std::get<1>(output)[i] =
                std::get<1>(output)[i - 1] +
                (std::get<2>(output)[i] - std::get<2>(output)[i - 1]) /
                    std::get<3>(output)[i - 1];
        }
        thread.join();
        std::get<4>(output) = Utility::derive(std::get<1>(output),
                                              std::get<3>(output), size / nderive);
        std::get<4>(output) = Utility::smooth(
            std::get<2>(output), std::get<4>(output),
            [](Type xi, Type xj) {
                return std::exp(-((xi - xj) * (xi - xj)) / twob2);
            },
            size / nsmooth);
        std::get<4>(output) = Utility::smooth(
            std::get<2>(output), std::get<4>(output),
            [](Type xi, Type xj) {
                return std::exp(-((xi - xj) * (xi - xj) * std::deca::num) / twob2);
            },
            size / (nsmooth * std::deca::num));
        std::get<3>(output) = Utility::integrate(
            std::get<1>(output), std::get<4>(output), std::get<3>(output)[0]);
        std::get<1>(output)[0] = zero;
        for (unsigned int i = 1; i < size; ++i) {
            std::get<1>(output)[i] =
                std::get<1>(output)[i - 1] +
                (std::get<2>(output)[i] - std::get<2>(output)[i - 1]) /
                    std::get<3>(output)[i - 1];
        }
        Utility::parallelize(size, [=, &output, &result](const unsigned int i) {
            std::get<0>(result)[i] = std::get<1>(output)[i];
            std::get<1>(result)[i] = std::get<2>(output)[i];
            std::get<2>(result)[i] = std::get<3>(output)[i];
            std::get<3>(result)[i] = std::get<4>(output)[i];
        });
        if (reverse) {
            std::get<0>(result) = Utility::reverse(
                std::get<0>(result), *std::max_element(std::get<0>(result).begin(),
                                                       std::get<0>(result).end()));
            Utility::parallelize(size, [=, &result](const unsigned int i) {
                std::get<0>(result)[i] = std::abs(std::get<0>(result)[i]);
                std::get<2>(result)[i] = -std::get<2>(result)[i];
            });
            Utility::parallelize(
                result.begin(), result.end(),
                [](std::vector<Type> &v) { std::reverse(v.begin(), v.end()); });
        }
    }
    // Finalization
    if (!((h < zero) || (h > zero)) || !((omegam < zero) || (omegam > zero)) ||
        !((lboxmpch < zero) || (lboxmpch > zero)) || (result.empty())) {
        h = zero;
        omegam = zero;
        lboxmpch = zero;
        result = empty;
    } else if (!outfile.empty()) {
        if (!magrathea::FileSystem::exist(outfile)) {
            outstream.open(outfile);
            if (outstream) {
                outstream << std::setprecision(precision);
                outstream << comment << space << h << space << omegam << space
                          << lboxmpch << space << (lboxmpch / h) << space
                          << parameters.mpc << space << hubble << space
                          << (one / hubble) << space
                          << (one / (hubble * std::giga::num * year)) << std::right
                          << std::endl;
                for (unsigned int i = 0; i < size; ++i) {
                    std::for_each(output.begin(), output.end(),
                                  [=, &outstream](const std::vector<Type> &v) {
                                      outstream << std::setw(width) << v[i];
                                  });
                    outstream << endl;
                }
                outstream.close();
            }
        }
    }
    return result;
}

// Parameter file parsing
/// \brief          Parameter file parsing.
/// \details        Parses the provided parameter file and returns a map
///                 of parameters.
/// \tparam         Container Output map type.
/// \tparam         Element Element type.
/// \param[in]      filename File name.
/// \param[in]      separator Separator string.
/// \param[in]      comment Comment string.
/// \return         Map of parameters.
template <class Container, class Element, class>
Container Input::parse(const std::string &filename,
                       const std::string &separator,
                       const std::string &comment) {
    // Initialization
    std::vector<std::string> line;
    std::vector<Element> mapping;
    Container result;
    unsigned int size = 0;
    std::ifstream stream(filename);

    // Parse file
    if (stream) {
        do {
            line.emplace_back();
        } while (std::getline(stream, line.back()));
        size = line.size();
        mapping.resize(size);
        Utility::parallelize(
            size, [=, &separator, &comment, &line, &mapping](const unsigned int i) {
                mapping[i] = partition(trim(line[i], comment), separator);
                std::get<0>(mapping[i]) = trim(std::get<0>(mapping[i]), comment);
                std::get<1>(mapping[i]) = trim(std::get<1>(mapping[i]), comment);
                if ((std::get<0>(mapping[i]).empty()) ||
                    (std::get<1>(mapping[i]).empty()))
                    mapping[i] = Element();
            });
        mapping.erase(std::remove(mapping.begin(), mapping.end(), Element()),
                      mapping.end());
        mapping.erase(std::unique(mapping.begin(), mapping.end(),
                                  [](const Element &x, const Element &y) {
                                      return std::get<0>(x) == std::get<0>(y);
                                  }),
                      mapping.end());
        std::for_each(mapping.begin(), mapping.end(),
                      [=, &result](const Element &e) {
                          result[std::get<0>(e)] = std::get<1>(e);
                      });
        stream.close();
    }

    // Finalization
    return result;
}
// -------------------------------------------------------------------------- //

// ---------------------------------- CONES --------------------------------- //
// Save temporary cone file
/// \brief          Save temporary cone file.
/// \details        Saves the octree in a temporary cone file.
/// \tparam         Octree Octree type.
/// \param[in,out]  octree Source octree.
/// \param[in]      filename File name.
/// \return         True on success, false otherwise.
template <class Octree, class>
bool Input::save(Octree &octree, const std::string &filename) {
    std::ofstream stream;
    bool ok = false;
    if (!magrathea::FileSystem::exist(filename)) {
        stream.open(filename);
        if (stream) {
            ok = magrathea::DataHandler::rwrite(stream, octree.data(),
                                                octree.data() + octree.size());
            stream.close();
        }
    }
    return ok;
}

// Load temporary cone file
/// \brief          Load temporary cone file.
/// \details        Loads a temporary cone file into an octree.
/// \tparam         Octree Octree type.
/// \param[in,out]  octree Destination octree.
/// \param[in]      filename File name.
/// \return         True on success, false otherwise.
template <class Octree, class>
bool Input::load(Octree &octree, const std::string &filename) {
    static const unsigned long long int factor =
        sizeof(unsigned long long int) * sizeof(unsigned long long int);
    const unsigned long long int original = octree.size();
    std::ifstream stream(filename);
    long long int size = magrathea::FileSystem::size(stream);
    bool ok = (stream) && (size >= 0);
    if (ok) {
        // Need to know the size of octree beforehand if reading octree with
        // different structure
        octree.reserve(original + (size / SIZEOFOCTREE) +
                       ((size / SIZEOFOCTREE) / factor));
        octree.resize(original + size / SIZEOFOCTREE);
        ok = magrathea::DataHandler::rread(stream, octree.data() + original,
                                           octree.data() + octree.size());
        stream.close();
    }
    return ok;
}
// -------------------------------------------------------------------------- //

// ------------------------------- CORRECTION ------------------------------- //
// Cosmology correction
/// \brief          Cosmology correction.
/// \details        Produces another cosmology based on interpolation of an
///                 homogeneous trajectory.
/// \tparam         Cosmology Cosmology container type.
/// \tparam         Trajectory Trajectory container.
/// \tparam         Type Data type.
/// \param[in]      cosmology Cosmology.
/// \param[in]      trajectory Homogeneous trajectory.
/// \return         Corrected cosmology.
template <class Cosmology, class Trajectory, class Type, class>
Cosmology Input::correct(const Cosmology &cosmology,
                         const Trajectory &trajectory) {
    // Initialization
    static const unsigned int npasses = 4;
    static const Type one = Type(1);
    static const Type two = one + one;
    static const unsigned int nsmooth = 500;
    const unsigned int size = std::get<0>(cosmology).size();
    const unsigned int length = trajectory.size();
    const Type twob2 =
        two * std::pow(((length > 0) ? (trajectory.back().t()) : (one)) /
                           Type(std::kilo::num),
                       two);
    std::vector<Type> t(length);
    std::vector<Type> a(length);
    Cosmology result = cosmology;

    // Apply correction
    if (length > 0) {
        Utility::parallelize(length,
                             [=, &trajectory, &t, &a](const unsigned int i) {
                                 t[i] = trajectory[i].t();
                                 a[i] = trajectory[i].a();
                             });
        // Extend cosmology with FLRW trajectory
        Utility::parallelize(size, [=, &t, &a, &result](const unsigned int i) {
            if (!((std::get<0>(result)[i] < t.front()) ||
                  (std::get<0>(result)[i] > t.back())))
                std::get<1>(result)[i] =
                    Utility::interpolate(std::get<0>(result)[i], t, a);
        });
        if (length < size) {
            for (unsigned int ipass = 0; ipass < npasses; ++ipass) {
                std::get<1>(result) = Utility::smooth(
                    std::get<0>(result), std::get<1>(result),
                    [=, &twob2](Type xi, Type xj) {
                        return std::exp(-((xi - xj) * (xi - xj)) / twob2);
                    },
                    size / nsmooth);
                std::get<1>(result) = Utility::smooth(
                    std::get<0>(result), std::get<1>(result),
                    [=, &twob2](Type xi, Type xj) {
                        return std::exp(-((xi - xj) * (xi - xj) * std::deca::num) /
                                        twob2);
                    },
                    size / (nsmooth * std::deca::num));
                Utility::parallelize(size, [=, &t, &a, &result](const unsigned int i) {
                    if (!((std::get<0>(result)[i] < t.front()) ||
                          (std::get<0>(result)[i] > t.back())))
                        std::get<1>(result)[i] =
                            Utility::interpolate(std::get<0>(result)[i], t, a);
                });
            }
        }
    }

    // Finalization
    return result;
}

// Octree correction
/// \brief          Octree correction.
/// \details        Corrects uncomplete tree and zones with empty rho.
/// \tparam         Parameter Parameter type.
/// \tparam         Check Check mode.
/// \tparam         Selection Index of the data to be corrected.
/// \tparam         Octree Octree type.
/// \tparam         Kind Kind of amin.
/// \tparam         Dimension Number of dimensions.
/// \tparam         Element Underlying element type.
/// \tparam         Index Index type.
/// \tparam         Data Data type.
/// \tparam         Type Selected data type.
/// \tparam         Position Position of the hyperoctree center.
/// \tparam         Extent Extent of the hyperoctree.
/// \param[in]      parameters Parameter structure.
/// \param[in,out]  octree Input octree.
/// \param[out]     amin Outputs the lowest value of a.
/// \return         Reference to the octree.
template <class Parameter, int Check, unsigned int Selection, class Octree,
          typename Kind, unsigned int Dimension, class Element, class Index,
          class Data, class Type, class Position, class Extent, class>
Octree &Input::correct(const Parameter &parameters, Octree &octree,
                       Kind &&amin) {
    // Initialization
    static const unsigned int reservation = 1 << 20;
    static const int zero = 0;
    static const int one = 1;
    unsigned int size = octree.update().size();
    unsigned int ncoarse =
        (!octree.empty())
            ? (std::get<0>(
                   *std::min_element(octree.begin(), octree.end(),
                                     [](const Element &x, const Element &y) {
                                         return std::get<0>(x).level() <
                                                std::get<0>(y).level();
                                     }))
                   .level())
            : (zero);
    unsigned int nmax =
        (!octree.empty())
            ? (std::get<0>(
                   *std::max_element(octree.begin(), octree.end(),
                                     [](const Element &x, const Element &y) {
                                         return std::get<0>(x).level() <
                                                std::get<0>(y).level();
                                     }))
                   .level())
            : (zero);
    unsigned int counter = zero;
    int asize = zero;
    std::vector<unsigned int> count(size);
    std::vector<unsigned int> index;
    std::vector<Type> a;
    std::mutex mutex;
    std::atomic<unsigned int> distance(static_cast<unsigned int>(zero));
    // Correct coarse level
    if (Check >= zero) {
        if (parameters.coarsecorrection) {
            // Check if we find a cell at coarse level for which the density is not
            // normal
            Utility::parallelize(
                size, [=, &ncoarse, &count, &octree](const unsigned int i) {
                    count[i] = ((std::get<0>(octree[i]).level() == ncoarse) &&
                                (!std::isnormal(std::get<1>(octree[i]).rho())));
                });
            std::for_each(count.begin(), count.end(), [=, &counter](unsigned int &i) {
                i = (i > zero) ? (++counter) : (zero);
            });
            index.resize(counter);
            Utility::parallelize(size, [=, &count, &index](const unsigned int i) {
                if (count[i] > zero) {
                    index[count[i] - one] = i;
                }
            });
            // For non-normal coarse densities, average over normal neighbours
            Utility::parallelize(
                counter, [=, &ncoarse, &index, &octree](const unsigned int i) {
                    std::get<1>(octree[index[i]]).template data<Selection>() =
                        mean(octree, octree[index[i]], ncoarse);
                });
            Utility::parallelize(count.begin(), count.end(),
                                 [](unsigned int &i) { i = zero; });
        }
    }
    // Correct refined levels
    if (Check >= zero) {
        if (parameters.correction) {
            // Check if the value in .dphidt() is nonnormal(zero, NaN or infinity)
            for (unsigned int n = ncoarse + 1; n <= nmax; ++n) {
                Utility::parallelize(size, [=, &count, &octree](const unsigned int i) {
                    count[i] = ((std::get<0>(octree[i]).level() == n) &&
                                (!std::isnormal(std::get<1>(octree[i]).dphidt())));
                });
                // If dphidt() is nonnormal, give the .dphidt() from the parent
                Utility::parallelize(size, [=, &count, &octree](const unsigned int i) {
                    if (count[i] > zero) {
                        std::get<1>(octree[i]).dphidt() =
                            std::get<1>(*octree.find(std::get<0>(octree[i]).parent()))
                                .dphidt();
                    }
                });
            }
            // Check if the value in .rho() is nonnormal(zero, NaN or infinity)
            for (unsigned int n = ncoarse + 1; n <= nmax; ++n) {
                Utility::parallelize(size, [=, &count, &octree](const unsigned int i) {
                    count[i] = ((std::get<0>(octree[i]).level() == n) &&
                                (!std::isnormal(std::get<1>(octree[i]).rho())));
                });
                // If rho() is nonnormal, give the .rho() from the parent
                Utility::parallelize(size, [=, &count, &octree](const unsigned int i) {
                    if (count[i] > zero) {
                        std::get<1>(octree[i]).rho() =
                            std::get<1>(*octree.find(std::get<0>(octree[i]).parent()))
                                .rho();
                    }
                });
            }
        } else {
            Utility::parallelize(
                size, [=, &ncoarse, &count, &octree](const unsigned int i) {
                    if ((std::get<0>(octree[i]).level() > ncoarse) &&
                        (!std::isnormal(
                            std::get<1>(octree[i]).template data<Selection>()))) {
                        count[i] =
                            std::distance(octree.begin(),
                                          octree.find(std::get<0>(octree[i]).parent())) +
                            one;
                    }
                });
            count.erase(std::remove(count.begin(), count.end(), zero), count.end());
            Utility::parallelize(count.begin(), count.end(),
                                 [](unsigned int &i) { --i; });
            std::sort(count.begin(), count.end());
            count.erase(std::unique(count.begin(), count.end()), count.end());
            size = count.size();
            for (unsigned int i = 0; i < size; ++i) {
                if (!std::get<0>(octree[count[i]]).invalidated()) {
                    octree.coarsen(octree.begin() + count[i]);
                }
            }
            octree.update();
            size = octree.size();
            count.resize(size);
        }
        Utility::parallelize(count.begin(), count.end(),
                             [](unsigned int &i) { i = zero; });
    }
    // Detect non complete zones
    if (parameters.correction) {
        // Check if we have 8 sons per cell
        Utility::parallelize(size, [=, &count, &octree](const unsigned int i) {
            if (!octree.leaf(octree.begin() + i)) {
                for (unsigned int j = 0; j < Index::sites(); ++j) {
                    count[i] +=
                        (octree.find(std::get<0>(octree[i]).child(j)) == octree.end());
                }
            }
        });
        // Loop over ceclls in octree
        for (unsigned int i = 0; i < size; ++i) { // Parallelize ?
            if (count[i] > zero) {
                // Get index and data from cell
                const Index idx = std::get<0>(*(octree.begin() + i));
                const Data dat = std::get<1>(*(octree.begin() + i));
                Data datinterp;
                Index idxchild;
                // Loop over sons
                for (unsigned int isite = 0; isite < Index::sites(); ++isite) {
                    // Create index of son
                    idxchild = idx.child(isite);
                    // Interpolate at son location.
                    if (ORDER == 0)
                        datinterp =
                            octree.ngp(idxchild.template center<Type, Position, Extent>(0),
                                       idxchild.template center<Type, Position, Extent>(1),
                                       idxchild.template center<Type, Position, Extent>(2));
                    else if (ORDER == 1)
                        datinterp =
                            octree.cic(idxchild.template center<Type, Position, Extent>(0),
                                       idxchild.template center<Type, Position, Extent>(1),
                                       idxchild.template center<Type, Position, Extent>(2));
                    else if (ORDER == 2)
                        datinterp =
                            octree.tsc(idxchild.template center<Type, Position, Extent>(0),
                                       idxchild.template center<Type, Position, Extent>(1),
                                       idxchild.template center<Type, Position, Extent>(2));
                    else {
                        std::cout << "# WARNING : in Makefile, ORDER must be 0, 1 or 2"
                                  << std::endl;
                        std::cout << "# Error at file " << __FILE__
                                  << ", line : " << __LINE__ << std::endl;
                        std::terminate();
                    }
                    octree.append((datinterp != Data()) ? Element(idxchild, datinterp)
                                                        : Element(idxchild, dat));
                }
            }
        }
        // If we do not want to correct uncomplete refinement levels, then erase
        // them
    } else {
        Utility::parallelize(size, [=, &count, &octree](const unsigned int i) {
            if (!octree.leaf(octree.begin() + i)) {
                for (unsigned int j = 0; j < Index::sites(); ++j) {
                    count[i] +=
                        (octree.find(std::get<0>(octree[i]).child(j)) == octree.end());
                }
            }
        });
        for (unsigned int i = 0; i < size; ++i) {
            if (count[i] > zero) {
                if (!std::get<0>(octree[i]).invalidated()) {
                    octree.coarsen(octree.begin() + i);
                }
            }
        }
    }
    // Correct a (specific to Full Universe Run simulations)
    if (Check >= zero) {
        if (parameters.acorrection) {
            std::vector<Type> data;
            a.reserve(reservation);
            Utility::parallelize(octree.begin(), octree.end(),
                                 [=, &a, &distance, &mutex](const Element &e) {
                                     if (!std::count(a.begin(), a.begin() + distance,
                                                     std::get<1>(e).a())) {
                                         mutex.lock();
                                         a.push_back(std::get<1>(e).a());
                                         distance = a.size();
                                         mutex.unlock();
                                     }
                                 });
            a.erase(std::remove_if(a.begin(), a.end(),
                                   [](const Type x) {
                                       return std::signbit(x) || !std::isnormal(x);
                                   }),
                    a.end());
            std::sort(a.begin(), a.end());
            a.erase(std::unique(a.begin(), a.end()), a.end());
            asize = a.size();
            data.resize(asize);
            std::iota(data.begin(), data.end(), one);
            if (asize > one) {
                data.insert(data.begin(),
                            a[0] - Utility::differentiate<1>(Type(), data, a, one));
                if (!std::signbit(data[zero]) && std::isnormal(data[zero]) &&
                    data[zero] < a[zero]) {
                    a.insert(a.begin(), data[zero]);
                    asize = a.size();
                }
            }
            if (asize > zero) {
                Utility::parallelize(octree.begin(), octree.end(), [=, &a](Element &e) {
                    std::get<1>(e).a() =
                        *(std::lower_bound(a.begin(), a.end(), std::get<1>(e).a()) - 1);
                });
                amin = a[zero];
            }
        }
    }
    // Finalization
    octree.update();
#ifdef VERBOSE
    std::cout << "# After correction, Octree size : " << octree.size()
              << " Octree capacity " << octree.capacity() << std::endl;
#endif
    return octree;
}

// -------------------------------------------------------------------------- //

// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Input.
/// \return         0 if no error.
/*int Input::example()
{
    // Initialize
    std::cout<<"BEGIN = Input::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    std::array<double, 3> first = std::array<double, 3>({{4., 8., 15.}});
    std::array<double, 3> second = std::array<double, 3>({{16., 23., 42.}});
    Gravity<float, 3> gravity;
    magrathea::SimpleHyperOctree<double,
magrathea::SimpleHyperOctreeIndex<unsigned long long int, 3>, unsigned int>
counter(0, 1); magrathea::SimpleHyperOctree<double,
magrathea::SimpleHyperOctreeIndex<unsigned long long int, 3>, std::string>
ftree(0, 2); magrathea::SimpleHyperOctree<double,
magrathea::SimpleHyperOctreeIndex<unsigned long long int, 3>, Gravity<float, 3>
> octree(0, 2); magrathea::HyperSphere<3> sphere =
magrathea::HyperSphere<3>::unit(); Cone<> cone(first, second, 0.42);
    std::vector<std::string> list;
    std::array<std::vector<double>, 4> cosmology;
    std::vector<Photon<double, 3> > trajectory;
    std::string string;

    // Construction
    Input input;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : " <<std::endl;
    std::cout<<std::setw(width)<<"Input() : "                                ;
Input(); std::cout<<std::endl; std::cout<<std::setw(width)<<"input = Input() : "
; input = Input(); std::cout<<std::endl;

    // Utilities
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Utilities : " <<std::endl;
    std::cout<<std::setw(width*3)<<"input.trim(string) : "
<<input.trim(string)<<std::endl;
    std::cout<<std::setw(width*3)<<"input.partition(string).first : "
<<input.partition(string).first<<std::endl;
    std::cout<<std::setw(width*3)<<"input.count(counter, ftree) : "
<<input.count(counter, ftree)<<std::endl;
    std::cout<<std::setw(width*3)<<"input.collide(octree,
std::get<0>(octree[0]), sphere, cone) : " <<input.collide(octree,
std::get<0>(octree[0]), sphere, cone)<<std::endl;
    std::cout<<std::setw(width*3)<<"input.mean(octree, octree[1]) : "
<<input.mean(octree, octree[1])<<std::endl;
    std::cout<<std::setw(width*3)<<"input.constantify(100, 0., 42., 1., 0.,
0.).size() : " <<input.constantify(100, 0., 42., 1., 0., 0.).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"input.sistemize(std::get<1>(octree[0]), 0.5,
0.72, 0.3, 21000.) : " <<input.sistemize(std::get<1>(octree[0]), 0.5, 0.72, 0.3,
21000.)<<std::endl; std::cout<<std::setw(width*3)<<"input.sistemize(octree,
0.72, 0.3, 21000.) : " <<input.sistemize(octree, 0.72, 0.3, 21000.)<<std::endl;
    std::cout<<std::setw(width*3)<<"input.homogenize(std::get<1>(octree[0])) : "
<<input.homogenize(std::get<1>(octree[0]))<<std::endl;
    std::cout<<std::setw(width*3)<<"input.homogenize(octree) : "
<<input.homogenize(octree)<<std::endl;
    std::cout<<std::setw(width*3)<<"input.schwarzschildify(std::get<1>(octree[0]),
first, second, 42, 1) : " <<input.schwarzschildify(std::get<1>(octree[0]),
first, second, 42, 1)<<std::endl;
    std::cout<<std::setw(width*3)<<"input.schwarzschildify(octree, first, 42.,
21000., [](Gravity<float, 3>&, unsigned int){return false;}) : "
<<input.schwarzschildify(octree, first, 42., 21000., [](Gravity<float, 3>&,
unsigned int){return false;})<<std::endl;

    // Files
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Files : " <<std::endl;
    std::cout<<std::setw(width*2)<<"input.filetree(ftree, \"/tmp/\",
\"file_%05d\") : "                                 <<input.filetree(ftree,
"/tmp/", "file_%05d")<<std::endl;
    std::cout<<std::setw(width*2)<<"input.prepare(list, ftree, sphere, cone) : "
<<input.prepare(list, ftree, sphere, cone)<<std::endl;

    // Data
    std::cout<<std::endl;
    std::cout<<std::setw(width*3)<<"Data : " <<std::endl;
    std::cout<<std::setw(width*3)<<"input.import(octree, \"/tmp/file_00001\",
13, 0, [](const decltype(octree.element())&){return true;}) : "
<<input.import(octree, "/tmp/file_00001", 13, 0, [](const
decltype(octree.element())&){return true;})<<std::endl;
    std::cout<<std::setw(width*3)<<"input.acquire(string, string, string,
first[0], first[1], first[2]).size() : " <<input.acquire(string, string, string,
first[0], first[1], first[2]).size()<<std::endl;
    std::cout<<std::setw(width*3)<<"input.parse(string, \"=\", \"#\").size() : "
<<input.parse(string, "=", "#").size()<<std::endl;

    // Cones
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Cones : " <<std::endl;
    std::cout<<std::setw(width)<<"input.save(counter, \"/tmp/file_00000\") :
"<<input.save(counter, "/tmp/file_00000")<<std::endl;
    std::cout<<std::setw(width)<<"input.load(counter, \"/tmp/file_00000\") :
"<<input.load(counter, "/tmp/file_00000")<<std::endl;

    // Correction
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Correction : " <<std::endl;
    std::cout<<std::setw(width*2)<<"input.correct(cosmology, trajectory).size()
: "                                     <<input.correct(cosmology,
trajectory).size()<<std::endl;
    std::cout<<std::setw(width*2)<<"input.correct<-1>(octree).size() : "
<<input.correct<-1>(octree).size()<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Input::example()"<<std::endl;
    return 0;
}*/
// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
#endif // INPUT_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
