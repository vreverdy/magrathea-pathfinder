/* ********************************** MISCELLANEOUS ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Miscellaneous
// DESCRIPTION :    Some miscellaneous function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           miscellaneous.h
/// \brief          Some miscellaneous function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef MISCELLANEOUS_H_INCLUDED
#define MISCELLANEOUS_H_INCLUDED

#include <ctime>
// Include C++
#include <algorithm>
#include <array>
#include <atomic>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#ifdef GCCBELOW7
#include <experimental/algorithm>
#endif
// Include libs
#include <mpi.h>
// Include project
#include "TReadHDF5.h"
#include "cone.h"
#include "input.h"
#include "magrathea/abstracthypersphere.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/hypercube.h"
#include "magrathea/hypersphere.h"
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/timer.h"
#include "output.h"
#include "utility.h"

using namespace magrathea;

using integer = int;
using uint = unsigned int;
using real = double;
using floating = float;
using point = std::array<real, 3>;
static constexpr uint zero = 0;
static constexpr uint one = 1;

class Miscellaneous {

    // Methodes
public:
    // Name of files in Directory
    static void getFilesinDir(const std::string dirName, std::vector<std::string> &fileNames);
    // Tokenize string
    static void Tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters = " ");

    // Clear and shrink
    template <class Vectored>
    static void clear_shrink(Vectored &vector);
    template <typename Type>
    static void fullclear_vector(std::vector<Type> &vector);
    // Ticket
    template <class Function, typename Integer>
    static void TicketizeFunction(const Integer rank, const Integer ntasks, Function &&function);
    // Get specs
    template <class Parameter, typename Scalar>
    static void get_narrow_specs(const Parameter &parameters, std::array<std::array<double, 3>, 3> &rotm1, Scalar &thetay, Scalar &thetaz);
    // Load & Correct octree
    template <class Octree, class Filelist, typename Integer>
    static void loadOctree(const Integer icone, Octree &octree, Filelist &conefile);
    template <class Octree, class Cosmology, class Parameters, typename Real>
    static void correctOctree(Octree &octree, const Cosmology &cosmology, Parameters &parameters, const Real h, const Real omegam, const Real lboxmpch, Real &amin);
    // Targets
    template <class Vector, typename Scalar, class Container>
    static std::vector<std::array<double, 8>> getTargets(const std::vector<std::array<double, 8>> &posTargets, const Cone<Vector, Scalar> &cone, const Container &cones);

    // Fill particles
    template <class Parameter, class Cone, typename Type1>
    static void fill_particles_vectors(const Parameter &parameters, const Cone &cone, const std::vector<std::string> &shellList, std::vector<Type1> &pos_part, std::vector<Type1> &force_part, std::vector<Type1> &potential_part, std::vector<Type1> &a_part, const double thetay, const double thetaz);

    // Visualise octree
    template <template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
    static void VizualizeOctree(const Octree<Type, Index, Data, Dimension, Position, Extent, Element, Container> &octree, const Type radius);

    // Read Angular position from previously computed catalog
    template <typename Integer, class Parameters>
    static void ReadFromCat(const Integer icone, const Parameters &parameters, std::vector<std::array<double, 18>> &catalogue);

    // Write and read cone orientation file
    template <class Cone, class Parameters>
    static void write_cone_orientation(const Cone &cones, const Cone &conesIfRot, const Parameters &parameters);
    template <class Cone, class Parameters>
    static void read_cone_orientation(Cone &cones, Cone &conesIfRot, const Parameters &parameters);
};

/// \brief          Get list of files.
/// \details        Get list of files and directories.
/// \param[in]      String Directory name.
/// \param[in,out]  vector<string> List of files and directories.
void Miscellaneous::getFilesinDir(const std::string dirName, std::vector<std::string> &fileNames) {

    DIR *pdir;
    struct dirent *pent;

    pdir = opendir(dirName.c_str()); //"." refers to the current dir
    if (!pdir) {
        std::cout << "opendir() failure; terminating" << std::endl;
        std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__ << " " << dirName << std::endl;
        std::terminate();
    }
    errno = 0;
    while ((pent = readdir(pdir))) {
        fileNames.push_back(pent->d_name);
    }
    if (errno) {
        std::cout << "readdir() failure; terminating" << std::endl;
        std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__ << " " << dirName << std::endl;
        std::terminate();
    }
    closedir(pdir);
}

// Clear vector
/// \brief          Tokenize string.
/// \details        Tokenize strings with delimiters.
/// \tparam         str String to tokenize
/// \param[in,out]  tokens Vector of strings.
/// \param[in]      delimiters Delimiters used to tokenize String.
/// \return         tokenized vector
void Miscellaneous::Tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters) {
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

// Clear vector
/// \brief          Clear and shrink vector.
/// \details        Clear and shrink vector.
/// \tparam         Vectored vector type
/// \param[in,out]  vector Vector to be cleared.
/// \return         Empty vector
template <class Vectored>
void Miscellaneous::clear_shrink(Vectored &vector) {

    vector.clear();
    vector.shrink_to_fit();
}
// Clear vector
/// \brief          Erase vector.
/// \details        Erase vector.
/// \tparam         Type Type type
/// \param[in,out]  vector Vector to be erased.
/// \return         Empty vector
template <typename Type>
void Miscellaneous::fullclear_vector(std::vector<Type> &vector) {
    std::vector<Type>().swap(vector);
}

// I/O Ticketing fuction
/// \brief          IO ticket funtion.
/// \details        IO ticket system for heavy IO computations.
/// \tparam         Function lambda function
/// \tparam         Integer integer type
/// \param[in]      rank Process rank .
/// \param[in]      ntasks Number of tasks
/// \param[in]      function Lambda function
/// \return         Operation using ticketing system.
template <class Function, typename Integer>
void Miscellaneous::TicketizeFunction(const Integer rank, const Integer ntasks, Function &&function) {

    Integer ioticket = 1;
    // If we use the ticket system
    if (IOGROUPSIZE > 1) {

        if (rank % IOGROUPSIZE > 0)
            MPI_Recv((void *)&ioticket, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        function();

        if ((rank + 1) % IOGROUPSIZE > 0 && rank + 1 < ntasks)
            MPI_Send((void *)&ioticket, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        // No ticket system
    } else {
        function();
    }
}

// Get narrow cone specs
/// \brief          Get narrow cone specs.
/// \details        Get narrow cone specs.
/// \tparam         Parameter Parameter type.
/// \tparam         Scalar scalar type
/// \param[in]      parameters Parameter structure
/// \param[in,out]  rotm1 Rotation matrix for narrow cone cells
/// \param[in,out]  thetay Semi-angle for solid angle in direction y
/// \param[in,out]  thetaz Semi-angle for solid angle in direction z
template <class Parameter, typename Scalar>
void Miscellaneous::get_narrow_specs(const Parameter &parameters, std::array<std::array<double, 3>, 3> &rotm1, Scalar &thetay, Scalar &thetaz) {

    std::size_t found;
    double theta_rot(0), phi_rot(0);
    std::array<std::array<double, 3>, 3> rotation = {{0}};
    std::vector<std::string> filelistingprior;
    std::string filelisting;

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.celldir, filelistingprior);
    for (uint ifiling = 0; ifiling < filelistingprior.size(); ++ifiling) {
        if (parameters.typefile == 1) {
            // Find any hdf5 file
            found = filelistingprior[ifiling].find(".h5");
            if (found != std::string::npos) {
                filelisting = parameters.celldir + filelistingprior[ifiling];
                break;
            }
        } else if (parameters.typefile == 2) {
            // For ASCII, search for the 'info_narrow_cone.txt' file
            if (filelistingprior[ifiling] == "info_narrow_cone.txt") {
                filelisting = parameters.celldir + filelistingprior[ifiling];
                break;
            }
        } else {
            std::cout << "# WARNING : Narrow cones can only be computed using HDF5 or ASCII input files" << std::endl;
            std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__ << std::endl;
            std::terminate();
        }
    }
    // Get informations from HDF5 files
    if (parameters.typefile == 1) {
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "phi", phi_rot);
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "theta", theta_rot);
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "thetay", thetay);
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "thetaz", thetaz);
        // Get informations from ASCII files
    } else if (parameters.typefile == 2) {
        std::map<std::string, std::string> parameterASCII;
        parameterASCII = Input::parse(filelisting);
        phi_rot = std::stoul(parameterASCII["phi"]);
        theta_rot = std::stod(parameterASCII["theta"]);
        thetay = std::stod(parameterASCII["thetay"]);
        thetaz = std::stod(parameterASCII["thetaz"]);
    } else {
        std::cout << "# WARNING : Narrow cones can only be computed using HDF5 or ASCII input files" << std::endl;
        std::cout << "# Error at file " << __FILE__ << ", line : " << __LINE__ << std::endl;
        std::terminate();
    }
    thetaz *= magrathea::Constants<double>::deg();
    thetay *= magrathea::Constants<double>::deg();
    phi_rot *= magrathea::Constants<double>::deg();
    theta_rot *= magrathea::Constants<double>::deg();

    // Compute the rotation matrix to rotate the cone
    rotation[0][0] = std::cos(theta_rot) * std::cos(phi_rot);
    rotation[0][1] = std::cos(theta_rot) * std::sin(phi_rot);
    rotation[0][2] = -std::sin(theta_rot);
    rotation[1][0] = -std::sin(phi_rot);
    rotation[1][1] = std::cos(phi_rot);
    rotation[1][2] = 0;
    rotation[2][0] = std::cos(phi_rot) * std::sin(theta_rot);
    rotation[2][1] = std::sin(phi_rot) * std::sin(theta_rot);
    rotation[2][2] = std::cos(theta_rot);
    // Inverse matrix
    rotm1 = Utility::invMatrix3d(rotation);
}

// Load  Octree
/// \brief          Load octree.
/// \details        Load binary octree from preparation phase.
/// \tparam         Octree octree type
/// \tparam         Filelist cone file type
/// \tparam         Integer integer type
/// \param[in]      icone cone number
/// \param[in,out]  octree Octree to be filled
/// \param[in]      conefile Cone names in conedir
template <class Octree, class Filelist, typename Integer>
void Miscellaneous::loadOctree(const Integer icone, Octree &octree, Filelist &conefile) {

    octree.fullclear();
#ifdef VERBOSE
    std::cout << "# Loading conefile[" << icone << "] " << conefile[icone] << std::endl;
#endif
    // Load octree
    Input::load(octree, conefile[icone]);
#ifdef VERBOSE
    std::cout << "# Conefile[" << icone << "] loaded. Size : " << octree.size() << std::endl;
#endif
#ifdef VELOCITYFIELD
    // When putting an octree generated with gravity.h in an octree with gravity2.h,
    // need to correct the position of data
    Utility::parallelize(octree.size(), [=, &octree](const uint i) {
        std::get<1>(octree[i]).rho() = std::get<1>(octree[i]).dphidy();
        std::get<1>(octree[i]).phi() = std::get<1>(octree[i]).dphidx();
        std::get<1>(octree[i]).dphidx() = std::get<1>(octree[i]).vz();
        std::get<1>(octree[i]).dphidy() = std::get<1>(octree[i]).dphidt();
        std::get<1>(octree[i]).dphidz() = std::get<1>(octree[i]).a();
        std::get<1>(octree[i]).a() = std::get<1>(octree[i]).vy();
        std::get<1>(octree[i]).dphidt() = std::get<1>(octree[i]).vx();
        std::get<1>(octree[i]).vxyz() = std::array<float, 3>();
    });
#endif
}

// Correct Octree
/// \brief          Correct octree.
/// \details        Process a loaded Octree.
/// \tparam         Octree octree type
/// \tparam         Cosmology cosmology type
/// \tparam         Parameters parameters type
/// \tparam         Real float/real type
/// \param[in,out]  octree Octree to be filled
/// \param[in]      cosmology Cosmological tables
/// \param[in]      conefile Cone names in conedir
/// \param[in]      parameters Parameters structure
/// \param[in]      h dimensionless Hubble parameter
/// \param[in]      omegam Matter density fraction
/// \param[in]      lboxmpch Size of simulation box
/// \param[in,out]  amin minimum value scale factor
template <class Octree, class Cosmology, class Parameters, typename Real>
void Miscellaneous::correctOctree(Octree &octree, const Cosmology &cosmology, Parameters &parameters, const Real h, const Real omegam, const Real lboxmpch, Real &amin) {

    // Convert from Ramses Units to SI
    Input::sistemize(parameters, octree, h, omegam, lboxmpch);
    // Initialise dphida to zero
    Utility::parallelize(octree.size(), [=, &octree](const uint i) { std::get<1>(octree[i]).dphidt() = 0; });
    // Apply correction to the octree,
    // also .update() is used twice (before and after corrections)
    Input::correct(parameters, octree, amin);
    // Convert from dphi/da to dphi/dt
    Utility::parallelize(octree.size(), [=, &octree](const uint i) { std::get<1>(octree[i]).dphidt() *= Utility::rinterpolate(std::get<1>(octree[i]).a(), std::get<1>(cosmology), std::get<2>(cosmology)); });
}

// Vizualize Octree
/// \brief          Write position of cells.
/// \details        Write position of cells under a given radius.
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in]      octree Octree to be filled
/// \param[in]      radius Maximum radius at which we write cells (in Ramses Units)
template <template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container>
void Miscellaneous::VizualizeOctree(const Octree<Type, Index, Data, Dimension, Position, Extent, Element, Container> &octree, const Type radius) {

    for (uint i = 0; i < octree.size(); i++) {
        double x = std::get<0>(octree[i]).template center<Type, Position, Extent>(0);
        double y = std::get<0>(octree[i]).template center<Type, Position, Extent>(1);
        double z = std::get<0>(octree[i]).template center<Type, Position, Extent>(2);
        if (x * x + y * y + z * z < radius * radius) {
            // Output result on terminal
            std::cout << x << " " << y << " " << z << " " << std::get<1>(octree[i]) << " #vizualize" << std::endl;
        }
    }
}

// ---------------------------------- TARGETS ---------------------------------- //

/// \brief          Get target.
/// \details        Get target if inside cone.
/// \tparam         Vector Position vector type.
/// \tparam         Scalar Scalar data type
/// \tparam         Container Container of cone type
/// \param[in]      posTargets vector of halo positions.
/// \param[in]      Cone current cone.
/// \param[in]	    Cones to compare distances.
/// return  	    Vector vector of targets inside cone
template <class Vector, typename Scalar, class Container>
std::vector<std::array<double, 8>> Miscellaneous::getTargets(const std::vector<std::array<double, 8>> &posTargets, const Cone<Vector, Scalar> &cone, const Container &cones) {
    // Fill selection vector with -1
    std::vector<int> selection(posTargets.size(), -1);

    // Loop over all targets
    Utility::parallelize(posTargets.size(), [=, &posTargets, &cone, &cones, &selection](const unsigned int ivec) {
        std::array<double, 3> position;
        double reference(0), length(0), distance(0);
        bool ok = false;
        position[0] = posTargets[ivec][0];
        position[1] = posTargets[ivec][1];
        position[2] = posTargets[ivec][2];
        reference = 0;
        length = 0;
        // Check if inside cone
        if (cone.inside(position)) {
            ok = true;
            // Compute scalar product of cone base and target direction
            for (unsigned int idim = 0; idim < 3; ++idim) {
                length += (cone.base(idim) - cone.vertex(idim)) * (position[idim] - cone.vertex(idim));
            }
            length /= cone.template pow<2>(cone.length());
            // Compute the ditance between cone base and target
            for (unsigned int idim = 0; idim < 3; ++idim) {
                reference += cone.template pow<2>(position[idim] - (cone.vertex(idim) + (cone.base(idim) - cone.vertex(idim)) * length));
            }
            // Loop over all the other cones
            for (unsigned int icone = 0; icone < cones.size(); ++icone) {
                length = 0;
                distance = 0;
                if (cones[icone] != cone) {
                    // Compute scalar product of other cone base and target direction
                    for (unsigned int idim = 0; idim < 3; ++idim) {
                        length += (cones[icone].base(idim) - cone.vertex(idim)) * (position[idim] - cones[icone].vertex(idim));
                    }
                    if (!(length < 0)) {
                        length /= cones[icone].template pow<2>(cones[icone].length());
                        // Compute the ditance between other cone base and target
                        for (unsigned int idim = 0; idim < 3; ++idim) {
                            distance += cones[icone].template pow<2>(position[idim] - (cones[icone].vertex(idim) + (cones[icone].base(idim) - cones[icone].vertex(idim)) * length));
                        }
                        // Check if other cone is closer to pixel
                        if (distance < reference) {
                            ok = false;
                            icone = cones.size();
                        }
                    } //  if length
                }     // if
            }         // for icone
            if (ok) {
                selection[ivec] = ivec;
            }
        } //  if
    });   //  while

    // Erase targets which are not inside the cone
    selection.erase(std::remove(std::begin(selection), std::end(selection), -1), std::end(selection));
    std::vector<std::array<double, 8>> pointsCible(selection.size());
    // Put targets in vector
    Utility::parallelize(selection.size(), [=, &posTargets, &pointsCible, &selection](const unsigned int ivec) {
        pointsCible[ivec] = posTargets[selection[ivec]];
    });
    return pointsCible;
}

// Fill vectors of particles
/// \brief          Rewrite particle in vector.
/// \details        Rewrite particle position/velocity from HDF5 in a vector.
/// \tparam         Parameter Parameter type
/// \tparam         Cone cone type
/// \tparam         Type1 type1 type
/// \param[in]      parameters Parameter structure
/// \param[in]      cone Cone parameters
/// \param[in]      shellList List of HDF5 files containing particles
/// \param[in,out]  pos_part Position of particles
/// \param[in,out]  force_part Force of particles
/// \param[in,out]  potential_part Potential of particles
/// \param[in,out]  a_part scale factor of particles
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
template <class Parameter, class Cone, typename Type1>
void Miscellaneous::fill_particles_vectors(const Parameter &parameters, const Cone &cone, const std::vector<std::string> &shellList, std::vector<Type1> &pos_part, std::vector<Type1> &force_part, std::vector<Type1> &potential_part, std::vector<Type1> &a_part, const double thetay, const double thetaz) {

    unsigned long long int marker1(0), marker2(0);

#ifdef VERBOSE
    std::cout << "# Now loading the position, force and potential of particles" << std::endl;
#endif
    Miscellaneous::fullclear_vector(pos_part);
    Miscellaneous::fullclear_vector(force_part);
    Miscellaneous::fullclear_vector(potential_part);
    // Loop over files
    for (uint ifiling = 0; ifiling < shellList.size(); ifiling++) {
        double aexp(0);
        // Get scale factor of shell
        TReadHDF5::getAttribute(shellList[ifiling], "metadata/cone_info", "aexp", aexp);
        // Get data from file, depending on whether it intersects the cone or not
        TReadHDF5::fillVectors_part(parameters, shellList[ifiling], thetay, thetaz, cone, "data", "position_part", pos_part, "gravitational_field_part", force_part, "potential_part", potential_part);
        float unit_l(0), unit_t(0);
        // Get conversion factor from Ramses Units to SI
        TReadHDF5::getAttribute(shellList[ifiling], "metadata/ramses_info", "unit_l", unit_l); // comobile
        TReadHDF5::getAttribute(shellList[ifiling], "metadata/ramses_info", "unit_t", unit_t); // superconformal unit
        // Factors to convert Potential and Force in SI
        const double factorpot = std::pow(unit_l * 1e-2 / unit_t, 2);
        const double factorforce = -aexp * unit_l * 1e-2 / (unit_t * unit_t);
        std::transform(std::begin(potential_part) + marker1, std::end(potential_part), std::begin(potential_part) + marker1, std::bind1st(std::multiplies<double>(), factorpot)); // SI units
        std::transform(std::begin(force_part) + marker2, std::end(force_part), std::begin(force_part) + marker2, std::bind1st(std::multiplies<double>(), factorforce));           // SI units
        std::vector<Type1> a_tmp(potential_part.size() - marker1);
        std::fill(a_tmp.begin(), a_tmp.end(), aexp);
        // Give the same scale factor value to all the particles in the same shell
        a_part.insert(std::end(a_part), std::begin(a_tmp), std::end(a_tmp));
        marker1 = potential_part.size();
        marker2 = force_part.size();
    }
}

// Read source properties from pre-computed catalogue
/// \brief          Read magrathea source catalogue.
/// \details        Read magrathea source catalogue.
/// \tparam         Integer Integer type
/// \tparam         Parameter Parameter type
/// \param[in]      parameters Parameter structure
/// \param[in]      icone cone number
/// \param[in,out]  Catalogue Vector containing the source catalogue
template <typename Integer, class Parameters>
void Miscellaneous::ReadFromCat(const Integer icone, const Parameters &parameters, std::vector<std::array<double, 18>> &catalogue) {
    // Create filename of catalogue
    std::string filename = parameters.outputdir + "/" + Output::name(parameters.base, "_", std::make_pair("%05d", icone), ".txt"); // Name of catalog, given icone, directory and base
#ifdef VERBOSE
    std::cout << "# Cone " << icone << " Read angular position from " << filename << std::endl;
#endif
    // Open filename
    std::ifstream streaming(filename.c_str());
    streaming.unsetf(std::ios_base::skipws);
    uint size = std::count(std::istream_iterator<char>(streaming), std::istream_iterator<char>(), '\n');
    streaming.close();
    std::ifstream stream(filename.c_str());
    catalogue.resize(size);
    // Read halo catalogue or particle catalogue (for the latter there is no 'npart' column)

    for (unsigned int i = 0; i < size; ++i) {
        stream >> catalogue[i][0] >> catalogue[i][1] >> catalogue[i][2] >> catalogue[i][3] >> catalogue[i][4] >> catalogue[i][5] >> catalogue[i][6] >> catalogue[i][7] >> catalogue[i][8] >> catalogue[i][9] >> catalogue[i][10] >> catalogue[i][11] >> catalogue[i][12] >> catalogue[i][13] >> catalogue[i][14] >> catalogue[i][15] >> catalogue[i][16] >> catalogue[i][17];
    }
}

// Write cone properties in ascii file
/// \brief          Write cone properties
/// \details        Write cone properties
/// \tparam         Cone cones type
/// \tparam         Parameter Parameter type
/// \param[in]      cones cone container
/// \param[in]      conesIfRot rotated cone container
/// \param[in]      parameters Parameter structure
template <class Cone, class Parameters>
void Miscellaneous::write_cone_orientation(const Cone &cones, const Cone &conesIfRot, const Parameters &parameters) {
    // Create filename of catalogue
    std::string filename, filename2;
    if (parameters.isfullsky) {
        filename = Output::name(parameters.conedir, "/cone_orientations_ncones_", std::make_pair("%05d", parameters.ncones), ".txt");
    } else {
        filename = Output::name(parameters.conedir, "/cone_orientations_ncones_", std::make_pair("%05d", parameters.ncones), "_buffer_", std::make_pair("%5.4f", parameters.buffer), ".txt");
        filename2 = Output::name(parameters.conedir, "/cone_orientations_rotated_ncones_", std::make_pair("%05d", parameters.ncones), "_buffer_", std::make_pair("%5.4f", parameters.buffer), ".txt");
    }
#ifdef VERBOSE
    std::cout << "# Writing in " << filename << std::endl;
#endif
    // Open filename
    std::ofstream streaming;
    streaming.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
    streaming.close();
    std::ofstream stream(filename.c_str(), std::ios::app);

    for (uint i = 0; i < cones.size(); i++) {
        stream << cones[i] << std::endl;
    }

    stream.close();

    if (!parameters.isfullsky) {
        std::ofstream streaming2;
        streaming2.open(filename2.c_str(), std::ofstream::out | std::ofstream::trunc);
        streaming2.close();
        std::ofstream stream2(filename2.c_str(), std::ios::app);

        for (uint i = 0; i < conesIfRot.size(); i++) {
            stream2 << conesIfRot[i] << std::endl;
        }
        stream2.close();
    }
}

// Read cone properties from ascii file
/// \brief          Read cone properties
/// \details        Read cone properties
/// \tparam         Cone cones type
/// \tparam         Parameter Parameter type
/// \param[in]      cones cone container
/// \param[in]      conesIfRot rotated cone container
/// \param[in]      parameters Parameter structure
template <class Cone, class Parameters>
void Miscellaneous::read_cone_orientation(Cone &cones, Cone &conesIfRot, const Parameters &parameters) {
    // Create filename of catalogue
    std::string filename, filename2;
    if (parameters.isfullsky) {
        filename = Output::name(parameters.conedir, "/cone_orientations_ncones_", std::make_pair("%05d", parameters.ncones), ".txt");
    } else {
        filename = Output::name(parameters.conedir, "/cone_orientations_ncones_", std::make_pair("%05d", parameters.ncones), "_buffer_", std::make_pair("%5.4f", parameters.buffer), ".txt");
        filename2 = Output::name(parameters.conedir, "/cone_orientations_rotated_ncones_", std::make_pair("%05d", parameters.ncones), "_buffer_", std::make_pair("%5.4f", parameters.buffer), ".txt");
    }
#ifdef VERBOSE
    std::cout << "# Reading " << filename << std::endl;
#endif
    // Open filename
    std::ifstream streaming(filename.c_str());
    streaming.unsetf(std::ios_base::skipws);
    uint size = std::count(std::istream_iterator<char>(streaming), std::istream_iterator<char>(), '\n');
    streaming.close();
    std::ifstream stream(filename.c_str());

    for (uint i = 0; i < size; ++i) {
        stream >> cones[i].vertex(0) >> cones[i].vertex(1) >> cones[i].vertex(2) >> cones[i].base(0) >> cones[i].base(1) >> cones[i].base(2) >> cones[i].angle();
    }

    stream.close();

    if (parameters.isfullsky) {
        conesIfRot = cones;
    } else {

        std::ifstream streaming2(filename2.c_str());
        streaming2.unsetf(std::ios_base::skipws);
        size = std::count(std::istream_iterator<char>(streaming2), std::istream_iterator<char>(), '\n');
        streaming2.close();
        std::ifstream stream2(filename2.c_str());

        for (uint i = 0; i < size; ++i) {
            stream2 >> conesIfRot[i].vertex(0) >> conesIfRot[i].vertex(1) >> conesIfRot[i].vertex(2) >> conesIfRot[i].base(0) >> conesIfRot[i].base(1) >> conesIfRot[i].base(2) >> conesIfRot[i].angle();
        }
        stream2.close();
    }
}

#endif // MISCELLANEOUS_H_INCLUDED
