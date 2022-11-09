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
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <utility>
#include <random>
#include <string>
#include <vector>
#include <deque>
#include <array>
#include <tuple>
#include <atomic>
#include <mutex>

#ifdef GCCBELOW7
#include <experimental/algorithm>
#endif
// Include libs
#include <mpi.h>
// Include project
#include "magrathea/simplehyperoctree.h"
#include "magrathea/simplehyperoctreeindex.h"
#include "magrathea/hypersphere.h"
#include "magrathea/abstracthypersphere.h"
#include "magrathea/hypercube.h"
#include "magrathea/constants.h"
#include "magrathea/evolution.h"
#include "magrathea/timer.h"
#include "cone.h"
#include "utility.h"
#include "input.h"
#include "output.h"
#include "TReadHDF5.h"


using namespace magrathea;

    using integer = int;
    using uint = unsigned int;
    using real = double;
    using floating = float;
    using point = std::array<real, 3>;
    static constexpr uint zero = 0;
    static constexpr uint one = 1;


class Miscellaneous {

	//Methodes
	public:
	
	// Name of files in Directory
	static void getFilesinDir(const std::string dirName, std::vector<std::string>& fileNames);
	// Tokenize string
	static void Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ");

	// Clear and shrink
	template <class Vectored> static void clear_shrink(Vectored& vector);
	template <typename Type> static void fullclear_vector(std::vector<Type>& vector);
	// Ticket
	template <class Function, typename Integer> static void TicketizeFunction(const Integer rank, const Integer ntasks, Function&& function);
	// Cones generation
	template <  class Cone, template <unsigned int, class, typename > class Sphere, unsigned int Dimension, class Vector, typename Scalar, typename Integer > static void GenerateFullskyCones(const Integer ncones, std::vector< Cone >& cone, std::vector< Cone >& coneIfRot, Sphere<Dimension, Vector, Scalar>& sphere);
	template <  class Parameter, class Cone, template <unsigned int, class, typename > class Sphere, unsigned int Dimension, class Vector, typename Scalar > static void GenerateNarrowCones(const Parameter& parameters, std::vector< Cone >& cone, std::vector< Cone >& coneIfRot, Sphere<Dimension, Vector, Scalar>& sphere,  std::array< std::array< double, 3 >, 3 >& rotm1, Scalar& thetay, Scalar& thetaz);
	// Load & Correct octree
	template < class Octree, class Filelist,  typename Integer > static void loadOctree(const Integer icone, Octree& octree, Filelist& conefile);
	template < class Octree, class Cosmology, class Parameters, typename Real > static void correctOctree(Octree& octree, const Cosmology& cosmology, Parameters& parameters, const Real h, const Real omegam, const Real lboxmpch, Real& amin);
	// Targets
	template < class Vector, typename Scalar,  class Container>  static std::vector< std::array<double,8> > getTargets(const std::vector< std::array<double,8> >& posTargets, const Cone<Vector, Scalar>& cone, const Container& cones);

	// Fill particles
	template < class Parameter, class Cone, typename Type1 > static void fill_particles_vectors(const Parameter& parameters, const Cone& cone, const std::vector< std::string >& shellList, std::vector<Type1>& pos_part, std::vector<Type1>& force_part, std::vector<Type1>& potential_part, std::vector<Type1>& a_part, const double thetay, const double thetaz);

	// Visualise octree
	template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > static void VizualizeOctree(const Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Type radius);

	// Read Angular position from previously computed catalog
	template <typename Integer, class Parameters  > static void ReadFromCat(const Integer icone, const Parameters& parameters, std::vector < std::array< double, 18 > >& catalogue);

};


/// \brief          Get list of files.
/// \details        Get list of files and directories.
/// \param[in]      String Directory name.
/// \param[in,out]  vector<string> List of files and directories.
void Miscellaneous::getFilesinDir(const std::string dirName, std::vector<std::string>& fileNames){

    DIR *pdir;
    struct dirent *pent;

    pdir=opendir(dirName.c_str()); //"." refers to the current dir
    if (!pdir){
	std::cout<<"opendir() failure; terminating"<<std::endl;
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<" "<<dirName<<std::endl;
	std::terminate();
    }
    errno=0;
    while ((pent=readdir(pdir))){
	fileNames.push_back(pent->d_name);
    }
    if (errno){
	std::cout<<"readdir() failure; terminating"<<std::endl;
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<" "<<dirName<<std::endl;
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
void Miscellaneous::Tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters){
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
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
void Miscellaneous::clear_shrink(Vectored& vector){

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
void Miscellaneous::fullclear_vector(std::vector<Type>& vector){
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
void Miscellaneous::TicketizeFunction(const Integer rank, const Integer ntasks, Function&& function){

    Integer ioticket = 1;
    // If we use the ticket system
    if (IOGROUPSIZE > 1){

        if (rank % IOGROUPSIZE > 0)
            MPI_Recv((void *) &ioticket, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	function();               

        if ((rank+1) % IOGROUPSIZE > 0 && rank+1 < ntasks)
                MPI_Send((void *) &ioticket, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    // No ticket system
    } else{
	function();
    }

}


// Fullsky cones generation
/// \brief          Generate fullsky cones.
/// \details        Generate fullsky cones.
/// \tparam         Cone cone type
/// \tparam         Sphere sphere type
/// \tparam         Dimension dimension
/// \tparam         Vector vector type
/// \tparam         Scalar scalar type
/// \tparam         Integer integer type
/// \param[in]      ncones Number of cones .
/// \param[in,out]  cone Cones generated
/// \param[in]      coneIfRot Cones generated
/// \param[in]      sphere Central sphere
/// \return         Generate the shape of fullsky cones.
template < class Cone, template <unsigned int, class, typename > class Sphere, unsigned int Dimension, class Vector, typename Scalar, typename Integer >
void Miscellaneous::GenerateFullskyCones(const Integer ncones, std::vector< Cone >& cone, std::vector< Cone >& coneIfRot, Sphere<Dimension, Vector, Scalar>& sphere){

    std::vector< std::array<double, 3> > tiling(ncones);
    // Cone angle from the maximum distance between points generated on the sphere. We multiply by an arbitrary factor which seems ideal to produce wide enough cones
    double alpha = 1.8*std::asin(sphere.template uniform<Dimension-1>(std::begin(tiling), std::end(tiling)).first/sphere.diameter());
    // Assign properties to each cone
    Utility::parallelize(ncones, [=, &tiling, &cone, &sphere, &alpha](const uint i){cone[i].assign(sphere.position(), tiling[i], alpha);});
    // No rotation for fullsky cones
    for(uint i = 0; i < ncones; i++){
	coneIfRot[i] = cone[i];
    }

}

// Narrow cones generation
/// \brief          Generate narrow cones.
/// \details        Generate narrow cones.
/// \tparam         Parameter Parameter type.
/// \tparam         Cone cone type
/// \tparam         Sphere sphere type
/// \tparam         Dimension dimension
/// \tparam         Vector vector type
/// \tparam         Scalar scalar type
/// \param[in]      parameters Parameter structure 
/// \param[in,out]  cone Cones generated
/// \param[in,out]      coneIfRot Cones generated
/// \param[in]      sphere Central sphere
/// \param[in,out]  rotm1 Rotation matrix for narrow cone cells
/// \param[in,out]  thetay Semi-angle for solid angle in direction y
/// \param[in,out]  thetaz Semi-angle for solid angle in direction z
template <  class Parameter, class Cone, template <unsigned int, class, typename > class Sphere, unsigned int Dimension, class Vector, typename Scalar >
void Miscellaneous::GenerateNarrowCones(const Parameter& parameters, std::vector< Cone >& cone, std::vector< Cone >& coneIfRot, Sphere<Dimension, Vector, Scalar>& sphere, std::array< std::array< double, 3 >, 3 >& rotm1, Scalar& thetay, Scalar& thetaz){

    static constexpr double pi = Constants<double>::pi();
    std::size_t found;
    std::vector< std::array<double, 3> > tiling(parameters.ncones), tilingbis(parameters.ncones);
    double theta_rot(0), phi_rot(0);
    std::array< std::array< double, 3 >, 3 > rotation = {{0}};
    std::vector<std::string> filelistingprior;
    std::string filelisting;
    const double eps = magrathea::Constants<double>::deg()*0.5; // 0.6 deg of buffer zone

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.celldir, filelistingprior);
    for(uint ifiling = 0; ifiling < filelistingprior.size(); ++ifiling){
        if(parameters.typefile == 1){
	    // Find any hdf5 file
            found = filelistingprior[ifiling].find(".h5");
            if(found!=std::string::npos){
	        filelisting = parameters.celldir+filelistingprior[ifiling];
	        break;
            }
        } else if(parameters.typefile == 2){
	    // For ASCII, search for the 'info_narrow_cone.txt' file
	    if(filelistingprior[ifiling] == "info_narrow_cone.txt"){
	        filelisting = parameters.celldir+filelistingprior[ifiling];
	        break;
	    }
        } else{
            std::cout<<"# WARNING : Narrow cones can only be computed using HDF5 or ASCII input files"<<std::endl;
	    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	    std::terminate();
        }
    }
    // Get informations from HDF5 files
    if(parameters.typefile == 1){
        TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "phi", phi_rot);
	TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "theta", theta_rot);
	TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "thetay", thetay);
	TReadHDF5::getAttribute(filelisting, "metadata/cone_info", "thetaz", thetaz);
    // Get informations from ASCII files
    } else if(parameters.typefile == 2){
	std::map<std::string, std::string> parameterASCII;
    	parameterASCII = Input::parse(filelisting);
    	phi_rot = std::stoul(parameterASCII["phi"]);
    	theta_rot = std::stod(parameterASCII["theta"]);
    	thetay = std::stod(parameterASCII["thetay"]);
    	thetaz = std::stod(parameterASCII["thetaz"]);
    } else{
	std::cout<<"# WARNING : Narrow cones can only be computed using HDF5 or ASCII input files"<<std::endl;
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	std::terminate();
    }
    thetaz *= magrathea::Constants<double>::deg();
    thetay *= magrathea::Constants<double>::deg();
    phi_rot *= magrathea::Constants<double>::deg();
    theta_rot *= magrathea::Constants<double>::deg();

    // Compute the rotation matrix to rotate the cone
    rotation[0][0] = std::cos(theta_rot)*std::cos(phi_rot);
    rotation[0][1] = std::cos(theta_rot)*std::sin(phi_rot);
    rotation[0][2] = - std::sin(theta_rot);
    rotation[1][0] = - std::sin(phi_rot);
    rotation[1][1] =  std::cos(phi_rot);
    rotation[1][2] = 0;
    rotation[2][0] = std::cos(phi_rot)*std::sin(theta_rot);
    rotation[2][1] = std::sin(phi_rot)*std::sin(theta_rot);
    rotation[2][2] =  std::cos(theta_rot);
    // Inverse matrix
    rotm1 = Utility::invMatrix3d(rotation);
    uint iloop(0);
    std::array<double, 3> ciblage;
    ciblage[0] = 1;
    ciblage[1] = 0;
    ciblage[2] = 0;
    double fullsky = 4*pi; 
    double portion = 2*thetay*(std::sin(thetaz) - std::cos(pi/2. + thetaz));
    // Inverse fraction of the sky
    uint fsp = fullsky/portion;
    tiling.resize(parameters.ncones*fsp);

    // Generate random points on the full sky. Then check is we have a number of points equal to 'ncones' in the area of interest
    sphere.template uniform<Dimension-1>(std::begin(tiling), std::end(tiling));
    for(uint i = 0; i < tiling.size(); ++i){
	double phi = std::atan2(tiling[i][1],tiling[i][0]);
	double theta = std::acos(tiling[i][2]/sphere.radius());
	// If point is inside the region of interest
	if(std::abs(phi) < thetay - eps && theta >  pi/2 - thetaz + eps && theta < pi/2 + thetaz - eps){
	    iloop++;
	}
    } 
    float irecuploop = 0;
    // Generally, we will not exactly have the good number of points, then need to iterate on the initial number of points on the full sky
    while(iloop != parameters.ncones){
        //std::cout<<irecuploop<<" "<<tiling.size()<<std::endl;
	if(irecuploop < 10) 
	    tiling.resize(tiling.size()+static_cast<int>((static_cast<float>(parameters.ncones)-static_cast<float>(iloop))*static_cast<float>(fsp)*(1/(1+irecuploop))));
	else 
	    tiling.resize(tiling.size()+parameters.ncones-iloop);
	sphere.template uniform<Dimension-1>(std::begin(tiling), std::end(tiling));
	iloop = 0;
	for(uint i = 0; i < tiling.size(); ++i){
	    double phi = std::atan2(tiling[i][1],tiling[i][0]);
	    double theta = std::acos(tiling[i][2]/sphere.radius());
	    if(std::abs(phi) < thetay - eps && theta >  pi/2 - thetaz + eps && theta < pi/2 + thetaz - eps){
	        iloop++;
	    }
	}
	irecuploop ++;
    }
    iloop = 0;
    // Put final result in tilingbis
    for(uint i = 0; i < tiling.size(); ++i){
	double phi = std::atan2(tiling[i][1],tiling[i][0]);
	double theta = std::acos(tiling[i][2]/sphere.radius());
	if(std::abs(phi) < thetay -eps && theta >  pi/2 - thetaz + eps && theta < pi/2 + thetaz - eps){
	    tilingbis[iloop][0] = tiling[i][0];
	    tilingbis[iloop][1] = tiling[i][1];
	    tilingbis[iloop][2] = tiling[i][2];
	    iloop++;
	}
    } 

    double resulting(0);
    resulting = sphere.diameter();
    // Check the maximum distance between two cones
    for (unsigned int i = 0; i < parameters.ncones; ++i) {
        for (unsigned int j = 0; j < parameters.ncones; ++j) {
            if (i != j) {
                double tmp = 0;
                for (unsigned int idim = 0; idim < 3; ++idim) {
                    tmp += std::pow(tilingbis[i][idim]-tilingbis[j][idim], 2);
                }
                tmp = std::sqrt(tmp);
                resulting = (resulting < tmp) ? (resulting) : (tmp);
            }
        }
    }

    // Assign an angle depending on the maximum distance between two cones. We multiply by an arbitrary factor which seems ideal to produce wide enough cones 
    double alpha = 1.8*std::asin(resulting/sphere.diameter());
    // Need a minimum angle to avoid thin cones. 0.1 = 6 degrees
    const double anglemin = 0.001;
    alpha = (anglemin > alpha) ? anglemin : alpha;
    std::cout<<"# Angle proposed : "<<1.8*std::asin(resulting/sphere.diameter())<<" angle chosen : "<<alpha<<std::endl;
    // Assign properties to cones
    Utility::parallelize(parameters.ncones, [=, &tilingbis, &cone](const uint i){cone[i].assign(sphere.position(), tilingbis[i], alpha);});
    // For narrow cones, need rotation
    for(uint i = 0; i < parameters.ncones; i++){
	coneIfRot[i] = cone[i];
	coneIfRot[i].base(0) = cone[i].base(0)*rotm1[0][0] + cone[i].base(1)*rotm1[0][1] + cone[i].base(2)*rotm1[0][2];
	coneIfRot[i].base(1) = cone[i].base(0)*rotm1[1][0] + cone[i].base(1)*rotm1[1][1] + cone[i].base(2)*rotm1[1][2];
	coneIfRot[i].base(2) = cone[i].base(0)*rotm1[2][0] + cone[i].base(1)*rotm1[2][1] + cone[i].base(2)*rotm1[2][2];
    }

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
template < class Octree, class Filelist,  typename Integer > 
void Miscellaneous::loadOctree(const Integer icone, Octree& octree, Filelist& conefile){
    
    octree.fullclear();
#ifdef VERBOSE
    std::cout<<"# Loading conefile["<<icone<<"] "<<conefile[icone]<<std::endl;
#endif
    // Load octree
    Input::load(octree, conefile[icone]);
#ifdef VERBOSE
    std::cout<<"# Conefile["<<icone<<"] loaded. Size : "<<octree.size()<<std::endl;
#endif
#ifdef VELOCITYFIELD
    // When putting an octree generated with gravity.h in an octree with gravity2.h, 
    // need to correct the position of data
    Utility::parallelize(octree.size(), [=, &octree](const uint i){
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
template < class Octree, class Cosmology, class Parameters, typename Real > 
void Miscellaneous::correctOctree(Octree& octree, const Cosmology& cosmology, Parameters& parameters, const Real h, const Real omegam, const Real lboxmpch, Real& amin){

    // Convert from Ramses Units to SI
    Input::sistemize(parameters, octree, h, omegam, lboxmpch);
    // Initialise dphida to zero
    Utility::parallelize(octree.size(), [=, &octree](const uint i){std::get<1>(octree[i]).dphidt() = 0;});
    // Apply correction to the octree,
    // also .update() is used twice (before and after corrections)
    Input::correct(parameters, octree, amin);
    // Convert from dphi/da to dphi/dt
    Utility::parallelize(octree.size(), [=, &octree](const uint i){std::get<1>(octree[i]).dphidt() *= Utility::interpolate2(std::get<1>(octree[i]).a(), std::get<1>(cosmology), std::get<2>(cosmology));});

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
template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > 
void Miscellaneous::VizualizeOctree(const Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Type radius){

    for(uint i = 0; i < octree.size(); i++){
        double x = std::get<0>(octree[i]).template center<Type,Position,Extent>(0);
	double y = std::get<0>(octree[i]).template center<Type,Position,Extent>(1);
	double z = std::get<0>(octree[i]).template center<Type,Position,Extent>(2);
        if(x*x+y*y+z*z < radius*radius){
	    // Output result on terminal
	    std::cout<<x<<" "<<y<<" "<<z<<" "<<std::get<1>(octree[i])<<" #vizualize"<<std::endl; 
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
template < class Vector, typename Scalar,  class Container>  
std::vector< std::array<double,8> > Miscellaneous::getTargets(const std::vector< std::array<double,8> >& posTargets, const Cone<Vector, Scalar>& cone, const Container& cones)
{   
    // Fill selection vector with -1
    std::vector<int> selection(posTargets.size(), -1);

    // Loop over all targets
    Utility::parallelize(posTargets.size(), [=, &posTargets, &cone, &cones, &selection](const unsigned int ivec){
        std::array<double,3> position;
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
                length += (cone.base(idim)-cone.vertex(idim))*(position[idim]-cone.vertex(idim));
            }
            length /= cone.template pow<2>(cone.length());
	    // Compute the ditance between cone base and target
            for (unsigned int idim = 0; idim < 3; ++idim) {
                reference += cone.template pow<2>(position[idim]-(cone.vertex(idim)+(cone.base(idim)-cone.vertex(idim))*length));
            }
	    // Loop over all the other cones	
            for (unsigned int icone = 0; icone < cones.size(); ++icone) {
		length = 0;
		distance = 0;
                if (cones[icone] != cone) {
		    // Compute scalar product of other cone base and target direction
                    for (unsigned int idim = 0; idim < 3; ++idim) {
                        length += (cones[icone].base(idim)-cone.vertex(idim))*(position[idim]-cones[icone].vertex(idim));
                    }
		    if(!(length < 0)){
                        length /= cones[icone].template pow<2>(cones[icone].length());
	    		// Compute the ditance between other cone base and target
                        for (unsigned int idim = 0; idim < 3; ++idim) {
                            distance += cones[icone].template pow<2>(position[idim]-(cones[icone].vertex(idim)+(cones[icone].base(idim)-cones[icone].vertex(idim))*length));
                        }
			// Check if other cone is closer to pixel
                        if (distance < reference) {
                            ok = false;
                            icone = cones.size();
                        }
		    } //  if length
                }  // if
            } // for icone
		if(ok){
		    selection[ivec] = ivec;
		}
        } //  if
    }); //  while

    // Erase targets which are not inside the cone
    selection.erase(std::remove(std::begin(selection), std::end(selection), -1), std::end(selection));
    std::vector< std::array<double,8> > pointsCible(selection.size());   
    // Put targets in vector
    Utility::parallelize(selection.size(), [=, &posTargets, &pointsCible, &selection](const unsigned int ivec){
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
template < class Parameter, class Cone, typename Type1>
void Miscellaneous::fill_particles_vectors(const Parameter& parameters, const Cone& cone, const std::vector< std::string >& shellList, std::vector<Type1>& pos_part, std::vector<Type1>& force_part, std::vector<Type1>& potential_part, std::vector<Type1>& a_part, const double thetay, const double thetaz){

    unsigned long long int marker1(0), marker2(0);

#ifdef VERBOSE
    std::cout<<"# Now loading the position, force and potential of particles"<<std::endl;
#endif
    Miscellaneous::fullclear_vector(pos_part);
    Miscellaneous::fullclear_vector(force_part);
    Miscellaneous::fullclear_vector(potential_part);
    // Loop over files
    for(uint ifiling = 0; ifiling < shellList.size(); ifiling++){
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
        const double factorpot = std::pow(unit_l*1e-2/unit_t, 2);
        const double factorforce = -aexp*unit_l*1e-2/(unit_t*unit_t);
        std::transform(std::begin(potential_part)+marker1, std::end(potential_part), std::begin(potential_part)+marker1, std::bind1st(std::multiplies<double>(), factorpot)); // SI units
        std::transform(std::begin(force_part)+marker2, std::end(force_part), std::begin(force_part)+marker2, std::bind1st(std::multiplies<double>(), factorforce)); // SI units
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
template <typename Integer, class Parameters >
void Miscellaneous::ReadFromCat(const Integer icone, const Parameters& parameters, std::vector < std::array< double, 18 > >& catalogue){
    // Create filename of catalogue
    std::string filename = parameters.outputdir + "/" + Output::name(parameters.base, "_", std::make_pair("%05d", icone), ".txt"); // Name of catalog, given icone, directory and base
#ifdef VERBOSE
    std::cout<<"# Cone "<<icone<<" Read angular position from "<<filename<<std::endl;
#endif
    // Open filename
    std::ifstream streaming(filename.c_str());
    streaming.unsetf(std::ios_base::skipws);
    uint size = std::count(std::istream_iterator<char>(streaming), std::istream_iterator<char>(), '\n');
    streaming.close();
    std::ifstream stream(filename.c_str());
    catalogue.resize(size);
    // Read halo catalogue or particle catalogue (for the latter there is no 'npart' column)

    for(unsigned int i = 0 ; i < size ; ++i){ 
	stream >> catalogue[i][0] >> catalogue[i][1] >> catalogue[i][2] >> catalogue[i][3] >> catalogue[i][4] >> catalogue[i][5] >> catalogue[i][6] >> catalogue[i][7] >> catalogue[i][8] >> catalogue[i][9] >> catalogue[i][10] >> catalogue[i][11] >> catalogue[i][12] >> catalogue[i][13] >> catalogue[i][14] >> catalogue[i][15] >> catalogue[i][16] >> catalogue[i][17];
    }   

}


#endif // MISCELLANEOUS_H_INCLUDED
