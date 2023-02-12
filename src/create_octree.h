/* ********************************** CREATE_OCTREE ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Create_octree
// DESCRIPTION :    Some create_octree function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           create_octree.h
/// \brief          Some create_octree function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef CREATE_OCTREE_H_INCLUDED
#define CREATE_OCTREE_H_INCLUDED


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


struct parameters_t {
    // Common
    real buffer;
    std::string celldir;
    std::string conedir;
    std::string conefmt;
    std::string evolfile;
    uint isfullsky;
    real mpc; 
    uint ncoarse;
    uint ncones;
    std::string paramfile;
    real rhoch2;
    uint seed;
    uint typefile;

    // Specific
    uint acorrection;
    uint allocation;
    std::string cellfmt;
    uint coarsecorrection;
    uint coarseonly;
    uint correction;
    std::string inputtype;
    uint microcoeff;
    std::string minicone;
    std::string partdir;
} parameters;


class Create_octree {

	//Methodes
	public:
	
	// Read parameter file
	template <class Parameters, class Map> static void ReadParamFile(Parameters& parameters, Map& parameter);

	// Preparation
	template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameters, class Cone, class FileList, class Sphere, typename Integer, typename Scalar, class Cosmology > static void PreparationHDF5_from_cells(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, const std::vector< Cone >& coneIfRot, const std::array< std::array< double, 3 >, 3 >& rotm1, const Scalar thetay, const Scalar thetaz, FileList& conefile, const Sphere& microsphere, const Scalar h,  const Scalar omegam, const Scalar lboxmpch, Scalar& amin, const Cosmology& cosmology);
	template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameters, class Cone, class FileList, class Sphere, typename Integer, typename Scalar > static void PreparationHDF5_from_particles(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, const std::vector< Cone >& coneIfRot, const Scalar thetay, const Scalar thetaz, FileList& conefile, const Sphere& microsphere);
	template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameters, class Cone, class Sphere, class FileList, class Index2, typename Integer, typename Scalar, class Cosmology > static void PreparationBinary(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, FileList& conefile, const Sphere& microsphere, Index2& filetree, const Scalar h,  const Scalar omegam, const Scalar lboxmpch, Scalar& amin, const Cosmology& cosmology);
	template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameters, class Cone, class FileList, class Sphere, typename Integer, typename Scalar, class Cosmology > static void PreparationASCII(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, const std::vector< Cone >& coneIfRot, const std::array< std::array< double, 3 >, 3 >& rotm1, FileList& conefile, const Sphere& microsphere, const Scalar h,  const Scalar omegam, const Scalar lboxmpch, Scalar& amin, const Cosmology& cosmology);


	// Create octree with particles
	template < class Parameters, typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > static void CreateOctreeWithCIC(bool& continue_refining, Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, std::vector<Type1>& pos_part, std::vector<Type1>& force_part, std::vector<Type1>& potential_part, std::vector<Type1>& a_part);
	template < class Parameters, typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > static void CreateOctreeWithTSC(bool& continue_refining, Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, std::vector<Type1>& pos_part, std::vector<Type1>& force_part, std::vector<Type1>& potential_part, std::vector<Type1>& a_part);

};



// Read parameter file
/// \brief          Read parameter file.
/// \details        Read and put in a structure the parameters.
/// \tparam         Parameters structure type
/// \tparam         Map map type
/// \param[in,out]  parameters Structure containing the parameters.
/// \param[in]      parameter Contains parameters to be rewritten
template <class Parameters, class Map> 
void Create_octree::ReadParamFile(Parameters& parameters, Map& parameter){

    parameters.paramfile = parameter["paramfile"];
    parameters.evolfile = parameter["evolfile"];
    parameters.celldir = parameter["celldir"];
    parameters.minicone = parameter["minicone"];
    parameters.conedir = parameter["conedir"];
    parameters.conefmt = parameter["conefmt"];
    parameters.cellfmt = parameter["cellfmt"];
    parameters.partdir = parameter["partdir"];
    parameters.inputtype = parameter["inputtype"];
    parameters.seed = std::stoul(parameter["seed"]);
    parameters.allocation = std::stoul(parameter["allocation"]);
    parameters.microcoeff = std::stoul(parameter["microcoeff"]);
    parameters.typefile = std::stoul(parameter["typefile"]);
    parameters.isfullsky = std::stoul(parameter["isfullsky"]);
    parameters.ncoarse = std::stoul(parameter["ncoarse"]);
    parameters.ncones = std::stoul(parameter["ncones"]);
    parameters.coarseonly = std::stoul(parameter["coarseonly"]);
    parameters.correction = std::stoul(parameter["correction"]);
    parameters.coarsecorrection = std::stoul(parameter["coarsecorrection"]);
    parameters.acorrection = std::stoul(parameter["acorrection"]);
    parameters.mpc = std::stod(parameter["mpc"]); 
    parameters.rhoch2 = std::stod(parameter["rhoch2"]);
    parameters.buffer = std::stod(parameter["buffer"]);

}



// Create Octree binary files
/// \brief          Create Octree binary files.
/// \details        Create Octree binary files from HDF5.
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \tparam         Parameters Parameters type
/// \tparam         Cone cone type
/// \tparam         Filelist filelist type
/// \tparam         Sphere sphere type
/// \tparam         Integer Integer type
/// \tparam         Scalar scalar type
/// \tparam         Cosmology cosmology type
/// \param[in,out]  octree Octree to be filled  
/// \param[in]      parameters Parameters structure
/// \param[in]      ntasks Number of tasks
/// \param[in]      rank Process rank
/// \param[in]      cone Cone properties
/// \param[in]      coneIfRot Cone with possibly rotation
/// \param[in]      rotm1 Rotation matrix for narrow cells
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
/// \param[in]      conefile Cone names in conedir
/// \param[in]      microsphere Central buffer zone for Octree
/// \param[in]      h Dimensionless Hubble parameter
/// \param[in]      omegam Matter density fraction
/// \param[in]      lboxmpch Size of simulation box
/// \param[in]      amin minimum value scale factor
/// \param[in]      cosmology Cosmological tables
template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameters, class Cone, class FileList, class Sphere, typename Integer, typename Scalar, class Cosmology > 
void Create_octree::PreparationHDF5_from_cells(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, const std::vector< Cone >& coneIfRot, const std::array< std::array< double, 3 >, 3 >& rotm1, const Scalar thetay, const Scalar thetaz, FileList& conefile, const Sphere& microsphere, const Scalar h,  const Scalar omegam, const Scalar lboxmpch, Scalar& amin, const Cosmology& cosmology){

    uint nfiles = zero;
    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());
    Cone conic, conicIfRot;
    std::size_t found;
    std::vector<std::string> filelistprior, filelistmini, filelist;

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.celldir, filelistprior);
    for(uint ifiling = 0; ifiling < filelistprior.size(); ++ifiling){
	// Check for all the .h5 files
	found = filelistprior[ifiling].find(".h5");  
	if(found!=std::string::npos){
	    double aexp(0);
	    TReadHDF5::getAttribute(parameters.celldir+filelistprior[ifiling], "metadata/cone_info","aexp",aexp);
	    // Check if the scale factor is normal. If it is then add the file to the list
	    if(std::isnormal(aexp)){
		filelist.push_back(parameters.celldir+filelistprior[ifiling]);
	    }
	}
    }
    Miscellaneous::fullclear_vector(filelistprior);
            
    // For narrow cones, also need to read the files from a 'minicone', i.e. a small spherical region around the observer (to allow for the first integrations steps)
    if(parameters.isfullsky == 0){ 
        std::vector<std::string> filelistpriormini;
	Miscellaneous::getFilesinDir(parameters.minicone, filelistpriormini);
	for(uint ifiling = 0; ifiling < filelistpriormini.size(); ++ifiling){
	    found = filelistpriormini[ifiling].find(".h5");  
	    if(found!=std::string::npos){
		double aexp(0);
		TReadHDF5::getAttribute(parameters.minicone+filelistpriormini[ifiling], "metadata/cone_info", "aexp", aexp);
		if(std::isnormal(aexp)){
		    filelistmini.push_back(parameters.minicone+filelistpriormini[ifiling]);
		}
	    }
	}	
        Miscellaneous::fullclear_vector(filelistpriormini);
    }

    // Pre-reserve some memory for the octree to prevent reallocation (slow)
    octree.reserve(parameters.allocation);
    for (uint icone = zero; icone < parameters.ncones; ++icone) {
        if (icone%static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
            octree.clear();
#ifdef VERBOSE
	    std::cout<<"# Creating cone "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
            conic = cone[icone];
	    conicIfRot = coneIfRot[icone]; 
	    // Shuffle the filenames so that different MPI tasks read different files
            std::shuffle(std::begin(filelist), std::end(filelist), engine1);
            nfiles = filelist.size();
	    // Read files and put cell information in octree
            for (uint ifile = zero; ifile < nfiles; ++ifile) {
	        Input::importhdf5<Position,Extent>(parameters, icone, rotm1, thetay, thetaz, conic, octree, filelist[ifile], [=, &octree](const Element& e){return Input::collide(octree, std::get<0>(e), microsphere, conicIfRot);});
#ifdef VERBOSE
		std::cout<<"# Cone : "<<icone<<" File "<<ifile+1<<"/"<<nfiles<<" Octree size : "<<octree.size()<<" Octree capacity : "<<octree.capacity()<<std::endl;
#endif
            }
	    // For narrow cones, read minicone files and put cell information in octree
	    if(parameters.isfullsky == 0){ 
                std::shuffle(std::begin(filelistmini), std::end(filelistmini), engine1);
                for (uint ifile = zero, nfiles = filelistmini.size(); ifile < nfiles; ++ifile) {
		    std::cout<<"# "<<filelistmini[ifile]<<std::endl;
		    Input::importfullhdf5<Position,Extent>(parameters, octree, filelistmini[ifile], [=, &octree](const Element& e){return Input::collide(octree, std::get<0>(e), microsphere, conicIfRot);});
#ifdef VERBOSE
		    std::cout<<"# MiniCone : "<<icone<<" File "<<ifile+1<<"/"<<nfiles<<" Octree size : "<<octree.size()<<std::endl;
#endif
                }
	    }
#ifdef VERBOSE
	    std::cout<<"# Before corrections, Cone : "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
	    // Apply correction on the octree
	    Miscellaneous::correctOctree(octree, cosmology, parameters, h,  omegam, lboxmpch, amin); // comment this line if you want to keep all cells in cones
#ifdef VERBOSE
	    std::cout<<"# Final Cone : "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
            Input::save(octree, conefile[icone]); 
        }
    }
    octree.fullclear();

}

// Create Octree binary files
/// \brief          Create Octree binary files.
/// \details        Create Octree binary files from HDF5.
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \tparam         Parameters Parameters type
/// \tparam         Cone cone type
/// \tparam         Filelist filelist type
/// \tparam         Sphere sphere type
/// \tparam         Integer Integer type
/// \tparam         Scalar scalar type
/// \tparam         Cosmology cosmology type
/// \param[in,out]  octree Octree to be filled  
/// \param[in]      parameters Parameters structure
/// \param[in]      ntasks Number of tasks
/// \param[in]      rank Process rank
/// \param[in]      cone Cone properties
/// \param[in]      coneIfRot Cone with possibly rotation
/// \param[in]      conefile Cone names in conedir
/// \param[in]      microsphere Central buffer zone for Octree
/// \param[in]      h Dimensionless Hubble parameter
/// \param[in]      omegam Matter density fraction
/// \param[in]      lboxmpch Size of simulation box
/// \param[in]      amin minimum value scale factor
/// \param[in]      cosmology Cosmological tables
template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameter, class Cone, class FileList, class Sphere, typename Integer, typename Scalar > 
void Create_octree::PreparationHDF5_from_particles(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameter& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, const std::vector< Cone >& coneIfRot, const Scalar thetay, const Scalar thetaz, FileList& conefile, const Sphere& microsphere){

    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());
    Cone conic, conicIfRot;
    std::size_t found;
    std::vector<std::string> filelistprior, filelistmini, filelist;
    // Get filenames in directory
    Miscellaneous::getFilesinDir(parameters.partdir, filelistprior);
    for(uint ifiling = 0; ifiling < filelistprior.size(); ++ifiling){
	// Check all files with suffix .h5
	found = filelistprior[ifiling].find(".h5");  
	if(found!=std::string::npos){
	    double aexp(0);
	    TReadHDF5::getAttribute(parameters.partdir+filelistprior[ifiling], "metadata/cone_info","aexp",aexp);
	    // Check if scale factor is normal
	    if(std::isnormal(aexp)){
		filelist.push_back(parameters.partdir+filelistprior[ifiling]);
	    }
	}
    }
    Miscellaneous::fullclear_vector(filelistprior);

    // Pre-reserve memory for the octree, to avoid reallocation (slow)
    octree.reserve(parameters.allocation);
    for (uint icone = zero; icone < parameters.ncones; ++icone) {
        if (icone%static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
	    bool continue_refining = true;
            octree.clear();
	    // Create a homogeneous octree at coarse level
	    Input::homogenize(octree.assign(parameters.ncoarse + std::log2(EXTENT), zero));
            conic = cone[icone];
	    conicIfRot = coneIfRot[icone]; 
	    std::cout<<"#size initial octree : "<<octree.size()<<std::endl;
	    // Only keep cells within the cone
    	    octree.resize(std::distance(std::begin(octree), std::remove_if(std::begin(octree), std::end(octree), [=, &octree](const Element& elem){return !(Input::collide(octree, std::get<0>(elem), microsphere, conicIfRot));})));
	    std::cout<<"#size cone octree : "<<octree.size()<<std::endl;
	    std::vector<float> pos_part, force_part, potential_part, a_part;
#ifdef VERBOSE
	    std::cout<<"# Creating cone "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
            std::shuffle(std::begin(filelist), std::end(filelist), engine1);
	    // Read all partiles within the cone
	    Miscellaneous::fill_particles_vectors(parameters, conicIfRot, filelist, pos_part, force_part, potential_part, a_part, thetay, thetaz);
	    // As long as we did not exceed the maximum refinement level, continue creating octree cells are finer levels
	    while(continue_refining){
		if(ORDER == 1){
		    Create_octree::CreateOctreeWithCIC(continue_refining, octree, parameters, pos_part, force_part, potential_part, a_part);
		} else if (ORDER == 2){
		    Create_octree::CreateOctreeWithTSC(continue_refining, octree, parameters, pos_part, force_part, potential_part, a_part);
		} else{
		    std::cout<<"# Please select -DORDER = 1 or 2 to compute the octree from particles"<<std::endl;
		    std::cout<<"# Error at file "<<__FILE__<<" line : "<<__LINE__<<std::endl;
		    std::terminate();
		}
	    }
#ifdef VERBOSE
	    std::cout<<"# Before Update , Cone : "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
	    octree.update();

	    // Cleaning the Octree at coarse level

	    // Get the coarse level from octree (should match the input coarse level + log2(EXTENT))
            const unsigned int ncoarse = (std::get<0>(*std::min_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level());
	    // Get the maximum level (At most, should be equal to ncoarse + max number of refinement, see inner functions) 
            const unsigned int nmax = (std::get<0>(*std::max_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level()); 
	    uint size = octree.size();
	    std::vector<uint> count(size), index;
	    uint counter(0);
	    // Count if for a given coarse cell the density is zero or NaN (in this case, also means that all the information is set to zero)
            Utility::parallelize(size, [=, &ncoarse, &count, &octree](const unsigned int i){count[i] = ((std::get<0>(octree[i]).level() == ncoarse) && (!std::isnormal(std::get<1>(octree[i]).rho())));});
	    // Total number of cells counted
            std::for_each(count.begin(), count.end(), [=, &counter](unsigned int& i){i = (i > zero) ? (++counter) : (zero);});
            index.resize(counter);
            Utility::parallelize(size, [=, &count, &index](const unsigned int i){if (count[i] > zero) {index[count[i]-one] = i;}});
	    // For these cells, average over all the closest cells with normal density
            Utility::parallelize(counter, [=, &ncoarse, &index, &octree](const unsigned int i){std::get<1>(octree[index[i]]) = Input::meanAll(octree, octree[index[i]], ncoarse);});
	    // Clear vector
            Utility::parallelize(count.begin(), count.end(), [](unsigned int& i){i = zero;});

	    // Cleaning the Octree at refined levels
            for (unsigned int n = ncoarse + 1; n <= nmax; ++n) {
	        // Count if for a given refined cell the density is zero or NaN (in this case, also means that all the information is set to zero)
		Utility::parallelize(size, [=, &count, &octree](const unsigned int i){count[i] = ((std::get<0>(octree[i]).level() == n) && (!std::isnormal(std::get<1>(octree[i]).rho())));});
		// If the density is not normal, then interpolate from coarser cell
                Utility::parallelize(size, [=, &count, &octree](const unsigned int i){if (count[i] > zero) {
		    if(ORDER == 1){
		        std::get<1>(octree[i]) = octree.cic(std::get<0>(octree[i]).template center<double, Position, Extent>(0), std::get<0>(octree[i]).template center<double, Position, Extent>(1), std::get<0>(octree[i]).template center<double, Position, Extent>(2));
		    } else if (ORDER == 2){
		        std::get<1>(octree[i]) = octree.tsc(std::get<0>(octree[i]).template center<double, Position, Extent>(0), std::get<0>(octree[i]).template center<double, Position, Extent>(1), std::get<0>(octree[i]).template center<double, Position, Extent>(2));
		    }		
		// If we could not interpolate, then give the value from parent cell
	        std::get<1>(octree[i]) = (std::get<1>(octree[i]) == Data()) ? std::get<1>(*octree.find(std::get<0>(octree[i]).parent())) : std::get<1>(octree[i]);}});
            }
#ifdef VERBOSE
	    std::cout<<"# Final Cone : "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
            Input::save(octree, conefile[icone]); 
        }
    }
    octree.fullclear();
}

// Create Octree binary files
/// \brief          Create Octree binary files.
/// \details        Create Octree binary files from Binary.
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \tparam         Parameters Parameters type
/// \tparam         Cone cone type
/// \tparam         Filelist filelist type
/// \tparam         Sphere sphere type
/// \tparam         Integer Integer type
/// \tparam         Scalar scalar type
/// \tparam         Cosmology cosmology type
/// \param[in,out]  octree Octree to be filled  
/// \param[in]      parameters Parameters structure
/// \param[in]      ntasks Number of tasks
/// \param[in]      rank Process rank
/// \param[in]      cone Cone properties
/// \param[in]      conefile Cone names in conedir
/// \param[in]      microsphere Central buffer zone for Octree
/// \param[in]      filetree Type of Octree containing strings
/// \param[in]      h dimensionless Hubble parameter
/// \param[in]      omegam Matter density fraction
/// \param[in]      lboxmpch Size of simulation box
/// \param[in]      amin minimum value scale factor
/// \param[in]      cosmology Cosmological tables
template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameter, class Cone, class Sphere, class FileList, class Index2, typename Integer, typename Scalar, class Cosmology > 
void Create_octree::PreparationBinary(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameter& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, FileList& conefile, const Sphere& microsphere, Index2& filetree, const Scalar h,  const Scalar omegam, const Scalar lboxmpch, Scalar& amin, const Cosmology& cosmology){

    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());
    std::vector<std::string> filelist;
    uint nfiles = zero;
    Cone conic;
    // Pre-reserve memory to avoid reallocation
    octree.reserve(parameters.allocation);
    // Get names from directory and put them in octree
    Input::filetree(filetree, parameters.celldir, parameters.cellfmt);
    for (uint icone = zero; icone < parameters.ncones; ++icone) {
        if (icone%static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
            octree.clear();
            filelist.clear();
            conic = cone[icone];
	    // Select filenames in cone
            Input::prepare(filelist, filetree, microsphere, conic);
	    // Shuffle so that different MPI tasks open different files
            std::shuffle(std::begin(filelist), std::end(filelist), engine1);
            nfiles = filelist.size();
	    // Fill octree with cells
            for (uint ifile = zero; ifile < nfiles; ++ifile) {
                Input::import(parameters, octree, filelist[ifile], [=, &octree](const Element& e){return Input::collide(octree, std::get<0>(e), microsphere, conic);});
            }
	    // Apply correction to octree
	    Miscellaneous::correctOctree(octree, cosmology, parameters, h,  omegam, lboxmpch, amin); // comment this line if you want to keep all cells in cones
            Input::save(octree, conefile[icone]);
        }
    }

}


// Create Octree binary files
/// \brief          Create Octree binary files.
/// \details        Create Octree binary files from ASCII iles.
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \tparam         Parameters Parameters type
/// \tparam         Cone cone type
/// \tparam         Filelist filelist type
/// \tparam         Sphere sphere type
/// \tparam         Integer Integer type
/// \tparam         Scalar scalar type
/// \tparam         Cosmology cosmology type
/// \param[in,out]  octree Octree to be filled  
/// \param[in]      parameters Parameters structure
/// \param[in]      ntasks Number of tasks
/// \param[in]      rank Process rank
/// \param[in]      cone Cone properties
/// \param[in]      coneIfRot Cone with possibly rotation
/// \param[in]      rotm1 Rotation matrix for narrow cells
/// \param[in]      conefile Cone names in conedir
/// \param[in]      microsphere Central buffer zone for Octree
/// \param[in]      h dimensionless Hubble parameter
/// \param[in]      omegam Matter density fraction
/// \param[in]      lboxmpch Size of simulation box
/// \param[in]      amin minimum value scale factor
/// \param[in]      cosmology Cosmological tables
template < template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container, class Parameter, class Cone, class FileList, class Sphere, typename Integer, typename Scalar, class Cosmology > 
void Create_octree::PreparationASCII(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameter& parameters, const Integer ntasks, const Integer rank, const std::vector< Cone >& cone, const std::vector< Cone >& coneIfRot, const std::array< std::array< double, 3 >, 3 >& rotm1, FileList& conefile, const Sphere& microsphere, const Scalar h,  const Scalar omegam, const Scalar lboxmpch, Scalar& amin, const Cosmology& cosmology){

    uint nfiles = zero;
    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());
    std::vector<std::string> filelist;
    Cone conic, conicIfRot;
    std::size_t found;
    std::vector<std::string> filelistprior;
    std::vector<std::string> filelistmini;

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.celldir, filelistprior);
    for(uint ifiling = 0; ifiling < filelistprior.size(); ++ifiling){
	// Select files with suffix .dat
	found = filelistprior[ifiling].find(".dat");  
	if(found!=std::string::npos){
	    filelist.push_back(parameters.celldir+filelistprior[ifiling]);
	}
    }
    Miscellaneous::fullclear_vector(filelistprior);
        
    // For narrow cones, need minicone directory
    if(parameters.isfullsky == 0){ 
        std::vector<std::string> filelistpriormini;
	Miscellaneous::getFilesinDir(parameters.minicone, filelistpriormini);
	for(uint ifiling = 0; ifiling < filelistpriormini.size(); ++ifiling){
	    found = filelistpriormini[ifiling].find(".dat");  
	    if(found!=std::string::npos){
		filelistmini.push_back(parameters.minicone+filelistpriormini[ifiling]);
	    }
	}	
        Miscellaneous::fullclear_vector(filelistpriormini);
    }
  
    octree.reserve(parameters.allocation);
    for (uint icone = zero; icone < parameters.ncones; ++icone) {
        if (icone%static_cast<uint>(ntasks) == static_cast<uint>(rank)) {
            octree.clear();
#ifdef VERBOSE
	    std::cout<<"# Creating cone "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
            conic = cone[icone];
	    conicIfRot = coneIfRot[icone]; 
	    // Shuffle names so that different MPI tasks do not read the same file
            std::shuffle(std::begin(filelist), std::end(filelist), engine1);
            nfiles = filelist.size();
	    // Fill octree with cell data
            for (uint ifile = zero; ifile < nfiles; ++ifile) {
	        Input::importascii<Position,Extent>(parameters, rotm1, octree, filelist[ifile], [=, &octree](const Element& e){return Input::collide(octree, std::get<0>(e), microsphere, conicIfRot);});
#ifdef VERBOSE
		std::cout<<"# Cone : "<<icone<<" File "<<ifile<<"/"<<nfiles<<" Octree size : "<<octree.size()<<std::endl;
#endif
            }
	    // For narrow cones, fill octree with minicone cell data
	    if(parameters.isfullsky == 0){ 
                std::shuffle(std::begin(filelistmini), std::end(filelistmini), engine1);
                for (uint ifile = zero, nfiles = filelistmini.size(); ifile < nfiles; ++ifile) {
		    Input::importascii<Position,Extent>(parameters, rotm1, octree, filelistmini[ifile], [=, &octree](const Element& e){return Input::collide(octree, std::get<0>(e), microsphere, conicIfRot);});
#ifdef VERBOSE
		    std::cout<<"# MiniCone : "<<icone<<" File "<<ifile<<"/"<<nfiles<<" Octree size : "<<octree.size()<<std::endl;
#endif
                }
	    }
#ifdef VERBOSE
	    std::cout<<"# Before corrections, Cone : "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
	    // Apply correction to octree
	    Miscellaneous::correctOctree(octree, cosmology, parameters, h,  omegam, lboxmpch, amin); // comment this line if you want to keep all cells in cones
#ifdef VERBOSE
	    std::cout<<"# Final Cone : "<<icone<<" Octree size : "<<octree.size()<<std::endl;
#endif
            Input::save(octree, conefile[icone]); 
        }
    }
    octree.fullclear();

}


// Compute CIC octree
/// \brief          Create Octree using CIC
/// \details        Create Octree using CIC from particles
/// \tparam         Parameters Parameters type
/// \tparam         Type1 Type1 type
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in,out]  continue_refining Bolean which controls if we continue refinement
/// \param[in,out]  octree Octree be to filled using particles 
/// \param[in]      parameters Parameters structure
/// \param[in]      pos_part Position of particles
/// \param[in]      force_part Force of particles
/// \param[in]      potential_part Potential of particles
/// \param[in]      a_part Scale factor of particles
template < class Parameters, typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > 
void Create_octree::CreateOctreeWithCIC(bool& continue_refining, Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, std::vector<Type1>& pos_part, std::vector<Type1>& force_part, std::vector<Type1>& potential_part, std::vector<Type1>& a_part){

    const unsigned int size = octree.size();
    std::vector<double> npart(size);
    const unsigned int n_mass_refine = 8; // Threshold of how many DM particles in a single cell to enable refinement
    const unsigned int n_leveldiff_refine = 8; // Threshold of how many refinement levels we can have at most

    const unsigned int lvlmax = (std::get<0>(*std::max_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level()); 
    const unsigned int leveldiff = lvlmax - (parameters.ncoarse + std::log2(EXTENT));

    // As long as the number of refinement levels is not superior or equal to 8, continue refining if possible 
    continue_refining = (leveldiff < n_leveldiff_refine);
    const float invextension = 1/(EXTENT/pow(2,lvlmax));  
    const float half = 0.5*EXTENT/pow(2,lvlmax); 
#ifdef VERBOSE
    std::cout<<"# Compute octree at level "<<lvlmax<<std::endl;
#endif
    // Loop over all the particles
    Utility::parallelize(pos_part.size()/3, [=, &octree, &npart, &pos_part, &force_part, &potential_part, &a_part, &half, &invextension, &lvlmax](const uint i){
        Index idxvertex;
        Data data; 
        unsigned long long int marker(0);
        double vratio(0);
	// Loop over the 8 neighboring cells of a given particle
        for(int iz = -1; iz <= 1; iz += 2){
	    for(int iy = -1; iy <= 1; iy += 2){
	        for(int ix = -1; ix <= 1; ix += 2){
		    // Create an index at the level of interest for the neighboring cell
		    idxvertex = idxvertex.template compute<Type,Position,Extent>(lvlmax, pos_part[3*i] + half*ix, pos_part[3*i+1] + half*iy, pos_part[3*i+2] + half*iz);
		    // Given an index in the octree that is consistent with the created index
       		    marker = std::distance(std::begin(octree), std::upper_bound(std::begin(octree), std::end(octree), Element(idxvertex, data), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
		    // If the index exists in the octree, compute CIC
		    if(std::get<0>(*(std::begin(octree)+marker-(marker > 0))) == idxvertex){
		        vratio = (1-std::abs(pos_part[3*i]-idxvertex.template center<Type,Position,Extent>(0))*invextension)
				*(1-std::abs(pos_part[3*i+1]-idxvertex.template center<Type,Position,Extent>(1))*invextension)
				*(1-std::abs(pos_part[3*i+2]-idxvertex.template center<Type,Position,Extent>(2))*invextension);
		        npart[marker - (marker > 0)] += vratio; 
		        std::get<1>(octree[marker - (marker > 0)]).rho() += vratio;
		        std::get<1>(octree[marker - (marker > 0)]).dphidx() += force_part[3*i]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).dphidy() += force_part[3*i+1]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).dphidz() += force_part[3*i+2]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).phi() += potential_part[i]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).a() += a_part[i]*vratio;
		    }		
	        }
	    }
	}
    });

    // Normalise the values by the mass
    Utility::parallelize(size, [=, &octree, &lvlmax](const uint i){
	if(std::get<0>(octree[i]).level() == lvlmax){
            double mass = std::get<1>(octree[i]).rho() + (std::get<1>(octree[i]).rho() == 0);
            std::get<1>(octree[i]).dphidx() /=  mass;   
            std::get<1>(octree[i]).dphidy() /=  mass;   
            std::get<1>(octree[i]).dphidz() /=  mass; 
            std::get<1>(octree[i]).phi() /=  mass; 
            std::get<1>(octree[i]).a() = (std::get<1>(octree[i]).a() - 1)/ mass; 
	}
    });

    // Count number of cells which contain a mass superior or equal to 8 DM particles
    const unsigned int num = std::count_if(npart.begin(), npart.end(), [=, &n_mass_refine](double& i){return i >= n_mass_refine;});
#ifdef VERBOSE
    std::cout<<"# Number of cells to refine : "<<num<<std::endl;
#endif
    // If no cell contains more than the threshold value, then stop the refinement
    if(num == 0){
#ifdef VERBOSE
	std::cout<<"# Cannot refine anymore "<<std::endl;
#endif
	continue_refining = false;
    }
    // If it is still possible to refine (density high enough and maximum refinement level not reached yet)
    if(continue_refining){
	for(unsigned int i = 0; i < size; i++){
	    // Refine cell (add 8 son cells in the octree) if the mass in one cell exceed the threshold
	    if(npart[i] >= n_mass_refine){
		octree.refine(octree.begin()+i);
	    }
	}
    }
    octree.update();
}

// Compute TSC octree
/// \brief          Create Octree using TSC
/// \details        Create Octree using TSC
/// \tparam         Parameters Parameters type
/// \tparam         Type1 Type1 type
/// \tparam         Octree octree type
/// \tparam         Type type type
/// \tparam         Index index type
/// \tparam         Data data type
/// \tparam         Dimension Number of dimensions
/// \tparam         Position position type
/// \tparam         Extent extent type
/// \tparam         Element element type
/// \tparam         Container container type
/// \param[in,out]  continue_refining Bolean which controls if we continue refinement
/// \param[in,out]  octree Octree be to filled using particles
/// \param[in]      parameters Parameters structure
/// \param[in]      pos_part Position of particles
/// \param[in]      force_part Force of particles
/// \param[in]      potential_part Potential of particles
/// \param[in]      a_part Scale factor of particles
template < class Parameters, typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > 
void Create_octree::CreateOctreeWithTSC(bool& continue_refining, Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, const Parameters& parameters, std::vector<Type1>& pos_part, std::vector<Type1>& force_part, std::vector<Type1>& potential_part, std::vector<Type1>& a_part){

    const unsigned int size = octree.size();
    std::vector<unsigned int> npart(size);
    const unsigned int n_mass_refine = 8; // Threshold of how many DM particles in a single cell to enable refinement
    const unsigned int n_leveldiff_refine = 8; // Threshold of how many refinement levels we can have at most

    const unsigned int lvlmax = (std::get<0>(*std::max_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level()); 
    const unsigned int leveldiff = lvlmax - (parameters.ncoarse + std::log2(EXTENT));
    // As long as the number of refinement levels is not superior or equal to 8, continue refining if possible 
    continue_refining = (leveldiff < n_leveldiff_refine);
    const double half = 0.5*EXTENT/pow(2,lvlmax); 
    const float twohalves = 2*half; 
#ifdef VERBOSE
    std::cout<<"# Compute octree at level "<<lvlmax<<std::endl;
#endif
    // Loop over all the particles
    Utility::parallelize(pos_part.size()/3, [=, &octree, &npart, &pos_part, &force_part, &potential_part, &a_part, &lvlmax, &twohalves](const uint i){
        Index idxvertex;
        Data data; 
	std::array<double, Index::dimension()> dist;
        unsigned long long int marker(0);
        double vratio(0), weightx(0), weighty(0), weightz(0);
	// Loop over the 27 neighboring cells of a given particle
        for(int iz = -1; iz <= 1; iz += 1){
	    const int aiz = abs(iz);
	    for(int iy = -1; iy <= 1; iy += 1){
		const int aiy = abs(iy);
	        for(int ix = -1; ix <= 1; ix += 1){
		    const int aix = abs(ix);
		    // Create an index at the level of interest for the neighboring cell
		    idxvertex = idxvertex.template compute<Type,Position,Extent>(lvlmax, pos_part[3*i] + twohalves*ix, pos_part[3*i+1] + twohalves*iy, pos_part[3*i+2] + twohalves*iz);
		    // Given an index in the octree that is consistent with the created index
       		    marker = std::distance(std::begin(octree), std::upper_bound(std::begin(octree), std::end(octree), Element(idxvertex, data), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
		    // If the index exists in the octree, compute TSC
		    if(std::get<0>(*(std::begin(octree)+marker-(marker > 0))) == idxvertex){
			for (unsigned int idim = 0; idim < Index::dimension(); ++idim) {
                    	    dist[idim] = std::abs((std::get<0>(*(std::begin(octree)+marker-(marker > 0))).template center<Type, Position, Extent>(idim)-pos_part[3*i+idim])/(twohalves));
			}
			weightx = aix*0.5*(1.5-dist[0])*(1.5-dist[0]) + (1-aix)*(0.75-dist[0]*dist[0]);
			weighty = aiy*0.5*(1.5-dist[1])*(1.5-dist[1]) + (1-aiy)*(0.75-dist[1]*dist[1]);
			weightz = aiz*0.5*(1.5-dist[2])*(1.5-dist[2]) + (1-aiz)*(0.75-dist[2]*dist[2]);
		        vratio = weightx*weighty*weightz;
			if(aiz+aiy+aix == 0){
		            npart[marker - (marker > 0)] += 1;
			}
		        std::get<1>(octree[marker - (marker > 0)]).rho() += vratio;
		        std::get<1>(octree[marker - (marker > 0)]).dphidx() += force_part[3*i]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).dphidy() += force_part[3*i+1]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).dphidz() += force_part[3*i+2]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).phi() += potential_part[i]*vratio;
		        std::get<1>(octree[marker - (marker > 0)]).a() += a_part[i]*vratio;
		    }		
	        }
	    }  
        }
    });

    // Normalise the values by the mass
    Utility::parallelize(size, [=, &octree, &lvlmax](const uint i){
	if(std::get<0>(octree[i]).level() == lvlmax){
            double mass = std::get<1>(octree[i]).rho() + (std::get<1>(octree[i]).rho() == 0);
            std::get<1>(octree[i]).dphidx() /=  mass;   
            std::get<1>(octree[i]).dphidy() /=  mass;   
            std::get<1>(octree[i]).dphidz() /=  mass; 
            std::get<1>(octree[i]).phi() /=  mass; 
            std::get<1>(octree[i]).a() = (std::get<1>(octree[i]).a() - 1)/ mass;
	}
    });

    // Count number of cells which contain a mass superior or equal to 8 DM particles
    unsigned int num = std::count_if(npart.begin(), npart.end(), [](unsigned int& i){return i >= n_mass_refine ;});
#ifdef VERBOSE
    std::cout<<"# Number of cells to refine : "<<num<<std::endl;
#endif
    if(num == 0){
#ifdef VERBOSE
	std::cout<<"# Cannot refine anymore "<<std::endl;
#endif
	continue_refining = false;
    }
    // If it is still possible to refine (density high enough and maximum refinement level not reached yet)
    if(continue_refining){ 
	for(unsigned int i = 0; i < size; i++){
	    // Refine cell (add 8 son cells in the octree) if the mass in one cell exceed the threshold
	    if(npart[i] >= n_mass_refine){
		octree.refine(octree.begin()+i);
	    }
	}
    }
    octree.update();
}

#endif // CREATE_OCTREE_H_INCLUDED
