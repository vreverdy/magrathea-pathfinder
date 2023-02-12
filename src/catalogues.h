/* ********************************** CATALOGUES ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Catalogues
// DESCRIPTION :    Some catalogues function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           catalogues.h
/// \brief          Some catalogues function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef CATALOGUES_H_INCLUDED
#define CATALOGUES_H_INCLUDED


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
#include "lensing.h"

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
    std::string base;
    real cat_accuracy;
    uint firstcone;
    uint halos;
    std::string beam;
    uint lastcone;
    uint npart;
    uint nsteps;
    real openingmin;
    std::string outputdir;
    std::string outputprefix;
    std::string plane;
    std::string sourcedir;
    std::string stop_bundle;
    uint use_previous_catalogues;
    real v0x;
    real v0y;
    real v0z;
    real zmin;
    real zmax;
} parameters;


class Catalogues {

	//Methodes
	public:

	template <class Parameters, class Map> static void ReadParamFile(Parameters& parameters, Map& parameter);	
	// Read particles
	template <typename Integer, class Parameter> static void ReadParticlesHDF5(const Integer rank, const Parameter& parameters, std::vector< std::array< double, 8 > >& caractVect_source);
	template <typename Integer, class Parameter> static void ReadParticlesASCII(const Integer rank, const Parameter& parameters, std::vector< std::array< double, 8 > >& caractVect_source);
	
	// Catalogues
	template < int Order = ORDER, bool RK4 = true, bool Verbose = false, class Point, class Cosmology, class Octree, class Type, class Parameter > static std::array< std::array<double, 2>, 2> newtonMethod2d(const Point& vobs, const Point& observer, const Point& trueTarget, const std::array<double, 2>& target, const Point& velocity, std::array< std::array<double, 2>, 2 >& jacobian, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h, std::vector<double>& redshifts, double& interpRef, const unsigned int iteration);
	template < int Order = ORDER, bool RK4 = true, bool Verbose = false, class Point, class Cosmology, class Octree, class Type, class Parameter > static std::array< std::array<double,2>, 2 > iterateNewtonMethod(const Point& vobs, const Point& observer, const Type phi, const Type theta, const Point& target, const Point& velocity, std::array< std::array<double, 2>, 2 >& jacobian, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h, std::vector<double>& redshifts, double& interpRef);
	template < class Point, class Cosmology, class Octree, class Type, class Parameter >  static void relCat(const Point& vobs, const std::array< std::array< double, 3 >, 3 >& rotm1, std::string& nomOutput, const Point& observer, const std::vector < std::array<double, 8> >& targets_position, const std::vector < std::array< double, 18 > >& previous_catalogue, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h); 
	template < class Point, class Cosmology, class Octree, class Type, class Parameter >  static void relCat_with_previous_cat(const Point& vobs, std::string& nomOutput, const Point& observer, std::vector < std::array< double, 18 > >& previous_catalogue, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h); 
	template < class Point, class Cosmology, class Octree, class Type, class Parameter >  static void relCat_with_previous_cat_forward(const Point& vobs, std::string& nomOutput, const Point& observer, std::vector < std::array< double, 18 > >& previous_catalogue, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h); 

};

// Read parameter file
/// \brief          Read parameter file.
/// \details        Read and put in a structure the parameters.
/// \tparam         Parameters structure type
/// \tparam         Map map type
/// \param[in,out]  parameters Structure containing the parameters.
/// \param[in]      parameter Contains parameters to be rewritten
/// \return         Filled parameters structure.
template <class Parameters, class Map> 
void Catalogues::ReadParamFile(Parameters& parameters, Map& parameter){

    parameters.npart = std::stoul(parameter["npart"]);
    parameters.zmin = std::stod(parameter["zmin"]);
    parameters.zmax = std::stod(parameter["zmax"]);
    parameters.firstcone = std::stoul(parameter["firstcone"]);
    parameters.lastcone = std::stoul(parameter["lastcone"]);
    parameters.cat_accuracy = std::stod(parameter["cat_accuracy"]);
    parameters.halos = std::stoul(parameter["halos"]);
    parameters.beam = parameter["beam"];
    parameters.conedir = parameter["conedir"];
    parameters.sourcedir = parameter["sourcedir"];
    parameters.outputdir = parameter["outputdir"];
    parameters.outputprefix = parameter["outputprefix"];
    parameters.paramfile = parameter["paramfile"];
    parameters.evolfile = parameter["evolfile"];
    parameters.mpc = std::stod(parameter["mpc"]); 
    parameters.rhoch2 = std::stod(parameter["rhoch2"]);
    parameters.typefile = std::stoul(parameter["typefile"]);
    parameters.isfullsky = std::stoul(parameter["isfullsky"]);
    parameters.v0x = std::stod(parameter["v0x"]);
    parameters.v0y = std::stod(parameter["v0y"]);
    parameters.v0z = std::stod(parameter["v0z"]);
    parameters.conefmt = parameter["conefmt"];
    parameters.stop_bundle = parameter["stop_bundle"];
    parameters.plane = parameter["plane"];
    parameters.base = parameter["base"];
    parameters.celldir = parameter["celldir"];
    parameters.ncoarse = std::stoul(parameter["ncoarse"]);
    parameters.ncones = std::stoul(parameter["ncones"]);
    parameters.use_previous_catalogues = std::stoul(parameter["use_previous_catalogues"]);
    parameters.seed = std::stoul(parameter["seed"]);
    parameters.openingmin = std::stod(parameter["openingmin"]);
    parameters.nsteps = std::stoul(parameter["nsteps"]);
    parameters.buffer = std::stod(parameter["buffer"]);
}


// Read Particles
/// \brief          Read Particles from HDF5 files
/// \details        Read particles from HDF5 files and put them in vectors
/// \tparam         Integer Integer type
/// \tparam         Parameter Parameter type
/// \param[in]      rank Rank
/// \param[in]      parameters Parameters structure
/// \param[in]      caractVect_source Source caracteristics
template <typename Integer, class Parameter> 
void Catalogues::ReadParticlesHDF5(const Integer rank, const Parameter& parameters, std::vector< std::array< double, 8 > >& caractVect_source){

    std::vector<std::string> filelistprior;
    std::size_t found;
    std::vector<std::string> partlist;
    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());

    double myamin = one/(one+parameters.zmin);
    double myamax = one/(one+parameters.zmax);

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.sourcedir, filelistprior);

    for(uint ifiling = 0; ifiling < filelistprior.size(); ++ifiling){
	if(parameters.halos){
	    // If sources are HDF5 haloes, we are looking for the file containing the keyword 'hfprop'
	    found = filelistprior[ifiling].find("hfprop"); 
	    if (found!=std::string::npos){
		partlist.push_back(parameters.sourcedir+filelistprior[ifiling]);
	    }
	} else{
	    // If sources are HDF5 DM particles, we are looking for all the files containing the keyword 'shell'
	    found = filelistprior[ifiling].find("shell");
	    if (found!=std::string::npos){
		if(!parameters.zmax)
		    partlist.push_back(parameters.sourcedir+filelistprior[ifiling]);
		else{
		    double amaxing(0), amining(0);
		    TReadHDF5::getAttribute(parameters.sourcedir+filelistprior[ifiling], "metadata/cone_info", "amax", amaxing); 
		    TReadHDF5::getAttribute(parameters.sourcedir+filelistprior[ifiling], "metadata/cone_info", "amin", amining); 
		    // Select shells within the redshifts (or scale factors) of interest
		    if(myamax <= amining &&  myamin >= amaxing)
			partlist.push_back(parameters.sourcedir+filelistprior[ifiling]);
		}
	    }
	}
    }
    if (rank == 0 && !(partlist.size())){
	std::cout<<"# WARNING:  No files detected: pleasure ensure that the halo file contains 'hfprop' or that particle files contain '.h5' in their name"<<std::endl;   
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	std::terminate();
    }
    // Random shuffling of the files (this ensures that different MPI tasks open different files)
    std::shuffle(std::begin(partlist), std::end(partlist), engine1);
    double fraction(0);
    double nparttot(0);
    if(!parameters.halos){
	// Get total number of particles in files within the redshift range
	for(uint ifiling = 0; ifiling < partlist.size(); ++ifiling){
	    std::vector<double> tmp;
	    TReadHDF5::fillVectors_part(partlist[ifiling], "metadata", "npart_file", tmp);
	    nparttot += tmp[0];
	}
	// Compute the fraction of particles that we want (given by npart) divided by the total number of particles within redshifts
	fraction = double(parameters.npart)/nparttot; 
    }

    std::array< double, 8 > caract_source; 
    if(parameters.halos){ 
	std::vector<double> pos_halos, vel_halos, id_halos; 
	std::vector<uint> npart_halos;
#ifdef VERBOSE
	if (rank == 0) 
	    std::cout<<"# Targets : Halos"<<std::endl;
#endif
	// Fill vectors from HDF5 files with halo properties
	TReadHDF5::fillVectors_part(partlist[0], "data", "position_halo", pos_halos, "velocity_halo", vel_halos, "identity_halo", id_halos, "npart_halo", npart_halos);
	for(uint i = 0; i < id_halos.size(); i++){
	    for(uint j = 0; j <  3; j++){
		caract_source[j] = pos_halos[3*i+j];
		caract_source[j+3] = vel_halos[3*i+j];
	    }
	    caract_source[6] = id_halos[i];
	    caract_source[7] = npart_halos[i];
	    caractVect_source.push_back(caract_source);
	}
    }
    else{ // if particles
#ifdef VERBOSE
	if (rank == 0){
	    std::cout<<"# Targets : Particles"<<std::endl;
	    std::cout<<"# Number of particles wanted : "<<parameters.npart<<" "<<" Total number of particles : "<<nparttot<<" Fraction "<<fraction <<" nfiles : "<<partlist.size()<<std::endl;
	}
#endif
	for(uint ifiling = 0; ifiling < partlist.size(); ++ifiling){
#ifdef VERBOSE
	    if(rank == 0) 
		std::cout<<"# File used : "<<partlist[ifiling]<<std::endl;
#endif
    	    std::vector<float> pos_part, vel_part; 
	    std::vector<double> id_part;
	    // Fill vectors from HDF5 files with DM particles properties (we take a random fraction given by 'fraction')
	    TReadHDF5::fillVectors_part(fraction, partlist[ifiling], "data", "position_part", pos_part, "velocity_part", vel_part, "identity_part_ramses", id_part);	
	    for(uint i = 0; i < id_part.size(); i++){
		for(uint j = 0; j <  3; j++){
		    caract_source[j] = pos_part[3*i+j];
		    caract_source[j+3] = vel_part[3*i+j];
		}
		caract_source[6] = id_part[i];
		caract_source[7] = 1;
		caractVect_source.push_back(caract_source);
	    }   
	}
    }
}

// Read Particles
/// \brief          Read Particles from ASCII files
/// \details        Read particles from ASCII files and put them in vectors
/// \tparam         Integer Integer type
/// \tparam         Parameter Parameter type
/// \param[in]      rank Rank
/// \param[in]      parameters Parameters structure
/// \param[in]      caractVect_source Source caracteristics
template <typename Integer, class Parameter> 
void Catalogues::ReadParticlesASCII(const Integer rank, const Parameter& parameters, std::vector< std::array< double, 8 > >& caractVect_source){
  
    std::vector<std::string> filelistprior;
    std::size_t found;
    std::vector<std::string> partlist;
    std::mt19937 engine1(parameters.seed > zero ? parameters.seed+rank : std::random_device()());

    // Get all filenames in directory
    Miscellaneous::getFilesinDir(parameters.sourcedir, filelistprior);


   for(uint ifiling = 0; ifiling < filelistprior.size(); ++ifiling){
	if(parameters.halos){
	    // If sources are ASCII haloes, we are looking for the file containing the keyword 'hfprop'
	    found = filelistprior[ifiling].find("hfprop"); 
	    if (found!=std::string::npos){
		partlist.push_back(parameters.sourcedir+filelistprior[ifiling]);
	    }
	} else{
	    // If sources are ASCII DM particles, we are looking for all the .dat files
	    found = filelistprior[ifiling].find(".dat");
	    if (found!=std::string::npos){
		partlist.push_back(parameters.sourcedir+filelistprior[ifiling]);
	    }
	}
    }
    if (rank == 0 && !(partlist.size())){
	std::cout<<"# WARNING:  No files detected: pleasure ensure that the halo file contains 'hfprop' or that particle files contain '.dat' in their name"<<std::endl; 
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	std::terminate();
    }  

    // Random shuffling of the files (this ensures that different MPI tasks open different files)
    std::shuffle(std::begin(partlist), std::end(partlist), engine1);

    double size(0);
    if(parameters.halos){ 
#ifdef VERBOSE
	if (rank == 0) 
	    std::cout<<"# Targets : Halos"<<std::endl;
#endif
	std::ifstream streaming(partlist[0].c_str());
	streaming.unsetf(std::ios_base::skipws);
	size = std::count(std::istream_iterator<char>(streaming), std::istream_iterator<char>(), '\n');
	streaming.close();
	std::ifstream stream(partlist[0].c_str());
	caractVect_source.resize(size);
	// Fill vectors from ASCII files with DM haloes properties (only 1 file)
	for(unsigned int i = 0 ; i < size ; ++i){ 
	    stream >> caractVect_source[i][0] >> caractVect_source[i][1] >> caractVect_source[i][2] >> caractVect_source[i][3] >> caractVect_source[i][4] >> caractVect_source[i][5] >> caractVect_source[i][6] >> caractVect_source[i][7];
	}   
    }
    else{ // if particles
	double caract_size(0);
	for(uint ifiling = 0; ifiling < partlist.size(); ++ifiling){
#ifdef VERBOSE
	    if(rank == 0) 
		std::cout<<"# File used : "<<partlist[ifiling]<<std::endl;
#endif
	    std::ifstream streaming(partlist[ifiling].c_str());
	    streaming.unsetf(std::ios_base::skipws);
	    size = std::count(std::istream_iterator<char>(streaming), std::istream_iterator<char>(), '\n');
	    streaming.close();
	    std::ifstream stream(partlist[ifiling].c_str());
	    caractVect_source.resize(caract_size + size);
	    // Fill vectors from ASCII files with DM particles properties
	    for(unsigned int i = 0 ; i < size ; ++i){ 
	        stream >> caractVect_source[caract_size+i][0] >> caractVect_source[caract_size+i][1] >> caractVect_source[caract_size+i][2] >> caractVect_source[caract_size+i][3] >> caractVect_source[caract_size+i][4] >> caractVect_source[caract_size+i][5] >> caractVect_source[caract_size+i][6] >> caractVect_source[caract_size+i][7];
	    }    
	    caract_size += size;
	}
    }
#ifdef VERBOSE
    if (rank == 0){
	std::cout<<"# Targets : Particles"<<std::endl;
	std::cout<<"# Number of particles wanted : "<<parameters.npart<<" "<<" Total number of particles : "<<caractVect_source.size()<<" nfiles : "<<partlist.size()<<std::endl;
    }
#endif

    // Only take a fraction of the DM particles if we want less that the total number within redshift range
    if(parameters.halos == 0 && parameters.npart < caractVect_source.size()){
	std::vector< std::array< double, 8 > > caractVect_source_tmp;
#ifdef GCCBELOW7
        std::experimental::sample(caractVect_source.begin(), caractVect_source.end(), std::back_inserter(caractVect_source_tmp), parameters.npart, std::mt19937{std::random_device{}()});
#else
        std::sample(caractVect_source.begin(), caractVect_source.end(), std::back_inserter(caractVect_source_tmp), parameters.npart, std::mt19937{std::random_device{}()});
#endif
	caractVect_source = caractVect_source_tmp;
    } 
}

// Methode de Newton 
/// \brief          Iterate once using Newton Method.
/// \details        Fit the better initial conditions for the photon 
///		    that intersects the source.
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2 for TSC
///                 or -1 for an homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam	    Point point type.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Parameter Parameter type
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      observer Observer cartesian position.
/// \param[in]      trueTarget True cartesian target position.
/// \param[in]      target Observed angular target position.
/// \param[in]      velocity Target velocity.
/// \param[in,out]  jacobian Jacobian matrix
/// \param[in]      parameters Parameters structure
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      h Dimensionless Hubble parameter
/// \param[in,out]  redshifts Vector containing the redshift decomposition for a given source.
/// \param[in,out]  interpRef value for interpolation of the bundle.
/// \param[in]      iteration Number of iterations for the root-finder
/// \return         2x2 array with NEW initial angles given by newton method
///		    and angle difference at the source between source and photon
template < int Order, bool RK4, bool Verbose, class Point, class Cosmology, class Octree, class Type, class Parameter > 
std::array< std::array<double, 2>, 2> Catalogues::newtonMethod2d(const Point& vobs, const Point& observer, const Point& trueTarget, const std::array<double, 2>& target, const Point& velocity, std::array< std::array<double, 2>, 2 >& jacobian, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h, std::vector<double>& redshifts, double& interpRef, const unsigned int iteration) 
{

    magrathea::Evolution<Photon<double, 3> > trajectory, trajectory_born;
    Photon<double, 3> photon;
    Point finalPoint;
    std::array<std::array<double, 2>, 2> result;
    double distTarget = std::sqrt(pow(trueTarget[0],2)+pow(trueTarget[1],2)+pow(trueTarget[2],2));
    std::array< std::array<double, 2>, 2 > jacobianinv;
    unsigned int firstid(0);
    static const double c = magrathea::Constants<double>::c();
    static const double c2 = magrathea::Constants<double>::c2();
    Point kiTarget;

    // Initialise photon
    const double phi = target[0];
    const double theta = target[1];
    photon = Integrator::launch(observer[0], observer[1], observer[2], phi, theta);	
    trajectory.append(photon);
    // Integrate on null geodesics until the photon reaches the source
    Integrator::integrate(trajectory, "radius", distTarget, cosmology, octree, vobs, length, parameters.nsteps);

    const unsigned int marked = trajectory.size() - 1;
    firstid = marked - (marked > 0);
    // Check if the photon reached the source
    if(marked == 0 || trajectory[marked].chi() < distTarget){
	result[0][0] = 42;
	result[0][1] = 42;
	result[1][0] = 0;	
	result[1][1] = firstid;
	return result;
    }
    // Prepare linar interpolation
    double previous = trajectory[firstid].chi();
    double next = trajectory[firstid+1].chi();
    double f = (next-distTarget)/(next-previous);

    // Estimate redshift, scale factor and position at the source
    // Redshift FLRW + Potential
    redshifts[1] = trajectory[firstid].redshift()*f + trajectory[firstid+1].redshift()*(1-f);
    double scale_factor =  trajectory[firstid].a()*f + trajectory[firstid+1].a()*(1-f);
    finalPoint[0] = trajectory[firstid].x()*f + trajectory[firstid+1].x()*(1-f);
    finalPoint[1] = trajectory[firstid].y()*f + trajectory[firstid+1].y()*(1-f);
    finalPoint[2] = trajectory[firstid].z()*f + trajectory[firstid+1].z()*(1-f);

    // If we compute the lensing matrix with a bundle
    if(parameters.beam == "bundle"){
	// Set plane
        if(parameters.plane == "sachs"){
	    kiTarget[0] = trajectory[firstid].dxdl()*f + trajectory[firstid+1].dxdl()*(1-f);
	    kiTarget[1] = trajectory[firstid].dydl()*f + trajectory[firstid+1].dydl()*(1-f);
	    kiTarget[2] = trajectory[firstid].dzdl()*f + trajectory[firstid+1].dzdl()*(1-f);	
        } else if (parameters.plane == "normal"){
	    kiTarget = finalPoint;
        } else if(parameters.plane == "exact"){
	    std::cout<<"# Jacobian 'exact' not yet implemented !"<<std::endl;
        } else{
	    std::cout<<"# WARNING : with beam = 'bundle', please choose plane = 'sachs', 'normal' or 'exact'"<<std::endl;
	    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	    std::terminate();
        }

        // Interpolation  for bundle
        if(parameters.stop_bundle == "redshift"){
	    interpRef = redshifts[1]; 
        } else if(parameters.stop_bundle == "a"){
	    interpRef = scale_factor; 	
        } else if((parameters.stop_bundle == "t") || (parameters.stop_bundle == "eta")){
	    interpRef = trajectory[firstid].t()*f + trajectory[firstid+1].t()*(1-f); 	
        } else if(parameters.stop_bundle == "lambda"){
	    interpRef = trajectory[firstid].lambda()*f + trajectory[firstid+1].lambda()*(1-f); 	
        } else if((parameters.stop_bundle == "r") || (parameters.stop_bundle == "radius")){
	    interpRef = distTarget;
        } else if(parameters.stop_bundle == "plane"){
	    interpRef = kiTarget[0]*finalPoint[0] + kiTarget[1]*finalPoint[1] + kiTarget[2]*finalPoint[2]; 	
        } else {
	    std::cout<<"# WARNING : Wrong stop criterion for integration"<<std::endl;
	    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	    std::terminate();
        }
    }

    // Distance between source and photon at the same comoving radius
    double dist_sep = std::sqrt(pow(trueTarget[0]-finalPoint[0],2)+pow(trueTarget[1]-finalPoint[1],2)+pow(trueTarget[2]-finalPoint[2],2));
    result[1][0] = dist_sep/distTarget;

    // Set angle for re-run rejected sources
    if(parameters.use_previous_catalogues == 2){
	result[0][0] = phi;
	result[0][1] = theta;
    }

    // Compute root-finder is we do not reach the desired accuracy
    if(result[1][0] > parameters.cat_accuracy){
	if(parameters.use_previous_catalogues == 0){ // run catalogues
	    // Compute jacobian for Newton's method
	    // First few rays using the Identity matrix
	    if(iteration < 3){
	        jacobianinv[0][0] = 1;
	        jacobianinv[0][1] = 0;
	        jacobianinv[1][0] = 0;
	        jacobianinv[1][1] = 1;
	    // After a few iterations, use the infinitesimal matrix
	    } else if( (iteration < 5) | ( (parameters.beam == "infinitesimal") | (parameters.beam == "infinitesimal_born") ) ){
	        jacobian = Lensing::dbetadtheta_infinitesimal(distTarget, trajectory, octree, length);
	        jacobianinv = Utility::invMatrix2d(jacobian); 
	    // If it still did not work and we want the jacobian with bundle, then use this method
	    } else{
	        jacobian = Lensing::dbetadtheta(parameters, kiTarget, interpRef, observer, phi, theta, distTarget, cosmology, octree, vobs, length);
	        if(jacobian[0][0] == 42 && jacobian[0][1] == 42){
	            result[0][0] = 42;
	            result[0][1] = 42;
	            result[1][0] = 22;
	            return result;
	        }
	        jacobianinv = Utility::invMatrix2d(jacobian); 
	    }
	// rerun rejected sources
	} else{
	    // Set the size of the grid used to try to find the null geodesic
	    const double err = 5*result[1][0]/(iteration+1);
	    // Store the angles found in previous iteration in temporary variables
	    double phi_tmp = phi;
	    double theta_tmp = theta;
	    // Store the distance between source and photon for the previous angles
	    double sep_tmp = result[1][0];
	    std::cout<<std::setprecision(17)<<"# Direction : "<<phi<<" "<<theta<<" result[1][0] : "<<result[1][0]<<" err : "<<err<<std::endl;
	    // Launch rays toward the direction of the grid, with the previous angle as the center, while size and resolution are functions of 'err'
	    for(double  i = -err; i <= err; i += 0.2*err){
		for(double j = -err; j <= err; j += 0.2*err){
		    magrathea::Evolution<Photon<double, 3> > trajectory_tmp;
		    photon = Integrator::launch(observer[0], observer[1], observer[2], phi + i, theta + j);	
		    trajectory_tmp.append(photon);
		    // Integrate new ray until the comoving radius of the source
		    Integrator::integrate(trajectory_tmp, "radius", distTarget, cosmology, octree, vobs, length, parameters.nsteps);
		    const double marked2 = trajectory_tmp.size()-1;
		    firstid = marked2 - (marked2 > 0);
		    previous = trajectory_tmp[firstid].chi();
		    next = trajectory_tmp[firstid+1].chi();
		    f = (next-distTarget)/(next-previous);
		    finalPoint[0] = trajectory_tmp[firstid].x()*f + trajectory_tmp[firstid+1].x()*(1-f);
		    finalPoint[1] = trajectory_tmp[firstid].y()*f + trajectory_tmp[firstid+1].y()*(1-f);
		    finalPoint[2] = trajectory_tmp[firstid].z()*f + trajectory_tmp[firstid+1].z()*(1-f);
		    // Get new distance between photon and source at the same comoving distance
		    sep_tmp = std::sqrt(pow(trueTarget[0]-finalPoint[0],2)+pow(trueTarget[1]-finalPoint[1],2)+pow(trueTarget[2]-finalPoint[2],2))/distTarget;
		    // If we find one trajectory that gives us the desired accuracy
		    if(sep_tmp < parameters.cat_accuracy){
			std::cout<<"Good, stop now "<<i<<" "<<j<<" delta : "<<err<<" separation : "<<sep_tmp<<" best for the moment : "<<result[1][0]<<" must be below : "<<parameters.cat_accuracy<<std::endl;
			result[1][0] = sep_tmp;
			phi_tmp = phi + i;
			theta_tmp = theta + j;
			trajectory = trajectory_tmp;
			// Prepare for the real jacobian bundle computation
			if(parameters.beam == "bundle"){
			    if(parameters.plane == "sachs"){
			        kiTarget[0] = trajectory[firstid].dxdl()*f + trajectory[firstid+1].dxdl()*(1-f);
			        kiTarget[1] = trajectory[firstid].dydl()*f + trajectory[firstid+1].dydl()*(1-f);
			        kiTarget[2] = trajectory[firstid].dzdl()*f + trajectory[firstid+1].dzdl()*(1-f);	
			    } else if (parameters.plane == "normal"){
			        kiTarget = finalPoint;
			    } else if(parameters.plane == "exact"){
	    			std::cout<<"# Jacobian 'exact' not yet implemented !"<<std::endl;
			    }else{
			        std::cout<<"# WARNING : Wrong plane, please choose 'sachs', 'normal' or 'exact'"<<std::endl;
			        std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
			        std::terminate();
			    }
			    // Interpolation 
			    if(parameters.stop_bundle == "redshift"){
			        interpRef = redshifts[1]; 
			    } else if(parameters.stop_bundle == "a"){
			        interpRef = scale_factor; 	
			    } else if((parameters.stop_bundle == "t") || (parameters.stop_bundle == "eta")){
			        interpRef = trajectory[firstid].t()*f + trajectory[firstid+1].t()*(1-f); 	
			    } else if(parameters.stop_bundle == "lambda"){
			        interpRef = trajectory[firstid].lambda()*f + trajectory[firstid+1].lambda()*(1-f); 	
			    }  else if((parameters.stop_bundle == "r") || (parameters.stop_bundle == "radius")){
			        interpRef = distTarget;
			    } else if(parameters.stop_bundle == "plane"){
			        interpRef = kiTarget[0]*finalPoint[0] + kiTarget[1]*finalPoint[1] + kiTarget[2]*finalPoint[2]; 	
			    }else {
			        std::cout<<"# WARNING : Wrong stop criterion for integration"<<std::endl;
			        std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
			        std::terminate();
			    }
			}
			// Break second loop
			break;
		    // If trajectory if better than previously but still not the desired accuracy
		    } else if (sep_tmp < result[1][0]){
			std::cout<<"# Better, continue "<<i<<" "<<j<<" delta : "<<err<<" separation : "<<sep_tmp<<" best for the moment : "<<result[1][0]<<" must be below : "<<parameters.cat_accuracy<<std::endl;
			result[1][0] = sep_tmp;
			phi_tmp = phi + i;
			theta_tmp = theta + j;
			trajectory = trajectory_tmp;
		    } else{
			// Break second loop
			continue;
		    }
	        }
		// Break first loop
	        if(sep_tmp < parameters.cat_accuracy){ // if good trajectory
		    break;
	        } else{
		    continue;
	        }
	    }
	    // Finalise
	    result[0][0] = phi_tmp;
	    result[0][1] = theta_tmp;
	    // Need to iterate on rejected sources
	    if(result[1][0] > parameters.cat_accuracy)
	        return result;
        }

	// Newton's method for root-finder (since it did not yet converge), using a 2d cartesian plane in the plane-parallel approximation
	if(parameters.use_previous_catalogues == 0){ 
            const double cp(std::cos(phi)), sp(std::sin(phi)), ct(std::cos(theta)), st(std::sin(theta));
            //Initialization plane. 
            Point  e1, e2, sep, pos, targeted;
            e1[0] = -sp;
            e1[1] = cp;
            e1[2] = 0;
            e2[0] = -cp*ct;
            e2[1] = -sp*ct;
            e2[2] = st;
            targeted[0] = distTarget*cp*st;
            targeted[1] = distTarget*sp*st;
            targeted[2] = distTarget*ct;

            for(unsigned int idim = 0; idim < 3; idim++){
                sep[idim] = finalPoint[idim] -trueTarget[idim];
                pos[idim] = targeted[idim] -trueTarget[idim];
            }
            double sepx = sep[0]*e1[0] + sep[1]*e1[1];
            double sepy = sep[0]*e2[0] + sep[1]*e2[1] + sep[2]*e2[2];
            double posx = pos[0]*e1[0] + pos[1]*e1[1];
            double posy = pos[0]*e2[0] + pos[1]*e2[1] + pos[2]*e2[2];
            double sepxnew = jacobianinv[0][0]*sepx + jacobianinv[0][1]*sepy;
            double sepynew = jacobianinv[1][0]*sepx + jacobianinv[1][1]*sepy;

	    double posxnew = posx - sepxnew;
	    double posynew = posy - sepynew;
            double xnew = trueTarget[0]+(posxnew*e1[0]+posynew*e2[0]);
            double ynew = trueTarget[1]+(posxnew*e1[1]+posynew*e2[1]);
            double znew = trueTarget[2]+posynew*e2[2];
            result[0][0] = std::atan2(ynew,xnew);
            result[0][1] = std::acos(znew/std::sqrt(xnew*xnew+ynew*ynew+znew*znew));
	    return result;
	}
    }


    // Compute jacobian
    if(parameters.beam == "bundle"){
	jacobian = Lensing::dbetadtheta(parameters, kiTarget, interpRef, observer, phi, theta, distTarget, cosmology, octree, vobs, length);
    } else if(parameters.beam == "infinitesimal"){
	jacobian = Lensing::dbetadtheta_infinitesimal(distTarget, trajectory, octree, length);
    } else if(parameters.beam == "infinitesimal_born"){
        trajectory_born.append(photon);
        Integrator::integrate<-1>(trajectory_born, "radius", distTarget, cosmology, octree, vobs, length, parameters.nsteps);
	jacobian = Lensing::dbetadtheta_infinitesimal(distTarget, trajectory_born, octree, length);
    } else{
	std::cout<<"# beam must be 'bundle' or 'infinitesimal'"<<std::endl;
	std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
	std::terminate();
    }
    // Check if the jacobian is correct 
    if( jacobian[0][0] == 42 && jacobian[0][1] == 42){
	result[0][0] = 42;
	result[0][1] = 42;
	result[1][0] = 22;
        return result;
    }

    // Compute redshift perturbations
    const double potential = trajectory[firstid].phi()*f + trajectory[firstid+1].phi()*(1-f);
    const double isw =  trajectory[firstid].isw()*f + trajectory[firstid+1].isw()*(1-f);
    // Redshift FLRW
    redshifts[0] = - 1 + 1./(scale_factor);

    // Compute conversion factor from Ramses Units to SI
    const double unit_t = scale_factor*scale_factor*parameters.mpc/(h*std::hecto::num*std::kilo::num);
    const double unit_l = scale_factor*length*100;
    const double velocityx = velocity[0]*unit_l*1e-2/(unit_t*c); // cm -> m
    const double velocityy = velocity[1]*unit_l*1e-2/(unit_t*c); // cm -> m
    const double velocityz = velocity[2]*unit_l*1e-2/(unit_t*c); // cm -> m
    // Compute the components for the exact definition of redshift at the observer
    const double gref = -trajectory.front().a()*c*trajectory[0].dtdl()*(1.+trajectory.front().phi()/c2 
			+ (vobs[0]*trajectory[0].dxdl() + vobs[1]*trajectory[0].dydl() + vobs[2]*trajectory[0].dzdl())/(c*trajectory[0].dtdl())
			+ 0.5*(vobs[0]*vobs[0] + vobs[1]*vobs[1] + vobs[2]*vobs[2]));
    const double dtdl = trajectory[firstid].dtdl()*f + trajectory[firstid+1].dtdl()*(1-f);
    const double dxdl = trajectory[firstid].dxdl()*f + trajectory[firstid+1].dxdl()*(1-f);
    const double dydl = trajectory[firstid].dydl()*f + trajectory[firstid+1].dydl()*(1-f);
    const double dzdl = trajectory[firstid].dzdl()*f + trajectory[firstid+1].dzdl()*(1-f);
    // Compute the components for the exact definition of redshift at the source
    const double gKmuUnu = -scale_factor*c*dtdl*(1. + potential/c2 
			+ (velocityx*dxdl + velocityy*dydl + velocityz*dzdl)/(c*dtdl) 
			+ 0.5*(velocityx*velocityx + velocityy*velocityy + velocityz*velocityz)); 
    // Redshift FLRW + Potential + Doppler
    redshifts[2] = redshifts[1] + (velocityx*dxdl + velocityy*dydl + velocityz*dzdl)/(scale_factor*c*dtdl) 
				       - (vobs[0]*trajectory[0].dxdl() + vobs[1]*trajectory[0].dydl() + vobs[2]*trajectory[0].dzdl())/(scale_factor*c*trajectory[0].dtdl()); 
    // Redshift FLRW + Potential + Doppler + Transverse Doppler
    redshifts[3] = redshifts[2] + 0.5*(velocityx*velocityx + velocityy*velocityy + velocityz*velocityz - vobs[0]*vobs[0] - vobs[1]*vobs[1] - vobs[2]*vobs[2])/scale_factor;
    // Redshift FLRW + Potential + Doppler + Transverse Doppler + ISW/RS
    redshifts[4] = redshifts[3] + isw/scale_factor;
    // Compute the exact definition of redshift in GR (at first order)
    redshifts[5] = gKmuUnu/gref - 1;

    if(parameters.use_previous_catalogues == 0){
	result[0][0] = phi;
	result[0][1] = theta; 
    }

    return result;
}


/// \brief          Iterate Newton Method.
/// \details        Iterate Newton Method to a certain precision 
/// \tparam         Order Octree interpolation order :  0 for NGP, 1 for CIC, 2 for TSC
///                 or -1 for a homogeneous universe.
/// \tparam         RK4 Runge-kutta of fourth order or euler.
/// \tparam         Verbose Verbose mode for debug purposes.
/// \tparam	    Point point type.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Parameter Parameter type.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      observer Observer position.
/// \param[in]      phi Source comoving angular coordinate phi.
/// \param[in]      theta Source comoving angular coordinate theta.
/// \param[in]      target Comoving position of the source
/// \param[in]      velocity Velocity of the target.
/// \param[in,out]  jacobian Jacobian matrix
/// \param[in]      parameters Parameter structure
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      h Dimensionless Hubble parameter
/// \param[in,out]  redshifts Vector containing the redshift decomposition for a given source.
/// \param[in]      interpRef value for interpolation of the bundle.
/// \return         2x2 array with NEW initial angles given after several iterations of newton method
///		    and errors in angles at same radius for the NEW initial angles.
template < int Order, bool RK4, bool Verbose, class Point, class Cosmology, class Octree, class Type, class Parameter > 
std::array< std::array<double,2>, 2 > Catalogues::iterateNewtonMethod(const Point& vobs, const Point& observer, const Type phi, const Type theta, const Point& target, const Point& velocity, std::array< std::array<double, 2>, 2 >& jacobian, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h, std::vector<double>& redshifts, double& interpRef) 
{

    std::array< std::array<double,2>, 2 > result;
    unsigned int iteration(0);
    jacobian = std::array< std::array<double, 2>, 2 >();
    // For regular run, maximum iterations at 10 (or else will be rejected). Increased to 100 to re-run rejected sources with grid method
    const uint nmaxiterations = (parameters.use_previous_catalogues == 0) ? 10 : 100;

    result[0][1] = theta;
    result[0][0] = phi;
    // Call Newton's method routine
    do {
	result = newtonMethod2d(vobs, observer, target, result[0], velocity, jacobian, parameters, cosmology, octree, length, h, redshifts, interpRef, iteration);
	iteration++;
	if(result[0][0] == 42 && result[0][1] == 42)
	    return result;
	if(iteration > nmaxiterations && result[1][0] > parameters.cat_accuracy){
	    result[1][1] = 42;
	    return result;
	}
    } while(result[1][0] > parameters.cat_accuracy);
    result[1][1] = iteration;

    return result; 
}

/// \brief          Give new position of points.
/// \details        Give lensed position of points given points.
/// \tparam	    Point point type.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam         Parameter Parameter type.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      rotm1 rotation matrix
/// \param[in]      filename output file name.
/// \param[in]      observer Observer position.
/// \param[in]      targets_position vector of target positions, velocities and radius.
/// \param[in]      previous_catalogue Previously computed catalogue, used to re-run rejected sources.
/// \param[in]      parameters Parameter structure
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      h Dimensionless Hubble parameter
template < class Point, class Cosmology, class Octree, class Type, class Parameter > 
void Catalogues::relCat(const Point& vobs, const std::array< std::array< double, 3 >, 3 >& rotm1, std::string& filename, const Point& observer, const std::vector < std::array<double, 8> >& targets_position, const std::vector < std::array< double, 18 > >& previous_catalogue, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h) 
{

    const unsigned int size = targets_position.size();
    std::vector<std::array<double, 16> > catalog(size);
    Utility::parallelize(size, [=, &catalog, &vobs, &rotm1, &observer, &targets_position, &previous_catalogue, &parameters, &cosmology, &octree, &length, &h](const uint i){

	Point trueTarget, velocityTarget;
        std::array< std::array<double, 2>, 2 > jacobian;
        std::array<std::array<double,2>, 2> result;
        double interpRef(0), phi(0), theta(0);
	std::vector<double> redshifts(6);
	// Fill position and velocity of sources (may need rotation for narrow cones)
	for(unsigned int j = 0; j < 3; j++){
	    if(parameters.isfullsky == 1){
		trueTarget[j] = targets_position[i][j];
		velocityTarget[j] = targets_position[i][j+3];	
	    }
	    else{
		trueTarget[j] = targets_position[i][0]*rotm1[j][0] +  targets_position[i][1]*rotm1[j][1] + targets_position[i][2]*rotm1[j][2];
	        velocityTarget[j] = targets_position[i][3]*rotm1[j][0] +  targets_position[i][4]*rotm1[j][1] + targets_position[i][5]*rotm1[j][2];;	
	    }
	}
	// If we re-run rejected sources, then start with the latest observed angle computed
	if(previous_catalogue.size() > 0){
	    theta = previous_catalogue[i][4];
	    phi = previous_catalogue[i][3];
	// If we compute the catalogue, the first guess to launch the ray is toward the comoving position of the source
	} else{
	    const double distTarget(std::sqrt(trueTarget[0]*trueTarget[0]+trueTarget[1]*trueTarget[1]+trueTarget[2]*trueTarget[2]));
	    theta = std::acos(trueTarget[2]/distTarget);
	    phi = std::atan2(trueTarget[1],trueTarget[0]);
	}
	result = Catalogues::iterateNewtonMethod(vobs, observer, phi, theta, trueTarget, velocityTarget, jacobian, parameters, cosmology, octree, length, h, redshifts, interpRef);
	if(previous_catalogue.size() > 0){
	    theta = previous_catalogue[i][2];
	    phi = previous_catalogue[i][1];
	} 
	// Put full results in array
	catalog[i][0] = phi;  // Comoving angle (phi)
	catalog[i][1] = theta; // Comoving angle (theta)
	catalog[i][2] = result[0][0]; // Observed angle (phi)
	catalog[i][3] = result[0][1]; // Observed angle (theta)
	catalog[i][4] = result[1][0]; // Error on angle at the source (phi)
	catalog[i][5] = result[1][1]; // Error on angle at the source (theta)
	catalog[i][6] = redshifts[0]; // Redshift FLRW
	catalog[i][7] = redshifts[1]; // Redshift FLRW + Potential
	catalog[i][8] = redshifts[2]; // Redshift FLRW + Potential + Doppler
	catalog[i][9] = redshifts[3]; // Redshift FLRW + Potential + Doppler + Transverse Doppler
	catalog[i][10] = redshifts[4]; // Redshift FLRW + Potential + Doppler + Transverse Doppler + ISW/RS
	catalog[i][11] = redshifts[5]; // Redshift GR (first order in metric perturbations)
	catalog[i][12] = jacobian[0][0]; // Lensing distortion matrix (a11)
	catalog[i][13] = jacobian[0][1]; // Lensing distortion matrix (a12)
	catalog[i][14] = jacobian[1][0]; // Lensing distortion matrix (a21)
	catalog[i][15] = jacobian[1][1]; // Lensing distortion matrix (a22)

    });

    // Output result in ASCII files
    std::string filenameError = Output::name(filename, "_err", ".txt"); 
    std::string filenameRej = Output::name(filename, "_reject", ".txt"); 
    filename = Output::name(filename, ".txt");

    std::ofstream ofst;
    ofst.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc); 
    ofst.close();
    std::ofstream monOutput(filename.c_str(), std::ios::app);

    std::ofstream ofsterr;
    ofsterr.open(filenameError.c_str(), std::ofstream::out | std::ofstream::trunc); 
    ofsterr.close();
    std::ofstream monOutputErr(filenameError.c_str(), std::ios::app);

    std::ofstream ofstrej;
    ofstrej.open(filenameRej.c_str(), std::ofstream::out | std::ofstream::trunc); 
    ofstrej.close();
    std::ofstream monOutputRej(filenameRej.c_str(), std::ios::app);

    	// Ecriture dans un fichier
    for(unsigned int i = 0; i < size; i++){ 
	// If everything is fine
	if(catalog[i][2] != 42 &&  catalog[i][3] != 42 && catalog[i][5] != 42){
	   if(monOutput)
	        monOutput <<std::setprecision(17)<<targets_position[i][6]<<" "<<catalog[i][0]<<" "<<catalog[i][1]<<" "<<catalog[i][2]<<" "<<catalog[i][3]<<" "<<catalog[i][4]<<" "<<catalog[i][5]<<" "<<catalog[i][6]<<" "<<catalog[i][7]<<" "<<catalog[i][8]<<" "<<catalog[i][9]<<" "<<catalog[i][10]<<" "<<catalog[i][11]<<" "<<catalog[i][12]<<" "<<catalog[i][13]<<" "<<catalog[i][14]<<" "<<catalog[i][15]<<" "<<targets_position[i][7]<<std::endl; 
	// If Source is outside of the cone
	} else if(catalog[i][2] == 42 &&  catalog[i][3] == 42){
	   if(monOutputErr)
		monOutputErr <<targets_position[i][6]<<" "<<catalog[i][0]<<" "<<catalog[i][1]<< " "<<catalog[i][2]<<" "<<catalog[i][3]<<" "<<catalog[i][4]<<" "<<catalog[i][4]<<std::endl;
	// If source could not converge (rejected)
	} else{
	   if(monOutputRej)
		monOutputRej <<std::setprecision(17)<<targets_position[i][6]<<" "<<catalog[i][0]<<" "<<catalog[i][1]<<" "<<catalog[i][2]<<" "<<catalog[i][3]<<" "<<catalog[i][4]<<" "<<catalog[i][5]<<" "<<catalog[i][6]<<" "<<catalog[i][7]<<" "<<catalog[i][8]<<" "<<catalog[i][9]<<" "<<catalog[i][10]<<" "<<catalog[i][11]<<" "<<catalog[i][12]<<" "<<catalog[i][13]<<" "<<catalog[i][14]<<" "<<catalog[i][15]<<" "<<targets_position[i][7]<<std::endl; 
	} 	
    } // i 
    // Close streams
    monOutput.close(); 
    monOutputErr.close(); 
    monOutputRej.close(); 
}


/// \brief          Give new position of points.
/// \details        Give lensed position of points given points.
/// \tparam	    Point point type.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam	    Parameter Parameter type.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      filename output file name.
/// \param[in]      observer Observer position.
/// \param[in]      previous_catalogue Previously computed catalogue, used to re-run rejected sources.
/// \param[in]      parameters Parameter structure
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      h Dimensionless Hubble parameter 
template < class Point, class Cosmology, class Octree, class Type, class Parameter > 
void Catalogues::relCat_with_previous_cat(const Point& vobs, std::string& filename, const Point& observer, std::vector < std::array< double, 18 > >& previous_catalogue, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h) 
{

    const unsigned int size = previous_catalogue.size();

    if(size > 0){

        Utility::parallelize(size, [=, &vobs, &observer, &previous_catalogue, &parameters, &cosmology, &octree, &length, &h](const uint i){
            std::array< std::array<double, 2>, 2 > jacobian;
	    magrathea::Evolution<Photon<double, 3> > trajectory, trajectory_born;
	    Photon<double, 3> photon;
	    unsigned int firstid(0);
	    Point kiTarget, finalPoint;
	    double interpRef(0);
	    const double scale_factor = 1./(1.+previous_catalogue[i][7]);
	    // Initialise photon
	    // If Born approximation, then launch toward the comoving position of the source. Otherwise, launch toward the observed position
	    const double phi = (parameters.beam == "infinitesimal_born") ? previous_catalogue[i][1] : previous_catalogue[i][3];
	    const double theta = (parameters.beam == "infinitesimal_born") ? previous_catalogue[i][2] :previous_catalogue[i][4];
	    // Launch photon
	    photon = Integrator::launch(observer[0], observer[1], observer[2], phi, theta);	
	    trajectory.append(photon);
	    // Propagate photon until it reaches the scale factor or the source
	    Integrator::integrate(trajectory, "a", scale_factor, cosmology, octree, vobs, length, parameters.nsteps);

	    const unsigned int marked = trajectory.size()-1;
	    firstid = marked - (marked > 0);
	    const double previous = trajectory[firstid].a();
	    const double next = trajectory[firstid+1].a();
	    double f = (next-scale_factor)/(next-previous);

	    finalPoint[0] = trajectory[firstid].x()*f + trajectory[firstid+1].x()*(1-f);
	    finalPoint[1] = trajectory[firstid].y()*f + trajectory[firstid+1].y()*(1-f);
	    finalPoint[2] = trajectory[firstid].z()*f + trajectory[firstid+1].z()*(1-f);
	    const double distTarget = trajectory[firstid].chi()*f + trajectory[firstid+1].chi()*(1-f);

	    // Compute the lensing distortion matrix
	    if(parameters.beam == "bundle"){
	        if(parameters.plane == "sachs"){
	            kiTarget[0] = trajectory[firstid].dxdl()*f + trajectory[firstid+1].dxdl()*(1-f);
	            kiTarget[1] = trajectory[firstid].dydl()*f + trajectory[firstid+1].dydl()*(1-f);
	            kiTarget[2] = trajectory[firstid].dzdl()*f + trajectory[firstid+1].dzdl()*(1-f);	
	        } else if (parameters.plane == "normal"){
	            kiTarget = finalPoint;
	        } else if (parameters.plane == "exact"){
	    	    std::cout<<"# Jacobian 'exact' not yet implemented !"<<std::endl;
	        } else{
	            std::cout<<"# WARNING : Wrong plane, please choose 'sachs', 'normal' or 'exact'"<<std::endl;
		    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		    std::terminate();
	        }

	        // Interpolation 
	        if(parameters.stop_bundle == "redshift"){
	            interpRef = trajectory[firstid].redshift()*f + trajectory[firstid+1].redshift()*(1-f);
	        } else if(parameters.stop_bundle == "a"){
	            interpRef = scale_factor; 	
	        } else if((parameters.stop_bundle == "t") || (parameters.stop_bundle == "eta")){
	            interpRef = trajectory[firstid].t()*f + trajectory[firstid+1].t()*(1-f); 	
	        } else if(parameters.stop_bundle == "lambda"){
	            interpRef = trajectory[firstid].lambda()*f + trajectory[firstid+1].lambda()*(1-f); 	
	        }  else if((parameters.stop_bundle == "r") || (parameters.stop_bundle == "radius")){
	            interpRef = distTarget;
	        } else if(parameters.stop_bundle == "plane"){
	            interpRef = kiTarget[0]*finalPoint[0] + kiTarget[1]*finalPoint[1] + kiTarget[2]*finalPoint[2]; 	
	        } else {
	            std::cout<<"# WARNING : Wrong stop criterion for integration"<<std::endl;
		    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		    std::terminate();
	        }
	        jacobian = Lensing::dbetadtheta(parameters, kiTarget, interpRef, observer, phi, theta, distTarget, cosmology, octree, vobs, length);
	    } else if(parameters.beam == "infinitesimal"){
	        jacobian = Lensing::dbetadtheta_infinitesimal(distTarget, trajectory, octree, length);
	    } else if(parameters.beam == "infinitesimal_born"){
	        trajectory_born.append(photon);
	        Integrator::integrate<-1>(trajectory_born, "a", scale_factor, cosmology, octree, vobs, length, parameters.nsteps);
	        jacobian = Lensing::dbetadtheta_infinitesimal(distTarget, trajectory_born, octree, length);
	    } else{
	        std::cout<<"# WARNING: beam must be 'bundle' or 'infinitesimal'"<<std::endl;
		std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		std::terminate();
	    }
	    // When re-running a previous catalogue, only modify the distortion matrix
	    previous_catalogue[i][13] = jacobian[0][0];
	    previous_catalogue[i][14] = jacobian[0][1];
	    previous_catalogue[i][15] = jacobian[1][0];
	    previous_catalogue[i][16] = jacobian[1][1];
        });

        // Clear the file
        std::string filenameError = Output::name(filename, "_err", ".txt"); 
        std::string filenameRej = Output::name(filename, "_reject", ".txt"); 
        filename = Output::name(filename, ".txt");

        std::ofstream ofst;
        ofst.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc); 
        ofst.close();
        std::ofstream monOutput(filename.c_str(), std::ios::app);

        std::ofstream ofsterr;
        ofsterr.open(filenameError.c_str(), std::ofstream::out | std::ofstream::trunc); 
        ofsterr.close();
        std::ofstream monOutputErr(filenameError.c_str(), std::ios::app);


    	// Write an ASCII file
        for(unsigned int i = 0; i < size; i++){
	    // If everything is fine, put in the output catalogue 
	    if(previous_catalogue[i][13] != 42 &&  previous_catalogue[i][14] != 42){
	       if(monOutput)
	            monOutput <<std::setprecision(17)<<previous_catalogue[i][0]<<" "<<previous_catalogue[i][1]<<" "<<previous_catalogue[i][2]<<" "<<previous_catalogue[i][3]<<" "<<previous_catalogue[i][4]<<" "<<previous_catalogue[i][5]<<" "<<previous_catalogue[i][6]<<" "<<previous_catalogue[i][7]<<" "<<previous_catalogue[i][8]<<" "<<previous_catalogue[i][9]<<" "<<previous_catalogue[i][10]<<" "<<previous_catalogue[i][11]<<" "<<previous_catalogue[i][12]<<" "<<previous_catalogue[i][13]<<" "<<previous_catalogue[i][14]<<" "<<previous_catalogue[i][15]<<" "<<previous_catalogue[i][16]<<" "<<previous_catalogue[i][17]<<std::endl; 
	    // If there is a problem (for example use of a very lare bundle method at the edge of a cone), then put in the rror file
	    } else{
	       if(monOutputErr)
	        monOutputErr <<std::setprecision(17)<<previous_catalogue[i][0]<<" "<<previous_catalogue[i][1]<<" "<<previous_catalogue[i][2]<<" "<<previous_catalogue[i][3]<<" "<<previous_catalogue[i][4]<<" "<<previous_catalogue[i][5]<<" "<<previous_catalogue[i][6]<<" "<<previous_catalogue[i][7]<<" "<<previous_catalogue[i][8]<<" "<<previous_catalogue[i][9]<<" "<<previous_catalogue[i][10]<<" "<<previous_catalogue[i][11]<<" "<<previous_catalogue[i][12]<<" "<<previous_catalogue[i][13]<<" "<<previous_catalogue[i][14]<<" "<<previous_catalogue[i][15]<<" "<<previous_catalogue[i][16]<<" "<<previous_catalogue[i][17]<<std::endl; 
	    } 
        } 	
	// Close streams
        monOutput.close(); 
        monOutputErr.close(); 
    }
}

/// \brief          Get Hessian lensing matrix
/// \details        Get Hessian lensing matrix for flexion
/// \tparam	    Point point type.
/// \tparam         Cosmology Cosmology evolution type.
/// \tparam         Octree Octree type.
/// \tparam         Type Scalar type.
/// \tparam	    Parameter Parameter type.
/// \param[in]      vobs Observer peculiar velocity, in SI
/// \param[in]      filename output file name.
/// \param[in]      observer Observer position.
/// \param[in]      previous_catalogue Previously computed catalogue, used to re-run rejected sources.
/// \param[in]      parameters Parameter structure
/// \param[in]      cosmology Cosmology evolution.
/// \param[in]      octree Octree.
/// \param[in]      length Spatial length in SI units.
/// \param[in]      h Dimensionless Hubble parameter 
template < class Point, class Cosmology, class Octree, class Type, class Parameter > 
void Catalogues::relCat_with_previous_cat_flexion(const Point& vobs, std::string& filename, const Point& observer, std::vector < std::array< double, 18 > >& previous_catalogue, const Parameter& parameters, const Cosmology& cosmology, const Octree& octree, const Type length, const Type h) 
{

    const unsigned int size = 1;//previous_catalogue.size();

    if(size > 0){

        Utility::parallelize(size, [=, &vobs, &observer, &previous_catalogue, &parameters, &cosmology, &octree, &length, &h](const uint i){
            std::array< std::array<double, 2>, 2 > jacobian;
	    magrathea::Evolution<Photon<double, 3> > trajectory;
	    Photon<double, 3> photon;
	    unsigned int firstid(0);
	    const double one(1);
	    Point kiTarget, finalPoint;
	    double interpRef(0);
	    const double aexp = 1./(1.+previous_catalogue[i][7]);
	    // Initialise photon
	    // If Born approximation, then launch toward the comoving position of the source. Otherwise, launch toward the observed position
	    const double phi = previous_catalogue[i][3];
	    const double theta = previous_catalogue[i][4];
	    // Launch photon
	    photon = Integrator::launch(observer[0], observer[1], observer[2], phi, theta);	
	    trajectory.append(photon);
	    // Propagate photon until it reaches the scale factor or the source
	    Integrator::integrate(trajectory, "a", aexp, cosmology, octree, vobs, length, parameters.nsteps);

	    const unsigned int marked = trajectory.size()-1;
	    firstid = marked - (marked > 0);
	    const double previous = trajectory[firstid].a();
	    const double next = trajectory[firstid+1].a();
	    double f = (next-aexp)/(next-previous);

	    finalPoint[0] = trajectory[firstid].x()*f + trajectory[firstid+1].x()*(1-f);
	    finalPoint[1] = trajectory[firstid].y()*f + trajectory[firstid+1].y()*(1-f);
	    finalPoint[2] = trajectory[firstid].z()*f + trajectory[firstid+1].z()*(1-f);
	    const double k0 = trajectory[firstid].dtdl()*f + trajectory[firstid+1].dtdl()*(1-f);
	    const double kx = trajectory[firstid].dxdl()*f + trajectory[firstid+1].dxdl()*(1-f);
	    const double ky = trajectory[firstid].dydl()*f + trajectory[firstid+1].dydl()*(1-f);
	    const double kz = trajectory[firstid].dzdl()*f + trajectory[firstid+1].dzdl()*(1-f);
	    const double distTarget = trajectory[firstid].chi()*f + trajectory[firstid+1].chi()*(1-f);

	    //std::terminate();
	    // Compute the lensing Hessian matrix
	    if(parameters.beam == "bundle"){
	        if(parameters.plane == "sachs"){
	            kiTarget[0] = trajectory[firstid].dxdl()*f + trajectory[firstid+1].dxdl()*(1-f);
	            kiTarget[1] = trajectory[firstid].dydl()*f + trajectory[firstid+1].dydl()*(1-f);
	            kiTarget[2] = trajectory[firstid].dzdl()*f + trajectory[firstid+1].dzdl()*(1-f);	
	        } else if (parameters.plane == "normal"){
	            kiTarget = finalPoint;
	        } else if (parameters.plane == "exact"){
	    	    std::cout<<"# Jacobian 'exact' not yet implemented !"<<std::endl;
	        } else{
	            std::cout<<"# WARNING : Wrong plane, please choose 'sachs', 'normal' or 'exact'"<<std::endl;
		    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		    std::terminate();
	        }

	        // Interpolation 
	        if(parameters.stop_bundle == "redshift"){
	            interpRef = trajectory[firstid].redshift()*f + trajectory[firstid+1].redshift()*(1-f);
	        } else if(parameters.stop_bundle == "a"){
	            interpRef = aexp; 	
	        } else if((parameters.stop_bundle == "t") || (parameters.stop_bundle == "eta")){
	            interpRef = trajectory[firstid].t()*f + trajectory[firstid+1].t()*(1-f); 	
	        } else if(parameters.stop_bundle == "lambda"){
	            interpRef = trajectory[firstid].lambda()*f + trajectory[firstid+1].lambda()*(1-f); 	
	        }  else if((parameters.stop_bundle == "r") || (parameters.stop_bundle == "radius")){
	            interpRef = distTarget;
	        } else if(parameters.stop_bundle == "plane"){
	            interpRef = kiTarget[0]*finalPoint[0] + kiTarget[1]*finalPoint[1] + kiTarget[2]*finalPoint[2]; 	
	        } else {
	            std::cout<<"# WARNING : Wrong stop criterion for integration"<<std::endl;
		    std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		    std::terminate();
	        }
	        jacobian = Lensing::dbetadtheta(parameters, kiTarget, interpRef, observer, phi, theta, distTarget, cosmology, octree, vobs, length);
	    } else{
	        std::cout<<"# WARNING: beam must be 'bundle' for flexion"<<std::endl;
		std::cout<<"# Error at file "<<__FILE__<<", line : "<<__LINE__<<std::endl;
		std::terminate();
	    }
	    // When re-running a previous catalogue, only modify the distortion matrix
	    previous_catalogue[i][13] = jacobian[0][0];
	    previous_catalogue[i][14] = jacobian[0][1];
	    previous_catalogue[i][15] = jacobian[1][0];
	    previous_catalogue[i][16] = jacobian[1][1];
        });

        // Clear the file
        std::string filenameError = Output::name(filename, "_err", ".txt"); 
        filename = Output::name(filename, ".txt");

        std::ofstream ofst;
        ofst.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc); 
        ofst.close();
        std::ofstream monOutput(filename.c_str(), std::ios::app);

        std::ofstream ofsterr;
        ofsterr.open(filenameError.c_str(), std::ofstream::out | std::ofstream::trunc); 
        ofsterr.close();
        std::ofstream monOutputErr(filenameError.c_str(), std::ios::app);


    	// Write an ASCII file
        for(unsigned int i = 0; i < size; i++){
	    // If everything is fine, put in the output catalogue 
	    if(previous_catalogue[i][13] != 42 &&  previous_catalogue[i][14] != 42){
	       if(monOutput)
	            monOutput <<std::setprecision(17)<<previous_catalogue[i][0]<<" "<<previous_catalogue[i][13]<<" "<<previous_catalogue[i][14]<<" "<<previous_catalogue[i][15]<<" "<<previous_catalogue[i][16]<<std::endl; 
	    // If there is a problem (for example use of a very lare bundle method at the edge of a cone), then put in the rror file
	    } else{
	       if(monOutputErr)
	        monOutputErr <<std::setprecision(17)<<previous_catalogue[i][0]<<" "<<previous_catalogue[i][13]<<" "<<previous_catalogue[i][14]<<" "<<previous_catalogue[i][15]<<" "<<previous_catalogue[i][16]<<std::endl; 
	    } 
        } 	
	// Close streams
        monOutput.close(); 
        monOutputErr.close(); 
    }
}


#endif // CATALOGUES_H_INCLUDED
