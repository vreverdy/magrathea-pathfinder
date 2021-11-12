/* ********************************** OBSERVER_VELOCITY ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Observer_velocity
// DESCRIPTION :    Some observer_velocity function
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           observer_velocity.h
/// \brief          Some observer_velocity function
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef OBSERVER_VELOCITY_H_INCLUDED
#define OBSEREVR_VELOCITY_H_INCLUDED


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
     uint ncoarse;
     std::string partdir;
     std::string velocity_field_v0;
} parameters;


class Observer_velocity {

	//Methodes
	public:
	
	// Read parameter file
	template <class Parameters, class Map> static void ReadParamFile(Parameters& parameters, Map& parameter);

	// Velocity field octree
	template < typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > static void CreateOctreeVelocityWithCIC(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, std::vector<Type1>& pos_part, std::vector<Type1>& vel_part);
	template < typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > static void CreateOctreeVelocityWithTSC(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, std::vector<Type1>& pos_part, std::vector<Type1>& vel_part);





};


// Read parameter file
/// \brief          Read parameter file.
/// \details        Read and put in a structure the parameters.
/// \tparam         Parameters structure type
/// \tparam         Map map type
/// \param[in,out]  parameters Structure containing the parameters.
/// \param[in]      parameter Contains parameters to be rewritten
template <class Parameters, class Map> 
void Observer_velocity::ReadParamFile(Parameters& parameters, Map& parameter){
    parameters.ncoarse = std::stoul(parameter["ncoarse"]);
    parameters.partdir = parameter["partdir"];
    parameters.velocity_field_v0 = parameter["velocity_field_v0"];
}

// Compute CIC velocity field
/// \brief          Add velocity to Octree using CIC
/// \details        Add velocity to Octree using CIC
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
/// \param[in,out]  octree Octree be to filled with velocity field  
/// \param[in]      pos_part Position of particles
/// \param[in]      vel_part Velocity of particles
template < typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > 
void Observer_velocity::CreateOctreeVelocityWithCIC(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, std::vector<Type1>& pos_part, std::vector<Type1>& vel_part){

    // Get levels at which we wish to compute the velocity field
    const unsigned int lvlmax = (std::get<0>(*std::max_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level()); 
    const unsigned int lvlmin = (std::get<0>(*std::min_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level());

    // Loop over levels
    for(uint ilvl = lvlmin; ilvl <= lvlmax; ilvl++){
        const double invextension = 1/(EXTENT/pow(2,ilvl));  
        const double half = 0.5*EXTENT/pow(2,ilvl); 
#ifdef VERBOSE
        std::cout<<"# Compute velocity field at level "<<ilvl<<std::endl;
#endif
	// Loop over all the particles
    	Utility::parallelize(pos_part.size()/3, [=, &octree, &pos_part, &vel_part, &invextension, &half](const uint i){
            Index idxvertex;
            Data data; 
            unsigned long long int marker(0);
            double vratio(0);
	    // Loop over the 8 neighboring cells of a given particle
            for(int iz = -1; iz <= 1; iz += 2){
	        for(int iy = -1; iy <= 1; iy += 2){
	            for(int ix = -1; ix <= 1; ix += 2){
		        // Create an index at the level of interest for the neighboring cell
		        idxvertex = idxvertex.template compute<Type,Position,Extent>(ilvl, pos_part[3*i] + half*ix, pos_part[3*i+1] + half*iy, pos_part[3*i+2] + half*iz);
		        // Given an index in the octree that is consistent with the created index
       		        marker = std::distance(std::begin(octree), std::upper_bound(std::begin(octree), std::end(octree), Element(idxvertex, data), [](const Element& first, const Element& second){return std::get<0>(first) < std::get<0>(second);}));
		        // If the index exists in the octree, compute CIC
		        if(std::get<0>(*(std::begin(octree)+marker-(marker > 0))) == idxvertex){
		            vratio = (1-std::abs(pos_part[3*i]-idxvertex.template center<Type,Position,Extent>(0))*invextension)*(1-std::abs(pos_part[3*i+1]-idxvertex.template center<Type,Position,Extent>(1))*invextension)*(1-std::abs(pos_part[3*i+2]-idxvertex.template center<Type,Position,Extent>(2))*invextension);
		            std::get<1>(octree[marker - (marker > 0)]).rho() += vratio;
		            std::get<1>(octree[marker - (marker > 0)]).dphidx() += vel_part[3*i]*vratio;
		            std::get<1>(octree[marker - (marker > 0)]).dphidy() += vel_part[3*i+1]*vratio;
		            std::get<1>(octree[marker - (marker > 0)]).dphidz() += vel_part[3*i+2]*vratio;
		        }		
	            }
	        }  
            }
	});
    }
}

// Compute TSC velocity field
/// \brief          Add velocity to Octree using TSC
/// \details        Add velocity to Octree using TSC
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
/// \param[in,out]  octree Octree be to filled with velocity field  
/// \param[in]      pos_part Position of particles
/// \param[in]      vel_part Velocity of particles
template < typename Type1, template <typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container> class Octree, typename Type, class Index, class Data, unsigned int Dimension, class Position, class Extent, class Element, class Container > 
void Observer_velocity::CreateOctreeVelocityWithTSC(Octree< Type, Index, Data, Dimension, Position, Extent, Element, Container >& octree, std::vector<Type1>& pos_part, std::vector<Type1>& vel_part){

    // Get levels at which we wish to compute the velocity field
    const unsigned int lvlmax = (std::get<0>(*std::max_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level()); 
    const unsigned int lvlmin = (std::get<0>(*std::min_element(std::begin(octree), std::end(octree), [](const Element& x, const Element& y){return std::get<0>(x).level() < std::get<0>(y).level();})).level());

    // Loop over levels
    for(uint ilvl = lvlmin; ilvl <= lvlmax; ilvl++){
        const double half = 0.5*EXTENT/pow(2,ilvl); 
        const double twohalves = 2*half; 
#ifdef VERBOSE
        std::cout<<"# Compute velocity field at level "<<ilvl<<std::endl;
#endif
	// Loop over all the particles
    	Utility::parallelize(pos_part.size()/3, [=, &octree, &pos_part, &vel_part, &half, &twohalves](const uint i){
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
		        idxvertex = idxvertex.template compute<Type,Position,Extent>(ilvl, pos_part[3*i] + twohalves*ix, pos_part[3*i+1] + twohalves*iy, pos_part[3*i+2] + twohalves*iz);
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
		            std::get<1>(octree[marker - (marker > 0)]).rho() += vratio;
		            std::get<1>(octree[marker - (marker > 0)]).dphidx() += vel_part[3*i]*vratio;
		            std::get<1>(octree[marker - (marker > 0)]).dphidy() += vel_part[3*i+1]*vratio;
		            std::get<1>(octree[marker - (marker > 0)]).dphidz() += vel_part[3*i+2]*vratio;
		        }		
	            }
	        }  
            }
	});
    }
}


#endif // OBSERVER_VELOCITY_H_INCLUDED
