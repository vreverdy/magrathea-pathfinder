/* ****************************** GRIDCATEGORY ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          GridCategory
// DESCRIPTION :    Category concept of grid : data related to the mesh
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           gridcategory.h
/// \brief          Category concept of grid : data related to the mesh
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef GRIDCATEGORY_H_INCLUDED
#define GRIDCATEGORY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
// Include libs
// Include project
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Category concept of grid : data related to the mesh
/// \brief          Category concept of grid : data related to the mesh.
/// \details        Categorizes data as a grid one. A grid data is a data 
///                 involved in the navigation process through the mesh applied
///                 to a numerical simulation.
class GridCategory final
{
    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of GridCategory.
/// \return         0 if no error.
int GridCategory::example()
{
    // Initialize
    std::cout<<"BEGIN = GridCategory::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    GridCategory category;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"GridCategory() : "                            ; GridCategory(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"category = GridCategory() : "                 ; category = GridCategory(); std::cout<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = GridCategory::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // GRIDCATEGORY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
