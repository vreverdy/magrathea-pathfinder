/* **************************** EULERIANCATEGORY **************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          EulerianCategory
// DESCRIPTION :    Category concept of eulerian : data at a fixed position
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           euleriancategory.h
/// \brief          Category concept of eulerian : data at a fixed position
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef EULERIANCATEGORY_H_INCLUDED
#define EULERIANCATEGORY_H_INCLUDED
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
// Category concept of eulerian : data at a fixed position
/// \brief          Category concept of eulerian : data at a fixed position.
/// \details        Categorizes data as an eulerian one. An eulerian data is 
///                 specified at a fixed location : this is a way of looking 
///                 at a fluid motion that focuses on specific locations in the 
///                 space through which the fluid flows as time passes.
class EulerianCategory final
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
/// \details        Tests and demonstrates the use of EulerianCategory.
/// \return         0 if no error.
int EulerianCategory::example()
{
    // Initialize
    std::cout<<"BEGIN = EulerianCategory::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    EulerianCategory category;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"EulerianCategory() : "                        ; EulerianCategory(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"category = EulerianCategory() : "             ; category = EulerianCategory(); std::cout<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = EulerianCategory::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // EULERIANCATEGORY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
