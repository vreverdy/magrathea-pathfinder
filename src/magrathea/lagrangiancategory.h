/* *************************** LAGRANGIANCATEGORY *************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA 
// TITLE :          LagrangianCategory
// DESCRIPTION :    Category concept of lagrangian : data moving with the flow
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com) 
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)] 
// LICENSE :        CECILL-B License 
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           lagrangiancategory.h
/// \brief          Category concept of lagrangian : data moving with the flow
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef LAGRANGIANCATEGORY_H_INCLUDED
#define LAGRANGIANCATEGORY_H_INCLUDED
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
// Category concept of lagrangian : data moving with the flow
/// \brief          Category concept of lagrangian : data moving with the flow.
/// \details        Categorizes data as an lagrangian one. A lagrangian data is 
///                 moving with the flow : this is a way of looking at fluid 
///                 motion where the observer follows an individual fluid parcel 
///                 as it moves through space and time.
class LagrangianCategory final
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
/// \details        Tests and demonstrates the use of LagrangianCategory.
/// \return         0 if no error.
int LagrangianCategory::example()
{
    // Initialize
    std::cout<<"BEGIN = LagrangianCategory::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    LagrangianCategory category;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"LagrangianCategory() : "                      ; LagrangianCategory(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"category = LagrangianCategory() : "           ; category = LagrangianCategory(); std::cout<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = LagrangianCategory::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // LAGRANGIANCATEGORY_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
