/* ********************************** SHAPE ********************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Shape
// DESCRIPTION :    Basic implementation of shape
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           shape.h
/// \brief          Basic implementation of shape
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef SHAPE_H_INCLUDED
#define SHAPE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
// Include libs
// Include project
#include "abstractshape.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Basic implementation of shape
/// \brief          Basic implementation of shape.
/// \details        This class is the direct derivation of AbstractShape. 
///                 It provides the most basic and generic shape object 
///                 without adding new functionalities to the abstract class.
class Shape final
: public AbstractShape
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
/// \details        Tests and demonstrates the use of Shape.
/// \return         0 if no error.
int Shape::example()
{
    // Initialize
    std::cout<<"BEGIN = Shape::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    const unsigned int width = 40;
    
    // Construction
    Shape shape;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"Shape() : "                                   ; Shape(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"shape = Shape() : "                           ; shape = Shape(); std::cout<<std::endl;

    // Mathematical functions
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Mathematical functions : "                                                          <<std::endl;    
    std::cout<<std::setw(width*2)<<"shape.pow<0>(42) : "                                                                <<shape.pow<0>(42)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.pow<4>(42) : "                                                                <<shape.pow<4>(42)<<std::endl;    
    std::cout<<std::setw(width*2)<<"shape.pow<4>(42.) : "                                                               <<shape.pow<4>(42.)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.pow<-4>(42.) : "                                                              <<shape.pow<-4>(42.)<<std::endl;    
    std::cout<<std::setw(width*2)<<"shape.rt<1>(42) : "                                                                 <<shape.rt<1>(42)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.rt<2>(42) : "                                                                 <<shape.rt<2>(42)<<std::endl;    
    std::cout<<std::setw(width*2)<<"shape.rt<3>(42.) : "                                                                <<shape.rt<3>(42.)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.rt<4>(42.) : "                                                                <<shape.rt<4>(42.)<<std::endl;    
    std::cout<<std::setw(width*2)<<"shape.rt<-1>(42) : "                                                                <<shape.rt<-1>(42)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.rt<-2>(42) : "                                                                <<shape.rt<-2>(42)<<std::endl;    
    std::cout<<std::setw(width*2)<<"shape.rt<-3>(42.) : "                                                               <<shape.rt<-3>(42.)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.rt<-4>(42.) : "                                                               <<shape.rt<-4>(42.)<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<0>() : "                                                            <<shape.factorial<0>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<1>() : "                                                            <<shape.factorial<1>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<4>() : "                                                            <<shape.factorial<4>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<8>() : "                                                            <<shape.factorial<8>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<20>() : "                                                           <<shape.factorial<20>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<0, false>() : "                                                     <<shape.factorial<0, false>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<1, false>() : "                                                     <<shape.factorial<1, false>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<8, false>() : "                                                     <<shape.factorial<8, false>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<15, false>() : "                                                    <<shape.factorial<15, false>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.factorial<0, true>() : "                                                      <<shape.factorial<0, true>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<1, true>() : "                                                      <<shape.factorial<1, true>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<8, true>() : "                                                      <<shape.factorial<8, true>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.factorial<15, true>() : "                                                     <<shape.factorial<15, true>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.binomial<0, 0, unsigned int>() : "                                            <<shape.binomial<0, 0, unsigned int>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.binomial<15, 8>() : "                                                         <<shape.binomial<15, 8>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.binomial<15, -8>() : "                                                        <<shape.binomial<15, -8>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.binomial<-15, 8>() : "                                                        <<shape.binomial<-15, 8>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.binomial<-15, -8>() : "                                                       <<shape.binomial<-15, -8>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.combinations<0, 0, false, unsigned int>() : "                                 <<shape.combinations<0, 0, false, unsigned int>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.combinations<15, 8, false>() : "                                              <<shape.combinations<15, 8, false>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.combinations<0, 0, true, unsigned int>() : "                                  <<shape.combinations<0, 0, true, unsigned int>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.combinations<15, 8, true>() : "                                               <<shape.combinations<15, 8, true>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.variations<0, 0, false, unsigned int>() : "                                   <<shape.variations<0, 0, false, unsigned int>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.variations<15, 8, false>() : "                                                <<shape.variations<15, 8, false>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.variations<0, 0, true, unsigned int>() : "                                    <<shape.variations<0, 0, true, unsigned int>()<<std::endl; 
    std::cout<<std::setw(width*2)<<"shape.variations<15, 8, true>() : "                                                 <<shape.variations<15, 8, true>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.alt<8>() : "                                                                  <<shape.alt<8>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.alt<15>() : "                                                                 <<shape.alt<15>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.alt<8, 0, 1>() : "                                                            <<shape.alt<8, 0, 1>()<<std::endl;
    std::cout<<std::setw(width*2)<<"shape.alt<15, 0, 1>() : "                                                           <<shape.alt<15, 0, 1>()<<std::endl;

    // Constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Constants : "                                 <<std::endl;    
    std::cout<<std::setw(width)<<"shape.golden() : "                            <<shape.golden()<<std::endl;   
    std::cout<<std::setw(width)<<"shape.golden<2>() : "                         <<shape.golden<2>()<<std::endl;  
    std::cout<<std::setw(width)<<"shape.pi() : "                                <<shape.pi()<<std::endl;    
    std::cout<<std::setw(width)<<"shape.pi<2>() : "                             <<shape.pi<2>()<<std::endl;
    std::cout<<std::setw(width)<<"shape.sqrtpi() : "                            <<shape.sqrtpi()<<std::endl;    
    std::cout<<std::setw(width)<<"shape.sqrtpi<2>() : "                         <<shape.sqrtpi<2>()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Shape::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // SHAPE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
