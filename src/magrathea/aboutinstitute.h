/* ***************************** ABOUTINSTITUTE ***************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          AboutInstitute
// DESCRIPTION :    Information about an institution or an organization
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           aboutinstitute.h
/// \brief          Information about an institution or an organization
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABOUTINSTITUTE_H_INCLUDED
#define ABOUTINSTITUTE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
// Include libs
// Include project
#include "abstractaboutobject.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Information about an institution or an organization
/// \brief          Information about an institution or an organization.
/// \details        Name, link, contact, address... of an institute.
class AboutInstitute final
: public AbstractAboutObject<AboutInstitute, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>
{
    // Setup
    public: using AbstractAboutObject<AboutInstitute, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>::operator=;
    
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        template <class... Misc> explicit inline AboutInstitute(Misc&&... misc);
    //@}
    
    // Accessors
    /// \name           Accessors 
    //@{
    public:
        inline std::string identification(const std::string& separator = " - ") const;
        inline std::string address(const std::string& separator = ", ") const;
    //@}

    // Mutators
    /// \name           Mutators 
    //@{
    public:
        inline AboutInstitute& identification(const std::string& stitle, const std::string& sname, const std::string& sextended);
        inline AboutInstitute& address(const std::string& sstreet, const std::string& szip, const std::string& scity, const std::string& sregion, const std::string& scountry);
    //@}
    
    // Getters
    /// \name           Getters 
    //@{
    public:
        inline const std::string& title() const;
        inline const std::string& name() const;
        inline const std::string& extended() const;
        inline const std::string& street() const;
        inline const std::string& zip() const;
        inline const std::string& city() const;
        inline const std::string& region() const;
        inline const std::string& country() const;
        inline const std::string& link() const;
        inline const std::string& contact() const;
    //@}
    
    // Setters
    /// \name           Setters 
    //@{
    public:
        inline AboutInstitute& title(const std::string& value);
        inline AboutInstitute& name(const std::string& value);
        inline AboutInstitute& extended(const std::string& value);
        inline AboutInstitute& street(const std::string& value);
        inline AboutInstitute& zip(const std::string& value);
        inline AboutInstitute& city(const std::string& value);
        inline AboutInstitute& region(const std::string& value);
        inline AboutInstitute& country(const std::string& value);
        inline AboutInstitute& link(const std::string& value);
        inline AboutInstitute& contact(const std::string& value);
    //@}

    // Predefined
    /// \name           Predefined 
    //@{
    public:
        static inline const AboutInstitute& luth();
        static inline const AboutInstitute& sap();
        static inline const AboutInstitute& ipag();
        static inline const AboutInstitute& obspm();
        static inline const AboutInstitute& iap();
        static inline const AboutInstitute& ias();
        static inline const AboutInstitute& cnrs();
        static inline const AboutInstitute& cnes();
        static inline const AboutInstitute& cea();
        static inline const AboutInstitute& esa();
    //@}
    
    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------------- LIFECYCLE -------------------------------- //
// Explicit generic constructor
/// \brief          Explicit generic constructor.
/// \details        Provides a generic interface to all constructors of the 
///                 base class.
/// \tparam         Misc (Miscellaneous types.)
/// \param[in]      misc Miscellaneous arguments.  
template <class... Misc> 
inline AboutInstitute::AboutInstitute(Misc&&... misc)
: AbstractAboutObject<AboutInstitute, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string>(std::forward<Misc>(misc)...)
{
    ;
}
// -------------------------------------------------------------------------- //



// -------------------------------- ACCESSORS ------------------------------- //
// Full identification accessor
/// \brief          Full identification accessor.
/// \details        Gets the value of the title, the name and the extended name 
///                 properties. If one of the field is empty, the separator is 
///                 not used.
/// \param[in]      separator String separator.
/// \return         Indentification.
inline std::string AboutInstitute::identification(const std::string& separator) const
{
    std::string result = std::get<0>(_data);
    result += ((result.empty()) || (std::get<1>(_data).empty())) ? (std::get<1>(_data)) : (separator+std::get<1>(_data));
    result += ((result.empty()) || (std::get<2>(_data).empty())) ? (std::get<2>(_data)) : (separator+std::get<2>(_data));
    return result;
}

// Address accessor
/// \brief          Address accessor.
/// \details        Gets the value of the street, zip, city, region and country 
///                 properties. If one of the field is empty, the separator is 
///                 not used.
/// \param[in]      separator String separator.
/// \return         Address.
inline std::string AboutInstitute::address(const std::string& separator) const
{
    std::string result = std::get<3>(_data);
    result += ((result.empty()) || (std::get<4>(_data).empty())) ? (std::get<4>(_data)) : (separator+std::get<4>(_data));
    result += ((result.empty()) || (std::get<5>(_data).empty())) ? (std::get<5>(_data)) : (separator+std::get<5>(_data));
    result += ((result.empty()) || (std::get<6>(_data).empty())) ? (std::get<6>(_data)) : (separator+std::get<6>(_data));
    result += ((result.empty()) || (std::get<7>(_data).empty())) ? (std::get<7>(_data)) : (separator+std::get<7>(_data));
    return result;
}
// -------------------------------------------------------------------------- //



// -------------------------------- MUTATORS -------------------------------- //
// Full identification mutator
/// \brief          Full identification mutator.
/// \details        Sets the value of the title, the name and the extended name 
///                 properties.
/// \param[in]      stitle Title or acronym.
/// \param[in]      sname Complete name.
/// \param[in]      sextended Extended name or information.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::identification(const std::string& stitle, const std::string& sname, const std::string& sextended)
{
    std::get<0>(_data) = stitle;
    std::get<1>(_data) = sname;
    std::get<2>(_data) = sextended;
    return *this;
}

// Address mutator
/// \brief          Address mutator.
/// \details        Sets the value of the street, zip, city, region and country 
///                 properties.
/// \param[in]      sstreet Street number.
/// \param[in]      szip Zip code.
/// \param[in]      scity City.
/// \param[in]      sregion Region.
/// \param[in]      scountry Country.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::address(const std::string& sstreet, const std::string& szip, const std::string& scity, const std::string& sregion, const std::string& scountry)
{
    std::get<3>(_data) = sstreet;
    std::get<4>(_data) = szip;
    std::get<5>(_data) = scity;
    std::get<6>(_data) = sregion;
    std::get<7>(_data) = scountry;
    return *this;
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// Title or acronym getter
/// \brief          Title or acronym getter.
/// \details        Gets the value of the title property.
/// \return         Title.
inline const std::string& AboutInstitute::title() const
{
    return std::get<0>(_data);
}

// Complete name getter
/// \brief          Complete name getter.
/// \details        Gets the value of the name property.
/// \return         Name.
inline const std::string& AboutInstitute::name() const
{
    return std::get<1>(_data);
}

// Extended name or information getter
/// \brief          Extended name or information getter.
/// \details        Gets the value of the extended property.
/// \return         Extended.
inline const std::string& AboutInstitute::extended() const
{
    return std::get<2>(_data);
}

// Street number getter
/// \brief          Street number getter.
/// \details        Gets the value of the street property.
/// \return         Street.
inline const std::string& AboutInstitute::street() const
{
    return std::get<3>(_data);
}

// Zip code getter
/// \brief          Zip code getter.
/// \details        Gets the value of the zip property.
/// \return         Zip.
inline const std::string& AboutInstitute::zip() const
{
    return std::get<4>(_data);
}

// City getter
/// \brief          City getter.
/// \details        Gets the value of the city property.
/// \return         City.
inline const std::string& AboutInstitute::city() const
{
    return std::get<5>(_data);
}

// Region getter
/// \brief          Region getter.
/// \details        Gets the value of the region property.
/// \return         Region.
inline const std::string& AboutInstitute::region() const
{
    return std::get<6>(_data);
}

// Country getter
/// \brief          Country getter.
/// \details        Gets the value of the country property.
/// \return         Country.
inline const std::string& AboutInstitute::country() const
{
    return std::get<7>(_data);
}

// Website or link getter
/// \brief          Website or link getter.
/// \details        Gets the value of the link property.
/// \return         Link.
inline const std::string& AboutInstitute::link() const
{
    return std::get<8>(_data);
}

// Contact information getter
/// \brief          Contact information getter.
/// \details        Gets the value of the contact property.
/// \return         Contact.
inline const std::string& AboutInstitute::contact() const
{
    return std::get<9>(_data);
}
// -------------------------------------------------------------------------- //



// --------------------------------- SETTERS -------------------------------- //
// Title or acronym setter
/// \brief          Title or acronym setter.
/// \details        Sets the value of the title property.
/// \param[in]      value Title.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::title(const std::string& value)
{
    std::get<0>(_data) = value;
    return *this;
}

// Complete name setter
/// \brief          Complete name setter.
/// \details        Sets the value of the name property.
/// \param[in]      value Name.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::name(const std::string& value)
{
    std::get<1>(_data) = value;
    return *this;
}

// Extended name or information setter
/// \brief          Extended name or information setter.
/// \details        Sets the value of the extended property.
/// \param[in]      value Extended.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::extended(const std::string& value)
{
    std::get<2>(_data) = value;
    return *this;
}

// Street number setter
/// \brief          Street number setter.
/// \details        Sets the value of the street property.
/// \param[in]      value Street.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::street(const std::string& value)
{
    std::get<3>(_data) = value;
    return *this;
}

// Zip code setter
/// \brief          Zip code setter.
/// \details        Sets the value of the zip property.
/// \param[in]      value Zip.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::zip(const std::string& value)
{
    std::get<4>(_data) = value;
    return *this;
}

// City setter
/// \brief          City setter.
/// \details        Sets the value of the city property.
/// \param[in]      value City.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::city(const std::string& value)
{
    std::get<5>(_data) = value;
    return *this;
}

// Region setter
/// \brief          Region setter.
/// \details        Sets the value of the region property.
/// \param[in]      value Region.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::region(const std::string& value)
{
    std::get<6>(_data) = value;
    return *this;
}

// Country setter
/// \brief          Country setter.
/// \details        Sets the value of the country property.
/// \param[in]      value Country.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::country(const std::string& value)
{
    std::get<7>(_data) = value;
    return *this;
}

// Website or link setter
/// \brief          Website or link setter.
/// \details        Sets the value of the link property.
/// \param[in]      value Link.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::link(const std::string& value)
{
    std::get<8>(_data) = value;
    return *this;
}

// Contact information setter
/// \brief          Contact information setter.
/// \details        Sets the value of the contact property.
/// \param[in]      value Contact.
/// \return         Self reference.
inline AboutInstitute& AboutInstitute::contact(const std::string& value)
{
    std::get<9>(_data) = value;
    return *this;
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Predefined LUTH
/// \brief          Predefined LUTH.
/// \details        Information on the Laboratoire Univers et Theories.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::luth()
{
    static const AboutInstitute singleton("LUTH", "Laboratoire Univers et Theories", "UMR 8102", "Observatoire de Paris, 5 place Jules Janssen", "92195", "Meudon", "Hauts-de-Seine", "France", "http://luth.obspm.fr", "");
    return singleton;
}

// Predefined SAp
/// \brief          Predefined SAp.
/// \details        Information on the Service d'Astrophysique du CEA Saclay.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::sap()
{
    static const AboutInstitute singleton("SAp", "Service d'Astrophysique du CEA Saclay", "IRFU", "Centre d'etudes de Saclay, Orme des Merisiers, batiment 709", "91191", "Gif-sur-Yvette", "Essonne", "France", "http://irfu.cea.fr/Sap", "astrophysique@cea.fr");
    return singleton;
}

// Predefined IPAG
/// \brief          Predefined IPAG.
/// \details        Information on the Institut de Planetologie et 
///                 d’Astrophysique de Grenoble
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::ipag()
{
    static const AboutInstitute singleton("IPAG", "Institut de Planetologie et d’Astrophysique de Grenoble", "UMR 5274", "414 rue de la Piscine, Domaine Universitaire", "38400", "St-Martin d’Heres", "Isere", "France", "http://ipag.osug.fr", "");
    return singleton;
}

// Predefined OBSPM
/// \brief          Predefined OBSPM.
/// \details        Information on the Observatoire de Paris.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::obspm()
{
    static const AboutInstitute singleton("OBSPM", "Observatoire de Paris", "", "77 avenue Denfert-Rochereau", "75014", "Paris", "Paris", "France", "http://www.obspm.fr", "");
    return singleton;
}

// Predefined IAP
/// \brief          Predefined IAP.
/// \details        Information on the Institut d'Astrophysique de Paris.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::iap()
{
    static const AboutInstitute singleton("IAP", "Institut d'Astrophysique de Paris", "UMR7095", "98 bis boulevard Arago", "75014", "Paris", "Paris", "France", "http://www.iap.fr", "");
    return singleton;
}

// Predefined IAS
/// \brief          Predefined IAS.
/// \details        Information on the Institut d'Astrophysique Spatiale.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::ias()
{
    static const AboutInstitute singleton("IAS", "Institut d'Astrophysique Spatiale", "UMR8617", "Centre universitaire d’Orsay, batiment 120 – 121", "91405 ", "Orsay", "Essonne", "France", "http://www.ias.u-psud.fr", "");
    return singleton;
}

// Predefined CNRS
/// \brief          Predefined CNRS.
/// \details        Information on the Centre National de la Recherche 
///                 Scientifique.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::cnrs()
{
    static const AboutInstitute singleton("CNRS", "Centre National de la Recherche Scientifique", "", "3 rue Michel-Ange", "75794", "Paris cedex 16", "Paris", "France", "http://www.cnrs.fr", "");
    return singleton;
}

// Predefined CNES
/// \brief          Predefined CNES.
/// \details        Information on the Centre National d'Etudes Spatiales.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::cnes()
{
    static const AboutInstitute singleton("CNES", "Centre National d'Etudes Spatiales", "", "2 place Maurice Quentin", "75039", "Paris cedex 01", "Paris", "France", "http://www.cnes.fr", "");
    return singleton;
}

// Predefined CEA
/// \brief          Predefined CEA.
/// \details        Information on the Commissariat a l'Energie Atomique et aux 
///                 Energies Alternatives.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::cea()
{
    static const AboutInstitute singleton("CEA", "Commissariat a l'Energie Atomique et aux Energies Alternatives", "", "", "", "", "", "France", "http://www.cea.fr", "");
    return singleton;
}

// Predefined ESA
/// \brief          Predefined ESA.
/// \details        Information on the European Space Agency.
/// \return         Immutable reference to institute singleton.
inline const AboutInstitute& AboutInstitute::esa()
{
    static const AboutInstitute singleton("ESA", "European Space Agency", "", "", "", "", "", "Europe", "http://www.esa.int", "");
    return singleton;
}
// -------------------------------------------------------------------------- //   



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of AboutInstitute.
/// \return         0 if no error.
int AboutInstitute::example()
{
    // Initialize
    std::cout<<"BEGIN = AboutInstitute::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    AboutInstitute institute("title", "name", "extended", "street", "zip", "city", "region", "country", "link", "contact");

    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"AboutInstitute() : "                          <<AboutInstitute()<<std::endl;
    std::cout<<std::setw(width)<<"AboutInstitute(\"institute\") : "             <<AboutInstitute("institute")<<std::endl;
    std::cout<<std::setw(width)<<"AboutInstitute(institute) : "                 <<AboutInstitute(institute)<<std::endl; 

    // Accessors
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Accessors : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"institute.identification() : "                <<institute.identification()<<std::endl;
    std::cout<<std::setw(width)<<"institute.address() : "                       <<institute.address()<<std::endl;
    
    // Mutators
    std::cout<<std::endl;
    std::cout<<std::setw(width*2)<<"Mutators : "                                                                        <<std::endl;
    std::cout<<std::setw(width*2)<<"institute.identification(\"t\", \"n\", \"e\") : "                                   <<institute.identification("t", "n", "e")<<std::endl;
    std::cout<<std::setw(width*2)<<"institute.address(\"s\", \"z\", \"c\", \"r\", \"c\") : "                            <<institute.address("s", "z", "c", "r", "c")<<std::endl;
    
    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"institute.title() : "                         <<institute.title()<<std::endl;
    std::cout<<std::setw(width)<<"institute.name() : "                          <<institute.name()<<std::endl;
    std::cout<<std::setw(width)<<"institute.extended() : "                      <<institute.extended()<<std::endl;
    std::cout<<std::setw(width)<<"institute.street() : "                        <<institute.street()<<std::endl;
    std::cout<<std::setw(width)<<"institute.zip() : "                           <<institute.zip()<<std::endl;
    std::cout<<std::setw(width)<<"institute.city() : "                          <<institute.city()<<std::endl;
    std::cout<<std::setw(width)<<"institute.region() : "                        <<institute.region()<<std::endl;
    std::cout<<std::setw(width)<<"institute.country() : "                       <<institute.country()<<std::endl;
    std::cout<<std::setw(width)<<"institute.link() : "                          <<institute.link()<<std::endl;
    std::cout<<std::setw(width)<<"institute.contact() : "                       <<institute.contact()<<std::endl;

    // Setters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Setters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"institute.title(\"Title\") : "                <<institute.title("Title")<<std::endl;
    std::cout<<std::setw(width)<<"institute.name(\"Name\") : "                  <<institute.name("Name")<<std::endl;
    std::cout<<std::setw(width)<<"institute.extended(\"Extended\") : "          <<institute.extended("Extended")<<std::endl;
    std::cout<<std::setw(width)<<"institute.street(\"Street\") : "              <<institute.street("Street")<<std::endl;
    std::cout<<std::setw(width)<<"institute.zip(\"Zip\") : "                    <<institute.zip("Zip")<<std::endl;
    std::cout<<std::setw(width)<<"institute.city(\"City\") : "                  <<institute.city("City")<<std::endl;
    std::cout<<std::setw(width)<<"institute.region(\"Region\") : "              <<institute.region("Region")<<std::endl;
    std::cout<<std::setw(width)<<"institute.country(\"Country\") : "            <<institute.country("Country")<<std::endl;
    std::cout<<std::setw(width)<<"institute.link(\"Link\") : "                  <<institute.link("Link")<<std::endl; 
    std::cout<<std::setw(width)<<"institute.contact(\"Contact\") : "            <<institute.contact("Contact")<<std::endl; 
    
    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Predefined : "                                <<std::endl;
    std::cout<<std::setw(width)<<"institute.luth() : "                          <<institute.luth()<<std::endl;
    std::cout<<std::setw(width)<<"institute.sap() : "                           <<institute.sap()<<std::endl;
    std::cout<<std::setw(width)<<"institute.ipag() : "                          <<institute.ipag()<<std::endl;
    std::cout<<std::setw(width)<<"institute.obspm() : "                         <<institute.obspm()<<std::endl;
    std::cout<<std::setw(width)<<"institute.iap() : "                           <<institute.iap()<<std::endl;
    std::cout<<std::setw(width)<<"institute.ias() : "                           <<institute.ias()<<std::endl;
    std::cout<<std::setw(width)<<"institute.cnrs() : "                          <<institute.cnrs()<<std::endl;
    std::cout<<std::setw(width)<<"institute.cnes() : "                          <<institute.cnes()<<std::endl;
    std::cout<<std::setw(width)<<"institute.cea() : "                           <<institute.cea()<<std::endl;
    std::cout<<std::setw(width)<<"institute.esa() : "                           <<institute.esa()<<std::endl;
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = AboutInstitute::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABOUTINSTITUTE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
