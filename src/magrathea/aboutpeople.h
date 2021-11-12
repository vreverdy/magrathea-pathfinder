/* ******************************* ABOUTPEOPLE ****************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          AboutPeople
// DESCRIPTION :    Information about a developer, an author, or a contributor
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           aboutpeople.h
/// \brief          Information about a developer, an author, or a contributor
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef ABOUTPEOPLE_H_INCLUDED
#define ABOUTPEOPLE_H_INCLUDED
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
// Information about a developer, an author, or a contributor
/// \brief          Information about a developer, an author, or a contributor.
/// \details        Name, status, contact, link, dates of contributions... for 
///                 authors and contributors.
class AboutPeople final
: public AbstractAboutObject<AboutPeople, std::string, std::string, int, int, std::string, std::string, std::string, std::string>
{
    // Setup
    public: using AbstractAboutObject<AboutPeople, std::string, std::string, int, int, std::string, std::string, std::string, std::string>::operator=;
    
    // Lifecycle
    /// \name           Lifecycle
    //@{
    public:
        template <class... Misc> explicit inline AboutPeople(Misc&&... misc);
    //@}
    
    // Accessors
    /// \name           Accessors 
    //@{
    public:
        inline std::string name(const std::string& separator = " ") const;
        inline std::string years(const std::string& separator = "-") const;
        inline std::string mails(const std::string& separator = " - ") const;
    //@}

    // Mutators
    /// \name           Mutators 
    //@{
    public:
        inline AboutPeople& name(const std::string& fname, const std::string& lname);
        inline AboutPeople& years(const int fyear, const int lyear);
        inline AboutPeople& mails(const std::string& fmail, const std::string& lmail);
    //@}
    
    // Getters
    /// \name           Getters 
    //@{
    public:
        inline const std::string& first() const;
        inline const std::string& last() const;
        inline int begin() const;
        inline int end() const;
        inline const std::string& mail() const;
        inline const std::string& altmail() const;
        inline const std::string& link() const;
        inline const std::string& contact() const;
    //@}

    // Setters
    /// \name           Setters 
    //@{
    public:
        inline AboutPeople& first(const std::string& value);
        inline AboutPeople& last(const std::string& value);
        inline AboutPeople& begin(const int value);
        inline AboutPeople& end(const int value);
        inline AboutPeople& mail(const std::string& value);
        inline AboutPeople& altmail(const std::string& value);
        inline AboutPeople& link(const std::string& value);
        inline AboutPeople& contact(const std::string& value);
    //@}

    // Predefined
    /// \name           Predefined 
    //@{
    public:
        static inline const AboutPeople& vreverdy();
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
inline AboutPeople::AboutPeople(Misc&&... misc)
: AbstractAboutObject<AboutPeople, std::string, std::string, int, int, std::string, std::string, std::string, std::string>(std::forward<Misc>(misc)...)
{
    ;
}
// -------------------------------------------------------------------------- //



// -------------------------------- ACCESSORS ------------------------------- //
// Full name accessor
/// \brief          Full name accessor.
/// \details        Gets the value of the first and last name properties. If one
///                 of the first or last name is empty, the separator is not 
///                 used.
/// \param[in]      separator String separator.
/// \return         Name.
inline std::string AboutPeople::name(const std::string& separator) const
{
    return std::get<0>(_data)+(((std::get<0>(_data).empty()) || (std::get<1>(_data).empty())) ? (std::string()) : (separator))+std::get<1>(_data);
}

// Year range accessor
/// \brief          Year range accessor.
/// \details        Gets the value of the begin and end year properties. If the
///                 end year is lower than the begin year, it is not displayed.
/// \param[in]      separator String separator.
/// \return         Years.
inline std::string AboutPeople::years(const std::string& separator) const
{
    return std::to_string(std::get<2>(_data))+((std::get<3>(_data) > std::get<2>(_data)) ? (separator+std::to_string(std::get<3>(_data))) : (std::string()));
}

// Mails accessor
/// \brief          Mails accessor.
/// \details        Gets the value of the main and alternative mails properties.
///                 If one of the mails is empty, the separator is not used.
/// \param[in]      separator String separator.
/// \return         Mails.
inline std::string AboutPeople::mails(const std::string& separator) const
{
    return std::get<4>(_data)+(((std::get<4>(_data).empty()) || (std::get<5>(_data).empty())) ? (std::string()) : (separator))+std::get<5>(_data);
}
// -------------------------------------------------------------------------- //



// -------------------------------- MUTATORS -------------------------------- //
// Full name mutator
/// \brief          Full name mutator.
/// \details        Sets the value of the first and last name properties.
/// \param[in]      fname First name.
/// \param[in]      lname Last name.
/// \return         Self reference.
inline AboutPeople& AboutPeople::name(const std::string& fname, const std::string& lname)
{
    std::get<0>(_data) = fname;
    std::get<1>(_data) = lname;
    return *this;
}

// Year range mutator
/// \brief          Year range mutator.
/// \details        Sets the value of the begin and end year properties.
/// \param[in]      fyear First year of contribution.
/// \param[in]      lyear Last year of contribution.
/// \return         Self reference.
inline AboutPeople& AboutPeople::years(const int fyear, const int lyear)
{
    std::get<2>(_data) = fyear;
    std::get<3>(_data) = lyear;
    return *this;
}

// Mails mutator
/// \brief          Mails mutator.
/// \details        Sets the value of the main and alternative mails properties.
/// \param[in]      fmail E-mail.
/// \param[in]      lmail Alternative e-mail.
/// \return         Mails.
inline AboutPeople& AboutPeople::mails(const std::string& fmail, const std::string& lmail)
{
    std::get<4>(_data) = fmail;
    std::get<5>(_data) = lmail;
    return *this;
}
// -------------------------------------------------------------------------- //



// --------------------------------- GETTERS -------------------------------- //
// First name getter
/// \brief          First name getter.
/// \details        Gets the value of the first property.
/// \return         First.
inline const std::string& AboutPeople::first() const
{
    return std::get<0>(_data);
}

// Last name getter
/// \brief          Last name getter.
/// \details        Gets the value of the last property.
/// \return         Last.
inline const std::string& AboutPeople::last() const
{
    return std::get<1>(_data);
}

// First year of contribution getter
/// \brief          First year of contribution getter.
/// \details        Gets the value of the begin property.
/// \return         Begin.
inline int AboutPeople::begin() const
{
    return std::get<2>(_data);
}

// Last year of contribution getter
/// \brief          Last year of contribution getter.
/// \details        Gets the value of the end property.
/// \return         End.
inline int AboutPeople::end() const
{
    return std::get<3>(_data);
}

// E-mail getter
/// \brief          E-mail getter.
/// \details        Gets the value of the mail property.
/// \return         Mail.
inline const std::string& AboutPeople::mail() const
{
    return std::get<4>(_data);
}

// Alternative e-mail getter
/// \brief          Alternative e-mail getter.
/// \details        Gets the value of the altmail property.
/// \return         Altmail.
inline const std::string& AboutPeople::altmail() const
{
    return std::get<5>(_data);
}

// Web link getter
/// \brief          Web link getter.
/// \details        Gets the value of the link property.
/// \return         Link.
inline const std::string& AboutPeople::link() const
{
    return std::get<6>(_data);
}

// Additional contact information getter
/// \brief          Additional contact information getter.
/// \details        Gets the value of the contact property.
/// \return         Contact.
inline const std::string& AboutPeople::contact() const
{
    return std::get<7>(_data);
}
// -------------------------------------------------------------------------- //



// --------------------------------- SETTERS -------------------------------- //
// First name setter
/// \brief          First name setter.
/// \details        Sets the value of the first property.
/// \param[in]      value First.
/// \return         Self reference.
inline AboutPeople& AboutPeople::first(const std::string& value)
{
    std::get<0>(_data) = value;
    return *this;
}

// Last name setter
/// \brief          Last name setter.
/// \details        Sets the value of the last property.
/// \param[in]      value Last.
/// \return         Self reference.
inline AboutPeople& AboutPeople::last(const std::string& value)
{
    std::get<1>(_data) = value;
    return *this;
}

// First year of contribution setter
/// \brief          First year of contribution setter.
/// \details        Sets the value of the begin property.
/// \param[in]      value Begin.
/// \return         Self reference.
inline AboutPeople& AboutPeople::begin(const int value)
{
    std::get<2>(_data) = value;
    return *this;
}

// Last year of contribution setter
/// \brief          Last year of contribution setter.
/// \details        Sets the value of the end property.
/// \param[in]      value End.
/// \return         Self reference.
inline AboutPeople& AboutPeople::end(const int value)
{
    std::get<3>(_data) = value;
    return *this;
}

// E-mail setter
/// \brief          E-mail setter.
/// \details        Sets the value of the mail property.
/// \param[in]      value Mail.
/// \return         Self reference.
inline AboutPeople& AboutPeople::mail(const std::string& value)
{
    std::get<4>(_data) = value;
    return *this;
}

// Alternative e-mail setter
/// \brief          Alternative e-mail setter.
/// \details        Sets the value of the altmail property.
/// \param[in]      value Altmail.
/// \return         Self reference.
inline AboutPeople& AboutPeople::altmail(const std::string& value)
{
    std::get<5>(_data) = value;
    return *this;
}

// Web link setter
/// \brief          Web link setter.
/// \details        Sets the value of the link property.
/// \param[in]      value Link.
/// \return         Self reference.
inline AboutPeople& AboutPeople::link(const std::string& value)
{
    std::get<6>(_data) = value;
    return *this;
}

// Additional contact information setter
/// \brief          Additional contact information setter.
/// \details        Sets the value of the contact property.
/// \param[in]      value Contact.
/// \return         Self reference.
inline AboutPeople& AboutPeople::contact(const std::string& value)
{
    std::get<7>(_data) = value;
    return *this;
}
// -------------------------------------------------------------------------- //



// ------------------------------- PREDEFINED ------------------------------- //
// Predefined Vincent Reverdy
/// \brief          Predefined Vincent Reverdy.
/// \details        Vincent Reverdy details.
/// \return         Immutable reference to people singleton.
inline const AboutPeople& AboutPeople::vreverdy()
{
    static const AboutPeople singleton("Vincent", "Reverdy", 2012, 2013, "vince.rev@gmail.com", "vincent.reverdy@obspm.fr", "http://vreverdy.fr", "");
    return singleton;
}
// -------------------------------------------------------------------------- //   



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of AboutPeople.
/// \return         0 if no error.
int AboutPeople::example()
{
    // Initialize
    std::cout<<"BEGIN = AboutPeople::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    AboutPeople people("first", "last", 0, 1, "mail", "altmail", "link", "contact");

    // Lifecycle
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"AboutPeople() : "                             <<AboutPeople()<<std::endl;
    std::cout<<std::setw(width)<<"AboutPeople(\"people\") : "                   <<AboutPeople("people")<<std::endl;
    std::cout<<std::setw(width)<<"AboutPeople(people) : "                       <<AboutPeople(people)<<std::endl; 

    // Accessors
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Accessors : "                                 <<std::endl;
    std::cout<<std::setw(width)<<"people.name() : "                             <<people.name()<<std::endl;
    std::cout<<std::setw(width)<<"people.years() : "                            <<people.years()<<std::endl;
    std::cout<<std::setw(width)<<"people.mails() : "                            <<people.mails()<<std::endl;

    // Mutators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Mutators : "                                  <<std::endl;
    std::cout<<std::setw(width)<<"people.name(\"f\", \"l\") : "                 <<people.name("f", "l")<<std::endl;
    std::cout<<std::setw(width)<<"people.years(2, 3) : "                        <<people.years(2, 3)<<std::endl;
    std::cout<<std::setw(width)<<"people.mails(\"m\", \"a\") : "                <<people.mails("m", "a")<<std::endl;
    
    // Getters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Getters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"people.first() : "                            <<people.first()<<std::endl;
    std::cout<<std::setw(width)<<"people.last() : "                             <<people.last()<<std::endl;
    std::cout<<std::setw(width)<<"people.begin() : "                            <<people.begin()<<std::endl;
    std::cout<<std::setw(width)<<"people.end() : "                              <<people.end()<<std::endl;
    std::cout<<std::setw(width)<<"people.mail() : "                             <<people.mail()<<std::endl;
    std::cout<<std::setw(width)<<"people.altmail() : "                          <<people.altmail()<<std::endl;
    std::cout<<std::setw(width)<<"people.link() : "                             <<people.link()<<std::endl;
    std::cout<<std::setw(width)<<"people.contact() : "                          <<people.contact()<<std::endl;

    // Setters
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Setters : "                                   <<std::endl;
    std::cout<<std::setw(width)<<"people.first(\"First\") : "                   <<people.first("First")<<std::endl;
    std::cout<<std::setw(width)<<"people.last(\"Last\") : "                     <<people.last("Last")<<std::endl;
    std::cout<<std::setw(width)<<"people.begin(4) : "                           <<people.begin(4)<<std::endl;
    std::cout<<std::setw(width)<<"people.end(5) : "                             <<people.end(5)<<std::endl;
    std::cout<<std::setw(width)<<"people.mail(\"Mail\") : "                     <<people.mail("Mail")<<std::endl;
    std::cout<<std::setw(width)<<"people.altmail(\"Altmail\") : "               <<people.altmail("Altmail")<<std::endl;
    std::cout<<std::setw(width)<<"people.link(\"Link\") : "                     <<people.link("Link")<<std::endl;
    std::cout<<std::setw(width)<<"people.contact(\"Contact\") : "               <<people.contact("Contact")<<std::endl;

    // Predefined
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Predefined : "                                <<std::endl;
    std::cout<<std::setw(width)<<"people.vreverdy() : "                         <<people.vreverdy()<<std::endl;   
    
    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = AboutPeople::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // ABOUTPEOPLE_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/


