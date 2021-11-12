/* ******************************** CONSTANTS ******************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        MAGRATHEA
// TITLE :          Constants
// DESCRIPTION :    Common mathematical and physical constants
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           constants.h
/// \brief          Common mathematical and physical constants
/// \author         Vincent Reverdy (vince.rev@gmail.com)
/// \date           2012-2013
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/



// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <iostream>
#include <iomanip>
// Include libs
// Include project
#include "constant.h"
// Misc
namespace magrathea { 
// -------------------------------------------------------------------------- //



// ---------------------------------- CLASS --------------------------------- //
// Common mathematical and physical constants
/// \brief          Common mathematical and physical constants.
/// \details        Provides an unified access to common mathematical and 
///                 physical constants in SI units. When the square of the 
///                 constant is widely used in common formulas, it is provided
///                 too. The values were taken in 2012 from : <ul>
///                 <li>Wikipedia : 
///                 http://en.wikipedia.org/wiki/Mathematical_constant
///                 http://en.wikipedia.org/wiki/Physical_constants</li>
///                 <li>CODATA : http://www.codata.org</li>
///                 <li>Particle Data Group : http://pdg.lbl.gov</li>
///                 <li>International Astronomical Union : http://www.iau.org
///                 </li></ul>
/// \tparam         Type Type of the constants.
template <typename Type = double>
class Constants final
{
    // Mathematical constants
    /// \name           Mathematical constants
    //@{
    public:
        static constexpr const Constant<Type> pi();
        static constexpr const Constant<Type> napier();
        static constexpr const Constant<Type> euler();
        static constexpr const Constant<Type> golden();
        static constexpr const Constant<Type> pythagoras();
        static constexpr const Constant<Type> theodorus();
        static constexpr const Constant<Type> glaisher();
        static constexpr const Constant<Type> khinchin();
        static constexpr const Constant<Type> soldner();
        static constexpr const Constant<Type> mertens();
        static constexpr const Constant<Type> plastic();
        static constexpr const Constant<Type> catalan();
        static constexpr const Constant<Type> feigenbaumd();
        static constexpr const Constant<Type> feigenbauma();
        static constexpr const Constant<Type> omega();
        static constexpr const Constant<Type> tanarcsec();
        static constexpr const Constant<Type> sqrtpi();
        static constexpr const Constant<Type> pi2();
    //@}

    // Universal constants
    /// \name           Universal constants
    //@{
    public:
        static constexpr const Constant<Type> c();
        static constexpr const Constant<Type> g();
        static constexpr const Constant<Type> h();
        static constexpr const Constant<Type> hbar();
        static constexpr const Constant<Type> c2();
        static constexpr const Constant<Type> g2();
        static constexpr const Constant<Type> h2();
        static constexpr const Constant<Type> hbar2();
        static constexpr const Constant<Type> c4();
    //@}

    // Units
    /// \name           Units
    //@{
    public:
        static constexpr const Constant<Type> celsius();
        static constexpr const Constant<Type> atomic();
        static constexpr const Constant<Type> atm();
        static constexpr const Constant<Type> gn();
        static constexpr const Constant<Type> au();
        static constexpr const Constant<Type> ev();
        static constexpr const Constant<Type> evc2();
        static constexpr const Constant<Type> deg();
        static constexpr const Constant<Type> arcmin();
        static constexpr const Constant<Type> arcsec();
        static constexpr const Constant<Type> deg2();
        static constexpr const Constant<Type> arcmin2();
        static constexpr const Constant<Type> arcsec2();
        static constexpr const Constant<Type> sphere();
        static constexpr const Constant<Type> second();
        static constexpr const Constant<Type> minute();
        static constexpr const Constant<Type> hour();
        static constexpr const Constant<Type> day();
        static constexpr const Constant<Type> year();
        static constexpr const Constant<Type> ly();
        static constexpr const Constant<Type> pc();
    //@}

    // Particle masses
    /// \name           Particle masses
    //@{
    public:
        static constexpr const Constant<Type> mquarku();
        static constexpr const Constant<Type> mquarkd();
        static constexpr const Constant<Type> mquarkc();
        static constexpr const Constant<Type> mquarks();
        static constexpr const Constant<Type> mquarkt();
        static constexpr const Constant<Type> mquarkb();
        static constexpr const Constant<Type> melectron();
        static constexpr const Constant<Type> mmuon();
        static constexpr const Constant<Type> mtau();
        static constexpr const Constant<Type> mnuelectron();
        static constexpr const Constant<Type> mnumuon();
        static constexpr const Constant<Type> mnutau();
        static constexpr const Constant<Type> mphoton();
        static constexpr const Constant<Type> mbosonw();
        static constexpr const Constant<Type> mbosonz();
        static constexpr const Constant<Type> mgluon();
        static constexpr const Constant<Type> mhiggs();
        static constexpr const Constant<Type> mgraviton();
        static constexpr const Constant<Type> mproton();
        static constexpr const Constant<Type> mneutron();
    //@}

    // Electromagnetic constants
    /// \name           Electromagnetic constants
    //@{
    public:
        static constexpr const Constant<Type> fcd();
        static constexpr const Constant<Type> kcd();
        static constexpr const Constant<Type> e();
        static constexpr const Constant<Type> e2();
        static constexpr const Constant<Type> mu0();
        static constexpr const Constant<Type> epsilon0();
        static constexpr const Constant<Type> z0();
        static constexpr const Constant<Type> ke();
        static constexpr const Constant<Type> magnetonb();
        static constexpr const Constant<Type> magnetonn();
        static constexpr const Constant<Type> quantumc();
        static constexpr const Constant<Type> quantumf();
        static constexpr const Constant<Type> josephson();
        static constexpr const Constant<Type> klitzing();
    //@}

    // Physico-chemical constants
    /// \name           Physico-chemical constants
    //@{
    public:
        static constexpr const Constant<Type> na();
        static constexpr const Constant<Type> wien();
        static constexpr const Constant<Type> kb();
        static constexpr const Constant<Type> kb2();
        static constexpr const Constant<Type> gas();
        static constexpr const Constant<Type> radiationf();
        static constexpr const Constant<Type> radiations();
        static constexpr const Constant<Type> faraday();
        static constexpr const Constant<Type> stefan();
    //@}

    // Atomic and nuclear constants
    /// \name           Atomic and nuclear constants
    //@{
    public:
        static constexpr const Constant<Type> hfs();
        static constexpr const Constant<Type> alpha();
        static constexpr const Constant<Type> alpha2();
        static constexpr const Constant<Type> rydberg();
        static constexpr const Constant<Type> relectron();
        static constexpr const Constant<Type> rbohr();
        static constexpr const Constant<Type> hartree();
        static constexpr const Constant<Type> thomson();
    //@}

    // Planck units
    /// \name           Planck units
    //@{
    public:
        static constexpr const Constant<Type> planckq();
        static constexpr const Constant<Type> planckl();
        static constexpr const Constant<Type> planckm();
        static constexpr const Constant<Type> planckt();
        static constexpr const Constant<Type> planckf();
        static constexpr const Constant<Type> plancke();
        static constexpr const Constant<Type> planckp();
        static constexpr const Constant<Type> plancktheta();
    //@}

    // Solar system masses
    /// \name           Solar system masses
    //@{
    public:
        static constexpr const Constant<Type> msun();
        static constexpr const Constant<Type> mearth();
        static constexpr const Constant<Type> mmercury();
        static constexpr const Constant<Type> mvenus();
        static constexpr const Constant<Type> mmars();
        static constexpr const Constant<Type> mjupiter();
        static constexpr const Constant<Type> msaturn();
        static constexpr const Constant<Type> muranus();
        static constexpr const Constant<Type> mneptune();
        static constexpr const Constant<Type> mmoon();
    //@}

    // Solar system radii
    /// \name           Solar system radii
    //@{
    public:
        static constexpr const Constant<Type> rsun();
        static constexpr const Constant<Type> rearth();
        static constexpr const Constant<Type> rmercury();
        static constexpr const Constant<Type> rvenus();
        static constexpr const Constant<Type> rmars();
        static constexpr const Constant<Type> rjupiter();
        static constexpr const Constant<Type> rsaturn();
        static constexpr const Constant<Type> ruranus();
        static constexpr const Constant<Type> rneptune();
        static constexpr const Constant<Type> rmoon();
    //@}

    // Test
    /// \name           Test
    //@{
    public:
        static int example();
    //@}
};
// -------------------------------------------------------------------------- //



// ------------------------- MATHEMATICAL CONSTANTS ------------------------- //
// Pi
/// \brief          Pi.
/// \details        Gets the value of \f$\pi\f$.
/// \return         \f$\pi\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::pi()
{
    return 3.1415926535897932384626433832795028841971693993751058209749445L;
}

// Napier constant
/// \brief          Napier constant.
/// \details        Gets the value of \f$e\f$.
/// \return         \f$e\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::napier()
{
    return 2.7182818284590452353602874713526624977572470936999595749669676L;
}

// Euler-Mascheroni constant
/// \brief          Euler-Mascheroni constant.
/// \details        Gets the value of \f$\gamma\f$.
/// \return         \f$\gamma\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::euler()
{
    return 0.5772156649015328606065120900824024310421593359399235988057672L;
}

// Golden ratio
/// \brief          Golden ratio.
/// \details        Gets the value of \f$\varphi = \frac{1+\sqrt{5}}{2}\f$.
/// \return         \f$\varphi = \frac{1+\sqrt{5}}{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::golden()
{
    return 1.6180339887498948482045868343656381177203091798057628621354486L;
}

// Pythagoras constant, square root of 2
/// \brief          Pythagoras constant, square root of 2.
/// \details        Gets the value of \f$\sqrt{2}\f$.
/// \return         \f$\sqrt{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::pythagoras()
{
    return 1.4142135623730950488016887242096980785696718753769480731766797L;
}

// Theodorus constant, square root of 3
/// \brief          Theodorus constant, square root of 3.
/// \details        Gets the value of \f$\sqrt{3}\f$.
/// \return         \f$\sqrt{3}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::theodorus()
{
    return 1.7320508075688772935274463415058723669428052538103806280558069L;
}

// Glaisher-Kinkelin constant
/// \brief          Glaisher-Kinkelin constant.
/// \details        Gets the value of \f$A\f$.
/// \return         \f$A\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::glaisher()
{
    return 1.2824271291006226368753425688697917277676889273250011920637400L;
}

// Khinchin constant
/// \brief          Khinchin constant.
/// \details        Gets the value of \f$K_{0}\f$.
/// \return         \f$K_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::khinchin()
{
    return 2.6854520010653064453097148354817956938203822939944629530511523L;
}

// Ramanujan-Soldner constant
/// \brief          Ramanujan-Soldner constant.
/// \details        Gets the value of \f$\mu\f$.
/// \return         \f$\mu\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::soldner()
{
    return 1.4513692348833810502839684858920274494930322836480158630930045L;
}

// Meissel-Mertens constant
/// \brief          Meissel-Mertens constant.
/// \details        Gets the value of \f$M_{1}\f$.
/// \return         \f$M_{1}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mertens()
{
    return 0.2614972128476427837554268386086958590515666482611992061920642L;
}

// Plastic number
/// \brief          Plastic number.
/// \details        Gets the value of \f$\rho\f$.
/// \return         \f$\rho\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::plastic()
{
    return 1.3247179572447460259609088544780973407344040569017333645340150L;
}

// Catalan constant
/// \brief          Catalan constant.
/// \details        Gets the value of \f$K\f$.
/// \return         \f$K\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::catalan()
{
    return 0.9159655941772190150546035149323841107741493742816721342664981L;
}

// First Feigenbaum constant
/// \brief          First Feigenbaum constant.
/// \details        Gets the value of \f$\delta\f$.
/// \return         \f$\delta\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::feigenbaumd()
{
    return 4.6692016091029906718532038204662016172581855774757686327456513L;
}

// Second Feigenbaum constant
/// \brief          Second Feigenbaum constant.
/// \details        Gets the value of \f$\alpha\f$.
/// \return         \f$\alpha\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::feigenbauma()
{
    return 2.5029078750958928222839028732182157863812713767271499773361920L;
}

// Omega constant
/// \brief          Omega constant.
/// \details        Gets the value of \f$\Omega\f$.
/// \return         \f$\Omega\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::omega()
{
    return 0.5671432904097838729999686622103555497538157871865125081351310L;
}

// Value of the tangent of one arcsecond
/// \brief          Value of the tangent of one arcsecond.
/// \details        Gets the value of \f$\tan\left(arcsec\right)\f$.
/// \return         \f$\tan\left(arcsec\right)\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::tanarcsec()
{
    return 4.8481368111333441675396429478852851658848753880815033842788396E-6L;
}

// Square root of pi
/// \brief          Square root of pi.
/// \details        Gets the value of \f$\sqrt{\pi}\f$.
/// \return         \f$\sqrt{\pi}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::sqrtpi()
{
    return 1.7724538509055160272981674833411451827975494561223871282138077L;
}

// Squared pi
/// \brief          Squared pi.
/// \details        Gets the value of \f$\pi^{2}\f$.
/// \return         \f$\pi^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::pi2()
{
    return pi()*pi();
}
// -------------------------------------------------------------------------- //



// --------------------------- UNIVERSAL CONSTANTS -------------------------- //
// Speed of light in vacuum
/// \brief          Speed of light in vacuum.
/// \details        Gets the value of \f$c\f$.
/// \return         \f$c\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::c()
{
    return 299792458;
}

// Newtonian constant of gravitation
/// \brief          Newtonian constant of gravitation.
/// \details        Gets the value of \f$G\f$.
/// \return         \f$G\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::g()
{
    return 6.67384E-11;
}

// Planck constant
/// \brief          Planck constant.
/// \details        Gets the value of \f$h\f$.
/// \return         \f$h\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::h()
{
    return 6.62606957E-34;
}

// Reduced Planck constant
/// \brief          Reduced Planck constant.
/// \details        Gets the value of \f$\bar{h}\f$.
/// \return         \f$\bar{h}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::hbar()
{
    return h()/(2*pi());
}

// Squared speed of light in vacuum
/// \brief          Squared speed of light in vacuum.
/// \details        Gets the value of \f$c^{2}\f$.
/// \return         \f$c^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::c2()
{
    return c()*c();
}

// Squared newtonian constant of gravitation
/// \brief          Squared newtonian constant of gravitation.
/// \details        Gets the value of \f$g^{2}\f$.
/// \return         \f$g^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::g2()
{
    return g()*g();
}

// Squared Planck constant
/// \brief          Squared Planck constant.
/// \details        Gets the value of \f$h^{2}\f$.
/// \return         \f$h^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::h2()
{
    return h()*h();
}

// Squared reduced Planck constant
/// \brief          Squared reduced Planck constant.
/// \details        Gets the value of \f$\bar{h}^{2}\f$.
/// \return         \f$\bar{h}^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::hbar2()
{
    return hbar()*hbar();
}

// Speed of light in vacuum to the fourth
/// \brief          Speed of light in vacuum to the fourth.
/// \details        Gets the value of \f$c^{4}\f$.
/// \return         \f$c^{4}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::c4()
{
    return c2()*c2();
}
// -------------------------------------------------------------------------- //



// ---------------------------------- UNITS --------------------------------- //
// Value of zero celsius degree in kelvin
/// \brief          Value of zero celsius degree in kelvin.
/// \details        Gets the value of \f$^{\circ}C\f$.
/// \return         \f$^{\circ}C\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::celsius()
{
    return 273.15;
}

// Value of the atomic mass in kilograms
/// \brief          Value of the atomic mass in kilograms.
/// \details        Gets the value of \f$u\f$.
/// \return         \f$u\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::atomic()
{
    return 1.660538921E-27;
}

// Value of one atmosphere in pascals
/// \brief          Value of one atmosphere in pascals.
/// \details        Gets the value of \f$atm\f$.
/// \return         \f$atm\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::atm()
{
    return 101325;
}

// Acceleration of gravity in meters per second squared
/// \brief          Acceleration of gravity in meters per second squared.
/// \details        Gets the value of \f$g_{n}\f$.
/// \return         \f$g_{n}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::gn()
{
    return 9.80665;
}

// Value of an astronomical unit in meters
/// \brief          Value of an astronomical unit in meters.
/// \details        Gets the value of \f$au\f$.
/// \return         \f$au\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::au()
{
    return 1.495978707E11;
}

// Value of an electron-volt in joules
/// \brief          Value of an electron-volt in joules.
/// \details        Gets the value of \f$eV\f$.
/// \return         \f$eV\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::ev()
{
    return 1.602176565E-19;
}

// Value of the electron-volt mass unit in kilograms
/// \brief          Value of the electron-volt mass unit in kilograms.
/// \details        Gets the value of \f$\frac{eV}{c^{2}}\f$.
/// \return         \f$\frac{eV}{c^{2}}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::evc2()
{
    return ev()/c2();
}

// Value of a degree in radians
/// \brief          Value of a degree in radians.
/// \details        Gets the value of \f$deg\f$.
/// \return         \f$deg\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::deg()
{
    return pi()/180;
}

// Value of an arcminute in radians
/// \brief          Value of an arcminute in radians.
/// \details        Gets the value of \f$arcmin\f$.
/// \return         \f$arcmin\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::arcmin()
{
    return deg()/60;
}

// Value of an arcsecond in radians
/// \brief          Value of an arcsecond in radians.
/// \details        Gets the value of \f$arcsec\f$.
/// \return         \f$arcsec\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::arcsec()
{
    return arcmin()/60;
}

// Value of a square degree in steradians
/// \brief          Value of a square degree in steradians.
/// \details        Gets the value of \f$deg^{2}\f$.
/// \return         \f$deg^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::deg2()
{
    return deg()*deg();
}

// Value of a square arcminute in steradians
/// \brief          Value of a square arcminute in steradians.
/// \details        Gets the value of \f$arcmin^{2}\f$.
/// \return         \f$arcmin^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::arcmin2()
{
    return arcmin()*arcmin();
}

// Value of a square arcsecond in steradians
/// \brief          Value of a square arcsecond in steradians.
/// \details        Gets the value of \f$arcsec^{2}\f$.
/// \return         \f$arcsec^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::arcsec2()
{
    return arcsec()*arcsec();
}

// Value of a sphere in steradians
/// \brief          Value of a sphere in steradians.
/// \details        Gets the value of \f$4 \pi sr\f$.
/// \return         \f$4 \pi sr\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::sphere()
{
    return 4*pi();
}

// Duration of a second in seconds
/// \brief          Duration of a second in seconds.
/// \details        Gets the value of \f$s\f$.
/// \return         \f$s\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::second()
{
    return 1;
}

// Duration of a minute in seconds
/// \brief          Duration of a minute in seconds.
/// \details        Gets the value of \f$min\f$.
/// \return         \f$min\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::minute()
{
    return 60*second();
}

// Duration of an hour in seconds
/// \brief          Duration of an hour in seconds.
/// \details        Gets the value of \f$h\f$.
/// \return         \f$h\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::hour()
{
    return 60*minute();
}

// Duration of a SI julian day in seconds
/// \brief          Duration of a SI julian day in seconds.
/// \details        Gets the value of \f$d\f$.
/// \return         \f$d\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::day()
{
    return 24*hour();
}

// Duration of a SI julian year in seconds
/// \brief          Duration of a SI julian year in seconds.
/// \details        Gets the value of \f$year\f$.
/// \return         \f$year\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::year()
{
    return 365*day()+6*hour();
}

// Value of a light-year in meters
/// \brief          Value of a light-year in meters.
/// \details        Gets the value of \f$ly\f$.
/// \return         \f$ly\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::ly()
{
    return c()*year();
}

// Value of a parsec in meters
/// \brief          Value of a parsec in meters.
/// \details        Gets the value of \f$pc\f$.
/// \return         \f$pc\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::pc()
{
    return au()/tanarcsec();
}
// -------------------------------------------------------------------------- //



// ----------------------------- PARTICLE MASSES ---------------------------- //
// Up quark mass
/// \brief          Up quark mass.
/// \details        Gets the value of \f$m_{u}\f$.
/// \return         \f$m_{u}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mquarku()
{
    return 2.3E6*evc2();
}

// Down quark mass
/// \brief          Down quark mass.
/// \details        Gets the value of \f$m_{d}\f$.
/// \return         \f$m_{d}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mquarkd()
{
    return 4.8E6*evc2();
}

// Charm quark mass
/// \brief          Charm quark mass.
/// \details        Gets the value of \f$m_{c}\f$.
/// \return         \f$m_{c}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mquarkc()
{
    return 1.275E9*evc2();
}

// Strange quark mass
/// \brief          Strange quark mass.
/// \details        Gets the value of \f$m_{s}\f$.
/// \return         \f$m_{s}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mquarks()
{
    return 95.E6*evc2();
}

// Top quark mass
/// \brief          Top quark mass.
/// \details        Gets the value of \f$m_{t}\f$.
/// \return         \f$m_{t}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mquarkt()
{
    return 173.500E9*evc2();
}

// Bottom quark mass
/// \brief          Bottom quark mass.
/// \details        Gets the value of \f$m_{b}\f$.
/// \return         \f$m_{b}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mquarkb()
{
    return 4.180E9*evc2();
}

// Electron mass
/// \brief          Electron mass.
/// \details        Gets the value of \f$m_{e}\f$.
/// \return         \f$m_{e}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::melectron()
{
    return 0.510998928E6*evc2();
}

// Muon mass
/// \brief          Muon mass.
/// \details        Gets the value of \f$m_{\mu}\f$.
/// \return         \f$m_{\mu}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mmuon()
{
    return 105.6583715E6*evc2();
}

// Tau mass
/// \brief          Tau mass.
/// \details        Gets the value of \f$m_{\tau}\f$.
/// \return         \f$m_{\tau}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mtau()
{
    return 1.77682E9*evc2();
}

// Electron neutrino mass
/// \brief          Electron neutrino mass.
/// \details        Gets the value of \f$m_{\nu_{e}}\f$.
/// \return         \f$m_{\nu_{e}}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mnuelectron()
{
    return 2.05*evc2();
}

// Muon neutrino mass
/// \brief          Muon neutrino mass.
/// \details        Gets the value of \f$m_{\nu_{\mu}}\f$.
/// \return         \f$m_{\nu_{\mu}}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mnumuon()
{
    return 0.19E6*evc2();
}

// Tau neutrino mass
/// \brief          Tau neutrino mass.
/// \details        Gets the value of \f$m_{\nu_{\tau}}\f$.
/// \return         \f$m_{\nu_{\tau}}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mnutau()
{
    return 18.2E6*evc2();
}

// Photon mass
/// \brief          Photon mass.
/// \details        Gets the value of \f$m_{\gamma}\f$.
/// \return         \f$m_{\gamma}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mphoton()
{
    return 0.*evc2();
}

// W boson mass
/// \brief          W boson mass.
/// \details        Gets the value of \f$m_{W}\f$.
/// \return         \f$m_{W}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mbosonw()
{
    return 80.385E9*evc2();
}

// Z boson mass
/// \brief          Z boson mass.
/// \details        Gets the value of \f$m_{Z}\f$.
/// \return         \f$m_{Z}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mbosonz()
{
    return 91.1876E9*evc2();
}

// Gluon mass
/// \brief          Gluon mass.
/// \details        Gets the value of \f$m_{g}\f$.
/// \return         \f$m_{g}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mgluon()
{
    return 0.*evc2();
}

// Higgs boson mass
/// \brief          Higgs boson mass.
/// \details        Gets the value of \f$m_{H}\f$.
/// \return         \f$m_{H}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mhiggs()
{
    return 125.300E9*evc2();
}

// Graviton mass
/// \brief          Graviton mass.
/// \details        Gets the value of \f$m_{graviton}\f$.
/// \return         \f$m_{graviton}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mgraviton()
{
    return 0.*evc2();
}

// Proton mass
/// \brief          Proton mass.
/// \details        Gets the value of \f$m_{p}\f$.
/// \return         \f$m_{p}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mproton()
{
    return 938.272046E6*evc2();
}

// Neutron mass
/// \brief          Neutron mass.
/// \details        Gets the value of \f$m_{n}\f$.
/// \return         \f$m_{n}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mneutron()
{
    return 939.565379E6*evc2();
}
// -------------------------------------------------------------------------- //



// ------------------------ ELECTROMAGNETIC CONSTANTS ----------------------- //
// Frequency of radiation used in candela definition
/// \brief          Frequency of radiation used in candela definition.
/// \details        Gets the value of \f$\nu_{cd}\f$.
/// \return         \f$\nu_{cd}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::fcd()
{
    return 5.4E14;
}

// Luminous efficacy used in candela definition
/// \brief          Luminous efficacy used in candela definition.
/// \details        Gets the value of \f$K_{cd}\f$.
/// \return         \f$K_{cd}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::kcd()
{
    return 683;
}

// Elementary charge
/// \brief          Elementary charge.
/// \details        Gets the value of \f$e\f$.
/// \return         \f$e\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::e()
{
    return ev();
}

// Squared elementary charge
/// \brief          Squared elementary charge.
/// \details        Gets the value of \f$e^{2}\f$.
/// \return         \f$e^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::e2()
{
    return e()*e();
}

// Vacuum permeability
/// \brief          Vacuum permeability.
/// \details        Gets the value of \f$\mu_{0}\f$.
/// \return         \f$\mu_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mu0()
{
    return 4.E-7*pi();
}

// Vacuum permittivity
/// \brief          Vacuum permittivity.
/// \details        Gets the value of \f$\epsilon_{0}\f$.
/// \return         \f$\epsilon_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::epsilon0()
{
    return 1/(mu0()*c2());
}

// Characteristic impedance of vacuum
/// \brief          Characteristic impedance of vacuum.
/// \details        Gets the value of \f$Z_{0}\f$.
/// \return         \f$Z_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::z0()
{
    return mu0()*c();
}

// Coulomb's constant
/// \brief          Coulomb's constant.
/// \details        Gets the value of \f$k_{e}\f$.
/// \return         \f$k_{e}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::ke()
{
    return 1/(4*pi()*epsilon0());
}

// Bohr magneton
/// \brief          Bohr magneton.
/// \details        Gets the value of \f$\mu_{B}\f$.
/// \return         \f$\mu_{B}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::magnetonb()
{
    return (e()*hbar())/(2*melectron());
}

// Nuclear magneton
/// \brief          Nuclear magneton.
/// \details        Gets the value of \f$\mu_{N}\f$.
/// \return         \f$\mu_{N}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::magnetonn()
{
    return (e()*hbar())/(2*mproton());
}

// Quantum conductance
/// \brief          Quantum conductance.
/// \details        Gets the value of \f$G_{0}\f$.
/// \return         \f$G_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::quantumc()
{
    return (2*e2())/h();
}

// Quantum flux
/// \brief          Quantum flux.
/// \details        Gets the value of \f$\phi_{0}\f$.
/// \return         \f$\phi_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::quantumf()
{
    return h()/(2*e());
}

// Josephson constant
/// \brief          Josephson constant.
/// \details        Gets the value of \f$K_{J}\f$.
/// \return         \f$K_{J}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::josephson()
{
    return (2*e())/h();
}

// Von Klitzing constant
/// \brief          Von Klitzing constant.
/// \details        Gets the value of \f$R_{K}\f$.
/// \return         \f$R_{K}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::klitzing()
{
    return h()/e2();
}
// -------------------------------------------------------------------------- //



// ----------------------- PHYSICO-CHEMICAL CONSTANTS ----------------------- //
// Avogadro number
/// \brief          Avogadro number.
/// \details        Gets the value of \f$N_{A}\f$.
/// \return         \f$N_{A}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::na()
{
    return 6.02214129E23;
}

// Wien constant
/// \brief          Wien constant.
/// \details        Gets the value of \f$b\f$.
/// \return         \f$b\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::wien()
{
    return 2.8977721E-3;
}

// Boltzmann constant
/// \brief          Boltzmann constant.
/// \details        Gets the value of \f$k_{B}\f$.
/// \return         \f$k_{B}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::kb()
{
    return 1.3806488E-23;
}

// Squared boltzmann constant
/// \brief          Squared boltzmann constant.
/// \details        Gets the value of \f$k_{B}^{2}\f$.
/// \return         \f$k_{B}^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::kb2()
{
    return kb()*kb();
}

// Ideal gas constant
/// \brief          Ideal gas constant.
/// \details        Gets the value of \f$R\f$.
/// \return         \f$R\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::gas()
{
    return na()*kb();
}

// First radiation constant
/// \brief          First radiation constant.
/// \details        Gets the value of \f$c_{1}\f$.
/// \return         \f$c_{1}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::radiationf()
{
    return 2*pi()*h()*c2();
}

// Second radiation constant
/// \brief          Second radiation constant.
/// \details        Gets the value of \f$c_{2}\f$.
/// \return         \f$c_{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::radiations()
{
    return (h()*c())/kb();
}

// Faraday constant
/// \brief          Faraday constant.
/// \details        Gets the value of \f$F\f$.
/// \return         \f$F\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::faraday()
{
    return na()*e();
}

// Stefan constant
/// \brief          Stefan constant.
/// \details        Gets the value of \f$\sigma\f$.
/// \return         \f$\sigma\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::stefan()
{
    return kb2()*kb2()*pi2()/(60*hbar2()*hbar()*c2());
}
// -------------------------------------------------------------------------- //



// ---------------------- ATOMIC AND NUCLEAR CONSTANTS ---------------------- //
// Ground state hyperfine splitting frequency of caesium 133
/// \brief          Ground state hyperfine splitting frequency of caesium 133.
/// \details        Gets the value of \f$\Delta\nu(^{133}Cs)_{hfs}\f$.
/// \return         \f$\Delta\nu(^{133}Cs)_{hfs}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::hfs()
{
    return 9.19263177E9;
}

// Fine structure constant
/// \brief          Fine structure constant.
/// \details        Gets the value of \f$\alpha\f$.
/// \return         \f$\alpha\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::alpha()
{
    return (mu0()*e2()*c())/(2*h());
}

// Squared fine structure constant
/// \brief          Squared fine structure constant.
/// \details        Gets the value of \f$\alpha^{2}\f$.
/// \return         \f$\alpha^{2}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::alpha2()
{
    return alpha()*alpha();
}

// Rydberg constant
/// \brief          Rydberg constant.
/// \details        Gets the value of \f$R_{\infty}\f$.
/// \return         \f$R_{\infty}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rydberg()
{
    return (alpha2()*melectron()*c())/(2*h());
}

// Classical electron radius
/// \brief          Classical electron radius.
/// \details        Gets the value of \f$r_{e}\f$.
/// \return         \f$r_{e}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::relectron()
{
    return (ke()*e2())/(melectron()*c2());
}

// Bohr radius
/// \brief          Bohr radius.
/// \details        Gets the value of \f$a_{0}\f$.
/// \return         \f$a_{0}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rbohr()
{
    return alpha()/(4*pi()*rydberg());
}

// Hartree energy
/// \brief          Hartree energy.
/// \details        Gets the value of \f$E_{h}\f$.
/// \return         \f$E_{h}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::hartree()
{
    return 2*rydberg()*h()*c();
}

// Thomson cross section
/// \brief          Thomson cross section.
/// \details        Gets the value of \f$\sigma_{T}\f$.
/// \return         \f$\sigma_{T}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::thomson()
{
    return (8*pi()*relectron()*relectron())/3;
}
// -------------------------------------------------------------------------- //



// ------------------------------ PLANCK UNITS ------------------------------ //
// Planck charge
/// \brief          Planck charge.
/// \details        Gets the value of \f$q_{P}\f$.
/// \return         \f$q_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::planckq()
{
    return 1.875545956;
}

// Planck length
/// \brief          Planck length.
/// \details        Gets the value of \f$l_{P}\f$.
/// \return         \f$l_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::planckl()
{
    return 1.616199;
}

// Planck mass
/// \brief          Planck mass.
/// \details        Gets the value of \f$m_{P}\f$.
/// \return         \f$m_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::planckm()
{
    return 2.17651;
}

// Planck time
/// \brief          Planck time.
/// \details        Gets the value of \f$t_{P}\f$.
/// \return         \f$t_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::planckt()
{
    return planckl()/c();
}

// Planck force
/// \brief          Planck force.
/// \details        Gets the value of \f$F_{P}\f$.
/// \return         \f$F_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::planckf()
{
    return (c2()*c2())/g();
}

// Planck energy
/// \brief          Planck energy.
/// \details        Gets the value of \f$E_{P}\f$.
/// \return         \f$E_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::plancke()
{
    return planckf()*planckl();
}

// Planck power
/// \brief          Planck power.
/// \details        Gets the value of \f$P_{P}\f$.
/// \return         \f$P_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::planckp()
{
    return plancke()/planckt();
}

// Planck temperature
/// \brief          Planck temperature.
/// \details        Gets the value of \f$\theta_{P}\f$.
/// \return         \f$\theta_{P}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::plancktheta()
{
    return (planckm()*c2())/kb();
}
// -------------------------------------------------------------------------- //



// --------------------------- SOLAR SYSTEM MASSES -------------------------- //
// Solar mass
/// \brief          Solar mass.
/// \details        Gets the value of \f$m_{\odot}\f$.
/// \return         \f$m_{\odot}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::msun()
{
    return 1.32712442099E20/g();
}

// Earth mass
/// \brief          Earth mass.
/// \details        Gets the value of \f$m_{E}\f$.
/// \return         \f$m_{E}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mearth()
{
    return 3.986004418E14/g();
}

// Mercury mass
/// \brief          Mercury mass.
/// \details        Gets the value of \f$m_{Me}\f$.
/// \return         \f$m_{Me}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mmercury()
{
    return msun()/6.0236E6;
}

// Venus mass
/// \brief          Venus mass.
/// \details        Gets the value of \f$m_{Ve}\f$.
/// \return         \f$m_{Ve}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mvenus()
{
    return msun()/4.08523719E5;
}

// Mars mass
/// \brief          Mars mass.
/// \details        Gets the value of \f$m_{Ma}\f$.
/// \return         \f$m_{Ma}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mmars()
{
    return msun()/3.09870359E6;
}

// Jupiter mass
/// \brief          Jupiter mass.
/// \details        Gets the value of \f$m_{J}\f$.
/// \return         \f$m_{J}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mjupiter()
{
    return msun()/1.047348644E3;
}

// Saturn mass
/// \brief          Saturn mass.
/// \details        Gets the value of \f$m_{Sa}\f$.
/// \return         \f$m_{Sa}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::msaturn()
{
    return msun()/3.4979018E3;
}

// Uranus mass
/// \brief          Uranus mass.
/// \details        Gets the value of \f$m_{U}\f$.
/// \return         \f$m_{U}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::muranus()
{
    return msun()/2.290298E4;
}

// Neptune mass
/// \brief          Neptune mass.
/// \details        Gets the value of \f$m_{N}\f$.
/// \return         \f$m_{N}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mneptune()
{
    return msun()/1.941226E4;
}

// Moon mass
/// \brief          Moon mass.
/// \details        Gets the value of \f$m_{Moon}\f$.
/// \return         \f$m_{Moon}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::mmoon()
{
    return 1.23000371E-2*mearth();
}
// -------------------------------------------------------------------------- //



// --------------------------- SOLAR SYSTEM RADII --------------------------- //
// Solar radius
/// \brief          Solar radius.
/// \details        Gets the value of \f$r_{\odot}\f$.
/// \return         \f$r_{\odot}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rsun()
{
    return 695990.0E3;
}

// Earth radius
/// \brief          Earth radius.
/// \details        Gets the value of \f$r_{E}\f$.
/// \return         \f$r_{E}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rmercury()
{
    return 2439.7E3;
}

// Mercury radius
/// \brief          Mercury radius.
/// \details        Gets the value of \f$r_{Me}\f$.
/// \return         \f$r_{Me}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rvenus()
{
    return 6051.8E3;
}

// Venus radius
/// \brief          Venus radius.
/// \details        Gets the value of \f$r_{Ve}\f$.
/// \return         \f$r_{Ve}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rearth()
{
    return 6371.00E3;
}

// Mars radius
/// \brief          Mars radius.
/// \details        Gets the value of \f$r_{Ma}\f$.
/// \return         \f$r_{Ma}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rmars()
{
    return 3389.50E3;
}

// Jupiter radius
/// \brief          Jupiter radius.
/// \details        Gets the value of \f$r_{J}\f$.
/// \return         \f$r_{J}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rjupiter()
{
    return 69911.0E3;
}

// Saturn radius
/// \brief          Saturn radius.
/// \details        Gets the value of \f$r_{Sa}\f$.
/// \return         \f$r_{Sa}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rsaturn()
{
    return 58232.0E3;
}

// Uranus radius
/// \brief          Uranus radius.
/// \details        Gets the value of \f$r_{U}\f$.
/// \return         \f$r_{U}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::ruranus()
{
    return 25362.0E3;
}

// Neptune radius
/// \brief          Neptune radius.
/// \details        Gets the value of \f$r_{N}\f$.
/// \return         \f$r_{N}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rneptune()
{
    return 24622.0E3;
}

// Moon radius
/// \brief          Moon radius.
/// \details        Gets the value of \f$r_{Moon}\f$.
/// \return         \f$r_{Moon}\f$.
template <typename Type>
constexpr const Constant<Type> Constants<Type>::rmoon()
{
    return 1737.5E3;
}
// -------------------------------------------------------------------------- //



// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Constants.
/// \return         0 if no error.
template <typename Type> 
int Constants<Type>::example()
{
    // Initialize
    std::cout<<"BEGIN = Constants::example()"<<std::endl;
    std::cout<<std::boolalpha<<std::left;
    static const unsigned int width = 40;
    
    // Construction
    Constants<double> constants;

    // Lifecycle and operators
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Lifecycle and operators : "                   <<std::endl;
    std::cout<<std::setw(width)<<"Constants<>() : "                             ; Constants<>(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Constants<double>() : "                       ; Constants<double>(); std::cout<<std::endl;
    std::cout<<std::setw(width)<<"constants = Constants<double>() : "           ; constants = Constants<double>(); std::cout<<std::endl;

    // Mathematical constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Mathematical constants : "                    <<std::endl;
    std::cout<<std::setw(width)<<"pi() : "                                      <<constants.pi()<<std::endl;
    std::cout<<std::setw(width)<<"napier() : "                                  <<constants.napier()<<std::endl;
    std::cout<<std::setw(width)<<"euler() : "                                   <<constants.euler()<<std::endl;
    std::cout<<std::setw(width)<<"golden() : "                                  <<constants.golden()<<std::endl;
    std::cout<<std::setw(width)<<"pythagoras() : "                              <<constants.pythagoras()<<std::endl;
    std::cout<<std::setw(width)<<"theodorus() : "                               <<constants.theodorus()<<std::endl;
    std::cout<<std::setw(width)<<"glaisher() : "                                <<constants.glaisher()<<std::endl;
    std::cout<<std::setw(width)<<"khinchin() : "                                <<constants.khinchin()<<std::endl;
    std::cout<<std::setw(width)<<"soldner() : "                                 <<constants.soldner()<<std::endl;
    std::cout<<std::setw(width)<<"mertens() : "                                 <<constants.mertens()<<std::endl;
    std::cout<<std::setw(width)<<"plastic() : "                                 <<constants.plastic()<<std::endl;
    std::cout<<std::setw(width)<<"catalan() : "                                 <<constants.catalan()<<std::endl;
    std::cout<<std::setw(width)<<"feigenbaumd() : "                             <<constants.feigenbaumd()<<std::endl;
    std::cout<<std::setw(width)<<"feigenbauma() : "                             <<constants.feigenbauma()<<std::endl;
    std::cout<<std::setw(width)<<"omega() : "                                   <<constants.omega()<<std::endl;
    std::cout<<std::setw(width)<<"tanarcsec() : "                               <<constants.tanarcsec()<<std::endl;
    std::cout<<std::setw(width)<<"sqrtpi() : "                                  <<constants.sqrtpi()<<std::endl;
    std::cout<<std::setw(width)<<"pi2() : "                                     <<constants.pi2()<<std::endl;

    // Universal constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Universal constants : "                       <<std::endl;
    std::cout<<std::setw(width)<<"c() : "                                       <<constants.c()<<std::endl;
    std::cout<<std::setw(width)<<"g() : "                                       <<constants.g()<<std::endl;
    std::cout<<std::setw(width)<<"h() : "                                       <<constants.h()<<std::endl;
    std::cout<<std::setw(width)<<"hbar() : "                                    <<constants.hbar()<<std::endl;
    std::cout<<std::setw(width)<<"c2() : "                                      <<constants.c2()<<std::endl;
    std::cout<<std::setw(width)<<"g2() : "                                      <<constants.g2()<<std::endl;
    std::cout<<std::setw(width)<<"h2() : "                                      <<constants.h2()<<std::endl;
    std::cout<<std::setw(width)<<"hbar2() : "                                   <<constants.hbar2()<<std::endl;
    std::cout<<std::setw(width)<<"c4() : "                                      <<constants.c4()<<std::endl;

    // Units
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Units : "                                     <<std::endl;
    std::cout<<std::setw(width)<<"celsius() : "                                 <<constants.celsius()<<std::endl;
    std::cout<<std::setw(width)<<"atomic() : "                                  <<constants.atomic()<<std::endl;
    std::cout<<std::setw(width)<<"atm() : "                                     <<constants.atm()<<std::endl;
    std::cout<<std::setw(width)<<"gn() : "                                      <<constants.gn()<<std::endl;
    std::cout<<std::setw(width)<<"au() : "                                      <<constants.au()<<std::endl;
    std::cout<<std::setw(width)<<"ev() : "                                      <<constants.ev()<<std::endl;
    std::cout<<std::setw(width)<<"evc2() : "                                    <<constants.evc2()<<std::endl;
    std::cout<<std::setw(width)<<"deg() : "                                     <<constants.deg()<<std::endl;
    std::cout<<std::setw(width)<<"arcmin() : "                                  <<constants.arcmin()<<std::endl;
    std::cout<<std::setw(width)<<"arcsec() : "                                  <<constants.arcsec()<<std::endl;
    std::cout<<std::setw(width)<<"deg2() : "                                    <<constants.deg2()<<std::endl;
    std::cout<<std::setw(width)<<"arcmin2() : "                                 <<constants.arcmin2()<<std::endl;
    std::cout<<std::setw(width)<<"arcsec2() : "                                 <<constants.arcsec2()<<std::endl;
    std::cout<<std::setw(width)<<"sphere() : "                                  <<constants.sphere()<<std::endl;
    std::cout<<std::setw(width)<<"second() : "                                  <<constants.second()<<std::endl;
    std::cout<<std::setw(width)<<"minute() : "                                  <<constants.minute()<<std::endl;
    std::cout<<std::setw(width)<<"hour() : "                                    <<constants.hour()<<std::endl;
    std::cout<<std::setw(width)<<"day() : "                                     <<constants.day()<<std::endl;
    std::cout<<std::setw(width)<<"year() : "                                    <<constants.year()<<std::endl;
    std::cout<<std::setw(width)<<"ly() : "                                      <<constants.ly()<<std::endl;
    std::cout<<std::setw(width)<<"pc() : "                                      <<constants.pc()<<std::endl;

    // Particle masses
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Particle masses : "                           <<std::endl;
    std::cout<<std::setw(width)<<"mquarku() : "                                 <<constants.mquarku()<<std::endl;
    std::cout<<std::setw(width)<<"mquarkd() : "                                 <<constants.mquarkd()<<std::endl;
    std::cout<<std::setw(width)<<"mquarkc() : "                                 <<constants.mquarkc()<<std::endl;
    std::cout<<std::setw(width)<<"mquarks() : "                                 <<constants.mquarks()<<std::endl;
    std::cout<<std::setw(width)<<"mquarkt() : "                                 <<constants.mquarkt()<<std::endl;
    std::cout<<std::setw(width)<<"mquarkb() : "                                 <<constants.mquarkb()<<std::endl;
    std::cout<<std::setw(width)<<"melectron() : "                               <<constants.melectron()<<std::endl;
    std::cout<<std::setw(width)<<"mmuon() : "                                   <<constants.mmuon()<<std::endl;
    std::cout<<std::setw(width)<<"mtau() : "                                    <<constants.mtau()<<std::endl;
    std::cout<<std::setw(width)<<"mnuelectron() : "                             <<constants.mnuelectron()<<std::endl;
    std::cout<<std::setw(width)<<"mnumuon() : "                                 <<constants.mnumuon()<<std::endl;
    std::cout<<std::setw(width)<<"mnutau() : "                                  <<constants.mnutau()<<std::endl;
    std::cout<<std::setw(width)<<"mphoton() : "                                 <<constants.mphoton()<<std::endl;
    std::cout<<std::setw(width)<<"mbosonw() : "                                 <<constants.mbosonw()<<std::endl;
    std::cout<<std::setw(width)<<"mbosonz() : "                                 <<constants.mbosonz()<<std::endl;
    std::cout<<std::setw(width)<<"mgluon() : "                                  <<constants.mgluon()<<std::endl;
    std::cout<<std::setw(width)<<"mhiggs() : "                                  <<constants.mhiggs()<<std::endl;
    std::cout<<std::setw(width)<<"mgraviton() : "                               <<constants.mgraviton()<<std::endl;
    std::cout<<std::setw(width)<<"mproton() : "                                 <<constants.mproton()<<std::endl;
    std::cout<<std::setw(width)<<"mneutron() : "                                <<constants.mneutron()<<std::endl;

    // Electromagnetic constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Electromagnetic constants : "                 <<std::endl;
    std::cout<<std::setw(width)<<"fcd() : "                                     <<constants.fcd()<<std::endl;
    std::cout<<std::setw(width)<<"kcd() : "                                     <<constants.kcd()<<std::endl;
    std::cout<<std::setw(width)<<"e() : "                                       <<constants.e()<<std::endl;
    std::cout<<std::setw(width)<<"e2() : "                                      <<constants.e2()<<std::endl;
    std::cout<<std::setw(width)<<"mu0() : "                                     <<constants.mu0()<<std::endl;
    std::cout<<std::setw(width)<<"epsilon0() : "                                <<constants.epsilon0()<<std::endl;
    std::cout<<std::setw(width)<<"z0() : "                                      <<constants.z0()<<std::endl;
    std::cout<<std::setw(width)<<"ke() : "                                      <<constants.ke()<<std::endl;
    std::cout<<std::setw(width)<<"magnetonb() : "                               <<constants.magnetonb()<<std::endl;
    std::cout<<std::setw(width)<<"magnetonn() : "                               <<constants.magnetonn()<<std::endl;
    std::cout<<std::setw(width)<<"quantumc() : "                                <<constants.quantumc()<<std::endl;
    std::cout<<std::setw(width)<<"quantumf() : "                                <<constants.quantumf()<<std::endl;
    std::cout<<std::setw(width)<<"josephson() : "                               <<constants.josephson()<<std::endl;
    std::cout<<std::setw(width)<<"klitzing() : "                                <<constants.klitzing()<<std::endl;

    // Physico-chemical constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Physico-chemical constants : "                <<std::endl;
    std::cout<<std::setw(width)<<"na() : "                                      <<constants.na()<<std::endl;
    std::cout<<std::setw(width)<<"wien() : "                                    <<constants.wien()<<std::endl;
    std::cout<<std::setw(width)<<"kb() : "                                      <<constants.kb()<<std::endl;
    std::cout<<std::setw(width)<<"kb2() : "                                     <<constants.kb2()<<std::endl;
    std::cout<<std::setw(width)<<"gas() : "                                     <<constants.gas()<<std::endl;
    std::cout<<std::setw(width)<<"radiationf() : "                              <<constants.radiationf()<<std::endl;
    std::cout<<std::setw(width)<<"radiations() : "                              <<constants.radiations()<<std::endl;
    std::cout<<std::setw(width)<<"faraday() : "                                 <<constants.faraday()<<std::endl;
    std::cout<<std::setw(width)<<"stefan() : "                                  <<constants.stefan()<<std::endl;

    // Atomic and nuclear constants
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Atomic and nuclear constants : "              <<std::endl;
    std::cout<<std::setw(width)<<"hfs() : "                                     <<constants.hfs()<<std::endl;
    std::cout<<std::setw(width)<<"alpha() : "                                   <<constants.alpha()<<std::endl;
    std::cout<<std::setw(width)<<"alpha2() : "                                  <<constants.alpha2()<<std::endl;
    std::cout<<std::setw(width)<<"rydberg() : "                                 <<constants.rydberg()<<std::endl;
    std::cout<<std::setw(width)<<"relectron() : "                               <<constants.relectron()<<std::endl;
    std::cout<<std::setw(width)<<"rbohr() : "                                   <<constants.rbohr()<<std::endl;
    std::cout<<std::setw(width)<<"hartree() : "                                 <<constants.hartree()<<std::endl;
    std::cout<<std::setw(width)<<"thomson() : "                                 <<constants.thomson()<<std::endl;

    // Planck units
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Planck units : "                              <<std::endl;
    std::cout<<std::setw(width)<<"planckq() : "                                 <<constants.planckq()<<std::endl;
    std::cout<<std::setw(width)<<"planckl() : "                                 <<constants.planckl()<<std::endl;
    std::cout<<std::setw(width)<<"planckm() : "                                 <<constants.planckm()<<std::endl;
    std::cout<<std::setw(width)<<"planckt() : "                                 <<constants.planckt()<<std::endl;
    std::cout<<std::setw(width)<<"planckf() : "                                 <<constants.planckf()<<std::endl;
    std::cout<<std::setw(width)<<"plancke() : "                                 <<constants.plancke()<<std::endl;
    std::cout<<std::setw(width)<<"planckp() : "                                 <<constants.planckp()<<std::endl;
    std::cout<<std::setw(width)<<"plancktheta() : "                             <<constants.plancktheta()<<std::endl;

    // Solar system masses
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Solar system masses : "                       <<std::endl;
    std::cout<<std::setw(width)<<"msun() : "                                    <<constants.msun()<<std::endl;
    std::cout<<std::setw(width)<<"mearth() : "                                  <<constants.mearth()<<std::endl;
    std::cout<<std::setw(width)<<"mmercury() : "                                <<constants.mmercury()<<std::endl;
    std::cout<<std::setw(width)<<"mvenus() : "                                  <<constants.mvenus()<<std::endl;
    std::cout<<std::setw(width)<<"mmars() : "                                   <<constants.mmars()<<std::endl;
    std::cout<<std::setw(width)<<"mjupiter() : "                                <<constants.mjupiter()<<std::endl;
    std::cout<<std::setw(width)<<"msaturn() : "                                 <<constants.msaturn()<<std::endl;
    std::cout<<std::setw(width)<<"muranus() : "                                 <<constants.muranus()<<std::endl;
    std::cout<<std::setw(width)<<"mneptune() : "                                <<constants.mneptune()<<std::endl;
    std::cout<<std::setw(width)<<"mmoon() : "                                   <<constants.mmoon()<<std::endl;

    // Solar system radii
    std::cout<<std::endl;
    std::cout<<std::setw(width)<<"Solar system radii : "                        <<std::endl;
    std::cout<<std::setw(width)<<"rsun() : "                                    <<constants.rsun()<<std::endl;
    std::cout<<std::setw(width)<<"rearth() : "                                  <<constants.rearth()<<std::endl;
    std::cout<<std::setw(width)<<"rmercury() : "                                <<constants.rmercury()<<std::endl;
    std::cout<<std::setw(width)<<"rvenus() : "                                  <<constants.rvenus()<<std::endl;
    std::cout<<std::setw(width)<<"rmars() : "                                   <<constants.rmars()<<std::endl;
    std::cout<<std::setw(width)<<"rjupiter() : "                                <<constants.rjupiter()<<std::endl;
    std::cout<<std::setw(width)<<"rsaturn() : "                                 <<constants.rsaturn()<<std::endl;
    std::cout<<std::setw(width)<<"ruranus() : "                                 <<constants.ruranus()<<std::endl;
    std::cout<<std::setw(width)<<"rneptune() : "                                <<constants.rneptune()<<std::endl;
    std::cout<<std::setw(width)<<"rmoon() : "                                   <<constants.rmoon()<<std::endl;

    // Finalize
    std::cout<<std::noboolalpha<<std::right<<std::endl;
    std::cout<<"END = Constants::example()"<<std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //



/*////////////////////////////////////////////////////////////////////////////*/
} // namespace
#endif // CONSTANTS_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
