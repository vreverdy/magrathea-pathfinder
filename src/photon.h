/* ********************************* PHOTON ********************************* */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          Photon
// DESCRIPTION :    Photon step implementation for raytracing
// AUTHOR(S) :      Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Vincent Reverdy (2012-2013), Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           photon.h
/// \brief          Photon step implementation for raytracing
/// \author         Vincent Reverdy (vince.rev@gmail.com), Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2012-2013, 2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
#ifndef PHOTON_H_INCLUDED
#define PHOTON_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/

// ------------------------------ PREPROCESSOR ------------------------------ //
// Include C++
#include <array>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <utility>
// Include libs
// Include project
#include "magrathea/abstractstep.h"
// Misc
// -------------------------------------------------------------------------- //

// ---------------------------------- CLASS --------------------------------- //
// Photon step implementation for raytracing
/// \brief          Photon step implementation for raytracing.
/// \details        A step of integration for a given photon.
/// \tparam         Type Data type.
/// \tparam         Dimension Number of space dimension.
template <typename Type = double, unsigned int Dimension = 3>
class Photon final
    : public magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>> {
    // Setup
public:
    using magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>::operator=;

    // Lifecycle
    /// \name           Lifecycle
    //@{
public:
    template <class... Misc>
    explicit inline Photon(Misc &&...misc);
    //@}

    // Data
    /// \name           Data
    //@{
public:
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template id<Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template index(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template id<Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template index(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template a(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template a(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template t(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template t(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<2, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template x(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<2, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template x(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<3, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template y(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<3, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template y(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<4, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template z(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<4, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template z(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<5, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dtdl(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<5, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dtdl(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<6, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dxdl(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<6, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dxdl(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<7, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dydl(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<7, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dydl(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<8, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dzdl(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template core<8, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dzdl(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template level(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<0, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template level(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template ah(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<1, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template ah(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<2, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template rho(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<2, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template rho(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<3, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template phi(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<3, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template phi(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<4, Values...>(0, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidx(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<4, Values...>(0, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidx(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<4, Values...>(1, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidy(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<4, Values...>(1, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidy(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<4, Values...>(2, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidz(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<4, Values...>(2, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidz(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<5, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidl(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<5, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidl(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<6, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template laplacian(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<6, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template laplacian(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<7, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template redshift(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<7, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template redshift(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<8, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dsdl2(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<8, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dsdl2(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<9, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template error(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<9, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template error(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<10, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template distance(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<10, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template distance(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<11, Values...>(0, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template isw(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<11, Values...>(0, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template isw(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<11, Values...>(1, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template iswold(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<11, Values...>(1, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template iswold(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<11, Values...>(2, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template chi(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<11, Values...>(2, std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template chi(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<12, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template lambda(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<12, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template lambda(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<13, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidt(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<13, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template dphidt(Misc &&...misc) const;
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<14, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template s(Misc &&...misc);
    template <unsigned int... Values, class... Misc, class Template = decltype(std::declval<const magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>>().template extra<14, Values...>(std::declval<Misc>()...)), class = typename std::enable_if<!std::is_void<Template>::value>::type>
    inline Template s(Misc &&...misc) const;

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
template <typename Type, unsigned int Dimension>
template <class... Misc>
inline Photon<Type, Dimension>::Photon(Misc &&...misc)
    : magrathea::AbstractStep<Photon<Type, Dimension>, unsigned int, std::array<Type, 1 + (1 + Dimension) * 2>, std::tuple<Type, Type, Type, Type, std::array<Type, Dimension>, Type, Type, Type, Type, Type, Type, std::array<Type, 3>, Type, Type, Type>>(std::forward<Misc>(misc)...) {
    ;
}
// -------------------------------------------------------------------------- //

// ---------------------------------- DATA ---------------------------------- //
// Access to the index data
/// \brief          Access to the index data.
/// \details        Provides an access to the index data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::index(Misc &&...misc) {
    return this->template id<Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the index data
/// \brief          Immutable access to the index data.
/// \details        Provides an immutable access to the index data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::index(Misc &&...misc) const {
    return this->template id<Values...>(std::forward<Misc>(misc)...);
}

// Access to the a data
/// \brief          Access to the a data.
/// \details        Provides an access to the a data by forwarding parameters
///                 to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::a(Misc &&...misc) {
    return this->template core<0, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the a data
/// \brief          Immutable access to the a data.
/// \details        Provides an immutable access to the a data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::a(Misc &&...misc) const {
    return this->template core<0, Values...>(std::forward<Misc>(misc)...);
}

// Access to the t data
/// \brief          Access to the t data.
/// \details        Provides an access to the t data by forwarding parameters
///                 to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::t(Misc &&...misc) {
    return this->template core<1, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the t data
/// \brief          Immutable access to the t data.
/// \details        Provides an immutable access to the t data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::t(Misc &&...misc) const {
    return this->template core<1, Values...>(std::forward<Misc>(misc)...);
}

// Access to the x data
/// \brief          Access to the x data.
/// \details        Provides an access to the x data by forwarding parameters
///                 to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::x(Misc &&...misc) {
    return this->template core<2, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the x data
/// \brief          Immutable access to the x data.
/// \details        Provides an immutable access to the x data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::x(Misc &&...misc) const {
    return this->template core<2, Values...>(std::forward<Misc>(misc)...);
}

// Access to the y data
/// \brief          Access to the y data.
/// \details        Provides an access to the y data by forwarding parameters
///                 to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::y(Misc &&...misc) {
    return this->template core<3, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the y data
/// \brief          Immutable access to the y data.
/// \details        Provides an immutable access to the y data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::y(Misc &&...misc) const {
    return this->template core<3, Values...>(std::forward<Misc>(misc)...);
}

// Access to the z data
/// \brief          Access to the z data.
/// \details        Provides an access to the z data by forwarding parameters
///                 to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::z(Misc &&...misc) {
    return this->template core<4, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the z data
/// \brief          Immutable access to the z data.
/// \details        Provides an immutable access to the z data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::z(Misc &&...misc) const {
    return this->template core<4, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dtdl data
/// \brief          Access to the dtdl data.
/// \details        Provides an access to the dtdl data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dtdl(Misc &&...misc) {
    return this->template core<5, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dtdl data
/// \brief          Immutable access to the dtdl data.
/// \details        Provides an immutable access to the dtdl data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dtdl(Misc &&...misc) const {
    return this->template core<5, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dxdl data
/// \brief          Access to the dxdl data.
/// \details        Provides an access to the dxdl data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dxdl(Misc &&...misc) {
    return this->template core<6, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dxdl data
/// \brief          Immutable access to the dxdl data.
/// \details        Provides an immutable access to the dxdl data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dxdl(Misc &&...misc) const {
    return this->template core<6, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dydl data
/// \brief          Access to the dydl data.
/// \details        Provides an access to the dydl data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dydl(Misc &&...misc) {
    return this->template core<7, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dydl data
/// \brief          Immutable access to the dydl data.
/// \details        Provides an immutable access to the dydl data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dydl(Misc &&...misc) const {
    return this->template core<7, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dzdl data
/// \brief          Access to the dzdl data.
/// \details        Provides an access to the dzdl data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dzdl(Misc &&...misc) {
    return this->template core<8, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dzdl data
/// \brief          Immutable access to the dzdl data.
/// \details        Provides an immutable access to the dzdl data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dzdl(Misc &&...misc) const {
    return this->template core<8, Values...>(std::forward<Misc>(misc)...);
}

// Access to the level data
/// \brief          Access to the level data.
/// \details        Provides an access to the level data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::level(Misc &&...misc) {
    return this->template extra<0, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the level data
/// \brief          Immutable access to the level data.
/// \details        Provides an immutable access to the level data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::level(Misc &&...misc) const {
    return this->template extra<0, Values...>(std::forward<Misc>(misc)...);
}

// Access to the ah data
/// \brief          Access to the ah data.
/// \details        Provides an access to the ah data by forwarding parameters
///                 to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::ah(Misc &&...misc) {
    return this->template extra<1, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the ah data
/// \brief          Immutable access to the ah data.
/// \details        Provides an immutable access to the ah data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::ah(Misc &&...misc) const {
    return this->template extra<1, Values...>(std::forward<Misc>(misc)...);
}

// Access to the rho data
/// \brief          Access to the rho data.
/// \details        Provides an access to the rho data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::rho(Misc &&...misc) {
    return this->template extra<2, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the rho data
/// \brief          Immutable access to the rho data.
/// \details        Provides an immutable access to the rho data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::rho(Misc &&...misc) const {
    return this->template extra<2, Values...>(std::forward<Misc>(misc)...);
}

// Access to the phi data
/// \brief          Access to the phi data.
/// \details        Provides an access to the phi data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::phi(Misc &&...misc) {
    return this->template extra<3, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the phi data
/// \brief          Immutable access to the phi data.
/// \details        Provides an immutable access to the phi data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::phi(Misc &&...misc) const {
    return this->template extra<3, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dphidx data
/// \brief          Access to the dphidx data.
/// \details        Provides an access to the dphidx data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidx(Misc &&...misc) {
    return this->template extra<4, Values...>(0, std::forward<Misc>(misc)...);
}

// Immutable access to the dphidx data
/// \brief          Immutable access to the dphidx data.
/// \details        Provides an immutable access to the dphidx data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidx(Misc &&...misc) const {
    return this->template extra<4, Values...>(0, std::forward<Misc>(misc)...);
}

// Access to the dphidy data
/// \brief          Access to the dphidy data.
/// \details        Provides an access to the dphidy data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidy(Misc &&...misc) {
    return this->template extra<4, Values...>(1, std::forward<Misc>(misc)...);
}

// Immutable access to the dphidy data
/// \brief          Immutable access to the dphidy data.
/// \details        Provides an immutable access to the dphidy data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidy(Misc &&...misc) const {
    return this->template extra<4, Values...>(1, std::forward<Misc>(misc)...);
}

// Access to the dphidz data
/// \brief          Access to the dphidz data.
/// \details        Provides an access to the dphidz data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidz(Misc &&...misc) {
    return this->template extra<4, Values...>(2, std::forward<Misc>(misc)...);
}

// Immutable access to the dphidz data
/// \brief          Immutable access to the dphidz data.
/// \details        Provides an immutable access to the dphidz data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidz(Misc &&...misc) const {
    return this->template extra<4, Values...>(2, std::forward<Misc>(misc)...);
}

// Access to the dphidl data
/// \brief          Access to the dphidl data.
/// \details        Provides an access to the dphidl data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidl(Misc &&...misc) {
    return this->template extra<5, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dphidl data
/// \brief          Immutable access to the dphidl data.
/// \details        Provides an immutable access to the dphidl data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidl(Misc &&...misc) const {
    return this->template extra<5, Values...>(std::forward<Misc>(misc)...);
}

// Access to the laplacian data
/// \brief          Access to the laplacian data.
/// \details        Provides an access to the laplacian data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::laplacian(Misc &&...misc) {
    return this->template extra<6, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the laplacian data
/// \brief          Immutable access to the laplacian data.
/// \details        Provides an immutable access to the laplacian data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::laplacian(Misc &&...misc) const {
    return this->template extra<6, Values...>(std::forward<Misc>(misc)...);
}

// Access to the redshift data
/// \brief          Access to the redshift data.
/// \details        Provides an access to the redshift data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::redshift(Misc &&...misc) {
    return this->template extra<7, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the redshift data
/// \brief          Immutable access to the redshift data.
/// \details        Provides an immutable access to the redshift data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::redshift(Misc &&...misc) const {
    return this->template extra<7, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dsdl2 data
/// \brief          Access to the dsdl2 data.
/// \details        Provides an access to the dsdl2 data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dsdl2(Misc &&...misc) {
    return this->template extra<8, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dsdl2 data
/// \brief          Immutable access to the dsdl2 data.
/// \details        Provides an immutable access to the dsdl2 data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dsdl2(Misc &&...misc) const {
    return this->template extra<8, Values...>(std::forward<Misc>(misc)...);
}

// Access to the error data
/// \brief          Access to the error data.
/// \details        Provides an access to the error data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::error(Misc &&...misc) {
    return this->template extra<9, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the error data
/// \brief          Immutable access to the error data.
/// \details        Provides an immutable access to the error data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::error(Misc &&...misc) const {
    return this->template extra<9, Values...>(std::forward<Misc>(misc)...);
}

// Access to the distance data
/// \brief          Access to the distance data.
/// \details        Provides an access to the distance data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::distance(Misc &&...misc) {
    return this->template extra<10, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the distance data
/// \brief          Immutable access to the distance data.
/// \details        Provides an immutable access to the distance data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::distance(Misc &&...misc) const {
    return this->template extra<10, Values...>(std::forward<Misc>(misc)...);
}

// Access to the isw data
/// \brief          Access to the isw data.
/// \details        Provides an access to the isw data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::isw(Misc &&...misc) {
    return this->template extra<11, Values...>(0, std::forward<Misc>(misc)...);
}

// Immutable access to the isw data
/// \brief          Immutable access to the isw data.
/// \details        Provides an immutable access to the isw data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::isw(Misc &&...misc) const {
    return this->template extra<11, Values...>(0, std::forward<Misc>(misc)...);
}

// Access to the iswold data
/// \brief          Access to the iswold data.
/// \details        Provides an access to the iswold data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::iswold(Misc &&...misc) {
    return this->template extra<11, Values...>(1, std::forward<Misc>(misc)...);
}

// Immutable access to the iswold data
/// \brief          Immutable access to the iswold data.
/// \details        Provides an immutable access to the iswold data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::iswold(Misc &&...misc) const {
    return this->template extra<11, Values...>(1, std::forward<Misc>(misc)...);
}

// Access to the chi data
/// \brief          Access to the chi data.
/// \details        Provides an access to the chi data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::chi(Misc &&...misc) {
    return this->template extra<11, Values...>(2, std::forward<Misc>(misc)...);
}

// Immutable access to the chi data
/// \brief          Immutable access to the chi data.
/// \details        Provides an immutable access to the chi data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::chi(Misc &&...misc) const {
    return this->template extra<11, Values...>(2, std::forward<Misc>(misc)...);
}

// Access to the lambda data
/// \brief          Access to the lambda data.
/// \details        Provides an access to the lambda data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::lambda(Misc &&...misc) {
    return this->template extra<12, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the lambda data
/// \brief          Immutable access to the lambda data.
/// \details        Provides an immutable access to the lambda data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::lambda(Misc &&...misc) const {
    return this->template extra<12, Values...>(std::forward<Misc>(misc)...);
}

// Access to the dphidt data
/// \brief          Access to the dphidt data.
/// \details        Provides an access to the dphidt data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidt(Misc &&...misc) {
    return this->template extra<13, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the dphidt data
/// \brief          Immutable access to the dphidt data.
/// \details        Provides an immutable access to the dphidt data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::dphidt(Misc &&...misc) const {
    return this->template extra<13, Values...>(std::forward<Misc>(misc)...);
}

// Access to the s data
/// \brief          Access to the s data.
/// \details        Provides an access to the s data by forwarding
///                 parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::s(Misc &&...misc) {
    return this->template extra<14, Values...>(std::forward<Misc>(misc)...);
}

// Immutable access to the s data
/// \brief          Immutable access to the s data.
/// \details        Provides an immutable access to the s data by
///                 forwarding parameters to the unified base accessor member.
/// \tparam         Values List of template values.
/// \tparam         Misc (Miscellaneous types.)
/// \tparam         Template (Deduced template type.)
/// \param[in]      misc Miscellaneous arguments.
/// \return         Forwarded result.
template <typename Type, unsigned int Dimension>
template <unsigned int... Values, class... Misc, class Template, class>
inline Template Photon<Type, Dimension>::s(Misc &&...misc) const {
    return this->template extra<14, Values...>(std::forward<Misc>(misc)...);
}

// -------------------------------------------------------------------------- //

// ---------------------------------- TEST ---------------------------------- //
// Example function
/// \brief          Example function.
/// \details        Tests and demonstrates the use of Photon.
/// \return         0 if no error.
template <typename Type, unsigned int Dimension>
int Photon<Type, Dimension>::example() {
    // Initialize
    std::cout << "BEGIN = Photon::example()" << std::endl;
    std::cout << std::boolalpha << std::left;
    const unsigned int width = 40;
    std::array<double, 9> array({{0, 1, 2, 3, 4, 5, 6, 7, 8}});
    std::array<double, 3> array2({{8, 9, 10}});

    // Construction
    Photon<double, 3> photon(42, array, std::make_tuple(11, 12, 13, 14, array2, 15, 16, 17, 18, 19, 20, array2));

    // Lifecycle
    std::cout << std::endl;
    std::cout << std::setw(width * 3) << "Lifecycle : " << std::endl;
    std::cout << std::setw(width * 3) << "Photon<>() : " << Photon<>() << std::endl;
    std::cout << std::setw(width * 3) << "Photon<double, 3>(42, array, std::make_tuple(11, 12, 13, 14, array2, 15, 16, 17, 18, 19, 20, array2)) : " << Photon<double, 3>(42, array, std::make_tuple(11, 12, 13, 14, array2, 15, 16, 17, 18, 19, 20, array2)) << std::endl;
    std::cout << std::setw(width * 3) << "Photon<double, 3>(photon) : " << Photon<double, 3>(photon) << std::endl;

    // Data
    std::cout << std::endl;
    std::cout << std::setw(width) << "Data : " << std::endl;
    std::cout << std::setw(width) << "photon.index() = 42 : " << (photon.index() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.index() : " << photon.index() << std::endl;
    std::cout << std::setw(width) << "photon.a() = 42 : " << (photon.a() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.a() : " << photon.a() << std::endl;
    std::cout << std::setw(width) << "photon.t() = 42 : " << (photon.t() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.t() : " << photon.t() << std::endl;
    std::cout << std::setw(width) << "photon.x() = 42 : " << (photon.x() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.x() : " << photon.x() << std::endl;
    std::cout << std::setw(width) << "photon.y() = 42 : " << (photon.y() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.y() : " << photon.y() << std::endl;
    std::cout << std::setw(width) << "photon.z() = 42 : " << (photon.z() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.z() : " << photon.z() << std::endl;
    std::cout << std::setw(width) << "photon.dtdl() = 42 : " << (photon.dtdl() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dtdl() : " << photon.dtdl() << std::endl;
    std::cout << std::setw(width) << "photon.dxdl() = 42 : " << (photon.dxdl() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dxdl() : " << photon.dxdl() << std::endl;
    std::cout << std::setw(width) << "photon.dydl() = 42 : " << (photon.dydl() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dydl() : " << photon.dydl() << std::endl;
    std::cout << std::setw(width) << "photon.dzdl() = 42 : " << (photon.dzdl() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dzdl() : " << photon.dzdl() << std::endl;
    std::cout << std::setw(width) << "photon.level() = 42 : " << (photon.level() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.level() : " << photon.level() << std::endl;
    std::cout << std::setw(width) << "photon.ah() = 42 : " << (photon.ah() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.ah() : " << photon.ah() << std::endl;
    std::cout << std::setw(width) << "photon.rho() = 42 : " << (photon.rho() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.rho() : " << photon.rho() << std::endl;
    std::cout << std::setw(width) << "photon.phi() = 42 : " << (photon.phi() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.phi() : " << photon.phi() << std::endl;
    std::cout << std::setw(width) << "photon.dphidx() = 42 : " << (photon.dphidx() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dphidx() : " << photon.dphidx() << std::endl;
    std::cout << std::setw(width) << "photon.dphidy() = 42 : " << (photon.dphidy() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dphidy() : " << photon.dphidy() << std::endl;
    std::cout << std::setw(width) << "photon.dphidz() = 42 : " << (photon.dphidz() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dphidz() : " << photon.dphidz() << std::endl;
    std::cout << std::setw(width) << "photon.dphidl() = 42 : " << (photon.dphidl() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dphidl() : " << photon.dphidl() << std::endl;
    std::cout << std::setw(width) << "photon.laplacian() = 42 : " << (photon.laplacian() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.laplacian() : " << photon.laplacian() << std::endl;
    std::cout << std::setw(width) << "photon.redshift() = 42 : " << (photon.redshift() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.redshift() : " << photon.redshift() << std::endl;
    std::cout << std::setw(width) << "photon.dsdl2() = 42 : " << (photon.dsdl2() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.dsdl2() : " << photon.dsdl2() << std::endl;
    std::cout << std::setw(width) << "photon.error() = 42 : " << (photon.error() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.error() : " << photon.error() << std::endl;
    std::cout << std::setw(width) << "photon.distance() = 42 : " << (photon.distance() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.distance() : " << photon.distance() << std::endl;
    std::cout << std::setw(width) << "photon.isw() = 42 : " << (photon.isw() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.isw() : " << photon.isw() << std::endl;
    std::cout << std::setw(width) << "photon.iswold() = 42 : " << (photon.iswold() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.iswold() : " << photon.iswold() << std::endl;
    std::cout << std::setw(width) << "photon.chi() = 42 : " << (photon.chi() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.chi() : " << photon.chi() << std::endl;
    std::cout << std::setw(width) << "photon.lambda() = 42 : " << (photon.lambda() = 42) << std::endl;
    std::cout << std::setw(width) << "photon.lambda() : " << photon.lambda() << std::endl;

    // Finalize
    std::cout << std::noboolalpha << std::right << std::endl;
    std::cout << "END = Photon::example()" << std::endl;
    return 0;
}
// -------------------------------------------------------------------------- //

/*////////////////////////////////////////////////////////////////////////////*/
#endif // PHOTON_H_INCLUDED
/*////////////////////////////////////////////////////////////////////////////*/
