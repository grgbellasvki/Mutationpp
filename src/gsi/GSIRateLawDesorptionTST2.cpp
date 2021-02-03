/**
 * @file GSIRateLawGammaT.cpp
 *
 * @brief Class which computes the reaction rate constant for a desorption
 *        surface reaction based on simple transition state theory.
 */

/*
 * Copyright 2018 von Karman Institute for Fluid Dynamics (VKI)
 *
 * This file is part of MUlticomponent Thermodynamic And Transport
 * properties for IONized gases in C++ (Mutation++) software package.
 *
 * Mutation++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Mutation++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Mutation++.  If not, see
 * <http://www.gnu.org/licenses/>.
 */


#include "Thermodynamics.h"
#include "Transport.h"

#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"
#include "SurfaceProperties.h"

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawDesorptionTST2 : public GSIRateLaw
{
public:
    GSIRateLawDesorptionTST2(ARGS args)
        : GSIRateLaw(args),
          m_surf_props(args.s_surf_props),
          mv_react(args.s_reactants),
          pos_T_trans(0),
          idx_react(0)
    {
        assert(args.s_node_rate_law.tag() == "desorption_tst2");

        args.s_node_rate_law.getAttribute("T", m_T_des,
            "Activation temperature should be provided for every desorption "
            "reaction.");

        m_mass_des = m_thermo.speciesMw(
            args.s_surf_props.surfaceToGasIndex(mv_react[idx_react])) / NA;
        m_site_categ = m_surf_props.siteSpeciesToSiteCategoryIndex(
            mv_react[idx_react]);
        m_n_sites = m_surf_props.nSiteDensityInCategory(m_site_categ);
        // m_n_sites = args.s_surf_props.nSiteDensity();
    }

//==============================================================================

    ~GSIRateLawDesorptionTST2( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Tsurf) const
    {
    	const double Tsurf = v_Tsurf(pos_T_trans);

    	const double pre_exp = 2 * PI * m_mass_des * KB * KB * Tsurf
            / (HP * HP * HP * m_n_sites);

        return pre_exp * exp(-m_T_des / Tsurf);
    }

private:
    const size_t idx_react;
    int m_site_categ;

    double m_T_des;

    double m_n_sites;
    double m_mass_des;

    const double pos_T_trans;

    const std::vector<int>& mv_react;
    const SurfaceProperties& m_surf_props;
};

ObjectProvider<
    GSIRateLawDesorptionTST2, GSIRateLaw>
    gsi_rate_law_desorption_tst2("desorption_tst2");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
