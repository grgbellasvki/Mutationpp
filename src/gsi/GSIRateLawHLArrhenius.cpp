/**
 * @file GSIRateLawGammaT.cpp
 *
 * @brief Class which computes the reaction rate constant for an HL
 *        surface reaction based on an Arrhenius formula.
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

class GSIRateLawLHArrhenius : public GSIRateLaw
{
public:
    GSIRateLawLHArrhenius(ARGS args)
        : GSIRateLaw(args),
          m_surf_props(args.s_surf_props),
          mv_react(args.s_reactants),
          pos_T_trans(0),
          pos_site_r(0)
    {
        assert(args.s_node_rate_law.tag() == "LH_arrhenius");

        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
            "The sticking coeffcient probability for the reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "T", m_T_act,
            "The activation temperature for the reaction "
            "should be provided for an adsorption reaction.");

        // For the sites
        m_idx_site = mv_react[pos_site_r];

        m_site_categ = m_surf_props.siteSpeciesToSiteCategoryIndex(m_idx_site);
        m_n_sites = m_surf_props.nSiteDensityInCategory(m_site_categ);

        int gas_index = m_surf_props.surfaceToGasIndex(m_idx_site);
        m_mw_r = m_thermo.speciesMw()[gas_index] / NA;
    }

//==============================================================================

    ~GSIRateLawLHArrhenius( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Tsurf) const
    {
    	const double Tsurf = v_Tsurf(pos_T_trans);

	return m_pre_exp *sqrt(1.0/m_n_sites) * sqrt(PI*KB*Tsurf/ (2.0*m_mw_r)) * exp(-m_T_act/ Tsurf);

    }

private:
    const size_t pos_T_trans;
    const size_t pos_site_r;

    int m_idx_site;
    int m_site_categ;
    double m_n_sites;

    double m_pre_exp;
    double m_T_act;
    
    double m_mw_r;

    const std::vector<int>& mv_react;
    const SurfaceProperties& m_surf_props;
};

ObjectProvider<
    GSIRateLawLHArrhenius, GSIRateLaw>
    gsi_rate_law_lh_arrhenius("LH_arrhenius");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
