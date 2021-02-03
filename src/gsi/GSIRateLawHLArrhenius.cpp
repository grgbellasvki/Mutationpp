/**
 * @file GSIRateLawGammaT.cpp
 *
 * @brief Class which computes the reaction rate constant for an adsorption
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
          pos_gas_r(0),
          pos_site_r(1)
    {
        assert(args.s_node_rate_law.tag() == "LH_arrhenius");

        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
            "The sticking coeffcient probability for the reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "T", m_T_act,
            "The activation temperature for the reaction "
            "should be provided for an adsorption reaction.");

        // For the gas in the reactants
        m_idx_gas = mv_react[pos_gas_r];
        //std::cout << "m_idx_gas is " << m_idx_gas << std::endl; //it is number 10
        // Error if m_idx_gas > ns

        // For the sites
        int idx_site = mv_react[pos_site_r];
        // Error if idx_site > ns

        m_site_categ = m_surf_props.siteSpeciesToSiteCategoryIndex(idx_site);
        m_n_sites = m_surf_props.nSiteDensityInCategory(m_site_categ);
    }

//==============================================================================

    ~GSIRateLawLHArrhenius( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Tsurf) const
    {
    	const double Tsurf = v_Tsurf(pos_T_trans);

    	const int set_state_with_rhoi_T = 1;
        m_thermo.setState(v_rhoi.data(), v_Tsurf.data(), set_state_with_rhoi_T);
        //const double thermal_speed =
        //    m_transport.speciesThermalSpeed(m_idx_gas);

	//std::cout << "th speed is" << thermal_speed << std::endl;
	double m_n = 0.0140067 / NA; //hardcoded for now
	//return m_pre_exp *sqrt(1.0/m_n_sites) * sqrt(PI*KB*v_Tsurf[0]/ (8.0*m_n)) * exp(-21000.0/ v_Tsurf[0]);
	//return m_pre_exp *sqrt(1.0/m_n_sites) * sqrt(PI*KB*v_Tsurf[0]/ (8.0*m_n)) * exp(-m_T_act/ v_Tsurf[0]);
	return m_pre_exp *sqrt(1.0/m_n_sites) * sqrt(PI*KB*v_Tsurf[0]/ (2.0*m_n)) * exp(-m_T_act/ v_Tsurf[0]);

        //return m_pre_exp * thermal_speed * sqrt(1.0/m_n_sites)/
        //    (4.) * exp(-m_T_act / Tsurf);
    }

private:
    const size_t pos_T_trans;
    const size_t pos_gas_r;
    const size_t pos_site_r;

    int m_idx_gas;
    int m_site_categ;
    double m_n_sites;

    double m_pre_exp;
    double m_T_act;

    const std::vector<int>& mv_react;
    const SurfaceProperties& m_surf_props;
};

ObjectProvider<
    GSIRateLawLHArrhenius, GSIRateLaw>
    gsi_rate_law_lh_arrhenius("LH_arrhenius");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
