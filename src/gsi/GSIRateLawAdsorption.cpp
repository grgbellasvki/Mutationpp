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

using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawAdsorption : public GSIRateLaw
{
public:
    GSIRateLawAdsorption(ARGS args)
        : GSIRateLaw(args),
          m_surf_props(args.s_surf_props),
          mv_react(args.s_reactants),
          pos_T_trans(0),
          idx_react(0)
    {
        assert(args.s_node_rate_law.tag() == "adsorption");

        const double huge_number = 1.e20;
        const double zero = 0.0;

        args.s_node_rate_law.getAttribute( "stick_coef", m_stick_coef,
            "The sticking coeffcient probability for the reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "T", m_T_act,
            "The activation temperature for the reaction "
            "should be provided for an adsorption reaction.");
        args.s_node_rate_law.getAttribute( "stick_coef_corr",
            m_beta_T_corr, zero);
        args.s_node_rate_law.getAttribute("T_thresh", m_T_thresh, huge_number);

        const int n_stick_sp = 1;
        m_stick_coef_power = mv_react.size() - n_stick_sp;
    }

//==============================================================================

    ~GSIRateLawAdsorption( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Tsurf) const
    {
/*    	double T_surf = v_Tsurf(pos_T_trans);

    	const int set_state_with_rhoi_T = 1;
        m_thermo.setState(v_rhoi.data(), v_Twall.data(), set_state_with_rhoi_T);
        double thermal_speed =
            m_transport.speciesThermalSpeed(mv_react[idx_react]);

        //  Sticking coefficient temperature correction
        double S_coef_T_corr = m_S_coeff;
        if (T_surf > m_T_thresh)
        	S_coef_T_corr *= exp(-m_beta_T_corr * (T_surf - m_T_thresh));

        double B = pow(
            m_surf_props.nSitesInCategory(site_categ), m_stick_coef_power);

       return S_coef_T_corr * thermal_speed / (4 * B) * exp(-m_T_act / T_surf)); */
       return 0.;
    }

private:
    const size_t pos_T_trans;
    const size_t idx_react;

    double m_stick_coef;
    double m_T_act;

    double m_T_thresh;
    double m_beta_T_corr;

    int m_stick_coef_power;

    const std::vector<int>& mv_react;
    const SurfaceProperties& m_surf_props;
};

ObjectProvider<
    GSIRateLawAdsorption, GSIRateLaw>
    gsi_rate_law_adsorption("adsorption");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
