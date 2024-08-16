/**
 * @file GSIRateLawGammaT.cpp
 *
 * @brief Class which computes the reaction rate constant for a gas 
 *        independent ablation surface reaction based on an Arrhenius formula.
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

class GSIRateLawLHAblaArrhenius : public GSIRateLaw
{
public:
    GSIRateLawLHAblaArrhenius(ARGS args)
        : GSIRateLaw(args),
          pos_T_trans(0)
    {
        assert(args.s_node_rate_law.tag() == "LH_abla_arrhenius");

        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
            "The sticking coeffcient probability for the reaction "
            "should be provided.");
        args.s_node_rate_law.getAttribute( "T", m_T_act,
            "The activation temperature for the reaction "
            "should be provided for an adsorption reaction.");
    }

//==============================================================================

    ~GSIRateLawLHAblaArrhenius( ){ }

//==============================================================================

    double forwardReactionRateCoefficient(
        const Eigen::VectorXd& v_rhoi, const Eigen::VectorXd& v_Tsurf) const
    {
    	const double Tsurf = v_Tsurf(pos_T_trans);

	return m_pre_exp*exp(-m_T_act / Tsurf); 
    }

private:
    const size_t pos_T_trans;

    double m_pre_exp;
    double m_T_act;
};

ObjectProvider<
    GSIRateLawLHAblaArrhenius, GSIRateLaw>
    gsi_rate_law_lh_abla_arrhenius("LH_abla_arrhenius");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
