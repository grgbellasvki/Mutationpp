/**
 * @file GSIRateLawPyrolysisGoldstein.cpp
 *
 * @brief
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


#include "AutoRegistration.h"
#include "Utilities.h"

#include "GSIRateLaw.h"
#include "SurfaceState.h"
#include "SolidProperties.h"

using namespace Eigen;

using namespace Mutation::Utilities;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIRateLawPyrolysisGoldstein: public GSIRateLaw
{
public:
    GSIRateLawPyrolysisGoldstein(ARGS args)
        : GSIRateLaw(args),
          m_surf_state(args.s_surf_state),
          pos_T_trans(0),
          pos_rho_s(args.s_reactants[0])
    {
        // assert(args.s_node_rate_law.tag() == "goldstein");
        args.s_node_rate_law.getAttribute( "pre_exp", m_pre_exp,
            "The pre-exponential coefficient for the reaction "
            "should be provided with goldstein pyrolysis.");
        args.s_node_rate_law.getAttribute( "T", m_T_act,
            "The activation temperature for the reaction "
            "should be provided with goldstein pyrolysis.");
        args.s_node_rate_law.getAttribute( "order", m_ord, 1);

        m_rho_i =
            m_surf_state.solidProps().
                getPyrolysingSolidInitialDensity(pos_rho_s);
        m_rho_f =
            m_surf_state.solidProps().
                getPyrolysingSolidFinalDensity(pos_rho_s);
    }

//==============================================================================

    ~GSIRateLawPyrolysisGoldstein(){}

//==============================================================================

    double forwardReactionRateCoefficient(
        const VectorXd& v_rhoi, const VectorXd& v_Tsurf) const
    {
        double rho_s =
            m_surf_state.solidProps().getPyrolysingSolidDensity(pos_rho_s);
        double frac = pow((rho_s-m_rho_f) / m_rho_i, m_ord);

        return -m_pre_exp * m_rho_i * exp(
            -m_T_act / v_Tsurf(pos_T_trans)) * frac;
    }

private:
    const SurfaceState& m_surf_state;

    const size_t pos_T_trans;
    const int pos_rho_s;

    double m_rho_i;
    double m_rho_f;

    int m_ord;
    double m_pre_exp;
    double m_T_act;
};

ObjectProvider<
    GSIRateLawPyrolysisGoldstein, GSIRateLaw>
    gsi_rate_law_pyrolysis_goldstein("goldstein");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
