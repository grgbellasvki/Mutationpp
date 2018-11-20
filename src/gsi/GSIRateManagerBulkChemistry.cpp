/**
 * @file GSIRateBulkChemistry.cpp
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


#include "Thermodynamics.h"
#include "Transport.h"

#include "GSIReaction.h"
#include "GSIRateLaw.h"
#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"
#include "SolidProperties.h"
#include "SurfaceState.h"

using namespace Eigen;

using namespace Mutation::Utilities::Config;
using namespace Mutation::Thermodynamics;

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerBulkChemistry : public GSIRateManager
{
public:
    GSIRateManagerBulkChemistry(DataGSIRateManager args)
        : GSIRateManager(args),
		  m_ns(args.s_thermo.nSpecies()),
		  m_nr(args.s_reactions.size()),
          mv_work(m_ns),
          mv_chem_rate(m_ns),
          mv_react_rate_const(m_nr),
          mv_r_to_ps(m_nr),
          mv_r_to_pg(m_nr)
    {
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            int pos_pg = args.s_reactions[i_r]->getProducts().size()-1;
            int pos_ps = 0;
            // Ensure that the pyro gas is in the last position
            mv_r_to_pg[i_r] = args.s_reactions[i_r]->getProducts()[pos_pg];
            mv_r_to_ps[i_r] = args.s_reactions[i_r]->getReactants()[pos_ps];
            // Here add the intermediate solid products!
        }
    }

//=============================================================================

    ~GSIRateManagerBulkChemistry(){}

//=============================================================================

    Eigen::VectorXd computeRates()
    {
        mv_chem_rate.setZero();
        double Tsurf = m_surf_state.getSurfaceT()(0);
        double Psurf = m_thermo.P();

        // Get reaction rate constant
        double react_rate_const;
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            react_rate_const =
                v_reactions[i_r]->getRateLaw()->
                    forwardReactionRateCoefficient(
                        m_surf_state.getSurfaceRhoi(),
                        m_surf_state.getSurfaceT());

            m_surf_state.solidProps().getPyrolysingGasEquilMassFrac(
                mv_r_to_pg[i_r], Psurf, Tsurf, mv_work);

            mv_chem_rate -= react_rate_const * mv_work;
        }
        return mv_chem_rate;
    }

//=============================================================================

    Eigen::VectorXd computeRatesPerReaction()
    {
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_react_rate_const(i_r) =
                v_reactions[i_r]->getRateLaw()->
                    forwardReactionRateCoefficient(
                        m_surf_state.getSurfaceRhoi(),
                        m_surf_state.getSurfaceT());
        }
        return -mv_react_rate_const;
    }

//=============================================================================

    void computeRatesGasAndSolid(VectorXd& v_chem_rate_per_gas_and_solid)
    {
        v_chem_rate_per_gas_and_solid.setZero();
        double Tsurf = m_surf_state.getSurfaceT()(0);
        double Psurf = m_thermo.P();

        // Get reaction rate constant
        double reac_rate_const;
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            reac_rate_const =
                v_reactions[i_r]->getRateLaw()->
                    forwardReactionRateCoefficient(
                        m_surf_state.getSurfaceRhoi(),
                        m_surf_state.getSurfaceT());

            m_surf_state.solidProps().getPyrolysingGasEquilMassFrac(
                mv_r_to_pg[i_r], Psurf, Tsurf, mv_work);

            v_chem_rate_per_gas_and_solid.head(m_ns) -= reac_rate_const * mv_work;
            v_chem_rate_per_gas_and_solid(mv_r_to_ps[i_r]) += reac_rate_const;
        }
    }

//=============================================================================

    int nSurfaceReactions(){ return m_nr; }

//=============================================================================
private:
    const size_t m_ns;
    const size_t m_nr;
    std::vector<int> mv_r_to_pg;
    std::vector<int> mv_r_to_ps;


    Eigen::VectorXd mv_work;
    Eigen::VectorXd mv_chem_rate;
    Eigen::VectorXd mv_react_rate_const;

    GSIStoichiometryManager m_reactants;
};

ObjectProvider<GSIRateManagerBulkChemistry, GSIRateManager>
    gsi_rate_manager_bulk_chemistry("bulk_chemistry");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
