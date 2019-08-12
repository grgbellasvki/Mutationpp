/**
 * @file GSIRateManagerPhenomenological.cpp
 *
 * @brief Class which computes the chemical production rate for each species
 *        based on detailed chemistry models for catalysis and ablation.
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


#include "NewtonSolver.h"
#include "Thermodynamics.h"
#include "Transport.h"

#include "GSIReaction.h"
#include "GSIRateLaw.h"
#include "GSIRateManager.h"
#include "GSIStoichiometryManager.h"
#include "SurfaceProperties.h"
#include "SurfaceState.h"

using namespace Eigen;

using namespace Mutation::Numerics;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

//=============================================================================

class GSIRateManagerDetailed :
    public GSIRateManager,
    public NewtonSolver<VectorXd, GSIRateManagerDetailed>
{
public:
    GSIRateManagerDetailed(DataGSIRateManager args)
        : GSIRateManager(args),
          m_surf_props(args.s_surf_state.getSurfaceProperties()),
		  m_ns(args.s_thermo.nSpecies()),
		  m_nr(args.s_reactions.size()),
          mv_react_rate_const(m_nr),
		  mv_work(m_ns),
          is_surf_steady_state(true)
    {
        for (int i_reac = 0; i_reac < m_nr; ++i_reac) {
            m_reactants.addReaction(
                i_reac, args.s_reactions[i_reac]->getReactants());
            m_irr_products.addReaction(
                i_reac, args.s_reactions[i_reac]->getProducts());
        }

        // Setup NewtonSolver
        setMaxIterations(5);
        setWriteConvergenceHistory(false);
        setEpsilon(1.e-13);
    }

//=============================================================================

    ~GSIRateManagerDetailed(){}

//=============================================================================

    Eigen::VectorXd computeRates()
    {
        // Get reaction rate constant
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_react_rate_const(i_r) =
                v_reactions[i_r]->getRateLaw()->
                    forwardReactionRateCoefficient(
                        m_surf_state.getSurfaceRhoi(),
                        m_surf_state.getSurfaceT());
        }

        if (is_surf_steady_state)
            computeSurfaceSteadyStateCoverage();

        // Getting all number densities
        // mv_nd.head(m_ns) = 0.;
        // mv_nd.tail(m_site_sp) = 0.;

        // Constant rate times densities of species
        mv_work.setZero();
        m_reactants.incrSpecies(mv_react_rate_const, mv_work);
        m_irr_products.decrSpecies(mv_react_rate_const, mv_work);

        // Multiply by molar mass
        return mv_work.cwiseProduct(m_thermo.speciesMw().matrix());

        /* lv_numb_dens.setZero();

    	// Getting the initial number densities
    	m_wall_state.getNdStateGasSurf(lv_numb_dens); // @todo return the value

    	// Getting the kfs with the initial conditions
        for (int i_reac = 0; i_reac < m_nr; ++i_reac) {
            v_kf(i_reac) =  v_reactions[i_reac]->
            		            getRateLaw()->forwardReactionRateCoefficient(
            		                m_wall_state.getWallRhoi(), m_wall_state.getWallT());
        }

        // Solving for the steady state for the Surface Species!
        computeSurfSpeciesBalance();

        v_RateofProduction = v_kf;
        m_reactants.multReactions(lv_numb_dens, v_RateofProduction);

        lv_numb_dens.setZero();
        m_reactants.decrSpecies(v_RateofProduction, lv_numb_dens);
        m_irr_products.incrSpecies(v_RateofProduction, lv_numb_dens);

        // Multiply by molar mass
        return (     -     lv_numb_dens.cwiseProduct(     // HERE IS THE FAMOUS "-" THAT I ONLY USE TO CONFUSE MYSELF AND A. TURCHI
        		Eigen::Map<const Eigen::VectorXd>(m_thermo.speciesMw(), m_ns_gas))/Mutation::NA).head(m_ns_gas); */
    }

//=============================================================================

    Eigen::VectorXd computeRatesPerReaction()
    {
        // IMPLIMENT ME!
    	// Getting the kfs with the initial conditions
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_react_rate_const(i_r) =
                v_reactions[i_r]->getRateLaw()->forwardReactionRateCoefficient(
            		m_surf_state.getSurfaceRhoi(), m_surf_state.getSurfaceT());
        }
        m_reactants.multReactions(
            m_surf_state.getSurfaceRhoi(), mv_react_rate_const);

        return mv_react_rate_const;
    }


//=============================================================================

    int nSurfaceReactions(){ return m_nr; }

//=============================================================================
private:
    void computeSurfaceSteadyStateCoverage(){

        // Initial coverage
        mv_X = m_surf_props.getSurfaceSiteCoverageFrac(); // Surface Coverage
        // mv_X = solve(mv_X);

        // Setting up the SurfaceSiteCoverage.
        // This is not essential for efficiency.
        m_surf_props.setSurfaceSiteCoverageFrac(mv_X);
    }
//=============================================================================

    void updateFunction(VectorXd& v_X){}
//=============================================================================

    void updateJacobian(VectorXd& v_X){}
//=============================================================================

    VectorXd& systemSolution(){}
//=============================================================================

    double norm(){return 0.;}
//=============================================================================
private:
    SurfaceProperties& m_surf_props;

    const size_t m_ns;
    const size_t m_nr;

    bool is_surf_steady_state;

    VectorXd mv_react_rate_const;
    VectorXd mv_work;

    // For the steady state coverage solver.
    VectorXd mv_X;
    VectorXd mv_dX;

    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;
};

ObjectProvider<GSIRateManagerDetailed, GSIRateManager>
    gsi_rate_manager_detailed_mass("detailed_mass");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
