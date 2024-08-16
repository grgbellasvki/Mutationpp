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
using namespace std;

using namespace Mutation::Numerics;
using namespace Mutation::Utilities::Config;
using namespace Mutation::Thermodynamics;

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
          mv_kf(m_nr),
          mv_rate(m_nr),
          mn_site_sp(m_surf_props.nSiteSpecies()),
          mn_site_cat(m_surf_props.nSiteCategories()),
          mv_sp_in_site(mn_site_cat), // Define
          mv_sigma(mn_site_cat), // Define
          mv_rhoi(m_ns),
          mv_nd(m_ns+mn_site_sp),
          is_surf_steady_state(true),
          m_tol(1.e-15),
          m_pert(1.e-1),
          mv_X(mn_site_sp),
          mv_dX(mn_site_sp),
          mv_f_unpert(m_ns+mn_site_sp),
          mv_f(m_ns+mn_site_sp),
          m_jac(mn_site_sp,mn_site_sp),
          mv_mass(m_thermo.speciesMw()/NA)
    {
        for (int i_reac = 0; i_reac < m_nr; ++i_reac) {
            m_reactants.addReaction(
                i_reac, args.s_reactions[i_reac]->getReactants());
            m_irr_products.addReaction(
                i_reac, args.s_reactions[i_reac]->getProducts());
        }

        for (int i = 0; i < mn_site_cat; i++) { // All these better be done in SurfaceProperties
            mv_sigma(i) = m_surf_props.nSiteDensityInCategory(i);        // Remove this function
            mv_sp_in_site[i] = m_surf_props.nSpeciesInSiteCategory(i);   // Remove this function
        }

        // Setup NewtonSolver
        setMaxIterations(50);
        setWriteConvergenceHistory(false);
        setEpsilon(m_tol);
    }

//=============================================================================

    ~GSIRateManagerDetailed(){}

//=============================================================================

    Eigen::VectorXd computeRates()
    {
        // Note mv_nd is density he
        mv_rhoi = m_surf_state.getSurfaceRhoi();

        // Get reaction rate constant
        for (int i_r = 0; i_r < m_nr; ++i_r)
            mv_kf(i_r) =
                v_reactions[i_r]->getRateLaw()->forwardReactionRateCoefficient(
                    mv_rhoi, m_surf_state.getSurfaceT());

        // Getting all number densities
        m_thermo.convert<RHO_TO_CONC>(mv_rhoi.data(), mv_nd.data());
        mv_nd.head(m_ns) *= NA;
        mv_nd.tail(mn_site_sp) = m_surf_props.getSurfaceSiteCoverageFrac();

        int k = 0;
        for (int i = 0; i < mn_site_cat; i++)
            for (int j = 0; j < mv_sp_in_site[i]; j++, k++)
                mv_nd(m_ns+k) *= mv_sigma(i);  // All these better be done in SurfaceProperties

        // In the case of surface at steady state is considered,
        // it updates the mv_nd.tail(mn_site_sp)
        bool is_surf_cov_steady_state = m_surf_props.isSurfaceCoverageSteady();
        if (is_surf_cov_steady_state)
            computeSurfaceSteadyStateCoverage();

        // Constant rate times densities of species
        mv_rate = mv_kf;
        m_reactants.multReactions(mv_nd, mv_rate);

        mv_nd.setZero();
        m_reactants.incrSpecies(mv_rate, mv_nd);
        m_irr_products.decrSpecies(mv_rate, mv_nd);

        // Multiply by molar mass
        return (mv_nd.cwiseProduct(mv_mass)).head(m_ns);
    }

//=============================================================================

    Eigen::VectorXd computeRatesPerReaction()
    {
        // Note mv_nd is density he
        mv_rhoi = m_surf_state.getSurfaceRhoi();

    	// Getting the kfs with the initial conditions
        for (int i_r = 0; i_r < m_nr; ++i_r) {
            mv_kf(i_r) =
                v_reactions[i_r]->getRateLaw()->forwardReactionRateCoefficient(
            		mv_rhoi, m_surf_state.getSurfaceT());
        }

        m_thermo.convert<RHO_TO_CONC>(mv_rhoi.data(), mv_nd.data());
        mv_nd.head(m_ns) *= NA;
        mv_nd.tail(mn_site_sp) = m_surf_props.getSurfaceSiteCoverageFrac();

        int k = 0;
        for (int i = 0; i < mn_site_cat; i++)
            for (int j = 0; j < mv_sp_in_site[i]; j++, k++)
                mv_nd(m_ns+k) *= mv_sigma(i);  // All these better be done in SurfaceProperties

        // In the case of surface at steady state,
        // it updates the mv_nd.tail(mn_site_sp)
        bool is_surf_cov_steady_state = m_surf_props.isSurfaceCoverageSteady();
        if (is_surf_cov_steady_state)
            computeSurfaceSteadyStateCoverage();

        mv_rate = mv_kf;
        m_reactants.multReactions(mv_nd, mv_rate);


        return mv_rate;
    }


//=============================================================================

    int nSurfaceReactions(){ return m_nr; }

//=============================================================================

    void updateFunction(VectorXd& v_X)
    {
        // std::cout << "v_X = " << v_X << std::endl;
        // mv_f.head(m_ns) = mv_nd.head(m_ns);
        mv_nd.tail(mn_site_sp) = v_X;

        // std::cout << "nd = " << mv_nd << std::endl;
        // std::cout << "kf[0] = " << mv_kf[0] << std::endl;
        // std::cout << "kf[1] = " << mv_kf[1] << std::endl;

        mv_rate = mv_kf;
        m_reactants.multReactions(mv_nd, mv_rate);
        // std::cout << "Rate[0] = " << mv_rate[0] << std::endl;
        // std::cout << "Rate[1] = " << mv_rate[1] << std::endl;

        mv_f.setZero();
    	m_reactants.incrSpecies(mv_rate, mv_f);
    	m_irr_products.decrSpecies(mv_rate, mv_f);
        // std::cout << "F = \n" << mv_f << std::endl;
    }
//=============================================================================

    void updateJacobian(VectorXd& v_X)
    {
    	mv_f_unpert = mv_f;
        // Make pert
        for (int i = 0; i < mn_site_sp; ++i) {
            double m_X_unpert = v_X(i);
            double pert = v_X(i) * m_pert;
    	    v_X(i) += pert;

            // Update Jacobian column
            updateFunction(v_X);

            m_jac.col(i) =
                (mv_f.tail(mn_site_sp)-mv_f_unpert.tail(mn_site_sp)) / pert;
            v_X(i) = m_X_unpert;
        }
    }

//=============================================================================

    VectorXd& systemSolution()
    {
        // std::cout << "X = \n" << mv_X << std::endl;

        // std::cout << "jac = \n" << m_jac << std::endl;

        double a = (m_jac.diagonal()).cwiseAbs().maxCoeff();
        // std::cout << "jac = \n" << m_jac + a*MatrixXd::Ones(mn_site_sp,mn_site_sp) << std::endl;
        // std::cout << "f = \n" << mv_f_unpert.tail(mn_site_sp) << std::endl;

        mv_dX = (m_jac + a*MatrixXd::Ones(mn_site_sp,mn_site_sp)).
            fullPivLu().solve(mv_f_unpert.tail(mn_site_sp));
        // m_jac(0,0) = 1;
        // m_jac(0,1) = 1;
        // mv_f_unpert(0) = mv_sigma(0);
        // mv_dX = (m_jac).
        //     fullPivLu().solve(mv_f_unpert.tail(mn_site_sp));

        // std::cout << "dx = \n" << mv_dX << std::endl;
        // double in; std::cin >> in;
        return mv_dX;
    }
//=============================================================================

    double norm()
    {
        //return mv_f_unpert.tail(mn_site_sp).lpNorm<Eigen::Infinity>();
        return mv_dX.lpNorm<Eigen::Infinity>(); 
    }
//=============================================================================
private:
    void computeSurfaceSteadyStateCoverage(){

        mv_X.setZero();
        for (int i = 0; i < mn_site_cat; ++i)
            mv_X.setConstant(mv_sigma(i) / mv_sp_in_site[i]);

        mv_X = solve(mv_X);

        applyTolerance(mv_X);

        // Setting up the SurfaceSiteCoverage.
        // This is not essential for efficiency.
        mv_nd.tail(mn_site_sp) = mv_X;
        m_surf_props.setSurfaceSiteCoverageFrac(mv_X/mv_sigma(0));
    }

//=============================================================================
    inline void applyTolerance(Eigen::VectorXd& v_x) const {
        for (int i = 0; i < v_x.size(); i++)
            if (std::abs(v_x(i) / mv_sigma(0)) < 1e-30) v_x(i) = 0.;
            //if (std::abs(v_x(i) / mv_sigma(0)) < m_tol) v_x(i) = 0.; //prevent system to converge
    }

//=============================================================================
private:
    SurfaceProperties& m_surf_props;

    const size_t m_ns;
    const size_t m_nr;

    bool is_surf_steady_state;

    VectorXd mv_kf;
    VectorXd mv_rate;

    // Surface properties
    const size_t mn_site_sp;
    const size_t mn_site_cat;
    vector<int> mv_sp_in_site;
    VectorXd mv_sigma;

    VectorXd mv_mass;
    VectorXd mv_rhoi;
    VectorXd mv_nd;

    // For the steady state coverage solver.
    const double m_tol;
    double m_pert;
    VectorXd mv_X;
    VectorXd mv_dX;
    VectorXd mv_work;
    VectorXd mv_f;
    VectorXd mv_f_unpert;
    Eigen::MatrixXd m_jac;

    GSIStoichiometryManager m_reactants;
    GSIStoichiometryManager m_irr_products;
};

ObjectProvider<GSIRateManagerDetailed, GSIRateManager>
    gsi_rate_manager_detailed_mass("detailed_mass");

ObjectProvider<GSIRateManagerDetailed, GSIRateManager>
    gsi_rate_manager_detailed_mass_energy("detailed_mass_energy");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
