/**
 * @file SurfaceBulkChemistry.cpp
 *
 * @brief .
 */

/*
 * Copyright 2014-2018 von Karman Institute for Fluid Dynamics (VKI)
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


#include "Errors.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "MassBlowingRate.h"
#include "SolidProperties.h"
#include "Surface.h"
#include "SurfaceChemistry.h"
#include "SurfaceState.h"

using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace GasSurfaceInteraction {

class SurfaceBulkChemistry : public Surface
{
public:
    SurfaceBulkChemistry(ARGS args)
        : m_thermo(args.s_thermo),
          m_surf_state(args.s_surf_state),
          mp_surf_chem(NULL),
          mp_mass_blowing_rate(NULL)
    {
        // Initializing surface chemistry
        if (args.xml_surf_chem.tag() == "surface_chemistry"){
            mp_surf_chem = new SurfaceChemistry(
                m_thermo,
                args.s_transport,
                args.s_gsi_mechanism,
                args.xml_surf_chem,
                m_surf_state);
        }

        // MassBlowingRate
        DataMassBlowingRate data_mass_blowing_rate = {m_thermo, *mp_surf_chem};
        const std::string s_mass_blowing = "isOn";
        mp_mass_blowing_rate = Factory<MassBlowingRate>::create(
            s_mass_blowing, data_mass_blowing_rate);
    }

//=============================================================================

    ~SurfaceBulkChemistry()
    {
        if (mp_surf_chem != NULL) { delete mp_surf_chem; }
        if (mp_mass_blowing_rate != NULL) { delete mp_mass_blowing_rate; }
    }

//=============================================================================

    void computeSurfaceReactionRates(Eigen::VectorXd& v_surf_reac_rates)
    {
        errorSurfaceStateNotSet();

        v_surf_reac_rates.setZero();
        if (mp_surf_chem != NULL)
            mp_surf_chem->surfaceReactionRates(v_surf_reac_rates);
    }

//=============================================================================

    Eigen::VectorXd computeSurfaceReactionRatesPerReaction()
    {
        const int nr = nSurfaceReactions();
        Eigen::VectorXd v_wrk(nr);

        if (mp_surf_chem != NULL){
            mp_surf_chem->surfaceReactionRatesPerReaction(v_wrk);
            return v_wrk;
        }
        throw LogicError()
            << "computeGSIReactionRatePerReaction cannot be invoked "
            << "without defining a surface_chemistry option in "
            << "Gas-Surface Interaction input file.";

        return v_wrk.setZero();
    }

//=============================================================================

    int nSurfaceReactions()
    {
        if (mp_surf_chem != NULL)
            return mp_surf_chem->nSurfaceReactions();

        return 0;
    }

//=============================================================================

    int nPyrolysingSolids() const {
        return m_surf_state.solidProps().nPyrolysingSolids();
    }

//=============================================================================

    void surfaceReactionRatesGasAndSolid(
        Eigen::VectorXd& v_surf_reac_rates_gas_solid)
    {
        errorSurfaceStateNotSet();

        if (mp_surf_chem != NULL)
            mp_surf_chem->surfaceReactionRatesGasAndSolid(
                v_surf_reac_rates_gas_solid);
    }

//=============================================================================

    void setPyrolysingSolidDensities(const Eigen::VectorXd& v_rho_pyro_solid) {
        m_surf_state.solidProps().setPyrolysingSolidDensities(
            v_rho_pyro_solid);
    }

//=============================================================================

    double massBlowingRate()
    {
        return mp_mass_blowing_rate->computeBlowingFlux();
    }

//==============================================================================
private:
    void errorSurfaceStateNotSet() const
    {
        if (!m_surf_state.isSurfaceStateSet()) {
            throw LogicError()
                << "The surface state must have been set!";
        }
    }

//==============================================================================
private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    SurfaceChemistry* mp_surf_chem;
    MassBlowingRate* mp_mass_blowing_rate;
    SurfaceState& m_surf_state;
};

ObjectProvider<
    SurfaceBulkChemistry, Surface>
    surface_bulk_chemistry("bulk_chemistry");

    } // namespace GasSurfaceInteraction
} // namespace Mutation
