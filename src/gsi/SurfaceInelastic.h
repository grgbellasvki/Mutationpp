/**
 * @file GasFourierHeatFluxCalculator.h
 *
 * @brief Declaration of GasFourierHeatFluxCalculator class.
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



#ifndef Surface_Inelastic
#define Surface_Inelastic

#include <eigen3/Eigen/Dense>

namespace Mutation { namespace Thermodynamics { class Thermodynamics; }}
namespace Mutation { namespace Utilities { namespace IO { class XmlElement; }}}

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 * Class responsible for computing the no-equilibrium term in the GSI.
 */

class SurfaceInelastic
{
public:
    /**
     *  Constructor
     */
    SurfaceInelastic(
    Mutation::Thermodynamics::Thermodynamics& thermo,
    const Mutation::Utilities::IO::XmlElement& xml_surf_inelastic);

//==============================================================================
    /**
     *  Destructor
     */
    ~SurfaceInelastic();

//==============================================================================
    /**
     *  Function which returns the inelastic term.
     */
    double surfaceInelasticTerm(const Eigen::VectorXd& v_X, const Eigen::VectorXd& v_h,  const Eigen::VectorXd& chem_souce);

//==============================================================================

private:
    Mutation::Thermodynamics::Thermodynamics& m_thermo;

    const int pos_E;
    double m_eff_coll;
    double m_beta;
    const double m_ns;
    const double m_index;

    Eigen::VectorXd m_speciesMw;
    Eigen::VectorXd m_therm_vel_over_T;

};

    } // namespace GasSurfaceInteraction
} // namespace Mutation

#endif // SURFACE_INELASTIC_H
