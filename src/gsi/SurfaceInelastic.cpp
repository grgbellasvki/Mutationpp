/**
 * @file GasFourierHeatFluxCalculator.cpp
 *
 * @brief Class which computes the gas heat flux needed by
 *        the surface energy balances.
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
#include "Utilities.h"

#include "SurfaceInelastic.h"
#include "SurfaceState.h"

using namespace Eigen;

using namespace Mutation;
using namespace Mutation::Utilities;

namespace Mutation {
    namespace GasSurfaceInteraction {

SurfaceInelastic::SurfaceInelastic(
    Mutation::Thermodynamics::Thermodynamics& thermo,
    const Mutation::Utilities::IO::XmlElement& xml_surf_inelastic)
        :  m_thermo(thermo),
           pos_E(thermo.nSpecies()),
	   m_ns(thermo.nSpecies()),
	   m_speciesMw(thermo.speciesMw()),
	   m_index(m_thermo.speciesIndex("CN")),
	   m_therm_vel_over_T(sqrt(RU/(2.*PI*(thermo.speciesMw()))))
{
    xml_surf_inelastic.getAttribute("effective_collisions", m_eff_coll, 1.);
    xml_surf_inelastic.getAttribute("accomodation_coef", m_beta, 1.);
}

//==============================================================================

SurfaceInelastic::~SurfaceInelastic(){}

//==============================================================================

double SurfaceInelastic::surfaceInelasticTerm(const VectorXd& v_X, const VectorXd& v_h, const VectorXd& chem_souce, const VectorXd& v_rhoi)
{

    //computing Vibrational traslational exchange
    double T_tra = v_X(pos_E);
    double T_vib = v_X(pos_E + 1);
    double inleastic_term = 0.;
    double thermal_speed;
    double num_dens_i;
    double one_over_tau;
    double h_VE;
    double h_VV_per_particle;

    Eigen::VectorXd h_tra(m_ns);
    Eigen::VectorXd h_vib(m_ns);
    Eigen::VectorXd h_el(m_ns);

    const int set_state_with_rhoi_T = 1;
    m_thermo.setState(
        v_rhoi.data(), v_X.tail(2).data(), set_state_with_rhoi_T);
    double number_density = m_thermo.numberDensity();

    m_thermo.speciesHOverRT(T_tra, T_tra, T_tra, T_tra, T_tra, NULL, h_tra.data(), NULL, h_vib.data(), h_el.data(), NULL);
    for(int i = 0; i < m_ns; ++i) {

	//thermal speed of species i [m/s]
	if (i < m_thermo.hasElectrons()) thermal_speed = sqrt(T_vib)*m_therm_vel_over_T(i);
	else thermal_speed = sqrt(T_tra)*m_therm_vel_over_T(i);

        //number denisty species
	num_dens_i = v_X(i)*number_density;

	//impinging particle flux [# / m^2 s ]
	one_over_tau = num_dens_i*thermal_speed;

	//vibronic enthaply at Teq of the single particles [J/ #]
	if (i < m_thermo.hasElectrons()) h_VE = h_tra[i]*T_tra*RU / NA;
	else h_VE = (h_vib(i) + h_el(i))*T_tra*RU / NA;

        //from mass entalphy to element enthalphy [J/ #]
        h_VV_per_particle =  v_h(m_ns + i) * m_speciesMw(i) / NA;

	// [  J / m^2 s ]
	inleastic_term += m_eff_coll*(h_VE - h_VV_per_particle)*one_over_tau;
	}

    //Compute vibrational chemical production  
    inleastic_term -= (1.-m_beta)*v_h(m_index)*chem_souce(m_index);

    return inleastic_term;
}
    } // namespace GasSurfaceInteraction
} // namespace Mutation
