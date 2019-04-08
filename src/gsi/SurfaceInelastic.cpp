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

double SurfaceInelastic::surfaceInelasticTerm(const VectorXd& v_X, const VectorXd& v_h, const VectorXd& chem_souce)
{

    //computing Vibrational traslational exchange
    double T_tra = v_X(pos_E);
    double T_vib = v_X(pos_E + 1);
    Eigen::VectorXd h_tra(m_ns);
    Eigen::VectorXd h_vib(m_ns);
    Eigen::VectorXd h_el(m_ns);
    m_thermo.speciesHOverRT(T_tra, T_tra, T_tra, T_tra, T_tra, NULL, h_tra.data(), NULL, h_vib.data(), h_el.data(), NULL);	
    double number_density = m_thermo.numberDensity();
    double inleastic_term = 0.;
    double tol = 1.E-19;
    if (m_eff_coll > tol) {
        for(int i = 0; i < m_ns; ++i) {
	    //thermal speed of species i
	    double thermal_speed;
	    //if (i < m_thermo.hasElectrons()) thermal_speed = sqrt(RU*T_vib/(2.*PI*m_speciesMw(i)));
	    //else thermal_speed = sqrt(RU*T_tra/(2.*PI*m_speciesMw(i))); 
	    if (i < m_thermo.hasElectrons()) thermal_speed = sqrt(T_vib)*m_therm_vel_over_T(i);
	    else thermal_speed = sqrt(T_tra)*m_therm_vel_over_T(i); 

	    //number denisty species
	    double num_dens_i = v_X(i)*number_density;

	    //characteristc time ^-1
	    double one_over_tau = num_dens_i*thermal_speed;
	    //std::cout << "one_over_tau is" << one_over_tau << std::endl;

	    //vibronic enthaply at Teq
	    double h_VE;
	    if (i < m_thermo.hasElectrons()) h_VE = (h_tra[i]*T_tra)*RU/m_speciesMw(i);  //check it
	    else h_VE = (h_vib(i) + h_el(i))*T_tra*RU/m_speciesMw(i);

	    inleastic_term += m_eff_coll*(h_VE-v_h(m_ns + i))*one_over_tau;	    	
	    }
    }

    //Compute vibrational chemical production  
    if (1.-m_beta > tol) inleastic_term += (1-m_beta)*v_h(m_index)*chem_souce(m_index);

    return inleastic_term;
}
    } // namespace GasSurfaceInteraction
} // namespace Mutation
