/**
 * @file SolidPropertiesBulkChemistry.cpp
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

#include "Composition.h"
#include "SolidProperties.h"
#include "Thermodynamics.h"

using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities::IO;
using namespace Mutation::Utilities::Config;

using namespace Eigen;

namespace Mutation {
    namespace GasSurfaceInteraction {
/**
 *
 */
class SolidPropertiesBulkChemistry : public SolidProperties {
public:
    /**
     * Default constructor
     */
    SolidPropertiesBulkChemistry(ARGS args);

//==============================================================================
    /**
     * Default destructor
     */
    ~SolidPropertiesBulkChemistry(){}

//==============================================================================

    int pyrolysisSpeciesIndex(const std::string& str_sp) const;

//==============================================================================

    int nPyrolysingSolids() const { return m_n_ps; }

//==============================================================================

    void setPyrolysingSolidDensities(
        const VectorXd& v_rho_pyro_solid) const {
        mv_rhops = v_rho_pyro_solid;
        for (int i = 0; i < m_n_ps; i++) {
            if (mv_rhops(i) > mv_rhops_i(i) || mv_rhops(i) < mv_rhops_f(i) ) {
                std::cerr << "ERROR THE SET DENSITIES ARE WRONG." << std::endl;
                exit(1);
            }
        }
    }

//==============================================================================

    double getPyrolysingSolidDensity(const int& sp) const;

//==============================================================================

    double getPyrolysingSolidInitialDensity(const int& sp) const;

//==============================================================================

    double getPyrolysingSolidFinalDensity(const int& sp) const;

//==============================================================================

    int nPyrolysingGases() const { return m_n_pg; }

//==============================================================================

    void getPyrolysingGasEquilMassFrac(
        const int& sp, const double& P, const double& T,
        Eigen::VectorXd& v_yi) const
    {
        int id = sp - m_n_g - m_n_ps;

        m_thermo.equilibriumComposition(
            T, P, mv_pg_comp[id].data(), v_yi.data());

        m_thermo.convert<X_TO_Y>(v_yi.data(), v_yi.data());
    }

//==============================================================================
private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;

    const int m_n_g;

    int m_n_ps;
    std::vector<std::string> mv_ps; // bulk stores chemical symbol
    Eigen::VectorXd mv_rhops_i;
    mutable Eigen::VectorXd mv_rhops;
    Eigen::VectorXd mv_rhops_f;

    int m_n_pg;
    std::vector<std::string> mv_pg; // gas stores chemical symbol
    std::vector<VectorXd> mv_pg_comp;
};

ObjectProvider<
    SolidPropertiesBulkChemistry, SolidProperties>
    solid_properties_bulk_chemistry("bulk_chemistry");

//==============================================================================

SolidPropertiesBulkChemistry::SolidPropertiesBulkChemistry(ARGS args)
    : SolidProperties(args),
      m_thermo(args.s_thermo),
      m_n_g(m_thermo.nSpecies())
{
    assert(args.s_node_solid_props.tag() == "solid_properties");

    // Helper vectors
    std::vector<double> v_dens_i;
    std::vector<double> v_dens_f;

    // Looping over all the children
    for (XmlElement::const_iterator iter_phase =
        args.s_node_solid_props.begin();
        iter_phase != args.s_node_solid_props.end();
        iter_phase++)
    {
        std::string option = iter_phase->tag();

        std::string label;
        std::string symbol;
        double wrk;

        if (option.compare("pyro_solid") == 0) {
            iter_phase->getAttribute(
                "label", label,
                "Error in surface_composition for the "
                "surface properties. A label should be provided.");

            iter_phase->getAttribute(
                "component", symbol,
                "Error in SolidPropertiesBulkChemistry species");

            mv_ps.push_back(symbol + '-' + label);

            iter_phase->getAttribute(
                "init_dens", wrk,
                "Error in solid_composition for the "
                "surface properties. An initial density should be provided.");
            v_dens_i.push_back(wrk);

            iter_phase->getAttribute(
                "final_dens", wrk,
                "Error in solid_composition for the "
                "solid properties. A final density should be provided.");
            v_dens_f.push_back(wrk);

        } else if (option.compare("pyro_gas") == 0) {
            iter_phase->getAttribute(
                 "label", label,
                "Error in solid_composition for the "
                "solid properties. A label should be provided.");

            iter_phase->getAttribute(
                "composition", symbol,
                "Error in SolidProperties, composition");

            symbol += '-' + label;
            mv_pg.push_back(symbol);

            std::string elem;
            iter_phase->getAttribute(
                "elem_comp", elem,
                "Error in SolidPropertiesBulkChemistry, elemental composition");

            std::map<std::string, int> map;
            for(int i= 0; i < args.s_thermo.nElements(); i++)
                map[args.s_thermo.elementName(i)] = i;

            VectorXd v_el_comp(args.s_thermo.nElements());
            Composition m_comp(symbol, elem);
            m_comp.getComposition(map, v_el_comp.data());

            mv_pg_comp.push_back(v_el_comp);

        } else {
            throw InvalidInputError("SolidProperties", option)
                << option << "is a wrong solid properties input "
                << "for bulk chemistry.";
        }
    }
    m_n_ps = mv_ps.size();
    m_n_pg = mv_pg_comp.size();
    // if (m_n_g = 0) error

    mv_rhops_i.resize(m_n_ps);
    mv_rhops_f.resize(m_n_ps);
    for (int i = 0; i < m_n_ps; i++){
        mv_rhops_i(i) = v_dens_i[i];
        mv_rhops_f(i) = v_dens_f[i];
    }

    // Initial Conditions:
    mv_rhops = mv_rhops_i;
}

//==============================================================================

    int SolidPropertiesBulkChemistry::pyrolysisSpeciesIndex(
        const std::string& str_sp) const
    {
        for (int i_ps = 0; i_ps < m_n_ps; i_ps++) {
            if (mv_ps[i_ps] == str_sp)
                return m_n_g + i_ps;
        }
        for (int i_pg = 0; i_pg < m_n_pg; i_pg++) {
            if (mv_pg[i_pg] == str_sp)
                return m_n_g + m_n_ps + i_pg;
        }
        return -1;
    }

//==============================================================================

    double SolidPropertiesBulkChemistry::getPyrolysingSolidDensity(
        const int& sp) const
    {
        int id = sp - m_n_g;
        if (id >= 0 && id < m_n_ps)
            return mv_rhops(id);

        return 0.;
    }

//==============================================================================

    double SolidPropertiesBulkChemistry::getPyrolysingSolidInitialDensity(
        const int& sp) const
    {
        int id = sp - m_n_g;
        if (id >= 0 && id < m_n_ps)
            return mv_rhops_i(id);

        return 0.;
    }

//==============================================================================

    double SolidPropertiesBulkChemistry::getPyrolysingSolidFinalDensity(
        const int& sp) const
    {
        int id = sp - m_n_g;
        if (id >= 0 && id < m_n_ps)
            return mv_rhops_f(id);

        return 0.;
    }


    } // namespace GasSurfaceInteraction
} // namespace Mutation
