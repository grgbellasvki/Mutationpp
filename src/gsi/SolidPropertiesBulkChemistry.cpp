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
#include "SurfaceState.h"
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
            if (mv_rhops_i(i) > mv_rhops_f(i) &&
                (mv_rhops(i) > mv_rhops_i(i) || mv_rhops(i) < mv_rhops_f(i))) {
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

    void solidEffectiveThermalConductivity(VectorXd v_solid_lambda) const;

//==============================================================================

    void solidHeatCapacity(VectorXd v_solid_cp) const;


//==============================================================================
private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const SurfaceState& m_surf_state;

    const int m_n_g;
    const size_t pos_T_trans;
    const double tol;

    int m_n_ps;
    std::vector<std::string> mv_ps; // bulk stores chemical symbol
    Eigen::VectorXd mv_rhops_i;
    mutable Eigen::VectorXd mv_rhops;
    Eigen::VectorXd mv_rhops_f;

    int m_n_pg;
    std::vector<std::string> mv_pg; // gas stores chemical symbol
    std::vector<VectorXd> mv_pg_comp;

    const size_t m_n_coeff;
    mutable Eigen::VectorXd mv_temp;

    bool is_isotropic;
    Eigen::VectorXd v_cond_dx;
    Eigen::VectorXd v_cond_dy;
    Eigen::VectorXd v_cond_dz;
    Eigen::VectorXd v_cond_fx;
    Eigen::VectorXd v_cond_fy;
    Eigen::VectorXd v_cond_fz;

    Eigen::VectorXd v_cp_d;
    Eigen::VectorXd v_cp_f;
};

ObjectProvider<
    SolidPropertiesBulkChemistry, SolidProperties>
    solid_properties_bulk_chemistry("bulk_chemistry");

//==============================================================================

SolidPropertiesBulkChemistry::SolidPropertiesBulkChemistry(ARGS args)
    : SolidProperties(args),
      m_thermo(args.s_thermo),
      m_surf_state(args.s_surf_state),
      m_n_g(m_thermo.nSpecies()),
      m_n_coeff(7),
      mv_temp(m_n_coeff),
      tol(1.e-6),
      pos_T_trans(0)
{
    assert(args.s_node_solid_props.tag() == "solid_properties");

    // Helper vectors
    std::vector<double> v_dens_i;
    std::vector<double> v_dens_f;

    // Looping over all the children
    for (XmlElement::const_iterator iter =
        args.s_node_solid_props.begin();
        iter != args.s_node_solid_props.end();
        iter++)
    {
        std::string option = iter->tag();

        std::string label;
        std::string symbol;
        double wrk;

        if (option.compare("effective_cond") == 0) {
            iter->getAttribute("isotropic", is_isotropic,
            "Isotropic property should be provided.");

            for (XmlElement::const_iterator mat = iter->begin();
                mat != iter->end();
                mat++)
            {
                std::string material = mat->tag();

                if (material.compare("decomposing") == 0) {

                    XmlElement::const_iterator dir = mat->findTag("x");
                    std::istringstream ssx(dir->text());

                    std::copy(
                        std::istream_iterator<double>(ssx),
                        std::istream_iterator<double>(),
                        v_cond_dx.data());

                    if (is_isotropic){
                        v_cond_dy = v_cond_dx;
                        v_cond_dz = v_cond_dx;
                    } else {
                        dir = iter->findTag("y");
                        std::istringstream ssy(dir->text());

                        std::copy(
                            std::istream_iterator<double>(ssy),
                            std::istream_iterator<double>(),
                            v_cond_dy.data());

                        dir = iter->findTag("z");
                        std::istringstream ssz(dir->text());

                        std::copy(
                            std::istream_iterator<double>(ssz),
                            std::istream_iterator<double>(),
                            v_cond_dz.data());
                    }

                } else if (material.compare("final") == 0) {
                    XmlElement::const_iterator dir = mat->findTag("x");
                    std::istringstream ssx(dir->text());

                    std::copy(
                        std::istream_iterator<double>(ssx),
                        std::istream_iterator<double>(),
                        v_cond_fx.data());

                    if (is_isotropic){
                        v_cond_fy = v_cond_fx;
                        v_cond_fz = v_cond_fx;
                    } else {
                        dir = iter->findTag("y");
                        std::istringstream ssy(dir->text());

                        std::copy(
                            std::istream_iterator<double>(ssy),
                            std::istream_iterator<double>(),
                            v_cond_fy.data());

                        dir = iter->findTag("z");
                        std::istringstream ssz(dir->text());

                        std::copy(
                            std::istream_iterator<double>(ssz),
                            std::istream_iterator<double>(),
                            v_cond_fz.data());
                    }
                } else {
                    // "ERROR"
                }

                // Assert = max_coeff
            }
        } else if (option.compare("heat_capacity") == 0) {
            for (XmlElement::const_iterator mat = iter->begin();
                mat != iter->end();
                mat++)
            {
                std::string material = mat->tag();

                if (material.compare("decomposing") == 0) {
                    std::istringstream ss(mat->text());

                    std::copy(
                        std::istream_iterator<double>(ss),
                        std::istream_iterator<double>(),
                        v_cp_d.data());
                } else if (material.compare("final") == 0) {
                    std::istringstream ss(mat->text());

                    std::copy(
                        std::istream_iterator<double>(ss),
                        std::istream_iterator<double>(),
                        v_cp_f.data());
                } else {
                    // "ERROR"
                }

                // Assert max_coeff
            }
        } else if (option.compare("pyro_solid") == 0) {
            iter->getAttribute(
                "label", label,
                "Error in surface_composition for the "
                "surface properties. A label should be provided.");

            iter->getAttribute(
                "component", symbol,
                "Error in SolidPropertiesBulkChemistry species");

            mv_ps.push_back(symbol + '-' + label);

            iter->getAttribute(
                "init_dens", wrk,
                "Error in solid_composition for the "
                "surface properties. An initial density should be provided.");
            v_dens_i.push_back(wrk);

            iter->getAttribute(
                "final_dens", wrk,
                "Error in solid_composition for the "
                "solid properties. A final density should be provided.");
            v_dens_f.push_back(wrk);

        } else if (option.compare("pyro_gas") == 0) {
            iter->getAttribute(
                 "label", label,
                "Error in solid_composition for the "
                "solid properties. A label should be provided.");

            iter->getAttribute(
                "composition", symbol,
                "Error in SolidProperties, composition");

            symbol += '-' + label;
            mv_pg.push_back(symbol);

            std::string elem;
            iter->getAttribute(
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

//==============================================================================

    void SolidPropertiesBulkChemistry::solidEffectiveThermalConductivity(
        VectorXd v_solid_lambda) const
    {
        double Tsolid1 = m_surf_state.getSurfaceT()(pos_T_trans);
        double Tsolid2 = Tsolid1  * Tsolid1;
        double Tsolid3 = Tsolid2 * Tsolid1;
        double Tsolid4 = Tsolid3 * Tsolid1;
        double Tsolidm1 = 1.       / Tsolid1;
        double Tsolidm2 = Tsolidm1 / Tsolid1;

        mv_temp(0) = Tsolid4;
        mv_temp(1) = Tsolid3;
        mv_temp(2) = Tsolid2;
        mv_temp(3) = Tsolid1;
        mv_temp(4) = 1.;
        mv_temp(5) = Tsolidm1;
        mv_temp(6) = Tsolidm2;

        if ( ((mv_rhops_f-mv_rhops).cwiseAbs()).maxCoeff() < tol ) {
                v_solid_lambda(0) = mv_temp.dot(v_cond_fx);
                v_solid_lambda(1) = mv_temp.dot(v_cond_fy);
                v_solid_lambda(2) = mv_temp.dot(v_cond_fz);
        } else {
                v_solid_lambda(0) = mv_temp.dot(v_cond_dx);
                v_solid_lambda(1) = mv_temp.dot(v_cond_dy);
                v_solid_lambda(2) = mv_temp.dot(v_cond_dz);
        }
    }

//==============================================================================

    void SolidPropertiesBulkChemistry::solidHeatCapacity(
        VectorXd v_solid_cp) const
        {

            //v_solid_cp = 
        }

//==============================================================================

    } // namespace GasSurfaceInteraction
} // namespace Mutation
