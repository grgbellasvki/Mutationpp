/**
 * @file SurfacePropertiesAblation.cpp
 *
 * @brief SurfaceProperties class when an ablative surface is modeled
 *        based on a gamma model.
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

#include <iterator>

#include "AutoRegistration.h"
#include "Thermodynamics.h"
#include "Utilities.h"

#include "SurfaceProperties.h"

using namespace std;

using namespace Mutation::Utilities::Config;
using namespace Mutation::Utilities::IO;

namespace Mutation {
    namespace GasSurfaceInteraction {

/**
 *
 */
class SurfacePropertiesDetailed : public SurfaceProperties
{
public:
	SurfacePropertiesDetailed(ARGS args)
        : SurfaceProperties(args),
          xml_surf_props(args.s_node_surf_props),
          m_thermo(args.s_thermo),
          n_gas_sp(m_thermo.nSpecies()),
          is_surface_set(false)
    {
        assert(xml_surf_props.tag() == "surface_properties");

        // For all bulk compositions
        for (XmlElement::const_iterator iter_phase =
                xml_surf_props.begin();
            iter_phase != xml_surf_props.end();
            iter_phase++)
        {
            string option = iter_phase->tag();

            if (option.compare("surface") == 0) {
                std::string label;
                iter_phase->getAttribute(
                    "label", label,
                    "Error in surface option for the "
                    "surface properties. A label should be provided for "
                    "this type of surface.");

                std::string species;
                iter_phase->getAttribute(
                    "species", species,
                    "Error in surface option for the "
                    "surface properties. Species should be provided for "
                    "this type of surface.");
                parseAblationSpecies(species, label);

                is_surface_set = true;
            } else if (option.compare("site") == 0) {
                std::string label;
                iter_phase->getAttribute("label", label,
                    "Error in sites option for the surface properties. "
                    "A label should be provided.");

                double sigma; // Number density of sites per m2
                iter_phase->getAttribute("sigma", sigma,
                    "Error in sites option for the surface properties. "
                    "The number of sites (sigma) should be provided.");
                v_sigma.push_back(sigma);

                // Inserting empty site
                v_site_sp.push_back(label);
                v_site_sp_to_gas_idx.push_back(-2);

                // Parsing species in site
                std::string species;
                iter_phase->getAttribute(
                    "species", species,
                    "Error in surface option for the "
                    "surface properties. Species should be provided for "
                    "this type of surface.");

                parseSiteSpecies(species, label);
            } else {
                throw InvalidInputError("SurfaceProperties", option)
                 << option << " is a wrong input for surface "
                << "properties.";
            }

            n_site_categ = v_sp_in_site.size();

            n_surf_comp_sp = v_surf_sp.size();
            n_site_sp = v_site_sp.size();
            n_tot_surf_sp = n_surf_comp_sp + n_site_sp;
        }

        for (int i = 0; i < n_site_categ; i++)
            for (int j = 0; j < v_sp_in_site[i]; j++)
                v_site_sp_to_site_cat.push_back(i);

        if (!is_surface_set && n_site_categ == 0 ){
            throw InvalidInputError("SurfaceProperties", xml_surf_props.tag())
            << "In the surface properties at least one type of surface or sites"
            << "should be provided.";
        }

    }

//==============================================================================
    /**
     * Destructor.
     */
    ~SurfacePropertiesDetailed(){ }

//==============================================================================

    /**
     * Returns the index of the surface species. Firstly the species in sites
     * and the empty sites according to the order they are declared in the
     * input file and then the surface species composing the
     * surface. All of them are following the gas phase species.
     */
    int surfaceSpeciesIndex(const std::string& str_sp) const {
        // Looping over sites
        for (int i = 0; i < n_site_sp; i++) {
            if (v_site_sp[i] == str_sp)
                return n_gas_sp + i;
        }
        // Looping over surface composition
        for (int i = 0; i < n_surf_comp_sp; i++) {
            if (v_surf_sp[i] == str_sp)
                return n_gas_sp + n_site_sp + i;
        }
        // Not found
        return -1;
    }

//==============================================================================
    /**
     * Returns the gas phase species associated with the surface species.
     * For empty sites a specific value is returned equal to -2.
     */
    int surfaceToGasIndex(const int& i_surf_sp) const {
        if(i_surf_sp > n_gas_sp + n_tot_surf_sp) return -1;

        if (i_surf_sp > n_gas_sp + n_site_sp - 1) {
            return v_surf_to_gas_idx[i_surf_sp - (n_gas_sp + n_site_sp)];}
        else if (i_surf_sp > n_gas_sp - 1)
            return v_site_sp_to_gas_idx[i_surf_sp - n_gas_sp];
        else
            return i_surf_sp;
    }

//==============================================================================
    /**
     * Returns the total number of surface species.
     */
    size_t nSurfaceSpecies() const { return n_tot_surf_sp; }

//==============================================================================
    /**
     * Returns to which site category the site species belong.
     */
    int siteSpeciesToSiteCategoryIndex(const int& i_site_sp) const {
        std::cout << "Site Species = " << i_site_sp << std::endl;
        if (i_site_sp > n_gas_sp - 1 &&
            i_site_sp < n_gas_sp + n_tot_surf_sp - 1){
            std::cout << "position in vector = " << i_site_sp - n_gas_sp << std::endl;
            std::cout << "Returning = " << v_site_sp_to_site_cat[i_site_sp - n_gas_sp] << std::endl;
            return v_site_sp_to_site_cat[i_site_sp - n_gas_sp];
        }
        return -1;
    }

//==============================================================================
    /**
     * Returns the total number of site species.
     */
    size_t nSiteSpecies() const { return n_site_sp; }

//==============================================================================
    /**
     * Returns the total number of different site categories.
     */
    size_t nSiteCategories() const { return n_site_categ; }

//==============================================================================
private:
    void parseAblationSpecies(
        const std::string& species,
        const std::string& label)
    {
        std::istringstream iss(species);
        vector<std::string> v_species;
        std::copy(
            std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(v_species));

        for (int i_sp = 0; i_sp < v_species.size(); ++i_sp) {
            int id_sp = m_thermo.speciesIndex(v_species[i_sp]);

            if (id_sp == -1) {
                throw InvalidInputError("SurfaceProperties",
                    v_species[i_sp]) << "Surface species " <<
                    v_species[i_sp] << " is not " <<
                    "a species of the gas mixture!";
            }

            v_surf_to_gas_idx.push_back(id_sp);
            v_surf_sp.push_back(v_species[i_sp] + '-' + label);
        }
    }

//==============================================================================
    void parseSiteSpecies(
        const std::string& species,
        const std::string& label)
    {
        std::istringstream iss(species);
        vector<std::string> v_species;
        std::copy(
            std::istream_iterator<std::string>(iss),
            std::istream_iterator<std::string>(),
            std::back_inserter(v_species));

        const int non_empty =  v_species.size();
        for (int i_sp = 0; i_sp < non_empty; ++i_sp) {
            int id_sp = m_thermo.speciesIndex(v_species[i_sp]);

            if (id_sp == -1) {
                throw InvalidInputError("SurfaceProperties",
                v_species[i_sp]) << "Site species " <<
                v_species[i_sp] << " is not " <<
                "a species of the gas mixture!";
            }

            v_site_sp_to_gas_idx.push_back(id_sp);
            v_site_sp.push_back(v_species[i_sp] + '-' + label);
        }
        size_t const empty_site = 1;
        size_t tot_sp_in_site = empty_site + non_empty;
        v_sp_in_site.push_back(tot_sp_in_site);
    }

private:
    const Mutation::Thermodynamics::Thermodynamics& m_thermo;
    const XmlElement& xml_surf_props;

    const size_t n_gas_sp;

//  Info about the sites
    int n_site_categ;
    vector<int> v_sp_in_site; // Number of different species in each site
    vector<double> v_sigma; // Surface site number density per m2
    int n_site_sp;
    vector<std::string> v_site_sp; // Contains all species chemical symbol
    vector<int> v_site_sp_to_gas_idx; // Maps sites to gas idx
    vector<int> v_site_sp_to_site_cat; // Maps sites to site category


//  Info about the surface composition
    bool is_surface_set;
    int n_surf_comp_sp;
    std::vector<std::string> v_surf_sp;
    std::vector<int> v_surf_to_gas_idx;

    int n_tot_surf_sp;
};

ObjectProvider<
    SurfacePropertiesDetailed, SurfaceProperties>
    surface_properties_detailed("detailed");

    } // namespace GasSurfaceInteraction
} // namespace Mutation