/**
 * @file GSIReactionAblation.cpp
 *
 * @brief Implements a surface reaction of ablation type.
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
#include "Utilities.h"

#include "GSIRateLaw.h"
#include "GSIReaction.h"
#include "SurfaceState.h"
#include "SurfaceProperties.h"


using namespace std;
using namespace Mutation::Utilities::Config;

namespace Mutation {
    namespace GasSurfaceInteraction {

class GSIReactionDetailed : public GSIReaction
{
public:
    GSIReactionDetailed(ARGS args)
        : GSIReaction(args)
    {
        assert(args.s_iter_reaction->tag() == "reaction");

        args.s_iter_reaction->getAttribute(
            "formula", m_formula,
            errorNoFormulainReaction());

        parseFormula(args.s_thermo, *(args.s_iter_reaction));

        const Mutation::Utilities::IO::XmlElement& node_rate_law =
            *(args.s_iter_reaction->begin());

        DataGSIRateLaw data_gsi_rate_law = {
            args.s_thermo,
            args.s_transport,
            m_surf_props,
            node_rate_law,
            m_reactants,
            m_products };

        mp_rate_law = Factory<GSIRateLaw>::create(
            node_rate_law.tag(), data_gsi_rate_law);

        if (mp_rate_law == NULL) {
            args.s_iter_reaction->parseError(
                "A rate law must be provided for this reaction!");
        }

        // Check for charge and mass conservation
        const size_t ns = args.s_thermo.nSpecies();
        const size_t ne = args.s_thermo.nElements();
        vector<int> v_sums(ne);
        const size_t n_site_c = m_surf_props.nSiteCategories();
        vector<int> v_sums_site(n_site_c);

        std::fill(v_sums.begin(), v_sums.end(), 0);
        // Reactants
        for (int ir = 0; ir < m_reactants.size(); ++ir) {
            int reactant = m_reactants[ir];
            for (int kr = 0; kr < ne; ++kr) {
                if (reactant < ns) {
                    v_sums[kr] += args.s_thermo.elementMatrix()(reactant, kr);
                } else {
                    int gas = m_surf_props.surfaceToGasIndex(reactant);
                    if (gas > -1)
                        v_sums[kr] += args.s_thermo.elementMatrix()(gas, kr);
                }
            }
            int site = m_surf_props.siteSpeciesToSiteCategoryIndex(reactant);
            if (site != -1) v_sums_site[site] += 1;
        }
        // Reactants in the surface phase
        for (int ir = 0; ir < m_reactants_surf.size(); ++ir)
            for (int kr = 0; kr < ne; ++kr)
                v_sums[kr] += args.s_thermo.elementMatrix()(
                    m_surf_props.surfaceToGasIndex(m_reactants_surf[ir]), kr);
        // Products
        for (int ip = 0; ip < m_products.size(); ++ip) {
            int product = m_products[ip];
            for (int kp = 0; kp < ne; ++kp) {
                if (product < ns) {
                    v_sums[kp] -= args.s_thermo.elementMatrix()(product, kp);
                } else {
                    int gas = m_surf_props.surfaceToGasIndex(product);
                    if (gas > -1)
                        v_sums[kp] -= args.s_thermo.elementMatrix()(gas, kp);
                }
            }
            int site = m_surf_props.siteSpeciesToSiteCategoryIndex(product);
            if (site != -1) v_sums_site[site] -= 1;
        }

        // Products in the surface phase
        for (int ip = 0; ip < m_products_surf.size(); ++ip){
            for (int kp = 0; kp < ne; ++kp) {
                v_sums[kp] -= args.s_thermo.elementMatrix()(
                    m_surf_props.surfaceToGasIndex(m_products_surf[ip]), kp);
            }
        }

        for (int i = 0; i < ne; ++i)
            m_conserves &= (v_sums[i] == 0);
        for (int i = 0; i < n_site_c; ++i)
            m_conserves &= (v_sums_site[i] == 0);

         if (m_conserves == false)
                throw InvalidInputError("formula", m_formula)
                    << "Reaction " << m_formula
                    << " does not conserve mass or sites.";
    }

//==============================================================================

    ~GSIReactionDetailed(){
        if (mp_rate_law != NULL){ delete mp_rate_law; }
    }

//==============================================================================

    void parseSpecies(
        std::vector<int>& species, std::vector<int>& species_surf,
        std::string& str_chem_species,
        const Mutation::Utilities::IO::XmlElement& node_reaction,
        const Mutation::Thermodynamics::Thermodynamics& thermo)
    {

        // State-Machine states for parsing a species formula
        enum ParseState {
            coefficient,
            name,
            plus
        } state = coefficient;

        size_t c = 0;
        size_t s = 0;
        size_t e = 0;
        int nu   = 1;
        bool add_species = false;

        Mutation::Utilities::String::eraseAll(str_chem_species, " ");

        // Loop over every character
        while (c != str_chem_species.size()) {
            switch(state) {
                case coefficient:
                    if (isdigit(str_chem_species[c])) {
                        nu = atoi(str_chem_species.substr(c, 1).c_str());
                        s = c + 1;
                    } else {
                        nu = 1;
                        s = c;
                    }
                    state = name;
                    break;
                case name:
                    if (str_chem_species[c] == '+')
                        state = plus;
                    break;
                case plus:
                    if (str_chem_species[c] != '+') {
                        e = c - 2;
                        c--;
                        add_species = true;
                        state = coefficient;
                    }
                    break;
            }
            if (c == str_chem_species.size() - 1) {
                add_species = true;
                e = c;
            }

            int index;
            if (add_species) {
                index = thermo.speciesIndex(
                    str_chem_species.substr(s, e-s+1));

                if (index == -1) {
                    index = m_surf_props.surfaceSpeciesIndex(
                        str_chem_species.substr(s, e-s+1));
                }

                if(index == -1) {
                    node_reaction.parseError((
                        "Species " + str_chem_species.substr(s, e-s+1)
                        + " is not in the mixture list or a species of the "
                        "surface!").c_str());
                }

                if (index < (thermo.nSpecies() + m_surf_props.nSiteSpecies())) {
                    for (int i_nu = 0; i_nu < nu; i_nu++)
                        species.push_back(index);
                } else {
                    for (int i_nu = 0; i_nu < nu; i_nu++)
                        species_surf.push_back(index);
                }

                add_species = false;
            }

            // Move on the next character
            c++;
        }

        // Sort the species vector
        std::sort(species.begin(), species.end());
    }
};

ObjectProvider<
    GSIReactionDetailed, GSIReaction>
    detailed_reaction("detailed");

    } // namespace GasSurfaceInteraction
}  // namespace Mutation
