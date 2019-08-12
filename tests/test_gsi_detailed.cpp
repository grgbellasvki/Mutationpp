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

#include "mutation++.h"
#include "Configuration.h"
#include "TestMacros.h"
#include <catch.hpp>
#include <Eigen/Dense>

#include "SurfaceProperties.h"

using namespace Mutation;
using namespace Catch;
using namespace Eigen;

TEST_CASE("Detailed surface chemictry tests.","[gsi]")
{
    const double tol = std::numeric_limits<double>::epsilon();
    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Surface Species and Coverage.")
    {
        // Setting up M++
        MixtureOptions opts("smb_detailed_coverage_NASA9_ChemNonEq1T");
        Mixture mix(opts);

        CHECK(mix.nSpecies() == 7);

        // Check global options
        CHECK(mix.nSurfaceReactions() == 0);
        CHECK(mix.getSurfaceProperties().nSurfaceSpecies() == 10);
        CHECK(mix.getSurfaceProperties().nSiteSpecies() == 9);

        // Check Species
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("N-s") == 8);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("O-s") == 9);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("N2-s") == 10);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("NO-c") == 12);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("C-b") == 16);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("O-c") == 13);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("C-p") == 15);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("A") == -1);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("s") == 7);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("c") == 11);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("p") == 14);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("b") == -1);

        // Check surface species association with gaseous species
        CHECK( mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("N-s")) == 0);
        CHECK( mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O-s")) == 1);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("N2-s")) == 3);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("NO-c")) == 2);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("C-b")) == 5);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O-c")) == 1);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("C-p")) == 5);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("s")) == -2);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("c")) == -2);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("p")) == -2);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("b")) == -1);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(100) == -1);

        // Check site species map correctly to the site category
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("N-s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O-c")) == 1);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O-s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("N2-s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("NO-c")) == 1);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("C-b")) == -1);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("C-p")) == 2);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("c")) == 1);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("p")) == 2);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("b")) == -1);

        CHECK(mix.getSurfaceProperties().nSiteDensityInCategory(0) == 3.e19);
        CHECK(mix.getSurfaceProperties().nSiteDensityInCategory(1) == 7.e19);
        CHECK(mix.getSurfaceProperties().nSiteDensityInCategory(2) == 1.e20);
        CHECK(mix.getSurfaceProperties().nSiteDensityInCategory(3) == -1);

    }

    SECTION("Adsorption-Desorption Equilibrium.")
    {
        // Setting up M++
        MixtureOptions opts("smb_ads_des_eq_NASA9_ChemNonEq1T");
        Mixture mix(opts);

        size_t ns = 4;
        size_t nr = 2;
        CHECK(mix.nSpecies() == ns);
        CHECK(mix.nSurfaceReactions() == 2);

        const size_t iO = 0;

        ArrayXd wdot(ns); ArrayXd wdotmpp(ns);
        ArrayXd rates(nr); ArrayXd ratesmpp(nr);
        wdot.setZero(); wdotmpp.setZero();

        ArrayXd mm = mix.speciesMw();

        // Gas conditions
        double P = 100000.;

        // Equilibrium Surface
        double T = 300.; // K
        double dT = 100.; // K
        for (int i = 0; i < 30; i++) {
            T += i * dT;

            mix.equilibrate(T, P);
            double nO = mix.X()[iO] * mix.numberDensity();

            mix.surfaceReactionRatesPerReaction(ratesmpp.data());

            const double B = 1.e-5 * NA;
            double kf1 = sqrt(RU * T / (2 * PI * mm(iO)));
            double kf2 = 2 * PI * mm(iO) * KB * KB / (B * HP * HP * HP);
            kf2 *= exp(-100000./T);

            // Changing the surface coverage
            const double a = 1.e-1;
            double s = B;

            for (int j = 0; j < 4; j++){
                s *= a;
                double Os = B - s;

                rates(0) = kf1 * nO * s;
                rates(1) = kf2 * Os;

                wdot(iO) = mm(iO) / NA *(-rates(0)+rates(1)); // O
            }

            // Equilibrium surface coverage
            s = B * kf2 / (nO * kf1 + kf2);
            double Os = B - s;

            rates(0) = kf1 * nO * s;
            rates(1) = kf2 * Os;

            wdot(0) = mm(0) / NA *(-rates(0)+rates(1)); // O
        }

    }

}