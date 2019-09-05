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
    const double tol = 100. * std::numeric_limits<double>::epsilon();
    const double tol_det = 1.e6 * std::numeric_limits<double>::epsilon();

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
        const int set_state_rhoi_T = 1;

        ArrayXd v_rhoi(ns);
        ArrayXd wdot(ns); ArrayXd wdotmpp(ns);
        ArrayXd rates(nr); ArrayXd ratesmpp(nr);
        wdot.setZero(); wdotmpp.setZero();

        ArrayXd mm = mix.speciesMw();

        CHECK(mix.getSurfaceProperties().isSurfaceCoverageSteady() == false);
        ArrayXd v_surf_cov_frac(mix.getSurfaceProperties().nSurfaceSpecies());
        ArrayXd v_surf_cov_ss_frac(mix.getSurfaceProperties().nSurfaceSpecies());

        // Equilibrium Surface
        double P = 1000.;
        double T; // K
        double dT = 500.; // K
        for (int i = 12; i < 30; i++) { // 30
            T = (i+1) * dT; // += i*dT;

            mix.equilibrate(T, P);
            mix.densities(v_rhoi.data());

            mix.setSurfaceState(v_rhoi.data(), &T, set_state_rhoi_T);
            double nO = mix.X()[iO] * mix.numberDensity();

            const double B = 1.e20;
            double kf1 = sqrt(RU * T / (2 * PI * mm(iO)));
            double kf2 = 2 * PI * mm(iO) / NA * KB * KB * T * T / (HP * HP * HP);
            kf2 *= exp(-100000./T);

            // Changing the surface coverage
            mix.getSurfaceProperties().setIsSurfaceCoverageSteady(false);
            v_surf_cov_frac(0) = 1.;
            const double a = 1.e-1;

            for (int j = 0; j < 0; j++){
                v_surf_cov_frac(0) *= a;                      // s
                v_surf_cov_frac(1) = 1. - v_surf_cov_frac(0); // Os

                mix.getSurfaceProperties().setSurfaceSiteCoverageFrac(v_surf_cov_frac);
                mix.surfaceReactionRatesPerReaction(ratesmpp.data());

                rates(0) = kf1 * nO * v_surf_cov_frac(0);
                rates(1) = kf2 * v_surf_cov_frac(1);
                CHECK(rates(0) == Approx(ratesmpp(0)).epsilon(tol));
                CHECK(rates(1) == Approx(ratesmpp(1)).epsilon(tol));

                wdot(iO) = mm(iO) / NA * (rates(0)-rates(1)); // From the surface point of view
                mix.surfaceReactionRates(wdotmpp.data());

                CHECK(wdot(0) == Approx(wdotmpp(0)).epsilon(tol));
                CHECK(wdot(1) == Approx(wdotmpp(1)).epsilon(tol));
                CHECK(wdot(2) == Approx(wdotmpp(2)).epsilon(tol));
                CHECK(wdot(3) == Approx(wdotmpp(3)).epsilon(tol));
            }

            // Equilibrium surface coverage
            mix.getSurfaceProperties().setIsSurfaceCoverageSteady(true);
            mix.surfaceReactionRatesPerReaction(ratesmpp.data());
            mix.surfaceReactionRates(wdotmpp.data());

            v_surf_cov_ss_frac = mix.getSurfaceProperties().getSurfaceSiteCoverageFrac();

            // THIS IS WRONG!
            v_surf_cov_frac(0) = kf2 / (nO * kf1 + kf2);
            // double asdf = 1.e0 - v_surf_cov_ss_frac(0);
            v_surf_cov_frac(1) = 1.e0 - v_surf_cov_frac(0);

            // CHECK(v_surf_cov_ss_frac(0) == Approx(v_surf_cov_frac(0)).epsilon(tol));
            // CHECK(v_surf_cov_ss_frac(1) == Approx(v_surf_cov_frac(1)).epsilon(tol));
            // std::cout << "i = " << i << std::endl;
            // std::cout << "    kf1 = " << kf1 << " kf2 = " << kf2 << std::endl;
            // std::cout << "HERE 1 = " << v_surf_cov_frac(0) << " HERE 2 = " << v_surf_cov_ss_frac(1) << std::endl;
            // std::cout << "Sum = " << v_surf_cov_frac(0) + v_surf_cov_ss_frac(1) << std::endl;
            // std::cout << "MPP  1 = " << v_surf_cov_ss_frac(0) << " MPP  2 = " << v_surf_cov_ss_frac(1) << std::endl;
            // std::cout << "Sum = " << v_surf_cov_ss_frac(0) + v_surf_cov_ss_frac(1) << std::endl;


            rates(0) = kf1 * nO * v_surf_cov_frac(0);
            rates(1) = kf2 * v_surf_cov_frac(1);
            // CHECK(rates(0) == Approx(ratesmpp(0)).epsilon(tol_det));
            CHECK(rates(1) == Approx(ratesmpp(1)).epsilon(tol_det));

            // std::cout << "is zero here= " << rates(0) << " " << rates(1)<< std::endl;
            // std::cout << "is zero mpp = " << ratesmpp(0) << " " << ratesmpp(1)<< std::endl;

            wdot(iO) = mm(iO) / NA * (rates(0)-rates(1)); // O
            //CHECK(wdot(0) == Approx(wdotmpp(0)).epsilon(tol));
            //CHECK(wdot(1) == Approx(wdotmpp(1)).epsilon(tol));
            //CHECK(wdot(2) == Approx(wdotmpp(2)).epsilon(tol));
            //CHECK(wdot(3) == Approx(wdotmpp(3)).epsilon(tol));

            // std::cout << "MPP = " << wdot(0) << " " << wdotmpp(0) << std::endl;
            // std::cout << "MPP = " << wdot(1) << " " << wdotmpp(1) << std::endl;
            // std::cout << "MPP = " << wdot(2) << " " << wdotmpp(2) << std::endl;
            // std::cout << "MPP = " << wdot(3) << " " << wdotmpp(3) << std::endl;
            // double in; std::cin >> in;

        }

    }

}