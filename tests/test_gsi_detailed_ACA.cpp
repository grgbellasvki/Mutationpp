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

void computeRates(double Rs) {
  Rs = 3.0; 
}

TEST_CASE("Detailed surface chemictry tests full ACA.","[gsi]")
{
    const double tol = 100. * std::numeric_limits<double>::epsilon();
    const double tol_det = 1.e2 * std::numeric_limits<double>::epsilon();

    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Surface Species and Coverage.")
    {
        // Setting up M++
        MixtureOptions opts("smb_full_ACA_NASA9_ChemNonEq1T");
        Mixture mix(opts);

        CHECK(mix.nSpecies() == 8);

        // Check global options
        CHECK(mix.nSurfaceReactions() == 20);
        CHECK(mix.getSurfaceProperties().nSurfaceSpecies() == 5);
        CHECK(mix.getSurfaceProperties().nSiteSpecies() == 4);

        // Check Species
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("s") == 8);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("N-s") == 9);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("O-s") == 10);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("O*-s") == 11);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("C-b") == 12);
        CHECK(mix.getSurfaceProperties().surfaceSpeciesIndex("A") == -1);

        // Check surface species association with gaseous species
        CHECK( mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("N-s")) == 2);
        CHECK( mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O-s")) == 0);
        CHECK( mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O*-s")) == 0);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("C-b")) == 4);

        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("s")) == -2);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("b")) == -1);
        CHECK(mix.getSurfaceProperties().surfaceToGasIndex(100) == -1);

        // Check site species map correctly to the site category
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("N-s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O-s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("O*-s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("C-b")) == -1);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("s")) == 0);
        CHECK(mix.getSurfaceProperties().siteSpeciesToSiteCategoryIndex(
            mix.getSurfaceProperties().surfaceSpeciesIndex("b")) == -1);

        CHECK(mix.getSurfaceProperties().nSiteDensityInCategory(0) == 6.022e18);
        CHECK(mix.getSurfaceProperties().nSiteDensityInCategory(3) == -1);
        
        // Check reaction rates
        double P = 1500; //Pa
        double T = 2000; // K
        
        int iO = 0;
        int iN = 2;
        int ns = mix.nSpecies();
        
        ArrayXd v_rhoi(ns);
        ArrayXd v_nd(ns);
        ArrayXd mm = mix.speciesMw();
        const int set_state_rhoi_T = 1;
        
        mix.equilibrate(T, P);
        mix.densities(v_rhoi.data());
        for (int i = 0; i < ns; ++i)
            v_nd(i) = v_rhoi(i) / mm(i);
        
        double B = 6.022e18;
        double FO = 1./B * sqrt(RU * T / (2 * PI * mm(iO)));
        double FO2 = 1./B * sqrt(RU * T / (2 * PI * 2 * mm(iO)));
        double FN = 1./B * sqrt(RU * T / (2 * PI * mm(iN)));
        
        //coefficients 
        double k1 = FO * 0.3;
        double k2 = 2.0 * PI * mm(iO) / NA * KB * KB * T * T / (HP * HP * HP) / B * exp(-44277.6/T);
        double k3 = FO * 100. * exp(-4000/T);
        double k4 = FO * exp(-500/T);
        double k5 = FO * 0.7;
        double k6 = 2.0 * PI * mm(iO) / NA * KB * KB * T * T / (HP * HP * HP) / B * exp(-96500.6/T);
        double k7 = FO * 1000. * exp(-4000/T);
        double k8 = sqrt(1.0/B) * sqrt(NA * PI * KB * T / (2.0 * mm(iO)))* 1.0e-3 * exp(-15000.0/ T);
        double k9 = sqrt(1.0/B) * sqrt(NA * PI * KB * T / (2.0 * mm(iO)))* 5.0e-5 * exp(-15000.0/ T);
        double k10 = FN * exp(-2500./T);
        double k11 = 2.0 * PI * mm(iN) / NA * KB * KB * T * T / (HP * HP * HP) / B * exp(-73971.6/T);
        double k12 = FN * 1.5 * exp(-7000./T);
        double k13 = FN * 0.5 * exp(-2000./T);
        double k14 = sqrt(1.0/(B)) * sqrt(NA*PI*KB*T/ (2.0 * mm(iN)))* 0.1 * exp(-21000.0/ T);
        double k15 = 1e8 * exp(-20676.0 / T);
        double k16 = FO2 / B * exp(-8000./T);
        double k17 = FO2 * 100. * exp(-4000./T);
        double k18 = FO2 * exp(-500./T);
        double k19 = FO2 / B * exp(-8000./T);
        double k20 = FO2 * 1000. * exp(-4000./T);
        
        double r;
        computeRates(r);
        
        //compute steady state
        mix.setSurfaceState(v_rhoi.data(), &T, set_state_rhoi_T);
        mix.getSurfaceProperties();
        
       

    }
    

}
