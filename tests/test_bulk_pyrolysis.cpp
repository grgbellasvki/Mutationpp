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
#include <catch/catch.hpp>
#include <eigen3/Eigen/Dense>

using namespace Mutation;
using namespace Catch;
using namespace Eigen;

TEST_CASE
(
    "Test for the bulk pyrolysis",
    "[gsi]"
)
{
    const double tol = std::numeric_limits<double>::epsilon();
    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Goldstein Pyrolysis.")
    {
        // Input Parameters
        double beta = .3; // K/s

        MixtureOptions opts("air5_carbon_gsi_bulk_pyro_goldstein_RRHO_ChemNonEq1T");
        Mixture mix(opts);

        const int npyrsol = 2;
        CHECK(npyrsol == mix.nPyrolysingSolids());

        const size_t with_rhoi_T = 1;
        const int ns = mix.nSpecies();
        const int nT = mix.nEnergyEqns();
        VectorXd v_rhoi(ns);
        VectorXd v_rho_ps(npyrsol);
        VectorXd v_wdot(ns+npyrsol);
        VectorXd v_T(nT);
        VectorXd v_Tb(nT);

        v_rhoi(0) = 6.23011e-15;
        v_rhoi(1) = 2.47146e-07;
        v_rhoi(2) = 0.000302155;
        v_rhoi(3) = 0.177305;
        v_rhoi(4) = 0.0537185;
        v_T(0) = 1500.;

        mix.setState(v_rhoi.data(), v_T.data(), with_rhoi_T);

        v_Tb(0) = 300.;
        double tlim = 2000.;
        double tstep = 1.;
        double t = 0.;
        v_rho_ps(0) = .8;
        v_rho_ps(1) = .8;
        while  (t < tlim)
        {
            mix.setSurfaceState(v_rhoi.data(), v_Tb.data(), with_rhoi_T);
            mix.setPyrolysingSolidDensities(v_rho_ps.data());

            mix.surfaceReactionRatesGasAndSolid(v_wdot.data());

        // Get Wall Production Rate

            v_Tb(0) += beta * tstep;
            t += tstep;
        }
    }
    SECTION("Zuram Pyrolysis.")
    {
    
    }

}