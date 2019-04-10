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
    "Tests for the bulk chemistry and thermochemical properties.",
    "[gsi]"
)
{
    const double tol = std::numeric_limits<double>::epsilon();
    Mutation::GlobalOptions::workingDirectory(TEST_DATA_FOLDER);

    SECTION("Testing the pyrolysis Goldstein rates.")
    {
        // Setting up M++
        MixtureOptions opts("bulk_air5_goldstein_RRHO_ChemNonEq1T");
        Mixture mix(opts);

        // Setting up
        const int set_state_with_rhoi_T = 1;
        const size_t ns = mix.nSpecies();
        ArrayXd rhoi(ns);
        const int nr = 2;
        const int nps = 2;

        // Setup checks
        CHECK(nr == mix.nSurfaceReactions());
        CHECK(nps == mix.nPyrolysingSolids());

        // Bulk Properties
        ArrayXd rho_ps(nps);
        ArrayXd rho_ps_i(nps); rho_ps_i << 4., 3.;
        ArrayXd rho_ps_f(nps); rho_ps_f << 2., 0.;

        // Necessary data
        ArrayXd wdot(nps);
        ArrayXd wdot_mpp(nps);
        ArrayXd wdot_tot(ns+nps);

        // Conditions Solids Initial Density
        size_t nT = 10;
        for (int iT = 0; iT < nT; iT++) {
            const double P = 10000.;
            double Ts = 500. + iT * 200;
            mix.equilibrate(Ts, P);
            mix.densities(rhoi.data());

            mix.setSurfaceState(rhoi.data(), &Ts, set_state_with_rhoi_T);

            // Getting production rates from Mutation++
            mix.setPyrolysingSolidDensities(rho_ps_i.data());
            mix.surfaceReactionRatesPerReaction(wdot_mpp.data());
            mix.surfaceReactionRatesGasAndSolid(wdot_tot.data());

            wdot(0) = 1. * (rho_ps_i(0) - rho_ps_f(0)) * exp(-500. / Ts);
            wdot(1) = 2. * rho_ps_i(1) * pow( (rho_ps_i(1) - rho_ps_f(1))
                / (rho_ps_i(1)), 2) * exp(-1200. / Ts);

            CHECK(wdot(0) == Approx(wdot_mpp(0)).epsilon(tol));
            CHECK(wdot(1) == Approx(wdot_mpp(1)).epsilon(tol));

            // Sum of mass gas = mass solid
            CHECK(wdot_tot.head(ns).sum() ==
                Approx(-wdot_tot.tail(nps).sum()).epsilon(tol));
        }

        // Conditions Solids Intermediate Density
        rho_ps << 3., .5;
        for (int iT = 0; iT < nT; iT++) {
            const double P = 10000.;
            double Ts = 500. + iT * 200;
            mix.equilibrate(Ts, P);
            mix.densities(rhoi.data());

            mix.setSurfaceState(rhoi.data(), &Ts, set_state_with_rhoi_T);

            // Getting it from Mutation++
            mix.setPyrolysingSolidDensities(rho_ps_i.data());
            mix.surfaceReactionRatesPerReaction(wdot_mpp.data());
            mix.surfaceReactionRatesGasAndSolid(wdot_tot.data());

            wdot(0) = 1. * (rho_ps_i(0) - rho_ps_f(0)) * exp(-500. / Ts);
            wdot(1) = 2. * rho_ps_i(1) * pow( (rho_ps_i(1) - rho_ps_f(1)) / (rho_ps_i(1)), 2) * exp(-1200. / Ts);

            CHECK(wdot(0) == Approx(wdot_mpp(0)).epsilon(tol));
            CHECK(wdot(1) == Approx(wdot_mpp(1)).epsilon(tol));

            // Sum of mass gas = mass solid
            CHECK(wdot_tot.head(ns).sum() ==
                Approx(-wdot_tot.tail(nps).sum()).epsilon(tol));
        }

        // Conditions Solids Final Density
        for (int iT = 0; iT < nT; iT++) {
            const double P = 10000.;
            double Ts = 500. + iT * 200;
            mix.equilibrate(Ts, P);
            mix.densities(rhoi.data());

            mix.setSurfaceState(rhoi.data(), &Ts, set_state_with_rhoi_T);

            // Getting it from Mutation++
            mix.setPyrolysingSolidDensities(rho_ps_f.data());
            mix.surfaceReactionRatesPerReaction(wdot_mpp.data());

            wdot(0) = 0.;
            wdot(1) = 0.;

            CHECK(wdot(0) == Approx(wdot_mpp(0)).epsilon(tol));
            CHECK(wdot(1) == Approx(wdot_mpp(1)).epsilon(tol));

            // Sum of mass gas = mass solid
            CHECK(wdot_tot.head(ns).sum() ==
                Approx(-wdot_tot.tail(nps).sum()).epsilon(tol));
        }
    }

        SECTION("Testing solid heat capacity and thermal conductivity.")
    {
        // Setting up M++
        /* MixtureOptions opts("bulk_air5_goldstein_RRHO_ChemNonEq1T");
        Mixture mix(opts);

        // Setting up
        const int set_state_with_rhoi_T = 1;
        const size_t ns = mix.nSpecies();
        ArrayXd rhoi(ns);
        const int nr = 2;
        const int nps = 2;

        // Setup checks
        CHECK(nr == mix.nSurfaceReactions());
        CHECK(nps == mix.nPyrolysingSolids());

        // Bulk Properties
        ArrayXd rho_ps(nps);
        ArrayXd rho_ps_i(nps); rho_ps_i << 4., 3.;
        ArrayXd rho_ps_f(nps); rho_ps_f << 2., 0.;

        // Necessary data
        ArrayXd wdot(nps);
        ArrayXd wdot_mpp(nps);
        ArrayXd wdot_tot(ns+nps);

        // Conditions Solids Initial Density
        size_t nT = 10;
        for (int iT = 0; iT < nT; iT++) {
        } */
    }
}
