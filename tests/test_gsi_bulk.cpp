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
        MixtureOptions opts("bulk_goldstein_NASA-9_ChemNonEq1T");
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

            // Getting production rates from Mutation++
            mix.setPyrolysingSolidDensities(rho_ps_i.data());
            mix.surfaceReactionRatesPerReaction(wdot_mpp.data());
            mix.surfaceReactionRatesGasAndSolid(wdot_tot.data());

            wdot(0) = 1. * (rho_ps_i(0) - rho_ps_f(0)) * exp(-500. / Ts);
            wdot(1) = 2. * rho_ps_i(1) * pow( (rho_ps_i(1) - rho_ps_f(1)) / (rho_ps_i(1)), 2) * exp(-1200. / Ts);

            CHECK(wdot(0) == Approx(wdot_mpp(0)).epsilon(tol));
            CHECK(wdot(1) == Approx(wdot_mpp(1)).epsilon(tol));

            // Making sure mass is conserved
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

            // Getting production rates from Mutation++
            mix.setPyrolysingSolidDensities(rho_ps_f.data());
            mix.surfaceReactionRatesPerReaction(wdot_mpp.data());

            wdot(0) = 0.;
            wdot(1) = 0.;

            CHECK(wdot(0) == Approx(wdot_mpp(0)).epsilon(tol));
            CHECK(wdot(1) == Approx(wdot_mpp(1)).epsilon(tol));

            // Making sure mass is conserved
            CHECK(wdot_tot.head(ns).sum() ==
                Approx(-wdot_tot.tail(nps).sum()).epsilon(tol));
        }
    }

        SECTION("Testing solid heat capacity and thermal conductivity.")
    {
        // Setting up M++
        MixtureOptions opts("bulk_thermal_props_NASA-9_ChemNonEq1T");
        Mixture mix(opts);

        // Setting up
        const int set_state_with_rhoi_T = 1;
        const size_t ns = mix.nSpecies();
        ArrayXd rhoi(ns);
        const int nr = 3;
        const int nps = 3;

        // Setup checks
        CHECK(nr == mix.nSurfaceReactions());
        CHECK(nps == mix.nPyrolysingSolids());

        // Bulk Properties
        // ArrayXd rho_ps(nps);
        ArrayXd rho_ps_i(nps); rho_ps_i << 15., 12., 60.;
        ArrayXd rho_ps_f(nps); rho_ps_f << 5., 0., 20.;

        // Necessary data
        double cp_mpp;
        const int ndim = 3;
        ArrayXd v_lambda_mpp(ndim);

        // Data for verification
        const int ncoeff = 7;
        const double Tref = 273.15;
        VectorXd v_cp_d(ncoeff); v_cp_d <<
            5.53472512800491e-9, 1.25740916978889e-05, -3.31957114005794e-2,
            1.7174710171051e1, 3.79730067921415e2, 0., 0.;
        VectorXd v_cp_f(ncoeff); v_cp_f <<
            0., 0., 0., 1305., 1933.8, -394500., 0.;

        VectorXd v_lambda_x_d(ncoeff); v_lambda_x_d <<
            0., 0., -6.79172921619856e-06, 0.0094013164, -1.3718533625, 0., 0.;
        VectorXd v_lambda_y_d(ncoeff); v_lambda_y_d <<
            0., 0., -2.35322821162944e-06, 0.0031852894, -0.5386417424, 0., 0.;
        VectorXd v_lambda_z_d(ncoeff); v_lambda_z_d <<
            0., 0., -2.35322821162944e-06, 0.0031852894, -0.5386417424, 0., 0.;

        VectorXd v_lambda_x_f(ncoeff); v_lambda_x_f <<
            0., 0., 8.01548924555013e-07, -0.0018332383, 2.2255293753, 0., 0.;
        VectorXd v_lambda_y_f(ncoeff); v_lambda_y_f <<
            0., 0., 6.57518729361747e-07, -0.0010739585, 0.8503031122, 0., 0.;
        VectorXd v_lambda_z_f(ncoeff); v_lambda_z_f <<
            0., 0., 6.57518729361747e-07, -0.0010739585, 0.8503031122, 0., 0.;

        VectorXd vT(ncoeff);
        ArrayXd v_lambda(ndim);
        double cp;

        // Conditions virgin solid
        size_t nT = 20;
        for (int iT = 0; iT < nT; iT++) {
            const double P = 10000.;
            double Ts = 500. + iT * 50;
            mix.equilibrate(Ts, P);
            mix.densities(rhoi.data());

            // Setting up surface conditions
            mix.setSurfaceState(rhoi.data(), &Ts, set_state_with_rhoi_T);
            mix.setPyrolysingSolidDensities(rho_ps_i.data());

            // Getting heat capacity and thermal conductivity
            mix.solidHeatCapacity(&cp_mpp);
            mix.solidEffectiveThermalConductivity(v_lambda_mpp.data());

            // Recomputing cp and lambda and checking
            vT(4) = 1.;
            vT(3) = Ts;
            vT(2) = vT(3) * Ts;
            vT(1) = vT(2) * Ts;
            vT(0) = vT(1) * Ts;
            vT(5) = 1/Ts;
            vT(6) = vT(5) / Ts;

            v_lambda(0) = v_lambda_x_d.dot(vT);
            v_lambda(1) = v_lambda_y_d.dot(vT);
            v_lambda(2) = v_lambda_z_d.dot(vT);

            CHECK(v_lambda(0) == Approx(v_lambda_mpp(0)).epsilon(tol));
            CHECK(v_lambda(1) == Approx(v_lambda_mpp(1)).epsilon(tol));
            CHECK(v_lambda(2) == Approx(v_lambda_mpp(2)).epsilon(tol));

            vT(4) = 1.;
            vT(3) = Ts-Tref;
            vT(2) = vT(3) * (Ts-Tref);
            vT(1) = vT(2) * (Ts-Tref);
            vT(0) = vT(1) * (Ts-Tref);
            vT(5) = 1. / (Ts-Tref);
            vT(6) = vT(5) / (Ts-Tref);

            cp = v_cp_d.dot(vT);
            CHECK(cp == Approx(cp_mpp).epsilon(tol));
        }

        // Conditions decomposed solid
        for (int iT = 0; iT < nT; iT++) {
            const double P = 10000.;
            double Ts = 500. + iT * 50;
            mix.equilibrate(Ts, P);
            mix.densities(rhoi.data());

            // Setting up surface conditions
            mix.setSurfaceState(rhoi.data(), &Ts, set_state_with_rhoi_T);
            mix.setPyrolysingSolidDensities(rho_ps_f.data());

            // Getting heat capacity and thermal conductivity
            mix.solidHeatCapacity(&cp_mpp);
            mix.solidEffectiveThermalConductivity(v_lambda_mpp.data());

            // Recomputing cp and lambda and checking
            vT(4) = 1.;
            vT(3) = Ts;
            vT(2) = vT(3) * Ts;
            vT(1) = vT(2) * Ts;
            vT(0) = vT(1) * Ts;
            vT(5) = 1/Ts;
            vT(6) = vT(5) / Ts;

            v_lambda(0) = v_lambda_x_f.dot(vT);
            v_lambda(1) = v_lambda_y_f.dot(vT);
            v_lambda(2) = v_lambda_z_f.dot(vT);

            CHECK(v_lambda(0) == Approx(v_lambda_mpp(0)).epsilon(tol));
            CHECK(v_lambda(1) == Approx(v_lambda_mpp(1)).epsilon(tol));
            CHECK(v_lambda(2) == Approx(v_lambda_mpp(2)).epsilon(tol));

            cp = v_cp_f.dot(vT);
            CHECK(cp == Approx(cp_mpp).epsilon(tol));
        }
    }
}
