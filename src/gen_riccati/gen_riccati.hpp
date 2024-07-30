
/*
Fatrop - A fast trajectory optimization solver
 Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.

This file is part of Fatrop.

Fatrop is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fatrop is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Fatrop.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include "blasfeo_wrapper.hpp"
#include "cocp.hpp"
#include "numeric_vector.hpp"
#include <algorithm>
#include <cmath>
namespace gen_riccati
{
    bool check_reg(const int m, MAT *sA, const int ai, const int aj);
    class OCPLSRiccati
    {
    public:
        OCPLSRiccati(const COCP &cocp);
        // solve a KKT system
        int solve_pd_sys(
            COCP &OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam);
        // solve a KKT system
        int solve_pd_sys_degenerate(
            COCP &OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam);
        // solve a KKT system
        int
        solve_pd_sys_normal(
            COCP &OCP,
            const double inertia_correction,
            const FatropVecBF &ux,
            const FatropVecBF &lam);
        int get_rhs(
            COCP &OCP,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g);
        int compute_pd_sys_times_vec(
            COCP &OCP,
            const double inertia_correction_w,
            const double inertia_correction_c,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g);
        int solve_rhs(
            COCP &OCP,
            const FatropVecBF &ux,
            const FatropVecBF &lam,
            const FatropVecBF &rhs_rq,
            const FatropVecBF &rhs_b,
            const FatropVecBF &rhs_g);
        FatropMemoryMatBF Ppt;
        FatropMemoryMatBF Hh;
        FatropMemoryMatBF AL;
        FatropMemoryMatBF RSQrqt_tilde;
        FatropMemoryMatBF Ggt_stripe;
        FatropMemoryMatBF Ggt_tilde;
        FatropMemoryMatBF GgLt;
        FatropMemoryMatBF RSQrqt_hat;
        FatropMemoryMatBF Llt;
        FatropMemoryMatBF Llt_shift; // needed because feature not implemented yet
        FatropMemoryMatBF GgIt_tilde;
        FatropMemoryMatBF GgLIt;
        FatropMemoryMatBF HhIt;
        FatropMemoryMatBF PpIt_hat;
        FatropMemoryMatBF LlIt;
        FatropMemoryVecBF rhs_rq;
        FatropMemoryVecBF rhs_b;
        FatropMemoryVecBF rhs_g;
        FatropMemoryVecBF rhs_rq2;
        FatropMemoryVecBF rhs_b2;
        FatropMemoryVecBF rhs_g2;
        FatropMemoryVecBF ux_test;
        FatropMemoryVecBF lam_test;
        FatropMemoryVecBF v_Ppt;
        FatropMemoryVecBF v_Hh;
        FatropMemoryVecBF v_AL;
        FatropMemoryVecBF v_RSQrqt_tilde;
        FatropMemoryVecBF v_Ggt_stripe;
        FatropMemoryVecBF v_Ggt_tilde;
        FatropMemoryVecBF v_GgLt;
        FatropMemoryVecBF v_RSQrqt_hat;
        FatropMemoryVecBF v_Llt;
        FatropMemoryVecBF v_Llt_shift;
        FatropMemoryVecBF v_GgIt_tilde;
        FatropMemoryVecBF v_GgLIt;
        FatropMemoryVecBF v_HhIt;
        FatropMemoryVecBF v_PpIt_hat;
        FatropMemoryVecBF v_LlIt;
        MemoryPermMat Pl;
        MemoryPermMat Pr;
        MemoryPermMat PlI;
        MemoryPermMat PrI;
        NumericVector gamma;
        NumericVector rho;
        int rankI = 0;
        struct LastUsed
        {
            int rankI = 0;
            double inertia_correction_w = 0;
            double inertia_correction_c = 0;
            double kappa_d = 0;
            double mu = 0;
        } lastused_;
        const NumericVector ux_offs;
        const NumericVector g_offs;
        const NumericVector dyn_offs;
        const NumericVector dyn_eq_offs;
        const int max_nu;
        const int max_nx;
        const int max_ng;
        bool it_ref = false;
    };
} // namespace gen_riccati