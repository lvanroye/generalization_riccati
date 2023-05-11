#include "gen_riccati.hpp"
using namespace std;
using namespace gen_riccati;
bool gen_riccati::check_reg(const int m, MAT *sA, const int ai, const int aj)
{
    for (int i = 0; i < m; i++)
    {
        if (MATEL(sA, ai + i, aj + i) < 1e-8)
            return false;
    }
    return true;
}
OCPLSRiccati::OCPLSRiccati(const COCP &cocp) : Ppt(cocp.nx + 1, cocp.nx, cocp.K),
                                               Hh(cocp.nx, cocp.nx + 1, cocp.K), // the number of eqs can never exceed nx
                                               AL(vector<int>(1, max(cocp.nu + cocp.nx + 1)), vector<int>(1, max(cocp.nx)), 1),
                                               RSQrqt_tilde(cocp.nu + cocp.nx + 1, cocp.nx + cocp.nu, cocp.K), // TODO, only save first rho rows (can never exceed nu)
                                               Ggt_stripe(vector<int>(1, max(cocp.nu + cocp.nx + 1)), vector<int>(1, max(cocp.nx + cocp.nu)), 1),
                                               Ggt_tilde(cocp.nu + cocp.nx + 1, cocp.nx + cocp.nu, cocp.K), // TODO, only save first rho rows (can never exceed nu)
                                               GgLt(vector<int>(1, max(cocp.nu + cocp.nx + 1)), vector<int>(1, max(cocp.nu + cocp.nx)), 1),
                                               RSQrqt_hat(vector<int>(1, max(cocp.nu + cocp.nx + 1)), vector<int>(1, max(cocp.nx + cocp.nu)), 1),
                                               Llt(cocp.nu + cocp.nx + 1, cocp.nu, cocp.K),
                                               Llt_shift(vector<int>(1, max(cocp.nu + cocp.nx + 1)), vector<int>(1, max(cocp.nu)), 1),
                                               GgIt_tilde(vector<int>(1, cocp.nx[0] + 1), vector<int>(1, cocp.nx[0]), 1),
                                               GgLIt(vector<int>(1, cocp.nx[0] + 1), vector<int>(1, cocp.nx[0]), 1),
                                               HhIt(vector<int>(1, cocp.nx[0] + 1), vector<int>(1, cocp.nx[0]), 1),
                                               PpIt_hat(vector<int>(1, cocp.nx[0] + 1), vector<int>(1, cocp.nx[0]), 1),
                                               LlIt(vector<int>(1, cocp.nx[0] + 1), vector<int>(1, cocp.nx[0]), 1),
                                               rhs_rq(sum(cocp.nu) + sum(cocp.nx), 1),
                                               rhs_b(sum(cocp.nx) - cocp.nx[0], 1),
                                               rhs_g(sum(cocp.ng), 1),
                                               rhs_rq2(sum(cocp.nu) + sum(cocp.nx), 1),
                                               rhs_b2(sum(cocp.nx) - cocp.nx[0], 1),
                                               rhs_g2(sum(cocp.ng), 1),
                                               ux_test(sum(cocp.nu) + sum(cocp.nx), 1),
                                               lam_test(sum(cocp.ng) + sum(cocp.nx) - cocp.nx[0], 1),
                                               v_Ppt(cocp.nx, cocp.K),
                                               v_Hh(cocp.nx, cocp.K),
                                               v_AL(vector<int>(1, max(cocp.nx)), 1),
                                               v_RSQrqt_tilde(cocp.nu + cocp.nx, cocp.K),
                                               v_Ggt_stripe(vector<int>(1, max(cocp.nx + cocp.nu)), 1),
                                               v_Ggt_tilde(cocp.nu + cocp.nx, cocp.K),
                                               v_GgLt(vector<int>(1, max(cocp.nu + cocp.nx)), 1),
                                               v_RSQrqt_hat(vector<int>(1, max(cocp.nx + cocp.nu)), 1),
                                               v_Llt(cocp.nu + cocp.nx, cocp.K),
                                               v_Llt_shift(vector<int>(1, max(cocp.nu)), 1),
                                               v_GgIt_tilde(vector<int>(1, cocp.nx[0]), 1),
                                               v_GgLIt(vector<int>(1, cocp.nx[0]), 1),
                                               v_HhIt(vector<int>(1, cocp.nx[0]), 1),
                                               v_PpIt_hat(vector<int>(1, cocp.nx[0]), 1),
                                               v_LlIt(vector<int>(1, cocp.nx[0]), 1),
                                               Pl(max(cocp.nu), cocp.K), // number of equations can never exceed nx
                                               Pr(max(cocp.nu), cocp.K),
                                               PlI(cocp.nx[0], 1),
                                               PrI(cocp.nx[0], 1),
                                               gamma(vector<int>(cocp.K, 0)),
                                               rho(vector<int>(cocp.K, 0)),
                                               ux_offs(offsets(cocp.nx + cocp.nu)),
                                               g_offs(offsets(cocp.ng)),
                                               dyn_offs(offsets(rotate(cocp.nx, 1))),
                                               dyn_eq_offs(offsets(rotate(cocp.nx, 1)) + sum(cocp.ng)),
                                               max_nu(max(cocp.nu)), max_nx(max(cocp.nx)),
                                               max_ng(max(cocp.ng)){};
int OCPLSRiccati::solve_pd_sys(
    COCP &OCP,
    const double inertia_correction_w,
    const double inertia_correction_c,
    const FatropVecBF &ux,
    const FatropVecBF &lam)
{
    lastused_.inertia_correction_c = inertia_correction_c;
    lastused_.inertia_correction_w = inertia_correction_w;
    if (inertia_correction_c == 0.0)
    {
        return solve_pd_sys_normal(OCP, inertia_correction_w, ux, lam);
    }
    else
    {
        return solve_pd_sys_degenerate(OCP, inertia_correction_w, inertia_correction_c, ux, lam);
    }
}
int OCPLSRiccati::solve_pd_sys_degenerate(
    COCP &OCP,
    const double inertia_correction_w,
    const double inertia_correction_c,
    const FatropVecBF &ux,
    const FatropVecBF &lam)
{
    // blasfeo_timer timer;
    // blasfeo_tic(&timer);
    // define compiler macros for notational convenience
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP.name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
    int K = OCP.K;
    // make variables local for efficiency
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    SOLVERMACRO(MAT *, Ppt, _p);
    SOLVERMACRO(MAT *, AL, _p);
    SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
    SOLVERMACRO(MAT *, Ggt_tilde, _p);
    SOLVERMACRO(MAT *, Llt, _p);
    SOLVERMACRO(MAT *, Llt_shift, _p);
    SOLVERMACRO(MAT *, LlIt, _p);
    SOLVERMACRO(VEC *, ux, _p);
    SOLVERMACRO(VEC *, lam, _p);
    int *nu_p = (int *)OCP.nu.data();
    int *nx_p = (int *)OCP.nx.data();
    int *ng_p = (int *)OCP.ng.data();
    MAT *RSQrq_hat_curr_p;
    double delta_cmin1 = 1 / inertia_correction_c;

    /////////////// recursion ///////////////

    // last stage
    {
        const int nx = nx_p[K - 1];
        const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
        const int ng = ng_p[K - 1];
        // const int offs_g_ineq_k = offs_g_ineq_p[K - 1];
        // Pp_Km1 <- Qq_Km1
        GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
        DIARE(nx, inertia_correction_w, Ppt_p + K - 1, 0, 0);
        GECP(nx + 1, ng, Ggt_p + K - 1, nu, 0, Ggt_tilde_p + K - 1, 0, 0); // needless operation because feature not implemented yet
        SYRK_LN_MN(nx + 1, nx, ng, delta_cmin1, Ggt_tilde_p + K - 1, 0, 0, Ggt_tilde_p + K - 1, 0, 0, 1.0, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
        TRTR_L(nx, Ppt_p + K - 1, 0, 0, Ppt_p + K - 1, 0, 0);
    }
    for (int k = K - 2; k >= 0; --k)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int nxp1 = nx_p[k + 1];
        const int ng = ng_p[k];
        // const int offs_g_ineq_k = offs_g_ineq_p[k];
        //////// SUBSDYN
        {
            // AL <- [BAb]^T_k P_kp1
            GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
            // AL[-1,:] <- AL[-1,:] + p_kp1^T
            GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
            // RSQrqt_stripe <- AL[BA] + RSQrqt
            SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
            // RSQrqt + 1/d_c At@A
            DIARE(nu + nx, inertia_correction_w, RSQrqt_tilde_p + k, 0, 0);
            SYRK_LN_MN(nu + nx + 1, nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, Ggt_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
        }
        //////// TRANSFORM_AND_SUBSEQ
        {
            RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
        }
        //////// SCHUR
        {
            // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
            POTRF_L_MN(nu + nx + 1, nu, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
            if (!check_reg(nu, Llt_p + k, 0, 0))
                return 1;
            // Pp_k = Qq_hatk - L_k^T @ Ll_k
            // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
            GECP(nx + 1, nu, Llt_p + k, nu, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
            SYRK_LN_MN(nx + 1, nx, nu, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu, nu, Ppt_p + k, 0, 0);
        }
        TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
    }
    //////// FIRST_STAGE
    {
        const int nx = nx_p[0];
        {
            POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
            if (!check_reg(nx, LlIt_p, 0, 0))
                return 2;
        }
    }
    ////// FORWARD_SUBSTITUTION:
    // first stage
    {
        const int nx = nx_p[0];
        const int nu = nu_p[0];
        // calculate xIb
        ROWEX(nx, -1.0, LlIt_p, nx, 0, ux_p, nu);
        // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        TRSV_LTN(nx, LlIt_p, 0, 0, ux_p, nu, ux_p, nu);
    }
    const int *offs_ux = ux_offs.data();
    const int *offs_g = g_offs.data();
    const int *offs_dyn_eq = dyn_eq_offs.data();
    // other stages
    // for (int k = 0; k < K - 1; k++)
    // int dyn_eqs_ofs = offs_g[K - 1] + ng_p[K - 1]; // this value is incremented at end of recursion
    for (int k = 0; k < K - 1; k++)
    {
        const int nx = nx_p[k];
        const int nu = nu_p[k];
        const int nxp1 = nx_p[k + 1];
        const int nup1 = nu_p[k + 1];
        const int offsp1 = offs_ux[k + 1];
        const int offs = offs_ux[k];
        const int offs_dyn_eq_k = offs_dyn_eq[k];
        /// calculate ukb_tilde
        // -Lkxk - lk
        ROWEX(nu, -1.0, Llt_p + k, nu + nx, 0, ux_p, offs);
        // assume aliasing of last two eliments is allowed
        GEMV_T(nx, nu, -1.0, Llt_p + k, nu, 0, ux_p, offs + nu, 1.0, ux_p, offs, ux_p, offs);
        TRSV_LTN(nu, Llt_p + k, 0, 0, ux_p, offs, ux_p, offs);
        // calculate xkp1
        ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_p, offsp1 + nup1);
        GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, 1.0, ux_p, offsp1 + nup1, ux_p, offsp1 + nup1);
        // calculate lam_dyn xp1
        ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
        GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
        // // calculate lam_eq xk
        // ROWEX(ng, -1.0, Ggt_p + k, 0, 0, lam_p, offs_g_k);
        // GEMV_T(nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
    }
    // calculate lam_eq xk
    for (int k = 0; k < K; k++)
    {
        const int nx = nx_p[k];
        const int nu = nu_p[k];
        const int ng = ng_p[k];
        const int offs = offs_ux[k];
        const int offs_g_k = offs_g[k];
        ROWEX(ng, delta_cmin1, Ggt_p + k, nu + nx, 0, lam_p, offs_g_k);
        GEMV_T(nu + nx, ng, delta_cmin1, Ggt_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
    }
    return 0;
}
int OCPLSRiccati::solve_pd_sys_normal(
    COCP &OCP,
    const double inertia_correction,
    const FatropVecBF &ux,
    const FatropVecBF &lam)
{
    bool increased_accuracy = true;
    // bool it_ref = true;
    // blasfeo_timer timer;
    // blasfeo_tic(&timer);
    // define compiler macros for notational convenience
#define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP.name)
#define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
#define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
    int K = OCP.K;
    // make variables local for efficiency
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    SOLVERMACRO(MAT *, Ppt, _p);
    SOLVERMACRO(MAT *, Hh, _p);
    SOLVERMACRO(MAT *, AL, _p);
    SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
    SOLVERMACRO(MAT *, Ggt_stripe, _p);
    SOLVERMACRO(MAT *, Ggt_tilde, _p);
    SOLVERMACRO(PMAT *, Pl, _p);
    SOLVERMACRO(PMAT *, Pr, _p);
    SOLVERMACRO(MAT *, GgLt, _p);
    SOLVERMACRO(MAT *, RSQrqt_hat, _p);
    SOLVERMACRO(MAT *, Llt, _p);
    SOLVERMACRO(MAT *, Llt_shift, _p);
    SOLVERMACRO(MAT *, GgIt_tilde, _p);
    SOLVERMACRO(MAT *, GgLIt, _p);
    SOLVERMACRO(MAT *, HhIt, _p);
    SOLVERMACRO(MAT *, PpIt_hat, _p);
    SOLVERMACRO(MAT *, LlIt, _p);
    SOLVERMACRO(PMAT *, PlI, _p);
    SOLVERMACRO(PMAT *, PrI, _p);
    SOLVERMACRO(VEC *, ux, _p);
    SOLVERMACRO(VEC *, lam, _p);
    int *nu_p = (int *)OCP.nu.data();
    int *nx_p = (int *)OCP.nx.data();
    int *ng_p = (int *)OCP.ng.data();
    int *gamma_p = (int *)gamma.data();
    int *rho_p = (int *)rho.data();
    MAT *RSQrq_hat_curr_p;
    int rank_k;

    /////////////// recursion ///////////////

    // last stage
    {
        const int nx = nx_p[K - 1];
        const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
        const int ng = ng_p[K - 1];
        // const int offs_g_ineq_k = offs_g_ineq_p[K - 1];
        // Pp_Km1 <- Qq_Km1
        GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
        DIARE(nx, inertia_correction, Ppt_p + K - 1, 0, 0);
        // Hh_Km1 <- Gg_Km1
        GETR(nx + 1, ng, Ggt_p + (K - 1), nu, 0, Hh_p + (K - 1), 0, 0);
        gamma_p[K - 1] = ng;
        rho_p[K - 1] = 0;
    }
    for (int k = K - 2; k >= 0; --k)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int nxp1 = nx_p[k + 1];
        const int ng = ng_p[k];
        // const int offs_g_k = offs_g_p[k];
        // const int offs_g_ineq_k = offs_g_ineq_p[k];
        // calculate the size of H_{k+1} matrix
        const int Hp1_size = gamma_p[k + 1] - rho_p[k + 1];
        if (Hp1_size + ng > nu + nx)
            return -1;
        // gamma_k <- number of eqs represented by Ggt_stripe
        const int gamma_k = Hp1_size + ng;
        // if(k==0) blasfeo_print_dmat(1, nxp1, Ppt_p+k+1, nx, 0);
        //////// SUBSDYN
        {
            // AL <- [BAb]^T_k P_kp1
            GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
            // AL[-1,:] <- AL[-1,:] + p_kp1^T
            GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
            // RSQrqt_stripe <- AL[BA] + RSQrqt
            SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
            // if(k==K-2) blasfeo_print_dmat(1, nu+nx, RSQrqt_tilde_p+k, nu+nx, 0);
            DIARE(nu + nx, inertia_correction, RSQrqt_tilde_p + k, 0, 0);
            gamma_p[k] = gamma_k;
            // if ng[k]>0
            if (gamma_k > 0)
            {
                // if Gk nonempty
                if (ng > 0)
                {
                    // Ggt_stripe  <- Ggt_k
                    GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_stripe_p, 0, 0);
                }
                // if Hkp1 nonempty
                if (Hp1_size > 0)
                {
                    // Ggt_stripe <- [Ggt_k [BAb_k^T]H_kp1]
                    GEMM_NT(nu + nx + 1, Hp1_size, nxp1, 1.0, BAbt_p + k, 0, 0, Hh_p + (k + 1), 0, 0, 0.0, Ggt_stripe_p, 0, ng, Ggt_stripe_p, 0, ng);
                    // Ggt_stripe[-1,ng:] <- Ggt_stripe[-1,ng:] + h_kp1^T
                    GEADTR(1, Hp1_size, 1.0, Hh_p + (k + 1), 0, nxp1, Ggt_stripe_p, nu + nx, ng);
                }
            }
            else
            {
                rho_p[k] = 0;
                rank_k = 0;
                RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
            }
        }
        //////// TRANSFORM_AND_SUBSEQ
        {
            // symmetric transformation, done a little different than in paper, in order to fuse LA operations
            // LU_FACT_TRANSPOSE(Ggtstripe[:gamma_k, nu+nx+1], nu max)
            // if(k==K-2) blasfeo_print_dmat(1, gamma_k, Ggt_stripe_p, nu+nx, 0);
            LU_FACT_transposed(gamma_k, nu + nx + 1, nu, rank_k, Ggt_stripe_p, Pl_p + k, Pr_p + k);

            rho_p[k] = rank_k;
            if (gamma_k - rank_k > 0)
            {
                // transfer eq's to next stage
                GETR(nx + 1, gamma_k - rank_k, Ggt_stripe_p, nu, rank_k, Hh_p + k, 0, 0);
            }
            if (rank_k > 0)
            {
                // Ggt_tilde_k <- Ggt_stripe[rho_k:nu+nx+1, :rho] L-T (note that this is slightly different from the implementation)
                TRSM_RLNN(nu - rank_k + nx + 1, rank_k, -1.0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, rank_k, 0, Ggt_tilde_p + k, 0, 0);
                // the following command copies the top block matrix (LU) to the bottom because it it needed later
                GECP(rank_k, gamma_k, Ggt_stripe_p, 0, 0, Ggt_tilde_p + k, nu - rank_k + nx + 1, 0);
                // permutations
                TRTR_L(nu + nx, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0); // copy lower part of RSQ to upper part
                (Pr_p + k)->PM(rank_k, RSQrqt_tilde_p + k);                          // TODO make use of symmetry
                (Pr_p + k)->MPt(rank_k, RSQrqt_tilde_p + k);
                // GL <- Ggt_tilde_k @ RSQ[:rho,:nu+nx] + RSQrqt[rho:nu+nx+1, rho:] (with RSQ[:rho,:nu+nx] = RSQrqt[:nu+nx,:rho]^T)
                // GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
                // GgLt[0].print();
                GECP(nu - rank_k + nx + 1, nu + nx, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, 0, GgLt_p, 0, 0);
                // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                SYRK_LN_MN(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                // GEMM_NT(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                RSQrq_hat_curr_p = RSQrqt_hat_p;
            }
            else
            {
                RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
            }
            // if(k==K-2) blasfeo_print_dmat(1, nu+nx-rank_k, RSQrq_hat_curr_p, nu+nx -rank_k, 0);
        }
        //////// SCHUR
        {
            if (nu - rank_k > 0)
            {
                // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
                POTRF_L_MN(nu - rank_k + nx + 1, nu - rank_k, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
                if (!check_reg(nu - rank_k, Llt_p + k, 0, 0))
                    return 1;
                // Pp_k = Qq_hatk - L_k^T @ Ll_k
                // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
                GECP(nx + 1, nu - rank_k, Llt_p + k, nu - rank_k, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
                // SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                GECP(nx + 1, nx, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
                // next steps are for better accuracy
                if (increased_accuracy)
                {
                    // copy eta
                    GETR(nu - rank_k, gamma_k - rank_k, Ggt_stripe_p, rank_k, rank_k, Ggt_stripe_p, 0, 0);
                    // blasfeo_print_dmat(gamma_k-rank_k, nu-rank_k, Ggt_stripe_p, 0,0);
                    // eta L^-T
                    TRSM_RLTN(gamma_k - rank_k, nu - rank_k, 1.0, Llt_p + k, 0, 0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, 0, 0);
                    // ([S^T \\ r^T] L^-T) @ (L^-1 eta^T)
                    // (eta L^-T) @ ([S^T \\ r^T] L^-T)^T
                    GEMM_NT(gamma_k - rank_k, nx + 1, nu - rank_k, -1.0, Ggt_stripe_p, 0, 0, Llt_p + k, nu - rank_k, 0, 1.0, Hh_p + k, 0, 0, Hh_p + k, 0, 0);
                    // keep (L^-1 eta^T) for forward recursion
                    GETR(gamma_k - rank_k, nu - rank_k, Ggt_stripe_p, 0, 0, Ggt_tilde_p + k, 0, rank_k);
                }
            }
            else
            {
                GECP(nx + 1, nx, RSQrq_hat_curr_p, 0, 0, Ppt_p + k, 0, 0);
            }
            TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
        }
    }
    rankI = 0;
    //////// FIRST_STAGE
    {
        const int nx = nx_p[0];
        int gamma_I = gamma_p[0] - rho_p[0];
        if (gamma_I > nx)
        {
            return -3;
        }
        if (gamma_I > 0)
        {
            GETR(gamma_I, nx + 1, Hh_p + 0, 0, 0, HhIt_p, 0, 0); // transposition may be avoided
            // HhIt[0].print();
            LU_FACT_transposed(gamma_I, nx + 1, nx, rankI, HhIt_p, PlI_p, PrI_p);
            if (rankI < gamma_I)
                return -2;
            // PpIt_tilde <- Ggt[rankI:nx+1, :rankI] L-T (note that this is slightly different from the implementation)
            TRSM_RLNN(nx - rankI + 1, rankI, -1.0, HhIt_p, 0, 0, HhIt_p, rankI, 0, GgIt_tilde_p, 0, 0);
            // permutations
            (PrI_p)->PM(rankI, Ppt_p); // TODO make use of symmetry
            (PrI_p)->MPt(rankI, Ppt_p);
            // // GL <- GgIt_tilde @ Pp[:rankI,:nx] + Ppt[rankI:nx+1, rankI:] (with Pp[:rankI,:nx] = Ppt[:nx,:rankI]^T)
            // GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
            // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
            GECP(nx - rankI + 1, nx, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
            GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, GgLIt_p, 0, 0, GgLIt_p, 0, 0);
            // // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
            SYRK_LN_MN(nx - rankI + 1, nx - rankI, rankI, 1.0, GgLIt_p, 0, 0, GgIt_tilde_p, 0, 0, 1.0, GgLIt_p, 0, rankI, PpIt_hat_p, 0, 0);
            // TODO skipped if nx-rankI = 0
            POTRF_L_MN(nx - rankI + 1, nx - rankI, PpIt_hat_p, 0, 0, LlIt_p, 0, 0);
            if (!check_reg(nx - rankI, LlIt_p, 0, 0))
                return 2;
        }
        else
        {
            rankI = 0;
            POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
            if (!check_reg(nx, LlIt_p, 0, 0))
                return 2;
        }
    }
    ////// FORWARD_SUBSTITUTION:
    // first stage
    {
        const int nx = nx_p[0];
        const int nu = nu_p[0];
        // calculate xIb
        ROWEX(nx - rankI, -1.0, LlIt_p, nx - rankI, 0, ux_p, nu + rankI);
        // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        TRSV_LTN(nx - rankI, LlIt_p, 0, 0, ux_p, nu + rankI, ux_p, nu + rankI);
        // calculate xIa
        ROWEX(rankI, 1.0, GgIt_tilde_p, nx - rankI, 0, ux_p, nu);
        // assume aliasing is possible for last two elements
        GEMV_T(nx - rankI, rankI, 1.0, GgIt_tilde_p, 0, 0, ux_p, nu + rankI, 1.0, ux_p, nu, ux_p, nu);
        //// lag
        ROWEX(rankI, -1.0, Ppt_p, nx, 0, lam_p, 0);
        // assume aliasing is possible for last two elements
        GEMV_T(nx, rankI, -1.0, Ppt_p, 0, 0, ux_p, nu, 1.0, lam_p, 0, lam_p, 0);
        // U^-T
        TRSV_LNN(rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
        // L^-T
        TRSV_UNU(rankI, rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
        (PlI_p)->PtV(rankI, lam_p, 0);
        (PrI_p)->PtV(rankI, ux_p, nu);
        // blasfeo_print_dvec(rankI, lam_p, 0);
    }
    int *offs_ux = (int *)ux_offs.data();
    int *offs_g = (int *)g_offs.data();
    int *offs_dyn_eq = (int *)dyn_eq_offs.data();
    // other stages
    // for (int k = 0; k < K - 1; k++)
    // int dyn_eqs_ofs = offs_g[K - 1] + ng_p[K - 1]; // this value is incremented at end of recursion
    for (int k = 0; k < K - 1; k++)
    {
        const int nx = nx_p[k];
        const int nu = nu_p[k];
        const int nxp1 = nx_p[k + 1];
        const int nup1 = nu_p[k + 1];
        const int offsp1 = offs_ux[k + 1];
        const int offs = offs_ux[k];
        const int rho_k = rho_p[k];
        const int numrho_k = nu - rho_k;
        const int offs_g_k = offs_g[k];
        const int offs_dyn_eq_k = offs_dyn_eq[k];
        const int offs_g_kp1 = offs_g[k + 1];
        const int gammamrho_k = gamma_p[k] - rho_p[k];
        const int gamma_k = gamma_p[k];
        const int gammamrho_kp1 = gamma_p[k + 1] - rho_p[k + 1];
        if (numrho_k > 0)
        {
            /// calculate ukb_tilde
            // -Lkxk - lk
            ROWEX(numrho_k, -1.0, Llt_p + k, numrho_k + nx, 0, ux_p, offs + rho_k);
            if (increased_accuracy)
            {
                GEMV_N(nu - rho_k, gamma_k - rho_k, -1.0, Ggt_tilde_p + k, 0, rho_k, lam_p, offs_g_k, 1.0, ux_p, offs + rho_k, ux_p, offs + rho_k);
            }
            // assume aliasing of last two eliments is allowed
            GEMV_T(nx, numrho_k, -1.0, Llt_p + k, numrho_k, 0, ux_p, offs + nu, 1.0, ux_p, offs + rho_k, ux_p, offs + rho_k);
            TRSV_LTN(numrho_k, Llt_p + k, 0, 0, ux_p, offs + rho_k, ux_p, offs + rho_k);
        }
        /// calcualate uka_tilde
        if (rho_k > 0)
        {
            ROWEX(rho_k, 1.0, Ggt_tilde_p + k, numrho_k + nx, 0, ux_p, offs);
            GEMV_T(nx + numrho_k, rho_k, 1.0, Ggt_tilde_p + k, 0, 0, ux_p, offs + rho_k, 1.0, ux_p, offs, ux_p, offs);
            // calculate lamda_tilde_k
            // copy vk to right location
            // we implemented a version of vector copy that starts with copy of last element, to avoid aliasing error
            VECCPR(gammamrho_k, lam_p, offs_g_k, lam_p, offs_g_k + rho_k);
            ROWEX(rho_k, -1.0, RSQrqt_tilde_p + k, nu + nx, 0, lam_p, offs_g_k);
            // assume aliasing of last two eliments is allowed
            GEMV_T(nu + nx, rho_k, -1.0, RSQrqt_tilde_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
            // nu-rank_k+nx,0
            // needless copy because feature not implemented yet in trsv_lnn
            GECP(rho_k, gamma_k, Ggt_tilde_p + k, nu - rho_k + nx + 1, 0, AL_p, 0, 0);
            // U^-T
            TRSV_LNN(rho_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
            // L^-T
            TRSV_UNU(rho_k, gamma_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
            (Pl_p + k)->PtV(rho_k, lam_p, offs_g_k);
            (Pr_p + k)->PtV(rho_k, ux_p, offs);
        }
        // calculate xkp1
        ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_p, offsp1 + nup1);
        GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, 1.0, ux_p, offsp1 + nup1, ux_p, offsp1 + nup1);
        // calculate lam_dyn xp1
        ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
        GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
        GEMV_T(gammamrho_kp1, nxp1, 1.0, Hh_p + (k + 1), 0, 0, lam_p, offs_g_kp1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
    }
    // double el = blasfeo_toc(&timer);
    // cout << "el time fact solve" << el << endl; //
    if (it_ref)
    {
        const int min_it_ref = 1;
        double err_curr = 0.0;
        // copy(ux, ux_test[0]);
        // copy(lam, lam_test[0]);
        // copy(delta_s, delta_s_test[0]);
        // blasfeo_tic(&timer);
        get_rhs(
            OCP,
            rhs_rq2[0],
            rhs_b2[0],
            rhs_g2[0]);
        // el = blasfeo_toc(&timer);
        // cout << "el time get rhs" << el << endl; //
        double max_norm = std::max(Linf(rhs_g2[0]), std::max(Linf(rhs_rq2[0]), Linf(rhs_b2[0])));
        max_norm = (max_norm == 0.0) ? 1.0 : max_norm;
        double error_prev = -1.0;
        for (int i = 0; i < 5; i++)
        {
            // blasfeo_tic(&timer);
            compute_pd_sys_times_vec(
                OCP,
                inertia_correction,
                0.0,
                ux,
                lam,
                rhs_rq[0],
                rhs_b[0],
                rhs_g[0]);
            // el = blasfeo_toc(&timer);
            // cout << "el time compute mv prod" << el << endl; //
            axpby(1.0, rhs_rq[0], 1.0, rhs_rq2[0], rhs_rq[0]);
            axpby(1.0, rhs_b[0], 1.0, rhs_b2[0], rhs_b[0]);
            axpby(1.0, rhs_g[0], 1.0, rhs_g2[0], rhs_g[0]);

            // cout << "residu rq:  " << Linf(rhs_rq[0]) / max_norm << "  ";
            // cout << "residu b:  " << Linf(rhs_b[0]) / max_norm << "  ";
            // cout << "residu g:  " << Linf(rhs_g[0]) / max_norm << "  ";
            // cout << "residu g_ineq:  " << Linf(rhs_g_ineq[0]) / max_norm << "  ";
            // cout << "residu gradb:  " << Linf(rhs_gradb[0]) / max_norm  << "  "<<endl;
            err_curr = std::max(Linf(rhs_g[0]), std::max(Linf(rhs_rq[0]), Linf(rhs_b[0]))) / max_norm;
            // cout << "residu:  " << err_curr << endl;
            if (i >= min_it_ref)
            {
                if (err_curr < 1e-8 || (error_prev > 0.0 && err_curr > 1.0 * error_prev))
                {
                    if (err_curr > 1e-8)
                    {
                        // cout << "stopped it_ref because insufficient decrease err_curr:  " << err_curr << endl;
                    }
                    return 0;
                }
            }
            // blasfeo_tic(&timer);
            solve_rhs(
                OCP,
                ux_test[0],
                lam_test[0],
                rhs_rq[0],
                rhs_b[0],
                rhs_g[0]);
            // el = blasfeo_toc(&timer);
            // cout << "el time solveRHS " << el << endl; //
            // el = blasfeo_toc(&timer);
            // cout << "el time solveRHS " << el << endl;
            axpby(1.0, ux_test[0], 1.0, ux, ux);
            axpby(1.0, lam_test[0], 1.0, lam, lam);

            // prepare next iteration
            error_prev = err_curr;
        }
    }
    return 0;
}
int OCPLSRiccati::get_rhs(
    COCP &OCP,
    const FatropVecBF &rhs_rq,
    const FatropVecBF &rhs_b,
    const FatropVecBF &rhs_g)
{
    int K = OCP.K;
    // make variables local for efficiency
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    SOLVERMACRO(VEC *, rhs_rq, _p);
    SOLVERMACRO(VEC *, rhs_b, _p);
    SOLVERMACRO(VEC *, rhs_g, _p);
    int *nu_p = (int *)OCP.nu.data();
    int *nx_p = (int *)OCP.nx.data();
    int *ng_p = (int *)OCP.ng.data();
    int *offs_ux = (int *)ux_offs.data();
    int *offs_g = (int *)g_offs.data();
    //////////////////////////////
    ////////////// rhs_rq
    //////////////////////////////
    for (int k = 0; k < K; k++)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int offs = offs_ux[k];
        ROWEX(nu + nx, 1.0, RSQrqt_p + k, nu + nx, 0, rhs_rq_p, offs);
    }
    //////////////////////////////
    ////////////// rhs_b
    //////////////////////////////
    {
        int offs_b = 0;
        for (int k = 0; k < K - 1; k++)
        {
            const int nu = nu_p[k];
            const int nx = nx_p[k];
            const int nxp1 = nx_p[k + 1];
            ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, rhs_b_p, offs_b);
            offs_b += nxp1;
        }
    }
    //////////////////////////////
    ////////////// rhs_g
    //////////////////////////////
    for (int k = 0; k < K; k++)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int ng = ng_p[k];
        const int offs_g_k = offs_g[k];
        ROWEX(ng, 1.0, Ggt_p + k, nu + nx, 0, rhs_g_p, offs_g_k);
    }
    return 0;
}
int OCPLSRiccati::compute_pd_sys_times_vec(
    COCP &OCP,
    const double inertia_correction_w,
    const double inertia_correction_c,
    const FatropVecBF &ux,
    const FatropVecBF &lam,
    const FatropVecBF &rhs_rq,
    const FatropVecBF &rhs_b,
    const FatropVecBF &rhs_g)
{
    int K = OCP.K;
    // make variables local for efficiency
    OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(MAT *, BAbt, _p);
    OCPMACRO(MAT *, Ggt, _p);
    SOLVERMACRO(VEC *, ux, _p);
    SOLVERMACRO(VEC *, lam, _p);
    SOLVERMACRO(VEC *, rhs_rq, _p);
    SOLVERMACRO(VEC *, rhs_b, _p);
    SOLVERMACRO(VEC *, rhs_g, _p);
    int *nu_p = (int *)OCP.nu.data();
    int *nx_p = (int *)OCP.nx.data();
    int *ng_p = (int *)OCP.ng.data();
    //////////////////////////////
    ////////////// rhs_rq
    //////////////////////////////
    int *offs_ux = (int *)ux_offs.data();
    int *offs_g = (int *)g_offs.data();
    int *offs_dyn_eq = (int *)dyn_eq_offs.data();
    // contribution of Hessian
    for (int k = 0; k < K; k++)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int offs = offs_ux[k];
        GEMV_N(nu + nx, nu + nx, 1.0, RSQrqt_p + k, 0, 0, ux_p, offs, 0.0, rhs_rq_p, offs, rhs_rq_p, offs);
        // contribution of inertia correction
        AXPY(nu + nx, inertia_correction_w, ux_p, offs, rhs_rq_p, offs, rhs_rq_p, offs);
    }
    // contribution of dynamics constraints
    for (int k = 0; k < K - 1; k++)
    {
        const int nu = nu_p[k];
        const int nup1 = nu_p[k + 1];
        const int nx = nx_p[k];
        const int nxp1 = nx_p[k + 1];
        const int offsp1 = offs_ux[k + 1];
        const int offs = offs_ux[k];
        const int offs_dyn_eq_k = offs_dyn_eq[k];
        GEMV_N(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, lam_p, offs_dyn_eq_k, 1.0, rhs_rq_p, offs, rhs_rq_p, offs);
        AXPY(nxp1, -1.0, lam_p, offs_dyn_eq_k, rhs_rq_p, offsp1 + nup1, rhs_rq_p, offsp1 + nup1);
    }
    // contribution of equality constraints
    for (int k = 0; k < K; k++)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int ng = ng_p[k];
        const int offs = offs_ux[k];
        const int offs_g_k = offs_g[k];
        GEMV_N(nu + nx, ng, 1.0, Ggt_p + k, 0, 0, lam_p, offs_g_k, 1.0, rhs_rq_p, offs, rhs_rq_p, offs);
    }

    //////////////////////////////
    ////////////// rhs_b
    //////////////////////////////
    {
        int offs_b = 0;
        for (int k = 0; k < K - 1; k++)
        {
            const int nu = nu_p[k];
            const int nup1 = nu_p[k + 1];
            const int nx = nx_p[k];
            const int nxp1 = nx_p[k + 1];
            const int offs = offs_ux[k];
            const int offsp1 = offs_ux[k + 1];
            VECCP(nxp1, ux_p, offsp1 + nup1, rhs_b_p, offs_b);
            GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, -1.0, rhs_b_p, offs_b, rhs_b_p, offs_b);
            offs_b += nxp1;
        }
    }

    //////////////////////////////
    ////////////// rhs_g
    //////////////////////////////
    // contribution of equality constraints
    for (int k = 0; k < K; k++)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int ng = ng_p[k];
        const int offs = offs_ux[k];
        const int offs_g_k = offs_g[k];
        GEMV_T(nu + nx, ng, 1.0, Ggt_p + k, 0, 0, ux_p, offs, 0.0, rhs_g_p, offs_g_k, rhs_g_p, offs_g_k);
    }
    if (inertia_correction_c != 0.0)
    {
        axpy(-inertia_correction_c, lam.block(offs_g[0], rhs_g.nels()), rhs_g, rhs_g);
    }
    return 0;
};
int OCPLSRiccati::solve_rhs(
    COCP &OCP,
    const FatropVecBF &ux,
    const FatropVecBF &lam,
    const FatropVecBF &rhs_rq,
    const FatropVecBF &rhs_b,
    const FatropVecBF &rhs_g)
{
    // double inertia_correction_w = lastused_.inertia_correction_w;
    double inertia_correction_c = lastused_.inertia_correction_c;
    assert(inertia_correction_c == 0.0); // not implemented yet
    bool increased_accuracy = true;
    //     // blasfeo_timer timer;
    //     // blasfeo_tic(&timer);
    //     // define compiler macros for notational convenience
    // #define OCPMACRO(type, name, suffix) type name##suffix = ((type)OCP->name)
    // #define AUXMACRO(type, name, suffix) type name##suffix = ((type)OCP->aux.name)
    // #define SOLVERMACRO(type, name, suffix) type name##suffix = ((type)name)
    int K = OCP.K;
    //     // make variables local for efficiency
    // OCPMACRO(MAT *, RSQrqt, _p);
    OCPMACRO(MAT *, BAbt, _p);
    // OCPMACRO(MAT *, Ggt, _p);
    SOLVERMACRO(MAT *, Ppt, _p);
    SOLVERMACRO(MAT *, Hh, _p);
    SOLVERMACRO(MAT *, AL, _p);
    SOLVERMACRO(MAT *, RSQrqt_tilde, _p);
    SOLVERMACRO(MAT *, Ggt_stripe, _p);
    SOLVERMACRO(MAT *, Ggt_tilde, _p);
    SOLVERMACRO(PMAT *, Pl, _p);
    SOLVERMACRO(PMAT *, Pr, _p);
    // SOLVERMACRO(MAT *, GgLt, _p);
    // SOLVERMACRO(MAT *, RSQrqt_hat, _p);
    SOLVERMACRO(MAT *, Llt, _p);
    SOLVERMACRO(MAT *, Llt_shift, _p);
    SOLVERMACRO(MAT *, GgIt_tilde, _p);
    // SOLVERMACRO(MAT *, GgLIt, _p);
    SOLVERMACRO(MAT *, HhIt, _p);
    // SOLVERMACRO(MAT *, PpIt_hat, _p);
    SOLVERMACRO(MAT *, LlIt, _p);
    // SOLVERMACRO(MAT *, Ggt_ineq_temp, _p);
    SOLVERMACRO(PMAT *, PlI, _p);
    SOLVERMACRO(PMAT *, PrI, _p);
    SOLVERMACRO(VEC *, ux, _p);
    SOLVERMACRO(VEC *, lam, _p);
    SOLVERMACRO(VEC *, rhs_rq, _p);
    SOLVERMACRO(VEC *, rhs_b, _p);
    SOLVERMACRO(VEC *, rhs_g, _p);
    VEC *v_RSQrq_hat_curr_p;
    // MAT *RSQrq_hat_curr_p;
    int rank_k;
    int *offs_ux = (int *)ux_offs.data();
    int *offs_g = (int *)g_offs.data();
    int *offs_dyn_eq = (int *)dyn_eq_offs.data();
    int *offs_dyn = (int *)dyn_offs.data();
    int *nu_p = (int *)OCP.nu.data();
    int *nx_p = (int *)OCP.nx.data();
    int *ng_p = (int *)OCP.ng.data();
    int *gamma_p = (int *)gamma.data();
    int *rho_p = (int *)rho.data();

    // todo make this member variables

    VEC *v_Ppt_p = (VEC *)v_Ppt[0];
    VEC *v_Hh_p = (VEC *)v_Hh[0];
    VEC *v_AL_p = (VEC *)v_AL[0];
    VEC *v_RSQrqt_tilde_p = (VEC *)v_RSQrqt_tilde[0];
    VEC *v_Ggt_stripe_p = (VEC *)v_Ggt_stripe[0];
    VEC *v_Ggt_tilde_p = (VEC *)v_Ggt_tilde[0];
    VEC *v_GgLt_p = (VEC *)v_GgLt[0];
    VEC *v_RSQrqt_hat_p = (VEC *)v_RSQrqt_hat[0];
    VEC *v_Llt_p = (VEC *)v_Llt[0];
    VEC *v_Llt_shift_p = (VEC *)v_Llt_shift[0];
    VEC *v_GgIt_tilde_p = (VEC *)v_GgIt_tilde[0];
    VEC *v_GgLIt_p = (VEC *)v_GgLIt[0];
    VEC *v_HhIt_p = (VEC *)v_HhIt[0];
    VEC *v_PpIt_hat_p = (VEC *)v_PpIt_hat[0];
    VEC *v_LlIt_p = (VEC *)v_LlIt[0];

    /////////////// recursion ///////////////

    // last stage
    {
        const int nx = nx_p[K - 1];
        const int nu = nu_p[K - 1]; // this should be zero but is included here in case of misuse
        const int ng = ng_p[K - 1];
        const int offs_g_k = offs_g[K - 1];
        const int offs = offs_ux[K - 1];
        //         GECP(nx + 1, nx, RSQrqt_p + (K - 1), nu, nu, Ppt_p + K - 1, 0, 0);
        VECCP(nx, rhs_rq_p, offs + nu, v_Ppt_p + K - 1, 0);
        //         DIARE(nx, inertia_correction, Ppt_p + K - 1, 0, 0);
        //         //// inequalities
        //         // Hh_Km1 <- Gg_Km1
        //         GETR(nx + 1, ng, Ggt_p + (K - 1), nu, 0, Hh_p + (K - 1), 0, 0);
        VECCP(ng, rhs_g_p, offs_g_k, v_Hh_p + (K - 1), 0);
        //         gamma_p[K - 1] = ng;
        //         rho_p[K - 1] = 0;
    }
    for (int k = K - 2; k >= 0; --k)
    {
        const int nu = nu_p[k];
        const int nx = nx_p[k];
        const int nxp1 = nx_p[k + 1];
        const int ng = ng_p[k];
        const int offs_dyn_k = offs_dyn[k];
        const int offs_g_k = offs_g[k];
        const int offs = offs_ux[k];
        // calculate the size of H_{k+1} matrix
        const int Hp1_size = gamma_p[k + 1] - rho_p[k + 1];
        // if (Hp1_size + ng > nu + nx)
        //     return -1;
        // gamma_k <- number of eqs represented by Ggt_stripe
        const int gamma_k = Hp1_size + ng;
        // if (k==0) blasfeo_print_dvec(nxp1, v_Ppt_p+k+1, 0);
        //         //////// SUBSDYN
        {
            //             // AL <- [BAb]^T_k P_kp1
            //             GEMM_NT(nu + nx + 1, nxp1, nxp1, 1.0, BAbt_p + k, 0, 0, Ppt_p + k + 1, 0, 0, 0.0, AL_p, 0, 0, AL_p, 0, 0);
            GEMV_N(nxp1, nxp1, 1.0, Ppt_p + k + 1, 0, 0, rhs_b_p, offs_dyn_k, 0.0, v_AL_p, 0, v_AL_p, 0);
            //             // AL[-1,:] <- AL[-1,:] + p_kp1^T
            //             GEAD(1, nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, AL_p, nx + nu, 0);
            AXPY(nxp1, 1.0, v_Ppt_p + (k + 1), 0, v_AL_p, 0, v_AL_p, 0);
            // if (k==K-2) blasfeo_print_dvec(nxp1, v_AL_p, 0);
            //             // RSQrqt_stripe <- AL[BA] + RSQrqt
            //             SYRK_LN_MN(nu + nx + 1, nu + nx, nxp1, 1.0, AL_p, 0, 0, BAbt_p + k, 0, 0, 1.0, RSQrqt_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0);
            GEMV_N(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, v_AL_p, 0, 1.0, rhs_rq_p, offs, v_RSQrqt_tilde_p + k, 0);
            // if (k==K-2) blasfeo_print_dvec(nu+nx, v_RSQrqt_tilde_p+k, 0);
            //             //// inequalities
            //             DIARE(nu + nx, inertia_correction, RSQrqt_tilde_p + k, 0, 0);
            //             gamma_p[k] = gamma_k;
            //             // if ng[k]>0
            if (gamma_k > 0)
            {
                //                 // if Gk nonempty
                if (ng > 0)
                {
                    //                     // Ggt_stripe  <- Ggt_k
                    //                     GECP(nu + nx + 1, ng, Ggt_p + k, 0, 0, Ggt_stripe_p, 0, 0);
                    VECCP(ng, rhs_g_p, offs_g_k, v_Ggt_stripe_p, 0);
                }
                //                 // if Hkp1 nonempty
                if (Hp1_size > 0)
                {
                    //                     // Ggt_stripe <- [Ggt_k [BAb_k^T]H_kp1]
                    //                     GEMM_NT(nu + nx + 1, Hp1_size, nxp1, 1.0, BAbt_p + k, 0, 0, Hh_p + (k + 1), 0, 0, 0.0, Ggt_stripe_p, 0, ng, Ggt_stripe_p, 0, ng);
                    GEMV_N(Hp1_size, nxp1, 1.0, Hh_p + (k + 1), 0, 0, rhs_b_p, offs_dyn_k, 0.0, v_Ggt_stripe_p, ng, v_Ggt_stripe_p, ng);
                    //                     // Ggt_stripe[-1,ng:] <- Ggt_stripe[-1,ng:] + h_kp1^T
                    //                     GEADTR(1, Hp1_size, 1.0, Hh_p + (k + 1), 0, nxp1, Ggt_stripe_p, nu + nx, ng);
                    AXPY(Hp1_size, 1.0, v_Hh_p + (k + 1), 0, v_Ggt_stripe_p, ng, v_Ggt_stripe_p, ng);
                }
            }
            else
            {
                //                 rho_p[k] = 0;
                rank_k = 0;
                v_RSQrq_hat_curr_p = v_RSQrqt_tilde_p + k;
            }
        }
        //         //////// TRANSFORM_AND_SUBSEQ
        {
            //             // symmetric transformation, done a little different than in paper, in order to fuse LA operations
            //             // LU_FACT_TRANSPOSE(Ggtstripe[:gamma_k, nu+nx+1], nu max)
            //             LU_FACT_transposed(gamma_k, nu + nx + 1, nu, rank_k, Ggt_stripe_p, Pl_p + k, Pr_p + k);
            // TODO GET RID OF THIS OPERATION
            rank_k = rho_p[k];
            // if (k==K-2) blasfeo_print_dvec(gamma_k, v_Ggt_stripe_p, 0);
            GECP(rank_k, gamma_k, Ggt_tilde_p + k, nu - rank_k + nx + 1, 0, Ggt_stripe_p, 0, 0);
            (Pl_p + k)->PV(rank_k, v_Ggt_stripe_p, 0);
            // L1^-1 g_stipe[:rho]
            TRSV_UTU(rank_k, Ggt_stripe_p, 0, 0, v_Ggt_stripe_p, 0, v_Ggt_stripe_p, 0);
            // -L2 L1^-1 g_stripe[:rho] + g_stripe[rho:]
            GEMV_T(rank_k, gamma_k - rank_k, -1.0, Ggt_stripe_p, 0, rank_k, v_Ggt_stripe_p, 0, 1.0, v_Ggt_stripe_p, rank_k, v_Ggt_stripe_p, rank_k);

            //             rho_p[k] = rank_k;
            if (gamma_k - rank_k > 0)
            {
                //                 // transfer eq's to next stage
                //                 GETR(nx + 1, gamma_k - rank_k, Ggt_stripe_p, nu, rank_k, Hh_p + k, 0, 0);
                VECCP(gamma_k - rank_k, v_Ggt_stripe_p, rank_k, v_Hh_p + k, 0);
            }
            if (rank_k > 0)
            {
                //                 // Ggt_tilde_k <- Ggt_stripe[rho_k:nu+nx+1, :rho] L-T (note that this is slightly different from the implementation)
                //                 TRSM_RLNN(nu - rank_k + nx + 1, rank_k, -1.0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, rank_k, 0, Ggt_tilde_p + k, 0, 0);
                VECCPSC(rank_k, -1.0, v_Ggt_stripe_p, 0, v_Ggt_tilde_p + k, 0);
                TRSV_LTN(rank_k, Ggt_stripe_p, 0, 0, v_Ggt_tilde_p + k, 0, v_Ggt_tilde_p + k, 0);
                //                 // the following command copies the top block matrix (LU) to the bottom because it it needed later
                //                 GECP(rank_k, gamma_k, Ggt_stripe_p, 0, 0, Ggt_tilde_p + k, nu - rank_k + nx + 1, 0);
                //                 // permutations
                //                 TRTR_L(nu + nx, RSQrqt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0); // copy lower part of RSQ to upper part
                //                 (Pr_p + k)->PM(rank_k, RSQrqt_tilde_p + k);                          // TODO make use of symmetry
                (Pr_p + k)->PV(rank_k, v_RSQrqt_tilde_p + k, 0);
                //                 (Pr_p + k)->MPt(rank_k, RSQrqt_tilde_p + k);
                //                 // GL <- Ggt_tilde_k @ RSQ[:rho,:nu+nx] + RSQrqt[rho:nu+nx+1, rho:] (with RSQ[:rho,:nu+nx] = RSQrqt[:nu+nx,:rho]^T)
                //                 // GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                //                 // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
                //                 // GgLt[0].print();
                //                 GECP(nu - rank_k + nx + 1, nu + nx, RSQrqt_tilde_p + k, rank_k, 0, GgLt_p, 0, 0);
                VECCP(nu + nx, v_RSQrqt_tilde_p + k, 0, v_GgLt_p, 0);
                //                 GEMM_NT(nu - rank_k + nx + 1, nu + nx, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, RSQrqt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, 0, GgLt_p, 0, 0);
                GEMV_N(nu + nx, rank_k, 1.0, RSQrqt_tilde_p + k, 0, 0, v_Ggt_tilde_p + k, 0, 1.0, v_GgLt_p, 0, v_GgLt_p, 0);
                //                 // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
                //                 SYRK_LN_MN(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                GEMV_N(nu + nx - rank_k, rank_k, 1.0, Ggt_tilde_p + k, 0, 0, v_GgLt_p, 0, 1.0, v_GgLt_p, rank_k, v_RSQrqt_hat_p, 0);
                //                 // GEMM_NT(nu - rank_k + nx + 1, nu + nx - rank_k, rank_k, 1.0, GgLt_p, 0, 0, Ggt_tilde_p + k, 0, 0, 1.0, GgLt_p, 0, rank_k, RSQrqt_hat_p, 0, 0);
                v_RSQrq_hat_curr_p = v_RSQrqt_hat_p;
                //    RSQrq_hat_curr_p = RSQrqt_hat_p;
            }
            else
            {
                v_RSQrq_hat_curr_p = v_RSQrqt_tilde_p + k;
                //    RSQrq_hat_curr_p = RSQrqt_tilde_p + k;
            }
            // if (k==K-2) blasfeo_print_dvec(nu+nx-rank_k, v_RSQrq_hat_curr_p, 0);
        }
        //         //////// SCHUR
        {
            if (nu - rank_k > 0)
            {
                //                 // DLlt_k = [chol(R_hatk); Llk@chol(R_hatk)^-T]
                //                 POTRF_L_MN(nu - rank_k + nx + 1, nu - rank_k, RSQrq_hat_curr_p, 0, 0, Llt_p + k, 0, 0);
                TRSV_LNN(nu - rank_k, Llt_p + k, 0, 0, v_RSQrq_hat_curr_p, 0, v_Llt_p + k, 0);
                //                 if (!check_reg(nu - rank_k, Llt_p + k, 0, 0))
                //                     return 1;
                //                 // Pp_k = Qq_hatk - L_k^T @ Ll_k
                //                 // SYRK_LN_MN(nx+1, nx, nu-rank_k, -1.0,Llt_p+k, nu-rank_k,0, Llt_p+k, nu-rank_k,0, 1.0, RSQrq_hat_curr_p, nu-rank_k, nu-rank_k,Pp+k,0,0); // feature not implmented yet
                ///// TODO get rid ot this operation
                GECP(nx + 1, nu - rank_k, Llt_p + k, nu - rank_k, 0, Llt_shift_p, 0, 0); // needless operation because feature not implemented yet
                VECCP(nu - rank_k, v_Llt_p + k, 0, v_Llt_shift_p, 0);
                //                 // SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                //                 GECP(nx + 1, nx, RSQrq_hat_curr_p, nu - rank_k, nu - rank_k, Ppt_p + k, 0, 0);
                VECCP(nx, v_RSQrq_hat_curr_p, nu - rank_k, v_Ppt_p + k, 0);
                //                 SYRK_LN_MN(nx + 1, nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, Llt_shift_p, 0, 0, 1.0, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
                GEMV_N(nx, nu - rank_k, -1.0, Llt_shift_p, 0, 0, v_Llt_shift_p, 0, 1.0, v_Ppt_p + k, 0, v_Ppt_p + k, 0);
                //                 // next steps are for better accuracy
                if (increased_accuracy)
                {
                    //                     // copy eta
                    //                     GETR(nu - rank_k, gamma_k - rank_k, Ggt_stripe_p, rank_k, rank_k, Ggt_stripe_p, 0, 0);
                    //                     // blasfeo_print_dmat(gamma_k-rank_k, nu-rank_k, Ggt_stripe_p, 0,0);
                    //                     // eta L^-T
                    //                     TRSM_RLTN(gamma_k - rank_k, nu - rank_k, 1.0, Llt_p + k, 0, 0, Ggt_stripe_p, 0, 0, Ggt_stripe_p, 0, 0);
                    //                     // ([S^T \\ r^T] L^-T) @ (L^-1 eta^T)
                    //                     // (eta L^-T) @ ([S^T \\ r^T] L^-T)^T
                    //                     GEMM_NT(gamma_k - rank_k, nx + 1, nu - rank_k, -1.0, Ggt_stripe_p, 0, 0, Llt_p + k, nu - rank_k, 0, 1.0, Hh_p + k, 0, 0, Hh_p + k, 0, 0);
                    //                     // keep (L^-1 eta^T) for forward recursion
                    //                     GETR(gamma_k - rank_k, nu - rank_k, Ggt_stripe_p, 0, 0, Ggt_tilde_p + k, 0, rank_k);
                    GEMV_T(nu - rank_k, gamma_k - rank_k, -1.0, Ggt_tilde_p + k, 0, rank_k, v_Llt_p + k, 0, 1.0, v_Hh_p + k, 0, v_Hh_p + k, 0);
                }
            }
            else
            {
                //                 GECP(nx + 1, nx, RSQrq_hat_curr_p, 0, 0, Ppt_p + k, 0, 0);
                VECCP(nx, v_RSQrq_hat_curr_p, 0, v_Ppt_p + k, 0);
            }
            //             TRTR_L(nx, Ppt_p + k, 0, 0, Ppt_p + k, 0, 0);
        }
    }
    //     rankI = 0;
    //     //////// FIRST_STAGE
    {
        const int nx = nx_p[0];
        int gamma_I = gamma_p[0] - rho_p[0];
        //         if (gamma_I > nx)
        //         {
        //             return -3;
        //         }
        if (gamma_I > 0)
        {
            //             GETR(gamma_I, nx + 1, Hh_p + 0, 0, 0, HhIt_p, 0, 0); // transposition may be avoided
            VECCP(gamma_I, v_Hh_p + 0, 0, v_HhIt_p, 0);

            //             // HhIt[0].print();
            //             LU_FACT_transposed(gamma_I, nx + 1, nx, rankI, HhIt_p, PlI_p, PrI_p);
            PlI_p->PV(rankI, v_HhIt_p, 0);
            // L1^-1 g_stipe[:rho]
            TRSV_UTU(rankI, HhIt_p, 0, 0, v_HhIt_p, 0, v_HhIt_p, 0);
            // -L2 L1^-1 g_stripe[:rho] + g_stripe[rho:]
            GEMV_T(rankI, gamma_I - rankI, -1.0, HhIt_p, 0, rankI, v_HhIt_p, 0, 1.0, v_HhIt_p, rankI, v_HhIt_p, rankI);
            //             if (rankI < gamma_I)
            //                 return -2;
            //             // PpIt_tilde <- Ggt[rankI:nx+1, :rankI] L-T (note that this is slightly different from the implementation)
            //             TRSM_RLNN(nx - rankI + 1, rankI, -1.0, HhIt_p, 0, 0, HhIt_p, rankI, 0, GgIt_tilde_p, 0, 0);
            VECCPSC(rankI, -1.0, v_HhIt_p, 0, v_GgIt_tilde_p, 0);
            TRSV_LTN(rankI, HhIt_p, 0, 0, v_GgIt_tilde_p, 0, v_GgIt_tilde_p, 0);
            //             // permutations
            //             (PrI_p)->PM(rankI, Ppt_p); // TODO make use of symmetry
            (PrI_p)->PV(rankI, v_Ppt_p, 0);
            //             (PrI_p)->MPt(rankI, Ppt_p);
            //             // // GL <- GgIt_tilde @ Pp[:rankI,:nx] + Ppt[rankI:nx+1, rankI:] (with Pp[:rankI,:nx] = Ppt[:nx,:rankI]^T)
            //             // GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
            //             // split up because valgrind was giving invalid read errors when C matrix has nonzero row offset
            //             GECP(nx - rankI + 1, nx, Ppt_p, rankI, 0, GgLIt_p, 0, 0);
            VECCP(nx, v_Ppt_p, 0, v_GgLIt_p, 0);
            //             GEMM_NT(nx - rankI + 1, nx, rankI, 1.0, GgIt_tilde_p, 0, 0, Ppt_p, 0, 0, 1.0, GgLIt_p, 0, 0, GgLIt_p, 0, 0);
            GEMV_N(nx, rankI, 1.0, Ppt_p, 0, 0, v_GgIt_tilde_p, 0, 1.0, v_GgLIt_p, 0, v_GgLIt_p, 0);
            //             // // RSQrqt_hat = GgLt[nu-rank_k + nx +1, :rank_k] * G[:rank_k, :nu+nx] + GgLt[rank_k:, :]  (with G[:rank_k,:nu+nx] = Gt[:nu+nx,:rank_k]^T)
            //             SYRK_LN_MN(nx - rankI + 1, nx - rankI, rankI, 1.0, GgLIt_p, 0, 0, GgIt_tilde_p, 0, 0, 1.0, GgLIt_p, 0, rankI, PpIt_hat_p, 0, 0);
            GEMV_N(nx - rankI, rankI, 1.0, GgIt_tilde_p, 0, 0, v_GgLIt_p, 0, 1.0, v_GgLIt_p, rankI, v_PpIt_hat_p, 0);
            //             // TODO skipped if nx-rankI = 0
            //             POTRF_L_MN(nx - rankI + 1, nx - rankI, PpIt_hat_p, 0, 0, LlIt_p, 0, 0);
            TRSV_LNN(nx - rankI, LlIt_p, 0, 0, v_PpIt_hat_p, 0, v_LlIt_p, 0);
            //             if (!check_reg(nx - rankI, LlIt_p, 0, 0))
            //                 return 2;
        }
        else
        {
            //             rankI = 0;
            //             POTRF_L_MN(nx + 1, nx, Ppt_p, 0, 0, LlIt_p, 0, 0);
            TRSV_LNN(nx, LlIt_p, 0, 0, v_LlIt_p, 0, v_LlIt_p, 0);
            //             if (!check_reg(nx, LlIt_p, 0, 0))
            //                 return 2;
        }
    }
    //     ////// FORWARD_SUBSTITUTION:
    //     // first stage
    {
        const int nx = nx_p[0];
        const int nu = nu_p[0];
        //         // calculate xIb
        //         ROWEX(nx - rankI, -1.0, LlIt_p, nx - rankI, 0, ux_p, nu + rankI);
        VECCPSC(nx - rankI, -1.0, v_LlIt_p, 0, ux_p, nu + rankI);
        //         // assume TRSV_LTN allows aliasing, this is the case in normal BLAS
        TRSV_LTN(nx - rankI, LlIt_p, 0, 0, ux_p, nu + rankI, ux_p, nu + rankI);
        //         // calculate xIa
        //         ROWEX(rankI, 1.0, GgIt_tilde_p, nx - rankI, 0, ux_p, nu);
        VECCP(rankI, v_GgIt_tilde_p, 0, ux_p, nu);
        //         // assume aliasing is possible for last two elements
        GEMV_T(nx - rankI, rankI, 1.0, GgIt_tilde_p, 0, 0, ux_p, nu + rankI, 1.0, ux_p, nu, ux_p, nu);
        //         //// lag
        // ROWEX(rankI, -1.0, Ppt_p, nx, 0, lam_p, 0);
        VECCPSC(rankI, -1.0, v_Ppt_p, 0, lam_p, 0);
        //         // assume aliasing is possible for last two elements
        GEMV_T(nx, rankI, -1.0, Ppt_p, 0, 0, ux_p, nu, 1.0, lam_p, 0, lam_p, 0);
        //         // U^-T
        TRSV_LNN(rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
        //         // L^-T
        TRSV_UNU(rankI, rankI, HhIt_p, 0, 0, lam_p, 0, lam_p, 0);
        (PlI_p)->PtV(rankI, lam_p, 0);
        (PrI_p)->PtV(rankI, ux_p, nu);
        // blasfeo_print_dvec(rankI, lam_p, 0);
    }
    //     int *offs_ux = (int *)OCP->aux.ux_offs.data();
    //     int *offs_g = (int *)OCP->aux.g_offs.data();
    //     int *offs_dyn_eq = (int *)OCP->aux.dyn_eq_offs.data();
    //     // other stages
    //     // for (int k = 0; k < K - 1; k++)
    //     // int dyn_eqs_ofs = offs_g[K - 1] + ng_p[K - 1]; // this value is incremented at end of recursion
    for (int k = 0; k < K - 1; k++)
    {
        const int nx = nx_p[k];
        const int nu = nu_p[k];
        const int nxp1 = nx_p[k + 1];
        const int nup1 = nu_p[k + 1];
        const int offsp1 = offs_ux[k + 1];
        const int offs = offs_ux[k];
        const int rho_k = rho_p[k];
        const int numrho_k = nu - rho_k;
        const int offs_g_k = offs_g[k];
        const int offs_dyn_eq_k = offs_dyn_eq[k];
        const int offs_dyn_k = offs_dyn[k];
        const int offs_g_kp1 = offs_g[k + 1];
        const int gammamrho_k = gamma_p[k] - rho_p[k];
        const int gamma_k = gamma_p[k];
        const int gammamrho_kp1 = gamma_p[k + 1] - rho_p[k + 1];
        if (numrho_k > 0)
        {
            //             /// calculate ukb_tilde
            //             // -Lkxk - lk
            //             ROWEX(numrho_k, -1.0, Llt_p + k, numrho_k + nx, 0, ux_p, offs + rho_k);
            VECCPSC(numrho_k, -1.0, v_Llt_p + k, 0, ux_p, offs + rho_k);
            if (increased_accuracy)
            {
                GEMV_N(nu - rho_k, gamma_k - rho_k, -1.0, Ggt_tilde_p + k, 0, rho_k, lam_p, offs_g_k, 1.0, ux_p, offs + rho_k, ux_p, offs + rho_k);
            }
            //             // assume aliasing of last two eliments is allowed
            GEMV_T(nx, numrho_k, -1.0, Llt_p + k, numrho_k, 0, ux_p, offs + nu, 1.0, ux_p, offs + rho_k, ux_p, offs + rho_k);
            TRSV_LTN(numrho_k, Llt_p + k, 0, 0, ux_p, offs + rho_k, ux_p, offs + rho_k);
        }
        //         /// calcualate uka_tilde
        if (rho_k > 0)
        {
            // ROWEX(rho_k, 1.0, Ggt_tilde_p + k, numrho_k + nx, 0, ux_p, offs);
            VECCP(rho_k, v_Ggt_tilde_p + k, 0, ux_p, offs);
            GEMV_T(nx + numrho_k, rho_k, 1.0, Ggt_tilde_p + k, 0, 0, ux_p, offs + rho_k, 1.0, ux_p, offs, ux_p, offs);
            //             // calculate lamda_tilde_k
            //             // copy vk to right location
            //             // we implemented a version of vector copy that starts with copy of last element, to avoid aliasing error
            VECCPR(gammamrho_k, lam_p, offs_g_k, lam_p, offs_g_k + rho_k);
            // ROWEX(rho_k, -1.0, RSQrqt_tilde_p + k, nu + nx, 0, lam_p, offs_g_k);
            VECCPSC(rho_k, -1.0, v_RSQrqt_tilde_p + k, 0, lam_p, offs_g_k);
            //             // assume aliasing of last two eliments is allowed
            GEMV_T(nu + nx, rho_k, -1.0, RSQrqt_tilde_p + k, 0, 0, ux_p, offs, 1.0, lam_p, offs_g_k, lam_p, offs_g_k);
            //             // nu-rank_k+nx,0
            //             // needless copy because feature not implemented yet in trsv_lnn
            GECP(rho_k, gamma_k, Ggt_tilde_p + k, nu - rho_k + nx + 1, 0, AL_p, 0, 0);
            //             // U^-T
            TRSV_LNN(rho_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
            //             // L^-T
            TRSV_UNU(rho_k, gamma_k, AL_p, 0, 0, lam_p, offs_g_k, lam_p, offs_g_k);
            (Pl_p + k)->PtV(rho_k, lam_p, offs_g_k);
            (Pr_p + k)->PtV(rho_k, ux_p, offs);
        }
        //         // calculate xkp1
        //         ROWEX(nxp1, 1.0, BAbt_p + k, nu + nx, 0, ux_p, offsp1 + nup1);
        VECCP(nxp1, rhs_b_p, offs_dyn_k, ux_p, offsp1 + nup1);
        GEMV_T(nu + nx, nxp1, 1.0, BAbt_p + k, 0, 0, ux_p, offs, 1.0, ux_p, offsp1 + nup1, ux_p, offsp1 + nup1);
        //         // calculate lam_dyn xp1
        //         ROWEX(nxp1, 1.0, Ppt_p + (k + 1), nxp1, 0, lam_p, offs_dyn_eq_k);
        VECCP(nxp1, v_Ppt_p + (k + 1), 0, lam_p, offs_dyn_eq_k);
        GEMV_T(nxp1, nxp1, 1.0, Ppt_p + (k + 1), 0, 0, ux_p, offsp1 + nup1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
        GEMV_T(gammamrho_kp1, nxp1, 1.0, Hh_p + (k + 1), 0, 0, lam_p, offs_g_kp1, 1.0, lam_p, offs_dyn_eq_k, lam_p, offs_dyn_eq_k);
    }
    return 0;
};
