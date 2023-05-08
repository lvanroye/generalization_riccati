#pragma once
#include <cocp.hpp>
#include <sparse_system.hpp>
namespace genriccati_benchmark
{
    class BlasfeoMatrixAdapter : public symspals::Matrix<double>
    {
    public:
        BlasfeoMatrixAdapter(const int m, const int n, MAT *mat, const int ai, const int aj) : symspals::Matrix<double>(m, n)
        {
            for (int i = 0; i < m; i++)
                for (int j = 0; j < n; j++)
                    this->operator()(i, j) = BLASFEO_DMATEL(mat, ai + i, aj + j);
        };
    };
    class BenchmarkSparse
    {
    public:
        BenchmarkSparse(const gen_riccati::COCP &cocp)
        {
            gen_riccati::NumericVector nu = cocp.nu;
            gen_riccati::NumericVector nx = cocp.nx;
            gen_riccati::NumericVector ng = cocp.ng;
            const int K = cocp.K;
            // add all variables
            for (int k = 0; k < K; k++)
            {
                RSQ_blocks.push_back(kkt_system_.parameter(nu[k] + nx[k], nu[k] + nx[k]));
                rq_grads.push_back(kkt_system_.parameter(nu[k] + nx[k], 1));
                ux_syms.push_back(kkt_system_.variable(RSQ_blocks[k], rq_grads[k]));
            }
            // add dynamics constraints
            for (int k = 0; k < K - 1; k++)
            {
                BA_blocks.push_back(kkt_system_.parameter(nx[k + 1], nu[k] + nx[k]));
                b_rhss.push_back(kkt_system_.parameter(nx[k + 1], 1));
                dyn_lags.push_back(kkt_system_.add_constraint(-ux_syms[k + 1].block(nu[k + 1], 0, nx[k + 1], 1) + (*BA_blocks[k]) * ux_syms[k], -b_rhss[k]));
            }
            // add the stagewise constraints
            for (int k = 0; k < K; k++)
            {
                G_blocks.push_back(kkt_system_.parameter(ng[k], nu[k] + nx[k]));
                g_rhss.push_back(kkt_system_.parameter(ng[k], 1));
                g_lags.push_back(kkt_system_.add_constraint((*G_blocks[k]) * ux_syms[k], -g_rhss[k]));
            }
            set_value(cocp);
            to_ux = kkt_system_.evaluator(vec(ux_syms));
            // concatenate the lagrangian multipliers
            vector<symspals::Matrix<symspals::Expression>> lam_syms = dyn_lags;
            lam_syms.insert(lam_syms.end(), g_lags.begin(), g_lags.end());
            to_lam = kkt_system_.evaluator(vec(lam_syms));
        }
        void set_value(const gen_riccati::COCP &cocp)
        {
            gen_riccati::NumericVector nu = cocp.nu;
            gen_riccati::NumericVector nx = cocp.nx;
            gen_riccati::NumericVector ng = cocp.ng;
            const int K = cocp.K;
            for (int k = 0; k < K; k++)
            {
                kkt_system_.set_value(RSQ_blocks[k], BlasfeoMatrixAdapter(nu[k] + nx[k], nu[k] + nx[k], (MAT *)cocp.RSQrqt[k], 0, 0));
                kkt_system_.set_value(rq_grads[k], BlasfeoMatrixAdapter(1, nu[k] + nx[k], (MAT *)cocp.RSQrqt[k], nu[k] + nx[k], 0).transpose());
                kkt_system_.set_value(G_blocks[k], BlasfeoMatrixAdapter(nu[k] + nx[k], ng[k], (MAT *)cocp.Ggt[k], 0, 0).transpose());
                kkt_system_.set_value(g_rhss[k], BlasfeoMatrixAdapter(1, ng[k], (MAT *)cocp.Ggt[k], nu[k] + nx[k], 0).transpose());
            }
            for (int k = 0; k < K - 1; k++)
            {
                kkt_system_.set_value(BA_blocks[k], BlasfeoMatrixAdapter(nu[k] + nx[k], nx[k + 1], (MAT *)cocp.BAbt[k], 0, 0).transpose());
                kkt_system_.set_value(b_rhss[k], BlasfeoMatrixAdapter(1, nx[k + 1], (MAT *)cocp.BAbt[k], nu[k] + nx[k], 0).transpose());
            }
        }
        void to_ux_lam(const vector<double> & sol, vector<double> &ux, vector<double> &lam)
        {
            ux = to_ux.Eval(sol);
            lam = to_lam.Eval(sol);
        };
        symspals::KKTSystem kkt_system_;
        vector<symspals::Matrix<symspals::Expression>> ux_syms;
        vector<symspals::Parameter> RSQ_blocks;
        vector<symspals::Parameter> rq_grads;
        vector<symspals::Parameter> BA_blocks;
        vector<symspals::Parameter> b_rhss;
        vector<symspals::Parameter> G_blocks;
        vector<symspals::Parameter> g_rhss;
        vector<symspals::Matrix<symspals::Expression>> dyn_lags;
        vector<symspals::Matrix<symspals::Expression>> g_lags;
        symspals::Function to_ux;
        symspals::Function to_lam;
    };
}