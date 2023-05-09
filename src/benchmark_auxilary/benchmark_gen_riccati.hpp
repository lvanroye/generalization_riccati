#pragma once
#include <gen_riccati.hpp>
#include <numeric_vector.hpp>
extern "C"
{
#include <timing.h>
}
namespace genriccati_benchmark
{
    class BenchmarkGenRiccati
    {
    public:
        BenchmarkGenRiccati(const gen_riccati::COCP &cocp) : ocp_ls_riccati(cocp), ux(sum(cocp.nu + cocp.nx), 1), lam(sum(cocp.ng + cocp.nx) - cocp.nx[0], 1), rhs_rq(sum(cocp.nu + cocp.nx), 2), rhs_b(sum(cocp.nx) - cocp.nx[0], 2), rhs_g(sum(cocp.ng), 2)
        {
        }
        void solve(gen_riccati::COCP &cocp)
        {
            blasfeo_timer timer;
            blasfeo_tic(&timer);
            ocp_ls_riccati.solve_pd_sys_normal(cocp, 0.0, ux[0], lam[0]);
            last_time = blasfeo_toc(&timer);
            ocp_ls_riccati.get_rhs(cocp, rhs_rq[1], rhs_b[1], rhs_g[1]);
            ocp_ls_riccati.compute_pd_sys_times_vec(cocp, 0.0, 0.0, ux[0], lam[0], rhs_rq[0], rhs_b[0], rhs_g[0]);
            axpby(1.0, rhs_rq[0], 1.0, rhs_rq[1], rhs_rq[0]);
            axpby(1.0, rhs_b[0], 1.0, rhs_b[1], rhs_b[0]);
            axpby(1.0, rhs_g[0], 1.0, rhs_g[1], rhs_g[0]);
            last_res = std::max(Linf(rhs_rq[0]), std::max(Linf(rhs_b[0]), Linf(rhs_g[0])));
        }
        void last_stats(double &time, double &res)
        {
            time = last_time;
            res = last_res;
        }
        void controle()
        {
            // first 10 numbers of ux[0]
            for (int i = 0; i < 10; i++)
                std::cout << VECEL((VEC*) ux[0], i) << " ";
        }
        gen_riccati::OCPLSRiccati ocp_ls_riccati;
        gen_riccati::FatropMemoryVecBF ux;
        gen_riccati::FatropMemoryVecBF lam;
        gen_riccati::FatropMemoryVecBF rhs_rq;
        gen_riccati::FatropMemoryVecBF rhs_b;
        gen_riccati::FatropMemoryVecBF rhs_g;
        double last_time;
        double last_res;
    };
} // namespace genriccati_benchmark