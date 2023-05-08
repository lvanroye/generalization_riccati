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
        BenchmarkGenRiccati(const gen_riccati::COCP &cocp) : ocp_ls_riccati(cocp), ux(sum(cocp.nu + cocp.nx), 1), lam(sum(cocp.ng + cocp.nx) - cocp.nx[0], 1)
        {
        }
        void solve(gen_riccati::COCP &cocp)
        {
            blasfeo_timer timer;
            blasfeo_tic(&timer);
            ocp_ls_riccati.solve_pd_sys_normal(cocp, 0.0, ux[0], lam[0]);
            blasfeo_toc(&timer);
            std::cout << "solve_pd_sys_normal: " << blasfeo_toc(&timer) << std::endl;
        }
        gen_riccati::OCPLSRiccati ocp_ls_riccati;
        gen_riccati::FatropMemoryVecBF ux;
        gen_riccati::FatropMemoryVecBF lam;
    };
} // namespace genriccati_benchmark