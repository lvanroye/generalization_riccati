#pragma once
#include <blasfeo.h>
#include "numeric_vector.hpp"
#include "blasfeo_wrapper.hpp"
namespace gen_riccati
{
    class COCP
    {
        public:
        COCP(const int K, const NumericVector &nu, const NumericVector &nx, const NumericVector &ng) : K(K), nu(nu), nx(nx), ng(ng), RSQrqt(nu + nx + 1, nu + nx, K), BAbt(nx + nu+ 1, rotate(nx, 1), K-1), Ggt(nx+nu+1, ng, K)
        {
        }
        COCP(COCP&& cpy): K(cpy.K), nu(cpy.nu), nx(cpy.nx), ng(cpy.ng), RSQrqt(std::move(cpy.RSQrqt)), BAbt(std::move(cpy.BAbt)), Ggt(std::move(cpy.Ggt))
        {
        }
        int K;
        /// number of controls
        const NumericVector nu;
        /// number of states
        const NumericVector nx;
        /// number of stagewise equality constraints
        const NumericVector ng;
        /// small-scale Hessian
        FatropMemoryMatBF RSQrqt;
        /// small-scale Jacobian dynamics
        FatropMemoryMatBF BAbt;
        /// small-scale Jacobian stagewise eq constraints
        FatropMemoryMatBF Ggt;
    };
} // namespace gen_riccati