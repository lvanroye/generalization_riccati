#pragma once
#include <ocp/StageOCPApplication.hpp>
#include <memory>
#include <numeric_vector.hpp>
#include <cocp.hpp>
#include <blasfeo.h>
namespace genriccati_benchmark
{
    /***
     * this class hacks into the fatrop internals to get access to the quantities needed to create a COCP
     * it is not meant to be used in production code, and solely serves the purpose of benchmarking for the paper
     */
    class FatropProblemWrap
    {
    public:
        FatropProblemWrap(const std::shared_ptr<fatrop::OCPAbstract> &ocp)
        {
            ocp_ = ocp;
            auto adapter = std::make_shared<fatrop::OCPAdapter>(ocp);
            std::shared_ptr<fatrop::FatropOptions> fatropoptions_ = std::make_shared<fatrop::FatropOptions>();
            std::shared_ptr<fatrop::FatropData> fatropdata_;
            std::shared_ptr<fatrop::FatropNLP> nlp_;
            std::shared_ptr<fatrop::FatropPrinter> printer_ = std::make_shared<fatrop::FatropPrinter>();
            fatropocp = fatrop::FatropOCPBuilder(ocp, fatropoptions_, printer_).build(adapter);
            std::shared_ptr<fatrop::FatropNLP> nlp = fatropocp;
            fatrop::AlgBuilder algbuilder;
            algbuilder.set_printer(printer_);
            std::shared_ptr<fatrop::Journaller> journaller_;
            algbuilder.build_fatrop_algorithm_objects(nlp, fatropoptions_, fatropdata_, journaller_);
            fatropoptions_->register_option(fatrop::IntegerOption::un_bounded("print_level", "print level", &printer_->print_level(), 10));
            fatropalg_ = algbuilder.build_algorithm();
        }
        void populate_cocp(gen_riccati::COCP &cocp)
        {
            eval_quantities();
            const int K = cocp.K;
            const gen_riccati::NumericVector nu = cocp.nu;
            const gen_riccati::NumericVector nx = cocp.nx;
            const gen_riccati::NumericVector ng = cocp.ng;
            fatrop::OCPKKTMemory &kkt_memory = fatropocp->ocpkktmemory_;
            for (auto k = 0; k < K; k++)
            {
                // populate the small-scale Lagrangian Hessians
                blasfeo_dgecp(nu[k] + nx[k] + 1, nu[k] + nx[k], (MAT *)kkt_memory.RSQrqt[k], 0, 0, (MAT *)cocp.RSQrqt[k], 0, 0);
                // populate the small-scale constraint Jacobians
                blasfeo_dgecp(nu[k] + nx[k] + 1, ng[k], (MAT *)kkt_memory.Ggt[k], 0, 0, (MAT *)cocp.Ggt[k], 0, 0);
            }
            for (auto k = 0; k < K - 1; k++)
            {
                // populate the discretized dynamics
                blasfeo_dgecp(nu[k] + nx[k] + 1, nx[k + 1], (MAT *)kkt_memory.BAbt[k], 0, 0, (MAT *)cocp.BAbt[k], 0, 0);
            }
        }
        void take_step(const vector<double> &delta_x, const vector<double> &delta_lam)
        {
            fatropalg_->fatropdata_->delta_x = delta_x;
            fatropalg_->fatropdata_->lam_calc = delta_lam;
            fatropalg_->fatropdata_->update_trial_step(1.0, 1.0);
            fatropalg_->fatropdata_->accept_trial_step();
        }
        gen_riccati::COCP create_cocp()
        {
            gen_riccati::NumericVector nu;
            gen_riccati::NumericVector nx;
            gen_riccati::NumericVector ng;
            const int K = ocp_->get_horizon_length();
            for (auto k = 0; k < K; k++)
            {
                nu.push_back(ocp_->get_nuk(k));
                nx.push_back(ocp_->get_nxk(k));
                ng.push_back(ocp_->get_ngk(k));
            }
            return gen_riccati::COCP(K, nu, nx, ng);
        }

        void controle()
        {
            fatropalg_->solve_pd_sys(0.0, 0.0, 0.0);
            // ux0
            fatropalg_->fatropdata_->delta_x.block(0, 10).print();
        }

    public:
        void least_squares_dual()
        {

            // initialize the lagrange multipliers
            fatropalg_->eval_obj_grad_curr();
            fatropalg_->eval_constr_jac();
            fatropalg_->perform_initializiation();
            fatropalg_->fatropdata_->accept_dual_initializiaton();
        }
        void eval_quantities()
        {
            fatropalg_->fatropdata_->reset();
            // least_squares_dual();
            // eval initial constraint jacobian
            fatropalg_->eval_constr_jac();
            // eval initial lagrangian hessian
            fatropalg_->eval_lag_hess();
        }
        std::shared_ptr<fatrop::FatropAlg> fatropalg_;
        std::shared_ptr<fatrop::OCPAbstract> ocp_;
        std::shared_ptr<fatrop::FatropOCP> fatropocp;
        // void eval_quantities()
    };
}