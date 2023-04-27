#pragma once
#include <ocp/StageOCPApplication.hpp>
#include <cocp.hpp>
#include <memory>
namespace genriccati_benchmark
{
    class FatropProblemWrap
    {
    public:
        FatropProblemWrap(const std::shared_ptr<fatrop::OCPAbstract> &ocp)
        {
            auto adapter = std::make_shared<fatrop::OCPAdapter>(ocp);
            std::shared_ptr<fatrop::FatropOptions> fatropoptions_ = std::make_shared<fatrop::FatropOptions>();
            std::shared_ptr<fatrop::FatropData> fatropdata_;
            std::shared_ptr<fatrop::FatropNLP> nlp_;
            std::shared_ptr<fatrop::FatropPrinter> printer_ = std::make_shared<fatrop::FatropPrinter>();
            std::shared_ptr<fatrop::FatropNLP> nlp(fatrop::FatropOCPBuilder(ocp, fatropoptions_, printer_).build(adapter));
            fatrop::AlgBuilder algbuilder;
            algbuilder.set_printer(printer_);
            std::shared_ptr<fatrop::Journaller> journaller_;
            algbuilder.build_fatrop_algorithm_objects(nlp, fatropoptions_, fatropdata_, journaller_);
            fatropoptions_->register_option(fatrop::IntegerOption::un_bounded("print_level", "print level", &printer_->print_level(), 10));
            fatropalg_ = algbuilder.build_algorithm();
        }
        void populate_cocp(const gen_riccati::COCP &cocp)
        {
        }
        // gen_riccati::COCP create_cocp()
        // {
        // }

    private:
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
            // eval initial constraint jacobian
            fatropalg_->eval_constr_jac();
            // eval initial lagrangian hessian
            fatropalg_->eval_lag_hess();
        }
        std::shared_ptr<fatrop::FatropAlg> fatropalg_;
        // void eval_quantities()
    };
}