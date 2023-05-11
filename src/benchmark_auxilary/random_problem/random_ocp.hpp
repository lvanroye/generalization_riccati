#pragma once
#include <random>
#include <ocp/OCPAbstract.hpp>
namespace genriccati_benchmark
{
    /***
     * this class implements all necesssary functions of the OCPAbstract class to create a random OCP to be tested by the recursion
     */
    class RandomOCP : public fatrop::OCPAbstract
    {
    public:
        RandomOCP(const int K, const int nx, const int nu, const int ne_init, const int ne_middle, const int ne_final)
            : K(K), nx(nx), nu(nu), ne_init(ne_init), ne_middle(ne_middle), ne_final(ne_final), gen(seed)
        {
        }
        void set_seed(const int seed)
        {
            this->seed = seed;
        }
        const int K; // horizon length
        const int nx;
        const int nu;
        const int ne_init;
        const int ne_middle;
        const int ne_final;
        int seed = 0;
        std::mt19937 gen;

        // OCPAbstract interface
        int get_nxk(const int k) const override
        {
            return nx;
        }
        int get_nuk(const int k) const override
        {
            if (k == K - 1)
                return 0;
            return nu;
        };
        int get_ngk(const int k) const override
        {
            if (k == 0)
            {
                return ne_init;
            }
            else if (k == K - 1)
            {
                return ne_final;
            }
            else
            {
                return ne_middle;
            }
        };
        int get_n_stage_params_k(const int k) const override
        {
            return 0;
        };
        int get_n_global_params() const override
        {
            return 0;
        };
        int get_default_stage_paramsk(double *stage_params, const int k) const override
        {
            return 0;
        };
        int get_default_global_params(double *global_params) const override
        {
            return 0;
        };
        int get_ng_ineq_k(const int k) const override
        {
            return 0;
        };
        int get_horizon_length() const override
        {
            return K;
        };
        int eval_BAbtk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res,
            const int k) override
        {
            random_matrix(nu + nx + 1, nx, res, 0, 0);
            return 0;
        };
        int eval_RSQrqtk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *lam_dyn_k,
            const double *lam_eq_k,
            const double *lam_eq_ineq_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res,
            const int k) override
        {
            int nu_k = get_nuk(k);
            random_posdef_matrix(nu_k + nx, res, 0, 0);
            random_matrix(1, nu_k + nx, res, nu_k + nx, 0);
            return 0;
        };
        int eval_Ggtk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res,
            const int k) override
        {
            int nu_k = get_nuk(k);
            if (k == 0)
            {
                random_matrix(nu_k + nx + 1, ne_init, res, 0, 0);
            }
            else if (k == K - 1)
            {
                random_matrix(nu_k + nx + 1, ne_final, res, 0, 0);
            }
            else
            {
                random_matrix(nu_k + nx + 1, ne_middle, res, 0, 0);
            }
            return 0;
        };
        int eval_Ggt_ineqk(
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            MAT *res,
            const int k) override
        {
            return 0;
        };
        int eval_bk(
            const double *states_kp1,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k)
        {
            return 0;
        };
        int eval_gk(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override
        {
            return 0;
        };
        int eval_gineqk(
            const double *states_k,
            const double *inputs_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override
        {
            return 0;
        };
        int eval_rqk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override
        {
            return 0;
        };
        int eval_Lk(
            const double *objective_scale,
            const double *inputs_k,
            const double *states_k,
            const double *stage_params_k,
            const double *global_params,
            double *res,
            const int k) override
        {
            return 0;
        };
        int get_boundsk(double *lower, double *upper, const int k) const override
        {
            return 0;
        };
        int get_initial_xk(double *xk, const int k) const override
        {
            return 0;
        };
        int get_initial_uk(double *uk, const int k) const override
        {
            return 0;
        };
        void random_matrix(int n, int m, MAT *A, int ai, int aj)
        {
            // fill A with zeros
            blasfeo_dgese(n, m, 0.0, A, ai, aj);
            std::uniform_real_distribution<> dis(1e-1, 1e0);
            // fill the matrix with random values
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    MATEL(A, ai + i, aj + j) = dis(gen);
                }
            }
        }
        void random_lower_matrix(int n, MAT *A, int ai, int aj)
        {
            // fill A with zeros
            blasfeo_dgese(n, n, 0.0, A, ai, aj);
            std::uniform_real_distribution<> dis(1e-1, 1e1);
            // fill the matrix with random values
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    MATEL(A, ai + i, aj + j) = dis(gen);
                }
                // MATEL(A, ai + i, aj + i) = 1.0;
            }
        }

        void random_posdef_matrix(int n, MAT *A, int ai, int aj)
        {
            // fill A with zeros
            blasfeo_dgese(n, n, 0.0, A, ai, aj);
            // create a temporary matrix B
            MAT B;
            blasfeo_allocate_dmat(n, n, &B);
            // make a lower triangular matrix
            random_lower_matrix(n, &B, 0, 0);
            // make it positive definite by multiplying with its transpose
            // multiply B with its transpose and save to A
            GEMM_NT(n, n, n, 1.0, &B, 0, 0, &B, 0, 0, 0.0, A, ai, aj, A, ai, aj);
            GETR(n, n, A, 0, 0, A, 0, 0);
            // free the temporary matrix
            blasfeo_free_dmat(&B);
        }
    };
}