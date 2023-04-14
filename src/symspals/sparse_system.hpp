#ifndef SPARSEMATRIXINCLUDED
#define SPARSEMATRIXINCLUDED
#include "expressions.hpp"
namespace symspals
{
    struct Parameter : Matrix<Expression>
    {
        Parameter(const string &name, const int m, const int n) : Matrix<Expression>(SymMatrix(name, m, n)){};
        void set_value(const Matrix<double> &value)
        {
            dirty = false;
            value_ = value;
        }
        const Matrix<double> &get_value()
        {
            if (dirty)
                throw runtime_error("Parameter value not set");
            return value_;
        }

    private:
        bool dirty = true;
        Matrix<double> value_;
    };

    class SparseLinearSystem
    {
    public:
        SymMatrix variable(int size)
        {
            SymMatrix sym("x", size, 1);
            variables.push_back(sym);
            dirty = true;
            return sym;
        }
        shared_ptr<Parameter> parameter(int m, int n)
        {
            auto sym = SymMatrix("p", m, n);
            parameters_syms.push_back(sym);
            parameters.push_back(make_shared<Parameter>("p", m, n));
            dirty = true;
            return parameters.back();
        }

        void add_equation(const Matrix<Expression> &lhs, const Matrix<Expression> &rhs)
        {
            equations.push_back(lhs);
            rhss.push_back(rhs);
            dirty = true;
        }
        void add_equation(const Matrix<Expression> &lhs, const Expression &rhs)
        {
            equations.push_back(lhs);
            rhss.push_back(Matrix<Expression>(rhs));
            dirty = true;
        }
        vector<Expression> get_coeffs()
        {
            if (dirty)
                make_clean();
            return coeffs;
        }
        vector<double> eval_coeffs()
        {
            if (dirty)
                make_clean();
            return coeffs_f.Eval(parameters_vals);
        }
        vector<double> eval_rhs()
        {
            if (dirty)
                make_clean();
            return rhs_f.Eval(parameters_vals);
        }
        vector<Index> get_sparsity()
        {
            if (dirty)
                make_clean();
            return sparsity;
        }
        // void solver(const )
        void solver(const string &solver_name)
        {
            dirty = true;
            if (solver_name == "ma57")
            {
                // ma57
            }
            else if (solver_name == "mumps")
            {
                // mumps
            }
            else if (solver_name == "pardiso")
            {
                // pardiso
            }
            else
            {
                throw runtime_error("Unknown solver");
            }
        }
        void solve(vector<double> &solution)
        {
        }
        vector<Matrix<Expression>> variables;
        vector<Matrix<Expression>> parameters_syms;
        vector<shared_ptr<Parameter>> parameters;
        vector<double> parameters_vals;
        vector<Matrix<Expression>> equations;
        vector<Matrix<Expression>> rhss;
        vector<Index> sparsity;
        vector<Expression> coeffs;
        Function coeffs_f; // computes coefficients from parameters
        Function rhs_f;
        string triplet_order = "ma57";

    private:
        void make_clean()
        {
            sparsity.resize(0);
            coeffs.resize(0);
            // coeffs.resize(0);
            auto eq_vec = vec(equations);
            auto var_vec = vec(variables);
            auto triplets = GetCoefficients(eq_vec, var_vec);
            const int n_eq = equations.size();
            const int n_var = variables.size();
            // go trough all triplets
            if (triplet_order == "ma57")
            {
                for (auto triplet : triplets)
                {
                    sparsity.push_back(triplet.index);
                    coeffs.push_back(triplet.value);
                }
            }

            // initialize the function that computes the coefficients
            coeffs_f = Function(vec(parameters_syms), coeffs);
            rhs_f = Function(vec(parameters_syms), vec(rhss));

            // evaluate the parameters
            for (auto p : parameters)
            {
                auto p_vals = vec(p->get_value());
                parameters_vals.insert(parameters_vals.end(), p_vals.begin(), p_vals.end());
            }
            dirty = false;
        }
        bool dirty = true;
    };
}

#endif // SPARSEMATRIXINCLUDED