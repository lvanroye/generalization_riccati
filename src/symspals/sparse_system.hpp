#ifndef SPARSEMATRIXINCLUDED
#define SPARSEMATRIXINCLUDED
#include "expressions.hpp"
#include <unordered_map>
namespace symspals
{
    class Parameter : public shared_ptr<Matrix<Expression>>
    {
    public:
        Parameter(const string &name, const int m, const int n) : shared_ptr<Matrix<Expression>>(make_shared<SymMatrix>(name, m, n)){};
        operator Matrix<Expression> &()
        {
            return *this->get();
        };
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
        Parameter parameter(int m, int n)
        {
            auto res = Parameter("p", m, n);
            parameters[res] = Matrix<double>(m, n);
            parameters_syms.push_back(res);
            dirty = true;
            return res;
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
        void set_value(const Parameter &p, const Matrix<double> &value)
        {
            parameters[p] = value;
        };
        void eval_params()
        {
            parameters_vals.resize(0);
            for (auto p : parameters_syms)
            {
                const Matrix<double> &p_vals = parameters[p];
                const vector<double> &p_vals_vec = vec(p_vals);
                parameters_vals.insert(parameters_vals.end(), p_vals_vec.begin(), p_vals_vec.end());
            }
        }
        vector<double> eval_coeffs()
        {
            if (dirty)
                make_clean();
            eval_params();
            return coeffs_f.Eval(parameters_vals);
        }
        vector<double> eval_rhs()
        {
            if (dirty)
                make_clean();
            eval_params();
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
        vector<Parameter> parameters_syms;
        unordered_map<shared_ptr<Matrix<Expression>>, Matrix<double>> parameters;
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
            vector<Expression> parameter_sym_vec;
            for (auto p : parameters_syms)
            {
                auto p_vec = vec(*p);
                parameter_sym_vec.insert(parameter_sym_vec.end(), p_vec.begin(), p_vec.end());
            }

            // initialize the function that computes the coefficients
            coeffs_f = Function(parameter_sym_vec, coeffs);
            rhs_f = Function(parameter_sym_vec, vec(rhss));
            dirty = false;
        }
        bool dirty = true;
    };
}

#endif // SPARSEMATRIXINCLUDED