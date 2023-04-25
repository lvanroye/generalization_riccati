#ifndef SPARSEMATRIXINCLUDED
#define SPARSEMATRIXINCLUDED
#include "expressions.hpp"
#include "interfaces/sparse_solver_interface.hpp"
#include "interfaces/mumps.hpp"
#include "interfaces/ma57.hpp"
#include <memory>
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
        // implementation assumes LOWER triangular matrix of symmetric system: TODO make this more general
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
        void solver(const string &solver_name)
        {
            dirty = true;
            if (solver_name == "ma57")
            {
                linear_solver = "ma57";
            }
            else if (solver_name == "mumps")
            {
                linear_solver = "mumps";
            }
            else if (solver_name == "pardiso")
            {
                linear_solver = "pardiso";
            }
            else
            {
                throw runtime_error("Unknown solver");
            }
        }
        Function evaluator(const Matrix<Expression> &expr)
        {
            return Function(vec(variables), vec(expr));
        }
        void solve(vector<double> &solution)
        {
            if (dirty)
                make_clean();
            eval_params();
            solver_ptr->set_coefficient_matrix(coeffs_f.Eval(parameters_vals));
            solution = rhs_f.Eval(parameters_vals);
            solver_ptr->solve(solution);
        }
        bool is_lower_triangular()
        {
            // iterate through triplets and check wether lower triangular
            return true;
        }
        void residu(const vector<double> &sol, vector<double> &res)
        {
            int n_vars = variables.size();
            res.resize(n_vars);
            vector<double> coeffs_vals = eval_coeffs();
            vector<double> rhs_vals = eval_rhs();
            // add rhs to res
            for (int i = 0; i < n_vars; i++)
            {
                res[i] = -rhs_vals[i];
            }
            // iterate over all non-zero coefficients
            for (size_t i = 0; i < coeffs_vals.size(); i++)
            {
                int row = sparsity[i].row;
                int col = sparsity[i].col;
                res[row] += coeffs_vals[i] * sol[col];
                // assume only lower triangular part of SYMMETRIC matrix is stored
                if (row != col)
                    res[col] += coeffs_vals[i] * sol[row];
            }
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
        string linear_solver = "ma57";
        unique_ptr<SparseSolverInterface> solver_ptr;

    private:
        void make_clean()
        {
            sparsity.resize(0);
            coeffs.resize(0);
            // coeffs.resize(0);
            auto eq_vec = vec(equations);
            auto var_vec = vec(variables);
            auto triplets = GetCoefficients(eq_vec, var_vec);
            // const int n_eq = equations.size();
            // const int n_var = variables.size();
            // go trough all triplets
            for (auto triplet : triplets)
            {
                sparsity.push_back(triplet.index);
                coeffs.push_back(triplet.value);
            }
            if (linear_solver == "ma57")
            {
                solver_ptr = make_unique<InterfaceMA57>(var_vec.size(), sparsity);
                solver_ptr->preprocess();
            }
            else if (linear_solver == "mumps")
            {
                solver_ptr = make_unique<InterfaceMUMPS>(var_vec.size(), sparsity);
                solver_ptr->preprocess();
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

    class KKTSystem : public SparseLinearSystem
    {
    public:
        // lower triangular
        SymMatrix variable(const Matrix<Expression> &hess_block, const Matrix<Expression> &grad)
        {
            int size = hess_block.n_cols();
            auto var = SparseLinearSystem::variable(size);
            add_equation(tril(hess_block) * var, -grad);
            return var;
        }
        SymMatrix add_constraint(const Matrix<Expression> &constraint, const Matrix<Expression> &rhs)
        {
            auto lag = SparseLinearSystem::variable(constraint.n_rows());
            lags.push_back(lag);
            add_equation(constraint, rhs);
            return lag;
        }
        vector<Matrix<Expression>> lags;
    };
}

#endif // SPARSEMATRIXINCLUDED