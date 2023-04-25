#pragma once
#include <vector>
#include "../common.hpp"
using namespace std;
namespace symspals
{
    // vector<double> SparseMV(const vector<Index> &ind, const vector<double> &x, const vector<double> &rhs, const int dim_);
    class SparseSolverInterface
    {
    public:
        virtual void preprocess() = 0;
        virtual void set_coefficient_matrix(const vector<double> &A) = 0;
        virtual void solve(vector<double> &rhs) = 0;
    };
} // namespace symspals