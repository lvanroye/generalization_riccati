#include "sparse_system.hpp"
using namespace symspals;
int main()
{
    SparseLinearSystem ls = SparseLinearSystem();
    auto x = ls.variable(1);
    auto y = ls.variable(1);
    ls.add_equation(69*x + 13*y, Matrix<Expression>(420));
    auto sps = ls.get_sparsity();
    auto coeffs = ls.eval_coeffs();
    auto rhs = ls.eval_rhs();
    
    cout << "Sparsity pattern: " << endl;
    for(size_t i = 0; i<sps.size(); i++)
    {
        auto sp = sps[i];
        auto coeff = coeffs[i];
        cout <<  "(" <<  sp.row << ", " << sp.col << ") --> " << coeff << endl;
    }
    cout << "rhs: " << endl;
    for(auto r: rhs)
        cout << r << endl;
    return 0;
}