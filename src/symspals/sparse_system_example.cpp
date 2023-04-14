#include "sparse_system.hpp"
using namespace symspals;
int main()
{
    SparseLinearSystem ls = SparseLinearSystem();
    auto x = ls.variable(1);
    auto y = ls.variable(1);
    ls.add_equation(x + y, Const(1));
    auto sps = ls.get_sparsity();
    auto coeffs = ls.eval_coeffs();
    
    cout << "Sparsity pattern: " << endl;
    for(size_t i = 0; i<sps.size(); i++)
    {
        auto sp = sps[i];
        auto coeff = coeffs[i];
        cout <<  "(" <<  sp.row << ", " << sp.col << ") --> " << coeff << endl;
    }
    cout << "rhs: " << endl;
    return 0;
}