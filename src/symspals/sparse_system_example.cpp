#include "sparse_system.hpp"
using namespace symspals;
int main()
{
    SparseLinearSystem ls = SparseLinearSystem();
    auto x = ls.variable(1);
    auto y = ls.variable(1);
    auto z = ls.variable(2);
    auto a = ls.parameter(1,1);
    ls.add_equation(69*x + (13. + a)*y, 420);
    ls.add_equation(y, 2 + a);
    ls.add_equation(123.*z, {101., 505.});
    ls.set_value(a, 10000.);
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