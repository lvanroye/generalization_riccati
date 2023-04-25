#include "sparse_system.hpp"
using namespace symspals;
int main()
{
    SparseLinearSystem ls = SparseLinearSystem();
    auto x = ls.variable(1);
    auto y = ls.variable(1);
    auto z = ls.variable(2);
    auto a = ls.parameter(1, 1);
    ls.add_equation(69 * x + (13. + a) * y, 420);
    ls.add_equation(y, 2 + a);
    ls.add_equation(123. * z, {101., 505.});
    ls.set_value(a, {10000.});
    auto sps = ls.get_sparsity();
    auto coeffs = ls.eval_coeffs();
    auto rhs = ls.eval_rhs();

    cout << "Sparsity pattern: " << endl;
    for (size_t i = 0; i < sps.size(); i++)
    {
        auto sp = sps[i];
        auto coeff = coeffs[i];
        cout << "(" << sp.row << ", " << sp.col << ") --> " << coeff << endl;
    }
    cout << "rhs: " << endl;
    for (auto r : rhs)
        cout << r << endl;

    cout << "sparse kkt system" << endl;
    KKTSystem kkt;
    auto x_kkt = kkt.variable({1}, {3});
    auto y_kkt = kkt.variable({0}, {0});
    auto z_kkt = kkt.variable({2}, {1});
    kkt.add_constraint(x_kkt + y_kkt + z_kkt, {0});
    kkt.add_constraint(x_kkt + y_kkt, {0});
    kkt.solver("mumps");
    sps = kkt.get_sparsity();
    coeffs = kkt.eval_coeffs();
    rhs = kkt.eval_rhs();
    // print the sparsity pattern and value
    cout << "Sparsity pattern: " << endl;
    for(size_t i = 0; i < sps.size(); i++)
    {
        auto sp = sps[i];
        auto coeff = coeffs[i];
        cout << "(" << sp.row << ", " << sp.col << ") --> " << coeff << endl;
    }

    cout << "rhs: " << endl;
    for(auto r : rhs)
        cout << r << endl;

    vector<double> res;
    kkt.solve(res);
    cout << kkt.evaluator(x_kkt).Eval(res).at(0) << endl;
    cout << "found solution: " << endl;
    for(auto r : res)
        cout << r << endl;
    // compute and print the residu
    cout << "residu: " << endl;
    vector<double> residu;
    kkt.residu(res, residu);
    for(auto r : residu)
        cout << r << endl;
    return 0;
}