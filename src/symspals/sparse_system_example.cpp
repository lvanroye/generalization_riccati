#include "sparse_system.hpp"
using namespace symspals;
int main()
{
    cout << "sparse kkt system" << endl;
    /* 
    specify the kkt system 
    */
    KKTSystem kkt;
    auto x_kkt = kkt.variable({69}, {3});
    auto y_kkt = kkt.variable({0}, {0});
    auto z_kkt = kkt.variable({2}, {1});
    auto p = kkt.parameter(1, 1);
    kkt.add_constraint(- y_kkt + 10*x_kkt  + (*p)*z_kkt, {0});
    kkt.add_constraint(-x_kkt + y_kkt, {0});
    /* 
    choose a solver
    */
    kkt.solver("pardiso");
    /* 
    set the parameters
    */
    kkt.set_value(p, {123});
    /* 
    solve the system
    */
    vector<double> res;
    kkt.solve(res);
    cout << kkt.evaluator(x_kkt).Eval(res).at(0) << endl;
    cout << "found solution: " << endl;
    for(auto r : res)
        cout << r << endl;
    /* 
    get some information about the system
    */
    auto sps = kkt.get_sparsity();
    auto coeffs = kkt.eval_coeffs();
    auto rhs = kkt.eval_rhs();
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
    /* 
    compute the residual
    */
    cout << "residu: " << endl;
    vector<double> residu;
    kkt.residu(res, residu);
    for(auto r : residu)
        cout << r << endl;
    return 0;
}