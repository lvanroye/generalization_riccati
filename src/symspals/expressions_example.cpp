#include <iostream>
#include "expressions.hpp"
using namespace SpaLS;
int main()
{
    Sym a("a");
    Sym b("b");
    Sym c("c");
    Sym x("x");
    Sym y("y");

    auto expr = (a + b + 1) * (x + y) - 5 * x - 0*x + (y+y);
    auto terms = GetTerms(expr, ExprVec({x, y}));
    // terms
    cout << "terms: "<< endl;
    for (auto term : terms)
        cout << term << endl;
    auto coeffs = GetCoefficients({expr}, {x, y});
    // coefficients
    cout << "coefficients: "<< endl;
    for (auto coeff : coeffs)
        // cout << coeff.row << endl;
        // cout << coeff.col << endl;
        cout << coeff.value << endl;
    return 0;
    // linear_system.add_equation(-X.at(k+1) + A.at(K) * X.at(k) + B.at(k) * U.at(k), 0);
};
