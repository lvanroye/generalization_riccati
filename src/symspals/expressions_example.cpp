#include <iostream>
#include "expressions.hpp"
using namespace symspals;
int main()
{
    Sym a("a");
    Sym b("b");
    Sym c("c");
    Sym x("x");
    Sym y("y");

    // auto expr = (a + b + 1) * (x + y) - 5 * x - 0*x + (y+y);
    auto expr = a*x + b*y;
    auto expr2 = -x + -y;
    auto terms = GetTerms(expr, ExprVec({x, y}));
    // terms
    cout << "terms: "<< endl;
    for (auto term : terms)
        cout << term << endl;
    auto coeffs = GetCoefficients({expr, expr2}, {x, y}, {a, b});
    // auto coeffs = CoefficientSimple()({expr, expr2}, {x, y});
    // coefficients
    cout << "coefficients: "<< endl;
    for (auto coeff : coeffs)
        cout << coeff.value << endl;

    
    auto coeffs_v = Function({a, b, c}, {coeffs.at(0).value, coeffs.at(1).value, coeffs.at(2).value, coeffs.at(3).value}).Eval({1, 2, 3});
    for (auto coeff: coeffs_v)
    {
        cout << coeff << endl;
    }
    return 0;
};
