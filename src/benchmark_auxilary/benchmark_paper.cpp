#include "random_problem/random_ocp.hpp"
#include "cart_pendulum_problem/cart_pendulum.hpp"
#include "fatrop_problem_wrap.hpp"
#include "benchmark_sparse.hpp"
#include <memory>
using namespace genriccati_benchmark;
using namespace std;
int main()
{
    // create a random OCP
    // shared_ptr<OCPAbstract> ocp = make_shared<RandomOCP>(100, 50, 25, 50, 10, 50);
    // creat a cart pendulum OCP
    shared_ptr<OCPAbstract> ocp = make_shared<StageOCP>(CartPendulumProblem());
    FatropProblemWrap fatrop_problem_wrap(ocp);
    auto cocp = fatrop_problem_wrap.create_cocp();
    fatrop_problem_wrap.populate_cocp(cocp);
    BenchmarkSparse benchmark_sparse(cocp);
    vector<double> sol;
    // cout << "solve 1" << endl;
    benchmark_sparse.kkt_system_.numeric_prune();
    benchmark_sparse.kkt_system_.print_coeff_values();

    benchmark_sparse.kkt_system_.solve(sol);
    fatrop_problem_wrap.controle();
    auto ux0 = benchmark_sparse.kkt_system_.evaluator(benchmark_sparse.ux_syms[0]).Eval(sol);
    // print ux0 with for loop
    cout << "ux0 = [";
    for (size_t i = 0; i < ux0.size(); i++)
    {
        cout << ux0[i] << " ";
    }
    cout << "]" << endl;
    return 0;
}