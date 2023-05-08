#pragma once
#include "random_problem/random_ocp.hpp"
#include "cart_pendulum_problem/cart_pendulum.hpp"
#include "fatrop_problem_wrap.hpp"
#include "benchmark_sparse.hpp"
#include "benchmark_gen_riccati.hpp"
#include <memory>
using namespace genriccati_benchmark;
using namespace std;
class BenchmarkHelperRandom
{
public:
    void gen_riccati(int K, int nx, int nu, int ne_init, int ne_middle, int ne_final, bool it_ref, bool incr_acc, vector<double> &res_time, vector<double> &res_acc)
    {
        shared_ptr<RandomOCP> ocp = make_shared<RandomOCP>(K, nx, nu, ne_init, ne_middle, ne_final);
        FatropProblemWrap fatrop_problem_wrap(ocp);
        auto cocp = fatrop_problem_wrap.create_cocp();
        BenchmarkGenRiccati benchmark_gen_riccati(cocp);
        for (int i = 0; i < 10; i++)
        {
            ocp->set_seed(i);
            // fatrop_problem_wrap.least_squares_dual();
            fatrop_problem_wrap.populate_cocp(cocp);
            benchmark_gen_riccati.ocp_ls_riccati.it_ref = it_ref;
            benchmark_gen_riccati.solve(cocp);
            double time;
            double residu;
            benchmark_gen_riccati.last_stats(time, residu);
            res_time.push_back(time);
            res_acc.push_back(residu);
        }
    }
    void sparse_solver(const string &solver_name, int K, int nx, int nu, int ne_init, int ne_middle, int ne_final, vector<double> &res_time, vector<double> &res_acc)
    {
        shared_ptr<RandomOCP> ocp = make_shared<RandomOCP>(K, nx, nu, ne_init, ne_middle, ne_final);
        FatropProblemWrap fatrop_problem_wrap(ocp);
        auto cocp = fatrop_problem_wrap.create_cocp();
        BenchmarkSparse benchmark_sparse(cocp);
        vector<double> sol;
        benchmark_sparse.kkt_system_.solver(solver_name);
        for (int i = 0; i < 10; i++)
        {
            ocp->set_seed(i);
            // fatrop_problem_wrap.least_squares_dual();
            fatrop_problem_wrap.populate_cocp(cocp);
            benchmark_sparse.set_value(cocp);
            benchmark_sparse.kkt_system_.solve(sol);
            res_time.push_back(benchmark_sparse.kkt_system_.last_time());
            vector<double> residu;
            benchmark_sparse.kkt_system_.residu(sol, residu);
            // get inf norm of residu
            double inf_norm = residu[0];
            for (double r : residu)
                inf_norm = std::max(inf_norm, std::abs(r));
            res_acc.push_back(inf_norm);
        }
    }
};
class BenchmarkHelperCP
{
public:
    void gen_riccati(bool it_ref, vector<double> &res_time, vector<double> &res_acc)
    {
        shared_ptr<OCPAbstract> ocp = make_shared<StageOCP>(CartPendulumProblem());
        FatropProblemWrap fatrop_problem_wrap(ocp);
        auto cocp = fatrop_problem_wrap.create_cocp();
        // fatrop_problem_wrap.least_squares_dual();
        BenchmarkGenRiccati benchmark_gen_riccati(cocp);
        for (int i = 0; i < 10; i++)
        {
            fatrop_problem_wrap.populate_cocp(cocp);
            benchmark_gen_riccati.ocp_ls_riccati.it_ref = it_ref;
            benchmark_gen_riccati.solve(cocp);
            double time;
            double residu;
            benchmark_gen_riccati.last_stats(time, residu);
            res_time.push_back(time);
            res_acc.push_back(residu);
        }
    }

    void sparse_solver(const string &solver_name, vector<double> &res_time, vector<double> &res_acc)
    {
        shared_ptr<OCPAbstract> ocp = make_shared<StageOCP>(CartPendulumProblem());
        FatropProblemWrap fatrop_problem_wrap(ocp);
        auto cocp = fatrop_problem_wrap.create_cocp();
        // fatrop_problem_wrap.least_squares_dual();
        BenchmarkSparse benchmark_sparse(cocp);
        benchmark_sparse.kkt_system_.solver(solver_name);
        fatrop_problem_wrap.populate_cocp(cocp);
        benchmark_sparse.set_value(cocp);
        benchmark_sparse.kkt_system_.numeric_prune();
        for (int i = 0; i < 10; i++)
        {
            vector<double> sol;
            fatrop_problem_wrap.populate_cocp(cocp);
            benchmark_sparse.set_value(cocp);
            benchmark_sparse.kkt_system_.solve(sol);
            res_time.push_back(benchmark_sparse.kkt_system_.last_time());
            vector<double> residu(sol.size(), 0);
            benchmark_sparse.kkt_system_.residu(sol, residu);
            // get inf norm of residu
            double inf_norm = residu[0];
            for (double r : residu)
                inf_norm = std::max(inf_norm, std::abs(r));
            res_acc.push_back(inf_norm);
        }
        };
    };