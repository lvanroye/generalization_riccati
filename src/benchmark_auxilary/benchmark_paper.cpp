#include "benchmark_helpers.hpp"
int main()
{
    vector<double> solve_time;
    vector<double> residu;
    // BenchmarkHelperQuadrotor benchmark_helper_cp;
    BenchmarkHelperRandom benchmark_helper_rand;
    for(int i =0; i<100; i++)
    benchmark_helper_rand.gen_riccati(50, 50, 25, 10, 10, 10, true, false, solve_time, residu);
    // benchmark_helper_cp.gen_riccati(true, solve_time, residu);
    // benchmark_helper_cp.sparse_solver("ma57", solve_time, residu);
    // benchmark_helper_rand.sparse_solver("ma57",50, 10, 5, 10, 0, 10, solve_time, residu);
    cout << "residu's: " << endl;
    for (double r : residu)
        cout << r << endl;
    cout << "timings: " << endl;
    for (double t : solve_time)
        cout << t << endl;
    return 0;
}