#include "benchmark_helpers.hpp"
int main()
{
    vector<double> solve_time;
    vector<double> residu;
    BenchmarkHelperCP benchmark_helper_cp;
    benchmark_helper_cp.gen_riccati(true, solve_time, residu);
    // benchmark_helper_cp.sparse_solver("ma57", solve_time, residu);
    // print the residu results
    cout << "residu's: " << endl;
    for (double r : residu)
        cout << r << endl;
    // print the timings
    cout << "timings: " << endl;
    for (double t : solve_time)
        cout << t << endl;
    return 0;
}