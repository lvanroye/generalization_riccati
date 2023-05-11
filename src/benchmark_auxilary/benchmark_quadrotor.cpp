
#include "benchmark_helpers.hpp"
void print_results(vector<double> &solve_time, vector<double> &residu, const string& solver_name)
{
    cout << "solver: " << solver_name << endl;
    cout << "residu's: " << endl;
    for (double r : residu)
        cout << r << endl;
    cout << "timings: " << endl;
    for (double t : solve_time)
        cout << t << endl;
    solve_time.clear();
    residu.clear();
}
int main()
{
    vector<double> solve_time;
    vector<double> residu;
    BenchmarkHelperQuadrotor benchmark_helper_cp;
    benchmark_helper_cp.gen_riccati(true, solve_time, residu);
    print_results(solve_time, residu, "Riccati with iterative refinement");
    benchmark_helper_cp.gen_riccati(false, solve_time, residu);
    print_results(solve_time, residu, "Riccati without iterative refinement");
    benchmark_helper_cp.sparse_solver("ma57", solve_time, residu);
    print_results(solve_time, residu, "Sparse solver ma57");
    benchmark_helper_cp.sparse_solver("mumps", solve_time, residu);
    print_results(solve_time, residu, "Sparse solver mumps");
    benchmark_helper_cp.sparse_solver("pardiso", solve_time, residu);
    print_results(solve_time, residu, "Sparse solver pardiso");
    return 0;
}