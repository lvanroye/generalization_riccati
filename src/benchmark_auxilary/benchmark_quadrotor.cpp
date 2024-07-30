/* Fatrop - A fast trajectory optimization solver
 Copyright (C) 2022 - 2024 Lander Vanroye, KU Leuven. All rights reserved.

This file is part of Fatrop.

Fatrop is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fatrop is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Fatrop.  If not, see <http://www.gnu.org/licenses/>.
*/
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