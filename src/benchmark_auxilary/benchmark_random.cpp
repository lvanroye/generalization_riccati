
#include "benchmark_helpers.hpp"
void print_results(vector<double> &solve_time, vector<double> &residu, const string &solver_name)
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
double mean(vector<double> &vec)
{
    double sum = 0;
    for (double v : vec)
        sum += v;
    return sum / vec.size();
}
double std_dev(vector<double> &vec)
{
    double m = mean(vec);
    double sum = 0;
    for (double v : vec)
        sum += (v - m) * (v - m);
    return sqrt(sum / vec.size());
}
void process_results(vector<double> &solve_time, vector<double> &residu, vector<double> &solve_time_mean, vector<double> &residu_mean, vector<double> &solve_time_std, vector<double> &residu_std)
{
    solve_time_mean.push_back(mean(solve_time));
    residu_mean.push_back(mean(residu));
    solve_time_std.push_back(std_dev(solve_time));
    residu_std.push_back(std_dev(residu));
    solve_time.clear();
    residu.clear();
}
void print_header(const string &var, const vector<string> &solver)
{
    cout << var << " & ";
    for (string s : solver)
    {
        cout << s << " & ";
    }
    cout << "\\\\" << endl;
}
void print_line(const int var, const vector<double> &quantity)
{
    cout << var << " & ";
    for (int i = 0; i < quantity.size(); i++)
    {
        cout << quantity[i] << " & ";
    }
    cout << "\\\\" << endl;
}
void print_matrix(const vector<int> var, const vector<vector<double>> &quantities)
{
    for (int i = 0; i < var.size(); i++)
    {
        print_line(var[i], quantities[i]);
    }
}
int main()
{

    BenchmarkHelperRandom benchmark_helper_cp;
    {
        print_header("K", {"Riccati with iterative refinement", "Riccati without iterative refinement", "Sparse solver ma57", "Sparse solver mumps", "Sparse solver pardiso"});
        vector<vector<double>> solve_time_mean_matrix;
        vector<vector<double>> residu_mean_matrix;
        vector<vector<double>> solve_time_std_matrix;
        vector<vector<double>> residu_std_matrix;
        vector<vector<double>> solve_time_matrix;
        vector<int> Karr = {10, 16, 27, 46, 77, 129, 215, 359, 599, 1000};
        for (int K : Karr)
        {
            vector<double> solve_time_mean;
            vector<double> residu_mean;
            vector<double> solve_time_std;
            vector<double> residu_std;
            vector<double> solve_time;
            vector<double> residu;
            // int K = 100;
            int nx = 10;
            int nu = 5;
            int ne_init = 10;
            int ne_middle = 0;
            int ne_final = 10;
            benchmark_helper_cp.gen_riccati(K, nx, nu, ne_init, ne_middle, ne_final, true, true, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.gen_riccati(K, nx, nu, ne_init, ne_middle, ne_final, false, true, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.sparse_solver("ma57", K, nx, nu, ne_init, ne_middle, ne_final, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.sparse_solver("mumps", K, nx, nu, ne_init, ne_middle, ne_final, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.sparse_solver("pardiso", K, nx, nu, ne_init, ne_middle, ne_final, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            // add to matrix
            solve_time_mean_matrix.push_back(solve_time_mean);
            residu_mean_matrix.push_back(residu_mean);
            solve_time_std_matrix.push_back(solve_time_std);
            residu_std_matrix.push_back(residu_std);
        }
        // print the results for the mean time
        cout << "mean time" << endl;
        print_header("K", {"Riccati with iterative refinement", "Riccati without iterative refinement", "Sparse solver ma57", "Sparse solver mumps", "Sparse solver pardiso"});
        print_matrix(Karr, solve_time_mean_matrix);
        cout << "mean residu" << endl;
        print_header("K", {"Riccati with iterative refinement", "Riccati without iterative refinement", "Sparse solver ma57", "Sparse solver mumps", "Sparse solver pardiso"});
        print_matrix(Karr, residu_mean_matrix);
    }
    {
        vector<vector<double>> solve_time_mean_matrix;
        vector<vector<double>> residu_mean_matrix;
        vector<vector<double>> solve_time_std_matrix;
        vector<vector<double>> residu_std_matrix;
        vector<vector<double>> solve_time_matrix;
        vector<int> nxarr = {4, 8, 12, 16, 24, 36, 54, 80, 120, 180};
        for (int nx : nxarr)
        {
            vector<double> solve_time_mean;
            vector<double> residu_mean;
            vector<double> solve_time_std;
            vector<double> residu_std;
            vector<double> solve_time;
            vector<double> residu;
            int K = 100;
            // int nx = 10;
            int nu = nx/2;
            int ne_init = nx/4;
            int ne_middle = nx/4;
            int ne_final = nx/4;
            benchmark_helper_cp.gen_riccati(K, nx, nu, ne_init, ne_middle, ne_final, true, true, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.gen_riccati(K, nx, nu, ne_init, ne_middle, ne_final, false, true, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.sparse_solver("ma57", K, nx, nu, ne_init, ne_middle, ne_final, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.sparse_solver("mumps", K, nx, nu, ne_init, ne_middle, ne_final, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            benchmark_helper_cp.sparse_solver("pardiso", K, nx, nu, ne_init, ne_middle, ne_final, solve_time, residu);
            process_results(solve_time, residu, solve_time_mean, residu_mean, solve_time_std, residu_std);
            // add to matrix
            solve_time_mean_matrix.push_back(solve_time_mean);
            residu_mean_matrix.push_back(residu_mean);
            solve_time_std_matrix.push_back(solve_time_std);
            residu_std_matrix.push_back(residu_std);
        }
        // print the results for the mean time
        cout << "mean time" << endl;
        print_header("nx", {"Riccati with iterative refinement", "Riccati without iterative refinement", "Sparse solver ma57", "Sparse solver mumps", "Sparse solver pardiso"});
        print_matrix(nxarr, solve_time_mean_matrix);
        cout << "mean residu" << endl;
        print_header("nx", {"Riccati with iterative refinement", "Riccati without iterative refinement", "Sparse solver ma57", "Sparse solver mumps", "Sparse solver pardiso"});
        print_matrix(nxarr, residu_mean_matrix);

    }
    return 0;
}