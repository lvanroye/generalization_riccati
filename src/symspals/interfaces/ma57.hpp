#pragma once
#include "sparse_solver_interface.hpp"
#include "../common.hpp"
#include <cmath>
#include <memory>
#include <vector>
#include <iostream>

using namespace std;
namespace symspals
{
#ifdef __cplusplus
    extern "C"
    {
#endif
        extern void ma57id_(double cntl[5], int icntl[20]);
        extern void ma57ad_(int *n, int *nz, int *irn, int *jcn,
                            int *lkeep, int *keep, int *iw,
                            int icntl[20], int info[40], double rinfo[20]);
        extern void ma57bd_(int *n, int *nz, double *a, double *fact,
                            int *lfact, int *ifact, int *lifact, int *lkeep,
                            int *keep, int *iw, int icntl[20], double cntl[5],
                            int info[40], double rinfo[20]);
        extern void ma57cd_(int *job, int *n, double *fact, int *lfact,
                            int *ifact, int *lifact, int *nrhs, double *rhs,
                            int *lrhs, double *w, int *lw, int *iw,
                            int icntl[20], int info[40]);
#ifdef __cplusplus
    }
#endif
    template <typename T_Vec, typename T_Mat, class aSolver>
    class linSolve
    {
    public:
        linSolve(const int n, const int nz, T_Mat A)
        {
            theSolver = new aSolver(n, nz, A);
        }
        void Factorize() { theSolver->Factorize(); }
        void Solve(T_Vec rhs) { theSolver->Solve(rhs); }
        T_Vec getSol() { return theSolver->getSol(); }

    private:
        unique_ptr<aSolver> theSolver;
    };
    class MA57
    {
    public:
        MA57(const int i_n, const vector<Index> &index)
        {
            n = i_n;
            nz = index.size();

            irn.resize(nz);
            jcn.resize(nz);

            a.resize(nz);
            for (int i = 0; i < nz; i++)
            {
                // !! ma57 requires upper matrix !!
                irn[i] = index.at(i).col + 1;
                jcn[i] = index.at(i).row + 1;
                a[i] = 0.0;
            }

            lkeep = (nz > n) ? (5 * n + 2 * nz + 42) : (6 * n + nz + 42);
            keep.resize(lkeep);
            iw.resize(5 * n);

            job = 1;

            Initialize();
            // number of iterative refinements
            icntl[8] = 0;
            ///// these are the default settings of IPOPT
            /* Custom settings for MA57. */
            icntl[0] = 0; /* Error stream */
            icntl[1] = 0; /* Warning stream. */

            icntl[3] = 1; /* Print statistics.  NOT Used. */
            icntl[4] = 0; /* Print error. */

            // icntl[14] = 0; // turn off automatic scaling
            Analyze();
            int nrhs = 1;            // number of right hand side being solved
            int lw = 1.2 * n * nrhs; // length of w; lw>=n*nrhs
            w.resize(lw);            // double workspace
            lfact = 1.2 * info[8];
            fact.resize(lfact);
            lifact = 1.2 * info[9];
            ifact.resize(lifact);
        }
        void Factorize()
        {
            ma57bd_(&n, &nz, a.data(), fact.data(), &lfact, ifact.data(), &lifact, &lkeep,
                    keep.data(), iw.data(), icntl, cntl, info, rinfo);
            // cout << ">>>>>>>>>>>>>>> finished factoring using ma57" << endl;
            if (info[0] != 0)
            {
                cout << "factorize ma57 error " << info[0] << endl;
            }
            // assert(info[0] == 0);
        }
        void SetCoefficientMatrix(const vector<double> &tripl)
        {
            a = tripl;
        }

        void Solve(double *rhs)
        {
            int nrhs = 1;            // number of right hand side being solved
            int lw = 1.2 * n * nrhs; // length of w; lw>=n*nrhs
            // double *w = new double[lw]; // double workspace
            int lrhs = n; // integer, lenght of rhs
            // cout << ">>>>>>>>>>>>>>> started solving using ma57" << endl;
            ma57cd_(&job, &n, fact.data(), &lfact, ifact.data(), &lifact, &nrhs, rhs,
                    &lrhs, w.data(), &lw, iw.data(), icntl, info);
            if (info[0] != 0)
            {
                cout << "solve ma57 error " << info[0] << endl;
            }
            cout << "ma57 number of delayed pivos " << info[22] << endl;
            // cout << ">>>>>>>>>>>>>>> finished solving using ma57" << endl;
            // delete[] w;
        }

    protected:
        int n;  // order of matrix, namely number of rows or cols
        int nz; // number of nonzeros entries

        double cntl[5]; // double array of length 5; contains double control values
        int icntl[20];  // integer array of length 20; contaitns integer control values

        vector<int> irn; // row index of input
        vector<int> jcn; // col index of input

        int lkeep;        // length of keep; lkeep ≥  5*n+nz+MAX(n,nz)+42
        vector<int> keep; // integer workspace of length lkeep
        vector<int> iw;   // integer array of 5*n; pivot sequence
        int info[40];     // integer array of length 40
        double rinfo[20]; // double array of length 20

        vector<double> a;    // double array of data
        vector<double> fact; // double array of length lfact; holds the entries of the factors
        int lfact;           // length of fact
        vector<int> ifact;   // integer array; integer indexing info on the matrix factors
        int lifact;          // length of ifact

        int job; // integer; job = 1 if solving system AX=B

        vector<double> sol; // Solution
        vector<double> w;   // double workspace

        // Internal functions
        void Initialize()
        {
            ma57id_(cntl, icntl);
            if (info[0] != 0)
            {
                cout << "init ma57 error " << info[0] << endl;
            }
            // icntl[4] = 4;
        }
        void Analyze()
        {
            ma57ad_(&n, &nz, irn.data(), jcn.data(), &lkeep, keep.data(), iw.data(),
                    icntl, info, rinfo);
            if (info[0] != 0)
            {
                cout << "analyze ma57 error " << info[0] << endl;
            }
        }
    };

    class InterfaceMA57 : public SparseSolverInterface
    {
    public:
        InterfaceMA57(const int dim, const vector<Index> &index) : ma57_(dim, index){};
        void preprocess()
        {
        }
        void set_coefficient_matrix(const vector<double> &A)
        {
            ma57_.SetCoefficientMatrix(A);
        }
        void solve(vector<double> &rhsvec)
        {
            ma57_.Factorize();
            ma57_.Solve(rhsvec.data());
        }
        ~InterfaceMA57()
        {
        }

    private:
        MA57 ma57_;
    };
} // namespace symspals