#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparse_solver_interface.hpp"
#include "../common.hpp"
#include <iostream>

using namespace std;

/* PARDISO prototype. */
extern "C" void pardisoinit(void *, int *, int *, int *, double *, int *);
extern "C" void pardiso(void *, int *, int *, int *, int *, int *,
                        double *, int *, int *, int *, int *, int *,
                        int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix(int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec(int *, int *, double *, int *);
extern "C" void pardiso_printstats(int *, int *, double *, int *, int *, int *,
                                   double *, int *);
namespace symspals
{
    class InterfacePardiso : public SparseSolverInterface
    {
    public:
        InterfacePardiso(const int dim, const vector<Index> &ind, const vector<double> &initial_coeff_est) : dim_(dim), nnz_(ind.size())
        {
            // !! assume ind is sorted by column major order and diagonal elements either zero or present
            /* -------------------------------------------------------------------- */
            /* ..  Setup Pardiso control parameters.                                */
            /* -------------------------------------------------------------------- */
            // print all triplets
            ia.resize(dim_ + 1);
            a.resize(nnz_);
            ja.resize(nnz_);
            int curr_col = -1;
            int curr_no_els = 0;
            for (size_t i = 0; i < ind.size(); i++)
            {
                while(curr_col != ind.at(i).col)
                {
                    // move to the next row
                    curr_col++;
                    if (curr_col != ind.at(i).col)
                    {
                        throw runtime_error("Pardiso requires column major order");
                    }
                    ia.at(curr_col) = (curr_no_els + 1);
                }
                ja.at(i) = (ind.at(i).row + 1);
                a.at(i) = initial_coeff_est.at(i);
                curr_no_els++;
            }
            ia.at(dim_) = (nnz_ + 1);

            error = 100;
            solver = 0; /* use sparse direct solver */
            for (int i = 0; i < 64; i++)
            {
                iparm[i] = 0;
            }
            pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
            // pardisoinit(pt, &mtype, iparm, dparm);

            if (error != 0)
            {
                cout << "error pardiso init " << error << endl;
                if (error == -10)
                    printf("No license file found \n");
                if (error == -11)
                    printf("License is expired \n");
                if (error == -12)
                    printf("Wrong username or hostname \n");
            }
            // else
            //     printf("[PARDISO]: License check was successful ... \n");

            /* Numbers of processors, value of OMP_NUM_THREADS */
            char *var = getenv("OMP_NUM_THREADS");
            if (var != NULL)
                sscanf(var, "%d", &num_procs);
            else
            {
                printf("Set environment OMP_NUM_THREADS to 1");
                exit(1);
            }
            iparm[2] = num_procs;
            // iparm[1] = 2; // METIS ordering

            // iparm[9] = 10; // no iterative refinement
            // iparm[7] = 0; /* Max numbers of iterative refinement steps. */
            // // iparm[10] = 2;
            iparm[12] = 1; // this is required for finding solutions
            // iparm[20] = 3; // Bunch-Kaufman pivoting + rescaling computed during factorization process

            // these are the default options for IPOPT

            // iparm[1] = 2;
            // // iparm[2] = num_procs; // Set the number of processors
            // // iparm[5] = 1;         // Overwrite right-hand side
            // iparm[7] = 0;

            // // // Options suggested by Olaf Schenk
            // iparm[9] = 12;
            // iparm[10] = 1; // Results in better scaling
            // // Matching information:  iparm[12] = 1 seems ok, but results in a
            // // large number of pivot perturbation
            // // Matching information:  iparm[12] = 2 robust,  but more  expensive method
            // iparm[12] = 2;

            // iparm[20] = 3; // Results in better accuracy
            // iparm[23] = 0; // parallel fac
            // iparm[24] = 0; // parallel solve
            // iparm[28] = 0; // 64-bit factorization (double precision)
            // // iparm[29] = 80; // we need this for IPOPT interface
            // //                  // iparm[33] = 1; // bit-by-bit identical results in parallel run
            // dparm[1] = 1e-6;
            // dparm[2] = 5000;
            // dparm[3] = 10;
            // dparm[4] = 0.5;
            // dparm[5] = 1e-1;
            // dparm[6] = 1e7;
            // dparm[7] = 5*1e6;

            maxfct = 1; /* Maximum number of numerical factorizations.  */
            mnum = 1;   /* Which factorization to use. */

            msglvl = 0; /* Print statistical information  */
            error = 0;  /* Initialize error flag */

            /* -------------------------------------------------------------------- */
            /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
            /*     notation.                                                        */
            /* -------------------------------------------------------------------- */

            // fortran indexing convention
            // ai = ai + 1;
            // aj = aj + 1;
            // SetCoefficientMatrix(tripl);

            /* -------------------------------------------------------------------- */
            /*  .. pardiso_chk_matrix(...)                                          */
            /*     Checks the consistency of the given matrix.                      */
            /*     Use this functionality only for debugging purposes               */
            /* -------------------------------------------------------------------- */
            int dimi = dim_;

            pardiso_chkmatrix(&mtype, &dimi, a.data(), ia.data(), ja.data(), &error);
            if (error != 0)
            {
                printf("\nERROR in consistency of matrix: %d", error);
                exit(1);
            }

            /* -------------------------------------------------------------------- */
            /* ..  pardiso_chkvec(...)                                              */
            /*     Checks the given vectors for infinite and NaN values             */
            /*     Input parameters (see PARDISO user manual for a description):    */
            /*     Use this functionality only for debugging purposes               */
            /* -------------------------------------------------------------------- */

            // pardiso_chkvec(&dimi, &nrhs, rhs, &error);
            // if (error != 0)
            // {
            //     printf("\nERROR  in right hand side: %d", error);
            //     exit(1);
            // }

            /* -------------------------------------------------------------------- */
            /* .. pardiso_printstats(...)                                           */
            /*    prints information on the matrix to STDOUT.                       */
            /*    Use this functionality only for debugging purposes                */
            /* -------------------------------------------------------------------- */

            // pardiso_printstats(&mtype, &n, a, ia, ja, &nrhs, b, &error);
            // if (error != 0)
            // {
            //     printf("\nERROR right hand side: %d", error);
            //     exit(1);
            // }
            // SetCoefficientMatrix(tripl);
        };
        void preprocess()
        {
            /* -------------------------------------------------------------------- */
            /* ..  Reordering and Symbolic Factorization.  This step also allocates */
            /*     all memory that is necessary for the factorization.              */
            /* -------------------------------------------------------------------- */
            phase = 11;
            int dimi = dim_;

            pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                    &dimi, a.data(), ia.data(), ja.data(), &idum, &nrhs,
                    iparm, &msglvl, &ddum, &ddum, &error, dparm);

            if (error != 0)
            {
                printf("\nERROR during symbolic factorization: %d", error);
                exit(1);
            }
            // printf("\nReordering completed ... ");
            // printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
            // printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
        };
        void set_coefficient_matrix(const vector<double> &A)
        {
            a = A;
        }
        void solve(vector<double> &rhsvec)
        {
            // double el_time = 0.0;

            /* -------------------------------------------------------------------- */
            /* ..  Back substitution and iterative refinement.                      */
            /* -------------------------------------------------------------------- */

            int dimi = dim_;
            vector<double> rhsvec_copy = rhsvec;

            phase = 22;
            pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                    &dimi, a.data(), ia.data(), ja.data(), &idum, &nrhs,
                    iparm, &msglvl, rhsvec_copy.data(), rhsvec.data(), &error, dparm);
            if (error != 0)
            {
                printf("\nERROR during factorization: %d", error);
                exit(3);
            }
            phase = 33;
            pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                    &dimi, a.data(), ia.data(), ja.data(), &idum, &nrhs,
                    iparm, &msglvl, rhsvec_copy.data(), rhsvec.data(), &error, dparm);

            if (error != 0)
            {
                printf("\nERROR during solution: %d", error);
                exit(3);
            }

        }
        ~InterfacePardiso()
        {
            /* -------------------------------------------------------------------- */
            /* ..  Termination and release of memory.                               */
            /* -------------------------------------------------------------------- */
            phase = -1; /* Release internal memory. */
            int dimi_ = dim_;

            pardiso(pt, &maxfct, &mnum, &mtype, &phase,
                    &dimi_, &ddum, ia.data(), ja.data(), &idum, &nrhs,
                    iparm, &msglvl, &ddum, &ddum, &error, dparm);
        }
        vector<double> a;
        vector<int> ia;
        vector<int> ja;
        vector<int> permutation;

    private:
        int mtype = -2; /* Real symmetric matrix */

        /* RHS and solution vectors. */
        int nrhs = 1; /* Number of right hand sides. */

        /* Internal solver memory pointer pt,                  */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
        /* or void *pt[64] should be OK on both architectures  */
        void *pt[64];

        /* Pardiso control parameters. */
        int iparm[64];
        double dparm[64];
        int maxfct, mnum, phase, error, msglvl, solver;

        /* Number of processors. */
        int num_procs;

        /* Auxiliary variables. */
        // char *var;
        int i;

        double ddum; /* Double dummy */
        int idum;    /* Integer dummy. */
        int dim_;
        int nnz_;
    };
} // namespace symspals