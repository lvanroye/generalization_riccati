#ifndef INTERFACEMUMPSINCLUDED
#define INTERFACEMUMPSINCLUDED
#include "sparse_solver_interface.hpp"
#include "../expressions.hpp"
extern "C"
{
#include "dmumps_c.h"
#include "mpi.h"
}
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

namespace symspals
{
    class InterfaceMUMPS : public SparseSolverInterface
    {
    public:
        InterfaceMUMPS(const int dim, const vector<Index> &index):dim_(dim), nnz_(index.size())
        {
            // fortran indexing convention
            for(auto index : index)
            {
                ai.push_back(index.row + 1);
                aj.push_back(index.col + 1);
                a.push_back(0.0);
            }
        };
        void preprocess()
        {
            MUMPS_INT n = dim_;
            MUMPS_INT8 nnz = nnz_;
            // mumps takes upper part of matrix!
            MUMPS_INT *irn = aj.data();
            MUMPS_INT *jcn = ai.data();

            int argc = 1;
            char name[] = "InterfaceMUMPS";
            char* name_p = (char*) name;
            char **argv;
            argv = &name_p;
            ierr = MPI_Init(&argc, &argv);
            ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

            /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
            id.comm_fortran = USE_COMM_WORLD;
            id.par = 1;
            id.sym = 2;

            id.job = JOB_INIT;
            dmumps_c(&id);
            /* Define the problem on the host */
            if (myid == 0)
            {
                id.n = n;
                id.nnz = nnz;
                id.irn = irn;
                id.jcn = jcn;
            }
#define ICNTL(I) icntl[(I)-1] /* macro sdsdfs.t. indices match documentation */
            /* No outputs */
            id.ICNTL(1) = 6; // std output stream printing
            id.ICNTL(2) = 6;
            id.ICNTL(3) = 6;
            id.ICNTL(4) = 0; // printing level
            // id.ICNTL(10) = 0; // max no iterative refinement steps
            id.job = 1;
            dmumps_c(&id);
            if (id.infog[0] < 0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                       myid, id.infog[0], id.infog[1]);
                error = 1;
            }
        };
        void set_coefficient_matrix(const vector<double>& A)
        {
            // MUMPS_INT n = dim_;
            MUMPS_INT8 nnz = nnz_;
            // // mumps takes upper part of matrix!
            MUMPS_INT *irn = aj.data();
            MUMPS_INT *jcn = ai.data();
            id.irn = irn;
            id.jcn = jcn;
            for (int i = 0; i < nnz; i++)
            {
                a.at(i) = A.at(i);
            }
        }
        void solve(vector<double> &rhsvec)
        {
            MUMPS_INT *irn = aj.data();
            MUMPS_INT *jcn = ai.data();
            id.irn = irn;
            id.jcn = jcn;
            id.a = a.data();
            id.rhs = rhsvec.data();

            id.job = 2;
            dmumps_c(&id);
            if (id.infog[0] < 0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                       myid, id.infog[0], id.infog[1]);
                error = 1;
            }
            id.job = 3;
            dmumps_c(&id);
            if (id.infog[0] < 0)
            {
                printf(" (PROC %d) ERROR RETURN: \tINFOG(1)= %d\n\t\t\t\tINFOG(2)= %d\n",
                       myid, id.infog[0], id.infog[1]);
                error = 1;
            }
        }
        ~InterfaceMUMPS()
        {
            /* Terminate instance. */
            id.job = JOB_END;
            dmumps_c(&id);
            if (myid == 0)
            {
                if (!error)
                {
                    // printf("Solution is : (%8.2f  %8.2f)\n", rhs[0], rhs[1]);
                }
                else
                {
                    printf("An error has occured, please check error code returned by MUMPS.\n");
                }
            }
            free(id.wk_user);
        }
        DMUMPS_STRUC_C id;
        int myid, ierr;
        int error = 0;
        vector<double> a;
        vector<int> ai;
        vector<int> aj;
        const int dim_;
        const int nnz_;
    };
} // namespace fatrop
#endif //INTERFACEMUMPSINCLUDED