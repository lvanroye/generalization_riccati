#include "blasfeo_wrapper.hpp"
namespace gen_riccati
{
    void fatrop_dcolsc(int kmax, double alpha, struct blasfeo_dmat *sA, int ai, int aj)
    {
        for (int k = 0; k < kmax; k++)
        {
            MATEL(sA, ai + k, aj) *= alpha;
        }
    }
    /** \brief copy elements from sx to sy but in reversed order to avoid aliasing issues in recursion */
    void fatrop_dveccp_reversed(int m, struct blasfeo_dvec *sx, int xi, struct blasfeo_dvec *sy, int yi)
    {
        for (int i = m - 1; i >= 0; i--)
        {
            VECEL(sy, yi + i) = VECEL(sx, xi + i);
        }
    }

    // void fatrop_potrf_l_mn(int m, int n, struct blasfeo_dmat *sC, int ci, int cj, struct blasfeo_dmat *sD, int di, int dj)
    // {
    //     blasfeo_dpotrf_l_mn(m, n, sC, ci, cj, sD, di, dj);
    //     int minmn = (m < n) ? m : n;
    //     for (int i =0; i<minmn; i++){
    //         assert(MATEL(sD, di+i,dj+i)>0);
    //     }
    // }
    /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    void fatrop_dtrsm_rlnn(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj)
    {
        sD->use_dA = 0;
        for (int aj = n - 1; aj >= 0; aj--)
        {
            double ajj = MATEL(sA, offs_ai + aj, aj + offs_aj);
            double inv_ajj = 1.0 / ajj;
            double scjj = alpha * inv_ajj;
            for (int k = 0; k < m; k++)
            {
                MATEL(sD, offs_di + k, offs_dj + aj) = scjj * MATEL(sB, offs_bi + k, offs_bj + aj);
            }
            for (int ai = aj + 1; ai < n; ai++)
            {
                double sc = -inv_ajj * MATEL(sA, offs_ai + ai, offs_aj + aj);
                for (int k = 0; k < m; k++)
                {
                    // this algorithm is "store bounded", the loops can be switched like in the alt version to make this more efficient
                    MATEL(sD, offs_di + k, offs_dj + aj) += sc * MATEL(sD, offs_di + k, offs_dj + ai);
                }
            }
        }
    }
    // /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    // // this is an experimental, more efficient, but the corner cases are not treated in unrolled loop!! it achieves a speed-up of about factor 3 w.r.t naive implementation
    void fatrop_dtrsm_rlnn_alt(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj)
    {
        for (int aj = n - 1; aj >= 0; aj--)
        {
            double ajj = MATEL(sA, offs_ai + aj, aj + offs_aj);
            double inv_ajj = 1.0 / ajj;
            double scjj = alpha * inv_ajj;
            for (int k = 0; k < m; k++)
            {
                // todo, check if possible to incude in main loop
                MATEL(sD, offs_di + k, offs_dj + aj) = scjj * MATEL(sB, offs_bi + k, offs_bj + aj);
            }
            for (int k = 0; k < m; k++)
            {
                double res = 0.0;
                // double res4 = 0.0;
                // double res5 = 0.0;
                // double res6 = 0.0;
                // double res7 = 0.0;
                for (int ai = aj + 1; ai < n; ai = ai + 1)
                {
                    // todo unroll loop -> more independent operations -> filled pipelines
                    double sc = -inv_ajj * MATEL(sA, offs_ai + ai, offs_aj + aj);
                    res += sc * MATEL(sD, offs_di + k, offs_dj + ai);
                }
                MATEL(sD, offs_di + k, offs_dj + aj) += res;
            }
        }
    }
    // /** \brief D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal, fatrop uses its own (naive) implementation since it  not implemented yet in blasfeo*/
    // // this is an experimental, more efficient, but the corner cases are not treated in unrolled loop!! it achieves a speed-up of about factor 3 w.r.t naive implementation
    // void fatrop_dtrsm_rlnn_alt(int m, int n, double alpha, MAT *sA, int offs_ai, int offs_aj, MAT *sB, int offs_bi, int offs_bj, MAT *sD, int offs_di, int offs_dj)
    // {
    //     for (int aj = n - 1; aj >= 0; aj--)
    //     {
    //         double ajj = MATEL(sA, offs_ai + aj, aj + offs_aj);
    //         double inv_ajj = 1.0 / ajj;
    //         double scjj = alpha * inv_ajj;
    //         for (int k = 0; k < m; k++)
    //         {
    //             // todo, check if possible to incude in main loop
    //             MATEL(sD, offs_di + k, offs_dj + aj) = scjj * MATEL(sB, offs_bi + k, offs_bj + aj);
    //         }
    //         for (int k = 0; k < m; k++)
    //         {
    //             double res = 0.0;
    //             double res1 = 0.0;
    //             double res2 = 0.0;
    //             double res3 = 0.0;
    //             // double res4 = 0.0;
    //             // double res5 = 0.0;
    //             // double res6 = 0.0;
    //             // double res7 = 0.0;
    //             for (int ai = aj + 1; ai < n; ai = ai + 4)
    //             {
    //                 // todo unroll loop -> more independent operations -> filled pipelines
    //                 double sc = -inv_ajj * MATEL(sA, offs_ai + ai, offs_aj + aj);
    //                 res += sc * MATEL(sD, offs_di + k, offs_dj + ai);
    //                 double sc1 = -inv_ajj * MATEL(sA, offs_ai + ai + 1, offs_aj + aj);
    //                 res1 += sc1 * MATEL(sD, offs_di + k, offs_dj + ai + 1);
    //                 double sc2 = -inv_ajj * MATEL(sA, offs_ai + ai + 2, offs_aj + aj);
    //                 res2 += sc2 * MATEL(sD, offs_di + k, offs_dj + ai + 2);
    //                 double sc3 = -inv_ajj * MATEL(sA, offs_ai + ai + 3, offs_aj + aj);
    //                 res3 += sc3 * MATEL(sD, offs_di + k, offs_dj + ai + 3);
    //                 // double sc4 = -inv_ajj * MATEL(sA, offs_ai + ai+4, offs_aj + aj);
    //                 // res4 += sc * MATEL(sD, offs_di + k, offs_dj + ai+4);
    //                 // double sc5 = -inv_ajj * MATEL(sA, offs_ai + ai+5, offs_aj + aj);
    //                 // res5 += sc * MATEL(sD, offs_di + k, offs_dj + ai+5);
    //                 // double sc6 = -inv_ajj * MATEL(sA, offs_ai + ai+6, offs_aj + aj);
    //                 // res6 += sc * MATEL(sD, offs_di + k, offs_dj + ai+6);
    //                 // double sc7 = -inv_ajj * MATEL(sA, offs_ai + ai+7, offs_aj + aj);
    //                 // res7 += sc * MATEL(sD, offs_di + k, offs_dj + ai+7);
    //             }
    //             // MATEL(sD, offs_di + k, offs_dj + aj) += (res+res1+res2+res3+res4+res5+res6+res7);
    //             MATEL(sD, offs_di + k, offs_dj + aj) += (res + res1) + (res2 + res3);
    //         }
    //     }
    // }
    // B <= B + alpha*A^T (B is mxn)
    void fatrop_dgead_transposed(int m, int n, double alpha, struct blasfeo_dmat *sA, int offs_ai, int offs_aj, struct blasfeo_dmat *sB, int offs_bi, int offs_bj)
    {
        for (int bj = 0; bj < n; bj++)
        {
            for (int bi = 0; bi < m; bi++)
            {
                MATEL(sB, offs_bi + bi, offs_bj + bj) += alpha * MATEL(sA, offs_ai + bj, offs_aj + bi);
            }
        }
    }

    /** \brief Returns the maximum element of a blasfeo matrix of size (m,n), starting at (ai,aj) */
    MatrixInd max_el(int m, int n, MAT *matr, int ai, int aj)
    {
        MatrixInd res;
        res.ai = ai;
        res.aj = aj;
        double valmax = 0.0;
        for (int j = aj; j < n; j++)
        {
            for (int i = ai; i < m; i++)
            {
                double valij = abs(MATEL(matr, i, j));
                if (valij >= valmax)
                {
                    valmax = valij;
                    res.ai = i;
                    res.aj = j;
                }
            }
        }
        return res;
    };
    /** \brief Function to calculate LU factorization result is saved in A, L is lower unitriangular */
    void LU_FACT(const int m, const int n, const int n_max, int &rank, MAT *A, PMAT *Pl_p, PMAT *Pr_p, double tol)
    {
        A->use_dA = 0;
        int *perm_left = (int *)(*Pl_p);
        int *perm_right = (int *)(*Pr_p);
        int minmn = MIN(m, n_max);
        int j = 0;
        for (int i = 0; i < minmn; i++)
        {
            MatrixInd max_curr = max_el(m, n_max, A, i, i);
            if (abs(MATEL(A, max_curr.ai, max_curr.aj)) < tol)
            {
                break;
            }
            // switch rows
            ROWSW(n, A, i, 0, A, max_curr.ai, 0);
            // save in permutation vector
            perm_left[i] = max_curr.ai;
            // switch cols
            COLSW(m, A, 0, i, A, 0, max_curr.aj);
            // save in permutation vector
            perm_right[i] = max_curr.aj;
            for (int j = i + 1; j < m; j++)
            {
                double Lji = MATEL(A, j, i) / MATEL(A, i, i);
                MATEL(A, j, i) = Lji;
                GEAD(1, n - (i + 1), -Lji, A, i, i + 1, A, j, i + 1);
            }
            j = i + 1;
        }
        rank = j;
    };
    /** \brief Function to calculate LU factorization but A, and result (L and U) are transposed, all indices refer to the dimensions of the original A matrix (and not the transposed one) */
    void LU_FACT_transposed(const int m, const int n, const int n_max, int &rank, MAT *At, PMAT *Pl_p, PMAT *Pr_p, double tol)
    {
        At->use_dA = 0;
        int *perm_left = (int *)(*Pl_p);
        int *perm_right = (int *)(*Pr_p);
        int minmn = MIN(m, n_max);
        int j = 0;
        for (int i = 0; i < minmn; i++)
        {
            MatrixInd max_curr = max_el(n_max, m, At, i, i);
            if (abs(MATEL(At, max_curr.ai, max_curr.aj)) < tol)
            {
                break;
            }
            // switch rows
            COLSW(n, At, 0, i, At, 0, max_curr.aj);
            // save in permutation vector
            perm_left[i] = max_curr.aj;
            // switch cols
            ROWSW(m, At, i, 0, At, max_curr.ai, 0);
            // save in permutation vector
            perm_right[i] = max_curr.ai;
            for (int j = i + 1; j < m; j++)
            {
                double Lji = MATEL(At, i, j) / MATEL(At, i, i);
                MATEL(At, i, j) = Lji;
                GEAD(n - (i + 1), 1, -Lji, At, i + 1, i, At, i + 1, j);
            }
            j = i + 1;
        }
        rank = j;
    };

    void fatrop_dtrsv_unu(const int m, const int n, blasfeo_dmat *sA, const int ai, const int aj, blasfeo_dvec *sx, const int xi, blasfeo_dvec *sz, const int zi)
    {
        for (int i = m; i < n; i++)
        {
            VECEL(sz, zi + i) = VECEL(sx, xi + i);
        }
        for (int i = m - 1; i >= 0; i--)
        {
            double res = VECEL(sx, xi + i);
            for (int j = i + 1; j < n; j++)
            {
                res -= MATEL(sA, ai + i, aj + j) * VECEL(sz, zi + j);
            }
            VECEL(sz, zi + i) = res;
        }
    }
    void fatrop_dtrsv_utu(const int m, blasfeo_dmat *sA, const int ai, const int aj, blasfeo_dvec *sx, const int xi, blasfeo_dvec *sz, const int zi)
    {
        for (int i = 0; i < m; i++)
        {
            double res = VECEL(sx, xi + i);
            for (int j = 0; j < i; j++)
            {
                res -= MATEL(sA, ai + j, aj + i) * VECEL(sz, zi + j);
            }
            VECEL(sz, zi + i) = res;
        }
    }
    void fatrop_identity(const int m, MAT *sA, const int ai, const int aj)
    {
        GESE(m, m, 0.0, sA, ai, aj);
        DIARE(m, 1.0, sA, ai, aj);
    }
    void fatrop_drowad(int kmax, double alpha, struct blasfeo_dvec *sx, int xi, struct blasfeo_dmat *sA, int ai, int aj)
    {
        for (int i = 0; i < kmax; i++)
        {
            MATEL(sA, ai, aj + i) += alpha * VECEL(sx, xi + i);
        }
    }
    void axpy(const double alpha, const FatropVecBF &va, const FatropVecBF &vb, const FatropVecBF &vc)
    {
        DBGASSERT(va.nels() == vb.nels());
        DBGASSERT(va.nels() == vc.nels());
        VEC *va_p = (VEC *)va;
        VEC *vb_p = (VEC *)vb;
        VEC *vc_p = (VEC *)vc;
        AXPY(va.nels(), alpha, va_p, va.offset(), vb_p, vb.offset(), vc_p, vc.offset());
    };
    void copy(const FatropVecBF &va, const FatropVecBF &vb)
    {
        DBGASSERT(va.nels() == vb.nels());
        VEC *va_p = (VEC *)va;
        VEC *vb_p = (VEC *)vb;
        VECCP(va.nels(), va_p, va.offset(), vb_p, vb.offset());
    };
    void axpby(const double alpha, const FatropVecBF &va, const double beta, const FatropVecBF &vb, const FatropVecBF &vc)
    {
        DBGASSERT(va.nels() == vb.nels());
        DBGASSERT(va.nels() == vc.nels());
        VEC *va_p = (VEC *)va;
        VEC *vb_p = (VEC *)vb;
        VEC *vc_p = (VEC *)vc;
        AXPBY(va.nels(), alpha, va_p, va.offset(), beta, vb_p, vb.offset(), vc_p, vc.offset());
    };
    double dot(const FatropVecBF &va, FatropVecBF &vb)
    {
        DBGASSERT(va.nels() == vb.nels());
        VEC *va_p = (VEC *)va;
        VEC *vb_p = (VEC *)vb;
        return DOT(va.nels(), va_p, va.offset(), vb_p, vb.offset());
    };
    double Linf(const FatropVecBF &va)
    {
        VEC *va_p = (VEC *)va;
        int nels = va.nels();
        int offset = va.offset();
        double res = 0.0;
        for (int i = offset; i < nels + offset; i++)
        {
            res = MAX(res, abs(VECEL(va_p, i)));
        }
        return res;
    };
    double LinfScaled(const FatropVecBF &va, const FatropVecBF &scales)
    {
        VEC *va_p = (VEC *)va;
        VEC *scales_p = (VEC *)scales;
        int nels = va.nels();
        int offset = va.offset();
        double res = 0.0;
        for (int i = offset; i < nels + offset; i++)
        {
            res = MAX(res, abs(VECEL(va_p, i)) / (1. + abs(VECEL(scales_p, i))));
        }
        return res;
    };
    double minabs(const FatropVecBF &va)
    {
        VEC *va_p = (VEC *)va;
        int nels = va.nels();
        int offset = va.offset();
        if (nels == 0)
        {
            return 0.0;
        }
        double res = abs(VECEL(va_p, offset));
        for (int i = offset + 1; i < nels + offset; i++)
        {
            res = MIN(res, abs(VECEL(va_p, i)));
        }
        return res;
    };
    double L1(const FatropVecBF &va)
    {
        VEC *va_p = (VEC *)va;
        int nels = va.nels();
        int offset = va.offset();
        double res = 0.0;
        for (int i = offset; i < nels + offset; i++)
        {
            res += abs(VECEL(va_p, i));
        }
        return res;
    };
} // namespace fatrop
using namespace gen_riccati;
FatropMatBF::FatropMatBF(const int nrows, const int ncols, const int row_offset, const int col_offset) : row_offset_(row_offset), col_offset_(col_offset), nrows_(nrows), ncols_(ncols)
{
}
FatropMatBF::FatropMatBF(const int nrows, const int ncols, const int row_offset, const int col_offset, MAT *matbf) : mat_(matbf), row_offset_(row_offset), col_offset_(col_offset), nrows_(nrows), ncols_(ncols)
{
}
FatropMatBF::FatropMatBF(MAT *matbf) : mat_(matbf), row_offset_(0), col_offset_(0), nrows_(matbf->m), ncols_(matbf->n)
{
}
void FatropMatBF::operator=(const FatropMat &fm)
{
    for (int ai = 0; ai < fm.nrows(); ai++)
    // for (int ai = 0; ai < nrows_; ai++)
    {
        for (int aj = 0; aj < fm.ncols(); aj++)
        // for (int aj = 0; aj < ncols_; aj++)
        {
            this->at(ai, aj) = fm.get_el(ai, aj);
        }
    }
}
FatropMemoryMatBF::FatropMemoryMatBF(const NumericVector &nrows, const NumericVector &ncols, int N) : N_(N), nrows_(nrows), ncols_(ncols)
// TODO: if rvalue-reference is used -> unecessary copy, use move sementics instead.
{
    set_up();
}
FatropMemoryMatBF::FatropMemoryMatBF(const int nrows, const int ncols, int N) : N_(N), nrows_(vector<int>(N, nrows)), ncols_(vector<int>(N, ncols))
{
    set_up();
}
int FatropMemoryMatBF::memory_size() const
{
    int result = 0;
    // size to store structs
    result += N_ * sizeof(MAT);
    // sufficient space for cache alignment
    result = (result + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE * CACHE_LINE_SIZE + CACHE_LINE_SIZE;
    // size to store date
    for (int i = 0; i < N_; i++)
    {
        // result += MEMSIZE_MAT(nrows_.at(i), ncols_.at(i));
        result += MEMSIZE_MAT(nrows_[i], ncols_[i]);
    }
    return result;
};
void FatropMemoryMatBF::set_up()
{
    free(mem);
    mem = malloc(this->memory_size());
    char *data_p = (char *)mem;
    MAT *bf_ptr = (MAT *)data_p;
    this->mat = bf_ptr;
    bf_ptr += N_;
    // align with cache line
    long long l_ptr = (long long)bf_ptr;
    l_ptr = (l_ptr + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE * CACHE_LINE_SIZE;
    data_p = (char *)l_ptr;
    double *d_ptr_begin = (double *)data_p;
    for (int i = 0; i < N_; i++)
    {
        CREATE_MAT(nrows_[i], ncols_.at(i), mat + i, data_p);
        data_p += MEMSIZE_MAT(nrows_.at(i), ncols_.at(i));
        // data_p += (mat + i)->memsize;
    }
    double *d_ptr_end = (double *)data_p;
    for (double *d_ptr = d_ptr_begin; d_ptr < d_ptr_end; d_ptr++)
    {
        *d_ptr = 0.0;
    }
    // cout << "allocated memory size " << this->memory_size()<< endl;
    // cout << "used memory size " << (unsigned long long) d_ptr_end - (unsigned long long) mem << endl;
    // cout << "difference " << this->memory_size()-((unsigned long long) d_ptr_end - (unsigned long long) mem)<< endl;
}
FatropMatBF FatropMemoryMatBF::operator[](const int N) const
{
#if DEBUG
    assert(N < N_);
#endif
    MAT *resmat = mat + N;
    FatropMatBF res(resmat->m, resmat->n, 0, 0, resmat);
    return res;
}
FatropMemoryMatBF::~FatropMemoryMatBF()
{
    free(mem);
}
FatropMemoryVecBF::FatropMemoryVecBF(const NumericVector &nels, int N) : N_(N), nels_(nels)
// TODO: if rvalue-reference is used -> unecessary copy, use move sementics instead.
{
    set_up();
}
FatropMemoryVecBF::FatropMemoryVecBF(const int nels, int N) : N_(N), nels_(vector<int>(N, nels))
{
    set_up();
}
int FatropMemoryVecBF::memory_size() const
{
    int result = 0;
    // size to store structs
    result += N_ * sizeof(VEC);
    // sufficient space for cache alignment
    result = (result + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE * CACHE_LINE_SIZE + CACHE_LINE_SIZE;
    // size to store date
    for (int i = 0; i < N_; i++)
    {
        result += MEMSIZE_VEC(nels_.at(i));
    }
    return result;
};
void FatropMemoryVecBF::set_up()
{
    free(mem);
    mem = malloc(this->memory_size());
    char *data_p = (char *)mem;
    VEC *bf_ptr = (VEC *)data_p;
    this->vec = bf_ptr;
    bf_ptr += N_;
    // align with cache line
    long long l_ptr = (long long)bf_ptr;
    l_ptr = (l_ptr + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE * CACHE_LINE_SIZE;
    data_p = (char *)l_ptr;
    double *d_ptr_begin = (double *)data_p;
    for (int i = 0; i < N_; i++)
    {
        CREATE_VEC(nels_.at(i), vec + i, data_p);
        data_p += MEMSIZE_VEC(nels_.at(i));
    }
    double *d_ptr_end = (double *)data_p;
    for (double *d_ptr = d_ptr_begin; d_ptr < d_ptr_end; d_ptr++)
    {
        *d_ptr = 0.0;
    }
}

FatropVecBF::FatropVecBF(const int nels, const int offset) : offset_(offset), nels_(nels)
{
}
FatropVecBF::FatropVecBF(const int nels, const int offset, VEC *vecbf) : vec_(vecbf), offset_(offset), nels_(nels)
{
}
FatropVecBF::operator VEC *() const
{
    return this->vec_;
}
double &FatropVecBF::at(const int ai) const
{
#if DEBUG
    assert(ai < nels_);
#endif
    return VECEL(vec_, ai + offset_);
};
double FatropVecBF::get_el(const int ai) const
{
    return this->at(ai);
}
int FatropVecBF::nels() const
{
    return nels_;
}
int FatropVecBF::offset() const
{
    return offset_;
};
void FatropVecBF::operator=(const FatropVec &fm)
{
    for (int ai = 0; ai < nels_; ai++)
    {
        this->at(ai) = fm.get_el(ai);
    }
}
void FatropVecBF::copy(const FatropVecBF &fm)
{
    DBGASSERT(fm.nels() == nels());
    VEC *fm_p = (VEC *)fm;
    VECCP(nels(), fm_p, fm.offset(), vec_, offset());
}
void FatropVecBF::copyto(vector<double> &fm) const
{
    fm.resize(nels_);
    for (int ai = 0; ai < nels_; ai++)
    {
        fm.at(ai) = this->at(ai);
    }
}
void FatropVecBF::operator=(const vector<double> &fm)
{
    for (int ai = 0; ai < nels_; ai++)
    {
        this->at(ai) = fm.at(ai);
    }
}
void FatropVecBF::set_datap(VEC *vecbf)
{
    vec_ = vecbf;
}
FatropVecBF FatropVecBF::block(const int i, const int p) const
{
    DBGASSERT(i + p <= nels_);
    return FatropVecBF(p, offset_ + i, this->vec_);
}
void FatropVecBF::SwapWith(FatropVecBF &vb)
{
    DBGASSERT(vb.offset_ == offset_);
    DBGASSERT(vb.nels_ == nels_);
    VEC *tmp = vec_;
    vec_ = vb.vec_;
    vb.vec_ = tmp;
}
void FatropVecBF::SetConstant(double constant) const
{
    VECSE(nels_, constant, vec_, offset_);
}
FatropVecBF FatropMemoryVecBF::operator[](const int N) const
{
#if DEBUG
    assert(N < N_);
#endif
    VEC *resvec = vec + N;
    FatropVecBF res(resvec->m, 0, resvec);
    return res;
}
FatropMemoryVecBF::~FatropMemoryVecBF()
{
    free(mem);
}

PermMat::PermMat(const int dim) : dim_(dim)
{
}
PermMat::PermMat(const int dim, int *data) : dim_(dim), data_(data)
{
}
double PermMat::get_el(const int ai, const int aj) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    int aj_one = data_[ai];
    int row_curr = ai - 1;
    while (row_curr >= 0)
    {
        if (aj_one == data_[row_curr])
        {
            aj_one = row_curr;
        }
        row_curr--;
    }
    if (aj == aj_one)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
};
void PermMat::print(const int kmax) const
{
    for (int k = 0; k < kmax; k++)
    {
        cout << k << " <-> " << data_[k] << endl;
    }
}
void PermMat::set_datap(int *data)
{
    data_ = data;
}
void PermMat::set_datap(const int i, const int val)
{
#if DEBUG
    assert(data_ != NULL);
    assert(i < dim_);
#endif
    data_[i] = val;
}
void PermMat::PM(const int kmax, MAT *M) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    ROWPE(kmax, data_, M);
}
void PermMat::PV(const int kmax, VEC *V, const int offs) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    VECPE(kmax, data_, V, offs);
}
void PermMat::PtV(const int kmax, VEC *V, const int offs) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    VECPEI(kmax, data_, V, offs);
}
void PermMat::PM(const int kmax, const int n, MAT *M, const int ai, const int aj) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    // invalidate stored inverse diagonal
    M->use_dA = 0;

    int ii;
    for (ii = 0; ii < kmax; ii++)
    {
        if (data_[ii] != ii)
            ROWSW(n, M, ai + ii, aj, M, ai + data_[ii], aj);
    }
    return;
}
void PermMat::PtM(const int kmax, MAT *M) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    ROWPEI(kmax, data_, M);
}
void PermMat::MP(const int kmax, MAT *M) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    COLPEI(kmax, data_, M);
}
void PermMat::MPt(const int kmax, MAT *M) const
{
#if DEBUG
    assert(data_ != NULL);
#endif
    COLPE(kmax, data_, M);
}

MemoryPermMat::MemoryPermMat(const int dim, const int N) : PermMat(dim), dim_(dim), N_(N)
{
    set_up();
};
int MemoryPermMat::memory_size() const
{
    int size = 0;
    size += N_ * sizeof(PermMat) + N_ * dim_ * sizeof(int);
    return size;
}
void MemoryPermMat::set_up()
{
    free(mem);
    mem = malloc(this->memory_size());
    char *char_p = (char *)mem;
    perm_p = (PermMat *)char_p;
    PermMat *perm_pp = perm_p;
    for (int i = 0; i < N_; i++)
    {
        new (perm_pp) PermMat(dim_);
        perm_pp++;
    }
    int *data_p = (int *)perm_pp;
    this->set_datap(data_p);
    for (int i = 0; i < N_; i++)
    {
        perm_p[i].set_datap(data_p);
        data_p += dim_;
    }
    char_p = (char *)data_p;
}
MemoryPermMat::~MemoryPermMat()
{
    free(mem);
}