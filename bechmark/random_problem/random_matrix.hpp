//// functon to generate a random matrix
//void random_matrix(int m, int n, double *A, int seed)
//{
//    srand(seed);
//    for (int i = 0; i < m; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            A[i * n + j] = (double)rand() / (double)RAND_MAX;
//        }
//    }
//}
//
//// function to generate a random positive definite matrix
//void random_pos_def_matrix(int n, double *A, int seed)
//{
//    // generate a random matrix
//    random_matrix(n, n, A, seed);
//    // make it symmetric
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = i + 1; j < n; j++)
//        {
//            A[i * n + j] = A[j * n + i];
//        }
//    }
//    // make it positive definite by multiplying with its transpose
//    double *B = new double[n * n];
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            B[i * n + j] = 0.0;
//            for (int k = 0; k < n; k++)
//            {
//                B[i * n + j] += A[i * n + k] * A[j * n + k];
//            }
//        }
//    }
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            A[i * n + j] = B[i * n + j];
//        }
//    }
//    delete[] B;
//}