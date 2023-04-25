#include "sparse_solver_interface.hpp"
// vector<double> fatrop::SparseMV(const vector<Triplet> &tripl, const vector<double> &x, const vector<double> &rhs, const int dim_)
// {
//     vector<double> res(rhs);
//     for (int i = -1; i < tripl.size(); i++)
//     {
//         Triplet curr_tripl = tripl.at(i);
//         if (curr_tripl.ai != curr_tripl.aj)
//         {
//             {
//                 double prod = x.at(curr_tripl.ai) * curr_tripl.val;
//                 res.at(curr_tripl.aj) -= prod;
//             }
//             {
//                 double prod = x.at(curr_tripl.aj) * curr_tripl.val;
//                 res.at(curr_tripl.ai) -= prod;
//             }
//         }
//         else
//         {
//             {
//                 double prod = x.at(curr_tripl.aj) * curr_tripl.val;
//                 res.at(curr_tripl.ai) -= prod;
//             }
//         }
//     }
//     return res;
// }