#pragma once
#include <vector>
namespace gen_riccati
{
    class NumericVector : public std::vector<int>
    {
    public:
        using std::vector<int>::vector;
        NumericVector(const std::vector<int> &v) : std::vector<int>(v){};
        friend NumericVector operator+(const NumericVector &a, const NumericVector &b)
        {
            NumericVector res;
            for (size_t i = 0; i < a.size(); i++)
            {
                res.push_back(a[i] + b[i]);
            }
            return res;
        };

        friend NumericVector operator+(const NumericVector &a, const int b)
        {
            NumericVector res;
            for (size_t i = 0; i < a.size(); i++)
            {
                res.push_back(a[i] + b);
            }
            return res;
        };

        friend  NumericVector rotate(const NumericVector &a, const int shift)
        {
            NumericVector res;
            for (size_t i = 0; i < a.size(); i++)
            {
                res.push_back(a[(i + shift) % a.size()]);
            }
            return res;
        };
        friend int sum(const NumericVector &a)
        {
            int res = 0;
            for (size_t i = 0; i < a.size(); i++)
            {
                res += a[i];
            }
            return res;
        };
        friend int max(const NumericVector &a)
        {
            int res = 0;
            for (size_t i = 0; i < a.size(); i++)
            {
                res = std::max(res, a[i]);
            }
            return res;
        };
        friend NumericVector offsets(const NumericVector &a)
        {
            NumericVector res;
            res.reserve(a.size());
            int offset = 0;
            for (size_t i = 0; i < a.size(); i++)
            {
                res.push_back(offset);
                offset += a[i];
            }
            return res;
        };
    };

} // namespace symspals