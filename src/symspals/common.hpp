#pragma once
#include <unordered_set>

struct Index
{
    int row = 0;
    int col = 0;
    bool operator==(const Index &other) const
    {
        return row == other.row && col == other.col;
    }
};
template <typename T>
struct Triplet
{
    Index index;
    T value;
    Triplet(int row, int col, T value) : index{row, col}, value(value){};
};

namespace std
{
    template <>
    struct hash<Index>
    {
        std::size_t operator()(const Index &k) const
        {
            using std::hash;
            using std::size_t;
            return (hash<int>()(k.row) ^ (hash<int>()(k.col) << 1));
        }
    };
};
