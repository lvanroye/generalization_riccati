#pragma once
struct Index
{
    int row = 0;
    int col = 0;
};
template <typename T>
struct Triplet
{
    Index index;
    T value;
    Triplet(int row, int col, T value) : index{row, col}, value(value){};
};