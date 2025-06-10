#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"
#include <cassert>

class Matrix {
private:
    int mNumRows;
    int mNumCols;
    double** mData;

public:
    Matrix(int numRows, int numCols);
    Matrix(const Matrix& other);
    ~Matrix();
    Matrix& operator=(const Matrix& other);
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;
    int getNumRows() const { return mNumRows; }
    int getNumCols() const { return mNumCols; }
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    Matrix operator*(const Matrix& other) const;
    Vector operator*(const Vector& v) const;
    Matrix transpose() const;
    double determinant() const;
    Matrix inverse() const;
    Matrix pseudoInverse() const;
    bool isSymmetric() const;
};

#endif