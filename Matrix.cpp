#include "Matrix.hpp"
#include <cmath>

Matrix::Matrix(int numRows, int numCols) : mNumRows(numRows), mNumCols(numCols) {
    assert(numRows > 0 && numCols > 0);
    mData = new double*[numRows];
    for (int i = 0; i < numRows; i++) {
        mData[i] = new double[numCols];
        for (int j = 0; j < numCols; j++) {
            mData[i][j] = 0.0;
        }
    }
}

Matrix::Matrix(const Matrix& other) : mNumRows(other.mNumRows), mNumCols(other.mNumCols) {
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; i++) {
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; j++) {
            mData[i][j] = other.mData[i][j];
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < mNumRows; i++) {
        delete[] mData[i];
    }
    delete[] mData;
}

Matrix& Matrix::operator=(const Matrix& other) {
    if (this != &other) {
        for (int i = 0; i < mNumRows; i++) {
            delete[] mData[i];
        }
        delete[] mData;
        mNumRows = other.mNumRows;
        mNumCols = other.mNumCols;
        mData = new double*[mNumRows];
        for (int i = 0; i < mNumRows; i++) {
            mData[i] = new double[mNumCols];
            for (int j = 0; j < mNumCols; j++) {
                mData[i][j] = other.mData[i][j];
            }
        }
    }
    return *this;
}

double& Matrix::operator()(int i, int j) {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i-1][j-1];
}

const double& Matrix::operator()(int i, int j) const {
    assert(i >= 1 && i <= mNumRows && j >= 1 && j <= mNumCols);
    return mData[i-1][j-1];
}

Matrix Matrix::operator+(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] + other.mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] - other.mData[i][j];
        }
    }
    return result;
}

Matrix Matrix::operator*(double scalar) const {
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[i][j] = mData[i][j] * scalar;
        }
    }
    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    assert(mNumCols == other.mNumRows);
    Matrix result(mNumRows, other.mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < other.mNumCols; j++) {
            for (int k = 0; k < mNumCols; k++) {
                result.mData[i][j] += mData[i][k] * other.mData[k][j];
            }
        }
    }
    return result;
}

Vector Matrix::operator*(const Vector& v) const {
    assert(mNumCols == v.getSize());
    Vector result(mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result[i] += mData[i][j] * v[j];
        }
    }
    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            result.mData[j][i] = mData[i][j];
        }
    }
    return result;
}

double Matrix::determinant() const {
    assert(mNumRows == mNumCols);
    if (mNumRows == 1) return mData[0][0];
    if (mNumRows == 2) return mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];
    double det = 0.0;
    for (int j = 0; j < mNumCols; j++) {
        Matrix minor(mNumRows - 1, mNumCols - 1);
        for (int i = 1; i < mNumRows; i++) {
            for (int k = 0; k < mNumCols; k++) {
                if (k < j) minor.mData[i-1][k] = mData[i][k];
                else if (k > j) minor.mData[i-1][k-1] = mData[i][k];
            }
        }
        det += (j % 2 == 0 ? 1 : -1) * mData[0][j] * minor.determinant();
    }
    return det;
}

Matrix Matrix::inverse() const {
    assert(mNumRows == mNumCols);
    int n = mNumRows;
    Matrix augmented(n, 2*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented(i+1, j+1) = mData[i][j];
        }
        augmented(i+1, n+i+1) = 1.0;
    }
    for (int k = 0; k < n; k++) {
        int pivot = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(augmented(i+1, k+1)) > std::abs(augmented(pivot+1, k+1))) {
                pivot = i;
            }
        }
        if (pivot != k) {
            for (int j = 0; j < 2*n; j++) {
                std::swap(augmented(k+1, j+1), augmented(pivot+1, j+1));
            }
        }
        double pivot_val = augmented(k+1, k+1);
        assert(std::abs(pivot_val) > 1e-10);
        for (int i = 0; i < n; i++) {
            if (i != k) {
                double factor = augmented(i+1, k+1) / pivot_val;
                for (int j = k + 1; j < 2*n; j++) {
                    augmented(i+1, j+1) -= factor * augmented(k+1, j+1);
                }
            }
        }
        for (int j = k + 1; j < 2*n; j++) {
            augmented(k+1, j+1) /= pivot_val;
        }
    }
    Matrix inv(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv(i+1, j+1) = augmented(i+1, n+j+1);
        }
    }
    return inv;
}

Matrix Matrix::pseudoInverse() const {
    if (mNumRows >= mNumCols) {
        Matrix At = transpose();
        Matrix AtA = At * (*this);
        Matrix AtA_inv = AtA.inverse();
        return AtA_inv * At;
    } else {
        Matrix At = transpose();
        Matrix AAt = (*this) * At;
        Matrix AAt_inv = AAt.inverse();
        return At * AAt_inv;
    }
}

bool Matrix::isSymmetric() const {
    if (mNumRows != mNumCols) return false;
    for (int i = 0; i < mNumRows; i++) {
        for (int j = i + 1; j < mNumCols; j++) {
            if (mData[i][j] != mData[j][i]) return false;
        }
    }
    return true;
}