#include "LinearSystem.hpp"
#include <algorithm>

LinearSystem::LinearSystem(const Matrix& A, const Vector& b) : mpA(&A), mpb(&b) {
    assert(A.getNumRows() == A.getNumCols());
    assert(A.getNumRows() == b.getSize());
    mSize = A.getNumRows();
}

Vector LinearSystem::Solve() const {
    Matrix A = *mpA;
    Vector b = *mpb;
    for (int k = 0; k < mSize - 1; k++) {
        int pivot = k;
        for (int i = k + 1; i < mSize; i++) {
            if (std::abs(A(i+1, k+1)) > std::abs(A(pivot+1, k+1))) {
                pivot = i;
            }
        }
        if (pivot != k) {
            for (int j = 0; j < mSize; j++) {
                std::swap(A(k+1, j+1), A(pivot+1, j+1));
            }
            std::swap(b(k+1), b(pivot+1));
        }
        for (int i = k + 1; i < mSize; i++) {
            double factor = A(i+1, k+1) / A(k+1, k+1);
            for (int j = k + 1; j < mSize; j++) {
                A(i+1, j+1) -= factor * A(k+1, j+1);
            }
            b(i+1) -= factor * b(k+1);
        }
    }
    Vector x(mSize);
    for (int i = mSize - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < mSize; j++) {
            sum += A(i+1, j+1) * x(j+1);
        }
        x(i+1) = (b(i+1) - sum) / A(i+1, i+1);
    }
    return x;
}