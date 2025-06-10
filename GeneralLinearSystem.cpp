#include "GeneralLinearSystem.hpp"

GeneralLinearSystem::GeneralLinearSystem(const Matrix& A, const Vector& b) : mpA(&A), mpb(&b) {
    assert(A.getNumRows() == b.getSize());
}

Vector GeneralLinearSystem::Solve() const {
    Matrix A_pinv = mpA->pseudoInverse();
    return A_pinv * (*mpb);
}