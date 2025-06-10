#include "PosSymLinSystem.hpp"
#include <cmath>

PosSymLinSystem::PosSymLinSystem(const Matrix& A, const Vector& b) : LinearSystem(A, b) {
    assert(A.isSymmetric());
}

Vector PosSymLinSystem::Solve() const {
    Vector x(mSize);
    Vector r = *mpb - (*mpA * x);
    Vector p = r;
    double rsold = r.dot(r);
    for (int iter = 0; iter < 1000; iter++) {
        Vector Ap = *mpA * p;
        double alpha = rsold / p.dot(Ap);
        x = x + p * alpha;
        r = r - Ap * alpha;
        double rsnew = r.dot(r);
        if (std::sqrt(rsnew) < 1e-6) break;
        p = r + p * (rsnew / rsold);
        rsold = rsnew;
    }
    return x;
}