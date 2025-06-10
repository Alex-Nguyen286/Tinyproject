#ifndef GENERALLINEARSYSTEM_HPP
#define GENERALLINEARSYSTEM_HPP

#include "Matrix.hpp"
#include "Vector.hpp"

class GeneralLinearSystem {
private:
    const Matrix* mpA;
    const Vector* mpb;

public:
    GeneralLinearSystem(const Matrix& A, const Vector& b);
    Vector Solve() const;
};

#endif