#ifndef POSSYMLINSYSTEM_HPP
#define POSSYMLINSYSTEM_HPP

#include "LinearSystem.hpp"

class PosSymLinSystem : public LinearSystem {
public:
    PosSymLinSystem(const Matrix& A, const Vector& b);
    Vector Solve() const override;
};

#endif