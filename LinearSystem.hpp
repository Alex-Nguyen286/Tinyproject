#ifndef LINEARSYSTEM_HPP
#define LINEARSYSTEM_HPP

#include "Matrix.hpp"
#include "Vector.hpp"

class LinearSystem {
protected:
 int mSize;
 const Matrix* mpA;
 const Vector* mpb;

public:
 LinearSystem(const Matrix& A, const Vector& b);
 virtual Vector Solve() const;
};

#endif