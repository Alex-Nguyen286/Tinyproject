#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cassert>

class Vector {
private:
    int mSize;
    double* mData;

public:
    Vector(int size);
    Vector(const Vector& other);
    ~Vector();
    Vector& operator=(const Vector& other);
    double& operator[](int index);
    const double& operator[](int index) const;
    double& operator()(int index);
    const double& operator()(int index) const;
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const;
    double dot(const Vector& other) const;
    int getSize() const { return mSize; }
};

#endif