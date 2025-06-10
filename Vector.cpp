#include "Vector.hpp"

Vector::Vector(int size) : mSize(size) {
    assert(size > 0);
    mData = new double[size];
    for (int i = 0; i < size; i++) {
        mData[i] = 0.0;
    }
}

Vector::Vector(const Vector& other) : mSize(other.mSize) {
    mData = new double[mSize];
    for (int i = 0; i < mSize; i++) {
        mData[i] = other.mData[i];
    }
}

Vector::~Vector() {
    delete[] mData;
}

Vector& Vector::operator=(const Vector& other) {
    if (this != &other) {
        delete[] mData;
        mSize = other.mSize;
        mData = new double[mSize];
        for (int i = 0; i < mSize; i++) {
            mData[i] = other.mData[i];
        }
    }
    return *this;
}

double& Vector::operator[](int index) {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

const double& Vector::operator[](int index) const {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

double& Vector::operator()(int index) {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

const double& Vector::operator()(int index) const {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

Vector Vector::operator+(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] + other.mData[i];
    }
    return result;
}

Vector Vector::operator-(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] - other.mData[i];
    }
    return result;
}

Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; i++) {
        result.mData[i] = mData[i] * scalar;
    }
    return result;
}

double Vector::dot(const Vector& other) const {
    assert(mSize == other.mSize);
    double sum = 0.0;
    for (int i = 0; i < mSize; i++) {
        sum += mData[i] * other.mData[i];
    }
    return sum;
}