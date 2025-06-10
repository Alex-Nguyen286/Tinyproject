#ifndef REGRESSION_H
#define REGRESSION_H

#include "Matrix.hpp"
#include "Vector.h"
#include <vector>
#include <string>

class Regression {
private:
    Matrix A;
    Vector b;
    Vector coefficients;

public:
    Regression(const std::string& csvFile);
    void Train();
    double Evaluate(const std::string& testCsvFile) const;
    void PrintCoefficients() const;
};

#endif
