#include "Regression.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

// Helper to parse CSV and create matrix A and vector b
Regression::Regression(const std::string& csvFile) {
    std::ifstream file(csvFile);
    std::string line;
    std::vector<std::vector<double>> features;
    std::vector<double> targets;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> row;
        int colIdx = 0;
        while (std::getline(ss, cell, ',')) {
            if (colIdx >= 2 && colIdx <= 7) row.push_back(std::stod(cell)); // features
            if (colIdx == 8) targets.push_back(std::stod(cell)); // PRP
            ++colIdx;
        }
        if (row.size() == 6) features.push_back(row);
    }

    int n = features.size();
    A = Matrix(n, 6);
    b = Vector(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 6; ++j)
            A(i + 1, j + 1) = features[i][j];
        b(i + 1) = targets[i];
    }
}

void Regression::Train() {
    Matrix At = Matrix(6, A.Rows());
    for (int i = 1; i <= A.Rows(); ++i)
        for (int j = 1; j <= A.Cols(); ++j)
            At(j, i) = A(i, j);

    Matrix AtA = At * A;
    Vector Atb = At * b;
    LinearSystem system(AtA, Atb);
    coefficients = system.Solve();
}

double Regression::Evaluate(const std::string& testCsvFile) const {
    std::ifstream file(testCsvFile);
    std::string line;
    std::vector<std::vector<double>> testX;
    std::vector<double> testY;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;
        std::vector<double> row;
        int colIdx = 0;
        while (std::getline(ss, cell, ',')) {
            if (colIdx >= 2 && colIdx <= 7) row.push_back(std::stod(cell));
            if (colIdx == 8) testY.push_back(std::stod(cell));
            ++colIdx;
        }
        if (row.size() == 6) testX.push_back(row);
    }

    double mse = 0.0;
    for (size_t i = 0; i < testX.size(); ++i) {
        double pred = 0.0;
        for (int j = 0; j < 6; ++j) pred += testX[i][j] * coefficients[j];
        double err = pred - testY[i];
        mse += err * err;
    }
    return std::sqrt(mse / testX.size());
}

void Regression::PrintCoefficients() const {
    std::cout << "Model Coefficients:\n";
    coefficients.Print();
}
