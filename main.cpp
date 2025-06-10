#include "Vector.hpp"
#include "Matrix.hpp"
#include "GeneralLinearSystem.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

// Function to standardize features (returns false if standardization fails)
bool standardize(Matrix& A, Vector& mean, Vector& std, int rows, int cols) {
    for (int j = 1; j <= cols; j++) {
        double sum = 0.0;
        for (int i = 1; i <= rows; i++) {
            sum += A(i, j);
        }
        mean(j) = sum / rows;
        double sum_sq = 0.0;
        for (int i = 1; i <= rows; i++) {
            sum_sq += (A(i, j) - mean(j)) * (A(i, j) - mean(j));
        }
        std(j) = std::sqrt(sum_sq / (rows - 1)); // Sample standard deviation
        if (std(j) < 1e-10) {
            std::cerr << "Warning: Standard deviation for feature " << j << " is near zero.\n";
            std(j) = 1.0; // Prevent division by zero
        }
        for (int i = 1; i <= rows; i++) {
            A(i, j) = (A(i, j) - mean(j)) / std(j);
        }
    }
    return true;
}

// Function to apply standardization to test set
void apply_standardization(Matrix& A, const Vector& mean, const Vector& std, int rows, int cols) {
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            A(i, j) = (A(i, j) - mean(j)) / std(j);
        }
    }
}

// Function to compute RMSE
double compute_rmse(const Matrix& A, const Vector& b, const Vector& x, int rows, double bias = 0.0) {
    Vector predicted = A * x;
    double sum_sq_error = 0.0;
    for (int i = 1; i <= rows; i++) {
        double error = b(i) - (predicted(i) + bias);
        sum_sq_error += error * error;
    }
    return std::sqrt(sum_sq_error / rows);
}

// Function to solve Ridge regression with bias term
std::pair<Vector, double> solve_ridge(const Matrix& A, const Vector& b, double lambda, int rows, int cols) {
    // Add bias term by including a column of ones in A
    Matrix A_with_bias(rows, cols + 1);
    for (int i = 1; i <= rows; i++) {
        for (int j = 1; j <= cols; j++) {
            A_with_bias(i, j) = A(i, j);
        }
        A_with_bias(i, cols + 1) = 1.0; // Bias column
    }

    Matrix At = A_with_bias.transpose();
    Matrix AtA = At * A_with_bias;
    Matrix I(cols + 1, cols + 1);
    for (int i = 1; i <= cols + 1; i++) I(i, i) = 1.0; // Identity matrix
    Matrix AtA_reg = AtA + I * lambda; // Add regularization term
    Matrix AtA_reg_inv = AtA_reg.inverse();
    Vector Atb = At * b;
    Vector x = AtA_reg_inv * Atb;

    // Extract bias (last element of x)
    double bias = x(cols + 1);
    Vector coefficients(cols);
    for (int i = 1; i <= cols; i++) {
        coefficients(i) = x(i);
    }
    return {coefficients, bias};
}

// Function to detect outliers (simple threshold-based)
bool is_outlier(const std::vector<double>& row, const Vector& mean, const Vector& std) {
    for (size_t j = 0; j < 6; j++) {
        double z_score = std::abs((row[j] - mean(j + 1)) / std(j + 1));
        if (z_score > 3.0) return true; // Z-score > 3 indicates an outlier
    }
    return false;
}

int main() {
    // Load dataset
    std::vector<std::vector<double>> data;
    std::ifstream file("machine.data");
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open machine.data\n";
        return 1;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<double> row;
        int count = 0;
        while (std::getline(ss, token, ',')) {
            if (count >= 2 && count <= 8) {
                try {
                    row.push_back(std::stod(token));
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing line: " << line << "\n";
                    file.close();
                    return 1;
                }
            }
            count++;
        }
        if (row.size() == 7) {
            data.push_back(row);
        }
    }
    file.close();
    if (data.size() != 209) {
        std::cerr << "Error: Expected 209 instances, got " << data.size() << "\n";
        return 1;
    }

    // Compute initial mean and std for outlier detection
    Vector mean(6), std(6);
    Matrix A_temp(data.size(), 6);
    for (size_t i = 0; i < data.size(); i++) {
        for (size_t j = 0; j < 6; j++) {
            A_temp(i + 1, j + 1) = data[i][j];
        }
    }
    standardize(A_temp, mean, std, data.size(), 6);

    // Remove outliers
    std::vector<std::vector<double>> clean_data;
    for (const auto& row : data) {
        if (!is_outlier(row, mean, std)) {
            clean_data.push_back(row);
        }
    }
    std::cout << "Removed " << (209 - clean_data.size()) << " outliers. New dataset size: " << clean_data.size() << "\n";

    // Shuffle indices for train-test split
    std::vector<int> indices(clean_data.size());
    for (size_t i = 0; i < clean_data.size(); i++) indices[i] = i;
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine{});

    int train_size = static_cast<int>(clean_data.size() * 0.8);
    int test_size = clean_data.size() - train_size;
    std::vector<int> training_indices(indices.begin(), indices.begin() + train_size);
    std::vector<int> testing_indices(indices.begin() + train_size, indices.end());

    // Populate training data
    Matrix A_train(train_size, 6);
    Vector b_train(train_size);
    for (int i = 0; i < train_size; i++) {
        int idx = training_indices[i];
        for (int j = 0; j < 6; j++) {
            A_train(i + 1, j + 1) = clean_data[idx][j];
        }
        b_train(i + 1) = clean_data[idx][6];
    }

    // Populate testing data
    Matrix A_test(test_size, 6);
    Vector b_test(test_size);
    for (int i = 0; i < test_size; i++) {
        int idx = testing_indices[i];
        for (int j = 0; j < 6; j++) {
            A_test(i + 1, j + 1) = clean_data[idx][j];
        }
        b_test(i + 1) = clean_data[idx][6];
    }

    // Standardize features
    standardize(A_train, mean, std, train_size, 6);
    apply_standardization(A_test, mean, std, test_size, 6);

    // Perform 10-fold cross-validation to select lambda
    double best_lambda = 0.0;
    double best_cv_rmse = std::numeric_limits<double>::max();
    std::vector<double> lambdas = {4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0};
    int k = 10; // k-fold cross-validation
    int fold_size = train_size / k;

    for (double lambda : lambdas) {
        double cv_rmse = 0.0;
        for (int fold = 0; fold < k; fold++) {
            Matrix A_cv_train(train_size - fold_size, 6);
            Vector b_cv_train(train_size - fold_size);
            Matrix A_cv_val(fold_size, 6);
            Vector b_cv_val(fold_size);

            int cv_train_idx = 0, cv_val_idx = 0;
            for (int i = 0; i < train_size; i++) {
                if (i >= fold * fold_size && i < (fold + 1) * fold_size) {
                    for (int j = 1; j <= 6; j++) {
                        A_cv_val(cv_val_idx + 1, j) = A_train(i + 1, j);
                    }
                    b_cv_val(cv_val_idx + 1) = b_train(i + 1);
                    cv_val_idx++;
                } else {
                    for (int j = 1; j <= 6; j++) {
                        A_cv_train(cv_train_idx + 1, j) = A_train(i + 1, j);
                    }
                    b_cv_train(cv_train_idx + 1) = b_train(i + 1);
                    cv_train_idx++;
                }
            }

            auto [x_cv, bias_cv] = solve_ridge(A_cv_train, b_cv_train, lambda, train_size - fold_size, 6);
            cv_rmse += compute_rmse(A_cv_val, b_cv_val, x_cv, fold_size, bias_cv);
        }
        cv_rmse /= k;
        if (cv_rmse < best_cv_rmse) {
            best_cv_rmse = cv_rmse;
            best_lambda = lambda;
        }
    }

    // Train final model with best lambda
    auto [x, bias] = solve_ridge(A_train, b_train, best_lambda, train_size, 6);

    // Output model parameters
    std::cout << "Model parameters (lambda = " << best_lambda << ", bias = " << bias << "):\n";
    for (int i = 1; i <= 6; i++) {
        std::cout << "x" << i << " = " << x(i) << "\n";
    }

    // Compute RMSE on test set
    double rmse = compute_rmse(A_test, b_test, x, test_size, bias);
    std::cout << "RMSE on testing set: " << rmse << "\n";

    return 0;
}