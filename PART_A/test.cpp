// This is a simple C++ program that tests matrix and vector operations,
// linear systems, positive definite symmetric linear systems, and Tikhonov regularization.
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Matrix.h"
#include "Vector.h"
#include "LinearSystem.h"
#include "PosSymLinSystem.h"
#include "Tikhonov.h"

void PrintMatrix(const Matrix& M) {
    for (int i = 1; i <= M.Rows(); ++i) {
        for (int j = 1; j <= M.Cols(); ++j) {
            std::cout << M(i, j) << " ";
        }
        std::cout << "\n";
    }
}

void PrintVector(const Vector& v) {
    for (int i = 1; i <= v.Size(); ++i) {
        std::cout << v(i) << " ";
    }
    std::cout << "\n";
}

int main() {
    std::srand(static_cast<unsigned>(std::time(0)));

    const int size = 3;

    Matrix A(size, size);
    Matrix B(size, size);
    Vector b(size);

    // Fill B and b with random values (1 to 10)
    for (int i = 1; i <= size; ++i) {
        b(i) = rand() % 10 + 1;
        for (int j = 1; j <= size; ++j) {
            B(i, j) = rand() % 10 + 1;
        }
    }

    // Fill A with random values and make it diagonally dominant
    for (int i = 1; i <= size; ++i) {
        double row_sum = 0;
        for (int j = 1; j <= size; ++j) {
            if (i != j) {
                A(i, j) = rand() % 10 + 1;
                row_sum += std::abs(A(i, j));
            }
        }
        A(i, i) = row_sum + (rand() % 10 + 1);
    }

    std::cout << "Matrix A:\n"; PrintMatrix(A);
    std::cout << "Matrix B:\n"; PrintMatrix(B);
    std::cout << "Vector b:\n"; PrintVector(b);

    std::cout << "\nA + B:\n"; PrintMatrix(A + B);
    std::cout << "A - B:\n"; PrintMatrix(A - B);
    std::cout << "A * B:\n"; PrintMatrix(A * B);
    std::cout << "Transpose of A:\n"; PrintMatrix(A.Transpose());
    std::cout << "A * 2:\n"; PrintMatrix(A * 2.0);

    if (A.Determinant() != 0) {
        std::cout << "Determinant of A: " << A.Determinant() << "\n";
        std::cout << "Inverse of A:\n"; PrintMatrix(A.Inverse());
    } else {
        std::cout << "Matrix A is singular, skipping inverse.\n";
    }

    // Test LinearSystem
    try {
        LinearSystem system(A, b);
        Vector x = system.Solve();
        std::cout << "\nSolution of Ax = b using LinearSystem:\n";
        PrintVector(x);
    } catch (...) {
        std::cout << "LinearSystem failed.\n";
    }

    // Test PosSymLinSystem
    try {
        std::cout << "\n=== Testing PosSymLinSystem ===\n";
        PosSymLinSystem psystem(A, b);
        Vector x_psd = psystem.Solve();
        std::cout << "Solution using PosSymLinSystem:\n";
        PrintVector(x_psd);
    } catch (const std::exception& e) {
        std::cout << "PosSymLinSystem failed: " << e.what() << "\n";
    } catch (...) {
        std::cout << "PosSymLinSystem failed with unknown error.\n";
    }


    // Test TikhonovSolver
    try {
        double lambda = 0.1;
        TikhonovSolver tikSolver(lambda);
        Vector x_tikh = tikSolver.Solve(A, b);
        std::cout << "\nSolution using Tikhonov Regularization (lambda = " << lambda << "):\n";
        PrintVector(x_tikh);
    } catch (...) {
        std::cout << "TikhonovSolver failed.\n";
    }

    return 0;
}