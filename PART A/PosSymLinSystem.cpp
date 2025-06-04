// PosSymLinSystem.cpp
#include "PosSymLinSystem.h"
#include <cmath>
#include <stdexcept>
#include <iostream> // for optional debug output

PosSymLinSystem::PosSymLinSystem(const Matrix& A, const Vector& b)
    : LinearSystem(A, b) {

    // Check that matrix is square
    if (A.Rows() != A.Cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }

    // Check that matrix and vector sizes match
    if (A.Rows() != b.Size()) {
        throw std::invalid_argument("Matrix row count must match vector size.");
    }

    // Check symmetry: A(i, j) == A(j, i)
    const double epsilon = 1e-8;
    for (int i = 1; i <= mSize; ++i) {
        for (int j = i + 1; j <= mSize; ++j) {
            if (std::abs((*mpA)(i, j) - (*mpA)(j, i)) > epsilon) {
                throw std::invalid_argument("Matrix is not symmetric.");
            }
        }
    }

    if (mSize == 0) {
        throw std::invalid_argument("Matrix size cannot be zero.");
    }
}

PosSymLinSystem::~PosSymLinSystem() {}

Vector PosSymLinSystem::Solve() const {
    const int n = mSize;
    const double tol = 1e-10;
    const int max_iter = 1000;

    Vector x(n);                   // Initial guess: zero vector
    Vector r = *mpb - (*mpA) * x;  // Initial residual
    Vector p = r;
    double rs_old = r * r;

    for (int i = 0; i < max_iter; ++i) {
        Vector Ap = (*mpA) * p;
        double alpha = rs_old / (p * Ap);

        x = x + alpha * p;
        r = r - alpha * Ap;

        double rs_new = r * r;
        if (std::sqrt(rs_new) < tol)
            break;

        p = r + (rs_new / rs_old) * p;
        rs_old = rs_new;
    }

    return x;
}
