// Tikhonov.cpp
#include "Tikhonov.h"

TikhonovSolver::TikhonovSolver(double lambda)
    : mLambda(lambda) {}

Vector TikhonovSolver::Solve(const Matrix& A, const Vector& b) {
    Matrix At = A.Transpose();
    Matrix AtA = At * A;

    int n = AtA.Rows();
    Matrix regularized = AtA;
    for (int i = 1; i <= n; ++i) {  //using 1-based indexing
        regularized(i, i) += mLambda * mLambda;
    }

    Matrix inv = regularized.Inverse();
    Vector Atb = At * b;
    Vector x = inv * Atb;
    return x;
}