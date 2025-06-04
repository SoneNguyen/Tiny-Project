#pragma once

#include "Matrix.h"
#include "Vector.h"

class TikhonovSolver {
private:
    double mLambda; // Regularization parameter

public:
    // Constructor
    TikhonovSolver(double lambda);

    // Solve Ax = b using Tikhonov regularization
    Vector Solve(const Matrix& A, const Vector& b);
};
