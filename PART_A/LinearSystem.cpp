// LinearSystem.cpp
#include "LinearSystem.h"
#include <cmath>
#include <stdexcept>

LinearSystem::LinearSystem(const Matrix& A, const Vector& b) {
    assert(A.Rows() == A.Cols());
    assert(A.Rows() == b.Size());

    mSize = A.Rows();
    mpA = new Matrix(A);
    mpb = new Vector(b);
}

LinearSystem::~LinearSystem() {
    delete mpA;
    delete mpb;
}

Vector LinearSystem::Solve() const {
    Matrix A = *mpA;
    Vector b = *mpb;
    int n = mSize;

    // Forward elimination with partial pivoting
    for (int k = 0; k < n; ++k) {
        int pivot = k;
        double maxVal = std::abs(A(k + 1, k + 1));
        for (int i = k + 1; i < n; ++i) {
            double val = std::abs(A(i + 1, k + 1));
            if (val > maxVal) {
                maxVal = val;
                pivot = i;
            }
        }

        if (pivot != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(A(k + 1, j + 1), A(pivot + 1, j + 1));
            }
            std::swap(b(k + 1), b(pivot + 1));
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = A(i + 1, k + 1) / A(k + 1, k + 1);
            A(i + 1, k + 1) = 0.0;
            for (int j = k + 1; j < n; ++j) {
                A(i + 1, j + 1) -= factor * A(k + 1, j + 1);
            }
            b(i + 1) -= factor * b(k + 1);
        }
    }

    // Back substitution
    Vector x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = b(i + 1);
        for (int j = i + 1; j < n; ++j) {
            sum -= A(i + 1, j + 1) * x(j + 1);
        }
        x(i + 1) = sum / A(i + 1, i + 1);
    }

    return x;
}