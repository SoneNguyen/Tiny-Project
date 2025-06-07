// Matrix.cpp
#include "Matrix.h"
#include <cassert>
#include <iostream>
#include <cmath>

// Constructor: initialize matrix with zeros
Matrix::Matrix(int rows, int cols) : mNumRows(rows), mNumCols(cols) {
    assert(rows > 0 && cols > 0);
    //std::cout << "[DEBUG] Creating matrix: " << rows << "x" << cols << std::endl;
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; ++i) {
        mData[i] = new double[mNumCols]();
    }
}

// Copy constructor
Matrix::Matrix(const Matrix& other) : mNumRows(other.mNumRows), mNumCols(other.mNumCols) {
    std::cout << "[DEBUG] Copy constructing matrix: " << mNumRows << "x" << mNumCols << std::endl;
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; ++i) {
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; ++j) {
            mData[i][j] = other.mData[i][j];
        }
    }
}

// Destructor
Matrix::~Matrix() {
    //std::cout << "[DEBUG] Destroying matrix: " << mNumRows << "x" << mNumCols << std::endl;
    for (int i = 0; i < mNumRows; ++i) {
        delete[] mData[i];
    }
    delete[] mData;
}

// Get number of rows
int Matrix::Rows() const {
    return mNumRows;
}

// Get number of cols
int Matrix::Cols() const {
    return mNumCols;
}

// Element access (modifiable)
double& Matrix::operator()(int i, int j) {
    assert(i >= 1 && i <= mNumRows);
    assert(j >= 1 && j <= mNumCols);
    return mData[i - 1][j - 1];
}

// Element access (const version)
const double& Matrix::operator()(int i, int j) const {
    assert(i >= 1 && i <= mNumRows);
    assert(j >= 1 && j <= mNumCols);
    return mData[i - 1][j - 1];
}

// Assignment operator
Matrix& Matrix::operator=(const Matrix& other) {
    if (this == &other) return *this;

    //std::cout << "[DEBUG] Assigning matrix: " << other.mNumRows << "x" << other.mNumCols << std::endl;

    // Clean old data
    for (int i = 0; i < mNumRows; ++i) {
        delete[] mData[i];
    }
    delete[] mData;

    // Copy sizes
    mNumRows = other.mNumRows;
    mNumCols = other.mNumCols;

    // Allocate new
    mData = new double*[mNumRows];
    for (int i = 0; i < mNumRows; ++i) {
        mData[i] = new double[mNumCols];
        for (int j = 0; j < mNumCols; ++j) {
            mData[i][j] = other.mData[i][j];
        }
    }
    return *this;
}

// Matrix + Matrix
Matrix Matrix::operator+(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            result.mData[i][j] = mData[i][j] + other.mData[i][j];
        }
    }
    return result;
}

// Matrix - Matrix
Matrix Matrix::operator-(const Matrix& other) const {
    assert(mNumRows == other.mNumRows && mNumCols == other.mNumCols);
    Matrix result(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            result.mData[i][j] = mData[i][j] - other.mData[i][j];
        }
    }
    return result;
}

// Matrix * Matrix
Matrix Matrix::operator*(const Matrix& other) const {
    assert(mNumCols == other.mNumRows);
    Matrix result(mNumRows, other.mNumCols);

    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < other.mNumCols; ++j) {
            double sum = 0.0;
            for (int k = 0; k < mNumCols; ++k) {
                sum += mData[i][k] * other.mData[k][j];
            }
            result.mData[i][j] = sum;
        }
    }
    return result;
}

// Transpose
Matrix Matrix::Transpose() const {
    //std::cout << "[DEBUG] Transposing " << mNumRows << "x" << mNumCols << " matrix" << std::endl;
    Matrix result(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; ++i) {
        for (int j = 0; j < mNumCols; ++j) {
            result.mData[j][i] = mData[i][j];
        }
    }
    std::cout << "[DEBUG] Transpose complete: " << result.mNumRows << "x" << result.mNumCols << std::endl;
    return result;
}

// Determinant (recursive implementation for square matrices)
double Matrix::Determinant() const {
    assert(mNumRows == mNumCols);
    int n = mNumRows;

    if (n == 1) {
        return mData[0][0];
    }
    if (n == 2) {
        return mData[0][0] * mData[1][1] - mData[0][1] * mData[1][0];
    }

    double det = 0.0;
    for (int p = 0; p < n; ++p) {
        Matrix subMatrix(n - 1, n - 1);
        for (int i = 1; i < n; ++i) {
            int colIndex = 0;
            for (int j = 0; j < n; ++j) {
                if (j == p) continue;
                subMatrix.mData[i - 1][colIndex] = mData[i][j];
                ++colIndex;
            }
        }
        det += (p % 2 == 0 ? 1 : -1) * mData[0][p] * subMatrix.Determinant();
    }
    return det;
}

// Inverse (using adjoint method)
Matrix Matrix::Inverse() const {
    assert(mNumRows == mNumCols);
    double det = Determinant();
    assert(det != 0);

    int n = mNumRows;
    Matrix adj(n, n);

    // Compute adjoint
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Matrix subMatrix(n - 1, n - 1);
            int subi = 0;
            for (int row = 0; row < n; ++row) {
                if (row == i) continue;
                int subj = 0;
                for (int col = 0; col < n; ++col) {
                    if (col == j) continue;
                    subMatrix.mData[subi][subj] = mData[row][col];
                    ++subj;
                }
                ++subi;
            }
            adj.mData[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * subMatrix.Determinant();
        }
    }

    Matrix inv(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inv.mData[i][j] = adj.mData[i][j] / det;
        }
    }
    return inv;
}

// PseudoInverse (using formula: (A^T A)^(-1) A^T )
Matrix Matrix::PseudoInverse() const {
    Matrix At = Transpose();
    Matrix AtA = At * (*this);
    Matrix AtA_inv = AtA.Inverse();
    return AtA_inv * At;
}

// === Free function implementations ===

// Matrix × Vector multiplication
Vector operator*(const Matrix& A, const Vector& v) {
    //std::cout << "[DEBUG] Matrix-vector mult: " << A.Rows() << "x" << A.Cols() << " * " << v.Size() << std::endl;
    assert(A.Cols() == v.Size());
    Vector result(A.Rows());
    for (int i = 1; i <= A.Rows(); ++i) {
        double sum = 0.0;
        for (int j = 1; j <= A.Cols(); ++j) {
            sum += A(i, j) * v(j);
        }
        result(i) = sum;
    }
    return result;
}

// Matrix × Scalar multiplication
Matrix operator*(const Matrix& A, double scalar) {
    Matrix result(A.Rows(), A.Cols());
    for (int i = 1; i <= A.Rows(); ++i) {
        for (int j = 1; j <= A.Cols(); ++j) {
            result(i, j) = A(i, j) * scalar;
        }
    }
    return result;
}

// Scalar × Matrix multiplication
Matrix operator*(double scalar, const Matrix& A) {
    return A * scalar;
}