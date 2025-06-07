#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.h"
#include <cassert>   // for assert()
#include <iostream>  // for printing

class Matrix {
private:
    int mNumRows;              // Number of rows
    int mNumCols;              // Number of columns
    double** mData;            // Pointer to 2D array storing matrix data

public:
    //initialize matrix with given dimensions, set all elements to 0
    Matrix(int rows, int cols);

    //create a deep copy of another matrix
    Matrix(const Matrix& other);

    //free memory
    ~Matrix();

    //return number of rows
    int Rows() const;

    //eturn number of columns
    int Cols() const;

    // element access by get/set value at position (i, j)
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    //assignment (no constant as we are modifying *this)
    Matrix& operator=(const Matrix& other);

    //add two matrices of the same size
    Matrix operator+(const Matrix& other) const;

    //subtract two matrices of the same size
    Matrix operator-(const Matrix& other) const;

    //multiply with another matrix (compatible dimensions)
    Matrix operator*(const Matrix& other) const;

    //determinant of the matrix (must be square)
    double Determinant() const;

    //Transpose (for pseudo Inverse)
    Matrix Transpose() const;
    
    //inverse of the matrix (if invertible)
    Matrix Inverse() const;

    //Moore-Penrose pseudoinverse (for non-square or ill-conditioned matrices)
    Matrix PseudoInverse() const;
};

// Matrix-vector multiplication
Vector operator*(const Matrix& A, const Vector& v);

// Matrix × scalar multiplication
Matrix operator*(const Matrix& A, double scalar);
Matrix operator*(double scalar, const Matrix& A);  // commutative scalar multiplication

#endif