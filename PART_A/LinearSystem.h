#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include "Matrix.h"
#include "Vector.h"
#include <cassert>

class LinearSystem {
protected:
    int mSize;
    Matrix* mpA;
    Vector* mpb;

public:
    // Constructor: initializes the system Ax = b
    LinearSystem(const Matrix& A, const Vector& b);

    // Virtual destructor (important for base class)
    virtual ~LinearSystem();

    // Disable copy constructor and assignment operator (Rule of Three)
    LinearSystem(const LinearSystem& other) = delete;
    LinearSystem& operator=(const LinearSystem& other) = delete;

    // Solve Ax = b using Gaussian Elimination with Pivoting
    virtual Vector Solve() const;
};

#endif // LINEARSYSTEM_H
