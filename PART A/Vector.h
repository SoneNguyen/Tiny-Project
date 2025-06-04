// Vector.h
#pragma once
#include <cassert>   // for assert()
#include <iostream>  // for printing

class Vector {
private:
    int mSize;        // Size of the vector
    double* mData;    // Pointer to dynamically allocated array of doubles

public:
    // Constructor: Allocate memory for vector
    Vector(int size);

    // Copy constructor: Used when copying another Vector
    Vector(const Vector& other);

    // Destructor: Clean up memory
    ~Vector();

    // Assignment operator: Copy data from one vector to another
    Vector& operator=(const Vector& other);

    // Element access
    double& operator[](int index) const;    // 0-based indexing with assert check
    double& operator()(int index) const;    // 1-based indexing (mathematical style)

    // Operator overloads
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const; // scalar multiplication

    // Utility functions
    int Size() const;       // Returns the size
    void Print() const;     // For printing the vector
};

// dot product
double operator*(const Vector& a, const Vector& b);
