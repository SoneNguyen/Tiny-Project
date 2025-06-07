// Vector.cpp
#include "Vector.h"

// Cónstructor
Vector::Vector(){
    mSize = 0;
    mData = nullptr;
}

// Constructor: create vector of given size and zero-initialize it
Vector::Vector(int size) : mSize(size) {
    //std::cout << "[DEBUG] Creating vector of size: " << size << std::endl;
    assert(size > 0);
    mData = new double[mSize]();
}

// Copy constructor: deep copy the contents of another vector
Vector::Vector(const Vector& other) : mSize(other.mSize) {
    //std::cout << "[DEBUG] Copy constructing vector of size: " << mSize << std::endl;
    mData = new double[mSize];
    for (int i = 0; i < mSize; ++i) {
        mData[i] = other.mData[i];
    }
}

// Destructor: release memory
Vector::~Vector() {
    //std::cout << "[DEBUG] Destroying vector of size: " << mSize << std::endl;
    delete[] mData;
}

// Assignment operator: deep copy
Vector& Vector::operator=(const Vector& other) {
    if (this == &other) return *this;

    //std::cout << "[DEBUG] Assigning vector of size: " << other.mSize << std::endl;
    delete[] mData;

    mSize = other.mSize;
    mData = new double[mSize];
    for (int i = 0; i < mSize; ++i) {
        mData[i] = other.mData[i];
    }

    return *this;
}

// [] operator — 0-based indexing
double& Vector::operator[](int index) const {
    assert(index >= 0 && index < mSize);
    return mData[index];
}

// () operator — 1-based indexing
double& Vector::operator()(int index) const {
    assert(index >= 1 && index <= mSize);
    return mData[index - 1];
}

// Vector addition
Vector Vector::operator+(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result[i] = mData[i] + other.mData[i];
    }
    return result;
}

// Vector subtraction
Vector Vector::operator-(const Vector& other) const {
    assert(mSize == other.mSize);
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result[i] = mData[i] - other.mData[i];
    }
    return result;
}

// Scalar multiplication
Vector Vector::operator*(double scalar) const {
    Vector result(mSize);
    for (int i = 0; i < mSize; ++i) {
        result[i] = mData[i] * scalar;
    }
    return result;
}

// dot product
double operator*(const Vector& a, const Vector& b) {
    assert(a.Size() == b.Size());
    double sum = 0.0;
    for (int i = 0; i < a.Size(); ++i) {
        sum += a(i + 1) * b(i + 1);
    }
    return sum;
}

// Return the size of the vector
int Vector::Size() const {
    return mSize;
}

// Print the vector
void Vector::Print() const {
    //std::cout << "[";
    for (int i = 0; i < mSize; ++i) {
        std::cout << mData[i];
        if (i < mSize - 1) std::cout << ", ";
    }
    //std::cout << "]\n";
}

void Vector::resize(int size){
    if (size == mSize) {return;}

    delete[] mData;
    mSize = size;
    mData = new double[mSize];
    for (int i = 0; i < size; i++){
        mData[i] = 0;
    }
}