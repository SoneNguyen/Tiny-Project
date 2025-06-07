#ifndef POSSYMLINSYSTEM_H
#define POSSYMLINSYSTEM_H

#include "LinearSystem.h"


class PosSymLinSystem : public LinearSystem {
public:
    PosSymLinSystem(const Matrix& A, const Vector& b);
    virtual ~PosSymLinSystem();

    // Override Solve method
    virtual Vector Solve() const override;
};

#endif // POSSYMLINSYSTEM_H


