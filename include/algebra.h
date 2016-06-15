#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "base.h"
#include <cassert>

namespace nurbs {
    
    /// Determinant of 3x3 matrix
    double det3x3(const DoubleVecVec& m);
    
    /// inverse of 3x3 matrix with given determinant
    DoubleVecVec inv3x3Mat(const DoubleVecVec& m,
                           const double det);
    
    /// Get inverse of 3x3 matrix
    DoubleVecVec inv3x3Mat(const DoubleVecVec& m);
    
    /// Get matrix of cofactors of matrix
    DoubleVecVec cofactorMat(const DoubleVecVec& m);
    
    /// Transpose a 3x3 matrix
    DoubleVecVec transpose3x3(const DoubleVecVec& m);
    
    /// Return coordinate system rotation matrix
    DoubleVecVec rotMatrix(const double theta);
    
    /// Cross product for two vectors
    template<typename T>
    std::vector<T> cross(const std::vector<T>& v1, const std::vector<T>& v2)
    {
        assert(v1.size() >= 3 && v2.size() >= 3);
        return {v1[1] * v2[2] - v1[2] * v2[1],
                v1[2] * v2[0] - v1[0] * v2[2],
                v1[0] * v2[1] - v1[1] * v2[0]};
    }
}
#endif