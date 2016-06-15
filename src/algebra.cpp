#include "algebra.h"
#include <cmath>

namespace nurbs {
    
    double det3x3(const DoubleVecVec& m)
    {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }
    
    DoubleVecVec inv3x3Mat(const DoubleVecVec& m,
                        const double det)
    {
        auto result = transpose3x3(cofactorMat(m));
        for(auto& r : result)
            for(auto& c : r)
                c *= 1 / det;
        return result;
    }
    
    DoubleVecVec inv3x3Mat(const DoubleVecVec& m)
    {
        return inv3x3Mat(m, det3x3(m));
    }
    
    DoubleVecVec cofactorMat(const DoubleVecVec& m)
    {
        DoubleVecVec result(3, DoubleVec(3, 0.0));
        return { {  m[1][1] * m[2][2] - m[1][2] * m[2][1],
                  -(m[1][0] * m[2][2] - m[1][2] * m[2][0]),
                    m[1][0] * m[2][1] - m[1][1] * m[2][0] },
                 {-(m[0][1] * m[2][2] - m[0][2] * m[2][1]),
                    m[0][0] * m[2][2] - m[0][2] * m[2][0],
                  -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) },
                 {  m[0][1] * m[1][2] - m[0][2] * m[1][1],
                  -(m[0][0] * m[1][2] - m[0][2] * m[1][0]),
                    m[0][0] * m[1][1] - m[0][1] * m[1][0] } };
    }
    
    DoubleVecVec transpose3x3(const DoubleVecVec& m)
    {
        DoubleVecVec r(3, DoubleVec(3, 0.0));
        for(uint i = 0; i < 3; ++i)
            for(uint j = 0; j < 3; ++j)
                r[i][j] = m[j][i];
        return r;
    }
    
    DoubleVecVec rotMatrix(const double theta)
    {
        return {{std::cos(theta), -std::sin(theta)}, {std::sin(theta), std::cos(theta)}};
    }
    
}