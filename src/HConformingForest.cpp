#include "HConformingForest.h"
#include "algebra.h"

namespace nurbs {
    
    DoubleVecVec HDivForest::transformBasis(const DoubleVecVec& basis,
                                            const DoubleVecVec& jacob,
                                            const double jdet) const
    {
        // apply the transformation 1 / det(J) J^T
        // jacob = {{dx/du, dy/du, dz/du}, {dx/dv, dy/dv, dz/dv}}
        // basis = {{N1_u, N1_v}, {N2_u, N2_v}, ...}
        // result = {{B1_x, B1_y, B1_z}, {B2_x, B2_y, B2_z}, ...}
        
        DoubleVecVec result(basis.size(), DoubleVec(3, 0.0));
        for(uint b = 0; b < basis.size(); ++b)
            for(uint i = 0; i < 3; ++i)
                for(uint j = 0; j < 2; ++j)
                    result[b][i] += 1.0 / jdet * jacob[j][i] * basis[b][j];
        return result;
    }
    
    DoubleVecVec HCurlForest::transformBasis(const DoubleVecVec& basis,
                                             const DoubleVecVec& jacob,
                                             const double jdet) const
    {
        // apply the transformation 1 / det(J) J^T
        // jacob = {{dx/du, dy/du, dz/du}, {dx/dv, dy/dv, dz/dv}}
        // basis = {{N1_u, N1_v}, {N2_u, N2_v}, ...}
        // result = {{B1_x, B1_y, B1_z}, {B2_x, B2_y, B2_z}, ...}
        
        DoubleVecVec result(basis.size(), DoubleVec(3, 0.0));
        DoubleVecVec j3x3 = jacob; // create 3x3 jacobian
        j3x3.push_back({1 / jdet * (jacob[0][1] * jacob[1][2] - jacob[0][2] * jacob[1][1]),
                        1 / jdet * (jacob[0][2] * jacob[1][0] - jacob[0][0] * jacob[1][2]),
                        1 / jdet * (jacob[0][0] * jacob[1][1] - jacob[0][1] * jacob[1][0])});
        const auto jinv = inv3x3Mat(j3x3);
        for(uint b = 0; b < basis.size(); ++b)
            for(uint i = 0; i < 3; ++i)
                for(uint j = 0; j < 2; ++j)
                    result[b][i] += jinv[i][j] * basis[b][j];
        return result;
    }
}