#ifndef ANALYSIS_ELEMENT_2D_H
#define ANALYSIS_ELEMENT_2D_H

#include <cassert>
#include <iostream>

#include "GeometryElement2D.h"
#include "base.h"

namespace nurbs {
    
    /// Forward declarations
    class Forest2D;
    class BSplineSpace2D;
    
    /// A class that represents an analysis element in 2D.
    /// This is the interface that is used to obtain relevant information
    /// required for assembly (e.g. normals, jacobian determinant, basis)
    
    class AnalysisElement2D : public GeometryElement2D {
        
    public:
        
        /// Construct with a forest, bspline space index and element index
        AnalysisElement2D(const Forest2D* f,
                          const uint ispace,
                          const uint iel);
        
        /// basis function number getter
        uint basisFuncN() const;
        
        /// Local connectivity (to this parameter space)
        UIntVec localBasisFuncI() const;
        
        /// Global connectivity (of nodes in forest).
        UIntVec globalBasisFuncI() const;
        
        /// Evaluate basis function given parent coordinate xi \in [-1,1]
        DoubleVec basis(const double xi) const;
        
        /// Degree of basis
        uint degree() const;
        
        /// Order of integration
        uint integrationOrder(const uint offset = 0)
        {
            assert(offset >= 0);
            return degree() + 1;
        }
        
        /// Get the analysis forest for this element
        const Forest2D* analysisForest() const;
        
        /// Get the analysis space
        const BSplineSpace2D* analysisSpace() const;
        
        /// print function
        void print(std::ostream& ost) const;
        
    protected:
        
        /// Get the span
        uint span() const;
        
        /// Element index getter
        uint elementI() const { return mElemI; }
        
    private:
        
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const AnalysisElement2D& e)
        { e.print(ost); return ost; }
        
        /// Member data
        
        /// Pointer to forest
        const Forest2D* mForest;
        
        /// Local element index (in bspline space)
        const uint mElemI;
        
        /// Cache the local basis connectivity
        mutable UIntVec mLocalBasisFuncI;
        
        /// Cache the global basis connectivity
        mutable UIntVec mGlobalBasisFuncI;
        
        /// cache the span
        mutable std::pair<bool, uint> mSpan;
        
    };
}
#endif