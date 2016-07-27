#ifndef HC_ELEMENT_H
#define HC_ELEMENT_H

#include "AnalysisElement.h"
#include "MultiForest.h"

#include <iostream>
#include <memory>
#include <cassert>
#include <tuple>

namespace nurbs {

    /// An interface for representing conforming elements with a vector basis.
    /// Specifically H(curl) and H(div) elements on 3d surface manifolds for
    /// boundary element analysis. These are used specifically for
    /// electromagnetic scattering using the EFIE.
    
    /// The appropriate construction of B-spline spaces is handled by the
    /// relevant multiforest which is passed to this element during construction.
    
    /// The interface of this class
    class HCElement : public VAnalysisElement {
        
    public:
        
        /// Constructor
        HCElement(const MultiForest* mforest,
                  const uint ispace,
                  const uint iel)
        :
        AnalysisElement(*(mforest->geometry()),
                          ispace,
                          mforest->knotIntervals(ispace, iel)),
        mMultiForest(mforest),
        mElemI(iel) {}
        
        /// Virtual destructor
        virtual ~HCElement() {}
        
        /// Override copy function
        std::unique_ptr<AnalysisElement> copy() const override
        {
            return make_unique<HCElement>(*this);
        }
        
        /// Specify number of components in basis
        uint componentN() const override { return 2; }
        
        /// No. of non-zero basis functions on this element
        uint basisFuncN() const override
        {
            return localBasisFuncI().size();
        }
        
        /// Get the non-zero local basis functions over this element
        UIntVec localBasisFuncI() const override;
        
        /// Get the global basis function indices that are non-zero over this element
        UIntVec globalBasisFuncI() const override;
        
        /// Get the vector basis given parent coordinates
        DoubleVecVec basis(const double u, const double v) const override;
        
        /// Get the basis before applying the Piola transform
        DoubleVecVec localBasis(const double u, const double v) const override;
        
        /// Get the local basis derivatives at the requested parent coordinate
        DoubleVecVec localBasisDers(const double u,
                                    const double v,
                                    const DerivType dtype) const override;
        
        /// Get the degree for the vector basis component and a given parametric direction
        uint degree(const ParamDir dir, const uint comp) const override
        {
            return space(ParamDirType(comp)).degree(dir);
        }
        
        /// Get the degree vector for the given vector basis component
        UIntVec degree(const uint comp) const override
        {
            return space(ParamDirType(comp)).degree();
        }
        
        /// Get the multiforest reference
        const MultiForest* multiForest() const { return mMultiForest; }
        
        /// Local element index getter
        uint elemI() const { return mElemI; }
        
        /// Get the global element index (in the multiforest) of this element
        uint globalElemI() const { return multiForest()->globalElI(spaceI(), elemI()); }
        
        /// Get the number of connected collocation points on this element
        uint collocPtN() const override
        {
            // TODO!
            return 0;
        }
        
        /// Get the global connected collocation indices
        UIntVec globalCollocConn() const override
        {
            // TODO!
            return UIntVec();
        }
        
        /// Get the global collocation index given a local index
        uint globalCollocI(const uint icpt) const override
        {
            // TODO!
            return 0;
        }
        
        /// Get the parent coordinate of a given collocation point index
        GPt2D collocParentCoord(const uint icpt) const override
        {
            // TODO!
            return GPt2D();
        }
        
        /// Print function
        void print(std::ostream& ost) const override
        {
            GeometryElement::print(ost);
            ost << "Multiforest: " << multiForest() << "\n";
            ost << "Space index: " << spaceI() << "\n";
            ost << "Local element index: " << elemI() << "\n";
            ost << "Local basis function indices: " << localBasisFuncI() << "\n";
            ost << "Global basis function indices: " << globalBasisFuncI() << "\n";
        }
            
        protected:
            
            /// Get the relevant BSpline space for this element and parametric direction
            const BSplineSpace& space(const ParamDir d) const
            {
                return multiForest()->space(spaceI(), d);
            }
            
        private:
            
        /// Get the indices that define the span of the present element in each
        /// B-spline space of the vector basis
        UIntVec span(const ParamDir d) const;
        
        /// Reference to multiforest
        const MultiForest* mMultiForest;
        
        /// Local element index in bspline space
        const uint mElemI;
        
        /// Vector of local basis function indices that are non-zero over this element
        mutable UIntVec mLocalBasisIVec;
        
        /// Vector of global basis function indices that are non-zero over this element
        mutable UIntVec mGlobalBasisIVec;
        
        /// Vector
        mutable std::vector<Sign> mGlobalSignVec;
        
        /// Span (as required by B-spline evaluation routines) for this element
        mutable std::map<ParamDir, UIntVec> mSpan;
                
        /// Overload output operator
        friend std::ostream& operator<<(std::ostream& ost, const HCElement& e)
        { e.print(ost); return ost; }
    };
}

#endif