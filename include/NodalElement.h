//
//  Element.h
//  nurbslib
//
//  Created by Robert Simpson on 15/07/2014.
//
//

#ifndef NURBS_ELEMENT_H
#define NURBS_ELEMENT_H

#include "AnalysisElement.h"
#include "BSplineSpace.h"
#include "Forest.h"
#include "base.h"

#include <fstream>
#include <algorithm>
#include <numeric>

namespace nurbs {

    /// An element class that represents an element with nodal basis functions.
    /// An equivalent class 'VectorElement' is used to represent an element
    /// with vector basis functions.
    
    /// The class inherits from GeomElement so that all geometry related
    /// functions can be used directly (e.g. normal(), tangent(), jacob()...
    
    class NodalElement : public NAnalysisElement {
    
    public:
        
        /// Construct a nodal element from a forest, space index and local
        /// element index
        NodalElement(const Forest* f, const uint ispace, const uint iel)
        : NAnalysisElement(*f->geometry(), ispace, f->knotIntervals(ispace, iel)),
          mForest(f),
          mSpace(&f->space(ispace)),
          mElemI(iel) {}
        
        /// copy this element returning as a unique ptr. Used by
        /// copy and assignment operator of Forest class.
        std::unique_ptr<AnalysisElement> copy() const override
        {
            return make_unique<NodalElement>(*this);
        }
        
        /// Scalar basis. Only one component.
        uint componentN() const override { return 1; }
        
        /// Number of non-zero basis functions over this element.
        uint basisFuncN() const override { return localBasisFuncI().size(); }
 
        UIntVec localBasisFuncI() const override
        {
            if(!mLocalBasisFuncI.empty())
                return mLocalBasisFuncI;
            ParamCoord c = paramCoord(lowerBound(S), lowerBound(T));
            mLocalBasisFuncI = space()->globalBasisFuncI(c.s, c.t);
            return mLocalBasisFuncI;
        }
        
        UIntVec globalBasisFuncI() const override
        {
            if(!mGlobalBasisFuncI.empty())
                return mGlobalBasisFuncI;
            for(const auto& i : localBasisFuncI())
                mGlobalBasisFuncI.emplace_back(forest()->globalI(spaceI(),i));
            return mGlobalBasisFuncI;
        }
        
        DoubleVec basis(const double u, const double v) const override
        {
            ParamCoord p = paramCoord(u,v);
            return space()->basis(p.s, p.t, span());
        }
        
        /// Get the local basis derivatives
        DoubleVec localBasisDers(const double u,
                                 const double v,
                                 const DerivType dtype) const override
        {
            ParamCoord p = paramCoord(u,v);
            return space()->basisDers(p.s, p.t, span(), dtype);
        }
        
        uint degree(const ParamDir dir, const uint comp = 0) const override
        {
            return space()->degree(dir);
        }
        
        UIntVec degree(const uint comp = 0) const override
        {
            return space()->degree();
        }
        
        /// Get the number of connected collocation points on this element
        uint collocPtN() const override
        {
            return forest()->connectedCollocPtN(spaceI(), localElementI());
        }
        
        /// Get the global connected collocation indices
        UIntVec globalCollocConn() const override
        {
            return forest()->connectedGlobalCollocI(spaceI(), localElementI());
        }
        
        /// Get the global collocation index given a local index
        uint globalCollocI(const uint icpt) const override
        {
            return forest()->connectedGlobalCollocI(spaceI(), localElementI(), icpt);
        }
        
        /// Get the parent coordinate of a given collocation point index
        GPt2D collocParentCoord(const uint icpt) const override
        {
            const uint ilocal = forest()->connectedLocalCollocI(spaceI(), localElementI(), icpt);
            return parentCoord(space()->grevilleAbscissaPt(ilocal));
        }
        
        void print(std::ostream& ost) const override
        {
            GeometryElement::print(ost);
            ost << "Forest: " << forest() << "\n";
            ost << "Space: " << space() << "\n";
            ost << "Local basis function indices: " << localBasisFuncI() << "\n";
            ost << "Global basis function indicies" << globalBasisFuncI() << "\n";
        }
 
    protected:
        
        /// Return a global basis function index given a local index
        uint globalBasisFuncI(const uint l_index) const
        {
            assert(l_index < basisFuncN());
            return globalBasisFuncI()[l_index];
        }
        
        /// Forest accessor
        const Forest* forest() const {return mForest;}
        
        /// Accessor for space
        const BSplineSpace* space() const{ return mSpace; }
        
        /// Get and cache the indices defining the span of this element
        UIntVec span() const
        {
            if(!mSpan.empty())
                return mSpan;
            mSpan = space()->span(lowerBound(S), lowerBound(T));
            return mSpan;
        }
        
    private:
        
        /// Local element index getter
        uint localElementI() const { return mElemI; }
                                                
        /// Pointer to the forest
        const Forest* mForest;
        
        /// Pointer to the space which contains this element
        const BSplineSpace* mSpace;
        
        /// Store the local (to the bspline space) index for the element
        const uint mElemI;
        
        /// Cache the knot span indices for evaluating basis functions
        mutable UIntVec mLocalBasisFuncI;
        
        /// A vector containing the global basis function indices. (i.e. the global forest basis index)
        mutable UIntVec mGlobalBasisFuncI;
        
        /// Vector defining the span of this element (used by BSplineSpace
        /// evaluation)
        mutable UIntVec mSpan;
        
        friend std::ostream& operator<<(std::ostream& ost, const NodalElement& e)
        { e.print(ost); return ost; }
        
    };
    
}

#endif
