#include <utility>

#include "AnalysisElement2D.h"
#include "Forest2D.h"
#include "BSplineSpace2D.h"

namespace nurbs {
    
    AnalysisElement2D::AnalysisElement2D(const Forest2D* f,
                                         const uint ispace,
                                         const uint iel)
    :
    GeometryElement2D(*f->geometry(),
                      ispace,
                      f->knotInterval(ispace, iel)),
    mForest(f),
    mElemI(iel),
    mSpan(std::make_pair(false, 0)) {}
    
    uint AnalysisElement2D::basisFuncN() const
    { return localBasisFuncI().size(); }
    
    UIntVec AnalysisElement2D::localBasisFuncI() const
    {
        if(!mLocalBasisFuncI.empty())
            return mLocalBasisFuncI;
        mLocalBasisFuncI = analysisSpace()->localBasisFuncI(lowerBound());
        return mLocalBasisFuncI;
    }
    
    UIntVec AnalysisElement2D::globalBasisFuncI() const
    {
        if(!mGlobalBasisFuncI.empty())
            return mGlobalBasisFuncI;
        for(const auto& i : localBasisFuncI())
            mGlobalBasisFuncI.emplace_back(analysisForest()->globalI(spaceI(), i));
        return mGlobalBasisFuncI;
    }
    
    DoubleVec AnalysisElement2D::basis(const double xi) const
    {
        return analysisSpace()->basis(paramCoord(xi), span());
    }
    
    /// Degree of basis
    uint AnalysisElement2D::degree() const
    { return analysisSpace()->degree(); }
    
    const Forest2D* AnalysisElement2D::analysisForest() const
    { return mForest; }
    
    const BSplineSpace2D* AnalysisElement2D::analysisSpace() const
    { return &analysisForest()->space(spaceI()); }
    
    uint AnalysisElement2D::span() const
    {
        if(mSpan.first)
            return mSpan.second;
        mSpan = std::make_pair(true, analysisSpace()->span(lowerBound()));
        return mSpan.second;
    }
    
    /// print function
    void AnalysisElement2D::print(std::ostream& ost) const
    {
        GeometryElement2D::print(ost);
        ost << "Forest: " << analysisForest() << "\n";
        ost << "Space: " << analysisSpace() << "\n";
        ost << "Local basis function indices: " << localBasisFuncI() << "\n";
        ost << "Global basis function indicies" << globalBasisFuncI() << "\n";
    }
    
}