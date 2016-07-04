#include <utility>

#include "BoundingBoxIterator.h"
#include "Forest.h"
#include "Geometry.h"
#include "Point3D.h"
#include "NodalElement.h"

namespace nurbs {
    
    BoundingBoxIterator::BoundingBoxIterator(const Forest& f)
    :
    mCurrentIndex(0),
    mForest(f)
    {
        // First compute and store collocation data
        const Forest& fref = forest();
        
        
//        mPtData.resize(fref.geometry()->controlPtN());
//        for(uint icpt = 0; icpt < fref.geometry()->controlPtN(); ++icpt)
//            mPtData[icpt] = fref.geometry()->controlPt(icpt).asCartesian();
        
        std::vector<Point3D> pdata;
        pdata.resize(fref.collocPtN());
        std::vector<bool> cached(fref.collocPtN(), false);
        for(uint ispace = 0; ispace < fref.spaceN(); ++ispace) {
            for(uint icpt = 0; icpt < fref.collocPtN(ispace); ++icpt) {
                const uint gindex = fref.globalCollocI(ispace, icpt);
                if(gindex > pdata.size())
                    error("Bad collocation data in BoundingBoxIterator");
                if(cached[gindex])
                    continue;
                pdata[gindex] = fref.collocPt(ispace, icpt);
                cached[gindex] = true;
            }
        }
        mPtData = pdata;
        
        mPtData.resize(fref.collocPtN());
        for(uint i = 0; i < pdata.size(); ++i)
            mPtData[i] = pdata[i];
        
        // Now compute bounding boxes of basis function spans
        const double minval = std::numeric_limits<double>::lowest();
        const double maxval = std::numeric_limits<double>::max();
        
        Point3D minpt(minval, minval,minval);
        Point3D maxpt(maxval, maxval, maxval);
        
        std::vector<std::pair<Point3D, Point3D>> bbdata(fref.globalDofN(),
                                                        std::make_pair(maxpt, minpt));
        
        for(uint ispace = 0; ispace < fref.spaceN(); ++ispace) {
            for(uint ielem = 0; ielem < fref.elemN(ispace); ++ielem) {
                const auto elem = fref.element(ispace, ielem);
                const Point3D upper = elem->approxUpperBound();
                const Point3D lower = elem->approxLowerBound();
                const auto globalBasisVec = elem->globalBasisFuncI();
                for(const auto& gindex : globalBasisVec) {
                    if(gindex > fref.globalDofN())
                        error("Bad basis function index in bounding box iterator");
                    auto& currentmin = bbdata[gindex].first;
                    auto& currentmax = bbdata[gindex].second;
                    currentmin = min(currentmin, lower);
                    currentmax = max(currentmax, upper);
                }
            }
        }
        mBBData = bbdata;
    }
    
    BoundingBoxIterator& BoundingBoxIterator::operator++()
    {
        ++mCurrentIndex;
        return *this;
    }
    
    bool BoundingBoxIterator::isDone() const
    {
        return mCurrentIndex >= forest().globalDofN();
    }
    
    Point3D BoundingBoxIterator::currentPt() const
    {
        return mPtData.at(currentIndex());
    }
    
    std::pair<Point3D, Point3D> BoundingBoxIterator::currentBoundingBox() const
    {
        return mBBData.at(currentIndex());
    }
    
    Point3D BoundingBoxIterator::currentLowerBound() const
    {
        return currentBoundingBox().first;
    }
    
    /// Get current upper bound
    Point3D BoundingBoxIterator::currentUpperBound() const
    {
        return currentBoundingBox().second;
    }
    
    
    
}